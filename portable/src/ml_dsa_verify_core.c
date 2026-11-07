// A module of internal functions used only in verification - ml_dsa_verify_core.c
#include "ml_dsa_core_header.h"

int ml_dsa_unpack_c_z(struct ml_dsa_verify_workspace *ctx, const u8 *sig_msg);
void ml_dsa_unpack_pk(s32 * t1, const u8 *pk_msg, enum ml_dsa_level_k k);
int ml_dsa_unpack_h(struct ml_dsa_verify_workspace *ctx, u8 *sign);
int ml_dsa_pack_w1(struct ml_dsa_verify_workspace *ctx);
int ml_dsa_get_c_prime(struct ml_dsa_verify_workspace *ctx);

static inline int bit_unpack_z_65_87(struct ml_dsa_verify_workspace *ctx, const u8 *ptr_z);
static inline int bit_unpack_z_44(struct ml_dsa_verify_workspace *ctx, const u8 *ptr_z);

int ml_dsa_unpack_c_z(struct ml_dsa_verify_workspace *ctx, const u8 *sig_msg)
{
	if(ctx->k == ML_DSA_44_K) 
	{ 
		memcpy(ctx->c_tilde, sig_msg, ML_DSA_32_BYTES);
		return bit_unpack_z_44(ctx, sig_msg + ML_DSA_32_BYTES);
	}else if(ctx->k == ML_DSA_65_K) 
	{ 
		memcpy(ctx->c_tilde, sig_msg, (ML_DSA_32_BYTES + 16));
		return bit_unpack_z_65_87(ctx, sig_msg + (ML_DSA_32_BYTES + 16));
	}else if(ctx->k == ML_DSA_87_K) 
	{ 
		memcpy(ctx->c_tilde, sig_msg, ML_DSA_64_BYTES);
		return bit_unpack_z_65_87(ctx, sig_msg + ML_DSA_64_BYTES);
	}else { return -ML_DSA_EINVAL; }
}

static inline int bit_unpack_z_44(struct ml_dsa_verify_workspace *ctx, const u8 *ptr_z)
{
    const u32 upper = 2 * ML_DSA_44_GAMMA1 - ML_DSA_44_BETA;

    for (u8 i = 0; i < ctx->l; i++)
    {
        size_t coef_ctr = 0;
        s32 *out = ctx->temp_vector_1 + i * ML_DSA_N;

        while (coef_ctr < ML_DSA_N)
        {
            u32 t[4];

            t[0] = (u32)ptr_z[0] | ((u32)ptr_z[1] << 8) | (((u32)ptr_z[2] & 0x03) << 16);

            t[1] = ((u32)ptr_z[2] >> 2) | ((u32)ptr_z[3] << 6) | (((u32)ptr_z[4] & 0x0F) << 14);

            t[2] = ((u32)ptr_z[4] >> 4) | ((u32)ptr_z[5] << 4) | (((u32)ptr_z[6] & 0x3F) << 12);

            t[3] = ((u32)ptr_z[6] >> 6) | ((u32)ptr_z[7] << 2) | ((u32)ptr_z[8] << 10);

            for (size_t j = 0; j < 4; j++)
            {
                if (t[j] <= ML_DSA_44_BETA || t[j] >= upper) { return -ML_DSA_EBADSIG; }
                out[coef_ctr + j] = (s32)ML_DSA_44_GAMMA1 - (s32)t[j];
            }

            coef_ctr += 4;
            ptr_z += 9;
        }
    }

    return 0;
}

static inline int bit_unpack_z_65_87(struct ml_dsa_verify_workspace *ctx, const u8 *ptr_z)
{
	s32 upper;
	if(ctx->k == ML_DSA_65_K) { upper = (s32)ML_DSA_65_87_GAMMA1 - (s32)ML_DSA_65_BETA; }
	else if(ctx->k == ML_DSA_87_K) { upper = (s32)ML_DSA_65_87_GAMMA1 - (s32)ML_DSA_87_BETA; }
	else { return -ML_DSA_EINVAL; }
	
	for(u8 i = 0; i < ctx->l; i++)
	{
		size_t coef_ctr = 0;
		s32 temp;
		while(coef_ctr < ML_DSA_N)
		{
			temp = ptr_z[0];
			temp |= (ptr_z[1] << 8);
			temp |= ((ptr_z[2] & 0x0F) << 16);
			
			ctx->temp_vector_1[coef_ctr + i * ML_DSA_N] = (s32)ML_DSA_65_87_GAMMA1 - temp;
			if((ctx->temp_vector_1[coef_ctr + i * ML_DSA_N] <= -upper) || (ctx->temp_vector_1[coef_ctr + i * ML_DSA_N] >= upper)) { return -ML_DSA_EBADSIG; }
			coef_ctr++;
			
			temp = (ptr_z[2] >> 4);
			temp |= (ptr_z[3] << 4);
			temp |= (ptr_z[4] << 12);
			
			ctx->temp_vector_1[coef_ctr + i * ML_DSA_N] = (s32)ML_DSA_65_87_GAMMA1 - temp;
			if((ctx->temp_vector_1[coef_ctr + i * ML_DSA_N] <= -upper) || (ctx->temp_vector_1[coef_ctr + i * ML_DSA_N] >= upper)) { return -ML_DSA_EBADSIG; }
			coef_ctr++;
			
			ptr_z += 5;
		}
	}
	return 0;
}

void ml_dsa_unpack_pk(s32 * t1, const u8 *pk_msg, enum ml_dsa_level_k k)
{
    // memcpy(ctx->rho, pk_msg, ML_DSA_32_BYTES);

    const u8 *ptr = pk_msg + ML_DSA_32_BYTES;

    for(u8 i = 0; i < k; i++)
    {
        size_t j = 0;
        while(j < ML_DSA_N)
        {
            t1[j + i * ML_DSA_N] = (s32)(*ptr);
            ptr++;
            t1[j + i * ML_DSA_N] |= (s32)(((*ptr) & 0x03) << 8);
            
            t1[j + i * ML_DSA_N] <<= ML_DSA_D;
            j++;

            t1[j + i * ML_DSA_N] = (s32)(*ptr >> 2);
            ptr++;
            t1[j + i * ML_DSA_N] |= (s32)(((*ptr) & 0x0F) << 6);
            
            t1[j + i * ML_DSA_N] <<= ML_DSA_D;
            j++;

            t1[j + i * ML_DSA_N] = (s32)(*ptr >> 4);
            ptr++;
            t1[j + i * ML_DSA_N] |= (s32)(((*ptr) & 0x3F) << 4);
            
            t1[j + i * ML_DSA_N] <<= ML_DSA_D;
            j++;

            t1[j + i * ML_DSA_N] = (s32)(*ptr >> 6);
            ptr++;
            t1[j + i * ML_DSA_N] |= (s32)((*ptr) << 2);
            
            t1[j + i * ML_DSA_N] <<= ML_DSA_D;
            j++;
            ptr++;
        }
    }
}

int ml_dsa_unpack_h(struct ml_dsa_verify_workspace *ctx, u8 *sign)
{
	if(!ctx || !sign) { return -ML_DSA_EINVAL; }
	
	u8 omega;
	u8 *ptr_h;
	u8 *k;
	if(ctx->k == ML_DSA_44_K) 
	{ 
		omega = ML_DSA_OMEGA_44; 
		ptr_h = sign + (ML_DSA_44_SIZE_SIG - (ML_DSA_OMEGA_44 + ctx->k)); 
		k = sign + (ML_DSA_44_SIZE_SIG - ctx->k); 
	}else if(ctx->k == ML_DSA_65_K) 
	{ 
		omega = ML_DSA_OMEGA_65;
		ptr_h = sign + (ML_DSA_65_SIZE_SIG - (ML_DSA_OMEGA_65 + ctx->k)); 
		k = sign + (ML_DSA_65_SIZE_SIG - ctx->k); 
	}else if(ctx->k == ML_DSA_87_K) 
	{ 
		omega = ML_DSA_OMEGA_87; 
		ptr_h = sign + (ML_DSA_87_SIZE_SIG - (ML_DSA_OMEGA_87 + ctx->k)); 
		k = sign + (ML_DSA_87_SIZE_SIG - ctx->k); 
	}else { return -ML_DSA_EINVAL; }
	
	ml_dsa_memzero(ctx->scratch_buffer, ML_DSA_SIZE_SCRATCH_BUFFER); // TODO - визначити розміри scratch_buffer при створенні модуля геенрації памяті
    u8 units_ctr = 0;
    for(u8 i = 0; i < ctx->k; i++)
    {
        if(k[i] > omega || k[i] < units_ctr) { return -ML_DSA_EINVAL; }

        u8 first = units_ctr;
        while(units_ctr < k[i])
        {
            if(units_ctr > first && ptr_h[units_ctr - 1] >= ptr_h[units_ctr]) { return -ML_DSA_EINVAL; }

            ctx->scratch_buffer[ptr_h[units_ctr] + i * ML_DSA_N] = 1;
            units_ctr++;
        }
    }
    
    while (units_ctr < omega)
    {
		if(ptr_h[units_ctr] != 0) { return -ML_DSA_EINVAL; }
		units_ctr++;
	}

    return 0;
}

int ml_dsa_pack_w1(struct ml_dsa_verify_workspace *ctx)
{
	if(!ctx) { return -ML_DSA_EINVAL; }
	
	if(ctx->k == ML_DSA_44_K)
	{
		size_t coef_counter = 0;
		size_t bytes_counter = 0;
		while(coef_counter < ctx->k * ML_DSA_N)
		{
			ctx->scratch_buffer[bytes_counter++] = (u8)(ctx->temp_vector_1[coef_counter] | (ctx->temp_vector_1[coef_counter + 1] << 6));
			coef_counter++;
			ctx->scratch_buffer[bytes_counter++] = (u8)((ctx->temp_vector_1[coef_counter] >> 2) | (ctx->temp_vector_1[coef_counter + 1] << 4));
			coef_counter++;
			ctx->scratch_buffer[bytes_counter++] = (u8)((ctx->temp_vector_1[coef_counter] >> 4) | (ctx->temp_vector_1[coef_counter + 1] << 2));
			coef_counter += 2;
		}
	}else if(ctx->k == ML_DSA_65_K || ctx->k == ML_DSA_87_K)
	{
		size_t coef_counter = 0;
		size_t bytes_counter = 0;
		while(coef_counter < ctx->k * ML_DSA_N)
		{
			ctx->scratch_buffer[bytes_counter++] = (u8)(ctx->temp_vector_1[coef_counter] | (ctx->temp_vector_1[coef_counter + 1] << 4));
			coef_counter += 2;
		}
	}else { return -ML_DSA_EINVAL; }
	
	return 0;
}

int ml_dsa_get_c_prime(struct ml_dsa_verify_workspace *ctx)
{
	if(!ctx) { return -ML_DSA_EINVAL; }
	
	size_t size_w1_bytes;
	size_t size_lambda;
    if(ctx->k == ML_DSA_44_K) { size_w1_bytes = 768; size_lambda = ML_DSA_32_BYTES; }
	else if(ctx->k == ML_DSA_65_K) { size_w1_bytes = 768; size_lambda = ML_DSA_32_BYTES + 16; }
	else if(ctx->k == ML_DSA_87_K) { size_w1_bytes = 1024; size_lambda = ML_DSA_64_BYTES; }
	else { return -ML_DSA_EINVAL; }
	
	u8 first_rate[136];
	memcpy(first_rate, ctx->mu, ML_DSA_64_BYTES);
	memcpy(first_rate + ML_DSA_64_BYTES, ctx->scratch_buffer, 136 - ML_DSA_64_BYTES);
	
	shake_ctx_zero(shake);
	shake_ctx_init(shake, ctx->mu, size_lambda, first_rate, 136);
	shake->rate = 136;
	shake_absorb(shake);
	shake->in = ctx->scratch_buffer + (136 - ML_DSA_64_BYTES);
	shake->inlen = size_w1_bytes - (136 - ML_DSA_64_BYTES);
	shake256(shake);
	
	for(u8 i = 0; i < size_lambda; i++)
	{
		if(ctx->mu[i] != ctx->c_tilde[i]) { return -ML_DSA_EBADSIG; }
	}
	
	return 0;
}
