#include "ml_dsa_core_header.h"

static const uint8_t oid_shake256[11] = {
    0x06, 0x09, 0x60, 0x86, 0x48, 0x01,
    0x65, 0x03, 0x04, 0x02, 0x0C
};

static inline void ml_dsa_decompose(s32 *r, s32 gamma2, s32 *r1, s32 *r0);
static inline size_t pack_for_c_44(struct ml_dsa_keys *ctx);
static inline size_t pack_for_c_65_87(struct ml_dsa_keys *ctx);
static inline void pack_z_44(s32 *coefs, u8 *buffer);
static inline void pack_z_65_87(s32 *coefs, u8 *buffer);
static inline void bit_unpack_z_44(struct ml_dsa_verify_ctx *ctx, const u8 *ptr_z);
static inline void bit_unpack_z_65_87(struct ml_dsa_verify_ctx *ctx, const u8 *ptr_z);
static inline int ml_dsa_unpack_h(struct ml_dsa_verify_ctx *ctx, const u8 *ptr_h, const u8 *k, u8 k_size, u8 omega);

int get_rho_double_prime(struct ml_dsa_keys *ctx, u8 deterministic, ml_dsa_entropy_fn entropy);
void ml_dsa_expand_mask(struct ml_dsa_keys *ctx, u32 kappa);
void ml_dsa_get_w1(struct ml_dsa_keys *ctx);
void ml_dsa_get_c(struct ml_dsa_keys *ctx, u8 *c);
void ml_dsa_get_poly_c(struct ml_dsa_keys *ctx, u8 *c);
bool ml_dsa_get_z(struct ml_dsa_keys *ctx, u8 *out);
bool ml_dsa_get_r0(struct ml_dsa_keys *ctx);
bool ml_dsa_get_h(struct ml_dsa_keys *ctx, u8 *out);
int ml_dsa_unpack_c_z(struct ml_dsa_verify_ctx *ctx, const u8 *sig_msg);

void get_mu(struct shake_ctx *shake, u8 *mu, u8 *tr, const uint8_t *msg, size_t msg_len, const uint8_t *context, size_t context_len, u8 prehash)
{
	// buffer for first three rate (136 * 3 = 408)
	u8 first_rate[408]; 
	ml_dsa_memzero(first_rate, 408);
	
	// copy hash pk, bytes prehash and size of contest
	memcpy(first_rate, tr, ML_DSA_64_BYTES); 
	first_rate[ML_DSA_64_BYTES] = prehash;        
	first_rate[ML_DSA_64_BYTES + 1] = (u8)context_len;
	
	// copy context
	size_t counter_ptr = ML_DSA_64_BYTES + 2;
	if(context_len) { memcpy(first_rate + counter_ptr, context, context_len); }
	counter_ptr += context_len;
	
	if(prehash == 1)
	{
		shake_ctx_zero(shake);
		shake_ctx_init(shake, mu, ML_DSA_64_BYTES, msg, msg_len);
		shake256(shake);
		memcpy(first_rate + counter_ptr, oid_shake256, ML_DSA_PREHASH_OID_SIZE);
		counter_ptr += ML_DSA_PREHASH_OID_SIZE;
		memcpy(first_rate + counter_ptr, mu, ML_DSA_64_BYTES);
		counter_ptr += ML_DSA_64_BYTES;
		
		shake_ctx_zero(shake);
		shake_ctx_init(shake, mu, ML_DSA_64_BYTES, first_rate, counter_ptr);
		shake256(shake);
		goto exit;
	}
	
	if(408 - counter_ptr >= msg_len) 
	{
		if(msg_len) { memcpy(first_rate + counter_ptr, msg, msg_len); }
		counter_ptr += msg_len;
		
		shake_ctx_zero(shake);
		shake_ctx_init(shake, mu, ML_DSA_64_BYTES, first_rate, counter_ptr);
		shake256(shake);
		goto exit;
	}
	
	memcpy(first_rate + counter_ptr, msg, 408 - counter_ptr); 
	shake_ctx_zero(shake);
	shake_ctx_init(shake, mu, ML_DSA_64_BYTES, first_rate, 408);
	shake->rate = 136;
	shake_absorb(shake);
	counter_ptr = 408 - counter_ptr;
	shake->in = msg + counter_ptr;
	shake->inlen = msg_len - counter_ptr;
	shake256(shake);
	
	exit:
		ml_dsa_memzero(first_rate, 408);
		shake_ctx_zero(shake);
		return;
}

int get_rho_double_prime(struct ml_dsa_keys *ctx, u8 deterministic, ml_dsa_entropy_fn entropy)
{
	// Copy K in first 32 bytes
	memcpy(ctx->workspace->scratch_buffer, ctx->K, ML_DSA_32_BYTES);
	
	// Get next 32 zeroies or random bytes
	if (deterministic == 1) { ml_dsa_memzero(ctx->workspace->scratch_buffer + ML_DSA_32_BYTES, ML_DSA_32_BYTES); }
	else if(entropy == NULL)
	{
		if (ml_dsa_entropy(ctx->workspace->scratch_buffer + ML_DSA_32_BYTES, ML_DSA_32_BYTES) != 0) { return ML_DSA_SYSTEM_ENTROPY_FAILED; }
	}else
	{
		if (entropy(ctx->workspace->scratch_buffer + ML_DSA_32_BYTES, ML_DSA_32_BYTES) != 0) { return ML_DSA_CALLBACK_ENTROPY_FAILED; }
	}
	
	// Copy finish 64 mu bytes
	memcpy(ctx->workspace->scratch_buffer + ML_DSA_64_BYTES, ctx->workspace->mu, ML_DSA_64_BYTES);
	
	// Get rho``
	shake_ctx_zero(ctx->workspace->shake);
	shake_ctx_init(ctx->workspace->shake, ctx->workspace->rho_double_prime, ML_DSA_64_BYTES, ctx->workspace->scratch_buffer, ML_DSA_64_BYTES * 2);
	shake256(ctx->workspace->shake);
	ml_dsa_memzero(ctx->workspace->scratch_buffer, ML_DSA_64_BYTES * 2);
	
	return 0;
}

void ml_dsa_expand_mask(struct ml_dsa_keys *ctx, u32 kappa)
{
	u32 gamma1;
    u32 coefficient_bits;
    size_t shake_outlen;
    
    if (ctx->l == ML_DSA_44_L) 
    { 
		gamma1 = ML_DSA_44_GAMMA1; 
		coefficient_bits = ML_DSA_44_GAMMA1_BITS; 
		shake_outlen = ML_DSA_44_POLYZ_PACKED_BYTES;
	}else if(ctx->l == ML_DSA_65_L || ctx->l == ML_DSA_87_L)
	{
		gamma1 = ML_DSA_65_87_GAMMA1;
		coefficient_bits = ML_DSA_65_87_GAMMA1_BITS;
		shake_outlen = ML_DSA_65_87_POLYZ_PACKED_BYTES;
	}
	
	u8 *buffer = ctx->workspace->scratch_buffer;
	for (u32 r = 0; r < (u32)ctx->l; ++r)
	{
		u32 nonce = kappa + r;
		for (size_t i = 0; i < ML_DSA_64_BYTES; ++i)
		{
			buffer[i] = ctx->workspace->rho_double_prime[i];
		}
		
		buffer[64] = (u8)nonce;
        buffer[65] = (u8)(nonce >> 8);
        
        shake_ctx_zero(ctx->workspace->shake);
        shake_ctx_init(ctx->workspace->shake, buffer, shake_outlen, buffer, 66);
        shake256(ctx->workspace->shake);
        
        u32 accumulator = 0;
        u32 bits_in_accumulator = 0;
        size_t byte_position = 0;
        u32 mask = (1U << coefficient_bits) - 1U;
        for (size_t i = 0; i < ML_DSA_N; ++i)
        {
			while (bits_in_accumulator < coefficient_bits)
			{
				accumulator |= (u32)buffer[byte_position++] << bits_in_accumulator;
				bits_in_accumulator += 8;
			}
			
			u32 value = accumulator & mask;
			accumulator >>= coefficient_bits;
            bits_in_accumulator -= coefficient_bits;
            ctx->workspace->vect_y[(size_t)r * ML_DSA_N + i] = (s32)gamma1 - (s32)value;
		}
	}
}


void ml_dsa_get_w1(const s32 *w, s32 *w1, const u8 *h, enum ml_dsa_level_k k)
{
	s32 m, gamma2;
	if(k == ML_DSA_44_K) { gamma2 = ML_DSA_44_GAMMA2; m = 44; } // For 44 level gamma2 and m = (q-1)/(2*𝛾2)
	else { gamma2 = ML_DSA_65_87_GAMMA2; m = 16; } // For 65 and 87 levels gamma2 and m = (q-1)/(2*𝛾2)

    for (size_t i = 0; i < k * ML_DSA_N; i++)
    {
        s32 r = w[i];
        s32 high, low;

        ml_dsa_decompose(&r, gamma2, &high, &low);

        if (h != NULL && h[i] == 1)
        {
            if (low > 0)
            {
                high = (high == m - 1) ? 0 : high + 1;
            }
            else
            {
                high = (high == 0) ? m - 1 : high - 1;
            }
        }

        w1[i] = high;
    }
}


static inline void ml_dsa_decompose(s32 *r, s32 gamma2, s32 *r1, s32 *r0)
{
    s32 high = (*r + 127) >> 7;

    if (gamma2 == ML_DSA_65_87_GAMMA2) 
    {
        high = (high * 1025 + (1 << 21)) >> 22;
        high &= 15;
    } else 
    {
        high = (high * 11275 + (1 << 23)) >> 24;
        high ^= ((43 - high) >> 31) & high;
    }

    s32 low = *r - high * 2 * gamma2;
    low -= ((ML_DSA_Q_HALF_MINUS_ONE - low) >> 31) & ML_DSA_Q;

    if (r1 != NULL) { *r1 = high; }
	if (r0 != NULL) { *r0 = low; }
}

void ml_dsa_get_c(struct ml_dsa_keys *ctx, u8 *c)
{
	size_t counter_bytes = 0;
	memcpy(ctx->workspace->scratch_buffer, ctx->workspace->mu, ML_DSA_64_BYTES);
	counter_bytes += ML_DSA_64_BYTES;
	
	if(ctx->k == ML_DSA_44_K) { counter_bytes += pack_for_c_44(ctx); }
	else { counter_bytes += pack_for_c_65_87(ctx); }
	
	size_t size_c = 0;
	if(ctx->k == ML_DSA_44_K) { size_c = ML_DSA_32_BYTES; }
	else if(ctx->k == ML_DSA_65_K) { size_c = ML_DSA_32_BYTES + 16; }
	else if(ctx->k == ML_DSA_87_K) { size_c = ML_DSA_64_BYTES; }
	
	shake_ctx_zero(ctx->workspace->shake);
	shake_ctx_init(ctx->workspace->shake, c, size_c, ctx->workspace->scratch_buffer, counter_bytes);
	shake256(ctx->workspace->shake);
}
	
			
static inline size_t pack_for_c_44(struct ml_dsa_keys *ctx)
{
    size_t counter_bytes = 0;
    u8 *temp_buff = ctx->workspace->scratch_buffer + ML_DSA_64_BYTES;
    const s32 *w = ctx->workspace->temp_vector_buffer;

    for (size_t i = 0; i < (size_t)ctx->k * ML_DSA_N; i += 4)
    {
        temp_buff[counter_bytes++] = (u8)(w[i] | (w[i + 1] << 6));
        temp_buff[counter_bytes++] = (u8)((w[i + 1] >> 2) | (w[i + 2] << 4));
        temp_buff[counter_bytes++] = (u8)((w[i + 2] >> 4) | (w[i + 3] << 2));
    }
    return counter_bytes;
}

static inline size_t pack_for_c_65_87(struct ml_dsa_keys *ctx)
{
    size_t counter_bytes = 0;
    u8 *temp_buff = ctx->workspace->scratch_buffer + ML_DSA_64_BYTES;
    const s32 *w = ctx->workspace->temp_vector_buffer;

    for (size_t i = 0; i < (size_t)ctx->k * ML_DSA_N; i += 2)
    {
        temp_buff[counter_bytes++] = (u8)(w[i] | (w[i + 1] << 4));
    }
    return counter_bytes;
}

//  SampleInBall
void ml_dsa_get_poly_c(s32 *poly_c, u8 *c, struct shake_ctx *shake, enum ml_dsa_level_k k)
{
	ml_dsa_memzero(poly_c, ML_DSA_N * sizeof(s32));
	
	// get τ
	u8 tau, size_c;
	if(k == ML_DSA_44_K) { tau = 39; size_c = 32; }
	else if(k == ML_DSA_65_K) { tau = 49; size_c = 48; }
	else if(k == ML_DSA_87_K) { tau = 60; size_c = 64; }

	// 136 - size of rate in shake256
	u8 temp_buff[136];
	shake_ctx_zero(shake);
	shake_ctx_init(shake, temp_buff, 136, c, size_c);
	shake256(shake);
	
	// First 8 bytes contain sign bits
	u8 bit_sing[8];
	memcpy(bit_sing, temp_buff, 8);
	
	size_t pos = 8;

	for (u16 i = 256 - tau; i < 256; i++)
	{
		u8 j;
		// Rejection sampling: 0 <= j <= i
        do {
            if (pos == sizeof(temp_buff)) { shake256(shake); pos = 0; }
            j = temp_buff[pos++];
        } while (j > i);

		// Fisher-Yates
		poly_c[i] = poly_c[j];

		// h[i + tau - 256]
		u8 h_index = i + tau - 256;

		u8 h = (bit_sing[h_index >> 3] >> (h_index & 7)) & 1;

		poly_c[j] = 1 - ((s32)h << 1);
	}
}

bool ml_dsa_get_z(struct ml_dsa_keys *ctx, u8 *out)
{
	bool ret = 0;
	
	s32 gamma;
	s32 beta;
	size_t bytes_counter;
	void (*pack_z)(s32 *, u8 *);
	s32 temp_mass[4];
	u8 *temp_out = out;
	if(ctx->k == ML_DSA_44_K) 
	{ 
		gamma = ML_DSA_44_GAMMA1;
		beta = ML_DSA_44_BETA;
		bytes_counter = 9;
		pack_z = pack_z_44;
		 
	}else 
	{ 
		gamma = ML_DSA_65_87_GAMMA1; 
		bytes_counter = 10;
		pack_z = pack_z_65_87;
		if(ctx->k == ML_DSA_65_K) { beta = ML_DSA_65_BETA; }
		else { beta = ML_DSA_87_BETA; }
	}
	
	for(size_t i = 0; i < ctx->l; i++)
		{
			for(size_t j = 0; j < ML_DSA_N; j++)
			{
				ctx->workspace->temp_vector_buffer[j + i * ML_DSA_N] = ml_dsa_montgomery_reduce((s64)ctx->workspace->poly_c[j] * (s64)ctx->s1[j + i * ML_DSA_N]);
			}
		}
		
		for(size_t i = 0; i < ctx->l; i++) { ml_dsa_intt(ctx->workspace->temp_vector_buffer + i * ML_DSA_N); }
		
		for(size_t i = 0; i < ctx->l; i++)
		{
			size_t j = 0;
			while(j < ML_DSA_N)
			{
				ctx->workspace->vect_y[j + i * ML_DSA_N] += ctx->workspace->temp_vector_buffer[j + i * ML_DSA_N];
				temp_mass[0] = gamma - ctx->workspace->vect_y[j + i * ML_DSA_N];
				ret |= (ctx->workspace->vect_y[j + i * ML_DSA_N] >= gamma - beta || ctx->workspace->vect_y[j + i * ML_DSA_N] <= -(gamma - beta));
				j++;
				
				ctx->workspace->vect_y[j + i * ML_DSA_N] += ctx->workspace->temp_vector_buffer[j + i * ML_DSA_N];
				temp_mass[1] = gamma - ctx->workspace->vect_y[j + i * ML_DSA_N];
				ret |= (ctx->workspace->vect_y[j + i * ML_DSA_N] >= gamma - beta || ctx->workspace->vect_y[j + i * ML_DSA_N] <= -(gamma - beta));
				j++;
				
				ctx->workspace->vect_y[j + i * ML_DSA_N] += ctx->workspace->temp_vector_buffer[j + i * ML_DSA_N];
				temp_mass[2] = gamma - ctx->workspace->vect_y[j + i * ML_DSA_N];
				ret |= (ctx->workspace->vect_y[j + i * ML_DSA_N] >= gamma - beta || ctx->workspace->vect_y[j + i * ML_DSA_N] <= -(gamma - beta));
				j++;
				
				ctx->workspace->vect_y[j + i * ML_DSA_N] += ctx->workspace->temp_vector_buffer[j + i * ML_DSA_N];
				temp_mass[3] = gamma - ctx->workspace->vect_y[j + i * ML_DSA_N];
				ret |= (ctx->workspace->vect_y[j + i * ML_DSA_N] >= gamma - beta || ctx->workspace->vect_y[j + i * ML_DSA_N] <= -(gamma - beta));
				j++;
				
				pack_z(temp_mass, temp_out);
				temp_out += bytes_counter;
			}
		}
		
		return ret;
}
				
static inline void pack_z_44(s32 *coefs, u8 *buffer)
{
	u8 i, j;
	i = 0;
	j = 0;
	buffer[i]   = (u8)coefs[j]; 
	buffer[++i] = (u8)((u32)coefs[j] >> 8); 
	buffer[++i] = (u8)((u32)coefs[j] >> 16);
	buffer[i]  |= (u8)((u32)coefs[++j] << 2);
	buffer[++i] = (u8)((u32)coefs[j] >> 6);
	buffer[++i] = (u8)((u32)coefs[j] >> 14);
	buffer[i]  |= (u8)((u32)coefs[++j] << 4);
	buffer[++i] = (u8)((u32)coefs[j] >> 4);
	buffer[++i] = (u8)((u32)coefs[j] >> 12);
	buffer[i]  |= (u8)((u32)coefs[++j] << 6);
	buffer[++i] = (u8)((u32)coefs[j] >> 2);
	buffer[++i] = (u8)((u32)coefs[j] >> 10);
}

static inline void pack_z_65_87(s32 *coefs, u8 *buffer)
{
	u8 i, j;
	i = 0;
	j = 0;
	buffer[i]   = (u8)coefs[j]; 
	buffer[++i] = (u8)((u32)coefs[j] >> 8);
	buffer[++i] = (u8)((u32)coefs[j] >> 16);
	buffer[i]  |= (u8)((u32)coefs[++j] << 4);
	buffer[++i] = (u8)((u32)coefs[j] >> 4);
	buffer[++i] = (u8)((u32)coefs[j] >> 12);
	buffer[++i] = (u8)coefs[++j]; 
	buffer[++i] = (u8)((u32)coefs[j] >> 8);
	buffer[++i] = (u8)((u32)coefs[j] >> 16);
	buffer[i]  |= (u8)((u32)coefs[++j] << 4);
	buffer[++i] = (u8)((u32)coefs[j] >> 4);
	buffer[++i] = (u8)((u32)coefs[j] >> 12);
}

bool ml_dsa_get_r0(struct ml_dsa_keys *ctx)
{
	bool ret = 0;
	
	s32 gamma;
	s32 beta;
	s32 temp = 0;
	if(ctx->k == ML_DSA_44_K)
	{
		gamma = ML_DSA_44_GAMMA2;
		beta = ML_DSA_44_BETA;
	}else if(ctx->k == ML_DSA_65_K)
	{
		gamma = ML_DSA_65_87_GAMMA2;
		beta = ML_DSA_65_BETA;
	}else if(ctx->k == ML_DSA_87_K)
	{
		gamma = ML_DSA_65_87_GAMMA2;
		beta = ML_DSA_87_BETA;
	}
	
	for(size_t i = 0; i < ctx->k; i++)
	{
		for(size_t j = 0; j < ML_DSA_N; j++)
		{
			ctx->workspace->temp_vector_buffer[j + i * ML_DSA_N] = ml_dsa_montgomery_reduce((s64)ctx->workspace->poly_c[j] * (s64)ctx->s2[j + i * ML_DSA_N]);
		}
	}
	
	for(size_t i = 0; i < ctx->k; i++) { ml_dsa_intt(ctx->workspace->temp_vector_buffer + i * ML_DSA_N); }
	
	for(size_t i = 0; i < ctx->k; i++)
	{
		for(size_t j = 0; j < ML_DSA_N; j++)
		{
			ctx->workspace->vect_w[j + i * ML_DSA_N] -= ctx->workspace->temp_vector_buffer[j + i * ML_DSA_N];
			ml_dsa_decompose(&ctx->workspace->vect_w[j + i * ML_DSA_N], gamma, NULL, &temp);
			ret |= ((temp >= gamma - beta) || (temp <= -(gamma - beta)));
		}
	}
	
	return ret;
}

// Check: mb cache leak
bool ml_dsa_get_h(struct ml_dsa_keys *ctx, u8 *out)
{
	bool ret = 0;
	u8 hint_end_offsets[8]; // 8 - max value for k, in ML-DSA-87
	
	// Mult t0 and c
	for(u8 i = 0; i < ctx->k; i++)
	{
		for(size_t j = 0; j < ML_DSA_N; j++)
		{
			ctx->workspace->temp_vector_buffer[j + i * ML_DSA_N] = ml_dsa_montgomery_reduce((s64)ctx->workspace->poly_c[j] * (s64)ctx->t0[j + i * ML_DSA_N]);
		}
	}
	
	// iNTT for c * t0
	for(u8 i = 0; i < ctx->k; i++) { ml_dsa_intt(ctx->workspace->temp_vector_buffer + i * ML_DSA_N); }
	
	s32 gamma2;
	u8 omega;
	if(ctx->k == ML_DSA_44_K) 
	{
		gamma2 = ML_DSA_44_GAMMA2;
		omega = ML_DSA_OMEGA_44;
	}else if(ctx->k == ML_DSA_65_K)
	{
		gamma2 = ML_DSA_65_87_GAMMA2;
		omega = ML_DSA_OMEGA_65;
	}else if(ctx->k == ML_DSA_87_K)
	{
		gamma2 = ML_DSA_65_87_GAMMA2;
		omega = ML_DSA_OMEGA_87;
	}
	
	ml_dsa_memzero(out, (size_t)omega + ctx->k);
	
	// MakeHint
	s32 z = 0;
	s32 r = 0;
	u8 omega_counter = 0;
	for(u8 i = 0; i < ctx->k; i++)
	{
		for(size_t j = 0; j < ML_DSA_N; j++)
		{
			r = ctx->workspace->vect_w[j + i * ML_DSA_N] + ctx->workspace->temp_vector_buffer[j + i * ML_DSA_N];
			z = r + ((-1) * ctx->workspace->temp_vector_buffer[j + i * ML_DSA_N]);
			
			r = ml_dsa_freeze(r);
			z = ml_dsa_freeze(z);
			
			ml_dsa_decompose(&r, gamma2, &r, NULL);
			ml_dsa_decompose(&z, gamma2, &z, NULL);
			
			// Check diapazone gamma
			ret |= ((ctx->workspace->temp_vector_buffer[j + i * ML_DSA_N] >= gamma2) || (ctx->workspace->temp_vector_buffer[j + i * ML_DSA_N] <= -gamma2));
			u8 pos_one = (r != z) * j; // value get 0 if r != z
			out[omega_counter] += pos_one; // value and counter == 0 both, so old value don`t because + 0
			omega_counter += (r != z); // counter don`t change if r != z
			ret |= omega_counter > omega;  // if counter more than omega ret get 1 and this functions return 1
			
			omega_counter *= !ret; // if counter >= omega, so it try is wrong. To avoid overflowing the provided buffer, the counter should always be multiplied by 0
		}
		hint_end_offsets[i] = omega_counter; // Write current value omega
	}
	
	memcpy(out + omega, hint_end_offsets, ctx->k);
	return ret;
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

 /* The alternative variours for unpack t1:
             t0 = (u32)in[0] | (((u32)in[1] & 0x03) << 8);
			 t1 = ((u32)in[1] >> 2) | (((u32)in[2] & 0x0F) << 6);
			 t2 = ((u32)in[2] >> 4) | (((u32)in[3] & 0x3F) << 4);
			 t3 = ((u32)in[3] >> 6) | ((u32)in[4] << 2);
			 */	

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
	
	ml_dsa_memzero(ctx->scratch_buffer, (size_t)ctx->k * ML_DSA_N); // TODO - визначити розміри scratch_buffer при створенні модуля геенрації памяті
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

static inline s32 ml_dsa_freeze(s32 val)
{
    val = ml_dsa_barrett_reduce(val);
    val += (val >> 31) & ML_DSA_Q;
    return val;
}
			

	
	
	
