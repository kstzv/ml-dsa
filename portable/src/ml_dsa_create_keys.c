// ML-DSA key generation
// This implementation is incomplete and under active development
#include "ml_dsa_core_header.h"

int ml_dsa_create_struct_keys(struct ml_dsa_keys *ctx, ml_dsa_entropy_fn entropy);
struct ml_dsa_keys *ml_dsa_alloc_struct_keys(ml_dsa_level_k k, ml_dsa_level_l l);
void ml_dsa_destroy_struct_keys(struct ml_dsa_keys *ctx);

// The main function of key generation
int ml_dsa_create_struct_keys(struct ml_dsa_keys *ctx, ml_dsa_entropy_fn entropy)
{
	if(!ctx) { return ML_DSA_EINVAL; }
	// TODO
	
}
	
// Memory allocation for the key structure and all its fields
struct ml_dsa_keys *ml_dsa_alloc_struct_keys(ml_dsa_level_k k, ml_dsa_level_l l)
{
	// Memory allocation for key structure
	struct ml_dsa_keys *ctx;
	ctx = ml_dsa_alloc(sizeof(struct ml_dsa_keys));
	if(!ctx) { goto err_1; }
	
	ctx->k = k;
	ctx->l = l;
	
	// Parameter η definition
	if(ctx->k == ML_DSA_44_K && ctx->l == ML_DSA_44_L) { ctx->eta = ML_DSA_44_87_ETA; }
	else if(ctx->k == ML_DSA_65_K && ctx->l == ML_DSA_65_L) { ctx->eta = ML_DSA_65_ETA; }
	else if(ctx->k == ML_DSA_87_K && ctx->l == ML_DSA_87_L) { ctx->eta = ML_DSA_44_87_ETA; }
	
	// Allocation for fields s1 and s2
	ctx->s1 = ml_dsa_alloc(l * ML_DSA_N + k * ML_DSA_N);
	if(!ctx->s1) { goto err_2; }
	ctx->s2 = ctx->s1 + l * ML_DSA_N;
	
	// Allocation for fields t0 and t1
	ctx->t0 = ml_dsa_alloc(2 * k * ML_DSA_N);
	if(!ctx->t0) { goto err_3; }
	ctx->t1 = ctx->t0 + k * ML_DSA_N;
	
	return ctx;
		
	err_3:
		ml_dsa_free(ctx->s1);
	err_2:
		ml_dsa_free(ctx);
	err_1:
		return NULL;
}

// Key destructor, secret reset, and memory cleanup
void ml_dsa_destroy_struct_keys(struct ml_dsa_keys *ctx)
{
	if(!ctx) { return; }
	
	if(ctx->s1) 
	{
		ml_dsa_memzero(ctx->s1, ctx->l * ML_DSA_N + ctx->k * ML_DSA_N);
		ml_dsa_free(ctx->s1); 
	}
	
	if(ctx->t0) 
	{ 
		ml_dsa_memzero(ctx->t0, 2 * ctx->k * ML_DSA_N);
		ml_dsa_free(ctx->t0);
	}
	
	 ml_dsa_memzero(ctx->K, ML_DSA_32_BYTES);
	 ml_dsa_memzero(ctx->rho, ML_DSA_32_BYTES);
	 ml_dsa_memzero(ctx->tr, ML_DSA_64_BYTES);
	 
	 ml_dsa_free(ctx); 
}

// The function for generating a primary seed and obtaining values 
// ​For the following key structure fields: ρ, K, s1 and s2
int get_rho_K_s1_s2(struct ml_dsa_keys *ctx, ml_dsa_entropy_fn entropy)
{
	// Buffers for first and expanded seeds 
	u8 first_seed[ML_DSA_32_BYTES + 2];
	u8 out_first_seed[ML_DSA_32_BYTES * 2 + ML_DSA_64_BYTES];
	
	// Check and use entropy
	if(entropy == NULL)
	{
		if(ml_dsa_entropy(first_seed, ML_DSA_32_BYTES) != 0) { return ML_DSA_SYSTEM_ENTROPY_FAILED; }
	}else{
		if(entropy(first_seed, ML_DSA_32_BYTES) != 0) { return ML_DSA_CALLBACK_ENTROPY_FAILED; }
	}
	
	// Set k and l
	first_seed[32] = ctx->k;
	first_seed[33] = ctx->l;
	
	// Get and write in struct K and ρ parameters
	ml_dsa_shake256(out_first_seed, ML_DSA_32_BYTES * 2 + ML_DSA_64_BYTES, first_seed, ML_DSA_32_BYTES + 2);
	memcpy(ctx->rho, out_first_seed, ML_DSA_32_BYTES);
	memcpy(ctx->K, out_first_seed + ML_DSA_32_BYTES + ML_DSA_64_BYTES, ML_DSA_32_BYTES);
	
	// Create and ready buffer for seed s1 and s2 vectors
	u8 s1_s2_seed[ML_DSA_64_BYTES + 2];
	s1_s2_seed[65] = 0;
	memcpy(s1_s2_seed, out_first_seed + ML_DSA_32_BYTES, ML_DSA_64_BYTES);
	
	// Determining the limit for selecting coefficients by parameter η
	// And introducing a buffer of 384 bytes, with the aim of a high 
	// Probability of obtaining the necessary bytes for one polynomial at a time
	u8 limit;
	u8 buffer_shake[384];
	if(ctx->eta == ML_DSA_44_87_ETA) { limit = 15; }
	else { limit = 9; }
	
	// Basic cycle for vectors s1 and s2
	for(u8 i = 0; i < ctx->k + ctx->l; i++)
	{
		// Gets bytes from SHAKE256
		s1_s2_seed[64] = i;
		ml_dsa_shake256(buffer_shake, 384, s1_s2_seed, ML_DSA_64_BYTES + 2);
		
		// Set pointer
		s8 *ptr;
		if(i < ctx->l) { ptr = ctx->s1 + ML_DSA_N * i; }
		else { ptr = ctx->s2 + (i - ctx->l) * ML_DSA_N; }
		
		// Half-byte selection according to algorithm 15 - CoeffFromHalfByte
		size_t count_coef, count_buff;
		count_coef = 0;
		count_buff = 0;
		while(count_coef < ML_DSA_N)
		{
			u8 val_0;
			u8 val_1;
			val_0 = buffer_shake[count_buff] & 0x0F;
			val_1 = buffer_shake[count_buff] >> 4;
			
			if(ctx->eta == ML_DSA_44_87_ETA && val_0 < limit)
			{
				val_0 = val_0 - (205 * val_0 >> 10) * 5;
				ptr[count_coef] = 2 - val_0;
				count_coef++;
			} else if(ctx->eta == ML_DSA_65_ETA && val_0 < limit)
			{
				ptr[count_coef] = 4 - val_0;
				count_coef++;
			}
			
			if(count_coef == ML_DSA_N) { break; }
			
			if(ctx->eta == ML_DSA_44_87_ETA && val_1 < limit)
			{
				val_1 = val_1 - (205 * val_1 >> 10) * 5;
				ptr[count_coef] = 2 - val_1;
				count_coef++;
			} else if(ctx->eta == ML_DSA_65_ETA && val_1 < limit)
			{
				ptr[count_coef] = 4 - val_1;
				count_coef++;
			}
			
			count_buff++;
			if(count_buff >= 384) { return ML_DSA_EAGAIN; }
		}
	}
	
	// Clearing static arrays of functions
	ml_dsa_memzero(first_seed, ML_DSA_32_BYTES + 2);
	ml_dsa_memzero(out_first_seed, ML_DSA_32_BYTES * 2 + ML_DSA_64_BYTES);
	ml_dsa_memzero(s1_s2_seed, ML_DSA_64_BYTES + 2);
	ml_dsa_memzero(buffer_shake, 384);
	
	return 0;
}
	
	
	
	
	
	













