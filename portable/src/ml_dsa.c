#include "ml_dsa_core_header.h"

struct ml_dsa_keys *ml_dsa_create_keys(u8 level, ml_dsa_entropy_fn entropy)
{
	enum ml_dsa_level_k k;
	enum ml_dsa_level_l l;
	int ret;
	if(level == 4 || level == 44)
	{
		k = ML_DSA_44_K; 
		l = ML_DSA_44_L;
	}else if(level == 6 || level == 65)
	{
		k = ML_DSA_65_K; 
		l = ML_DSA_65_L;
	}else if(level == 8 || level == 87)
	{
		k = ML_DSA_87_K; 
		l = ML_DSA_87_L;
	}else { ret = ML_DSA_EINVAL; goto err_1; }
	
	s32 *t;
	t = ml_dsa_alloc(ML_DSA_N * k * sizeof(s32));
	if(!t) { ret = ML_DSA_ENOMEM; goto err_1; }
	ml_dsa_memzero(t, ML_DSA_N * k * sizeof(s32));
	
	struct ml_dsa_keys *ctx;
	ctx = ml_dsa_alloc_struct_keys(k, l);
	if(!ctx) { ret = ML_DSA_ENOMEM; goto err_2; }
	
	ret = get_rho_K_s1_s2(ctx, entropy);
	if(ret != 0) { goto err_3; }
	
	#if defined(ML_DSA_FULL_MATRIX_BUFFER)
	
	ret = get_full_matrix(ctx);
	if(ret != 0) { goto err_3; }
	
	#endif
	
	for (size_t i = 0; i < l; i++) { ml_dsa_ntt(ctx->s1 + i * ML_DSA_N); }
	
	ret = mult_matrix(ctx, ctx->s1, t);
	if(ret != 0) { goto err_3; }
	
	// TODO: слід вирішити домен для s1(мб і не тільки)
	for (size_t i = 0; i < k; i++) { ml_dsa_intt(t + i * ML_DSA_N); }
	
	for(size_t i = 0; i < k; i++)
	{
		for(size_t j = 0; j < ML_DSA_N; j++)
		{
			t[ML_DSA_N * i + j] += ctx->s2[ML_DSA_N * i + j];
		}
	}
	
	ret = create_t0_t1(ctx, t);
	if(ret != 0) { goto err_3; }
	
	ret = ml_dsa_create_public_data(ctx);
	if(ret != 0) { goto err_3; }
	
	ml_dsa_memzero(t, ML_DSA_N * k * sizeof(s32));
	ml_dsa_free(t);
	
	return ctx;
		
	err_3:
		ml_dsa_destroy_struct_keys(ctx);
	err_2:
		ml_dsa_memzero(t, ML_DSA_N * k * sizeof(s32));
		ml_dsa_free(t);
	err_1:
		return err_create_keys(ret);
}
	
	

