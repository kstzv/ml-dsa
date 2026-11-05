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
	
	struct ml_dsa_keys *ctx;
	ctx = ml_dsa_alloc_struct_keys(k, l);
	if(!ctx) { ret = ML_DSA_ENOMEM; goto err_1; }
	
	ret = get_rho_K_s1_s2(ctx, entropy);
	if(ret != 0) { goto err_2; }
	
	#if defined(ML_DSA_FULL_MATRIX_BUFFER)
	
	ret = get_full_matrix(ctx);
	if(ret != 0) { goto err_2; }
	
	#endif

	// Here should be copy s1 before NTT(it is can do in vect_y - temp)
	for (size_t i = 0; i < l; i++) { ml_dsa_ntt(ctx->s1 + i * ML_DSA_N); }
	
	mult_matrix(ctx->rho, ctx->s1, ctx->workspace->temp_vector_buffer, ctx->workspace->shake, ctx->workspace->scratch_buffer, k, l);
	
	// TODO: слід вирішити домен для s1(мб і не тільки)
	for (size_t i = 0; i < k; i++) { ml_dsa_intt(ctx->workspace->temp_vector_buffer + i * ML_DSA_N); }
	
	for(size_t i = 0; i < k; i++)
	{
		for(size_t j = 0; j < ML_DSA_N; j++)
		{
			ctx->workspace->temp_vector_buffer[ML_DSA_N * i + j] += ctx->s2[ML_DSA_N * i + j];
		}
	}
	
	ret = create_t0_t1(ctx, ctx->workspace->temp_vector_buffer);
	if(ret != 0) { goto err_2; }
	
	ret = ml_dsa_create_public_data(ctx);
	if(ret != 0) { goto err_2; }

	// Here - create functions for format pk messege and sk messege
	
	for (size_t i = 0; i < k; i++) 
	{ 
		ml_dsa_ntt(ctx->s2 + i * ML_DSA_N);
		ml_dsa_ntt(ctx->t0 + i * ML_DSA_N); 
	}
	ml_dsa_memzero(ctx->workspace->vect_y, l * ML_DSA_N * sizeof(s32));
	
	return ctx;
		
	err_2:
		ml_dsa_destroy_struct_keys(ctx);
	err_1:
		return err_create_keys(ret);
}
	
	
int ml_dsa_sign(struct ml_dsa_keys *ctx, const uint8_t *msg, size_t msg_len, const uint8_t *context, size_t context_len, uint8_t *sig, size_t sig_len, u8 prehash, u8 deterministic, ml_dsa_entropy_fn entropy)
{
	if(!ctx || !ctx->workspace || !ctx->workspace->shake) { return ML_DSA_EINVAL; }
	else if((!msg && msg_len) || !sig) { return ML_DSA_EINVAL; }
	else if(context_len > ML_DSA_SIZE_MAX_CONTEXT) { return ML_DSA_EINVAL; }
	else if(context_len != 0 && !context) { return ML_DSA_EINVAL; }
	else if(ctx->k == ML_DSA_44_K && sig_len < ML_DSA_44_SIZE_SIG) { return ML_DSA_EINVAL; }
	else if(ctx->k == ML_DSA_65_K && sig_len < ML_DSA_65_SIZE_SIG) { return ML_DSA_EINVAL; }
	else if(ctx->k == ML_DSA_87_K && sig_len < ML_DSA_87_SIZE_SIG) { return ML_DSA_EINVAL; }
	else if(deterministic > 1) { return ML_DSA_EINVAL; }
	else if(prehash > 1) { return ML_DSA_EINVAL; }
	
	int ret = 0;
	
	// Create mu
	get_mu(ctx->workspace->shake, ctx->workspace->mu, ctx->tr, msg, msg_len, context, context_len, prehash);
	ret = get_rho_double_prime(ctx, deterministic, entropy);
	if(ret != 0) { return ret; }
	
	
	u32 kappa = 0;
	u16 counter = 0;
	u16 size_z;
	if(ctx->k == ML_DSA_44_K) { size_z = ctx->l * ML_DSA_44_POLYZ_PACKED_BYTES; }
	else { size_z = ctx->l * ML_DSA_65_87_POLYZ_PACKED_BYTES; }
	while(counter <= 814)
	{
		// Get vector y
		ml_dsa_expand_mask(ctx, kappa);
		memcpy(ctx->workspace->temp_vector_buffer, ctx->workspace->vect_y, sizeof(s32) * ctx->l * ML_DSA_N);
		
		// Get vectors w and w1
		for(size_t i = 0; i < ctx->l; i++) { ml_dsa_ntt(ctx->workspace->temp_vector_buffer + i * ML_DSA_N); }
		mult_matrix(ctx->rho, ctx->workspace->temp_vector_buffer, ctx->workspace->vect_w, ctx->workspace->shake, ctx->workspace->scratch_buffer, ctx->k, ctx->l);
		for(size_t i = 0; i < ctx->k; i++) { ml_dsa_intt(ctx->workspace->vect_w + i * ML_DSA_N); }
		for(size_t i = 0; i < ctx->k; i++) { ml_dsa_canonicalize(ctx->workspace->vect_w + i * ML_DSA_N); }
		ml_dsa_get_w1(ctx->workspace->vect_w, ctx->workspace->temp_vector_buffer, NULL, ctx->k);
		
		// Get polynomial c
		u8 c[ML_DSA_64_BYTES];
		ml_dsa_get_c(ctx, c);
		ml_dsa_get_poly_c(ctx->workspace->poly_c, c, ctx->workspace->shake, ctx->k);
		ml_dsa_ntt(ctx->workspace->poly_c);
		
		// Get z
		// Обов'язково перевірити, чи не треба занулювати певні поля робочих структур для початку нової спроби
		if(ml_dsa_get_z(ctx, sig + ctx->c_tilde_bytes)) { kappa += ctx->l; counter++; continue; } 
		if(ml_dsa_get_r0(ctx)) { kappa += ctx->l; counter++; continue; }
		if(ml_dsa_get_h(ctx, sig + ctx->c_tilde_bytes + size_z)) { kappa += ctx->l; counter++; continue; }
		memcpy(sig, c, ctx->c_tilde_bytes);
		break;
	}
	
	if(counter > 814) { ret = ML_DSA_EAGAIN; }
	
	return ret;
}

int ml_dsa_verify(struct ml_dsa_verify_ctx ctx, const u8 *pk, size_t pk_len, const u8 *msg, size_t msg_len, const u8 *context, size_t context_len, const u8 *sig, size_t sig_len, u8 prehash)
{
	if(!pk || (!msg && msg_len) || (!context && context_len) || !sig || context_len > ML_DSA_SIZE_MAX_CONTEXT || prehash > 1) { return ML_DSA_EINVAL; }
	
	enum ml_dsa_level_k k;
	enum ml_dsa_level_l l;
	if(pk_len == ML_DSA_44_SIZE_PK && sig_len == ML_DSA_44_SIZE_SIG)
	{
		k = ML_DSA_44_K;
		l = ML_DSA_44_L;
	}else if(pk_len == ML_DSA_65_SIZE_PK && sig_len == ML_DSA_65_SIZE_SIG)
	{
		k = ML_DSA_65_K;
		l = ML_DSA_65_L;
	}else if(pk_len == ML_DSA_87_SIZE_PK && sig_len == ML_DSA_87_SIZE_SIG)
	{
		k = ML_DSA_87_K;
		l = ML_DSA_87_L;
	}else { return ML_DSA_EINVAL; }
	
	struct ml_dsa_verify_workspace inside_ctx;
	bool fixed_pk = 0;
	if(!ctx) 
	{ 
		inside_ctx = create_ctx_verify(enum ml_dsa_level_k k, enum ml_dsa_level_l l);
		if(!inside_ctx) { return -ML_DSA_ENOMEM; }
	}else if(ctx->size_pool == 0 || !ctx->work) { return -ML_DSA_EINVAL; }
	
	// TODO тут плануються атомарні операції із підписом, це в процесі, буде присвоєно захоплену робочу структуру із пулу для inside_ctx і якщо ключ відомий то fixed_pk = 1;
	
	int ret = ml_dsa_unpack_c_z(inside_ctx, sig);
	if(ret != 0) { goto exit; }
	for(size_t i = 0; i < ctx->l; i++) { ml_dsa_ntt(inside_ctx->temp_vector_1 + i * ML_DSA_N); }
	
	if(!fixed_pk) 
	{ 
		shake_ctx_zero(inside_ctx->shake);
		shake_ctx_init(inside_ctx->shake, inside_ctx->tr, ML_DSA_64_BYTES, pk, pk_len);
		shake256(inside_ctx->shake);
	}
	
	get_mu(inside_ctx->shake, inside_ctx->mu, inside_ctx->tr, msg, msg_len, context, context_len, prehash);
	
	#if ML_DSA_FULL_MATRIX_BUFFER
	
	if(!fixed_pk) 
	{
		mult_matrix(pk, inside_ctx->temp_vector_1, inside_ctx->temp_vector_2, inside_ctx->shake, inside_ctx->scratch_buffer, inside_ctx->k, inside_ctx->l);
	}else
	{
		ret = mult_matrix(inside_ctx->matrix_buffer, inside_ctx->temp_vector_1, inside_ctx->temp_vector_2, inside_ctx->k, inside_ctx->l);
		if(ret != 0) { goto exit; }
	}
	
	#else
	
	mult_matrix(pk, inside_ctx->temp_vector_1, inside_ctx->temp_vector_2, inside_ctx->shake, inside_ctx->scratch_buffer, inside_ctx->k, inside_ctx->l);
	
	#endif
	
	if(!fixed_pk)
	{
		ml_dsa_unpack_pk(inside_ctx->temp_vector_1, pk, inside_ctx->k);
		for(size_t i = 0; i < inside_ctx->k; i++) { ml_dsa_ntt(inside_ctx->temp_vector_1 + i * ML_DSA_N); }
	}
	
	ml_dsa_get_poly_c(inside_ctx->poly_c, inside_ctx->c_tilde, inside_ctx->shake, inside_ctx->k);
	ml_dsa_ntt(inside_ctx->poly_c);
	
	for(size_t i = 0; i < inside_ctx->k; i++)
	{
		for(size_t j = 0; j < ML_DSA_N; j++)
		{
			inside_ctx->temp_vector_1[j + i * ML_DSA_N] = ml_dsa_montgomery_reduce((s64)inside_ctx->poly_c[j] * inside_ctx->temp_vector_1[j + i * ML_DSA_N]);
			inside_ctx->temp_vector_2[j + i * ML_DSA_N] -= inside_ctx->temp_vector_1[j + i * ML_DSA_N];
		}
	}
	
	for(size_t i = 0; i < inside_ctx->k; i++) { ml_dsa_intt(inside_ctx->temp_vector_2 + i * ML_DSA_N); } // get 𝐰Approx
	ret = ml_dsa_unpack_h(inside_ctx, sig);
	if(ret != 0) { goto exit; }
	for(size_t i = 0; i < ML_DSA_N * inside_ctx->k; i++) { inside_ctx->temp_vector_2[i] = ml_dsa_freeze(inside_ctx->temp_vector_2[i]); }
	ml_dsa_get_w1(inside_ctx->temp_vector_2, inside_ctx->temp_vector_1, inside_ctx->scratch_buffer, inside_ctx->k);
	ret = ml_dsa_pack_w1(inside_ctx);
	if(ret != 0) { goto exit; }
	
	ret = ml_dsa_get_c_prime(inside_ctx);
		
	exit:
		if(!ctx) { remove_ctx_verify(inside_ctx); }
		else
		{
			// TODO - тут буде функція атомарного звільнення слоту
		}
		
		return ret;
}
	
	
	
	
	
	
	
	










