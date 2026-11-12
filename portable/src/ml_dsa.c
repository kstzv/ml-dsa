#include "ml_dsa_core_header.h"

struct ml_dsa_keys *ml_dsa_create_keys(u8 level, ml_dsa_entropy_fn entropy);
struct ml_dsa_keys *ml_dsa_create_keys_object_from_imported_sk(u8 *sk_buffer, size_t size_buffer, u8 level);
void ml_dsa_destroy_struct_keys(struct ml_dsa_keys *ctx);

int ml_dsa_sign(struct ml_dsa_keys *ctx, const u8 *msg, size_t msg_len, const u8 *context, size_t context_len, u8 *sig, size_t sig_len, u8 prehash, u8 deterministic, ml_dsa_entropy_fn entropy);
int ml_dsa_verify(struct ml_dsa_verify_ctx *ctx, const u8 *pk, size_t pk_len, const u8 *msg, size_t msg_len, const u8 *context, size_t context_len, const u8 *sig, size_t sig_len, u8 prehash);

int ml_dsa_get_pk(struct ml_dsa_keys *ctx, u8 *buffer, size_t size_buffer);
int ml_dsa_get_sk(struct ml_dsa_keys *ctx, u8 *buffer, size_t size_buffer);

struct ml_dsa_verify_ctx *ml_dsa_alloc_verify_pool(u8 level, size_t number_slots, u8 *pk, size_t size_pk);
void destructor_verify_pool(struct ml_dsa_verify_ctx *ctx);

static inline size_t unpack_t0(s32 *t0, const u8 *buffer, u8 k);
static inline size_t pack_t0(const s32 *t0, u8 *buffer, u8 k);
static inline size_t pack_s1_and_s2_65(const s32 *s, u8 *buffer, u8 lng_loops);
static inline size_t pack_s1_and_s2_44_87(const s32 *s, u8 *buffer, u8 lng_loops);
static inline size_t unpack_s1_and_s2_65(s32 *s, const u8 *buffer, u8 lng_loops);
static inline size_t unpack_s1_and_s2_44_87(s32 *s, const u8 *buffer, u8 lng_loops);

static inline int get_ctx_keys_from_seed(struct ml_dsa_keys *ctx);
static inline int get_ctx_keys_from_sk(struct ml_dsa_keys *ctx, u8 *sk_buffer);

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
	
	ret = get_rho_K_s1_s2(ctx, entropy, 0, NULL, 0);
	if(ret != 0) { goto err_2; }

	ret = get_ctx_keys_from_seed(ctx);
	if(ret != 0) {  goto err_2; }
	
	ml_dsa_memzero(ctx->workspace->temp_vector_buffer, (ctx->l + 2 * ctx->k + 1) * ML_DSA_N * sizeof(s32));
	
	return ctx;
		
	err_2:
		ml_dsa_destroy_struct_keys(ctx);
	err_1:
		return err_create_keys(-ret);
}

// Key destructor, secret reset, and memory cleanup
void ml_dsa_destroy_struct_keys(struct ml_dsa_keys *ctx)
{
	if(!ctx) { return; }
	
	if(ctx->s1) 
	{
		ml_dsa_memzero(ctx->s1, (ctx->l + ctx->k) * ML_DSA_N * sizeof(s32));
		ml_dsa_free(ctx->s1); 
	}
	
	if(ctx->t0) 
	{ 
		ml_dsa_memzero(ctx->t0, (ctx->k * ML_DSA_N) * sizeof(s32));
		ml_dsa_free(ctx->t0);
	}
	
	if(ctx->pk) { ml_dsa_free(ctx->pk); }
	
	if(ctx->workspace && ctx->workspace->temp_vector_buffer)
	{
		size_t size_s32 = ctx->k * ML_DSA_N * sizeof(s32); 		  
		size_s32 += ctx->k * ML_DSA_N * sizeof(s32); 		 
		size_s32 += ctx->l * ML_DSA_N * sizeof(s32); 		 
		size_s32 += ML_DSA_N * sizeof(s32); 		
		ml_dsa_memzero(ctx->workspace->temp_vector_buffer, size_s32);
		ml_dsa_free(ctx->workspace->temp_vector_buffer);
	}
	
	if(ctx->workspace && ctx->workspace->scratch_buffer)
	{
		ml_dsa_memzero(ctx->workspace->scratch_buffer, ML_DSA_SIZE_SCRATCH_BUFFER);
		ml_dsa_free(ctx->workspace->scratch_buffer);
	}
	
	if(ctx->workspace && ctx->workspace->shake)
	{
		if(ctx->workspace->shake->state)
		{
			shake_ctx_zero(ctx->workspace->shake);
			ml_dsa_free(ctx->workspace->shake->state);
		}
		ml_dsa_free(ctx->workspace->shake);
	}
	 
	 ml_dsa_memzero(ctx->workspace, sizeof(struct ml_dsa_workspace));
	 ml_dsa_free(ctx->workspace);
	 ml_dsa_memzero(ctx, sizeof(struct ml_dsa_keys));
	 ml_dsa_free(ctx); 
}
	
	
int ml_dsa_sign(struct ml_dsa_keys *ctx, const u8 *msg, size_t msg_len, const u8 *context, size_t context_len, u8 *sig, size_t sig_len, u8 prehash, u8 deterministic, ml_dsa_entropy_fn entropy)
{
	if(!ctx || !ctx->workspace || !ctx->workspace->shake) { return err_ml_dsa_sign(-ML_DSA_EINVAL); }
	else if((!msg && msg_len) || !sig) { return err_ml_dsa_sign(-ML_DSA_EINVAL); }
	else if(context_len > ML_DSA_SIZE_MAX_CONTEXT) { return err_ml_dsa_sign(-ML_DSA_EINVAL); }
	else if(context_len != 0 && !context) { return err_ml_dsa_sign(-ML_DSA_EINVAL); }
	else if(ctx->k == ML_DSA_44_K && sig_len < ML_DSA_44_SIZE_SIG) { return err_ml_dsa_sign(-ML_DSA_EINVAL); }
	else if(ctx->k == ML_DSA_65_K && sig_len < ML_DSA_65_SIZE_SIG) { return err_ml_dsa_sign(-ML_DSA_EINVAL); }
	else if(ctx->k == ML_DSA_87_K && sig_len < ML_DSA_87_SIZE_SIG) { return err_ml_dsa_sign(-ML_DSA_EINVAL); }
	else if(deterministic > 1) { return err_ml_dsa_sign(-ML_DSA_EINVAL); }
	else if(prehash > 1) { return err_ml_dsa_sign(-ML_DSA_EINVAL); }
	
	int ret = 0;
	
	// Create mu
	get_mu(ctx->workspace->shake, ctx->workspace->mu, ctx->tr, msg, msg_len, context, context_len, prehash);
	ret = get_rho_double_prime(ctx, deterministic, entropy);
	if(ret != 0) { return err_ml_dsa_sign(-ret); }
	
	
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
		if(ml_dsa_get_z(ctx, sig + ctx->c_tilde_bytes)) { kappa += ctx->l; counter++; ml_dsa_memzero(c, ML_DSA_64_BYTES); continue; } 
		if(ml_dsa_get_r0(ctx)) { kappa += ctx->l; counter++; ml_dsa_memzero(c, ML_DSA_64_BYTES); continue; }
		if(ml_dsa_get_h(ctx, sig + ctx->c_tilde_bytes + size_z)) { kappa += ctx->l; counter++; ml_dsa_memzero(c, ML_DSA_64_BYTES); continue; }
		memcpy(sig, c, ctx->c_tilde_bytes);
		ml_dsa_memzero(c, ML_DSA_64_BYTES);
		break;
	}
	
	if(counter > 814) 
	{ 
		ml_dsa_wipper_sign_struct(ctx->workspace, ctx->k, ctx->l); 
		ml_dsa_memzero(sig, sig_len);
		return err_ml_dsa_sign(-ML_DSA_EAGAIN); 
	}
	
	return 0;
}

int ml_dsa_verify(struct ml_dsa_verify_ctx *ctx, const u8 *pk, size_t pk_len, const u8 *msg, size_t msg_len, const u8 *context, size_t context_len, const u8 *sig, size_t sig_len, u8 prehash)
{
	if((!msg && msg_len) || (!context && context_len) || !sig || context_len > ML_DSA_SIZE_MAX_CONTEXT || prehash > 1) { return err_ml_dsa_verify(-ML_DSA_EINVAL); }
	
	if(sig_len != ML_DSA_44_SIZE_SIG && sig_len != ML_DSA_65_SIZE_SIG && sig_len != ML_DSA_87_SIZE_SIG) { return err_ml_dsa_verify(-ML_DSA_EINVAL); }
	
	enum ml_dsa_level_k k;
	enum ml_dsa_level_l l;
	if((!ctx && pk && pk_len != 0) || (ctx && !ctx->fixed_pk && pk && pk_len !=0))
	{
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
		}else { return err_ml_dsa_verify(-ML_DSA_EINVAL); }
	}else if(ctx && ctx->fixed_pk && ctx->work && !pk && pk_len == 0)
	{
		k = ctx->work->k;
		l = ctx->work->l;
		
		if(sig_len == ML_DSA_44_SIZE_SIG && (k != ML_DSA_44_K || l != ML_DSA_44_L)) { return err_ml_dsa_verify(-ML_DSA_EINVAL); }
		if(sig_len == ML_DSA_65_SIZE_SIG && (k != ML_DSA_65_K || l != ML_DSA_65_L)) { return err_ml_dsa_verify(-ML_DSA_EINVAL); }
		if(sig_len == ML_DSA_87_SIZE_SIG && (k != ML_DSA_87_K || l != ML_DSA_87_L)) { return err_ml_dsa_verify(-ML_DSA_EINVAL); }
		
	}else { return err_ml_dsa_verify(-ML_DSA_EINVAL); }
	
	// Validate verification pool
	if(ctx)
	{
		if(ctx->size_pool == 0 || !ctx->work) { return err_ml_dsa_verify(-ML_DSA_EINVAL); }
		if(ctx->work[0].k != k || ctx->work[0].l != l) { return err_ml_dsa_verify(-ML_DSA_EINVAL); }
	}
	
	// Initialize local verification state
	struct ml_dsa_verify_workspace *inside_ctx = NULL;
	bool fixed_pk = 0;

	if(!ctx)
	{
		inside_ctx = ml_dsa_alloc_verify_works_struct(k, l, 1);
		if(!inside_ctx) { return err_ml_dsa_verify(-ML_DSA_ENOMEM); }
	}else
	{
		for(size_t i = 0; i < ctx->size_pool; i++)
		{
			if(ml_dsa_try_acquire_slot(&ctx->work[i].is_free))
			{
				inside_ctx = &ctx->work[i];
				if(ctx->fixed_pk) { fixed_pk = 1; } 
				break;
			}
		}
	}

    if(!inside_ctx) { return err_ml_dsa_verify(-ML_DSA_EAGAIN); }
		
	int ret = ml_dsa_unpack_c_z(inside_ctx, sig);
	if(ret != 0) { goto err; }
	for(size_t i = 0; i < l; i++) { ml_dsa_ntt(inside_ctx->temp_vector_1 + i * ML_DSA_N); }
	
	if(!fixed_pk) 
	{ 
		shake_ctx_zero(inside_ctx->shake);
		shake_ctx_init(inside_ctx->shake, inside_ctx->tr, ML_DSA_64_BYTES, pk, pk_len);
		shake256(inside_ctx->shake);
	}
	
	get_mu(inside_ctx->shake, inside_ctx->mu, inside_ctx->tr, msg, msg_len, context, context_len, prehash);
	
	#if ML_DSA_FULL_MATRIX_BUFFER

	if(fixed_pk)
	{
		ret = mult_full_matrix(ctx->matrix_buffer, inside_ctx->temp_vector_1, inside_ctx->temp_vector_2, inside_ctx->k, inside_ctx->l);
		if(ret != 0) { goto err; }
	}else
	{
		mult_matrix((u8 *)pk, inside_ctx->temp_vector_1, inside_ctx->temp_vector_2, inside_ctx->shake, inside_ctx->scratch_buffer, inside_ctx->k, inside_ctx->l);
	}

	#else
	
	const u8 *matrix_seed;
	if(fixed_pk) { matrix_seed = inside_ctx->rho; } 
	else { matrix_seed = pk; }

	mult_matrix(matrix_seed, inside_ctx->temp_vector_1, inside_ctx->temp_vector_2, inside_ctx->shake, inside_ctx->scratch_buffer, inside_ctx->k, inside_ctx->l);

	#endif
	
	ml_dsa_get_poly_c(inside_ctx->poly_c, inside_ctx->c_tilde, inside_ctx->shake, inside_ctx->k);
	ml_dsa_ntt(inside_ctx->poly_c);
	
	if(!fixed_pk)
	{
		ml_dsa_unpack_pk(inside_ctx->temp_vector_1, pk, inside_ctx->k);
		for(size_t i = 0; i < inside_ctx->k; i++) { ml_dsa_ntt(inside_ctx->temp_vector_1 + i * ML_DSA_N); }
	
		for(size_t i = 0; i < inside_ctx->k; i++)
		{
			for(size_t j = 0; j < ML_DSA_N; j++)
			{
				inside_ctx->temp_vector_1[j + i * ML_DSA_N] = ml_dsa_montgomery_reduce((s64)inside_ctx->poly_c[j] * inside_ctx->temp_vector_1[j + i * ML_DSA_N]);
				inside_ctx->temp_vector_2[j + i * ML_DSA_N] -= inside_ctx->temp_vector_1[j + i * ML_DSA_N];
			}
		}
	
	}else
	{
		for(size_t i = 0; i < inside_ctx->k; i++)
		{
			for(size_t j = 0; j < ML_DSA_N; j++)
			{
				inside_ctx->temp_vector_1[j + i * ML_DSA_N] = ml_dsa_montgomery_reduce((s64)inside_ctx->poly_c[j] * ctx->t1[j + i * ML_DSA_N]);
				inside_ctx->temp_vector_2[j + i * ML_DSA_N] -= inside_ctx->temp_vector_1[j + i * ML_DSA_N];
			}
		}
	}
	
	for(size_t i = 0; i < inside_ctx->k; i++) { ml_dsa_intt(inside_ctx->temp_vector_2 + i * ML_DSA_N); } // get 𝐰Approx
	ret = ml_dsa_unpack_h(inside_ctx, sig);
	if(ret != 0) { goto err; }
	for(size_t i = 0; i < ML_DSA_N * inside_ctx->k; i++) { inside_ctx->temp_vector_2[i] = ml_dsa_freeze(inside_ctx->temp_vector_2[i]); }
	ml_dsa_get_w1(inside_ctx->temp_vector_2, inside_ctx->temp_vector_1, inside_ctx->scratch_buffer, inside_ctx->k);
	ret = ml_dsa_pack_w1(inside_ctx);
	if(ret != 0) { goto err; }
	
	// Тут слід буде продумати вихід та позначення 0 чи 1 в ситуації коли був отримано результат збігу чи ні
	ret = ml_dsa_get_c_prime(inside_ctx);
		
	err:
		if(!ctx) { destructor_verify_works_struct(inside_ctx, 1); }
		else { ml_dsa_release_slot(&inside_ctx->is_free); }
		
		return err_ml_dsa_verify(ret);
}
	
// Створення пулу верифікацій
struct ml_dsa_verify_ctx *ml_dsa_alloc_verify_pool(u8 level, size_t number_slots, u8 *pk, size_t size_pk)
{
	enum ml_dsa_level_k k;
	enum ml_dsa_level_l l;
	int ret = 0;
	
	if(level == 44 || level == 4) { k = ML_DSA_44_K; l = ML_DSA_44_L; }
	else if(level == 65 || level == 6) { k = ML_DSA_65_K; l = ML_DSA_65_L; }
	else if(level == 87 || level == 8) { k = ML_DSA_87_K; l = ML_DSA_87_L; }
	else { return err_ml_dsa_alloc_verify_pool(-ML_DSA_EINVAL); }
	
	if(k == ML_DSA_44_K && pk && size_pk != ML_DSA_44_SIZE_PK) { return err_ml_dsa_alloc_verify_pool(-ML_DSA_EINVAL); }
	if(k == ML_DSA_65_K && pk && size_pk != ML_DSA_65_SIZE_PK) { return err_ml_dsa_alloc_verify_pool(-ML_DSA_EINVAL); }
	if(k == ML_DSA_87_K && pk && size_pk != ML_DSA_87_SIZE_PK) { return err_ml_dsa_alloc_verify_pool(-ML_DSA_EINVAL); }
	
	if(number_slots == 0 || number_slots > MAX_NUMBER_SLOTS) { return err_ml_dsa_alloc_verify_pool(-ML_DSA_EINVAL); }
	
	struct ml_dsa_verify_workspace *ctx;
	ctx = ml_dsa_alloc_verify_works_struct(k, l, number_slots); 
	if(!ctx) { return err_ml_dsa_alloc_verify_pool(-ML_DSA_ENOMEM); }
	
	struct ml_dsa_verify_ctx *pool_ctx;
	pool_ctx = ml_dsa_alloc(sizeof(struct ml_dsa_verify_ctx));
	if(!pool_ctx) { ret = ML_DSA_ENOMEM; goto err_1; }
	ml_dsa_memzero(pool_ctx, sizeof(struct ml_dsa_verify_ctx));
	
	pool_ctx->work = ctx;
	pool_ctx->size_pool = number_slots;
	
	if(pk)
	{
		if(k == ML_DSA_44_K && size_pk != ML_DSA_44_SIZE_PK) { ret = ML_DSA_EINVAL; goto err_2; }
		if(k == ML_DSA_65_K && size_pk != ML_DSA_65_SIZE_PK) { ret = ML_DSA_EINVAL; goto err_2; }
		if(k == ML_DSA_87_K && size_pk != ML_DSA_87_SIZE_PK) { ret = ML_DSA_EINVAL; goto err_2; }
		
		size_t size_mem_alloc = k * ML_DSA_N * sizeof(s32);
		#if ML_DSA_FULL_MATRIX_BUFFER
		size_mem_alloc += k * l * ML_DSA_N * sizeof(s32);
		#endif
		
		pool_ctx->t1 = ml_dsa_alloc(size_mem_alloc);
		if(!pool_ctx->t1) { ret = ML_DSA_ENOMEM; goto err_2; }
		ml_dsa_unpack_pk(pool_ctx->t1, pk, k);
		for(size_t i = 0; i < k; i++) { ml_dsa_ntt(pool_ctx->t1 + i * ML_DSA_N); }
		
		#if ML_DSA_FULL_MATRIX_BUFFER
		pool_ctx->matrix_buffer = pool_ctx->t1 + k * ML_DSA_N;
		get_full_matrix(pk, k, l, pool_ctx->work->shake, pool_ctx->matrix_buffer, pool_ctx->work->scratch_buffer);
		#else
		for(size_t i = 0; i < number_slots; i++) { memcpy(pool_ctx->work[i].rho, pk, ML_DSA_32_BYTES); }
		#endif
		
		shake_ctx_zero(pool_ctx->work->shake);
		shake_ctx_init(pool_ctx->work->shake, pool_ctx->work->tr, ML_DSA_64_BYTES, pk, size_pk);
		shake256(pool_ctx->work->shake);
		shake_ctx_zero(pool_ctx->work->shake);
		
		for(size_t i = 1; i < number_slots; i++) { memcpy(pool_ctx->work[i].tr, pool_ctx->work->tr, ML_DSA_64_BYTES); }
		
		pool_ctx->fixed_pk = 1;
	}else { pool_ctx->t1 = NULL; pool_ctx->fixed_pk = 0; }
	
	for(size_t i = 0; i < number_slots; i++) { ml_dsa_atomic_init_slot(&pool_ctx->work[i].is_free); }
	
	return pool_ctx;
	
	err_2:
		ml_dsa_free(pool_ctx);
	err_1:
		destructor_verify_works_struct(ctx, number_slots);
		return err_ml_dsa_alloc_verify_pool(-ret);
}

// Знищення пулу верифікацій - може приберу другйи параметр і очищення матриці залежатиме від значення вказівника t1 в структурі
void destructor_verify_pool(struct ml_dsa_verify_ctx *ctx)
{
	if(!ctx) { return; }
	
	if(ctx->work)
	{
		destructor_verify_works_struct(ctx->work, ctx->size_pool);
	}
	
	if(ctx->t1) { ml_dsa_free(ctx->t1); }
	
	ml_dsa_free(ctx);
}

// Отримання закритого ключ
int ml_dsa_get_sk(struct ml_dsa_keys *ctx, u8 *buffer, size_t size_buffer)
{
	if(!ctx || !ctx->workspace || !buffer) { return err_ml_dsa_get_sk(-ML_DSA_EINVAL); }
	
	if(!ctx->s1 || !ctx->s2 || !ctx->t0 || !ctx->workspace->vect_y || !ctx->workspace->vect_w || !ctx->workspace->temp_vector_buffer)
	{
		return err_ml_dsa_get_sk(-ML_DSA_EINVAL);
	}
	
	u8 eta;
	enum ml_dsa_level_l l;
	enum ml_dsa_level_k k;
	if(size_buffer == ML_DSA_44_SIZE_SK) { eta = 2; k = ML_DSA_44_K; l = ML_DSA_44_L; }
	else if(size_buffer == ML_DSA_65_SIZE_SK) { eta = 4; k = ML_DSA_65_K; l = ML_DSA_65_L; }
	else if(size_buffer == ML_DSA_87_SIZE_SK) { eta = 2; k = ML_DSA_87_K; l = ML_DSA_87_L; }
	else { return err_ml_dsa_get_sk(-ML_DSA_EINVAL); }
	
	if(k != ctx->k || l !=ctx->l) { return err_ml_dsa_get_sk(-ML_DSA_EINVAL); }
	
	size_t buffer_counter = 0;
	memcpy(buffer, ctx->rho, ML_DSA_32_BYTES);
	buffer_counter += ML_DSA_32_BYTES;
	
	memcpy(buffer + buffer_counter, ctx->K, ML_DSA_32_BYTES);
	buffer_counter += ML_DSA_32_BYTES;
	
	memcpy(buffer + buffer_counter, ctx->tr, ML_DSA_64_BYTES);
	buffer_counter += ML_DSA_64_BYTES;
	
	memcpy(ctx->workspace->vect_y, ctx->s1, ML_DSA_N * ctx->l * sizeof(s32));
	memcpy(ctx->workspace->vect_w, ctx->s2, ML_DSA_N * ctx->k * sizeof(s32));
	memcpy(ctx->workspace->temp_vector_buffer, ctx->t0, ML_DSA_N * ctx->k * sizeof(s32));
	
	for (size_t i = 0; i < l; i++) 
	{ 
		ml_dsa_intt(ctx->workspace->vect_y + i * ML_DSA_N);  
		for(size_t j = 0; j < ML_DSA_N; j++)
		{
			ctx->workspace->vect_y[i * ML_DSA_N + j] = ml_dsa_montgomery_reduce((s64)ctx->workspace->vect_y[i * ML_DSA_N + j]);
		}
	}
	
	for (size_t i = 0; i < k; i++) 
	{ 
		ml_dsa_intt(ctx->workspace->vect_w + i * ML_DSA_N);
		ml_dsa_intt(ctx->workspace->temp_vector_buffer + i * ML_DSA_N); 
		for(size_t j = 0; j < ML_DSA_N; j++)
		{
			ctx->workspace->vect_w[i * ML_DSA_N + j] = ml_dsa_montgomery_reduce((s64)ctx->workspace->vect_w[i * ML_DSA_N + j]);
			ctx->workspace->temp_vector_buffer[i * ML_DSA_N + j] = ml_dsa_montgomery_reduce((s64)ctx->workspace->temp_vector_buffer[i * ML_DSA_N + j]);
		}
	}
	
	if(eta == 2)
	{
		buffer_counter += pack_s1_and_s2_44_87(ctx->workspace->vect_y, buffer + buffer_counter, (u8)l);
		buffer_counter += pack_s1_and_s2_44_87(ctx->workspace->vect_w, buffer + buffer_counter, (u8)k);
	}else
	{
		buffer_counter += pack_s1_and_s2_65(ctx->workspace->vect_y, buffer + buffer_counter, (u8)l);
		buffer_counter += pack_s1_and_s2_65(ctx->workspace->vect_w, buffer + buffer_counter, (u8)k);
	}
	
	buffer_counter += pack_t0(ctx->workspace->temp_vector_buffer, buffer + buffer_counter, (u8)k);
	
	ml_dsa_memzero(ctx->workspace->temp_vector_buffer,(2 * k + l) * ML_DSA_N * sizeof(s32));
		
	if(buffer_counter != size_buffer) 
	{ 
		ml_dsa_memzero(buffer, size_buffer);
		return err_ml_dsa_get_sk(-ML_DSA_EINVAL); 
	}
	
	return (int)buffer_counter;
}

// Отримання відкритого ключа
int ml_dsa_get_pk(struct ml_dsa_keys *ctx, u8 *buffer, size_t size_buffer)
{
	if(!ctx || !buffer) { return err_ml_dsa_get_pk(-ML_DSA_EINVAL); }
	
	if(ctx->k == ML_DSA_44_K && ctx->l != ML_DSA_44_L) { return err_ml_dsa_get_pk(-ML_DSA_EINVAL); }
	if(ctx->k == ML_DSA_65_K && ctx->l != ML_DSA_65_L) { return err_ml_dsa_get_pk(-ML_DSA_EINVAL); }
	if(ctx->k == ML_DSA_87_K && ctx->l != ML_DSA_87_L) { return err_ml_dsa_get_pk(-ML_DSA_EINVAL); }
	
	if(ctx->k == ML_DSA_44_K && size_buffer != ML_DSA_44_SIZE_PK) { return err_ml_dsa_get_pk(-ML_DSA_EINVAL); }
	if(ctx->k == ML_DSA_65_K && size_buffer != ML_DSA_65_SIZE_PK) { return err_ml_dsa_get_pk(-ML_DSA_EINVAL); }
	if(ctx->k == ML_DSA_87_K && size_buffer != ML_DSA_87_SIZE_PK) { return err_ml_dsa_get_pk(-ML_DSA_EINVAL); }
	
	memcpy(buffer, ctx->pk, size_buffer);
	
	return (int)size_buffer;
}
	
	
//---------------------------------------------------------------------------------------------------------------------------------
// Створення контексту на основі закритого ключа-повідомлення згідно стандарту 
struct ml_dsa_keys *ml_dsa_create_keys_object_from_imported_sk(u8 *sk_buffer, size_t size_buffer, u8 level)
{
	if(!sk_buffer) { return err_ml_dsa_create_keys_object_from_imported_sk(-ML_DSA_EINVAL); }
	
	enum ml_dsa_level_l l;
	enum ml_dsa_level_k k;
	int ret;
	if(size_buffer == ML_DSA_32_BYTES && (level == 44 || level == 4))
	{
		k = ML_DSA_44_K; 
		l = ML_DSA_44_L;
	}else if(size_buffer == ML_DSA_32_BYTES && (level == 65 || level == 6))
	{
		k = ML_DSA_65_K; 
		l = ML_DSA_65_L;
	}else if(size_buffer == ML_DSA_32_BYTES && (level == 87 || level == 8))
	{
		k = ML_DSA_87_K; 
		l = ML_DSA_87_L;
	}else if(size_buffer == ML_DSA_44_SIZE_SK && (level == 44 || level == 4))
	{
		k = ML_DSA_44_K; 
		l = ML_DSA_44_L;
	}else if(size_buffer == ML_DSA_65_SIZE_SK && (level == 65 || level == 6))
	{
		k = ML_DSA_65_K; 
		l = ML_DSA_65_L;
	}else if(size_buffer == ML_DSA_87_SIZE_SK && (level == 87 || level == 8))
	{
		k = ML_DSA_87_K; 
		l = ML_DSA_87_L;
	}else { return err_ml_dsa_create_keys_object_from_imported_sk(-ML_DSA_EINVAL); }
	
	struct ml_dsa_keys *ctx;
	ctx = ml_dsa_alloc_struct_keys(k, l);
	if(!ctx) { return err_ml_dsa_create_keys_object_from_imported_sk(-ML_DSA_ENOMEM); }
	
	if(size_buffer == ML_DSA_32_BYTES)
	{
		ret = get_rho_K_s1_s2(ctx, NULL, 1, sk_buffer, ML_DSA_32_BYTES);
		if(ret != 0) { goto err; }
		ret = get_ctx_keys_from_seed(ctx);
		if(ret != 0) { goto err; }
	}else
	{
		ret = get_ctx_keys_from_sk(ctx, sk_buffer);
		if(ret != 0) { goto err; }
	}
	
	ml_dsa_memzero(ctx->workspace->temp_vector_buffer, (ctx->l + 2 * ctx->k + 1) * ML_DSA_N * sizeof(s32));
	
	return ctx;
	
	err:
		ml_dsa_destroy_struct_keys(ctx);
		return err_ml_dsa_create_keys_object_from_imported_sk(-ret);
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

static inline int get_ctx_keys_from_sk(struct ml_dsa_keys *ctx, u8 *sk_buffer)
{
	if(!ctx || !ctx->workspace || !sk_buffer) { return ML_DSA_EINVAL; }
	if(!ctx->s1 || !ctx->s2 || !ctx->t0 || !ctx->workspace) { return ML_DSA_EINVAL; }
	if(!ctx->workspace->temp_vector_buffer || !ctx->workspace->vect_w || !ctx->workspace->vect_y 
		|| !ctx->workspace->poly_c || !ctx->workspace->scratch_buffer || !ctx->workspace->shake) { return ML_DSA_EINVAL; }
	
	if(ctx->k == ML_DSA_44_K && ctx->l != ML_DSA_44_L) { return ML_DSA_EINVAL; }
	if(ctx->k == ML_DSA_65_K && ctx->l != ML_DSA_65_L) { return ML_DSA_EINVAL; }
	if(ctx->k == ML_DSA_87_K && ctx->l != ML_DSA_87_L) { return ML_DSA_EINVAL; }
	if(ctx->k != ML_DSA_44_K && ctx->k != ML_DSA_65_K && ctx->k != ML_DSA_87_K) { return ML_DSA_EINVAL; }
	
	size_t buffer_counter = 0;
	
	memcpy(ctx->rho, sk_buffer, ML_DSA_32_BYTES);
	buffer_counter += ML_DSA_32_BYTES;
	
	memcpy(ctx->K, sk_buffer + buffer_counter, ML_DSA_32_BYTES);
	buffer_counter += ML_DSA_32_BYTES;
	
	memcpy(ctx->tr, sk_buffer + buffer_counter, ML_DSA_64_BYTES);
	buffer_counter += ML_DSA_64_BYTES;
	
	// The magic value 128 and 96 it is number bytes for one coefficient in bytes stream
	size_t temp_val = 0;
	if(ctx->k == ML_DSA_44_K || ctx->k == ML_DSA_87_K)
	{
		temp_val = unpack_s1_and_s2_44_87(ctx->s1, sk_buffer + buffer_counter, ctx->l);
		if(temp_val != 96 * ctx->l) { return ML_DSA_EINVAL; }
		else { buffer_counter += temp_val; }
		
		temp_val = unpack_s1_and_s2_44_87(ctx->s2, sk_buffer + buffer_counter, ctx->k);
		if(temp_val != 96 * ctx->k) { return ML_DSA_EINVAL; }
		else { buffer_counter += temp_val; }
	}else
	{
		temp_val = unpack_s1_and_s2_65(ctx->s1, sk_buffer + buffer_counter, ctx->l);
		if(temp_val != 128 * ctx->l) { return ML_DSA_EINVAL; }
		else { buffer_counter += temp_val; }
		
		temp_val = unpack_s1_and_s2_65(ctx->s2, sk_buffer + buffer_counter, ctx->k);
		if(temp_val != 128 * ctx->k) { return ML_DSA_EINVAL; }
		else { buffer_counter += temp_val; }
	}
	
	buffer_counter += unpack_t0(ctx->t0, sk_buffer + buffer_counter, ctx->k);
	
	if(ctx->k == ML_DSA_44_K) { if(buffer_counter != ML_DSA_44_SIZE_SK) { return ML_DSA_EINVAL; } }
	else if(ctx->k == ML_DSA_65_K) { if(buffer_counter != ML_DSA_65_SIZE_SK) { return ML_DSA_EINVAL; } }
	else if(ctx->k == ML_DSA_87_K) { if(buffer_counter != ML_DSA_87_SIZE_SK) { return ML_DSA_EINVAL; } }
	else { return ML_DSA_EINVAL; }
	
	// Копіювання перших поліномів в vect_y та останнтього(якщо не 44 рівень) в поліном poly_c - з метою збереження розпакованого t0
	for (size_t i = 0; i < ctx->k; i++)  { ml_dsa_ntt(ctx->t0 + i * ML_DSA_N); }
	memcpy(ctx->workspace->vect_y, ctx->t0, ctx->l * ML_DSA_N * sizeof(s32));
	if(ctx->k != ML_DSA_44_K)
	{
		memcpy(ctx->workspace->poly_c, ctx->t0 + ctx->l * ML_DSA_N, ML_DSA_N * sizeof(s32));
	}
	
	if(get_ctx_keys_from_seed(ctx)) { return ML_DSA_EINVAL; }
	
	bool val = 1;
	size_t i = 0;
	while(i < ML_DSA_N * ctx->l)
	{
		val &= ctx->workspace->vect_y[i] == ctx->t0[i];
		i++;
	}
	
	if(ctx->k != ML_DSA_44_K)
	{
		while(i < ML_DSA_N * ctx->k)
		{
			val &= ctx->workspace->poly_c[i - ctx->l * ML_DSA_N] == ctx->t0[i];
			i++;
		}
	}
	
	u8 *temp_tr = sk_buffer + ML_DSA_64_BYTES;
	for(size_t i = 0; i < ML_DSA_64_BYTES; i++) { val &= ctx->tr[i] == temp_tr[i]; }
	
	if(!val) { return ML_DSA_EINVAL; }
	
	return 0;
}

static inline int get_ctx_keys_from_seed(struct ml_dsa_keys *ctx)
{
	if(!ctx) { return ML_DSA_EINVAL; }
	int ret = 0;
	
	// Here should be copy s1 before NTT(it is can do in vect_y - temp)
	for (size_t i = 0; i < ctx->l; i++) { ml_dsa_ntt(ctx->s1 + i * ML_DSA_N); }
	
	mult_matrix(ctx->rho, ctx->s1, ctx->workspace->temp_vector_buffer, ctx->workspace->shake, ctx->workspace->scratch_buffer, ctx->k, ctx->l);
	
	for (size_t i = 0; i < ctx->k; i++) { ml_dsa_intt(ctx->workspace->temp_vector_buffer + i * ML_DSA_N); }
	
	for(size_t i = 0; i < ctx->k; i++)
	{
		for(size_t j = 0; j < ML_DSA_N; j++)
		{
			ctx->workspace->temp_vector_buffer[ML_DSA_N * i + j] += ctx->s2[ML_DSA_N * i + j];
		}
	}
	
	ret = create_t0_t1(ctx, ctx->workspace->temp_vector_buffer);
	if(ret != 0) { return ret; }
	
	ret = ml_dsa_create_public_data(ctx);
	if(ret != 0) { return ret; }

	// Here - create functions for format pk messege and sk messege
	
	for (size_t i = 0; i < ctx->k; i++) 
	{ 
		ml_dsa_ntt(ctx->s2 + i * ML_DSA_N);
		ml_dsa_ntt(ctx->t0 + i * ML_DSA_N); 
	}
	
	return 0;
}

//---------------------------------------------------------------------------------------------------------------------------------

static inline size_t unpack_s1_and_s2_44_87(s32 *s, const u8 *buffer, u8 lng_loops)
{
	size_t one_coef_ctr = 0;
	size_t buffer_ctr = 0;
	bool invalid = 1;
	while(one_coef_ctr < ML_DSA_N * lng_loops)
	{
		s[one_coef_ctr] = 2 - ((s32)(buffer[buffer_ctr] & 0x7));
		invalid &= s[one_coef_ctr] <= 2 && s[one_coef_ctr] >= -2;
		one_coef_ctr++;
		
		s[one_coef_ctr] = 2 - ((s32)((buffer[buffer_ctr] >> 3) & 0x7));
		invalid &= s[one_coef_ctr] <= 2 && s[one_coef_ctr] >= -2;
		one_coef_ctr++;
		
		s[one_coef_ctr] = 2 - ((s32)((buffer[buffer_ctr] >> 6) & 0x3) | ((buffer[buffer_ctr+1] & 0x1) << 2));
		invalid &= s[one_coef_ctr] <= 2 && s[one_coef_ctr] >= -2;
		one_coef_ctr++;
		buffer_ctr++;
		
		s[one_coef_ctr] = 2 - ((s32)((buffer[buffer_ctr] >> 1) & 0x7));
		invalid &= s[one_coef_ctr] <= 2 && s[one_coef_ctr] >= -2;
		one_coef_ctr++;
		
		s[one_coef_ctr] = 2 - ((s32)((buffer[buffer_ctr] >> 4) & 0x7));
		invalid &= s[one_coef_ctr] <= 2 && s[one_coef_ctr] >= -2;
		one_coef_ctr++;
		
		s[one_coef_ctr] = 2 - ((s32)((buffer[buffer_ctr] >> 7) & 0x1) | ((buffer[buffer_ctr+1] & 0x3) << 1));
		invalid &= s[one_coef_ctr] <= 2 && s[one_coef_ctr] >= -2;
		one_coef_ctr++;
		buffer_ctr++;
		
		s[one_coef_ctr] = 2 - ((s32)((buffer[buffer_ctr] >> 2) & 0x7));
		invalid &= s[one_coef_ctr] <= 2 && s[one_coef_ctr] >= -2;
		one_coef_ctr++;
		
		s[one_coef_ctr] = 2 - ((s32)((buffer[buffer_ctr] >> 5) & 0x7));
		invalid &= s[one_coef_ctr] <= 2 && s[one_coef_ctr] >= -2;
		one_coef_ctr++;
		buffer_ctr++;
	}
	
	return buffer_ctr * (size_t)invalid;
}
		
static inline size_t unpack_s1_and_s2_65(s32 *s, const u8 *buffer, u8 lng_loops)
{
	size_t one_coef_ctr = 0;
	size_t buffer_ctr = 0;
	bool invalid = 1;
	while(one_coef_ctr < ML_DSA_N * lng_loops)
	{
		s[one_coef_ctr] = 4 - (s32)(buffer[buffer_ctr] & 0xF);
		invalid &= s[one_coef_ctr] <= 4 && s[one_coef_ctr] >= -4;
		one_coef_ctr++;
		
		s[one_coef_ctr] = 4 - (s32)(buffer[buffer_ctr] >> 4);
		invalid &= s[one_coef_ctr] <= 4 && s[one_coef_ctr] >= -4;
		one_coef_ctr++;
		buffer_ctr++;
	}
	
	return buffer_ctr * (size_t)invalid;
}		
	
static inline size_t pack_s1_and_s2_44_87(const s32 *s, u8 *buffer, u8 lng_loops)
{
	size_t one_coef_ctr = 0;
	size_t buffer_ctr = 0;
	s32 temp;
	while(one_coef_ctr < ML_DSA_N * lng_loops)
	{
		temp = 2 - s[one_coef_ctr];
		buffer[buffer_ctr] = (u8)temp;
		one_coef_ctr++;
		
		temp = 2 - s[one_coef_ctr];
		buffer[buffer_ctr] |= (u8)(temp << 3);
		one_coef_ctr++;
		
		temp = 2 - s[one_coef_ctr];
		buffer[buffer_ctr] |= (u8)(temp << 6);
		one_coef_ctr++;
		
		buffer_ctr++;
		buffer[buffer_ctr] = (u8)(temp >> 2);
		
		temp = 2 - s[one_coef_ctr];
		buffer[buffer_ctr] |= (u8)(temp << 1);
		one_coef_ctr++;
		
		temp = 2 - s[one_coef_ctr];
		buffer[buffer_ctr] |= (u8)(temp << 4);
		one_coef_ctr++;
		
		temp = 2 - s[one_coef_ctr];
		buffer[buffer_ctr] |= (u8)(temp << 7);
		one_coef_ctr++;
		
		buffer_ctr++;
		buffer[buffer_ctr] = (u8)(temp >> 1);
		
		temp = 2 - s[one_coef_ctr];
		buffer[buffer_ctr] |= (u8)(temp << 2);
		one_coef_ctr++;
		
		temp = 2 - s[one_coef_ctr];
		buffer[buffer_ctr] |= (u8)(temp << 5);
		one_coef_ctr++;
		buffer_ctr++;
	}
	
	return buffer_ctr;
}

static inline size_t pack_s1_and_s2_65(const s32 *s, u8 *buffer, u8 lng_loops)
{
	size_t one_coef_ctr = 0;
	size_t buffer_ctr = 0;
	s32 temp;
	while(one_coef_ctr < ML_DSA_N * lng_loops)
	{
		temp = 4 - s[one_coef_ctr];
		buffer[buffer_ctr] = (u8)temp;
		one_coef_ctr++;
		
		temp = 4 - s[one_coef_ctr];
		buffer[buffer_ctr] |= (u8)(temp << 4);
		one_coef_ctr++;
		buffer_ctr++;
	}
	
	return buffer_ctr;
}

//==========================================================================================================================================================

static inline size_t pack_t0(const s32 *t0, u8 *buffer, u8 k)
{
    size_t buffer_ctr = 0;
    u32 accumulator = 0;
    u32 bits = 0;

    for(size_t i = 0; i < ML_DSA_N * k; i++)
    {
        u32 value = (u32)(4096 - t0[i]);

        accumulator |= value << bits;
        bits += 13;

        while(bits >= 8)
        {
            buffer[buffer_ctr++] = (u8)accumulator;
            accumulator >>= 8;
            bits -= 8;
        }
    }
    
    return buffer_ctr;
}

static inline size_t unpack_t0(s32 *t0, const u8 *buffer, u8 k)
{
    size_t buffer_ctr = 0;
    u32 accumulator = 0;
    u32 bits = 0;

    for(size_t i = 0; i < ML_DSA_N * (size_t)k; i++)
    {
        while(bits < 13)
        {
            accumulator |= (u32)buffer[buffer_ctr++] << bits;
            bits += 8;
        }

        t0[i] = 4096 - (s32)(accumulator & 0x1FFFU);

        accumulator >>= 13;
        bits -= 13;
    }

    return buffer_ctr;
}
		
	
	
	
	










