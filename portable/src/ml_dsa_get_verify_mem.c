
#include "ml_dsa_core_header.h"

struct ml_dsa_verify_workspace *ml_dsa_alloc_verify_works_struct(enum ml_dsa_level_k k, enum ml_dsa_level_l l, size_t number_slots)
{
	if(number_slots == 0 || number_slots > MAX_NUMBER_SLOTS) { return NULL; }
	
	size_t size_scratch_buffer;
	if(k == ML_DSA_44_K) { size_scratch_buffer = ML_DSA_SIZE_SCRATCH_BUFFER; }
	else if(k == ML_DSA_65_K) { size_scratch_buffer = ML_DSA_SIZE_SCRATCH_BUFFER + 512; }
	else if(k == ML_DSA_87_K) { size_scratch_buffer = ML_DSA_SIZE_SCRATCH_BUFFER * 2; }
	else { return NULL; }
	
	if(l != ML_DSA_44_L && l != ML_DSA_65_L && l != ML_DSA_87_L) { return NULL; }
	
	if(k == ML_DSA_44_K && l != ML_DSA_44_L) { return NULL; }
	if(k == ML_DSA_65_K && l != ML_DSA_65_L) { return NULL; }
	if(k == ML_DSA_87_K && l != ML_DSA_87_L) { return NULL; }
	
	struct ml_dsa_verify_workspace *ctx;
	ctx = ml_dsa_alloc(sizeof(struct ml_dsa_verify_workspace) * number_slots);
	if(!ctx) { return NULL; }
	ml_dsa_memzero(ctx, sizeof(struct ml_dsa_verify_workspace) * number_slots);
	
	ctx->poly_c = ml_dsa_alloc(((ML_DSA_N + (ML_DSA_N * k * 2)) * sizeof(s32)) * number_slots);
	if(!ctx->poly_c) { goto err_1; }
	ml_dsa_memzero(ctx->poly_c, ((ML_DSA_N + (ML_DSA_N * k * 2)) * sizeof(s32)) * number_slots);
	ctx->temp_vector_1 = ctx->poly_c + ML_DSA_N;
	ctx->temp_vector_2 = ctx->temp_vector_1 + ML_DSA_N * k;
	
	ctx->scratch_buffer = ml_dsa_alloc(size_scratch_buffer * number_slots);
	if(!ctx->scratch_buffer) { goto err_2; }
	ml_dsa_memzero(ctx->scratch_buffer, size_scratch_buffer * number_slots);
	
	ctx->shake = ml_dsa_alloc(sizeof(struct shake_ctx) * number_slots);
	if(!ctx->shake) { goto err_3; }
	ml_dsa_memzero(ctx->shake, sizeof(struct shake_ctx) * number_slots);
	
	ctx->shake->state = ml_dsa_alloc(ML_DSA_SIZE_MEM_KECCAK * sizeof(u64) * number_slots);
	if(!ctx->shake->state) { goto err_4; }
	ml_dsa_memzero(ctx->shake->state, ML_DSA_SIZE_MEM_KECCAK * sizeof(u64) * number_slots);
	int ret = shake_ctx_get_mem(ctx->shake, ctx->shake->state, ML_DSA_SIZE_MEM_KECCAK);
	if(ret != 0) { goto err_5; }
	
	#if ML_DSA_FULL_MATRIX_BUFFER
	ctx->matrix_buffer = ml_dsa_alloc((ML_DSA_N * k * k) * sizeof(s32));
	if(!ctx->matrix_buffer) { goto err_5; }
	ml_dsa_memzero(ctx->matrix_buffer, (ML_DSA_N * k * k) * sizeof(s32));
	#endif
	
	ctx->k = k;
	ctx->l = l;
	ctx->next = NULL;
	// Додати ще атомарну ініціалізації ml_dsa_atomic_t is_free;
	
	if(number_slots == 1) { return ctx; }
	
	s32 s32_conter = (ML_DSA_N + (ML_DSA_N * k * 2));
	size_t scratch_buffer_conter = size_scratch_buffer;
	u8 shake_counter = 1;
	u64 shake_mem_conter = ML_DSA_SIZE_MEM_KECCAK;
	
	for(size_t i = 1; i < number_slots; i++)
	{
		struct ml_dsa_verify_workspace *inside_ctx = &ctx[i];
		
		inside_ctx->poly_c = ctx->poly_c + s32_conter;
		s32_conter += ML_DSA_N;
		
		inside_ctx->temp_vector_1 = ctx->poly_c + s32_conter;
		s32_conter += ML_DSA_N * k;
		
		inside_ctx->temp_vector_2 = ctx->poly_c + s32_conter;
		s32_conter += ML_DSA_N * k;
		
		inside_ctx->scratch_buffer = ctx->scratch_buffer + scratch_buffer_conter;
		scratch_buffer_conter += size_scratch_buffer;
		
		inside_ctx->shake = ctx->shake + shake_conter;
		shake_conter++;
		
		(void)shake_ctx_get_mem(inside_ctx->shake, ctx->shake->state + shake_mem_conter, ML_DSA_SIZE_MEM_KECCAK);
		shake_mem_conter += ML_DSA_SIZE_MEM_KECCAK;
		
		#if ML_DSA_FULL_MATRIX_BUFFER
		inside_ctx->matrix_buffer = ctx->matrix_buffer;
		#endif
		
		inside_ctx->k = k;
		inside_ctx->l = l;
		// Додати ще атомарну ініціалізації ml_dsa_atomic_t is_free;
	}
	
	for(size_t i = 0; i + 1 < number_slots; i++)
	{
		ctx[i].next = &ctx[i + 1];
	}
	ctx[number_slots - 1].next = NULL;
	
	return ctx;
	
	err_5:
		ml_dsa_free(ctx->shake->state);
	err_4:
		ml_dsa_free(ctx->shake);
	err_3:
		ml_dsa_free(ctx->scratch_buffer);
	err_2:
		ml_dsa_free(ctx->poly_c);
	err_1:
		ml_dsa_free(ctx);
	return NULL;
}


struct ml_dsa_verify_ctx *ml_dsa_alloc_verify_pool(enum ml_dsa_level_k k, enum ml_dsa_level_l l, size_t number_slots, u8 *pk)
{
	if(number_slots == 0 || number_slots > MAX_NUMBER_SLOTS) { return NULL; }
	struct ml_dsa_verify_workspace *ctx;
	if(k == ML_DSA_44_K && l != ML_DSA_44_L) { return NULL; }
	else if(k == ML_DSA_65_K && l != ML_DSA_65_L) { return NULL; }
	else if(k == ML_DSA_87_K && l != ML_DSA_87_L) { return NULL; }
	else if(k != ML_DSA_44_K && k != ML_DSA_65_K && k != ML_DSA_87_K) { return NULL; } 
	else{ ctx = ml_dsa_alloc_verify_works_struct(k, l, number_slots); if(!ctx) { return NULL; } }
	
	struct ml_dsa_verify_ctx *pool_ctx;
	pool_ctx = ml_dsa_alloc(sizeof(struct ml_dsa_verify_ctx));
	if(!pool_ctx) { destructor_verify_works_struct(ctx, number_slots); return NULL; }
	
	pool_ctx->work = ctx;
	pool_ctx->size_pool = number_slots;
	
	if(pk)
	{
		// TODO - тут плануєтсья логіка розпакування відкритого ключа та заповнення всіх елементів пулу значенням готових t0 та ін.)
	}
	
	return pool_ctx;
}

void destructor_verify_works_struct(struct ml_dsa_verify_workspace *ctx, size_t size_pool)
{
	if(!ctx) { return; }
	
	#if ML_DSA_FULL_MATRIX_BUFFER
	if(ctx->matrix_buffer)
	{
		ml_dsa_memzero(ctx->matrix_buffer, (ML_DSA_N * ctx->k * ctx->k) * sizeof(s32));
		ml_dsa_free(ctx->matrix_buffer);
	}
	#endif
	
	if(ctx->poly_c)
	{
		ml_dsa_memzero(ctx->poly_c, ((ML_DSA_N + (ML_DSA_N * ctx->k * 2)) * sizeof(s32)) * size_pool);
		ml_dsa_free(ctx->poly_c);
	}
	
	if(ctx->scratch_buffer)
	{
		size_t size_scratch_buffer = ML_DSA_SIZE_SCRATCH_BUFFER;
		if(ctx->k == ML_DSA_65_K) { size_scratch_buffer = ML_DSA_SIZE_SCRATCH_BUFFER + 512; }
		else if(ctx->k == ML_DSA_87_K) { size_scratch_buffer = ML_DSA_SIZE_SCRATCH_BUFFER * 2; }
		ml_dsa_memzero(ctx->scratch_buffer, size_scratch_buffer * size_pool);
		ml_dsa_free(ctx->scratch_buffer);
	}
	
	if(ctx->shake->state)
	{
		ml_dsa_memzero(ctx->shake->state, (ML_DSA_SIZE_MEM_KECCAK * sizeof(u64)) * size_pool);
		ml_dsa_free(ctx->shake->state);
	}
	
	if(ctx->shake)
	{
		ml_dsa_memzero(ctx->shake, sizeof(struct shake_ctx) * size_pool);
		ml_dsa_free(ctx->shake);
	}
	
	ml_dsa_memzero(ctx, sizeof(struct ml_dsa_verify_workspace) * size_pool);
	ml_dsa_free(ctx);
}

void destructor_verify_pool(struct ml_dsa_verify_ctx *ctx)
{
	if(!ctx || !ctx->work) { return; }
	
	destructor_verify_works_struct(ctx->work, ctx->size_pool);
	
	ml_dsa_free(ctx);
}
	







