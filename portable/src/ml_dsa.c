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
	
	ret = mult_matrix(ctx, ctx->s1, ctx->workspace->temp_vector_buffer);
	if(ret != 0) { goto err_2; }
	
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
	else if(!msg || !sig) { return ML_DSA_EINVAL; }
	else if(context_len > ML_DSA_SIZE_MAX_CONTEXT) { return ML_DSA_EINVAL; }
	else if(context_len != 0 && !context) { return ML_DSA_EINVAL; }
	else if(ctx->k == ML_DSA_44_K && sig_len < ML_DSA_44_SIZE_SIG) { return ML_DSA_EINVAL; }
	else if(ctx->k == ML_DSA_65_K && sig_len < ML_DSA_65_SIZE_SIG) { return ML_DSA_EINVAL; }
	else if(ctx->k == ML_DSA_87_K && sig_len < ML_DSA_87_SIZE_SIG) { return ML_DSA_EINVAL; }
	else if(deterministic > 1) { return ML_DSA_EINVAL; }
	else if(prehash > 1) { return ML_DSA_EINVAL; }
	
	// Якщо режим тільки буфера, то тільки перевіряю відповідність розміра, якщо не збіг - повертаю помилку
	#if ML_DSA_MEM_MODE == ML_DSA_MEM_BUFFER 
	
	if(prehash == 0 && msg_len > ML_DSA_BUFFER_SIZE) { return ML_DSA_ENOMEM; }
	else if(prehash == 1 && (2 + context_len + ML_DSA_PREHASH_OID_SIZE + 2 * ML_DSA_64_BYTES) > ML_DSA_MAX_SIZE_FOR_FORMAT_M) { return ML_DSA_ENOMEM; }
	
	// Якщо гібридний режим, тоді .....
	#elif ML_DSA_MEM_MODE == ML_DSA_MEM_HYBRID
	
	// Створюю резервний тимчасовий вказівник для буфера та значення для кількості аллокаованної пам'яті
	u8 *buffer_ptr = NULL;
	size_t size_mem = 0;
	
	// Перевіряю, чи вистачить місце в буферу, якщо умова не виконана, то вистачить і аллокація не відбуваєстяь і тимчасові змінні не потрібні
	if((prehash == 0 && msg_len > ML_DSA_BUFFER_SIZE) ||
		(prehash == 1 && (2 + context_len + ML_DSA_PREHASH_OID_SIZE + 2 * ML_DSA_64_BYTES) > ML_DSA_MAX_SIZE_FOR_FORMAT_M))
		{
			// Якщо умова виконана, то щоб не було витоку - буфер зберігаю в тимчасовому вказівнику, аллокую пам'ять у відповідне
			// Поле структуру, бо працюю тільки з ним
			buffer_ptr = ctx->workspace->messege;
			if(prehash == 1) { size_mem = ML_DSA_64_BYTES + 2 + context_len + ML_DSA_PREHASH_OID_SIZE + ML_DSA_64_BYTES; }
			else { size_mem = ML_DSA_64_BYTES + 2 + context_len + msg_len; }
			ctx->workspace->messege = ml_dsa_alloc(size_mem);
			if(!ctx->workspace->messege) { return ML_DSA_ENOMEM; }
		}
	
	// Якщо режим тільки аллокації - то просто вираховую стільки треба, запитую в системи - якщо ні, то помилка, не вистачає памя'ті
	#elif ML_DSA_MEM_MODE == ML_DSA_MEM_ALLOC
	
	size_t size_mem = 0;
	if(prehash == 1) { size_mem = ML_DSA_64_BYTES + 2 + context_len + ML_DSA_PREHASH_OID_SIZE + ML_DSA_64_BYTES; }
	else { size_mem = ML_DSA_64_BYTES + 2 + context_len + msg_len; }
	
	ctx->workspace->messege = ml_dsa_alloc(size_mem);
	if(!ctx->workspace->messege) { return ML_DSA_ENOMEM; }
		
	
	#endif // Завершена гра з режимами - далі працює алгоритм
	
	
	
	
	
	
	// Create mu
	get_mu(ctx, msg, msg_len, context, context_len, prehash);
	int ret = get_rho_double_prime(ctx, deterministic, entropy);
	if(ret != 0) { return ret; }
	
	
	u32 kappa = 0;
	u16 counter = 0;
	while(counter <= 814)
	{
		// Get vector y
		ml_dsa_expand_mask(ctx, kappa);
		memcpy(ctx->workspace->temp_vector_buffer, ctx->workspace->vect_y, sizeof(s32) * l * ML_DSA_N);
		
		// Get vectors w and w1
		for(size_t i = 0; i < ctx->l; i++) { ml_dsa_ntt(ctx->workspace->temp_vector_buffer + i * ML_DSA_N); }
		mult_matrix(ctx, ctx->workspace->temp_vector_buffer, ctx->workspace->vect_w);
		for(size_t i = 0; i < ctx->k; i++) { ml_dsa_intt(ctx->workspace->vect_w + i * ML_DSA_N); }
		for(size_t i = 0; i < ctx->k; i++) { ml_dsa_canonicalize(ctx->workspace->vect_w + i * ML_DSA_N); }
		ml_dsa_get_w1(ctx);
		
		// Get polynomial c
		u8 c[ML_DSA_64_BYTES];
		ml_dsa_get_c(ctx, c);
		ml_dsa_get_poly_c(ctx, c);
		
	
	
	
	
	
	
	
	
	
	
	
	
	
	// Якщо тільки буфер, то просто занулюю його і все
	#if ML_DSA_MEM_MODE == ML_DSA_MEM_BUFFER
	
	ml_dsa_memzero(ctx->workspace->messege, ML_DSA_MAX_SIZE_FOR_FORMAT_M);
	
	// При гібридному преевіряю .....
	#elif ML_DSA_MEM_MODE == ML_DSA_MEM_HYBRID
	
	// Якщо вказівник тимчасового зберігання пустий - то працювало тільки через буфер, і аллокації не було, то тільки занулення
	if(!buffer_ptr) { ml_dsa_memzero(ctx->workspace->messege, ML_DSA_MAX_SIZE_FOR_FORMAT_M); }
	else // Якщо ні - то була аллокація, занулюю, очищаю, та перевизначю вказівники
	{
		ml_dsa_memzero(ctx->workspace->messege, size_mem);
		ml_dsa_free(ctx->workspace->messege);
		ctx->workspace->messege = buffer_ptr;
		buffer_ptr = NULL;
	}
	
	// Якщо режим тільки аллокації, то просто занулюю і очищаю
	#elif ML_DSA_MEM_MODE == ML_DSA_MEM_ALLOC
	
	ml_dsa_memzero(ctx->workspace->messege, size_mem);
	ml_dsa_free(ctx->workspace->messege);
	ctx->workspace->messege = NULL;
	
	#endif
	
	return 0;
}
