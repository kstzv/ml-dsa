
static const uint8_t oid_shake256[11] = {
    0x06, 0x09, 0x60, 0x86, 0x48, 0x01,
    0x65, 0x03, 0x04, 0x02, 0x0C
};

void get_mu(struct ml_dsa_keys *ctx, const uint8_t *msg, size_t msg_len, const uint8_t *context, size_t context_len, u8 prehash)
{
	u8 *temp_ptr = ctx->workspace->messege;
	size_t counter_ptr = 0;
	memcpy(temp_ptr, ctx->tr, ML_DSA_64_BYTES);
	counter_ptr += ML_DSA_64_BYTES;
	temp_ptr += ML_DSA_64_BYTES;
	
	*temp_ptr = prehash;
	temp_ptr++;
	counter_ptr++;
	
	*temp_ptr = (u8)context_len;
	temp_ptr++;
	counter_ptr++;
	
	memcpy(temp_ptr, context, context_len);
	counter_ptr += context_len;
	temp_ptr += context_len;
	
	if(prehash == 0)
	{
		memcpy(temp_ptr, msg, msg_len);
		counter_ptr += msg_len;
		temp_ptr += msg_len;
	}else
	{
		shake_ctx_zero(ctx->workspace->shake);
		shake_ctx_init(ctx->workspace->shake, ctx->workspace->mu, ML_DSA_64_BYTES, msg, msg_len);
		shake256(ctx->workspace->shake);
		memcpy(temp_ptr, oid_shake256, ML_DSA_PREHASH_OID_SIZE);
		counter_ptr += ML_DSA_PREHASH_OID_SIZE;
		temp_ptr += ML_DSA_PREHASH_OID_SIZE;
		memcpy(temp_ptr, ctx->workspace->mu, ML_DSA_64_BYTES);
		counter_ptr += ML_DSA_64_BYTES;
		temp_ptr += ML_DSA_64_BYTES;
	}
	
	shake_ctx_zero(ctx->workspace->shake);
	shake_ctx_init(ctx->workspace->shake, ctx->workspace->mu, ML_DSA_64_BYTES, ctx->workspace->messege, counter_ptr);
	shake256(ctx->workspace->shake);
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

void ml_dsa_expand_mask(struct ml_dsa_keys *ctx, s32 *y, u32 kappa)
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
            y[(size_t)r * ML_DSA_N + i] = (s32)gamma1 - (s32)value;
		}
	}
}
	
		
	
	
	
	
