
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
	
		
	
	
	
	
