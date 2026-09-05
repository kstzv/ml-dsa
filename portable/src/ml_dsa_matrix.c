// ML-DSA create and use matrix A
// This implementation is incomplete and under active development
#include "ml_dsa_core_header.h"

// Multiply matrix A by a vector in the NTT domain
// The result buffer must be zero-initialized by the caller

int mult_matrix(struct ml_dsa_keys *ctx, s32 *vect, s32 *result);

#if defined(ML_DSA_FULL_MATRIX_BUFFER)

// Expand the full matrix A from rho and store it in matrix_buffer
// Matrix entries are generated as RejNTTPoly(rho || column || row)
int get_full_matrix(struct ml_dsa_keys *ctx)
{
	if(!ctx) { return ML_DSA_EINVAL; }
	
	u8 curr_seed[ML_DSA_32_BYTES + 2];
	memcpy(curr_seed, ctx->rho, ML_DSA_32_BYTES);
	curr_seed[32] = 0;
	curr_seed[33] = 0;
	
	for(size_t i = 0; i < ctx->k; i++)
	{
		for(size_t j = 0; j < ctx->l; j++)
		{
			curr_seed[32] = j;
			curr_seed[33] = i;
			ml_dsa_shake128(ctx->scratch_buffer, 894, curr_seed, ML_DSA_32_BYTES + 2);
			
			s32 *poly = ctx->matrix_buffer + ((i * ctx->l + j) * ML_DSA_N);
			
			// Rejection-sample 23-bit values below q
			size_t counter_in = 0;
            size_t counter_out = 0;
			while(counter_in < ML_DSA_N && counter_out + 2 < 894)
			{
				u8 a = ctx->scratch_buffer[counter_out++];
                u8 b = ctx->scratch_buffer[counter_out++];
                u8 c = ctx->scratch_buffer[counter_out++];
				s32 z = a | ((s32)b << 8) | ((s32)(c & 0x7F) << 16);
				
				if(z < ML_DSA_Q)
				{
					poly[counter_in++] = z;
				}
			}
			if (counter_in != ML_DSA_N) { return ML_DSA_EAGAIN; }
		}
	}
	return 0;
}

// Multiply the pre-expanded matrix A by vect
// Both operands are interpreted in the NTT domain
int mult_matrix(struct ml_dsa_keys *ctx, s32 *vect, s32 *result)
{
	if(!ctx || !vect || !result) { return ML_DSA_EINVAL; }
	
	for (size_t i = 0; i < ctx->k; i++)
	{
		for (size_t j = 0; j < ctx->l; j++)
		{
			s32 *poly_m = ctx->matrix_buffer + ((i * ctx->l + j) * ML_DSA_N);
			s32 *poly_v = vect + (j * ML_DSA_N);
			s32 *poly_r = result + (i * ML_DSA_N);

			// Pointwise multiplication in the NTT domain
			for (int c = 0; c < ML_DSA_N; c++)
			{
				poly_r[c] += ml_dsa_montgomery_reduce((s64)poly_m[c] * poly_v[c]);
			}
		}
	}
	return 0;
}

#else

// Stream matrix A from rho and multiply entries immediately
// No full matrix is materialized in memory
int mult_matrix(struct ml_dsa_keys *ctx, s32 *vect, s32 *result)
{
	if(!ctx || !vect || !result) { return ML_DSA_EINVAL; }
	
	u8 curr_seed[ML_DSA_32_BYTES + 2];
	memcpy(curr_seed, ctx->rho, ML_DSA_32_BYTES);
	
	for(size_t i = 0; i < ctx->k; i++)
	{
		s32 *poly_r = result + i * ML_DSA_N;
		for(size_t j = 0; j < ctx->l; j++)
		{
			const s32 *poly_v = vect + j * ML_DSA_N;
			curr_seed[32] = (u8)j;
			curr_seed[33] = (u8)i;
			ml_dsa_shake128(ctx->scratch_buffer, 894, curr_seed, ML_DSA_32_BYTES + 2);
			
			// Rejection-sample matrix coefficients and consume them immediately in the pointwise product
			size_t counter_in = 0;
            size_t counter_out = 0;
            while(counter_in < ML_DSA_N && counter_out + 2 < 894)
            {
				u8 a = ctx->scratch_buffer[counter_out++];
                u8 b = ctx->scratch_buffer[counter_out++];
                u8 c = ctx->scratch_buffer[counter_out++];
				s32 z = a | ((s32)b << 8) | ((s32)(c & 0x7F) << 16);
				
				if(z < ML_DSA_Q)
				{
					poly_r[counter_in] += ml_dsa_montgomery_reduce((s64)z * poly_v[counter_in]);
					counter_in++;
				}
			}
			if (counter_in != ML_DSA_N) { return ML_DSA_EAGAIN; }
		}
	}
	return 0;
}			

#endif
