// ML-DSA create and use matrix A
// This implementation is incomplete and under active development
#include "ml_dsa_core_header.h"

// Multiply matrix A by a vector in the NTT domain
// The result buffer must be zero-initialized by the caller

void mult_matrix(u8 *seed_matrix, s32 *vect, s32 *result, struct shake_ctx *shake, u8 *scratch_buffer, enum ml_dsa_level_k k, enum ml_dsa_level_l l);

// Stream matrix A from rho and multiply entries immediately
// No full matrix is materialized in memory
// (works with secrets or data dependent on a secret)
void mult_matrix(u8 *seed_matrix, s32 *vect, s32 *result, struct shake_ctx *shake, u8 *scratch_buffer, enum ml_dsa_level_k k, enum ml_dsa_level_l l)
{
	ml_dsa_memzero(result, k * ML_DSA_N * sizeof(s32));
	
	u8 curr_seed[ML_DSA_32_BYTES + 2];
	memcpy(curr_seed, seed_matrix, ML_DSA_32_BYTES);
	
	for(size_t i = 0; i < k; i++)
	{
		s32 *poly_r = result + i * ML_DSA_N;
		for(size_t j = 0; j < l; j++)
		{
			const s32 *poly_v = vect + j * ML_DSA_N;
			curr_seed[32] = (u8)j;
			curr_seed[33] = (u8)i;
			shake_ctx_zero(shake);
			shake_ctx_init(shake, scratch_buffer, ML_DSA_SIZE_SCRATCH_BUFFER, curr_seed, ML_DSA_32_BYTES + 2);
			shake128(shake);
			
			// Rejection-sample matrix coefficients and consume them immediately in the pointwise product
			size_t counter_in = 0;
            size_t counter_out = 0;
            while(counter_in < ML_DSA_N)
            {
				if(counter_out == ML_DSA_SIZE_SCRATCH_BUFFER) { shake128(shake); counter_out = 0; }
				u8 a = scratch_buffer[counter_out++]; 
				if(counter_out == ML_DSA_SIZE_SCRATCH_BUFFER) { shake128(shake); counter_out = 0; }
                u8 b = scratch_buffer[counter_out++]; 
                if(counter_out == ML_DSA_SIZE_SCRATCH_BUFFER) { shake128(shake); counter_out = 0; }
                u8 c = scratch_buffer[counter_out++]; 
                s32 z = a | ((s32)b << 8) | ((s32)(c & 0x7F) << 16);
				
				if(z < ML_DSA_Q)
				{
					poly_r[counter_in] += ml_dsa_montgomery_reduce((s64)z * poly_v[counter_in]);
					counter_in++;
				}
			}
		}
	}
}	

#if ML_DSA_FULL_MATRIX_BUFFER

void get_full_matrix(u8 *rho, enum ml_dsa_level_k k, enum ml_dsa_level_l l, struct shake_ctx *shake, s32 *matrix_buffer, u8 *scratch_buffer);
int mult_full_matrix(s32 *matrix_buffer, s32 *vect, s32 *result, enum ml_dsa_level_k k, enum ml_dsa_level_l l);

// Expand the full matrix A from rho and store it in matrix_buffer
// Matrix entries are generated as RejNTTPoly(rho || column || row)
void get_full_matrix(u8 *rho, enum ml_dsa_level_k k, enum ml_dsa_level_l l, struct shake_ctx *shake, s32 *matrix_buffer, u8 *scratch_buffer)
{
	u8 curr_seed[ML_DSA_32_BYTES + 2];
	memcpy(curr_seed, rho, ML_DSA_32_BYTES);
	curr_seed[32] = 0;
	curr_seed[33] = 0;
	
	for(size_t i = 0; i < k; i++)
	{
		for(size_t j = 0; j < l; j++)
		{
			curr_seed[32] = j;
			curr_seed[33] = i;
			shake_ctx_zero(shake);
			shake_ctx_init(shake, scratch_buffer, ML_DSA_SIZE_SCRATCH_BUFFER, curr_seed, ML_DSA_32_BYTES + 2);
			shake128(shake);
			
			s32 *poly = matrix_buffer + ((i * l + j) * ML_DSA_N);
			
			// Rejection-sample 23-bit values below q
			size_t counter_in = 0;
            size_t counter_out = 0;
			while(counter_in < ML_DSA_N)
			{
				if(counter_out == ML_DSA_SIZE_SCRATCH_BUFFER) { shake128(shake); counter_out = 0; }
				u8 a = scratch_buffer[counter_out++];
				if(counter_out == ML_DSA_SIZE_SCRATCH_BUFFER) { shake128(shake); counter_out = 0; }
                u8 b = scratch_buffer[counter_out++];
                if(counter_out == ML_DSA_SIZE_SCRATCH_BUFFER) { shake128(shake); counter_out = 0; }
                u8 c = scratch_buffer[counter_out++];
				s32 z = a | ((s32)b << 8) | ((s32)(c & 0x7F) << 16);
				
				if(z < ML_DSA_Q)
				{
					poly[counter_in++] = z;
				}	
					
			}
		}
	}
	shake_ctx_zero(shake);
	ml_dsa_memzero(scratch_buffer, ML_DSA_SIZE_SCRATCH_BUFFER);
}

// Multiply the pre-expanded matrix A by vect
// Both operands are interpreted in the NTT domain
// (works with secrets or data dependent on a secret)
int mult_full_matrix(s32 *matrix_buffer, s32 *vect, s32 *result, enum ml_dsa_level_k k, enum ml_dsa_level_l l)
{
	if(!matrix_buffer || !vect || !result) { return -ML_DSA_EINVAL; }
	
	if(k == ML_DSA_44_K && l != ML_DSA_44_L) { return -ML_DSA_EINVAL; }
	if(k == ML_DSA_65_K && l != ML_DSA_65_L) { return -ML_DSA_EINVAL; }
	if(k == ML_DSA_87_K && l != ML_DSA_87_L) { return -ML_DSA_EINVAL; }
	
	ml_dsa_memzero(result, k * ML_DSA_N * sizeof(s32));
	
	for (size_t i = 0; i < k; i++)
	{
		for (size_t j = 0; j < l; j++)
		{
			s32 *poly_m = matrix_buffer + ((i * l + j) * ML_DSA_N);
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
#endif //ML_DSA_FULL_MATRIX_BUFFER

