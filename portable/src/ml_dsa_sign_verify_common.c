// A module of internal functions used in both signing and verification - ml_dsa_sign_verify_common.c
#include "ml_dsa_core_header.h"

void get_mu(struct shake_ctx *shake, u8 *mu, u8 *tr, const u8 *msg, size_t msg_len, const u8 *context, size_t context_len, u8 prehash);
void ml_dsa_get_w1(const s32 *w, s32 *w1, const u8 *h, enum ml_dsa_level_k k);
void ml_dsa_get_poly_c(s32 *poly_c, u8 *c, struct shake_ctx *shake, enum ml_dsa_level_k k);

void get_mu(struct shake_ctx *shake, u8 *mu, u8 *tr, const u8 *msg, size_t msg_len, const u8 *context, size_t context_len, u8 prehash)
{
	// buffer for first three rate (136 * 3 = 408)
	u8 first_rate[408]; 
	ml_dsa_memzero(first_rate, 408);
	
	// copy hash pk, bytes prehash and size of contest
	memcpy(first_rate, tr, ML_DSA_64_BYTES); 
	first_rate[ML_DSA_64_BYTES] = (bool)prehash;        
	first_rate[ML_DSA_64_BYTES + 1] = (u8)context_len;
	
	// copy context
	size_t counter_ptr = ML_DSA_64_BYTES + 2;
	if(context_len) { memcpy(first_rate + counter_ptr, context, context_len); }
	counter_ptr += context_len;
	
	if(prehash != 0)
	{
		u8 last_byte;
		size_t hash_len;

		shake_ctx_zero(shake);
		
		if(prehash == 1) { last_byte = 0x0B; hash_len = 32; } 
		else if(prehash == 2) { last_byte = 0x0C; hash_len = 64; }
		else if(prehash == 3) { last_byte = 0x07; hash_len = 28; }
		else if(prehash == 4) { last_byte = 0x08; hash_len = 32; }
		else if(prehash == 5) { last_byte = 0x09; hash_len = 48; }
		else if(prehash == 6) { last_byte = 0x0A; hash_len = 64; }
		else if(prehash == 7) { last_byte = 0x04; hash_len = 28; }
		else if(prehash == 8) { last_byte = 0x01; hash_len = 32; }
		else { goto exit; } // Все перевірено до, тому тут функція void, вона внутрішня, це щоб компілятор попередженнями не смердів

		shake_ctx_init(shake, mu, hash_len, msg, msg_len);
		
		#if defined(SHA2_32)
		struct sha2_32_ctx sha2_ctx;
		#endif

		switch(prehash)
		{
			case 1: shake128(shake); break;
			case 2: shake256(shake); break;
			case 3: sha3_224(shake); break;
			case 4: sha3_256(shake); break;
			case 5: sha3_384(shake); break;
			case 6: sha3_512(shake); break;
			#if defined(SHA2_32)
			case 7: sha224(&sha2_ctx, msg, msg_len, mu); break;
			case 8: sha256(&sha2_ctx, msg, msg_len, mu); break;
			#endif
		}

		memcpy(first_rate + counter_ptr, oid_hash_base, ML_DSA_PREHASH_OID_SIZE);
		counter_ptr += ML_DSA_PREHASH_OID_SIZE;
		first_rate[counter_ptr++] = last_byte;

		memcpy(first_rate + counter_ptr, mu, hash_len);
		counter_ptr += hash_len;

		shake_ctx_zero(shake);
		shake_ctx_init(shake, mu, ML_DSA_64_BYTES, first_rate, counter_ptr);
		shake256(shake);
		
		#if defined(SHA2_32)
		ml_dsa_memzero(&sha2_ctx, sizeof(struct sha2_32_ctx));
		#endif
		
		goto exit;
	}
	
	if(408 - counter_ptr >= msg_len) 
	{
		if(msg_len) { memcpy(first_rate + counter_ptr, msg, msg_len); }
		counter_ptr += msg_len;
		
		shake_ctx_zero(shake);
		shake_ctx_init(shake, mu, ML_DSA_64_BYTES, first_rate, counter_ptr);
		shake256(shake);
		goto exit;
	}
	
	memcpy(first_rate + counter_ptr, msg, 408 - counter_ptr); 
	shake_ctx_zero(shake);
	shake_ctx_init(shake, mu, ML_DSA_64_BYTES, first_rate, 408);
	shake->rate = 136;
	shake_absorb(shake);
	counter_ptr = 408 - counter_ptr;
	shake->in = msg + counter_ptr;
	shake->inlen = msg_len - counter_ptr;
	shake256(shake);
	
	exit:
		ml_dsa_memzero(first_rate, 408);
		shake_ctx_zero(shake);
		return;
}

// (works with secrets or data dependent on a secret)
void ml_dsa_get_w1(const s32 *w, s32 *w1, const u8 *h, enum ml_dsa_level_k k)
{
	s32 m, gamma2;
	if(k == ML_DSA_44_K) { gamma2 = ML_DSA_44_GAMMA2; m = 44; } // For 44 level gamma2 and m = (q-1)/(2*𝛾2)
	else { gamma2 = ML_DSA_65_87_GAMMA2; m = 16; } // For 65 and 87 levels gamma2 and m = (q-1)/(2*𝛾2)

    for (size_t i = 0; i < k * ML_DSA_N; i++)
    {
        s32 r = w[i];
        s32 high, low;

        ml_dsa_decompose(&r, gamma2, &high, &low);

        if (h != NULL && h[i] == 1)
        {
            if (low > 0)
            {
                high = (high == m - 1) ? 0 : high + 1;
            }
            else
            {
                high = (high == 0) ? m - 1 : high - 1;
            }
        }

        w1[i] = high;
    }
}

//  SampleInBall
void ml_dsa_get_poly_c(s32 *poly_c, u8 *c, struct shake_ctx *shake, enum ml_dsa_level_k k)
{
	ml_dsa_memzero(poly_c, ML_DSA_N * sizeof(s32));
	
	// get τ
	u8 tau, size_c;
	if(k == ML_DSA_44_K) { tau = 39; size_c = 32; }
	else if(k == ML_DSA_65_K) { tau = 49; size_c = 48; }
	else if(k == ML_DSA_87_K) { tau = 60; size_c = 64; }
	else { return; }

	// 136 - size of rate in shake256
	u8 temp_buff[136];
	shake_ctx_zero(shake);
	shake_ctx_init(shake, temp_buff, 136, c, size_c);
	shake256(shake);
	
	// First 8 bytes contain sign bits
	u8 bit_sing[8];
	memcpy(bit_sing, temp_buff, 8);
	
	size_t pos = 8;

	for (u16 i = 256 - tau; i < 256; i++)
	{
		u8 j;
		// Rejection sampling: 0 <= j <= i
        do {
            if (pos == sizeof(temp_buff)) { shake256(shake); pos = 0; }
            j = temp_buff[pos++];
        } while (j > i);

		// Fisher-Yates
		poly_c[i] = poly_c[j];

		// h[i + tau - 256]
		u8 h_index = i + tau - 256;

		u8 h = (bit_sing[h_index >> 3] >> (h_index & 7)) & 1;

		poly_c[j] = 1 - ((s32)h << 1);
	}
}
