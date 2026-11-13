// SPDX-License-Identifier: GPL-2.0 OR MIT
// Copyright (c) 2026 K.S.Zavertailo
#include <stdio.h>
#include "ml_dsa_userspace_port.h"
#include "ml_dsa.h"

#define SIZE_BUFFER 16384
#define NUMBERS_VECTORS 75

u8 mass_seed[ML_DSA_32_BYTES];
u8 mass_pk[NUMBERS_VECTORS][ML_DSA_87_SIZE_PK];
u8 mass_sk[NUMBERS_VECTORS][ML_DSA_87_SIZE_SK];

u8 temp_mass_pk[ML_DSA_87_SIZE_PK];
u8 temp_mass_sk[ML_DSA_87_SIZE_SK];

static int hexval(int c);
int read_hex_byte(const char *in, u8 *out, size_t size_out_mass);

int main()
{
	// Buffer for reading lines from JSON file
	char buffer[SIZE_BUFFER];
	size_t ret;
	
	// Open NIST test vector file (KeyGen vectors).
	FILE *f = fopen("expectedResults.json", "r");
	if(!f) { return -1; }
	
	size_t size_pk = ML_DSA_44_SIZE_PK;
	size_t size_sk = ML_DSA_44_SIZE_SK;
	u32 i = 0;
	
	while(fgets(buffer, sizeof(buffer), f))
	{
		if(i == 25)
		{
			size_pk = ML_DSA_65_SIZE_PK;
			size_sk = ML_DSA_65_SIZE_SK;
		}
		
		if(i == 50)
		{
			size_pk = ML_DSA_87_SIZE_PK;
			size_sk = ML_DSA_87_SIZE_SK;
		}
		
		if(i >= SIZE_BUFFER) { break; }
		
		if(strstr(buffer, "\"pk\": \""))
		{
			ret = read_hex_byte(buffer, mass_pk[i], size_pk);
			if(ret != size_pk) { return -1; }
		}
		
		if(strstr(buffer, "\"sk\": \""))
		{
			ret = read_hex_byte(buffer, mass_sk[i], size_sk);
			if(ret != size_sk) { return -1; }
			i++;
		}
	}
	
	fclose(f);
	i = 0;
	size_pk = ML_DSA_44_SIZE_PK;
	size_sk = ML_DSA_44_SIZE_SK;
	u8 level = 44;
	
	// Open NIST test vector file (KeyGen vectors).
	f = fopen("prompt.json", "r");
	if(!f) { return -1; }
	
	// Skip initial metadata until first "parameterSet" entry.
	while(fgets(buffer, sizeof(buffer), f))
	{
		if(strstr(buffer, "parameterSet")) { break; }
	}
	
	printf("<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>\n");
	while(fgets(buffer, sizeof(buffer), f))
	{
		if(i >= NUMBERS_VECTORS) { break; }
		else if(strstr(buffer, "\"parameterSet\": \"ML-DSA-65\""))
		{
			size_pk = ML_DSA_65_SIZE_PK;
			size_sk = ML_DSA_65_SIZE_SK;
			level = 65;
			printf("<><><><><><><><><><><><><><><><><><><><><><><>><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>\n");
		}else if(strstr(buffer, "\"parameterSet\": \"ML-DSA-87\""))
		{
			size_pk = ML_DSA_87_SIZE_PK;
			size_sk = ML_DSA_87_SIZE_SK;
			level = 87;
			printf("<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><<><><><><><><><><><><><><><><><><><><><><><><><><><>\n");
		}else if(strstr(buffer, "\"seed\": \""))
		{
			ret = read_hex_byte(buffer, mass_seed, ML_DSA_32_BYTES);
			if(ret != ML_DSA_32_BYTES) { return -1; }
			
			struct ml_dsa_keys *ctx = ml_dsa_create_keys_object_from_imported_sk(mass_seed, ML_DSA_32_BYTES, level);
			if(!ctx) { printf("error in ml_dsa_create_keys_object_from_imported_sk()\n"); return -1; }
			
			ret = ml_dsa_get_pk(ctx, temp_mass_pk, size_pk);
			if(ret != size_pk) { ml_dsa_destroy_struct_keys(ctx); printf("error in ml_dsa_get_pk()\n"); return -1; }
			
			ret = ml_dsa_get_sk(ctx, temp_mass_sk, size_sk);
			if(ret != size_sk) { ml_dsa_destroy_struct_keys(ctx); printf("error in ml_dsa_get_sk()\n"); return -1; }
			
			if(memcmp(mass_pk[i], temp_mass_pk, size_pk) != 0)
			{
				printf("\nWrong! In %u pk vector, in ML-DSA-%u level sequrity!", i, level);
			}else
			{
				printf("\nDone.... In %u pk vector, in ML-DSA-%u level sequrity is identical with vector NIST....", i, level);
			}
			
			if(memcmp(mass_sk[i], temp_mass_sk, size_sk) != 0)
			{
				printf("\nWrong! In %u sk vector, in ML-DSA-%u level sequrity!", i, level);
			}else
			{
				printf("\nDone.... In %u sk vector, in ML-DSA-%u level sequrity is identical with vector NIST....", i, level);
			}
			
			ml_dsa_destroy_struct_keys(ctx);
			i++;
			printf("\n");
		}
	}
	fclose(f);
	printf("<><><><><><><><><><><><><><><><><><><><><><><><><><><<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>\n");
	
	return 0;
}

// Convert single HEX character to numeric value. 
// Returns:  0..15 on success or -1 on invalid input
static int hexval(int c)
{
	if(c >= '0' && c <= '9') { return c - '0'; }
	else if(c >= 'a' && c <= 'f') { return c - 'a' + 10; }
	else if(c >= 'A' && c <= 'F') { return c - 'A' + 10; }
	else { return -1; }
}

// Parse HEX string from JSON line into byte array.
// Expected format: "key": "A1B2C3..."
// Behavior: 1. Skips until third quote (start of HEX data)
// 2. Converts each pair of HEX chars into one byte
// 3. Performs validation on input format
// Returns: number of bytes written or -1 on error	
int read_hex_byte(const char *in, u8 *out, size_t size_out_mass)
{
	if(!in || !out) { return -1; }
	
	int symb;
	size_t in_counter = 0;
	
	// Skip until HEX data start (after third quote)
	size_t dq_counter = 0;
	while(dq_counter != 3 && in_counter < SIZE_BUFFER)
	{
		symb = in[in_counter];
		if(symb == '"') { dq_counter++; }
		in_counter++;
	}
	
	size_t out_counter = 0;
	while(out_counter < size_out_mass && in_counter < SIZE_BUFFER)
	{
		// High nibble
		symb = in[in_counter];
		if(symb == '"') { break; }
		if(symb == '\0' || symb == ' ' || symb == '\n') { return -1; }
		symb = hexval(symb);
		if(symb < 0) { return symb; }
		
		out[out_counter] = symb;
		out[out_counter] <<= 4;
		in_counter++;
		
		// Low nibble
		symb = in[in_counter];
		if(symb == '\0' || symb == ' ' || symb == '\n' || symb == '"') { return -1; }
		symb = hexval(symb);
		if(symb < 0) { return symb; }
		
		out[out_counter] |= symb;
		
		out_counter++;
		in_counter++;
	}
	return out_counter;
}
