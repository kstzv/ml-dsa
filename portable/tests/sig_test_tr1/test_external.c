// SPDX-License-Identifier: GPL-2.0 OR MIT
// Copyright (c) 2026 K.S.Zavertailo
#include <stdio.h>
#include "ml_dsa_userspace_port.h"
#include "ml_dsa.h"
#include <jansson.h>
#include <limits.h>
#include <string.h>

size_t buffer_size = 10000000;
#define MSG_MAX_SIZE    8192
#define NUMBERS_VECTORS 720
#define NUMBERS_BLOCKS  48

u8 mes[MSG_MAX_SIZE];
u8 seed_sk[ML_DSA_87_SIZE_SK];
u8 context_m[ML_DSA_SIZE_MAX_CONTEXT];
u8 temp_buff_sig[ML_DSA_87_SIZE_SIG];
u8 temp_buff_sig_expected[ML_DSA_87_SIZE_SIG];

u8 test_rnd[ML_DSA_32_BYTES];

static int hexval(int c);
int read_hex_byte(const char *in, u8 *out, size_t size_out_mass);

int main()
{
    json_error_t error;

    json_t *prompt = json_load_file("prompt.json", 0, &error);
    if (!prompt) { return -1; }

    json_t *expected = json_load_file("expectedResults.json", 0, &error);
    if (!expected) { json_decref(prompt); return -1; }

    // Кореневі поля - поки закоментовані, мб пізніше знадобляться при переробці коду тесту
    // const char *algorithm = json_string_value(json_object_get(prompt, "algorithm"));
    // const char *mode = json_string_value(json_object_get(prompt, "mode"));
    // const char *revision = json_string_value(json_object_get(prompt, "revision"));
    // json_int_t vsId = json_integer_value(json_object_get(prompt, "vsId"));

    // Масиви груп обох файлів
    json_t *groups = json_object_get(prompt, "testGroups");
    json_t *expected_groups = json_object_get(expected, "testGroups");

    size_t group_count = json_array_size(groups);

    for (size_t i = 0; i < group_count; i++) 
    {
        // Поточна група в кожному файлі
        json_t *group = json_array_get(groups, i);
        json_t *expected_group = json_array_get(expected_groups, i);

        // Поля поточної групи
        json_int_t tgId = json_integer_value(json_object_get(group, "tgId"));

        // const char *testType = json_string_value(json_object_get(group, "testType"));
        const char *parameterSet = json_string_value(json_object_get(group, "parameterSet"));
        // const char *keyFormat = json_string_value(json_object_get(group, "keyFormat"));
        const char *signatureInterface = json_string_value(json_object_get(group, "signatureInterface"));
        // const char *preHash = json_string_value(json_object_get(group, "preHash"));

        bool deterministic = json_is_true(json_object_get(group, "deterministic"));
        // bool externalMu = json_is_true(json_object_get(group, "externalMu"));
        
        if(strstr(signatureInterface, "internal")) { continue; }
        
        u8 level;
        size_t len_mes, len_context, len_sig, len_sk;
        if(strstr(parameterSet, "ML-DSA-44")) { level = 44; len_sig = ML_DSA_44_SIZE_SIG; len_sk = ML_DSA_44_SIZE_SK; }
        else if(strstr(parameterSet, "ML-DSA-65")) { level = 65; len_sig = ML_DSA_65_SIZE_SIG; len_sk = ML_DSA_65_SIZE_SK; }
        else if(strstr(parameterSet, "ML-DSA-87")) { level = 87; len_sig = ML_DSA_87_SIZE_SIG; len_sk = ML_DSA_87_SIZE_SK; }
        else { json_decref(prompt); json_decref(expected); return -1; }

        // Масиви тестів поточної групи
        json_t *tests = json_object_get(group, "tests");
        json_t *expected_tests = json_object_get(expected_group, "tests");

        size_t test_count = json_array_size(tests);

        for (size_t j = 0; j < test_count; j++) 
        {
            // Поточний тест в обох файлах
            json_t *test = json_array_get(tests, j);
            json_t *expected_test = json_array_get(expected_tests, j);

            // Ідентифікатори тесту
            json_int_t tcId = json_integer_value(json_object_get(test, "tcId"));
            // json_int_t expected_tcId = json_integer_value(json_object_get(expected_test, "tcId"));

            // Вхідні дані
            const char *message = json_string_value(json_object_get(test, "message"));
            const char *seed = json_string_value(json_object_get(test, "seed"));
            const char *sk = json_string_value(json_object_get(test, "sk"));
            const char *context = json_string_value(json_object_get(test, "context"));
            const char *hashAlg = json_string_value(json_object_get(test, "hashAlg"));
            // const char *mu = json_string_value(json_object_get(test, "mu"));
            const char *rnd = json_string_value(json_object_get(test, "rnd"));
            
            u8 prehash;
            if(hashAlg && strstr(hashAlg, "SHA2-512")) { continue; }
            else if(hashAlg && strstr(hashAlg, "SHA2-512/224")) { continue; }
            else if(hashAlg && strstr(hashAlg, "SHA2-256")) { continue; }
            else if(hashAlg && strstr(hashAlg, "SHA2-224")) { continue; }
            else if(hashAlg && strstr(hashAlg, "SHA2-512/256")) { continue; }
            else if(hashAlg && strstr(hashAlg, "SHA2-384")) { continue; } 
			else if(hashAlg && strstr(hashAlg, "SHAKE-128")) { prehash = 1; }
            else if(hashAlg && strstr(hashAlg, "SHAKE-256")) { prehash = 2; }
            else if(hashAlg && strstr(hashAlg, "SHA3-224")) { prehash = 3; }
            else if(hashAlg && strstr(hashAlg, "SHA3-256")) { prehash = 4; }
            else if(hashAlg && strstr(hashAlg, "SHA3-384")) { prehash = 5; }
            else if(hashAlg && strstr(hashAlg, "SHA3-512")) { prehash = 6; }
            else { prehash = 0; }

            // Очікуваний результат
            const char *signature = json_string_value(json_object_get(expected_test, "signature"));
            
            int ret;
            
            if(seed) 
            { 
				len_sk = ML_DSA_32_BYTES; 
				ret = read_hex_byte(seed, seed_sk, len_sk);
				if(ret != (int)len_sk) { printf("Err in read_hex_byte()"); json_decref(prompt); json_decref(expected); return -1; }
			}else if(sk)
			{
				ret = read_hex_byte(sk, seed_sk, len_sk);
				if(ret != (int)len_sk) { printf("Err in read_hex_byte()"); json_decref(prompt); json_decref(expected); return -1; }
			}else { printf("Err in json-file"); json_decref(prompt); json_decref(expected); return -1; }
			
			if(!rnd) { ml_dsa_memzero(test_rnd, ML_DSA_32_BYTES); }
			else 
			{ 
				ret = read_hex_byte(rnd, test_rnd, ML_DSA_32_BYTES);
				if(ret != ML_DSA_32_BYTES) { printf("Err in read_hex_byte()"); json_decref(prompt); json_decref(expected); return -1; }
			}
            
			ret = read_hex_byte(message, mes, MSG_MAX_SIZE);
			if(ret < 0) { printf("Err in read_hex_byte()"); json_decref(prompt); json_decref(expected); return -1; }
			else { len_mes = (size_t)ret; }
			
			ret = read_hex_byte(context, context_m, ML_DSA_SIZE_MAX_CONTEXT);
			if(ret < 0) { printf("Err in read_hex_byte()"); json_decref(prompt); json_decref(expected); return -1; }
			else { len_context = (size_t)ret; }
				
			struct ml_dsa_keys *ctx = ml_dsa_create_keys_object_from_imported_sk(seed_sk, len_sk, level);
			if(!ctx) { printf("Err in ml_dsa_create_keys_object_from_imported_sk()"); json_decref(prompt); json_decref(expected); return -1; }
				
			ret = ml_dsa_sign(ctx, mes, len_mes, context_m, len_context, temp_buff_sig, len_sig, prehash, deterministic, NULL);
			if(ret != 0) { printf("Err in ml_dsa_sign()"); json_decref(prompt); json_decref(expected); ml_dsa_destroy_struct_keys(ctx); return -1; }
				
			ret = read_hex_byte(signature, temp_buff_sig_expected, len_sig);
			if(ret < 0 || (size_t)ret != len_sig) { json_decref(prompt); json_decref(expected); ml_dsa_destroy_struct_keys(ctx); return -1; }
				
			// Compare generated public key with NIST reference. Full byte-wise comparison ensures deterministic correctness.
			if(memcmp(temp_buff_sig_expected, temp_buff_sig, len_sig) != 0) { printf("\nWrong! In %zu vector in %zu block, in ML-DSA-%u level sequrity!", i, j, level); }
			else{ printf("\nDone.... In %zu vector in %zu block, in ML-DSA-%u level sequrity is identical with vector NIST....", j, i, level); }
				
			ml_dsa_destroy_struct_keys(ctx);
			

            // Тут усі дані поточного вектора доступні.
            printf("tgId = %lld, tcId = %lld\n", (long long)tgId, (long long)tcId);

        }
    }

    json_decref(prompt);
    json_decref(expected);

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
    if (!in || !out) { return -1; }

    size_t len = strlen(in);

    if (len % 2 != 0) { return -1; }
    if (len / 2 > size_out_mass) { return -1; }
    if (len / 2 > INT_MAX) { return -1; }

    size_t out_counter = 0;

    while (out_counter < len / 2)
    {
        int high = hexval((unsigned char)in[out_counter * 2]);
        int low  = hexval((unsigned char)in[out_counter * 2 + 1]);

        if (high < 0 || low < 0) { return -1; }

        out[out_counter] = (u8)((high << 4) | low);

        out_counter++;
    }

    return (int)out_counter;
}
