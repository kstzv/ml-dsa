// ML_DSA_H
#ifndef ML_DSA_H
#define ML_DSA_H

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

#endif // ML_DSA_H

