#ifndef SHAKE_H
#define SHAKE_H

#include "ml_dsa_core_header.h"

struct shake_ctx {
    u64 *state;
    size_t rate;
    size_t pos;

    const u8 *in;
    size_t inlen;

    u8 *out;
    size_t outlen;
};

int shake_ctx_get_mem(struct shake_ctx *ctx, u64 *state, size_t size_mem);
void shake_ctx_init(struct shake_ctx *ctx, u8 *out, size_t outlen, const u8 *in, size_t inlen);
void shake_ctx_zero(struct shake_ctx *ctx);
void shake128(struct shake_ctx *ctx);
void shake256(struct shake_ctx *ctx);

void shake_absorb(struct shake_ctx *ctx);
void shake_finalize(struct shake_ctx *ctx);
void shake_squeeze(struct shake_ctx *ctx);

void keccak_f1600(u64 s[25]);

#endif
