#ifndef SHAKE_H
#define SHAKE_H

#include <stdint.h>
#include <stddef.h>

struct shake_ctx {
    uint64_t *state;
    size_t rate;
    size_t pos;

    const uint8_t *in;
    size_t inlen;

    uint8_t *out;
    size_t outlen;
};

int shake_ctx_get_mem(struct shake_ctx *ctx, uint64_t *state, size_t size_mem)
void shake_ctx_init(struct shake_ctx *ctx,uint8_t *out, size_t outlen, const uint8_t *in, size_t inlen);
void shake_ctx_zero(struct shake_ctx *ctx);
void shake128(struct shake_ctx *ctx);
void shake256(struct shake_ctx *ctx);

void keccak_f1600(uint64_t s[25]);

#endif
