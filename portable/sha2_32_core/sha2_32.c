// SPDX-License-Identifier: GPL-2.0 OR MIT
#include "sha2_32.h"

static const u32 sha224_iv[8] = {
    UINT32_C(0xc1059ed8), UINT32_C(0x367cd507),
    UINT32_C(0x3070dd17), UINT32_C(0xf70e5939),
    UINT32_C(0xffc00b31), UINT32_C(0x68581511),
    UINT32_C(0x64f98fa7), UINT32_C(0xbefa4fa4)
};

static const u32 sha256_iv[8] = {
    UINT32_C(0x6a09e667), UINT32_C(0xbb67ae85),
    UINT32_C(0x3c6ef372), UINT32_C(0xa54ff53a),
    UINT32_C(0x510e527f), UINT32_C(0x9b05688c),
    UINT32_C(0x1f83d9ab), UINT32_C(0x5be0cd19)
};

static void wipe(void *p, size_t len)
{
    volatile u8 *v = (volatile u8 *)p;
    while (len--) *v++ = 0;
}

static void init(struct sha2_32_ctx *ctx, const u32 iv[8], u8 digest_len)
{
    if (!ctx) return;
    memset(ctx, 0, sizeof(*ctx));
    memcpy(ctx->state, iv, sizeof(ctx->state));
    ctx->digest_len = digest_len;
}

void sha224_init(struct sha2_32_ctx *ctx) { init(ctx, sha224_iv, SHA224_DIGEST_SIZE); }
void sha256_init(struct sha2_32_ctx *ctx) { init(ctx, sha256_iv, SHA256_DIGEST_SIZE); }

int sha2_32_update(struct sha2_32_ctx *ctx, const u8 *input, size_t len)
{
    if (!ctx || (!input && len) || (ctx->digest_len != SHA224_DIGEST_SIZE && ctx->digest_len != SHA256_DIGEST_SIZE)) { return -1; }
    if (len > UINT64_MAX / UINT64_C(8) - ctx->total_len) { return -1; }

    ctx->total_len += (u64)len;
    if (ctx->buffered) 
    {
        size_t take = SHA2_32_BLOCK_SIZE - ctx->buffered;
        if(take > len) { take = len; }
        if(take) { memcpy(ctx->buffer + ctx->buffered, input, take); }
        ctx->buffered += take;
        if(take) { input += take; }
        len -= take;
        if(ctx->buffered == SHA2_32_BLOCK_SIZE)
        {
            sha2_32_compress(ctx->state, ctx->buffer, ctx->schedule);
            ctx->buffered = 0;
        }
    }
    while (len >= SHA2_32_BLOCK_SIZE) 
    {
        sha2_32_compress(ctx->state, input, ctx->schedule);
        input += SHA2_32_BLOCK_SIZE;
        len -= SHA2_32_BLOCK_SIZE;
    }
    if(len) 
    {
        memcpy(ctx->buffer, input, len);
        ctx->buffered = len;
    }
    return 0;
}

int sha2_32_final(struct sha2_32_ctx *ctx, u8 *digest, size_t digest_size)
{
    if (!ctx || !digest || (ctx->digest_len != SHA224_DIGEST_SIZE && ctx->digest_len != SHA256_DIGEST_SIZE) 
		|| digest_size < ctx->digest_len || ctx->buffered >= SHA2_32_BLOCK_SIZE) { return -1; }

    const u64 bit_len = ctx->total_len * UINT64_C(8);
    size_t n = ctx->buffered;
    ctx->buffer[n++] = UINT8_C(0x80);
    if(n > 56u) 
    {
        memset(ctx->buffer + n, 0, SHA2_32_BLOCK_SIZE - n);
        sha2_32_compress(ctx->state, ctx->buffer, ctx->schedule);
        n = 0;
    }
    
    memset(ctx->buffer + n, 0, 56u - n);
    
    for (unsigned j = 0; j < 8; ++j) { ctx->buffer[56u + j] = (u8)(bit_len >> (56u - 8u * j)); }
    
    sha2_32_compress(ctx->state, ctx->buffer, ctx->schedule);

    for (size_t i = 0; i < ctx->digest_len / 4u; ++i)
    {
        const u32 w = ctx->state[i];
        digest[4u * i]     = (u8)(w >> 24);
        digest[4u * i + 1] = (u8)(w >> 16);
        digest[4u * i + 2] = (u8)(w >> 8);
        digest[4u * i + 3] = (u8)w;
    }
    wipe(ctx, sizeof(*ctx));
    return 0;
}

int sha224(struct sha2_32_ctx *ctx, const u8 *input, size_t len, u8 digest[SHA224_DIGEST_SIZE])
{
    if(!ctx || !digest || (!input && len)) { return -1; }
    sha224_init(ctx);
    if (sha2_32_update(ctx, input, len) != 0) { wipe(ctx, sizeof(*ctx)); return -1; }
    return sha2_32_final(ctx, digest, SHA224_DIGEST_SIZE);
}

int sha256(struct sha2_32_ctx *ctx, const u8 *input, size_t len, u8 digest[SHA256_DIGEST_SIZE])
{
    if(!ctx || !digest || (!input && len)) { return -1; }
    sha256_init(ctx);
    if(sha2_32_update(ctx, input, len) != 0) { wipe(ctx, sizeof(*ctx)); return -1; }
    return sha2_32_final(ctx, digest, SHA256_DIGEST_SIZE);
}
