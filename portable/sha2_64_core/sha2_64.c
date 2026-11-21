// SPDX-License-Identifier: GPL-2.0 OR MIT
#include "sha2_64.h"

static const u64 sha384_iv[8] = {
    UINT64_C(0xcbbb9d5dc1059ed8), UINT64_C(0x629a292a367cd507),
    UINT64_C(0x9159015a3070dd17), UINT64_C(0x152fecd8f70e5939),
    UINT64_C(0x67332667ffc00b31), UINT64_C(0x8eb44a8768581511),
    UINT64_C(0xdb0c2e0d64f98fa7), UINT64_C(0x47b5481dbefa4fa4)
};

static const u64 sha512_iv[8] = {
    UINT64_C(0x6a09e667f3bcc908), UINT64_C(0xbb67ae8584caa73b),
    UINT64_C(0x3c6ef372fe94f82b), UINT64_C(0xa54ff53a5f1d36f1),
    UINT64_C(0x510e527fade682d1), UINT64_C(0x9b05688c2b3e6c1f),
    UINT64_C(0x1f83d9abfb41bd6b), UINT64_C(0x5be0cd19137e2179)
};

// SHA-512/224 and SHA-512/256 have their own IVs; neither is just a
// truncation of ordinary SHA-512.
static const u64 sha512_224_iv[8] = {
    UINT64_C(0x8c3d37c819544da2), UINT64_C(0x73e1996689dcd4d6),
    UINT64_C(0x1dfab7ae32ff9c82), UINT64_C(0x679dd514582f9fcf),
    UINT64_C(0x0f6d2b697bd44da8), UINT64_C(0x77e36f7304c48942),
    UINT64_C(0x3f9d85a86a1d36c8), UINT64_C(0x1112e6ad91d692a1)
};

static const u64 sha512_256_iv[8] = {
    UINT64_C(0x22312194fc2bf72c), UINT64_C(0x9f555fa3c84c64c2),
    UINT64_C(0x2393b86b6f53b151), UINT64_C(0x963877195940eabd),
    UINT64_C(0x96283ee2a88effe3), UINT64_C(0xbe5e1e2553863992),
    UINT64_C(0x2b0199fc2c85b8aa), UINT64_C(0x0eb72ddc81c52ca2)
};

static void wipe(void *p, size_t len)
{
    volatile u8 *v = (volatile u8 *)p;
    while (len--) *v++ = 0;
}

static void init(struct sha2_64_ctx *ctx, const u64 iv[8], u8 digest_len)
{
    if (!ctx) return;
    memset(ctx, 0, sizeof(*ctx));
    memcpy(ctx->state, iv, sizeof(ctx->state));
    ctx->digest_len = digest_len;
}

void sha384_init(struct sha2_64_ctx *ctx) { init(ctx, sha384_iv, SHA384_DIGEST_SIZE); }
void sha512_init(struct sha2_64_ctx *ctx) { init(ctx, sha512_iv, SHA512_DIGEST_SIZE); }
void sha512_224_init(struct sha2_64_ctx *ctx) { init(ctx, sha512_224_iv, SHA512_224_DIGEST_SIZE); }
void sha512_256_init(struct sha2_64_ctx *ctx) { init(ctx, sha512_256_iv, SHA512_256_DIGEST_SIZE); }

static int valid_digest_len(u8 n)
{
    return n == SHA384_DIGEST_SIZE || n == SHA512_DIGEST_SIZE ||
           n == SHA512_224_DIGEST_SIZE || n == SHA512_256_DIGEST_SIZE;
}

int sha2_64_update(struct sha2_64_ctx *ctx, const u8 *input, size_t len)
{
    if (!ctx || (!input && len) || !valid_digest_len(ctx->digest_len) ||
        ctx->buffered >= SHA2_64_BLOCK_SIZE) { return -1; }

    // Track byte count as 128 bits; make sure multiplication by eight at
    // finalization fits in the 128-bit FIPS 180-4 bit-length field.
    const u64 add = (u64)len;
    const u64 next_lo = ctx->total_len_lo + add;
    const u64 carry = (next_lo < ctx->total_len_lo);
    if (ctx->total_len_hi >= (UINT64_C(1) << 61) ||
        ctx->total_len_hi + carry >= (UINT64_C(1) << 61)) { return -1; }

    ctx->total_len_lo = next_lo;
    ctx->total_len_hi += carry;

    if (ctx->buffered)
    {
        size_t take = SHA2_64_BLOCK_SIZE - ctx->buffered;
        if (take > len) { take = len; }
        if (take) { memcpy(ctx->buffer + ctx->buffered, input, take); }
        ctx->buffered += take;
        if (take) { input += take; }
        len -= take;
        if (ctx->buffered == SHA2_64_BLOCK_SIZE)
        {
            sha2_64_compress(ctx->state, ctx->buffer, ctx->schedule);
            ctx->buffered = 0;
        }
    }
    while (len >= SHA2_64_BLOCK_SIZE)
    {
        sha2_64_compress(ctx->state, input, ctx->schedule);
        input += SHA2_64_BLOCK_SIZE;
        len -= SHA2_64_BLOCK_SIZE;
    }
    if (len)
    {
        memcpy(ctx->buffer, input, len);
        ctx->buffered = len;
    }
    return 0;
}

int sha2_64_final(struct sha2_64_ctx *ctx, u8 *digest, size_t digest_size)
{
    if (!ctx || !digest || !valid_digest_len(ctx->digest_len) ||
        digest_size < ctx->digest_len || ctx->buffered >= SHA2_64_BLOCK_SIZE ||
        ctx->total_len_hi >= (UINT64_C(1) << 61)) { return -1; }

    const u64 bits_hi = (ctx->total_len_hi << 3) | (ctx->total_len_lo >> 61);
    const u64 bits_lo = ctx->total_len_lo << 3;
    size_t n = ctx->buffered;
    ctx->buffer[n++] = UINT8_C(0x80);
    if (n > 112u)
    {
        memset(ctx->buffer + n, 0, SHA2_64_BLOCK_SIZE - n);
        sha2_64_compress(ctx->state, ctx->buffer, ctx->schedule);
        n = 0;
    }

    memset(ctx->buffer + n, 0, 112u - n);
    for (unsigned j = 0; j < 8; ++j)
    {
        ctx->buffer[112u + j] = (u8)(bits_hi >> (56u - 8u * j));
        ctx->buffer[120u + j] = (u8)(bits_lo >> (56u - 8u * j));
    }
    sha2_64_compress(ctx->state, ctx->buffer, ctx->schedule);

    // Byte-wise output also handles the partial last word of SHA-512/224.
    for (size_t i = 0; i < ctx->digest_len; ++i)
    {
        digest[i] = (u8)(ctx->state[i / 8u] >> (56u - 8u * (i & 7u)));
    }
    wipe(ctx, sizeof(*ctx));
    return 0;
}

int sha384(struct sha2_64_ctx *ctx, const u8 *input, size_t len, u8 digest[SHA384_DIGEST_SIZE])
{
    if (!ctx || !digest || (!input && len)) { return -1; }
    sha384_init(ctx);
    if (sha2_64_update(ctx, input, len) != 0) { wipe(ctx, sizeof(*ctx)); return -1; }
    return sha2_64_final(ctx, digest, SHA384_DIGEST_SIZE);
}

int sha512(struct sha2_64_ctx *ctx, const u8 *input, size_t len, u8 digest[SHA512_DIGEST_SIZE])
{
    if (!ctx || !digest || (!input && len)) { return -1; }
    sha512_init(ctx);
    if (sha2_64_update(ctx, input, len) != 0) { wipe(ctx, sizeof(*ctx)); return -1; }
    return sha2_64_final(ctx, digest, SHA512_DIGEST_SIZE);
}

int sha512_224(struct sha2_64_ctx *ctx, const u8 *input, size_t len, u8 digest[SHA512_224_DIGEST_SIZE])
{
    if (!ctx || !digest || (!input && len)) { return -1; }
    sha512_224_init(ctx);
    if (sha2_64_update(ctx, input, len) != 0) { wipe(ctx, sizeof(*ctx)); return -1; }
    return sha2_64_final(ctx, digest, SHA512_224_DIGEST_SIZE);
}

int sha512_256(struct sha2_64_ctx *ctx, const u8 *input, size_t len, u8 digest[SHA512_256_DIGEST_SIZE])
{
    if (!ctx || !digest || (!input && len)) { return -1; }
    sha512_256_init(ctx);
    if (sha2_64_update(ctx, input, len) != 0) { wipe(ctx, sizeof(*ctx)); return -1; }
    return sha2_64_final(ctx, digest, SHA512_256_DIGEST_SIZE);
}
