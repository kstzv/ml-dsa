// SPDX-License-Identifier: GPL-2.0 OR MIT
#ifndef SHA2_64_H
#define SHA2_64_H

#include "ml_dsa_core_header.h"

#define SHA384_DIGEST_SIZE     48u
#define SHA512_DIGEST_SIZE     64u
#define SHA512_224_DIGEST_SIZE 28u
#define SHA512_256_DIGEST_SIZE 32u
#define SHA2_64_BLOCK_SIZE   128u

// Caller owns the context; no heap allocation is performed.
struct sha2_64_ctx {
    u64 state[8];
    u64 schedule[16];
    u8 buffer[SHA2_64_BLOCK_SIZE];
    u64 total_len_hi; // High 64 bits of the number of message bytes.
    u64 total_len_lo; // Low 64 bits of the number of message bytes.
    size_t buffered;
    u8 digest_len; // 48 (384), 64 (512), 28 (512/224), or 32 (512/256).
};

// SHA-384/SHA-512/SHA-512/224/SHA-512/256 compression (FIPS 180-4).
// state:    eight 64-bit chaining words, updated in place.
// block:    exactly 128 input bytes (1024 bits), in standard byte order.
// schedule: caller-owned scratch space of sixteen u64 words (128 bytes).
// The caller sets the algorithm-specific IV, handles padding and the 128-bit
// message length, and serializes the final digest in big-endian byte order.
// This function processes one complete block only.
// state and schedule must not overlap. All pointers must be valid.
// No heap allocation, no global mutable state, no secret-dependent branches.
void sha2_64_compress(u64 state[8], const u8 block[SHA2_64_BLOCK_SIZE], u64 schedule[16]);

// Return 0 on success, -1 for invalid parameters or exceeded bit length.
// update() permits input=NULL only when len=0.
// The maximum supported input length is 2^125 - 1 bytes (FIPS 180-4).
// After final() succeeds the context is wiped and must be reinitialized.
// Input cannot overlap writable parts of ctx; output cannot overlap ctx.
// A context cannot be used concurrently by multiple threads.
void sha384_init(struct sha2_64_ctx *ctx);
void sha512_init(struct sha2_64_ctx *ctx);
void sha512_224_init(struct sha2_64_ctx *ctx);
void sha512_256_init(struct sha2_64_ctx *ctx);
int sha2_64_update(struct sha2_64_ctx *ctx, const u8 *input, size_t len);
int sha2_64_final(struct sha2_64_ctx *ctx, u8 *digest, size_t digest_size);

// One-shot: ctx is caller-owned scratch memory, wiped on success and on an
// update() failure. Output size is fixed by the chosen hash function.
int sha384(struct sha2_64_ctx *ctx, const u8 *input, size_t len, u8 digest[SHA384_DIGEST_SIZE]);
int sha512(struct sha2_64_ctx *ctx, const u8 *input, size_t len, u8 digest[SHA512_DIGEST_SIZE]);
int sha512_224(struct sha2_64_ctx *ctx, const u8 *input, size_t len, u8 digest[SHA512_224_DIGEST_SIZE]);
int sha512_256(struct sha2_64_ctx *ctx, const u8 *input, size_t len, u8 digest[SHA512_256_DIGEST_SIZE]);

#endif // SHA2_64_H
