// SPDX-License-Identifier: GPL-2.0 OR MIT
#ifndef SHA2_32_H
#define SHA2_32_H

#include "ml_dsa_core_header.h"

#define SHA224_DIGEST_SIZE 28u
#define SHA256_DIGEST_SIZE 32u
#define SHA2_32_BLOCK_SIZE 64u

// Caller owns the context; no heap allocation is performed
struct sha2_32_ctx {
    u32 state[8];
    u32 schedule[16];
    u8 buffer[64];
    u64 total_len; // Number of message bytes fed to update()
    size_t buffered;
    u8 digest_len; // 28 for SHA-224, 32 for SHA-256
};

//  SHA-224/SHA-256 32-bit compression function (FIPS 180-4, sec. 6.2/6.3).
//  state:    eight 32-bit chaining words, updated in place.
//  block:    exactly 64 input bytes (512 bits), in standard byte order.
//  schedule: caller-owned scratch space of 16 uint32_t words (64 bytes).
//  The caller initializes state to the algorithm-specific IV, handles padding
//  and message length, and serializes the final digest in big-endian order.
//  This function processes one complete block only.
//  state and schedule must not overlap. All pointers must be valid.
//  No heap allocations, no global mutable state, no secret-dependent branches.
void sha2_32_compress(u32 state[8], const u8 block[64], u32 schedule[16]);

// Return 0 on success, -1 for invalid pointers/length or exceeded bit length.
// update() accepts input=NULL only when len=0.
// Once final() succeeds, the context is wiped and must be reinitialized.
// Input cannot overlap writable parts of ctx; output cannot overlap ctx.
// A context cannot be used concurrently by multiple threads.
void sha224_init(struct sha2_32_ctx *ctx);
void sha256_init(struct sha2_32_ctx *ctx);
int sha2_32_update(struct sha2_32_ctx *ctx, const u8 *input, size_t len);
int sha2_32_final(struct sha2_32_ctx *ctx, u8 *digest, size_t digest_size);

// One-shot: ctx is caller-owned scratch memory and wiped on success
int sha224(struct sha2_32_ctx *ctx, const u8 *input, size_t len, u8 digest[SHA224_DIGEST_SIZE]);
int sha256(struct sha2_32_ctx *ctx, const u8 *input, size_t len, u8 digest[SHA256_DIGEST_SIZE]);

#endif
