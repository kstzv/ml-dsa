// SPDX-License-Identifier: GPL-2.0 OR MIT
#include "sha2_32.h"

static const u32 sha2_32_k[64] = {
    UINT32_C(0x428a2f98), UINT32_C(0x71374491), UINT32_C(0xb5c0fbcf), UINT32_C(0xe9b5dba5),
    UINT32_C(0x3956c25b), UINT32_C(0x59f111f1), UINT32_C(0x923f82a4), UINT32_C(0xab1c5ed5),
    UINT32_C(0xd807aa98), UINT32_C(0x12835b01), UINT32_C(0x243185be), UINT32_C(0x550c7dc3),
    UINT32_C(0x72be5d74), UINT32_C(0x80deb1fe), UINT32_C(0x9bdc06a7), UINT32_C(0xc19bf174),
    UINT32_C(0xe49b69c1), UINT32_C(0xefbe4786), UINT32_C(0x0fc19dc6), UINT32_C(0x240ca1cc),
    UINT32_C(0x2de92c6f), UINT32_C(0x4a7484aa), UINT32_C(0x5cb0a9dc), UINT32_C(0x76f988da),
    UINT32_C(0x983e5152), UINT32_C(0xa831c66d), UINT32_C(0xb00327c8), UINT32_C(0xbf597fc7),
    UINT32_C(0xc6e00bf3), UINT32_C(0xd5a79147), UINT32_C(0x06ca6351), UINT32_C(0x14292967),
    UINT32_C(0x27b70a85), UINT32_C(0x2e1b2138), UINT32_C(0x4d2c6dfc), UINT32_C(0x53380d13),
    UINT32_C(0x650a7354), UINT32_C(0x766a0abb), UINT32_C(0x81c2c92e), UINT32_C(0x92722c85),
    UINT32_C(0xa2bfe8a1), UINT32_C(0xa81a664b), UINT32_C(0xc24b8b70), UINT32_C(0xc76c51a3),
    UINT32_C(0xd192e819), UINT32_C(0xd6990624), UINT32_C(0xf40e3585), UINT32_C(0x106aa070),
    UINT32_C(0x19a4c116), UINT32_C(0x1e376c08), UINT32_C(0x2748774c), UINT32_C(0x34b0bcb5),
    UINT32_C(0x391c0cb3), UINT32_C(0x4ed8aa4a), UINT32_C(0x5b9cca4f), UINT32_C(0x682e6ff3),
    UINT32_C(0x748f82ee), UINT32_C(0x78a5636f), UINT32_C(0x84c87814), UINT32_C(0x8cc70208),
    UINT32_C(0x90befffa), UINT32_C(0xa4506ceb), UINT32_C(0xbef9a3f7), UINT32_C(0xc67178f2)
};

static inline u32 rotr32(u32 x, unsigned n)
{
    // Used only with fixed n in [2, 25].
    return (x >> n) | (x << (32u - n));
}

static inline u32 load_be32(const u8 p[4])
{
    return ((u32)p[0] << 24) | ((u32)p[1] << 16) | ((u32)p[2] << 8)  |  (u32)p[3];
}

static inline u32 choose32(u32 x, u32 y, u32 z)
{
    return (x & y) ^ (~x & z);
}

static inline u32 majority32(u32 x, u32 y, u32 z)
{
    return (x & y) ^ (x & z) ^ (y & z);
}

static inline u32 big_sigma0(u32 x)
{
    return rotr32(x, 2) ^ rotr32(x, 13) ^ rotr32(x, 22);
}

static inline u32 big_sigma1(u32 x)
{
    return rotr32(x, 6) ^ rotr32(x, 11) ^ rotr32(x, 25);
}

static inline u32 small_sigma0(u32 x)
{
    return rotr32(x, 7) ^ rotr32(x, 18) ^ (x >> 3);
}

static inline u32 small_sigma1(u32 x)
{
    return rotr32(x, 17) ^ rotr32(x, 19) ^ (x >> 10);
}

void sha2_32_compress(u32 state[8], const u8 block[64], u32 schedule[16])
{
    u32 a = state[0], b = state[1], c = state[2], d = state[3];
    u32 e = state[4], f = state[5], g = state[6], h = state[7];

    /* Read the 512-bit block as sixteen big-endian 32-bit words. */
    for (unsigned t = 0; t < 16; ++t) { schedule[t] = load_be32(block + 4u * t); }

    /* Compute W[t] on demand using a 16-word ring buffer. */
    for (unsigned t = 0; t < 64; ++t) {
        u32 wt;
        if (t < 16) {
            wt = schedule[t];
        } else {
            const unsigned i = t & 15u;
            wt = schedule[i]
               + small_sigma0(schedule[(i + 1u) & 15u])
               + schedule[(i + 9u) & 15u]
               + small_sigma1(schedule[(i + 14u) & 15u]);
            schedule[i] = wt;
        }

        const u32 t1 = h + big_sigma1(e) + choose32(e, f, g)
                          + sha2_32_k[t] + wt;
        const u32 t2 = big_sigma0(a) + majority32(a, b, c);

        h = g;
        g = f;
        f = e;
        e = d + t1;
        d = c;
        c = b;
        b = a;
        a = t1 + t2;
    }

    /* Feed-forward: update the chaining state modulo 2^32. */
    state[0] += a;
    state[1] += b;
    state[2] += c;
    state[3] += d;
    state[4] += e;
    state[5] += f;
    state[6] += g;
    state[7] += h;
}
