// SPDX-License-Identifier: GPL-2.0 OR MIT
#include "sha2_64.h"

static const u64 sha2_64_k[80] = {
    UINT64_C(0x428a2f98d728ae22), UINT64_C(0x7137449123ef65cd),
    UINT64_C(0xb5c0fbcfec4d3b2f), UINT64_C(0xe9b5dba58189dbbc),
    UINT64_C(0x3956c25bf348b538), UINT64_C(0x59f111f1b605d019),
    UINT64_C(0x923f82a4af194f9b), UINT64_C(0xab1c5ed5da6d8118),
    UINT64_C(0xd807aa98a3030242), UINT64_C(0x12835b0145706fbe),
    UINT64_C(0x243185be4ee4b28c), UINT64_C(0x550c7dc3d5ffb4e2),
    UINT64_C(0x72be5d74f27b896f), UINT64_C(0x80deb1fe3b1696b1),
    UINT64_C(0x9bdc06a725c71235), UINT64_C(0xc19bf174cf692694),
    UINT64_C(0xe49b69c19ef14ad2), UINT64_C(0xefbe4786384f25e3),
    UINT64_C(0x0fc19dc68b8cd5b5), UINT64_C(0x240ca1cc77ac9c65),
    UINT64_C(0x2de92c6f592b0275), UINT64_C(0x4a7484aa6ea6e483),
    UINT64_C(0x5cb0a9dcbd41fbd4), UINT64_C(0x76f988da831153b5),
    UINT64_C(0x983e5152ee66dfab), UINT64_C(0xa831c66d2db43210),
    UINT64_C(0xb00327c898fb213f), UINT64_C(0xbf597fc7beef0ee4),
    UINT64_C(0xc6e00bf33da88fc2), UINT64_C(0xd5a79147930aa725),
    UINT64_C(0x06ca6351e003826f), UINT64_C(0x142929670a0e6e70),
    UINT64_C(0x27b70a8546d22ffc), UINT64_C(0x2e1b21385c26c926),
    UINT64_C(0x4d2c6dfc5ac42aed), UINT64_C(0x53380d139d95b3df),
    UINT64_C(0x650a73548baf63de), UINT64_C(0x766a0abb3c77b2a8),
    UINT64_C(0x81c2c92e47edaee6), UINT64_C(0x92722c851482353b),
    UINT64_C(0xa2bfe8a14cf10364), UINT64_C(0xa81a664bbc423001),
    UINT64_C(0xc24b8b70d0f89791), UINT64_C(0xc76c51a30654be30),
    UINT64_C(0xd192e819d6ef5218), UINT64_C(0xd69906245565a910),
    UINT64_C(0xf40e35855771202a), UINT64_C(0x106aa07032bbd1b8),
    UINT64_C(0x19a4c116b8d2d0c8), UINT64_C(0x1e376c085141ab53),
    UINT64_C(0x2748774cdf8eeb99), UINT64_C(0x34b0bcb5e19b48a8),
    UINT64_C(0x391c0cb3c5c95a63), UINT64_C(0x4ed8aa4ae3418acb),
    UINT64_C(0x5b9cca4f7763e373), UINT64_C(0x682e6ff3d6b2b8a3),
    UINT64_C(0x748f82ee5defb2fc), UINT64_C(0x78a5636f43172f60),
    UINT64_C(0x84c87814a1f0ab72), UINT64_C(0x8cc702081a6439ec),
    UINT64_C(0x90befffa23631e28), UINT64_C(0xa4506cebde82bde9),
    UINT64_C(0xbef9a3f7b2c67915), UINT64_C(0xc67178f2e372532b),
    UINT64_C(0xca273eceea26619c), UINT64_C(0xd186b8c721c0c207),
    UINT64_C(0xeada7dd6cde0eb1e), UINT64_C(0xf57d4f7fee6ed178),
    UINT64_C(0x06f067aa72176fba), UINT64_C(0x0a637dc5a2c898a6),
    UINT64_C(0x113f9804bef90dae), UINT64_C(0x1b710b35131c471b),
    UINT64_C(0x28db77f523047d84), UINT64_C(0x32caab7b40c72493),
    UINT64_C(0x3c9ebe0a15c9bebc), UINT64_C(0x431d67c49c100d4c),
    UINT64_C(0x4cc5d4becb3e42b6), UINT64_C(0x597f299cfc657e2a),
    UINT64_C(0x5fcb6fab3ad6faec), UINT64_C(0x6c44198c4a475817)
};

static inline u64 rotr64(u64 x, unsigned n)
{
    // Called only with fixed n in [1, 61].
    return (x >> n) | (x << (64u - n));
}

static inline u64 load_be64(const u8 p[8])
{
    return ((u64)p[0] << 56) | ((u64)p[1] << 48) |
           ((u64)p[2] << 40) | ((u64)p[3] << 32) |
           ((u64)p[4] << 24) | ((u64)p[5] << 16) |
           ((u64)p[6] << 8)  | (u64)p[7];
}

static inline u64 choose64(u64 x, u64 y, u64 z)
{
    return (x & y) ^ (~x & z);
}

static inline u64 majority64(u64 x, u64 y, u64 z)
{
    return (x & y) ^ (x & z) ^ (y & z);
}

static inline u64 big_sigma0(u64 x)
{
    return rotr64(x, 28) ^ rotr64(x, 34) ^ rotr64(x, 39);
}

static inline u64 big_sigma1(u64 x)
{
    return rotr64(x, 14) ^ rotr64(x, 18) ^ rotr64(x, 41);
}

static inline u64 small_sigma0(u64 x)
{
    return rotr64(x, 1) ^ rotr64(x, 8) ^ (x >> 7);
}

static inline u64 small_sigma1(u64 x)
{
    return rotr64(x, 19) ^ rotr64(x, 61) ^ (x >> 6);
}

void sha2_64_compress(u64 state[8], const u8 block[SHA2_64_BLOCK_SIZE], u64 schedule[16])
{
    u64 a = state[0], b = state[1], c = state[2], d = state[3];
    u64 e = state[4], f = state[5], g = state[6], h = state[7];

    /* Read the 1024-bit block as sixteen big-endian 64-bit words. */
    for (unsigned t = 0; t < 16; ++t) { schedule[t] = load_be64(block + 8u * t); }

    /* Compute W[t] on demand using a sixteen-word ring buffer. */
    for (unsigned t = 0; t < 80; ++t)
    {
        u64 wt;
        if (t < 16)
        {
            wt = schedule[t];
        }
        else
        {
            const unsigned i = t & 15u;
            wt = schedule[i]
               + small_sigma0(schedule[(i + 1u) & 15u])
               + schedule[(i + 9u) & 15u]
               + small_sigma1(schedule[(i + 14u) & 15u]);
            schedule[i] = wt;
        }

        const u64 t1 = h + big_sigma1(e) + choose64(e, f, g)
                       + sha2_64_k[t] + wt;
        const u64 t2 = big_sigma0(a) + majority64(a, b, c);

        h = g;
        g = f;
        f = e;
        e = d + t1;
        d = c;
        c = b;
        b = a;
        a = t1 + t2;
    }

    /* Feed-forward: update the chaining state modulo 2^64. */
    state[0] += a;
    state[1] += b;
    state[2] += c;
    state[3] += d;
    state[4] += e;
    state[5] += f;
    state[6] += g;
    state[7] += h;
}
