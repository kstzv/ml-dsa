#include "shake.h"

static const uint64_t keccak_rc[24] = {
    0x0000000000000001ULL,
    0x0000000000008082ULL,
    0x800000000000808AULL,
    0x8000000080008000ULL,
    0x000000000000808BULL,
    0x0000000080000001ULL,
    0x8000000080008081ULL,
    0x8000000000008009ULL,
    0x000000000000008AULL,
    0x0000000000000088ULL,
    0x0000000080008009ULL,
    0x000000008000000AULL,
    0x000000008000808BULL,
    0x800000000000008BULL,
    0x8000000000008089ULL,
    0x8000000000008003ULL,
    0x8000000000008002ULL,
    0x8000000000000080ULL,
    0x000000000000800AULL,
    0x800000008000000AULL,
    0x8000000080008081ULL,
    0x8000000000008080ULL,
    0x0000000080000001ULL,
    0x8000000080008008ULL
};
	
static inline uint64_t rotl64(uint64_t x, unsigned n)
{
    return (x << n) | (x >> (64 - n));
}

static inline void keccak_theta(uint64_t s[25]);
static inline void keccak_rho(uint64_t s[25]);
static inline void keccak_pi(uint64_t s[25]);
static inline void keccak_chi(uint64_t s[25]);
static inline void keccak_iota(uint64_t s[25], unsigned round);

void keccak_f1600(uint64_t s[25])
{
    for (unsigned round = 0; round < 24; round++) 
    {
        keccak_theta(s);
        keccak_rho(s);
        keccak_pi(s);
        keccak_chi(s);
        keccak_iota(s, round);
    }
}

static inline void keccak_theta(uint64_t s[25])
{
    uint64_t c[5];
    uint64_t d[5];

    c[0] = s[0] ^ s[5] ^ s[10] ^ s[15] ^ s[20];
    c[1] = s[1] ^ s[6] ^ s[11] ^ s[16] ^ s[21];
    c[2] = s[2] ^ s[7] ^ s[12] ^ s[17] ^ s[22];
    c[3] = s[3] ^ s[8] ^ s[13] ^ s[18] ^ s[23];
    c[4] = s[4] ^ s[9] ^ s[14] ^ s[19] ^ s[24];

    d[0] = c[4] ^ rotl64(c[1], 1);
    d[1] = c[0] ^ rotl64(c[2], 1);
    d[2] = c[1] ^ rotl64(c[3], 1);
    d[3] = c[2] ^ rotl64(c[4], 1);
    d[4] = c[3] ^ rotl64(c[0], 1);

    s[0]  ^= d[0];
    s[5]  ^= d[0];
    s[10] ^= d[0];
    s[15] ^= d[0];
    s[20] ^= d[0];

    s[1]  ^= d[1];
    s[6]  ^= d[1];
    s[11] ^= d[1];
    s[16] ^= d[1];
    s[21] ^= d[1];

    s[2]  ^= d[2];
    s[7]  ^= d[2];
    s[12] ^= d[2];
    s[17] ^= d[2];
    s[22] ^= d[2];

    s[3]  ^= d[3];
    s[8]  ^= d[3];
    s[13] ^= d[3];
    s[18] ^= d[3];
    s[23] ^= d[3];

    s[4]  ^= d[4];
    s[9]  ^= d[4];
    s[14] ^= d[4];
    s[19] ^= d[4];
    s[24] ^= d[4];
}

static inline void keccak_rho(uint64_t s[25])
{
    // s[0] rotate 0 -- nothing to do
    s[1] = rotl64(s[1], 1);
    s[2] = rotl64(s[2], 62);
    s[3] = rotl64(s[3], 28);
    s[4] = rotl64(s[4], 27);

    s[5] = rotl64(s[5], 36);
    s[6] = rotl64(s[6], 44);
    s[7] = rotl64(s[7], 6);
    s[8] = rotl64(s[8], 55);
    s[9] = rotl64(s[9], 20);

    s[10] = rotl64(s[10], 3);
    s[11] = rotl64(s[11], 10);
    s[12] = rotl64(s[12], 43);
    s[13] = rotl64(s[13], 25);
    s[14] = rotl64(s[14], 39);

    s[15] = rotl64(s[15], 41);
    s[16] = rotl64(s[16], 45);
    s[17] = rotl64(s[17], 15);
    s[18] = rotl64(s[18], 21);
    s[19] = rotl64(s[19], 8);

    s[20] = rotl64(s[20], 18);
    s[21] = rotl64(s[21], 2);
    s[22] = rotl64(s[22], 61);
    s[23] = rotl64(s[23], 56);
    s[24] = rotl64(s[24], 14);
}

static inline void keccak_pi(uint64_t s[25])
{
    uint64_t t = s[1];

    s[1]  = s[6];
    s[6]  = s[9];
    s[9]  = s[22];
    s[22] = s[14];
    s[14] = s[20];
    s[20] = s[2];
    s[2]  = s[12];
    s[12] = s[13];
    s[13] = s[19];
    s[19] = s[23];
    s[23] = s[15];
    s[15] = s[4];
    s[4]  = s[24];
    s[24] = s[21];
    s[21] = s[8];
    s[8]  = s[16];
    s[16] = s[5];
    s[5]  = s[3];
    s[3]  = s[18];
    s[18] = s[17];
    s[17] = s[11];
    s[11] = s[7];
    s[7]  = s[10];
    s[10] = t;

    // s[0] remains in place
}

static inline void keccak_chi(uint64_t s[25])
{
    uint64_t a0, a1, a2, a3, a4;


    a0 = s[0];
    a1 = s[1];
    a2 = s[2];
    a3 = s[3];
    a4 = s[4];

    s[0] = a0 ^ ((~a1) & a2);
    s[1] = a1 ^ ((~a2) & a3);
    s[2] = a2 ^ ((~a3) & a4);
    s[3] = a3 ^ ((~a4) & a0);
    s[4] = a4 ^ ((~a0) & a1);

    
    a0 = s[5];
    a1 = s[6];
    a2 = s[7];
    a3 = s[8];
    a4 = s[9];

    s[5] = a0 ^ ((~a1) & a2);
    s[6] = a1 ^ ((~a2) & a3);
    s[7] = a2 ^ ((~a3) & a4);
    s[8] = a3 ^ ((~a4) & a0);
    s[9] = a4 ^ ((~a0) & a1);

    
    a0 = s[10];
    a1 = s[11];
    a2 = s[12];
    a3 = s[13];
    a4 = s[14];

    s[10] = a0 ^ ((~a1) & a2);
    s[11] = a1 ^ ((~a2) & a3);
    s[12] = a2 ^ ((~a3) & a4);
    s[13] = a3 ^ ((~a4) & a0);
    s[14] = a4 ^ ((~a0) & a1);

    
    a0 = s[15];
    a1 = s[16];
    a2 = s[17];
    a3 = s[18];
    a4 = s[19];

    s[15] = a0 ^ ((~a1) & a2);
    s[16] = a1 ^ ((~a2) & a3);
    s[17] = a2 ^ ((~a3) & a4);
    s[18] = a3 ^ ((~a4) & a0);
    s[19] = a4 ^ ((~a0) & a1);

    
    a0 = s[20];
    a1 = s[21];
    a2 = s[22];
    a3 = s[23];
    a4 = s[24];

    s[20] = a0 ^ ((~a1) & a2);
    s[21] = a1 ^ ((~a2) & a3);
    s[22] = a2 ^ ((~a3) & a4);
    s[23] = a3 ^ ((~a4) & a0);
    s[24] = a4 ^ ((~a0) & a1);
}

static inline void keccak_iota(uint64_t s[25], unsigned round)
{
    s[0] ^= keccak_rc[round];
}		
