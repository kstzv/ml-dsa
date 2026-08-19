// ML_DSA_CORE_HEADER_H
#ifndef ML_DSA_CORE_HEADER_H
#define ML_DSA_CORE_HEADER_H

// The header is under development, parts of the code will be written, edited, changed and deleted

// Algorithm constants
#define ML_DSA_32_BYTES 32
#define ML_DSA_64_BYTES 64

#define ML_DSA_44_87_ETA 2
#define ML_DSA_65_ETA    4

#define ML_DSA_Q 8380417
#define ML_DSA_N 256

#define ML_DSA_D                  13
#define ML_DSA_POWER2ROUND_BASE   (1 << ML_DSA_D)
#define ML_DSA_POWER2ROUND_HALF   (1 << (ML_DSA_D - 1))

// Error code for ML-DSA
#define ML_DSA_EINVAL                  22
#define ML_DSA_ENOMEM                  12
#define ML_DSA_SYSTEM_ENTROPY_FAILED   134
#define ML_DSA_CALLBACK_ENTROPY_FAILED 135

// HEADERS
#if defined(USERSPACE)

#include <ml_dsa_userspace_port.h>

#elif defined(LINUX_KERNEL) // TODO

#include "ml_dsa_linux_kernel_port.h"

#elif defined(FREERTOS) //TODO

#include "ml_dsa_freertos_port.h"

#else
#error "No ML-DSA port selected."
#endif // HEADERS

// Callback used to obtain cryptographically secure random bytes.
// Parameters:
//   buf - destination buffer
//   len - number of bytes to generate
// Returns:
//   0 on success
//   non-zero on failure
typedef int (*ml_dsa_entropy_fn)(void *buf, size_t len);

extern void *ml_dsa_alloc(size_t size_alloc);
extern void ml_dsa_free(void *ptr);

extern void ml_dsa_memzero(void *ptr, size_t len);
extern int ml_dsa_entropy(void *buf, size_t len);


// Paramentr for level sequrity
enum ml_dsa_level_k {
	ML_DSA_44_K = 4,
	ML_DSA_65_K = 6,
	ML_DSA_87_K = 8,
};

enum ml_dsa_level_l {
	ML_DSA_44_L = 4,
	ML_DSA_65_L = 5,
	ML_DSA_87_L = 7,
};

// Struct for save keys parametrs
struct ml_dsa_keys {
	enum ml_dsa_level_k k;   // level k sequrity for key
	enum ml_dsa_level_l l;   // level l sequrity for key
	u8 rho[ML_DSA_32_BYTES]; // ρ - seed for Matrix
	u8 K[ML_DSA_32_BYTES];   // secret signing seed, for create ρ′′
	u8 tr[ML_DSA_64_BYTES];  // hash pk, for create μ
	
	s8 *s1;
	s8 *s2;
	
	s16 *t0;
	s16 *t1;
};

// ---------------------------Functions of Barrett Reductions--------------------------------------------------------------------------------
// Special constants for Barrett reduction in ML-DSA.
// The ML-DSA standard (FIPS 204) defines the modulus q = 8380417,
// but does not mandate a specific reduction algorithm

#define ML_DSA_BARRETT_ROUND 4194304
// #define ML_DSA_QINV    58728449
#define CONST_K   23
static inline s32 ml_dsa_barrett_reduce(s32 val)
{
	s32 temp = val + ML_DSA_BARRETT_ROUND;
	temp >>= CONST_K;
	temp *= ML_DSA_Q;
	return val - temp;
}

// Special table for twiddle-factor in NTT
static const s32 ml_dsa_zetas_mont[ML_DSA_N] = {
    0,    25847, -2608894, -518909,   237124, -777960, -876248,   466468,
    1826347,  2353451, -359251, -2091905,  3119733, -2884855,  3111497,  2680103,
    2725464,  1024112, -1079900,  3585928, -549488, -1119584,  2619752, -2108549,
    -2118186, -3859737, -1399561, -3277672,  1757237, -19422,  4010497,   280005,
    2706023,    95776,  3077325,  3530437, -1661693, -3592148, -2537516,  3915439,
    -3861115, -3043716,  3574422, -2867647,  3539968, -300467,  2348700, -539299,
    -1699267, -1643818,  3505694, -3821735,  3507263, -2140649, -1600420,  3699596,
    811944,   531354,   954230,  3881043,  3900724, -2556880,  2071892, -2797779,
    -3930395, -1528703, -3677745, -3041255, -1452451,  3475950,  2176455, -1585221,
    -1257611,  1939314, -4083598, -1000202, -3190144, -3157330, -3632928,   126922,
    3412210, -983419,  2147896,  2715295, -2967645, -3693493, -411027, -2477047,
    -671102, -1228525, -22981, -1308169, -381987,  1349076,  1852771, -1430430,
    -3343383,   264944,   508951,  3097992,    44288, -1100098,   904516,  3958618,
    -3724342, -8578,  1653064, -3249728,  2389356, -210977,   759969, -1316856,
    189548, -3553272,  3159746, -1851402, -2409325, -177440,  1315589,  1341330,
    1285669, -1584928, -812732, -1439742, -3019102, -3881060, -3628969,  3839961,
    2091667,  3407706,  2316500,  3817976, -3342478,  2244091, -2446433, -3562462,
    266997,  2434439, -1235728,  3513181, -3520352, -3759364, -1197226, -3193378,
    900702,  1859098,   909542,   819034,   495491, -1613174, -43260, -522500,
    -655327, -3122442,  2031748,  3207046, -3556995, -525098, -768622, -3595838,
    342297,   286988, -2437823,  4108315,  3437287, -3342277,  1735879,   203044,
    2842341,  2691481, -2590150,  1265009,  4055324,  1247620,  2486353,  1595974,
    -3767016,  1250494,  2635921, -3548272, -2994039,  1869119,  1903435, -1050970,
    -1333058,  1237275, -3318210, -1430225, -451100,  1312455,  3306115, -1962642,
    -1279661,  1917081, -2546312, -1374803,  1500165,   777191,  2235880,  3406031,
    -542412, -2831860, -1671176, -1846953, -2584293, -3724270,   594136, -3776993,
    -2013608,  2432395,  2454455, -164721,  1957272,  3369112,   185531, -1207385,
    -3183426,   162844,  1616392,  3014001,   810149,  1652634, -3694233, -1799107,
    -3038916,  3523897,  3866901,   269760,  2213111, -975884,  1717735,   472078,
    -426683,  1723600, -1803090,  1910376, -1667432, -1104333, -260646, -3833893,
    -2939036, -2235985, -420899, -2286327,   183443, -976891,  1612842, -3545687,
    -554416,  3919660, -48306, -1362209,  3937738,  1400424, -846154,  1976782
};

// Montgomery reduction for multiplication
static inline s32 ml_dsa_montgomery_reduce(s64 val)
{
	s32 t;

	t = (s32)((s64)(s32)val * INT64_C(58728449));
	return (s32)((val - (s64)t * INT64_C(8380417)) >> 32);
}

// NTT - for representing polynomials in point format
static inline void ml_dsa_ntt(s32 w[ML_DSA_N])
{
    size_t m = 0;
    s32 t;
    
    for (size_t len = 128; len > 0; len >>= 1)
    {
		for (size_t start = 0; start < ML_DSA_N; start += 2 * len)
		{
			 m++;
			 s32 z = ml_dsa_zetas_mont[m];
			 for (size_t j = start; j < start + len; ++j)
			 {
				 t = ml_dsa_montgomery_reduce((s64)z * (s64)w[j + len]);
				 w[j + len] = w[j] - t;
				 w[j] = w[j] + t;
			 }
		 }
	 }
}
 
// iNTT - for representing polynomials in coefficient format
static inline void ml_dsa_intt(s32 w[ML_DSA_N])
{
    size_t m = 256;
    s32 t;
    
    for (size_t len = 1; len < 256; len <<= 1)
    {
		for (size_t start = 0; start < ML_DSA_N; start += 2 * len)
		{
			 m--;
			 s32 z = -ml_dsa_zetas_mont[m];
			 for (size_t j = start; j < start + len; ++j)
			 {
				 t = w[j];
				 w[j] = t + w[j + len];
				 w[j + len] = t - w[j + len];
				 w[j + len] = ml_dsa_montgomery_reduce((s64)z * (s64)w[j + len]);
			 }
		 }
	 }
	 
	 for(size_t i = 0; i < ML_DSA_N; i++)
	 {
		 w[i] = ml_dsa_montgomery_reduce(INT64_C(41978) *  (s64)w[i]);
	 }
 }
				 
				 

#endif // ML_DSA_CORE_HEADER_H
