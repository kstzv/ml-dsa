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

#define ML_DSA_44_BETA 78
#define ML_DSA_65_BETA 196
#define ML_DSA_87_BETA 120

#define ML_DSA_OMEGA_44 80
#define ML_DSA_OMEGA_65 55
#define ML_DSA_OMEGA_87 75

#define ML_DSA_44_GAMMA1          (1U << 17)
#define ML_DSA_65_87_GAMMA1       (1U << 19)

#define ML_DSA_44_GAMMA1_BITS         18
#define ML_DSA_65_87_GAMMA1_BITS      20

#define ML_DSA_44_POLYZ_PACKED_BYTES    576
#define ML_DSA_65_87_POLYZ_PACKED_BYTES 640

#define ML_DSA_44_GAMMA2    	 ((ML_DSA_Q - 1) / 88)
#define ML_DSA_65_87_GAMMA2 	 ((ML_DSA_Q - 1) / 32)
#define ML_DSA_Q_HALF_MINUS_ONE  ((ML_DSA_Q - 1) / 2) // 4190208

#define ML_DSA_44_SIZE_SIG 2420
#define ML_DSA_65_SIZE_SIG 3309
#define ML_DSA_87_SIZE_SIG 4627

#define ML_DSA_44_SIZE_PK 1312
#define ML_DSA_65_SIZE_PK 1952
#define ML_DSA_87_SIZE_PK 2592

#define ML_DSA_SIZE_SCRATCH_BUFFER 1088
#define ML_DSA_SIZE_MAX_CONTEXT    255
#define ML_DSA_PREHASH_OID_SIZE    11

// Error code for ML-DSA
#define ML_DSA_EINVAL                  22
#define ML_DSA_ENOMEM                  12
#define ML_DSA_EAGAIN          		   11
#define ML_DSA_SYSTEM_ENTROPY_FAILED   134
#define ML_DSA_CALLBACK_ENTROPY_FAILED 135
#define ML_DSA_EBADSIG                 136

#define MAX_NUMBER_SLOTS 128

#ifndef ML_DSA_FULL_MATRIX_BUFFER
#define ML_DSA_FULL_MATRIX_BUFFER 0
#endif

// for SIMD-frendly
#ifndef ML_DSA_SIMD
#error "ML_DSA_SIMD must be defined"
#endif

#define ML_DSA_SIZE_KECCAK_SIMD_1 25
#define ML_DSA_SIZE_KECCAK_SIMD_2 50
#define ML_DSA_SIZE_KECCAK_SIMD_4 100

#if ML_DSA_SIMD == 1
#define ML_DSA_SIZE_MEM_KECCAK ML_DSA_SIZE_KECCAK_SIMD_1
#include "shake.h"
#elif ML_DSA_SIMD == 2
#define ML_DSA_SIZE_MEM_KECCAK ML_DSA_SIZE_KECCAK_SIMD_2
#include "shake_x2.h"
#elif ML_DSA_SIMD == 4
#define ML_DSA_SIZE_MEM_KECCAK ML_DSA_SIZE_KECCAK_SIMD_4
#include "shake_x4.h"
#else 
#error "You need to specify 1, 2, or 4 for SIMD"
#endif

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

static const u8 oid_shake256[11] = {
    0x06, 0x09, 0x60, 0x86, 0x48, 0x01,
    0x65, 0x03, 0x04, 0x02, 0x0C
};

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

struct ml_dsa_workspace {
	s32 *matrix_buffer;
	s32 *temp_vector_buffer;
	s32 *vect_w;
	s32 *vect_y;
	s32 *poly_c;
	u8 deterministic;
	u8 rho_double_prime[ML_DSA_64_BYTES];
	u8 mu[ML_DSA_64_BYTES];
	u8 *scratch_buffer; // [ML_DSA_SIZE_SCRATCH_BUFFER];
	struct shake_ctx *shake;
	size_t len_messege;
};

// Struct for save keys parametrs
struct ml_dsa_keys {
	enum ml_dsa_level_k k;   // level k sequrity for key
	enum ml_dsa_level_l l;   // level l sequrity for key
	u8 eta;                  // Secret coefficient bound η
	u8 c_tilde_bytes;        // size of c, first part sign
	u8 rho[ML_DSA_32_BYTES]; // ρ - seed for Matrix
	u8 K[ML_DSA_32_BYTES];   // secret signing seed, for create ρ′′
	u8 tr[ML_DSA_64_BYTES];  // hash pk, for create μ
	u8 *pk;      			 // messege pk
	
	s32 *s1;                 
	s32 *s2;
	
	s32 *t0;
	s32 *t1;
	
	struct ml_dsa_workspace *workspace; // Data processing structure
};

struct ml_dsa_verify_workspace {
	ml_dsa_atomic_t is_free; // атомарні операції перевизначатимуться в портах
	struct ml_dsa_verify_workspace *next;
    enum ml_dsa_level_k k;
    enum ml_dsa_level_l l;

    s32 *poly_c;
    s32 *temp_vector_1;
    s32 *temp_vector_2;
 
    #if ML_DSA_FULL_MATRIX_BUFFER
    s32 *matrix_buffer;
    #else
    u8 rho[ML_DSA_32_BYTES];
    #endif
    
    u8 tr[ML_DSA_64_BYTES];
    u8 mu[ML_DSA_64_BYTES];
    u8 c_tilde[ML_DSA_64_BYTES];
    u8 *scratch_buffer;
    
    struct shake_ctx *shake;
};

struct ml_dsa_verify_ctx {
	struct ml_dsa_verify_workspace *work;
	size_t size_pool;
	bool fixed_pk;
};
	

// Internal functions from module create keys
extern struct ml_dsa_keys *ml_dsa_alloc_struct_keys(enum ml_dsa_level_k k, enum ml_dsa_level_l l);
extern void ml_dsa_destroy_struct_keys(struct ml_dsa_keys *ctx);
extern int get_rho_K_s1_s2(struct ml_dsa_keys *ctx, ml_dsa_entropy_fn entropy);
extern int create_t0_t1(struct ml_dsa_keys *ctx, s32 *t);
extern int ml_dsa_create_public_data(struct ml_dsa_keys *ctx);

// Internal function from matrix module
#if ML_DSA_FULL_MATRIX_BUFFER
extern int get_full_matrix(struct ml_dsa_keys *ctx);
extern int mult_matrix(s32 *matrix_buffer, s32 *vect, s32 *result, enum ml_dsa_level_k k, enum ml_dsa_level_l l);
#else
extern void mult_matrix(u8 *seed_matrix, s32 *vect, s32 *result, struct shake_ctx *shake, u8 *scratch_buffer, enum ml_dsa_level_k k, enum ml_dsa_level_l l);
#endif

// from ml_dsa_sign_verify_common.c module
extern void get_mu(struct shake_ctx *shake, u8 *mu, u8 *tr, const uint8_t *msg, size_t msg_len, const uint8_t *context, size_t context_len, u8 prehash);
extern void ml_dsa_get_w1(const s32 *w, s32 *w1, const u8 *h, enum ml_dsa_level_k k);
extern void ml_dsa_get_poly_c(s32 *poly_c, u8 *c, struct shake_ctx *shake, enum ml_dsa_level_k k);

// from ml_dsa_sign_core.c module
extern int get_rho_double_prime(struct ml_dsa_keys *ctx, u8 deterministic, ml_dsa_entropy_fn entropy);
extern void ml_dsa_expand_mask(struct ml_dsa_keys *ctx, u32 kappa);
extern bool ml_dsa_get_z(struct ml_dsa_keys *ctx, u8 *out);
extern bool ml_dsa_get_r0(struct ml_dsa_keys *ctx);
extern bool ml_dsa_get_h(struct ml_dsa_keys *ctx, u8 *out);
extern void ml_dsa_get_c(struct ml_dsa_keys *ctx, u8 *c);

// from ml_dsa_verify_core.c module
extern int ml_dsa_unpack_c_z(struct ml_dsa_verify_workspace *ctx, const u8 *sig_msg);
extern void ml_dsa_unpack_pk(s32 * t1, const u8 *pk_msg, enum ml_dsa_level_k k);
extern int ml_dsa_unpack_h(struct ml_dsa_verify_workspace *ctx, u8 *sign);
extern int ml_dsa_pack_w1(struct ml_dsa_verify_workspace *ctx);
extern int ml_dsa_get_c_prime(struct ml_dsa_verify_workspace *ctx);

// ---------------------------Functions of Barrett Reductions--------------------------------------------------------------------------------
// Special constants for Barrett reduction in ML-DSA.
// The ML-DSA standard (FIPS 204) defines the modulus q = 8380417,
// but does not mandate a specific reduction algorithm
// (works with secrets or data dependent on a secret)
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

// Montgomery reduction for multiplication (works with secrets or data dependent on a secret)
static inline s32 ml_dsa_montgomery_reduce(s64 val)
{
	s32 t;

	t = (s32)((s64)(s32)val * INT64_C(58728449));
	return (s32)((val - (s64)t * INT64_C(8380417)) >> 32);
}

// NTT - for representing polynomials in point format (works with secrets or data dependent on a secret)
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
 
// iNTT - for representing polynomials in coefficient format (works with secrets or data dependent on a secret)
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

// (works with secrets or data dependent on a secret)
static inline void ml_dsa_canonicalize(s32 w[ML_DSA_N])
{
    for (size_t i = 0; i < ML_DSA_N; i++)
    {
        w[i] += (w[i] >> 31) & ML_DSA_Q;
    }
}

// Algorith 36 (works with secrets or data dependent on a secret)
inline void ml_dsa_decompose(s32 *r, s32 gamma2, s32 *r1, s32 *r0)
{
    s32 high = (*r + 127) >> 7;

    if (gamma2 == ML_DSA_65_87_GAMMA2) 
    {
        high = (high * 1025 + (1 << 21)) >> 22;
        high &= 15;
    } else 
    {
        high = (high * 11275 + (1 << 23)) >> 24;
        high ^= ((43 - high) >> 31) & high;
    }

    s32 low = *r - high * 2 * gamma2;
    low -= ((ML_DSA_Q_HALF_MINUS_ONE - low) >> 31) & ML_DSA_Q;

    if (r1 != NULL) { *r1 = high; }
	if (r0 != NULL) { *r0 = low; }
}				 

// Приведення до модуля q
static inline s32 ml_dsa_freeze(s32 val)
{
    val = ml_dsa_barrett_reduce(val);
    val += (val >> 31) & ML_DSA_Q;
    return val;
}				 

#endif // ML_DSA_CORE_HEADER_H
