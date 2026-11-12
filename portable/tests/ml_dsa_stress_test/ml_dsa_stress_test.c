/* ML-DSA userspace end-to-end stress test (C11).
 * Tests keygen -> sign -> verify, negative checks, and fixed/dynamic verification pools.
 */
#include "ml_dsa.h"

#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAX_MESSAGE_LEN 4096u
#define VERIFY_SLOTS 2u

static u64 rng_state = UINT64_C(0x4d4c445341535452);

/* Non-cryptographic PRNG: only test messages/context use it.
 * Key generation and hedged signing use the userspace entropy port. */
static u64 next_random(void)
{
    u64 x = rng_state;
    x ^= x << 13;
    x ^= x >> 7;
    x ^= x << 17;
    return rng_state = x;
}

static void fill_bytes(u8 *dst, size_t len)
{
    for (size_t i = 0; i < len; ++i) {
        dst[i] = (u8)next_random();
    }
}

static int checked_verify(const char *which,
                          struct ml_dsa_verify_ctx *pool,
                          const u8 *pk, size_t pk_len,
                          const u8 *msg, size_t msg_len,
                          const u8 *context, size_t context_len,
                          const u8 *sig, size_t sig_len,
                          u8 prehash, int should_accept,
                          u8 level, size_t iter)
{
    const int rc = ml_dsa_verify(pool, pk, pk_len, msg, msg_len,
                                 context, context_len, sig, sig_len, prehash);
    /* A tampered message/signature/context must fail as a bad signature.
     * EINVAL/EAGAIN are not considered a successful negative test. */
    const int expected = should_accept ? 0 : -ML_DSA_EBADSIG;
    if (rc != expected) {
        fprintf(stderr,
                "[FAIL] ML-DSA-%u iter=%zu check=%s: got %d, expected %d\n",
                (unsigned)level, iter, which, rc, expected);
        return -1;
    }
    return 0;
}

static int run_level(u8 level, size_t iterations)
{
    static const size_t msg_lengths[] = {
        0, 1, 31, 32, 63, 64, 135, 136, 137, 255,
        256, 407, 408, 409, 1024, 2048, MAX_MESSAGE_LEN
    };
    static const size_t context_lengths[] = {0, 1, 32, 255};

    const size_t pk_len = level == 44 ? ML_DSA_44_SIZE_PK :
                          level == 65 ? ML_DSA_65_SIZE_PK : ML_DSA_87_SIZE_PK;
    const size_t sig_len = level == 44 ? ML_DSA_44_SIZE_SIG :
                           level == 65 ? ML_DSA_65_SIZE_SIG : ML_DSA_87_SIZE_SIG;

    u8 pk[ML_DSA_87_SIZE_PK];
    u8 sig[ML_DSA_87_SIZE_SIG];
    u8 sig_again[ML_DSA_87_SIZE_SIG];
    u8 msg[MAX_MESSAGE_LEN];
    u8 context[ML_DSA_SIZE_MAX_CONTEXT];
    struct ml_dsa_keys *keys = NULL;
    struct ml_dsa_verify_ctx *dynamic_pool = NULL;
    struct ml_dsa_verify_ctx *fixed_pool = NULL;
    int result = -1;

    printf("[START] ML-DSA-%u, iterations=%zu, verify slots=%u\n",
           (unsigned)level, iterations, VERIFY_SLOTS);
    fflush(stdout);

    keys = ml_dsa_create_keys(level, NULL);
    if (!keys) {
        fprintf(stderr, "[FAIL] ML-DSA-%u: key generation returned NULL\n", level);
        goto cleanup;
    }
    const int pk_result = ml_dsa_get_pk(keys, pk, pk_len);
    if (pk_result != (int)pk_len) {
        fprintf(stderr, "[FAIL] ML-DSA-%u: get_pk returned %d\n", level, pk_result);
        goto cleanup;
    }

    dynamic_pool = ml_dsa_alloc_verify_pool(level, VERIFY_SLOTS, NULL, 0);
    if (!dynamic_pool) {
        fprintf(stderr, "[FAIL] ML-DSA-%u: dynamic verify pool allocation\n", level);
        goto cleanup;
    }
    fixed_pool = ml_dsa_alloc_verify_pool(level, VERIFY_SLOTS, pk, pk_len);
    if (!fixed_pool) {
        fprintf(stderr, "[FAIL] ML-DSA-%u: fixed verify pool allocation\n", level);
        goto cleanup;
    }

    for (size_t i = 0; i < iterations; ++i) {
        const size_t msg_len = msg_lengths[i % (sizeof(msg_lengths) / sizeof(msg_lengths[0]))];
        const size_t ctx_len = context_lengths[(i / 3) % (sizeof(context_lengths) / sizeof(context_lengths[0]))];
        const u8 *ctx_ptr = ctx_len ? context : NULL;
        const u8 prehash = (u8)((i >> 1) & 1u);
        const u8 deterministic = (u8)(i & 1u);
        fill_bytes(msg, msg_len ? msg_len : 1u);
        fill_bytes(context, ctx_len);

        const int sign_rc = ml_dsa_sign(keys, msg, msg_len, ctx_ptr, ctx_len,
                                        sig, sig_len, prehash, deterministic, NULL);
        if (sign_rc != 0) {
            fprintf(stderr,
                    "[FAIL] ML-DSA-%u iter=%zu sign=%d msg_len=%zu ctx_len=%zu prehash=%u deterministic=%u\n",
                    level, i, sign_rc, msg_len, ctx_len, prehash, deterministic);
            goto cleanup;
        }

#define VERIFY(name, pool, m, mlen, c, clen, s, accept)                       \
        do {                                                                  \
            const u8 *verify_pk = ((pool) == fixed_pool) ? NULL : pk;        \
            size_t verify_pk_len = ((pool) == fixed_pool) ? 0u : pk_len;    \
            if (checked_verify((name), (pool), verify_pk, verify_pk_len,    \
                               (m), (mlen), (c), (clen), (s), sig_len,      \
                               prehash, (accept), level, i) != 0) {        \
                goto cleanup;                                                \
            }                                                                 \
        } while (0)

        VERIFY("dynamic pool / valid", dynamic_pool, msg, msg_len, ctx_ptr, ctx_len, sig, 1);
        VERIFY("fixed pool / valid", fixed_pool, msg, msg_len, ctx_ptr, ctx_len, sig, 1);

        /* Slow path without a pool, checked occasionally. */
        if (i % 17u == 0) {
            VERIFY("no pool / valid", NULL, msg, msg_len, ctx_ptr, ctx_len, sig, 1);
        }

        /* A changed message must invalidate an otherwise valid signature.
         * For an empty message use a different (one-byte) message. */
        if (msg_len == 0) {
            msg[0] ^= 1u;
            VERIFY("changed message", fixed_pool, msg, 1u, ctx_ptr, ctx_len, sig, 0);
            msg[0] ^= 1u;
        } else {
            msg[msg_len / 2u] ^= 1u;
            VERIFY("changed message", fixed_pool, msg, msg_len, ctx_ptr, ctx_len, sig, 0);
            msg[msg_len / 2u] ^= 1u;
        }

        sig[0] ^= 1u;
        VERIFY("changed signature", dynamic_pool, msg, msg_len, ctx_ptr, ctx_len, sig, 0);
        sig[0] ^= 1u;

        if (ctx_len == 0) {
            context[0] = 0x53;
            VERIFY("changed context", fixed_pool, msg, msg_len, context, 1u, sig, 0);
        } else {
            context[0] ^= 1u;
            VERIFY("changed context", fixed_pool, msg, msg_len, context, ctx_len, sig, 0);
            context[0] ^= 1u;
        }

        /* Signing must be repeatable in deterministic mode for identical input. */
        if (deterministic && i % 31u == 1u) {
            const int again_rc = ml_dsa_sign(keys, msg, msg_len, ctx_ptr, ctx_len,
                                               sig_again, sig_len, prehash, 1, NULL);
            if (again_rc != 0 || memcmp(sig, sig_again, sig_len) != 0) {
                fprintf(stderr, "[FAIL] ML-DSA-%u iter=%zu deterministic repeat (rc=%d)\n",
                        level, i, again_rc);
                goto cleanup;
            }
        }
#undef VERIFY
        if ((i + 1) % 100u == 0 || i + 1 == iterations) {
            printf("[PROGRESS] ML-DSA-%u %zu/%zu\n", level, i + 1, iterations);
            fflush(stdout);
        }
    }
    printf("[PASS] ML-DSA-%u (%zu iterations)\n", level, iterations);
    result = 0;

cleanup:
    /* Fixed pool owns the t1 allocation; dynamic pool does not. */
    destructor_verify_pool(fixed_pool);
    destructor_verify_pool(dynamic_pool);
    ml_dsa_destroy_struct_keys(keys);
    return result;
}

static int parse_count(const char *s, size_t *out)
{
    char *end = NULL;
    errno = 0;
    unsigned long long number = strtoull(s, &end, 10);
    if (errno != 0 || end == s || *end != '\0' || !number || number > SIZE_MAX) {
        return -1;
    }
    *out = (size_t)number;
    return 0;
}

int main(int argc, char **argv)
{
    size_t iterations = 100;
    const char *level = argc > 1 ? argv[1] : "all";
    if (argc > 3 || (argc > 2 && parse_count(argv[2], &iterations) != 0)) {
        fprintf(stderr, "Usage: %s [44|65|87|all] [positive_iterations]\n", argv[0]);
        return EXIT_FAILURE;
    }
    if (strcmp(level, "all") == 0) {
        const u8 levels[] = {44, 65, 87};
        for (size_t i = 0; i < sizeof(levels) / sizeof(levels[0]); ++i) {
            if (run_level(levels[i], iterations) != 0) return EXIT_FAILURE;
        }
    } else if (strcmp(level, "44") == 0 || strcmp(level, "65") == 0 || strcmp(level, "87") == 0) {
        if (run_level((u8)atoi(level), iterations) != 0) return EXIT_FAILURE;
    } else {
        fprintf(stderr, "Unknown ML-DSA level: %s (expected 44, 65, 87 or all)\n", level);
        return EXIT_FAILURE;
    }
    puts("[PASS] All requested ML-DSA stress checks passed.");
    return EXIT_SUCCESS;
}
