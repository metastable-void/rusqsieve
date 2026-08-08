#ifndef RUSQSIEVE_H
#define RUSQSIEVE_H

#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/*
 * Completely opaque Rust-owned factorization result. Factors are exposed as
 * borrowed NUL-terminated decimal strings, sorted ascending and repeated
 * according to multiplicity.
 */
typedef struct rusqsieve_factors rusqsieve_factors;

enum rusqsieve_status {
    RUSQSIEVE_OK = 0,
    RUSQSIEVE_INVALID_ARGUMENT = 1,
    RUSQSIEVE_INVALID_DECIMAL = 2,
    RUSQSIEVE_INPUT_OUT_OF_RANGE = 3,
    RUSQSIEVE_FACTORIZATION_FAILED = 4,
    RUSQSIEVE_INTERNAL_ERROR = 5,
    RUSQSIEVE_CANCELLED = 6
};

struct rusqsieve_progress {
    uint32_t phase;
    uint64_t completed;
    uint64_t total;
    uint32_t total_kind;
    uint32_t unit;
};

/* Return 0 to continue or nonzero to request cooperative cancellation. */
typedef int (*rusqsieve_progress_callback)(
    const struct rusqsieve_progress *progress,
    void *context
);

/* Return the ABI version implemented by this library (currently 3). */
uint32_t rusqsieve_abi_version(void);

/* Return a process-lifetime static description for a status code. */
const char *rusqsieve_strerror(int status);

/* Allocate a new empty result. */
rusqsieve_factors *rusqsieve_factors_new(void);

/*
 * Destroy a result and all strings it owns. A NULL pointer is ignored.
 * Every non-NULL result must be freed exactly once.
 */
void rusqsieve_factors_free(rusqsieve_factors *factors);

/*
 * Return the number of prime factors, including multiplicity.
 *
 * NULL and empty results have length zero.
 */
size_t rusqsieve_factors_len(const rusqsieve_factors *factors);

/*
 * Return a borrowed NUL-terminated decimal factor.
 *
 * Returns NULL when factors is NULL or index is out of bounds. The returned
 * pointer remains valid until the next rusqsieve_factor() call on this result
 * or until rusqsieve_factors_free().
 */
const char *rusqsieve_factors_get(
    const rusqsieve_factors *factors,
    size_t index
);

/*
 * Factor the positive base-10 integer n ("1" has an empty factorization).
 *
 * threads == 0 uses available parallelism capped at 48. A nonzero value is
 * capped at 256; the engine may cap smaller inputs further.
 *
 * The width of n is not itself limited: small factors are removed by trial
 * division, perfect-power detection, and Pollard-Brent, none of which is width
 * gated. RUSQSIEVE_INPUT_OUT_OF_RANGE is returned when n exceeds the library's
 * integer capacity, or when a composite wider than the quadratic sieve's
 * 400-bit range survives those stages and would have to be sieved.
 *
 * On success, factors is replaced with the complete sorted factorization.
 * Prime factors are repeated according to multiplicity; n == "1" succeeds
 * with rusqsieve_factors_len(factors) == 0. On failure, factors is left empty.
 *
 * n must point to a valid NUL-terminated string. factors must be a live object
 * returned by rusqsieve_factors_new().
 *
 * Calls that read, factor into, or free the same result must not overlap across
 * threads. Independent result objects may be used concurrently.
 */
enum rusqsieve_status rusqsieve_factor(
    const char *n,
    size_t threads,
    rusqsieve_factors *factors
);

/*
 * Optional behavior for rusqsieve_factor_ex().
 *
 * RUSQSIEVE_FLAG_ENABLE_ECM runs the elliptic curve method on composites the
 * quadratic sieve would otherwise handle by itself. ECM finds a medium-size
 * factor -- roughly 20 to 30 digits -- in time governed by the size of that
 * factor rather than of the input, which is the one shape neither Pollard-Brent
 * nor the sieve handles well. It is off by default because a curve cannot
 * succeed on a balanced semiprime and is pure overhead ahead of the sieve.
 *
 * The flag governs a narrower range than it looks like. Curves already run
 * without it on a composite wider than the sieve accepts, where the alternative
 * is refusing the number, and on one already known to be unbalanced because
 * trial division peeled a factor or Pollard-Brent split an ancestor. Neither
 * can happen for a balanced semiprime, which is what this flag is about.
 */
#define RUSQSIEVE_FLAG_ENABLE_ECM 1u

/*
 * Factor with optional behavior selected by a flag word.
 *
 * Unknown flag bits are ignored, so a caller built against a later flag set
 * still runs here. flags == 0 behaves exactly like rusqsieve_factor().
 */
enum rusqsieve_status rusqsieve_factor_ex(
    const char *n,
    size_t threads,
    uint32_t flags,
    rusqsieve_factors *factors
);

/*
 * Factor with progress and cancellation. callback may be NULL. The callback
 * runs on the calling thread; context is passed through unchanged.
 *
 * Stable phase codes are preprocessing=0, factor-base=1, sieving=2,
 * linear-algebra=6, extraction=7, complete=9.
 *
 * total_kind is unknown=0, exact=1, estimated=2. Unit codes are candidates=0,
 * primes=1, sieve-positions=3, relations=4, matrix-rows=5,
 * matrix-columns=6, matrix-nonzeros=7, iterations=8, matrix-products=9,
 * tasks=10.
 */
enum rusqsieve_status rusqsieve_factor_with_progress(
    const char *n,
    size_t threads,
    rusqsieve_factors *factors,
    rusqsieve_progress_callback callback,
    void *context
);

/*
 * Factor with both a flag word and progress reporting. The two are orthogonal,
 * so a caller never has to choose between selecting optional behavior and being
 * able to cancel. rusqsieve_factor_with_progress() is this with flags == 0, and
 * rusqsieve_factor_ex() is this with a NULL callback.
 */
enum rusqsieve_status rusqsieve_factor_ex_with_progress(
    const char *n,
    size_t threads,
    uint32_t flags,
    rusqsieve_factors *factors,
    rusqsieve_progress_callback callback,
    void *context
);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* RUSQSIEVE_H */
