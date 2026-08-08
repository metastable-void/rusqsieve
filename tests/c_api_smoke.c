#include "../rusqsieve.h"

#include <assert.h>
#include <stddef.h>
#include <string.h>

static int cancel_immediately(
    const struct rusqsieve_progress *progress,
    void *context
) {
    unsigned *calls = context;
    assert(progress != NULL);
    ++*calls;
    return 1;
}

int main(void) {
    rusqsieve_factors *factors = rusqsieve_factors_new();
    assert(factors != NULL);
    assert(rusqsieve_abi_version() == 3);
    assert(strcmp(rusqsieve_strerror(RUSQSIEVE_OK), "success") == 0);

    assert(rusqsieve_factor("360", 1, factors) == RUSQSIEVE_OK);
    assert(rusqsieve_factors_len(factors) == 6);
    assert(strcmp(rusqsieve_factors_get(factors, 0), "2") == 0);
    assert(strcmp(rusqsieve_factors_get(factors, 1), "2") == 0);
    assert(strcmp(rusqsieve_factors_get(factors, 2), "2") == 0);
    assert(strcmp(rusqsieve_factors_get(factors, 3), "3") == 0);
    assert(strcmp(rusqsieve_factors_get(factors, 4), "3") == 0);
    assert(strcmp(rusqsieve_factors_get(factors, 5), "5") == 0);
    assert(rusqsieve_factors_get(factors, 6) == NULL);

    /* A result object is reusable; factor() releases its previous contents. */
    assert(rusqsieve_factor("1", 0, factors) == RUSQSIEVE_OK);
    assert(rusqsieve_factors_len(factors) == 0);
    assert(rusqsieve_factors_get(factors, 0) == NULL);

    assert(rusqsieve_factor("invalid", 1, factors) ==
           RUSQSIEVE_INVALID_DECIMAL);
    assert(rusqsieve_factors_len(factors) == 0);

    /* The flag word is the only difference between the two entry points, and an unknown bit must
       not change the answer. */
    assert(rusqsieve_factor_ex("1000036000099", 1, 0, factors) == RUSQSIEVE_OK);
    assert(rusqsieve_factors_len(factors) == 2);
    assert(rusqsieve_factor_ex("1000036000099", 1, RUSQSIEVE_FLAG_ENABLE_ECM,
                               factors) == RUSQSIEVE_OK);
    assert(rusqsieve_factors_len(factors) == 2);
    assert(rusqsieve_factor_ex("1000036000099", 1, 0x8000u, factors) ==
           RUSQSIEVE_OK);
    assert(rusqsieve_factors_len(factors) == 2);
    assert(rusqsieve_factor_ex(NULL, 1, 0, factors) ==
           RUSQSIEVE_INVALID_ARGUMENT);

    /* Flags and progress are orthogonal, so the combined entry point must honour both. */
    unsigned ex_calls = 0;
    assert(rusqsieve_factor_ex_with_progress("1000036000099", 1,
                                             RUSQSIEVE_FLAG_ENABLE_ECM, factors,
                                             cancel_immediately, &ex_calls) ==
           RUSQSIEVE_CANCELLED);
    assert(ex_calls == 1);
    assert(rusqsieve_factors_len(factors) == 0);
    assert(rusqsieve_factor_ex_with_progress("1000036000099", 1, 0, factors,
                                             NULL, NULL) == RUSQSIEVE_OK);
    assert(rusqsieve_factors_len(factors) == 2);

    unsigned calls = 0;
    assert(rusqsieve_factor_with_progress(
               "1000036000099", 1, factors, cancel_immediately, &calls) ==
           RUSQSIEVE_CANCELLED);
    assert(calls == 1);
    assert(rusqsieve_factors_len(factors) == 0);
    rusqsieve_factors_free(factors);
    rusqsieve_factors_free(NULL);
    return 0;
}
