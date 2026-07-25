#include "../rusqsieve.h"

#include <assert.h>
#include <stddef.h>
#include <string.h>

int main(void) {
    rusqsieve_factors *factors = rusqsieve_factors_new();
    assert(factors != NULL);

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
    rusqsieve_factors_free(factors);
    rusqsieve_factors_free(NULL);
    return 0;
}
