#include "stc/rng.h"

#include <math.h>
#include <stddef.h>

void stc_rng32_seed(StcRng32 *rng, uint32_t seed)
{
    if (rng == NULL) {
        return;
    }
    rng->state = seed;
    rng->has_spare = 0;
    rng->spare = 0.0;
}

uint32_t stc_rng32_next(StcRng32 *rng)
{
    /* mulberry32 */
    uint32_t z;
    if (rng == NULL) {
        return 0u;
    }
    rng->state += 0x6D2B79F5u;
    z = rng->state;
    z = (z ^ (z >> 15)) * (z | 1u);
    z ^= z + (z ^ (z >> 7)) * (z | 61u);
    return z ^ (z >> 14);
}

void stc_rng64_seed(StcRng64 *rng, uint64_t seed)
{
    if (rng == NULL) {
        return;
    }
    rng->state = seed;
    rng->has_spare = 0;
    rng->spare = 0.0;
}

uint64_t stc_rng64_next(StcRng64 *rng)
{
    /* splitmix64 */
    uint64_t z;
    if (rng == NULL) {
        return 0u;
    }
    rng->state += 0x9E3779B97F4A7C15ULL;
    z = rng->state;
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
    return z ^ (z >> 31);
}

double stc_rng32_uniform(StcRng32 *rng)
{
    return (double)stc_rng32_next(rng) / 4294967296.0;
}

double stc_rng64_uniform(StcRng64 *rng)
{
    /* Top 53 bits -> [0,1). */
    uint64_t v = stc_rng64_next(rng) >> 11;
    return (double)v / 9007199254740992.0;
}

double stc_rng32_uniform_range(StcRng32 *rng, double lo, double hi)
{
    return lo + (hi - lo) * stc_rng32_uniform(rng);
}

double stc_rng64_uniform_range(StcRng64 *rng, double lo, double hi)
{
    return lo + (hi - lo) * stc_rng64_uniform(rng);
}

double stc_rng32_normal(StcRng32 *rng)
{
    /* Box-Muller with cached spare. Deterministic per seed. */
    double u1, u2, mag;
    if (rng == NULL) {
        return 0.0;
    }
    if (rng->has_spare) {
        rng->has_spare = 0;
        return rng->spare;
    }
    do {
        u1 = stc_rng32_uniform(rng);
    } while (u1 <= 0.0);
    u2 = stc_rng32_uniform(rng);
    mag = sqrt(-2.0 * log(u1));
    rng->spare = mag * sin(2.0 * 3.14159265358979323846 * u2);
    rng->has_spare = 1;
    return mag * cos(2.0 * 3.14159265358979323846 * u2);
}

double stc_rng64_normal(StcRng64 *rng)
{
    double u1, u2, mag;
    if (rng == NULL) {
        return 0.0;
    }
    if (rng->has_spare) {
        rng->has_spare = 0;
        return rng->spare;
    }
    do {
        u1 = stc_rng64_uniform(rng);
    } while (u1 <= 0.0);
    u2 = stc_rng64_uniform(rng);
    mag = sqrt(-2.0 * log(u1));
    rng->spare = mag * sin(2.0 * 3.14159265358979323846 * u2);
    rng->has_spare = 1;
    return mag * cos(2.0 * 3.14159265358979323846 * u2);
}
