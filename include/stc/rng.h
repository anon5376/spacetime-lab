#ifndef STC_RNG_H
#define STC_RNG_H

/* SPACETIME LAB — deterministic seeded RNGs.
 * mulberry32 (32-bit) + splitmix64 (64-bit). Uniform [0,1) + N(0,1).
 * No rand()/time() anywhere in sim code; seed explicitly.
 */

#include <stdint.h>

typedef struct {
    uint32_t state;
    int has_spare;
    double spare;
} StcRng32;

typedef struct {
    uint64_t state;
    int has_spare;
    double spare;
} StcRng64;

void stc_rng32_seed(StcRng32 *rng, uint32_t seed);
uint32_t stc_rng32_next(StcRng32 *rng);

void stc_rng64_seed(StcRng64 *rng, uint64_t seed);
uint64_t stc_rng64_next(StcRng64 *rng);

double stc_rng32_uniform(StcRng32 *rng); /* [0, 1) */
double stc_rng64_uniform(StcRng64 *rng); /* [0, 1) */
double stc_rng32_uniform_range(StcRng32 *rng, double lo, double hi);
double stc_rng64_uniform_range(StcRng64 *rng, double lo, double hi);
double stc_rng32_normal(StcRng32 *rng); /* N(0,1), Box-Muller */
double stc_rng64_normal(StcRng64 *rng); /* N(0,1), Box-Muller */

#endif /* STC_RNG_H */
