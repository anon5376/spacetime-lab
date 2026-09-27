/* SPACETIME LAB — fuzz (section 13): seeded randomized N-body finiteness sweep.
 *
 * Own file: src/verify/fuzz.c -> build/stfuzz (own main, no library).
 * Includes only stc/body.h + stc core (vec via body.h, rng, units) + C stdlib.
 *
 * Method: 200 deterministic cases (StcRng64 splitmix64, seed = base + case).
 * Each case builds a random sim with N = 2..12 bodies, log-uniform masses
 * (1e20..1e30 kg), radii (1e5..1e8 m), positions (+/-1e6..1e11 m),
 * velocities (+/-1e1..1e5 m/s) and timestep (1e0..1e5 s), then integrates
 * 50 velocity-Verlet steps (stc_sim_step) with collisions on (merge mode on
 * even cases, bounce mode with e = 0.6 on odd cases — both are the body's
 * collision pass, alternating only to cover both paths). After every step
 * each live body's position/velocity components and the total sim energy
 * must be finite; any violation prints a failure line and the run exits 1.
 * Prints ALL_CLEAN iff all 200 cases are clean. SI doubles, no I/O besides
 * the verdict lines, no malloc, no rand()/time().
 */

#include <math.h>
#include <stdint.h>
#include <stdio.h>

#include "stc/body.h"
#include "stc/rng.h"
#include "stc/units.h"

/* Sweep shape (all fixed for determinism). */
#define STF_NCASES 200
#define STF_NSTEPS 50
#define STF_NMIN 2
#define STF_NMAX 12

/* Base seed; case i uses STF_SEED_BASE + i. */
#define STF_SEED_BASE 20260913ULL

/* Log-uniform bands (base-10 exponents). */
#define STF_M_LO 20.0
#define STF_M_HI 30.0
#define STF_R_LO 5.0
#define STF_R_HI 8.0
#define STF_P_LO 6.0
#define STF_P_HI 11.0
#define STF_V_LO 1.0
#define STF_V_HI 5.0
#define STF_DT_LO 0.0
#define STF_DT_HI 5.0

/* Fixed sim knobs: Plummer softening (m) + bounce restitution. */
#define STF_SOFT 1e3
#define STF_RESTITUTION 0.6

/* Equation: log-uniform draw — 10^U(lo,hi); uniform in the exponent so each
 * decade is equally likely. */
static double stf_log_uniform(StcRng64 *rng, double lo_exp, double hi_exp)
{
    double u = stc_rng64_uniform_range(rng, lo_exp, hi_exp);

    return pow(10.0, u);
}

/* Equation: signed log-uniform draw — random sign times log-uniform. */
static double stf_signed_log_uniform(StcRng64 *rng, double lo_exp,
                                     double hi_exp)
{
    double mag = stf_log_uniform(rng, lo_exp, hi_exp);

    return (stc_rng64_uniform(rng) < 0.5) ? -mag : mag;
}

/* Check one case; prints FAIL lines, returns count of failing steps. */
static int stf_run_case(int icase, uint64_t seed)
{
    StcRng64 rng;
    StcSim s;
    char name[32];
    int i, step, n, fails = 0;
    double dt;

    stc_rng64_seed(&rng, seed);

    /* Equation: velocity-Verlet N-body sim (stc_sim_step); collisions on:
     * merge on even cases, bounce (e = 0.6) on odd cases. */
    stc_sim_init(&s, STC_G, STF_SOFT, STF_RESTITUTION, (icase % 2) == 0);

    n = STF_NMIN + (int)(stc_rng64_next(&rng) % (uint64_t)(STF_NMAX - STF_NMIN + 1));
    for (i = 0; i < n; i++) {
        StcVec3 p, v;
        double m = stf_log_uniform(&rng, STF_M_LO, STF_M_HI);
        double r = stf_log_uniform(&rng, STF_R_LO, STF_R_HI);

        p.x = stf_signed_log_uniform(&rng, STF_P_LO, STF_P_HI);
        p.y = stf_signed_log_uniform(&rng, STF_P_LO, STF_P_HI);
        p.z = stf_signed_log_uniform(&rng, STF_P_LO, STF_P_HI);
        v.x = stf_signed_log_uniform(&rng, STF_V_LO, STF_V_HI);
        v.y = stf_signed_log_uniform(&rng, STF_V_LO, STF_V_HI);
        v.z = stf_signed_log_uniform(&rng, STF_V_LO, STF_V_HI);

        (void)snprintf(name, sizeof(name), "f%d_%d", icase, i);
        if (stc_sim_add(&s, name, m, r, p, v, 0) < 0) {
            (void)printf("FAIL case=%d seed=%llu step=-1 body=%d what=add_rejected\n",
                         icase, (unsigned long long)seed, i);
            return 1;
        }
    }

    dt = stf_log_uniform(&rng, STF_DT_LO, STF_DT_HI);

    for (step = 0; step < STF_NSTEPS; step++) {
        double e;

        stc_sim_step(&s, dt);

        for (i = 0; i < s.n; i++) {
            const StcVec3 *p = &s.b[i].p;
            const StcVec3 *v = &s.b[i].v;

            if (!isfinite(p->x) || !isfinite(p->y) || !isfinite(p->z)) {
                (void)printf("FAIL case=%d seed=%llu step=%d body=%d what=pos_nonfinite\n",
                             icase, (unsigned long long)seed, step, i);
                fails++;
                break;
            }
            if (!isfinite(v->x) || !isfinite(v->y) || !isfinite(v->z)) {
                (void)printf("FAIL case=%d seed=%llu step=%d body=%d what=vel_nonfinite\n",
                             icase, (unsigned long long)seed, step, i);
                fails++;
                break;
            }
        }
        if (fails > 0)
            break;

        /* Equation: total energy — kinetic + pairwise potential. */
        e = stc_sim_energy(&s);
        if (!isfinite(e)) {
            (void)printf("FAIL case=%d seed=%llu step=%d body=-1 what=energy_nonfinite\n",
                         icase, (unsigned long long)seed, step);
            fails++;
            break;
        }
    }

    return fails;
}

/* Equation: none (harness main) — runs the 200-case sweep, exit 0 iff all
 * positions/velocities/energies stayed finite. */
int main(void)
{
    int icase, fails = 0;

    for (icase = 0; icase < STF_NCASES; icase++)
        fails += stf_run_case(icase, STF_SEED_BASE + (uint64_t)icase);

    if (fails == 0)
        (void)printf("ALL_CLEAN cases=%d steps_each=%d\n", STF_NCASES, STF_NSTEPS);

    return (fails == 0) ? 0 : 1;
}
