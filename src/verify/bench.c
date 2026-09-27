/* SPACETIME LAB — bench: N-body step throughput benchmark.
 *
 * Own file: src/verify/bench.c -> build/stbench (own main, no library).
 * May include stc/body.h + stc/core (vec/units) + C stdlib only.
 *
 * Method: analytic circular ring — one fixed central mass plus N-1
 * tracers on a circle (XZ plane, y-up), each at circular velocity
 * v = sqrt(G*M/R) (cite: circular velocity); steps the Verlet core
 * stc_sim_step 200 times per config and reports wall-clock throughput.
 * Info only: always exits 0.
 *
 * SI doubles throughout; analytic layout (no RNG); timing via
 * clock() is measurement only, not sim input.
 */

#include <math.h>
#include <stdio.h>
#include <time.h>

#include "stc/body.h"
#include "stc/units.h"
#include "stc/vec.h"

/* Steps per ring config (spec: 200). */
#define STC_BENCH_STEPS 200

/* Ring radius in meters. */
#define STC_BENCH_RING_R 1.0e9

/* Central fixed mass in kg. */
#define STC_BENCH_CENTER_M 1.0e30

/* Tracer mass in kg (negligible vs center, keeps ring stable). */
#define STC_BENCH_RING_M 1.0e20

/* Method: circular velocity v = sqrt(mu/R). */
static double bench_circular_v(double mu, double r)
{
    return sqrt(mu / r);
}

/* Method: Kepler-III period T = 2*pi*sqrt(R^3/mu). */
static double bench_period(double mu, double r)
{
    return 2.0 * STC_PI * sqrt(r * r * r / mu);
}

/* Method: analytic ring build — fixed center + (n-1) circular tracers.
 * Returns dt for the run, or a negative value on build failure. */
static double bench_build(StcSim *s, int n)
{
    StcVec3 pc = { 0.0, 0.0, 0.0 };
    StcVec3 vc = { 0.0, 0.0, 0.0 };
    double mu = STC_G * STC_BENCH_CENTER_M;
    double vcirc = bench_circular_v(mu, STC_BENCH_RING_R);
    double period = bench_period(mu, STC_BENCH_RING_R);
    int k;

    stc_sim_init(s, STC_G, 0.0, 0.0, 0);
    s->vmax = 0.0;
    if (stc_sim_add(s, "center", STC_BENCH_CENTER_M, 0.0, pc, vc, 1) != 0) {
        return -1.0;
    }
    for (k = 0; k < n - 1; k++) {
        double th = 2.0 * STC_PI * (double)k / (double)(n - 1);
        StcVec3 p = { 0.0, 0.0, 0.0 };
        StcVec3 v = { 0.0, 0.0, 0.0 };

        /* XZ plane, y-up; prograde (+Y angular momentum). */
        p.x = STC_BENCH_RING_R * cos(th);
        p.z = STC_BENCH_RING_R * sin(th);
        v.x = vcirc * sin(th);
        v.z = -vcirc * cos(th);
        if (stc_sim_add(s, "ring", STC_BENCH_RING_M, 0.0, p, v, 0) < 0) {
            return -1.0;
        }
    }
    /* Small stable step: 1/2000 of an orbit (200 steps = 0.1 orbit). */
    return period / 2000.0;
}

/* Method: clock-tick delta to seconds (one-tick floor when too fast). */
static double bench_elapsed(clock_t t0, clock_t t1)
{
    double dt;

    if (t1 <= t0) {
        return 1.0 / (double)CLOCKS_PER_SEC;
    }
    dt = (double)(t1 - t0) / (double)CLOCKS_PER_SEC;
    if (!(dt > 0.0)) {
        return 1.0 / (double)CLOCKS_PER_SEC;
    }
    return dt;
}

/* Method: timed Verlet run — 200 stc_sim_step calls, throughput row. */
static void bench_run(int n)
{
    StcSim s;
    double dt = bench_build(&s, n);
    clock_t t0;
    clock_t t1;
    double el;
    double steps_per_s;
    double pairs_per_s;
    double pairs;
    int i;

    if (!(dt > 0.0)) {
        printf("%3d %5d %12s %14s %10s\n", n, STC_BENCH_STEPS,
               "BUILD_FAIL", "BUILD_FAIL", "n/a");
        return;
    }
    t0 = clock();
    for (i = 0; i < STC_BENCH_STEPS; i++) {
        stc_sim_step(&s, dt);
    }
    t1 = clock();
    el = bench_elapsed(t0, t1);
    steps_per_s = (double)STC_BENCH_STEPS / el;
    pairs = (double)n * (double)(n - 1) / 2.0;
    pairs_per_s = pairs * (double)STC_BENCH_STEPS / el;
    printf("%3d %5d %12.1f %14.1f %10.4f\n", n, STC_BENCH_STEPS,
           steps_per_s, pairs_per_s, el);
}

/* Method: bench driver — ring configs N = 8/16/32/64, info table. */
int main(void)
{
    static const int ns[] = { 8, 16, 32, 64 };
    unsigned int i;

    printf("stbench: ring N-body benchmark (%d steps/config, info only)\n",
           STC_BENCH_STEPS);
    printf("%3s %5s %12s %14s %10s\n", "N", "steps", "steps/s",
           "pairs/s", "elapsed_s");
    for (i = 0; i < sizeof(ns) / sizeof(ns[0]); i++) {
        bench_run(ns[i]);
    }
    return 0;
}
