/* SPACETIME LAB -- nbody-test (section 11): body-core API/determinism tests.
 *
 * Own file: src/verify/nbody_test.c -> build/sttest (own main, no library).
 * May include body.h + core (vec/units/rng) + C stdlib only.
 * SI doubles; fixed seeds; no rand()/time(); no malloc; tiny and fast.
 *
 * Equation index (one per check):
 *   add_capacity   spawn bound: 64 slots ok, 65th returns -1 (STC_BODY_MAX).
 *   invalid_args   spawn guards: m > 0, r >= 0 else -1, n unchanged.
 *   remove_compact swap-with-last removal; out-of-range returns -1.
 *   presets        seeded presets 0..3 load (n > 0), bad index fails.
 *   determinism    same preset + 200 velocity-Verlet steps twice ->
 *                  bitwise-identical t and positions (memcmp + max diff).
 *   energies       total energy KE + PE finite before/after 200 steps.
 */

#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>

#include "stc/body.h"
#include "stc/units.h"
#include "stc/vec.h"

/* Fixed seed for all preset loads (determinism: same argv+seed = same sim). */
#define STC_TEST_SEED ((uint64_t)0x12345678u)

/* Determinism/energy step count (normative: 200 steps twice). */
#define STC_TEST_STEPS 200

/* Equation: per-preset stable dt — small fraction of each preset's
 * dynamical time (solar/binary ~days -> 1 hr; hole_disk sub-second orbits
 * -> 1e-4 s; slingshot cruise ~1e6 s -> 1e3 s). Determinism holds for any
 * dt; these keep the orbits finite so the finite-energy check is sane. */
static double stc_test_dt(int preset)
{
    switch (preset) {
    case 0:
        return 3600.0;
    case 1:
        return 3600.0;
    case 2:
        return 1.0e-4;
    case 3:
        return 1000.0;
    default:
        return 1.0;
    }
}

/* Equation: spawn bound — append STC_BODY_MAX bodies, each returns its
 * index; the 65th append returns -1 and n stays at 64. */
static int check_add_capacity(void)
{
    StcSim s;
    int ok = 1;
    int idx63 = -2;
    int rc65 = 0;
    int i;

    stc_sim_init(&s, STC_G, 0.0, 0.0, 1);
    for (i = 0; i < STC_BODY_MAX; i++) {
        StcVec3 p = { (double)i * 1.0e6, 0.0, 0.0 };
        StcVec3 v = { 0.0, 0.0, 0.0 };
        int rc = stc_sim_add(&s, "b", 1.0e10, 1.0, p, v, 0);

        if (i == STC_BODY_MAX - 1) {
            idx63 = rc;
        }
        ok = ok && (rc == i);
    }
    ok = ok && (s.n == STC_BODY_MAX);
    if (s.n == STC_BODY_MAX) {
        StcVec3 p = { 0.0, 0.0, 0.0 };
        StcVec3 v = { 0.0, 0.0, 0.0 };

        rc65 = stc_sim_add(&s, "overflow", 1.0e10, 1.0, p, v, 0);
        ok = ok && (rc65 == -1);
        ok = ok && (s.n == STC_BODY_MAX);
    } else {
        rc65 = -2;
        ok = 0;
    }
    printf("%s add_capacity: n=%d idx63=%d add65=%d expected=64 63 -1\n",
           ok ? "PASS" : "FAIL", s.n, idx63, rc65);
    return ok;
}

/* Equation: spawn guards — m > 0, r >= 0 required; NULL name allowed.
 * Rejections return -1 and leave n unchanged; a later valid add works. */
static int check_invalid_args(void)
{
    StcSim s;
    StcVec3 p = { 0.0, 0.0, 0.0 };
    StcVec3 v = { 0.0, 0.0, 0.0 };
    int ok = 1;
    int rc_m0;
    int rc_mneg;
    int rc_rneg;
    int rc_valid;

    stc_sim_init(&s, STC_G, 0.0, 0.0, 1);
    rc_m0 = stc_sim_add(&s, "m0", 0.0, 1.0, p, v, 0);
    rc_mneg = stc_sim_add(&s, "mneg", -1.0e10, 1.0, p, v, 0);
    rc_rneg = stc_sim_add(&s, "rneg", 1.0e10, -1.0, p, v, 0);
    ok = ok && (rc_m0 == -1);
    ok = ok && (rc_mneg == -1);
    ok = ok && (rc_rneg == -1);
    ok = ok && (s.n == 0);
    rc_valid = stc_sim_add(&s, "ok", 1.0e10, 0.0, p, v, 0);
    ok = ok && (rc_valid == 0);
    ok = ok && (s.n == 1);
    printf("%s invalid_args: m0=%d mneg=%d rneg=%d valid=%d n=%d "
           "expected=-1 -1 -1 0 1\n",
           ok ? "PASS" : "FAIL", rc_m0, rc_mneg, rc_rneg, rc_valid, s.n);
    return ok;
}

/* Equation: swap-with-last removal — removing idx 1 of [a,b,c] leaves
 * n == 2 with former last (c) moved into slot 1; bad idx returns -1. */
static int check_remove_compacts(void)
{
    StcSim s;
    StcVec3 p = { 0.0, 0.0, 0.0 };
    StcVec3 v = { 0.0, 0.0, 0.0 };
    int ok = 1;
    int rc_rm;
    int rc_bad_lo;
    int rc_bad_hi;
    double moved_m = NAN;

    stc_sim_init(&s, STC_G, 0.0, 0.0, 1);
    ok = ok && (stc_sim_add(&s, "a", 1.0e10, 1.0, p, v, 0) == 0);
    ok = ok && (stc_sim_add(&s, "b", 2.0e10, 1.0, p, v, 0) == 1);
    ok = ok && (stc_sim_add(&s, "c", 3.0e10, 1.0, p, v, 0) == 2);
    ok = ok && (s.n == 3);
    rc_rm = stc_sim_remove(&s, 1);
    ok = ok && (rc_rm == 0);
    ok = ok && (s.n == 2);
    if (s.n == 2) {
        moved_m = s.b[1].m;
        ok = ok && (s.b[1].m == 3.0e10);
        ok = ok && (strcmp(s.b[1].name, "c") == 0);
        ok = ok && (s.b[0].m == 1.0e10);
    } else {
        ok = 0;
    }
    rc_bad_lo = stc_sim_remove(&s, -1);
    rc_bad_hi = stc_sim_remove(&s, 5);
    ok = ok && (rc_bad_lo == -1);
    ok = ok && (rc_bad_hi == -1);
    ok = ok && (s.n == 2);
    printf("%s remove_compacts: n=%d moved_m=%.3g bad_lo=%d bad_hi=%d "
           "expected=2 3e+10 -1 -1\n",
           ok ? "PASS" : "FAIL", s.n, moved_m, rc_bad_lo, rc_bad_hi);
    return ok;
}

/* Equation: seeded presets — stc_preset_count() == 4, indices 0..3 load
 * with rc == 0 and n > 0; indices -1/4/99 fail nonzero. */
static int check_presets(void)
{
    int ok = 1;
    int count;
    int ns[4] = { -1, -1, -1, -1 };
    int rc_bad[3] = { 0, 0, 0 };
    const int bad_idx[3] = { -1, 4, 99 };
    int i;

    count = stc_preset_count();
    ok = ok && (count == 4);
    for (i = 0; i < 4; i++) {
        StcSim s;
        const char *nm;

        stc_sim_init(&s, STC_G, 1000.0, 0.4, 1);
        ok = ok && (stc_preset_load(&s, i, STC_TEST_SEED) == 0);
        ns[i] = s.n;
        ok = ok && (s.n > 0 && s.n <= STC_BODY_MAX);
        nm = stc_preset_name(i);
        ok = ok && (nm != NULL && nm[0] != '\0');
    }
    for (i = 0; i < 3; i++) {
        StcSim s;

        stc_sim_init(&s, STC_G, 1000.0, 0.4, 1);
        rc_bad[i] = stc_preset_load(&s, bad_idx[i], STC_TEST_SEED);
        ok = ok && (rc_bad[i] != 0);
    }
    printf("%s presets: count=%d n=[%d %d %d %d] bad=[%d %d %d] "
           "expected=4 >0 >0 >0 >0 nonzero nonzero nonzero\n",
           ok ? "PASS" : "FAIL", count, ns[0], ns[1], ns[2], ns[3],
           rc_bad[0], rc_bad[1], rc_bad[2]);
    return ok;
}

/* Equation: velocity-Verlet determinism — same preset + seed stepped 200
 * times twice gives bitwise-identical t and positions (per-double memcmp);
 * maxdiff is the max |component| difference (0.0 when identical). */
static int check_determinism(void)
{
    int ok = 1;
    int memcmp_ok = 1;
    double maxdiff = 0.0;
    int i;

    for (i = 0; i < 4; i++) {
        StcSim a;
        StcSim b;
        double dt = stc_test_dt(i);
        int k;
        int j;

        stc_sim_init(&a, STC_G, 1000.0, 0.4, 1);
        stc_sim_init(&b, STC_G, 1000.0, 0.4, 1);
        ok = ok && (stc_preset_load(&a, i, STC_TEST_SEED) == 0);
        ok = ok && (stc_preset_load(&b, i, STC_TEST_SEED) == 0);
        ok = ok && (a.n == b.n);
        for (k = 0; k < STC_TEST_STEPS; k++) {
            stc_sim_step(&a, dt);
            stc_sim_step(&b, dt);
        }
        ok = ok && (a.n == b.n);
        if (a.n != b.n) {
            memcmp_ok = 0;
            maxdiff = HUGE_VAL;
            continue;
        }
        if (memcmp(&a.t, &b.t, sizeof(double)) != 0) {
            memcmp_ok = 0;
        }
        {
            double d = fabs(a.t - b.t);

            if (!(d <= maxdiff)) {
                maxdiff = d;
            }
        }
        for (j = 0; j < a.n; j++) {
            double dx;
            double dy;
            double dz;

            if (memcmp(&a.b[j].p.x, &b.b[j].p.x, sizeof(double)) != 0 ||
                memcmp(&a.b[j].p.y, &b.b[j].p.y, sizeof(double)) != 0 ||
                memcmp(&a.b[j].p.z, &b.b[j].p.z, sizeof(double)) != 0) {
                memcmp_ok = 0;
            }
            dx = fabs(a.b[j].p.x - b.b[j].p.x);
            dy = fabs(a.b[j].p.y - b.b[j].p.y);
            dz = fabs(a.b[j].p.z - b.b[j].p.z);
            if (!(dx <= maxdiff)) {
                maxdiff = dx;
            }
            if (!(dy <= maxdiff)) {
                maxdiff = dy;
            }
            if (!(dz <= maxdiff)) {
                maxdiff = dz;
            }
        }
    }
    ok = ok && memcmp_ok && (maxdiff == 0.0);
    printf("%s determinism: steps=%d maxdiff=%.3e memcmp=%s "
           "expected=0.000e+00 identical\n",
           ok ? "PASS" : "FAIL", STC_TEST_STEPS, maxdiff,
           memcmp_ok ? "identical" : "DIFF");
    return ok;
}

/* Equation: total energy KE (free bodies) + pairwise PE — finite for
 * every preset before and after 200 steps. */
static int check_energies_finite(void)
{
    int ok = 1;
    double e_before[4] = { NAN, NAN, NAN, NAN };
    double e_after[4] = { NAN, NAN, NAN, NAN };
    int i;

    for (i = 0; i < 4; i++) {
        StcSim s;
        double dt = stc_test_dt(i);
        int k;

        stc_sim_init(&s, STC_G, 1000.0, 0.4, 1);
        ok = ok && (stc_preset_load(&s, i, STC_TEST_SEED) == 0);
        e_before[i] = stc_sim_energy(&s);
        ok = ok && isfinite(e_before[i]);
        for (k = 0; k < STC_TEST_STEPS; k++) {
            stc_sim_step(&s, dt);
        }
        e_after[i] = stc_sim_energy(&s);
        ok = ok && isfinite(e_after[i]);
    }
    printf("%s energies_finite: before=[%.3g %.3g %.3g %.3g] "
           "after=[%.3g %.3g %.3g %.3g] expected=all finite\n",
           ok ? "PASS" : "FAIL", e_before[0], e_before[1], e_before[2],
           e_before[3], e_after[0], e_after[1], e_after[2], e_after[3]);
    return ok;
}

/* Equation: test dispatch — run the 6 API/determinism checks, exit 0 iff
 * all pass. */
int main(void)
{
    int passed = 0;

    passed += check_add_capacity();
    passed += check_invalid_args();
    passed += check_remove_compacts();
    passed += check_presets();
    passed += check_determinism();
    passed += check_energies_finite();
    printf("sttest: %d/6 PASS\n", passed);
    return (passed == 6) ? 0 : 1;
}
