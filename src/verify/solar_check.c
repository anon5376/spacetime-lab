/* SPACETIME LAB — solar_check (sandbox swarm section 12).
 *
 * Own file: src/verify/solar_check.c -> build/stsolar (OWN main, no library).
 * Includes only stc/body.h + stc/core (vec/units) + C stdlib.
 *
 * Method: load solar preset (index 0), integrate dt = 3600 s for 8760 steps
 * (~1 yr) with the shipped velocity-Verlet core, then check:
 *   radius drift  d = |r1 - r0| / r0 < 2% per planet, where
 *                   r = |p_planet - p_sun| (Euclidean norm, SI meters);
 *                   r0 equals the initial a for the circular preset starts.
 *   Earth angle   cumulative swept angle S = sum_k atan2(|u_k x v_k|,
 *                   u_k . v_k) over successive Sun-relative Earth vectors
 *                   (plane-independent pairwise angle), |S - 2pi| / 2pi < 2%.
 *   energy drift  |E1 - E0| / |E0| < 0.5% via stc_sim_energy.
 * Prints a per-planet table + verdict lines. Exit 0 iff every check passes.
 * Deterministic: fixed preset index + fixed seed, no RNG, no wall clock.
 */

#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>

#include "stc/body.h"
#include "stc/units.h"
#include "stc/vec.h"

/* Integration + tolerance parameters (spec section 12). */
#define STSOLAR_DT 3600.0
#define STSOLAR_STEPS 8760
#define STSOLAR_SEED 0ULL
#define STSOLAR_RADIUS_TOL 0.02
#define STSOLAR_ANGLE_TOL 0.02
#define STSOLAR_ENERGY_TOL 0.005

/* Equation: ASCII lowercase (locale-independent helper for name matching). */
static char stsolar_ascii_lower(char c)
{
    if (c >= 'A' && c <= 'Z') {
        return (char)(c - 'A' + 'a');
    }
    return c;
}

/* Equation: bounded case-insensitive name equality (names are char[32]). */
static int stsolar_name_eq_ci(const char *a, const char *b)
{
    int i;

    for (i = 0; i < 32; i++) {
        char ca = stsolar_ascii_lower(a[i]);
        char cb = stsolar_ascii_lower(b[i]);

        if (ca != cb) {
            return 0;
        }
        if (ca == '\0') {
            return 1;
        }
    }
    return 1;
}

/* Equation: bounded exact-name lookup over live bodies; index or -1. */
static int stsolar_find(const StcSim *s, const char *name)
{
    int i;

    for (i = 0; i < s->n; i++) {
        if (strncmp(s->b[i].name, name, 32) == 0) {
            return i;
        }
    }
    return -1;
}

/* Equation: pairwise sweep angle dtheta = atan2(|u x v|, u . v) in [0, pi]. */
static double stsolar_step_angle(StcVec3 u, StcVec3 v)
{
    double lu = stc_vec3_len(u);
    double lv = stc_vec3_len(v);
    StcVec3 c;

    if (!(lu > 0.0) || !(lv > 0.0)) {
        return 0.0;
    }
    if (!isfinite(lu) || !isfinite(lv)) {
        return 0.0;
    }
    c = stc_vec3_cross(u, v);
    return atan2(stc_vec3_len(c), stc_vec3_dot(u, v));
}

/* Equation: 1-yr solar-system acceptance run (see file header for methods). */
int main(void)
{
    StcSim s;
    char sun_name[32];
    char names[STC_BODY_MAX][32];
    char earth_name[32];
    double r0[STC_BODY_MAX];
    const char *pname;
    double e0;
    double e1;
    double swept;
    double want;
    double aerr;
    double eerr;
    int n0;
    int i;
    int ok;
    int nplanet;
    int nplanet_ok;
    int angle_ok;
    int energy_ok;
    int all_ok;

    /* Mirror the shipped sandbox defaults, then load (app init-load pattern). */
    stc_sim_init(&s, STC_G, 1000.0, 0.4, 1);
    if (stc_preset_count() < 1 || stc_preset_load(&s, 0, STSOLAR_SEED) != 0) {
        printf("FAIL solar_load: preset 0 unavailable\n");
        printf("stsolar: FAIL\n");
        return 1;
    }
    n0 = s.n;
    pname = stc_preset_name(0);
    printf("stsolar: preset \"%s\" n=%d dt=%.0fs steps=%d seed=%llu\n",
        pname ? pname : "?", n0, STSOLAR_DT, STSOLAR_STEPS,
        (unsigned long long)STSOLAR_SEED);
    if (n0 < 2 || n0 > STC_BODY_MAX) {
        printf("FAIL solar_load: n=%d, need 2..%d\n", n0, STC_BODY_MAX);
        printf("stsolar: FAIL\n");
        return 1;
    }
    memcpy(sun_name, s.b[0].name, 32);
    for (i = 1; i < n0; i++) {
        memcpy(names[i], s.b[i].name, 32);
        r0[i] = stc_vec3_dist(s.b[i].p, s.b[0].p);
    }
    /* Earth by name; fallback to index 3 (Sun + Mercury, Venus, Earth). */
    {
        int ei = -1;

        for (i = 1; i < n0; i++) {
            if (stsolar_name_eq_ci(names[i], "earth")) {
                ei = i;
                break;
            }
        }
        if (ei < 0) {
            ei = (n0 > 3) ? 3 : 1;
        }
        memcpy(earth_name, names[ei], 32);
    }
    e0 = stc_sim_energy(&s);

    swept = 0.0;
    ok = 1;
    ok = ok && isfinite(e0);
    for (i = 0; i < STSOLAR_STEPS; i++) {
        int si = stsolar_find(&s, sun_name);
        int ei = stsolar_find(&s, earth_name);
        StcVec3 prev;
        StcVec3 cur;

        if (si < 0 || ei < 0) {
            ok = 0;
            break;
        }
        prev = stc_vec3_sub(s.b[ei].p, s.b[si].p);
        stc_sim_step(&s, STSOLAR_DT);
        si = stsolar_find(&s, sun_name);
        ei = stsolar_find(&s, earth_name);
        if (si < 0 || ei < 0) {
            ok = 0;
            break;
        }
        cur = stc_vec3_sub(s.b[ei].p, s.b[si].p);
        swept += stsolar_step_angle(prev, cur);
        if (!isfinite(swept)) {
            ok = 0;
            break;
        }
    }
    if (s.n != n0) {
        printf("note: n changed %d -> %d during run (collision merge?)\n",
            n0, s.n);
        ok = 0;
    }
    e1 = stc_sim_energy(&s);

    printf("%-8s %12s %12s %9s %s\n", "planet", "a0[AU]", "r1[AU]",
        "drift[%]", "verdict");
    nplanet = 0;
    nplanet_ok = 0;
    for (i = 1; i < n0; i++) {
        int si = stsolar_find(&s, sun_name);
        int bi = stsolar_find(&s, names[i]);
        double r1 = NAN;
        double drift = NAN;
        int rok = 0;

        if (si >= 0 && bi >= 0) {
            r1 = stc_vec3_dist(s.b[bi].p, s.b[si].p);
        }
        if (isfinite(r1) && isfinite(r0[i]) && r0[i] > 0.0) {
            drift = fabs(r1 - r0[i]) / r0[i];
            rok = (drift < STSOLAR_RADIUS_TOL);
        }
        nplanet++;
        nplanet_ok += rok ? 1 : 0;
        printf("%-8.32s %12.6f %12.6f %9.4f %s\n", names[i],
            r0[i] / STC_AU, r1 / STC_AU, drift * 100.0,
            rok ? "PASS" : "FAIL");
    }

    want = 2.0 * STC_PI;
    aerr = fabs(swept - want) / want;
    angle_ok = ok && isfinite(swept) && (aerr < STSOLAR_ANGLE_TOL);
    printf("%s earth_angle: swept=%.6f rad expected=%.6f err=%.4f%% tol=2%%\n",
        angle_ok ? "PASS" : "FAIL", swept, want, aerr * 100.0);

    eerr = 1.0;
    energy_ok = 0;
    if (isfinite(e0) && isfinite(e1) && e0 != 0.0) {
        eerr = fabs(e1 - e0) / fabs(e0);
        energy_ok = ok && (eerr < STSOLAR_ENERGY_TOL);
    }
    printf("%s energy: E0=%.6e E1=%.6e |dE/E|=%.4f%% tol=0.5%%\n",
        energy_ok ? "PASS" : "FAIL", e0, e1, eerr * 100.0);

    all_ok = ok && (nplanet_ok == nplanet) && angle_ok && energy_ok;
    printf("stsolar: %s (%d/%d planets radius<2%%, angle %s, energy %s)\n",
        all_ok ? "PASS" : "FAIL", nplanet_ok, nplanet,
        angle_ok ? "PASS" : "FAIL", energy_ok ? "PASS" : "FAIL");
    return all_ok ? 0 : 1;
}
