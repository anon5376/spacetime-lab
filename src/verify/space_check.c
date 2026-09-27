/* SPACETIME LAB — space_check (SWE2 worker swe2-space): SPACE-CORE
 * acceptance checks -> build/stspace.
 *
 * Every PASS/FAIL line consumes the shipped stc_space_* API (R1).
 * Oracles (analytic Kepler period, Einstein perihelion formula) are
 * second opinions only.
 *
 * Checks:
 *   kepler      circular 1 AU orbit, h = 1 d, one period: position
 *               error < 1e-9 AU, relative energy error < 1e-12.
 *   order       e = 0.6 two-body at h = P/50..P/400: observed
 *               convergence order in [5.5, 6.5].
 *   reversible  10-body seeded bound system, pn=0, collisions=0:
 *               +100 orbits then -100 -> back to start, < 1e-9 rel.
 *   mercury     Sun (STAR) + Mercury, pn=1: perihelion advance over
 *               100 Julian years = 42.98" +- 1.5%; pn=0 -> < 0.05".
 *   merge       head-on contact: 1 survivor, mass/momentum conserved
 *               to 1e-12 rel, MERGE logged.
 *   blackhole   10 Msun BH + plunging 1 Msun star: DISRUPT spawns
 *               disrupt_tracers tracers; aimed tracer -> CAPTURE.
 *   seek        merge run, ckpt every 5 d: seek(10 d)/seek(50 d)
 *               bitwise-equal the recorded states (memcmp).
 *   j2_node     Io-like moon about an oblate planet (J2 0.0147,
 *               R 71492 km, i 0.5 deg): measured nodal regression over
 *               200 d = -(3/2) n J2 (R/a)^2 cos i within 1%.
 *   j2_momentum oblate planet + massive moons: total momentum drift
 *               < 1e-13 relative over 30 d (zonal reaction included).
 *   pn_step     Sun+Mercury+Venus+Earth, pn=1, 10 yr at h=6 h vs
 *               h=1 h: Mercury positions agree < 1 km (implicit-
 *               midpoint kick keeps the Yoshida order).
 *   threads     same 2000-body tree run under STC_PAR_THREADS=1
 *               (re-exec'd child) and default threads -> same hash.
 *   tree        5000-body Plummer sphere: theta 0.5 vs direct, median
 *               rel force error < 1e-3, p99 < 1e-2.
 *   perf        INFO-only substep timings.
 *
 * Claim class of every result here: MEASURED (analytic orbit formulae
 * are closed-form expectations used only as toleranced targets).
 */

#define _POSIX_C_SOURCE 200809L

#include <math.h>
#include <stdarg.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "stc/rng.h"
#include "stc/space.h"
#include "stc/units.h"

static int g_fails;

static void report(int ok, const char *name, double computed,
                   double expected, double tol)
{
    printf("%s %s: computed=%.10g expected=%.10g tol=%.6g\n",
           ok ? "PASS" : "FAIL", name, computed, expected, tol);
    if (!ok) {
        g_fails++;
    }
}

static void info(const char *name, const char *fmt, ...)
{
    char buf[256];
    va_list ap;
    va_start(ap, fmt);
    vsnprintf(buf, sizeof(buf), fmt, ap);
    va_end(ap);
    printf("INFO %s: %s\n", name, buf);
}

static double now_s(void)
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + 1e-9 * (double)ts.tv_nsec;
}

/* --------------------------------------------------------------- */
/* shared two-body setup: returns relative-orbit barycentric state   */

static double two_body_period(double a, double mu)
{
    return 2.0 * STC_PI * sqrt(a * a * a / mu);
}

/* Put (primary, orbiter) into barycentric orbit: separation `sep` on
 * +x at periapsis with the given eccentricity (e=0 -> circular). */
static void add_two_body(StcSpace *s, double gm1, double gm2, double a,
                         double e)
{
    double mu = gm1 + gm2;
    double r = a * (1.0 - e);
    double vrel = sqrt(mu * (1.0 + e) / (a * (1.0 - e)));
    double f1 = gm2 / mu; /* primary sits at -f1 * r_rel */
    double f2 = gm1 / mu;
    double p1[3] = {-r * f1, 0.0, 0.0};
    double v1[3] = {0.0, -vrel * f1, 0.0};
    double p2[3] = {r * f2, 0.0, 0.0};
    double v2[3] = {0.0, vrel * f2, 0.0};
    stc_space_add(s, "prim", STC_SK_STAR, -1, gm1, 0.0, p1, v1);
    stc_space_add(s, "orb", STC_SK_PLANET, 0, gm2, 0.0, p2, v2);
}

static double pos_err_rel(const StcSpace *s, int i,
                          double x0, double y0, double z0, double scale)
{
    double dx = s->x[i] - x0;
    double dy = s->y[i] - y0;
    double dz = s->z[i] - z0;
    return sqrt(dx * dx + dy * dy + dz * dz) / scale;
}

/* --------------------------------------------------------------- */
/* 1: kepler                                                        */

static void check_kepler(void)
{
    StcSpace s;
    double a = STC_AU;
    double gm1 = STC_GM_SUN, gm2 = STC_GM_EARTH;
    double mu = gm1 + gm2;
    double P = two_body_period(a, mu);
    double x0, y0, z0, e0, e1, perr, derr;
    int ok;

    stc_space_init(&s, 4, 0);
    s.pn = 0;
    s.h = STC_DAY;
    add_two_body(&s, gm1, gm2, a, 0.0);
    x0 = s.x[1]; y0 = s.y[1]; z0 = s.z[1];
    e0 = stc_space_energy(&s);
    stc_space_advance(&s, P, 0);
    e1 = stc_space_energy(&s);
    perr = pos_err_rel(&s, 1, x0, y0, z0, STC_AU);
    derr = fabs(e1 - e0) / fabs(e0);
    ok = (perr < 1e-9) && (derr < 1e-12);
    report(ok, "kepler", perr, 0.0, 1e-9);
    info("kepler", "pos_err=%.3e AU  dE/E=%.3e  substeps=%ld",
         perr, derr, s.steps);
    stc_space_free(&s);
}

/* --------------------------------------------------------------- */
/* 2: order                                                         */

static void check_order(void)
{
    static const int div[4] = {50, 100, 200, 400};
    double a = STC_AU;
    double gm1 = STC_GM_SUN, gm2 = STC_GM_EARTH;
    double mu = gm1 + gm2;
    double P = two_body_period(a, mu);
    double err[4];
    double x0, y0, z0;
    int k;
    double orders[3];

    for (k = 0; k < 4; k++) {
        StcSpace s;
        stc_space_init(&s, 4, 0);
        s.pn = 0;
        s.h = P / (double)div[k];
        s.ckpt_max = 0; /* no ckpt splitting: pure step study */
        add_two_body(&s, gm1, gm2, a, 0.6);
        x0 = s.x[1]; y0 = s.y[1]; z0 = s.z[1];
        stc_space_advance(&s, P, 0);
        err[k] = pos_err_rel(&s, 1, x0, y0, z0, a);
        stc_space_free(&s);
    }
    for (k = 0; k < 4; k++) {
        info("order", "h=P/%-4d err=%.6e", div[k], err[k]);
    }
    for (k = 0; k < 3; k++) {
        orders[k] = err[k + 1] > 0.0 ? log(err[k] / err[k + 1]) / log(2.0)
                                     : NAN;
    }
    /* middle pair usually sits between coarse noise and the roundoff
     * floor; report the median of the non-floored estimates. */
    {
        double use = NAN;
        if (isfinite(orders[1])) {
            use = orders[1];
        } else if (isfinite(orders[0])) {
            use = orders[0];
        } else if (isfinite(orders[2])) {
            use = orders[2];
        }
        report(isfinite(use) && use >= 5.5 && use <= 6.5, "order",
               use, 6.0, 0.5);
        info("order", "observed orders: %.3f %.3f %.3f (used %.3f)",
             orders[0], orders[1], orders[2], use);
    }
}

/* --------------------------------------------------------------- */
/* 3: reversible                                                    */

static void check_reversible(void)
{
    StcSpace s;
    StcRng64 rng;
    double gm_c = 1.0e20;
    double rmin = 0.0, rsys;
    double *x0, *y0, *z0;
    double t_fwd, dev, maxdev = 0.0;
    int i, n = 10;

    stc_space_init(&s, 16, 0);
    s.pn = 0;
    s.collisions = 0;
    stc_rng64_seed(&rng, 0x5eed1234ULL);
    {
        double p[3] = {0, 0, 0}, v[3] = {0, 0, 0};
        stc_space_add(&s, "prim", STC_SK_STAR, -1, gm_c, 0.0, p, v);
    }
    for (i = 1; i < n; i++) {
        double r = stc_rng64_uniform_range(&rng, 5e9, 8e10);
        double phi = stc_rng64_uniform_range(&rng, 0.0, 2.0 * STC_PI);
        double inc = stc_rng64_uniform_range(&rng, -0.2, 0.2);
        double vc = sqrt(gm_c / r) *
                    stc_rng64_uniform_range(&rng, 0.95, 1.05);
        double gm = stc_rng64_uniform_range(&rng, 1e10, 1e16);
        double p[3], v[3];
        p[0] = r * cos(phi);
        p[1] = r * sin(phi);
        p[2] = r * inc;
        /* velocity ~ perpendicular to radius in a tilted plane */
        v[0] = -sin(phi) * vc;
        v[1] = cos(phi) * vc;
        v[2] = 0.3 * inc * vc;
        stc_space_add(&s, "sat", STC_SK_PLANET, 0, gm, 0.0, p, v);
        if (i == 1 || r < rmin) {
            rmin = r;
        }
    }
    rsys = 8e10;
    s.h = two_body_period(rmin, gm_c) / 64.0;
    t_fwd = 100.0 * two_body_period(rmin, gm_c);
    x0 = malloc(sizeof(double) * (size_t)s.nm * 6u);
    y0 = x0 + s.nm;
    z0 = y0 + s.nm;
    memcpy(x0, s.x, sizeof(double) * (size_t)s.nm);
    memcpy(y0, s.y, sizeof(double) * (size_t)s.nm);
    memcpy(z0, s.z, sizeof(double) * (size_t)s.nm);
    stc_space_advance(&s, t_fwd, 0);
    stc_space_advance(&s, -t_fwd, 0);
    for (i = 0; i < s.nm; i++) {
        double dx = s.x[i] - x0[i];
        double dy = s.y[i] - y0[i];
        double dz = s.z[i] - z0[i];
        dev = sqrt(dx * dx + dy * dy + dz * dz) / rsys;
        if (dev > maxdev) {
            maxdev = dev;
        }
    }
    free(x0);
    report(maxdev < 1e-9, "reversible", maxdev, 0.0, 1e-9);
    stc_space_free(&s);
}

/* --------------------------------------------------------------- */
/* 4: mercury                                                       */

static double perihelion_longitude(const StcSpace *s, int i, int j)
{
    StcOrbit o;
    if (stc_space_orbit(s, i, j, &o) != 0) {
        return 0.0;
    }
    return o.node + o.argp;
}

/* Integrate Mercury-like orbit for ~100 Julian years sampling the
 * perihelion longitude once per orbit (phase-matched, so periodic
 * oscillations cancel). Returns the advance in arcsec per century. */
static double measure_advance_h(int pn, double h, int years)
{
    StcSpace s;
    double gm_sun = 1.32712440041e20; /* IAU nominal GM_sun, m^3/s^2 */
    double gm_merc = 2.2032e13;
    double a = 0.387098 * STC_AU;
    double e = 0.205630;
    double mu = gm_sun + gm_merc;
    double P = two_body_period(a, mu);
    double arcsec = STC_PI / (180.0 * 3600.0);
    int n_orb = (int)floor((double)years * STC_YEAR / P + 0.5);
    double wprev, total = 0.0, out;
    int k;

    stc_space_init(&s, 4, 0);
    s.pn = pn;
    s.h = h;
    s.ckpt_interval = 0.0; /* no ckpt grid: pure integration */
    s.ckpt_max = 0;
    add_two_body(&s, gm_sun, gm_merc, a, e);
    wprev = perihelion_longitude(&s, 1, 0);
    for (k = 0; k < n_orb; k++) {
        double w, dw;
        stc_space_advance(&s, P, 0);
        w = perihelion_longitude(&s, 1, 0);
        dw = w - wprev;
        while (dw > STC_PI) dw -= 2.0 * STC_PI;
        while (dw < -STC_PI) dw += 2.0 * STC_PI;
        total += dw;
        wprev = w;
    }
    stc_space_free(&s);
    out = total * (100.0 * STC_YEAR / ((double)n_orb * P)) / arcsec;
    return out;
}

static void check_mercury(void)
{
    double gm_sun = 1.32712440041e20;
    double a = 0.387098 * STC_AU;
    double e = 0.205630;
    double P = two_body_period(a, gm_sun + 2.2032e13);
    double arcsec = STC_PI / (180.0 * 3600.0);
    double adv_pn, adv_n, theory;

    adv_pn = measure_advance_h(1, 10800.0, 100);
    adv_n = measure_advance_h(0, 10800.0, 100);
    /* Einstein: 6 pi GM / (c^2 a (1-e^2)) per orbit (oracle) */
    theory = 6.0 * STC_PI * gm_sun / (STC_C * STC_C * a * (1 - e * e)) *
             (100.0 * STC_YEAR / P) / arcsec;
    info("mercury", "measured %.4f\"/cy, newtonian %.4g\"/cy, "
         "formula %.4f\"/cy", adv_pn, adv_n, theory);
    report(adv_pn > 42.98 * 0.985 && adv_pn < 42.98 * 1.015,
           "mercury", adv_pn, 42.98, 42.98 * 0.015);
    report(fabs(adv_n) < 0.05, "mercury_newtonian", adv_n, 0.0, 0.05);
}

/* --------------------------------------------------------------- */
/* 5: merge                                                         */

static void check_merge(void)
{
    StcSpace s;
    double p0[3], p1[3], dp;
    double m0, m1;
    int i, saw = 0;

    stc_space_init(&s, 4, 0);
    s.pn = 0;
    s.h = 3600.0;
    {
        double pa[3] = {-5e9, 0, 0}, va[3] = {1e4, 0, 0};
        double pb[3] = {5e9, 0, 0}, vb[3] = {-1e4, 0, 0};
        stc_space_add(&s, "a", STC_SK_ROCK, -1, 1e15, 1e8, pa, va);
        stc_space_add(&s, "b", STC_SK_ROCK, -1, 2e15, 1e8, pb, vb);
    }
    m0 = (s.gm[0] + s.gm[1]) / STC_G;
    stc_space_momentum(&s, p0);
    stc_space_advance(&s, 10.0 * STC_DAY, 0);
    stc_space_momentum(&s, p1);
    m1 = 0.0;
    {
        int j;
        for (j = 0; j < s.nm; j++) {
            if (s.alive[j]) {
                m1 += s.gm[j];
            }
        }
        m1 /= STC_G;
    }
    for (i = 0; i < s.nev && i < STC_SPACE_EVLOG; i++) {
        if (s.ev[i].kind == STC_SEV_MERGE) {
            saw = 1;
        }
    }
    dp = 0.0;
    for (i = 0; i < 3; i++) {
        double d = fabs(p1[i] - p0[i]);
        if (fabs(p0[i]) > 0) d /= fabs(p0[i]);
        if (d > dp) dp = d;
    }
    report(stc_space_live_massive(&s) == 1 &&
               fabs(m1 - m0) / m0 < 1e-12 && dp < 1e-12 && saw,
           "merge", (double)stc_space_live_massive(&s), 1.0, 0.0);
    info("merge", "survivor gm=%.6g dp=%.3e nev=%d", m1, dp, s.nev);
    stc_space_free(&s);
}

/* --------------------------------------------------------------- */
/* 6: blackhole                                                     */

static void check_blackhole(void)
{
    StcSpace s;
    double gm_bh = 10.0 * 1.32712440041e20;
    double gm_star = 1.32712440041e20;
    double r_star = 7e8;
    double rs_bh = 2.0 * gm_bh / (STC_C * STC_C);
    double r_p = 8e8; /* aim inside r_t ~ 1.5e9, outside capture */
    double r0 = 1e11;
    double vt = sqrt(2.0 * gm_bh * r_p) / r0;
    int i, saw_d = 0, saw_c = 0, tidx;
    int nt0;

    stc_space_init(&s, 8, 1024);
    s.h = 3600.0;
    s.disrupt_tracers = 64;
    {
        double p0[3] = {0, 0, 0}, v0[3] = {0, 0, 0};
        double ps[3] = {r0, 0, 0}, vs[3] = {0, vt, 0};
        double pt[3] = {-5e9, 0, 0}, vtx[3] = {1e5, 0, 0};
        stc_space_add(&s, "bh", STC_SK_BLACKHOLE, -1, gm_bh, rs_bh,
                      p0, v0);
        stc_space_add(&s, "star", STC_SK_STAR, 0, gm_star, r_star,
                      ps, vs);
        tidx = stc_space_add_tracer(&s, pt, vtx);
    }
    nt0 = s.nt;
    stc_space_advance(&s, 30.0 * STC_DAY, 0);
    for (i = 0; i < s.nev && i < STC_SPACE_EVLOG; i++) {
        if (s.ev[i].kind == STC_SEV_DISRUPT) saw_d = 1;
        if (s.ev[i].kind == STC_SEV_CAPTURE && s.ev[i].b == -1 - tidx)
            saw_c = 1;
    }
    report(saw_d && s.nt - nt0 == s.disrupt_tracers,
           "blackhole_disrupt", (double)(s.nt - nt0),
           (double)s.disrupt_tracers, 0.0);
    info("blackhole", "disrupt=%d capture_tracer=%d nev=%d "
         "nt_delta=%d live_t=%d live_m=%d",
         saw_d, saw_c, s.nev, s.nt - nt0,
         stc_space_live_tracers(&s), stc_space_live_massive(&s));
    report(saw_c, "blackhole_capture", saw_c ? 1.0 : 0.0, 1.0, 0.0);
    stc_space_free(&s);
}

/* --------------------------------------------------------------- */
/* 7: seek                                                          */

typedef struct {
    int nm, nt, nev;
    long steps;
    double t;
    unsigned char *blob;
    size_t blob_n;
} Snap;

static size_t snap_bytes(const StcSpace *s)
{
    return (size_t)s->nm * (14 * sizeof(double) + 2 * sizeof(int) +
                          1 + STC_SPACE_NAME) +
           (size_t)s->nt * (6 * sizeof(double) + 1);
}

static void snap_take(const StcSpace *s, Snap *sn)
{
    unsigned char *p;
    size_t nm = (size_t)s->nm, nt = (size_t)s->nt;
    sn->nm = s->nm; sn->nt = s->nt; sn->nev = s->nev;
    sn->steps = s->steps; sn->t = s->t;
    sn->blob_n = snap_bytes(s);
    sn->blob = malloc(sn->blob_n);
    p = sn->blob;
#define PACK(arr, sz) memcpy(p, (arr), (sz)); p += (sz)
    PACK(s->x, nm * sizeof(double));
    PACK(s->y, nm * sizeof(double));
    PACK(s->z, nm * sizeof(double));
    PACK(s->vx, nm * sizeof(double));
    PACK(s->vy, nm * sizeof(double));
    PACK(s->vz, nm * sizeof(double));
    PACK(s->gm, nm * sizeof(double));
    PACK(s->radius, nm * sizeof(double));
    PACK(s->j2, nm * sizeof(double));
    PACK(s->j4, nm * sizeof(double));
    PACK(s->rref, nm * sizeof(double));
    PACK(s->polx, nm * sizeof(double));
    PACK(s->poly, nm * sizeof(double));
    PACK(s->polz, nm * sizeof(double));
    PACK(s->kind, nm * sizeof(int));
    PACK(s->parent, nm * sizeof(int));
    PACK(s->alive, nm);
    PACK(s->name, nm * STC_SPACE_NAME);
    PACK(s->tx, nt * sizeof(double));
    PACK(s->ty, nt * sizeof(double));
    PACK(s->tz, nt * sizeof(double));
    PACK(s->tvx, nt * sizeof(double));
    PACK(s->tvy, nt * sizeof(double));
    PACK(s->tvz, nt * sizeof(double));
    PACK(s->talive, nt);
#undef PACK
}

static int snap_same(const StcSpace *s, const Snap *sn)
{
    unsigned char *p;
    size_t nm = (size_t)sn->nm, nt = (size_t)sn->nt;
    if (s->nm != sn->nm || s->nt != sn->nt || s->nev != sn->nev ||
        s->steps != sn->steps || s->t != sn->t) {
        return 0;
    }
    if (s->nm == 0) {
        return 1;
    }
    p = sn->blob;
#define CMP(arr, sz) if (memcmp(p, (arr), (sz)) != 0) return 0; p += (sz)
    CMP(s->x, nm * sizeof(double));
    CMP(s->y, nm * sizeof(double));
    CMP(s->z, nm * sizeof(double));
    CMP(s->vx, nm * sizeof(double));
    CMP(s->vy, nm * sizeof(double));
    CMP(s->vz, nm * sizeof(double));
    CMP(s->gm, nm * sizeof(double));
    CMP(s->radius, nm * sizeof(double));
    CMP(s->j2, nm * sizeof(double));
    CMP(s->j4, nm * sizeof(double));
    CMP(s->rref, nm * sizeof(double));
    CMP(s->polx, nm * sizeof(double));
    CMP(s->poly, nm * sizeof(double));
    CMP(s->polz, nm * sizeof(double));
    CMP(s->kind, nm * sizeof(int));
    CMP(s->parent, nm * sizeof(int));
    CMP(s->alive, nm);
    CMP(s->name, nm * STC_SPACE_NAME);
    CMP(s->tx, nt * sizeof(double));
    CMP(s->ty, nt * sizeof(double));
    CMP(s->tz, nt * sizeof(double));
    CMP(s->tvx, nt * sizeof(double));
    CMP(s->tvy, nt * sizeof(double));
    CMP(s->tvz, nt * sizeof(double));
    CMP(s->talive, nt);
#undef CMP
    return 1;
}

static void check_seek(void)
{
    StcSpace s;
    Snap s10, s50;
    double gap = 3.46e10; /* contact at ~20.02 d with 2e4 m/s closing */
    int ok_back, ok_fwd, rc;

    stc_space_init(&s, 4, 0);
    s.pn = 0;
    s.h = 3600.0;
    s.ckpt_interval = 5.0 * STC_DAY;
    s.ckpt_max = 16;
    {
        double pa[3] = {-gap * 0.5, 0, 0}, va[3] = {1e4, 0, 0};
        double pb[3] = {gap * 0.5, 0, 0}, vb[3] = {-1e4, 0, 0};
        stc_space_add(&s, "a", STC_SK_ROCK, -1, 1e15, 1e7, pa, va);
        stc_space_add(&s, "b", STC_SK_ROCK, -1, 2e15, 1e7, pb, vb);
    }
    stc_space_advance(&s, 10.0 * STC_DAY, 0);
    snap_take(&s, &s10);
    stc_space_advance(&s, 40.0 * STC_DAY, 0);
    snap_take(&s, &s50);
    /* rewind across the merge to day 10 */
    rc = stc_space_seek(&s, 10.0 * STC_DAY);
    ok_back = (rc == 0) && snap_same(&s, &s10) &&
              stc_space_live_massive(&s) == 2;
    report(ok_back, "seek_back", (double)rc, 0.0, 0.0);
    rc = stc_space_seek(&s, 50.0 * STC_DAY);
    ok_fwd = (rc == 0) && snap_same(&s, &s50) &&
             stc_space_live_massive(&s) == 1;
    report(ok_fwd, "seek_fwd", (double)rc, 0.0, 0.0);
    info("seek", "nev=%d steps=%ld live=%d", s.nev, s.steps,
         stc_space_live_massive(&s));
    free(s10.blob);
    free(s50.blob);
    stc_space_free(&s);
}

/* --------------------------------------------------------------- */
/* 8: j2_node — zonal J2 nodal regression                            */

/* An orbiter at periapsis: apsidal line along azimuth phi in the
 * ecliptic plane, orbit plane tilted by inc about that line. pos/vel
 * are scaled by frac (use gm1/mu-style fractions to place the pair
 * barycentrically; frac = 1 gives the orbiter relative to a fixed
 * centre at the origin). */
static void orbiter_state(double gm_c, double gm, double a, double e,
                          double inc, double phi, double frac,
                          double pos[3], double vel[3])
{
    double mu = gm_c + gm;
    double rp = a * (1.0 - e);
    double vp = sqrt(mu * (1.0 + e) / rp);
    double ci = cos(inc), si = sin(inc);
    pos[0] = frac * rp * cos(phi);
    pos[1] = frac * rp * sin(phi);
    pos[2] = 0.0;
    vel[0] = -frac * vp * sin(phi) * ci;
    vel[1] = frac * vp * cos(phi) * ci;
    vel[2] = frac * vp * si;
}

static void check_j2_node(void)
{
    StcSpace s;
    StcOrbit o;
    double gm_p = 1.26686534e17; /* Jupiter, m^3/s^2 */
    double gm_m = 5.9599e12;     /* Io-like moon */
    double a = 4.217e8;          /* 421700 km */
    double j2 = 0.0147, rref = 7.1492e7;
    double inc = 0.5 * STC_PI / 180.0;
    double mu = gm_p + gm_m;
    double n = sqrt(mu / (a * a * a));
    double P = 2.0 * STC_PI / n;
    double f1 = gm_m / mu, f2 = gm_p / mu;
    double pole[3] = {0.0, 0.0, 1.0};
    int norb = (int)floor(200.0 * STC_DAY / P);
    double prev, total = 0.0, meas, an, relerr;
    int k, ip, im;

    stc_space_init(&s, 8, 0);
    s.pn = 0;
    s.collisions = 0;
    s.h = 3600.0;
    s.ckpt_max = 0;
    {
        double pp[3], vp[3], pm[3], vm[3];
        orbiter_state(gm_p, gm_m, a, 0.0, inc, 0.0, -f1, pp, vp);
        orbiter_state(gm_p, gm_m, a, 0.0, inc, 0.0, f2, pm, vm);
        ip = stc_space_add(&s, "jup", STC_SK_PLANET, -1, gm_p, 0.0,
                           pp, vp);
        im = stc_space_add(&s, "moon", STC_SK_MOON, ip, gm_m, 0.0,
                           pm, vm);
    }
    stc_space_set_oblateness(&s, ip, j2, 0.0, rref, pole);
    stc_space_orbit(&s, im, ip, &o);
    info("j2_node", "a=%.6e m  e=%.4f  inc=%.4f deg  P=%.4f d",
         o.a, o.e, o.inc * 180.0 / STC_PI, o.period / STC_DAY);
    prev = o.node;
    for (k = 0; k < norb; k++) {
        double d;
        stc_space_advance(&s, P, 0);
        stc_space_orbit(&s, im, ip, &o);
        d = o.node - prev;
        while (d > STC_PI) d -= 2.0 * STC_PI;
        while (d < -STC_PI) d += 2.0 * STC_PI;
        total += d;
        prev = o.node;
    }
    meas = total / ((double)norb * P); /* rad/s */
    an = -1.5 * n * j2 * (rref / a) * (rref / a) * cos(inc); /* oracle */
    relerr = fabs(meas - an) / fabs(an);
    info("j2_node", "regression: measured %.6e rad/s, analytic %.6e, "
         "rel err %.4e over %d orbits", meas, an, relerr, norb);
    report(relerr < 0.01, "j2_node", meas, an, 0.01 * fabs(an));
    stc_space_free(&s);
}

/* --------------------------------------------------------------- */
/* 9: j2_momentum — zonal reaction conserves total momentum          */

static void check_j2_momentum(void)
{
    StcSpace s;
    double gm_p = 1.26686534e17;
    double p0[3], p1[3];
    double scale = 0.0, dp = 0.0;
    double pole[3] = {0.03, -0.02, 1.0}; /* tilted; normalised inside */
    int i, ip, im1, im2, it;

    stc_space_init(&s, 8, 8);
    s.pn = 0;
    s.collisions = 0;
    s.h = 3600.0;
    s.ckpt_max = 0;
    {
        double pp[3] = {0.0, 0.0, 0.0}, vp[3] = {0.0, 0.0, 0.0};
        double d50 = 0.5 * STC_PI / 180.0;
        ip = stc_space_add(&s, "jup", STC_SK_PLANET, -1, gm_p, 0.0,
                           pp, vp);
        /* two massive moons: reactions must cancel their zonal pulls */
        {
            double p[3], v[3];
            orbiter_state(gm_p, 5.9599e12, 4.217e8, 0.004, d50, 0.0,
                          1.0, p, v);
            im1 = stc_space_add(&s, "io", STC_SK_MOON, ip, 5.9599e12,
                                0.0, p, v);
            orbiter_state(gm_p, 7.1793e12, 1.8827e9, 0.0074,
                          0.3 * STC_PI / 180.0, 1.2, 1.0, p, v);
            im2 = stc_space_add(&s, "cal", STC_SK_MOON, ip, 7.1793e12,
                                0.0, p, v);
            /* tracer at 6e8 m: feels J2, must NOT get a reaction */
            orbiter_state(gm_p, 0.0, 6.0e8, 0.01,
                          1.0 * STC_PI / 180.0, 2.5, 1.0, p, v);
            it = stc_space_add_tracer(&s, p, v);
        }
    }
    (void)im1; (void)im2; (void)it;
    stc_space_set_oblateness(&s, ip, 0.0147, -5.9e-4, 7.1492e7, pole);
    stc_space_momentum(&s, p0);
    stc_space_advance(&s, 30.0 * STC_DAY, 0);
    stc_space_momentum(&s, p1);
    for (i = 0; i < s.nm; i++) {
        if (s.alive[i]) {
            double vv = sqrt(s.vx[i] * s.vx[i] + s.vy[i] * s.vy[i] +
                             s.vz[i] * s.vz[i]);
            scale += s.gm[i] * vv / STC_G;
        }
    }
    for (i = 0; i < 3; i++) {
        double d = fabs(p1[i] - p0[i]);
        if (d > dp) dp = d;
    }
    dp /= scale;
    info("j2_momentum", "|dp|/sum m|v| = %.4e over 30 d", dp);
    report(dp < 1e-13, "j2_momentum", dp, 0.0, 1e-13);
    stc_space_free(&s);
}

/* --------------------------------------------------------------- */
/* 10: pn_step — implicit-midpoint kick restores step independence   */

static void run_inner_pn(double h, double out[3])
{
    StcSpace s;
    double gm_sun = 1.32712440041e20;
    static const double pl[3][4] = {
        /* gm, a (AU), e, apsidal azimuth phi (rad) */
        {2.2032e13, 0.387098, 0.205630, 0.0},     /* Mercury */
        {3.24858592e14, 0.723332, 0.006772, 2.1}, /* Venus */
        {4.03503e14, 1.000000, 0.016709, 4.0},    /* Earth+Moon */
    };
    int k;
    stc_space_init(&s, 8, 0);
    s.pn = 1;
    s.collisions = 0;
    s.h = h;
    s.ckpt_max = 0;
    {
        double p[3] = {0.0, 0.0, 0.0}, v[3] = {0.0, 0.0, 0.0};
        stc_space_add(&s, "Sun", STC_SK_STAR, -1, gm_sun, 0.0, p, v);
    }
    for (k = 0; k < 3; k++) {
        double p[3], v[3];
        orbiter_state(gm_sun, pl[k][0], pl[k][1] * STC_AU, pl[k][2],
                      0.0, pl[k][3], 1.0, p, v);
        stc_space_add(&s, "pl", STC_SK_PLANET, 0, pl[k][0], 0.0, p, v);
    }
    stc_space_advance(&s, 10.0 * STC_YEAR, 0);
    out[0] = s.x[1];
    out[1] = s.y[1];
    out[2] = s.z[1];
    stc_space_free(&s);
}

static void check_pn_step(void)
{
    double a6[3], a1[3], d;
    run_inner_pn(21600.0, a6);
    run_inner_pn(3600.0, a1);
    d = sqrt((a6[0] - a1[0]) * (a6[0] - a1[0]) +
             (a6[1] - a1[1]) * (a6[1] - a1[1]) +
             (a6[2] - a1[2]) * (a6[2] - a1[2]));
    info("pn_step", "Mercury 10 yr: |x(h=6h) - x(h=1h)| = %.3f km",
         d / 1000.0);
    report(d < 1000.0, "pn_step_independence", d, 0.0, 1000.0);
}

/* --------------------------------------------------------------- */
/* 11: threads                                                      */

static uint64_t fnv1a(const void *buf, size_t n, uint64_t h)
{
    const unsigned char *p = (const unsigned char *)buf;
    size_t i;
    for (i = 0; i < n; i++) {
        h ^= p[i];
        h *= 0x100000001b3ULL;
    }
    return h;
}

static uint64_t hash_state(const StcSpace *s)
{
    uint64_t h = 1469598103934665603ULL;
    size_t nm = (size_t)s->nm, nt = (size_t)s->nt;
    h = fnv1a(&s->nm, sizeof(s->nm), h);
    h = fnv1a(&s->nt, sizeof(s->nt), h);
    h = fnv1a(&s->t, sizeof(s->t), h);
    h = fnv1a(&s->steps, sizeof(s->steps), h);
    h = fnv1a(&s->nev, sizeof(s->nev), h);
    h = fnv1a(s->x, nm * sizeof(double), h);
    h = fnv1a(s->y, nm * sizeof(double), h);
    h = fnv1a(s->z, nm * sizeof(double), h);
    h = fnv1a(s->vx, nm * sizeof(double), h);
    h = fnv1a(s->vy, nm * sizeof(double), h);
    h = fnv1a(s->vz, nm * sizeof(double), h);
    h = fnv1a(s->gm, nm * sizeof(double), h);
    h = fnv1a(s->alive, nm, h);
    h = fnv1a(s->tx, nt * sizeof(double), h);
    h = fnv1a(s->ty, nt * sizeof(double), h);
    h = fnv1a(s->tz, nt * sizeof(double), h);
    h = fnv1a(s->tvx, nt * sizeof(double), h);
    h = fnv1a(s->tvy, nt * sizeof(double), h);
    h = fnv1a(s->tvz, nt * sizeof(double), h);
    h = fnv1a(s->talive, nt, h);
    return h;
}

static uint64_t run_threads_scenario(void)
{
    StcSpace s;
    StcRng64 rng;
    int i;
    stc_space_init(&s, 2048, 512);
    s.pn = 1;
    s.h = 3600.0;
    s.theta = 0.5;
    s.tree_threshold = 512;
    stc_rng64_seed(&rng, 0xdecafULL);
    {
        double p[3] = {0, 0, 0}, v[3] = {0, 0, 0};
        stc_space_add(&s, "prim", STC_SK_STAR, -1, 1e21, 0.0, p, v);
    }
    for (i = 1; i < 2000; i++) {
        double r = stc_rng64_uniform_range(&rng, 1e9, 1e11);
        double phi = stc_rng64_uniform_range(&rng, 0.0, 2.0 * STC_PI);
        double zz = stc_rng64_uniform_range(&rng, -0.05, 0.05);
        double vc = sqrt(1e21 / r);
        double gm = stc_rng64_uniform_range(&rng, 1e10, 1e17);
        double p[3], v[3];
        p[0] = r * cos(phi);
        p[1] = r * sin(phi);
        p[2] = r * zz;
        v[0] = -sin(phi) * vc;
        v[1] = cos(phi) * vc;
        v[2] = 0.1 * zz * vc;
        stc_space_add(&s, "b", STC_SK_PLANET, 0, gm, 0.0, p, v);
    }
    for (i = 0; i < 500; i++) {
        double r = stc_rng64_uniform_range(&rng, 1e9, 1e11);
        double phi = stc_rng64_uniform_range(&rng, 0.0, 2.0 * STC_PI);
        double vc = sqrt(1e21 / r);
        double p[3], v[3];
        p[0] = r * cos(phi);
        p[1] = r * sin(phi);
        p[2] = stc_rng64_uniform_range(&rng, -1e9, 1e9);
        v[0] = -sin(phi) * vc;
        v[1] = cos(phi) * vc;
        v[2] = 0.0;
        stc_space_add_tracer(&s, p, v);
    }
    stc_space_advance(&s, 2.0 * STC_DAY, 0);
    {
        uint64_t h = hash_state(&s);
        stc_space_free(&s);
        return h;
    }
}

static void check_threads(const char *argv0)
{
    uint64_t h_par = run_threads_scenario();
    char cmd[1024];
    char line[128];
    uint64_t h_ser = 0;
    int got = 0;
    FILE *fp;

    snprintf(cmd, sizeof(cmd), "\"%s\" threads-child", argv0);
    fp = popen(cmd, "r");
    if (fp != NULL) {
        if (fgets(line, sizeof(line), fp) != NULL &&
            sscanf(line, "HASH %llx", (unsigned long long *)&h_ser) == 1) {
            got = 1;
        }
        pclose(fp);
    }
    info("threads", "par=%016llx ser=%016llx",
         (unsigned long long)h_par, (unsigned long long)h_ser);
    report(got && h_par == h_ser, "threads", (double)(h_par == h_ser),
           1.0, 0.0);
}

/* --------------------------------------------------------------- */
/* 9: tree accuracy                                                  */

static int dbl_cmp(const void *a, const void *b)
{
    double x = *(const double *)a, y = *(const double *)b;
    return (x > y) - (x < y);
}

static void check_tree(void)
{
    StcSpace s;
    StcRng64 rng;
    int i, n = 5000;
    double *ax, *ay, *az, *bx, *by, *bz, *err;
    double med, p99;
    double a_pl = 3e10;

    stc_space_init(&s, n, 0);
    s.pn = 0;
    s.ckpt_max = 0;
    stc_rng64_seed(&rng, 0x91a7bULL);
    for (i = 0; i < n; i++) {
        double u = stc_rng64_uniform(&rng);
        double rr = a_pl / sqrt(pow(u, -2.0 / 3.0) - 1.0);
        double dx = stc_rng64_normal(&rng);
        double dy = stc_rng64_normal(&rng);
        double dz = stc_rng64_normal(&rng);
        double dn = sqrt(dx * dx + dy * dy + dz * dz);
        double p[3], v[3] = {0, 0, 0};
        if (!(dn > 0)) dn = 1.0;
        p[0] = rr * dx / dn;
        p[1] = rr * dy / dn;
        p[2] = rr * dz / dn;
        stc_space_add(&s, "p", STC_SK_ASTEROID, 0, STC_GM_EARTH, 0.0,
                      p, v);
    }
    ax = malloc(sizeof(double) * (size_t)n * 6u);
    ay = ax + n; az = ay + n;
    bx = az + n; by = bx + n; bz = by + n;
    err = malloc(sizeof(double) * (size_t)n);

    s.tree_threshold = 1000000; /* direct */
    stc_space_accel(&s, ax, ay, az);
    s.tree_threshold = 1;       /* tree, theta 0.5 */
    s.theta = 0.5;
    stc_space_accel(&s, bx, by, bz);
    for (i = 0; i < n; i++) {
        double dx = bx[i] - ax[i], dy = by[i] - ay[i], dz = bz[i] - az[i];
        double da = sqrt(dx * dx + dy * dy + dz * dz);
        double aa = sqrt(ax[i] * ax[i] + ay[i] * ay[i] + az[i] * az[i]);
        err[i] = aa > 0 ? da / aa : da;
    }
    qsort(err, (size_t)n, sizeof(double), dbl_cmp);
    med = err[n / 2];
    p99 = err[(int)(0.99 * (n - 1))];
    info("tree", "median=%.4e p99=%.4e", med, p99);
    report(med < 1e-3 && p99 < 1e-2, "tree", med, 1e-3, 1e-3);
    free(ax);
    free(err);
    stc_space_free(&s);
}

/* --------------------------------------------------------------- */
/* 10: perf                                                         */

static void perf_case(const char *tag, int nm, int nt, int tree, int pn,
                      double h, double days)
{
    StcSpace s;
    StcRng64 rng;
    int i;
    double t0, t1;
    long took;

    stc_space_init(&s, nm, nt);
    s.pn = pn;
    s.h = h;
    s.ckpt_max = 0;
    if (tree) {
        s.tree_threshold = 64;
    } else {
        s.tree_threshold = 1 << 28;
    }
    stc_rng64_seed(&rng, 0xbeefULL);
    {
        double p[3] = {0, 0, 0}, v[3] = {0, 0, 0};
        stc_space_add(&s, "prim", STC_SK_STAR, -1, 1e21, 0.0, p, v);
    }
    for (i = 1; i < nm; i++) {
        double r = stc_rng64_uniform_range(&rng, 1e9, 1e11);
        double phi = stc_rng64_uniform_range(&rng, 0.0, 2.0 * STC_PI);
        double vc = sqrt(1e21 / r);
        double p[3], v[3];
        p[0] = r * cos(phi);
        p[1] = r * sin(phi);
        p[2] = stc_rng64_uniform_range(&rng, -1e9, 1e9);
        v[0] = -sin(phi) * vc;
        v[1] = cos(phi) * vc;
        v[2] = 0.0;
        stc_space_add(&s, "b", STC_SK_PLANET, 0, 1e15, 0.0, p, v);
    }
    for (i = 0; i < nt; i++) {
        double r = stc_rng64_uniform_range(&rng, 1e9, 1e11);
        double phi = stc_rng64_uniform_range(&rng, 0.0, 2.0 * STC_PI);
        double vc = sqrt(1e21 / r);
        double p[3], v[3];
        p[0] = r * cos(phi);
        p[1] = r * sin(phi);
        p[2] = stc_rng64_uniform_range(&rng, -1e9, 1e9);
        v[0] = -sin(phi) * vc;
        v[1] = cos(phi) * vc;
        v[2] = 0.0;
        stc_space_add_tracer(&s, p, v);
    }
    t0 = now_s();
    took = stc_space_advance(&s, days * STC_DAY, 0);
    t1 = now_s();
    info("perf",
         "%s nm=%d nt=%d tree=%d pn=%d: %.2f ms/substep (%ld in %.3f s)",
         tag, nm, nt, tree, pn,
         1e3 * (t1 - t0) / (double)(took > 0 ? took : 1),
         took, t1 - t0);
    stc_space_free(&s);
}

/* --------------------------------------------------------------- */

int main(int argc, char **argv)
{
    if (argc > 1 && strcmp(argv[1], "threads-child") == 0) {
        /* forced-serial child: setenv BEFORE the first pool use */
        setenv("STC_PAR_THREADS", "1", 1);
        printf("HASH %016llx\n",
               (unsigned long long)run_threads_scenario());
        return 0;
    }
    printf("stspace: SPACE-CORE verification\n");
    check_kepler();
    check_order();
    check_reversible();
    check_mercury();
    check_merge();
    check_blackhole();
    check_seek();
    check_j2_node();
    check_j2_momentum();
    check_pn_step();
    check_threads(argv[0]);
    check_tree();
    perf_case("small+tracers", 45, 20000, 1, 0, 3600.0, 0.05);
    perf_case("direct2000", 2000, 0, 0, 0, 3600.0, 0.05);
    perf_case("tree20000", 20000, 0, 1, 0, 3600.0, 3600.0 / STC_DAY);
    perf_case("bigtracers", 20, 1000000, 0, 0, 3600.0, 3600.0 / STC_DAY);
    perf_case("bigtracers-pn", 20, 1000000, 0, 1, 3600.0,
              3600.0 / STC_DAY);
    printf("RESULT: %s (%d failures)\n", g_fails ? "FAIL" : "PASS",
           g_fails);
    return g_fails ? 1 : 0;
}
