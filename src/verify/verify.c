/* SPACETIME LAB — verify (builder 12): classic-GR acceptance checks.
 *
 * Own file: src/verify/verify.c -> build/stverify (own main, no library).
 * May include all sibling stc headers; includes shipped headers (core +
 * sr/geodesic/clocks/lensing/kerr/waves + body/rk4) plus C stdlib.
 * Gate-0 R1 repair: every check's PASS/FAIL consumes the shipped stc_*
 * API. The independent reference implementations of the CONTRACT orbit
 * equations (photon u''+u = 1.5*rs*u^2, timelike u''+u = GM/L^2 +
 * 1.5*rs*u^2, Lorentz gamma/interval) are kept as `oracle_*` second
 * opinions; their agreement with the shipped calls is asserted at the
 * stated tolerance and printed in the detail line.
 *
 * Equation index (one per check):
 *   sr_invariance    Lorentz factor gamma = 1/sqrt(1-v^2/c^2) (stc_sr_*);
 *                    interval s^2 = (c dt)^2 - dx^2 Lorentz invariance
 *                    (stc_sr_interval via stc_sr_boost); velocity-addition
 *                    round-trip (stc_sr_add_velocity); twin tau = T/gamma
 *                    (stc_sr_proper_time).
 *   photon_bending   Einstein deflection alpha = 4GM/(c^2 b)
 *                    (stc_lensing_deflection); geodesic leg through
 *                    shipped stc_geodesic_photon_init + stc_geodesic_
 *                    integrate/status; critical impact b_crit =
 *                    3*sqrt(3)*M_geom.
 *   mercury_precession Einstein shift dphi = 6*pi*GM/(c^2 a (1-e^2))
 *                    (stc_clocks_perihelion_shift / _mercury_arcsec_per_
 *                    century); timelike orbit leg through shipped
 *                    stc_geodesic_orbit_init + integrate.
 *   gps_offset       SR transverse Doppler dF = (vg^2-vs^2)/(2c^2);
 *                    Schwarzschild ratio dF = gsat/ggnd - 1;
 *                    circular velocity v = sqrt(GM/r).
 *   shapiro_delay    Grazing dt = (K GM/c^3) ln(4 r1 r2/b^2), K = 2/4;
 *                    PPN triangle dt = (2GM/c^3) ln((s+d)/(s-d)).
 *   gw_chirp         Quadrupole h = 4 G^2 m1 m2/(c^4 a R);
 *                    chirp f(tau) = (1/pi)(5/256tau)^(3/8)(GMc/c^3)^(-5/8);
 *                    inverse tau(f) = (5/256)(GMc/c^3)^(-5/3)(pi f)^(-8/3).
 *   kerr_isco        Bardeen-Press-Teukolsky ISCO (6M / 1M / 9M limits).
 *   nbody_orbit      Verlet circular orbit + energy conservation.
 *   momentum_merge   Verlet step + poststep merge conserves momentum.
 *   rk4_poststep     RK4 step + stc_sim_poststep merge: same policy as
 *                    the Verlet path (momentum-conserving merge, n
 *                    compaction).
 *   bbh_timeline     TaylorT4 3.5PN inspiral (arXiv:0907.0700 Eq. 6)
 *                    internal consistency: dphi/dt == omega = v^3,
 *                    ISCO boundary values, seam continuity, remnant
 *                    fits, determinism (stc_bbh_init/eval).
 *   bbh_t2_crosscheck shipped T4 t_m vs TaylorT2 t(v) closed form
 *                    (arXiv:0907.0700 Eq. 8b oracle), 35 Hz -> ISCO,
 *                    q = 1.24, 65 Msun.
 *   bbh_gw150914_context  INFO only: production 35 Hz -> t_peak
 *                    duration and GW cycle count next to the LIGO
 *                    ~0.2 s / ~8 cycle quote (arXiv:1602.03837).
 *
 * SI doubles throughout; fixed RK4 step counts (documented below); no
 * malloc, no RNG in physics (one fixed-seed determinism probe in gw_chirp).
 */

#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>

#include "stc/body.h"
#include "stc/clocks.h"
#include "stc/geodesic.h"
#include "stc/kerr.h"
#include "stc/lensing.h"
#include "stc/merger.h"
#include "stc/rk4.h"
#include "stc/rng.h"
#include "stc/sr.h"
#include "stc/units.h"
#include "stc/waves.h"

/* Mars semi-major axis in AU (circular approx, cf. Venus in clocks.h). */
#define STC_VERIFY_MARS_AU 1.52371

/* Fixed RK4 step counts (determinism: no adaptive stepping). */
#define STC_VERIFY_PHOTON_STEPS 40000  /* oracle photon integrator */
#define STC_VERIFY_MERCURY_STEPS 40000 /* oracle + shipped orbit legs */
#define STC_VERIFY_FATE_STEPS 60000    /* photon fate probes, 6*pi sweep */
#define STC_VERIFY_BEND_STEPS 80000    /* shipped bending leg */
#define STC_VERIFY_BEND_DPHI 8.0e-5    /* shipped bending leg step, rad */

/* Shared post-step radius/azimuth histories for the shipped
 * stc_geodesic_integrate legs (sized to the largest fixed step count;
 * static, no malloc). */
#define STC_VERIFY_GEO_HIST STC_VERIFY_BEND_STEPS
static double stc_verify_rhist[STC_VERIFY_GEO_HIST];
static double stc_verify_phist[STC_VERIFY_GEO_HIST];

/* Relative difference |a-b|/|b| (b != 0 at all call sites). */
static double stc_verify_reldiff(double a, double b)
{
    return fabs(a - b) / fabs(b);
}

/* ORACLE (R1 second opinion, not the checked path): Lorentz factor
 * gamma = 1/sqrt(1-v^2/c^2). v in m/s. */
static double oracle_gamma(double v_m_s)
{
    double b = v_m_s / STC_C;

    return 1.0 / sqrt(1.0 - b * b);
}

/* ORACLE (R1 second opinion, not the checked path): photon orbit
 * equation u'' + u = 1.5*rs*u^2, fixed-step RK4. Incoming photon with
 * impact parameter b: u(0) = 0, u'(0) = 1/b. Returns total deflection
 * (outgoing asymptote minus pi) in radians, or NaN if no outgoing
 * asymptote is found in range. */
static double oracle_photon_deflection(double rs_m, double b_m)
{
    const double phi_max = STC_PI + 0.02;
    const double h = phi_max / (double)STC_VERIFY_PHOTON_STEPS;
    double u = 0.0;
    double w = 1.0 / b_m;
    double phi = 0.0;
    double prev_u = 0.0;
    double prev_phi = 0.0;
    int i;

    for (i = 0; i < STC_VERIFY_PHOTON_STEPS; i++) {
        double k1u = w;
        double k1w = -u + 1.5 * rs_m * u * u;
        double u2 = u + 0.5 * h * k1u;
        double w2 = w + 0.5 * h * k1w;
        double k2u = w2;
        double k2w = -u2 + 1.5 * rs_m * u2 * u2;
        double u3 = u + 0.5 * h * k2u;
        double w3 = w + 0.5 * h * k2w;
        double k3u = w3;
        double k3w = -u3 + 1.5 * rs_m * u3 * u3;
        double u4 = u + h * k3u;
        double w4 = w + h * k3w;
        double k4u = w4;
        double k4w = -u4 + 1.5 * rs_m * u4 * u4;

        prev_u = u;
        prev_phi = phi;
        u += h * (k1u + 2.0 * k2u + 2.0 * k3u + k4u) / 6.0;
        w += h * (k1w + 2.0 * k2w + 2.0 * k3w + k4w) / 6.0;
        phi += h;
        if (phi > 1.0 && u < 0.0) {
            double f = prev_u / (prev_u - u);

            return (prev_phi + f * (phi - prev_phi)) - STC_PI;
        }
    }
    return NAN;
}

/* ORACLE (R1 second opinion): photon fate through the same orbit
 * equation: 1 = escaped (u returns to 0), 2 = captured (r < rs),
 * 0 = neither within 6*pi (trapped/winding). rmin_out gets the smallest
 * positive r seen (meters, HUGE if none). */
static int oracle_photon_fate(double rs_m, double b_m, double *rmin_out)
{
    const double phi_max = 6.0 * STC_PI;
    const double h = phi_max / (double)STC_VERIFY_FATE_STEPS;
    double u = 0.0;
    double w = 1.0 / b_m;
    double phi = 0.0;
    double rmin = HUGE_VAL;
    int i;

    for (i = 0; i < STC_VERIFY_FATE_STEPS; i++) {
        double k1u = w;
        double k1w = -u + 1.5 * rs_m * u * u;
        double u2 = u + 0.5 * h * k1u;
        double w2 = w + 0.5 * h * k1w;
        double k2u = w2;
        double k2w = -u2 + 1.5 * rs_m * u2 * u2;
        double u3 = u + 0.5 * h * k2u;
        double w3 = w + 0.5 * h * k2w;
        double k3u = w3;
        double k3w = -u3 + 1.5 * rs_m * u3 * u3;
        double u4 = u + h * k3u;
        double w4 = w + h * k3w;
        double k4u = w4;
        double k4w = -u4 + 1.5 * rs_m * u4 * u4;

        u += h * (k1u + 2.0 * k2u + 2.0 * k3u + k4u) / 6.0;
        w += h * (k1w + 2.0 * k2w + 2.0 * k3w + k4w) / 6.0;
        phi += h;
        if (u > 0.0) {
            double r = 1.0 / u;

            if (r < rmin) {
                rmin = r;
            }
        }
        if (u > 1.0 / rs_m) {
            *rmin_out = rmin;
            return 2;
        }
        if (phi > 0.5 && u < 0.0) {
            *rmin_out = rmin;
            return 1;
        }
    }
    *rmin_out = rmin;
    return 0;
}

/* ORACLE (R1 second opinion): timelike orbit equation in dimensionless
 * x = u*p, where p = a(1-e^2): x'' + x = 1 + K*x^2 with K = 1.5*rs/p,
 * fixed-step RK4. Starts at perihelion (x = 1+e, x' = 0); returns the
 * next perihelion's advance past 2*pi in radians (NaN if not found). */
static double oracle_mercury_dphi(double kk, double ecc)
{
    const double phi_max = 2.0 * STC_PI + 0.002;
    const double h = phi_max / (double)STC_VERIFY_MERCURY_STEPS;
    double x = 1.0 + ecc;
    double w = 0.0;
    double phi = 0.0;
    double px2 = 0.0;
    double px1 = 0.0;
    double pp1 = 0.0;
    int i;

    for (i = 0; i < STC_VERIFY_MERCURY_STEPS; i++) {
        double k1x = w;
        double k1w = 1.0 - x + kk * x * x;
        double x2 = x + 0.5 * h * k1x;
        double w2 = w + 0.5 * h * k1w;
        double k2x = w2;
        double k2w = 1.0 - x2 + kk * x2 * x2;
        double x3 = x + 0.5 * h * k2x;
        double w3 = w + 0.5 * h * k2w;
        double k3x = w3;
        double k3w = 1.0 - x3 + kk * x3 * x3;
        double x4 = x + h * k3x;
        double w4 = w + h * k3w;
        double k4x = w4;
        double k4w = 1.0 - x4 + kk * x4 * x4;

        px2 = px1;
        px1 = x;
        pp1 = phi;
        x += h * (k1x + 2.0 * k2x + 2.0 * k3x + k4x) / 6.0;
        w += h * (k1w + 2.0 * k2w + 2.0 * k3w + k4w) / 6.0;
        phi += h;
        if (phi > STC_PI && px1 > px2 && px1 > x) {
            /* Parabolic vertex refine around pp1. */
            double denom = px2 - 2.0 * px1 + x;
            double off = 0.0;

            if (denom != 0.0) {
                off = 0.5 * h * (px2 - x) / denom;
            }
            return (pp1 + off) - 2.0 * STC_PI;
        }
    }
    return NAN;
}

/* SHIPPED leg (the checked path): photon deflection through the public
 * geodesic API — stc_geodesic_photon_init inbound at r_start, fixed-step
 * stc_geodesic_integrate with the post-step radius history, escape
 * classified by the returned StcGeoStatus. The u = 1/r_start crossing is
 * linearly interpolated between the last two written samples (the same
 * scheme stc_geodesic_photon_deflection uses internally), and the
 * undeflected finite-start sweep 2*acos(b/r_start) is subtracted.
 * Returns total deflection in radians, NaN on failure. */
static double shipped_photon_bend(double rs_m, double b_m, double r_start_m)
{
    StcGeodesic g;
    StcGeoStatus st;
    double u_esc, u_prev, u_now, frac, swept;
    int n;

    if (stc_geodesic_photon_init(&g, rs_m, b_m, r_start_m, 0.0) != 0) {
        return NAN;
    }
    st = stc_geodesic_integrate(&g, STC_VERIFY_BEND_DPHI,
                              STC_VERIFY_BEND_STEPS, r_start_m,
                              stc_verify_rhist, stc_verify_phist);
    if (st != STC_GEO_ESCAPED) {
        return NAN;
    }
    /* integrate advanced exactly n = (g.phi - phi_start)/dphi steps. */
    n = (int)(g.phi / STC_VERIFY_BEND_DPHI + 0.5);
    if (n < 2 || n > STC_VERIFY_BEND_STEPS) {
        return NAN;
    }
    u_esc = 1.0 / r_start_m;
    u_prev = 1.0 / stc_verify_rhist[n - 2]; /* last sub-escape sample */
    u_now = 1.0 / stc_verify_rhist[n - 1];  /* escaping sample (~0) */
    frac = (u_now != u_prev) ? (u_prev - u_esc) / (u_prev - u_now) : 1.0;
    if (frac < 0.0) {
        frac = 0.0;
    }
    if (frac > 1.0) {
        frac = 1.0;
    }
    /* Accumulated per-sample phi from phi_hist (matches the shipped
     * one-shot deflection's interpolation to ~machine precision). */
    swept = stc_verify_phist[n - 2] +
        frac * (stc_verify_phist[n - 1] - stc_verify_phist[n - 2]);
    return swept - 2.0 * acos(b_m / r_start_m);
}

/* SHIPPED leg (the checked path): photon fate via stc_geodesic_photon_
 * init + stc_geodesic_integrate terminal StcGeoStatus — maps to the same
 * convention as the oracle: 2 = captured, 1 = escaped, 0 = still
 * orbiting at the fixed-step cap, -1 = bad input. rmin_out gets the
 * smallest post-step radius written by integrate (meters). */
static int shipped_photon_fate(double rs_m, double b_m, double r_start_m,
                               double *rmin_out)
{
    StcGeodesic g;
    StcGeoStatus st;
    const double dphi = 6.0 * STC_PI / (double)STC_VERIFY_FATE_STEPS;
    double rmin = HUGE_VAL;
    int i, n;

    if (stc_geodesic_photon_init(&g, rs_m, b_m, r_start_m, 0.0) != 0) {
        *rmin_out = rmin;
        return -1;
    }
    st = stc_geodesic_integrate(&g, dphi, STC_VERIFY_FATE_STEPS, 0.0,
                              stc_verify_rhist, NULL);
    /* The terminal status must agree with an explicit status query. */
    if (stc_geodesic_status(&g, 0.0) != st) {
        *rmin_out = rmin;
        return -1;
    }
    n = (int)(g.phi / dphi + 0.5); /* samples actually written */
    if (n > STC_VERIFY_FATE_STEPS) {
        n = STC_VERIFY_FATE_STEPS;
    }
    for (i = 0; i < n; i++) {
        double r = stc_verify_rhist[i];

        if (r > 0.0 && r < rmin) {
            rmin = r;
        }
    }
    *rmin_out = rmin;
    if (st == STC_GEO_CAPTURED) {
        return 2;
    }
    if (st == STC_GEO_ESCAPED) {
        return 1;
    }
    return 0;
}

/* SHIPPED leg (the checked path): timelike orbit through stc_geodesic_
 * orbit_init (periapsis rp = a(1-e), ecc) + fixed-step integrate; the
 * next periapsis is the smallest post-step radius past phi = pi, with a
 * parabolic vertex refine on u = 1/r (same refine the oracle applies to
 * its own samples). Returns the perihelion advance past 2*pi in
 * radians, NaN if not found. */
static double shipped_mercury_dphi(double gm_si, double rs_m, double a_m,
                                   double ecc)
{
    StcGeodesic g;
    const double span = 2.0 * STC_PI + 0.002;
    const double dphi = span / (double)STC_VERIFY_MERCURY_STEPS;
    double r_min = HUGE_VAL;
    double u_prev, u_max, u_next, denom, off;
    int i, k = -1, n;

    if (stc_geodesic_orbit_init(&g, gm_si, rs_m, a_m * (1.0 - ecc), ecc,
                              0.0) != 0) {
        return NAN;
    }
    stc_geodesic_integrate(&g, dphi, STC_VERIFY_MERCURY_STEPS, 0.0,
                           stc_verify_rhist, NULL);
    n = (int)(g.phi / dphi + 0.5); /* samples actually written */
    if (n > STC_VERIFY_MERCURY_STEPS) {
        n = STC_VERIFY_MERCURY_STEPS;
    }
    for (i = 0; i < n; i++) {
        double r = stc_verify_rhist[i];

        if ((double)(i + 1) * dphi > STC_PI && r > 0.0 && r < r_min) {
            r_min = r;
            k = i;
        }
    }
    if (k < 1 || k + 1 >= n) {
        return NAN;
    }
    /* Parabolic vertex on u = 1/r through samples k-1, k, k+1; sample i
     * sits at phi = (i+1)*dphi. */
    u_prev = 1.0 / stc_verify_rhist[k - 1];
    u_max = 1.0 / stc_verify_rhist[k];
    u_next = 1.0 / stc_verify_rhist[k + 1];
    denom = u_prev - 2.0 * u_max + u_next;
    off = (denom != 0.0) ? 0.5 * dphi * (u_prev - u_next) / denom : 0.0;
    return ((double)(k + 1) * dphi + off) - 2.0 * STC_PI;
}

/* Check 1. Equations: gamma = 1/sqrt(1-v^2/c^2); s^2 = (c dt)^2 - dx^2
 * Lorentz invariance; velocity addition round-trip; twin tau = T/gamma.
 * Every measured quantity below comes from the shipped stc_sr_* API. */
static int check_sr_invariance(void)
{
    const double v = 0.6 * STC_C;
    /* Shipped API: gamma, boost, interval, velocity add, proper time. */
    const double g = stc_sr_gamma(v);
    const double g08 = stc_sr_gamma(0.8 * STC_C);
    /* Interval endpoints: A at origin, B at (2 s, 4e8 m). Boost B by
     * v = 0.6c through the shipped Lorentz transform; s^2 must hold. */
    const StcSrEvent ev_a = { 0.0, 0.0 };
    const StcSrEvent ev_b = { 2.0, 4.0e8 };
    const StcSrEvent ev_bp = stc_sr_boost(ev_b, v);
    const double s2 = stc_sr_interval_events(ev_b, ev_a);
    const double s2_p = stc_sr_interval_events(ev_bp, ev_a);
    const double inv_err = stc_verify_reldiff(s2_p, s2);
    /* Velocity-addition round-trip: w = u(+)v, then u' = w(+)(-v). */
    const double u = 0.5 * STC_C;
    const double w = stc_sr_add_velocity(u, v);
    const double u_rt = stc_sr_add_velocity(w, -v);
    const double vadd_err = stc_verify_reldiff(u_rt, u);
    /* Twin: 10 s Earth frame at 0.6c out-and-back, instant turnaround. */
    const double t_earth = 10.0;
    const double t_twin = stc_sr_proper_time(t_earth, v);
    /* Oracle second opinion on gamma (independent re-implementation). */
    const double o_err = stc_verify_reldiff(oracle_gamma(v), g);
    int ok = 1;

    ok = ok && isfinite(w) && isfinite(u_rt) && isfinite(s2_p);
    ok = ok && (stc_verify_reldiff(g, 1.25) < 1.0e-12);
    ok = ok && (stc_verify_reldiff(g08, 5.0 / 3.0) < 1.0e-12);
    ok = ok && (stc_sr_gamma(0.0) == 1.0);
    ok = ok && (inv_err < 1.0e-9);
    ok = ok && (vadd_err < 1.0e-12);
    ok = ok && (t_twin < t_earth);
    ok = ok && (stc_verify_reldiff(t_twin, 8.0) < 1.0e-12);
    ok = ok && (o_err < 1.0e-12); /* oracle agrees with shipped gamma */
    printf("%s sr_invariance: computed=gamma(0.6c) %.12f interval-relerr %.3e "
           "vadd-relerr %.3e twin %.6f/%.6fs oracle-gamma-relerr %.3e "
           "expected=1.25 0.0 0.0 8.0/10.0s tol=1e-12 1e-9 1e-12 1e-12\n",
           ok ? "PASS" : "FAIL", g, inv_err, vadd_err, t_twin, t_earth,
           o_err);
    return ok;
}

/* Check 2. Equations: alpha = 4GM/(c^2 b); photon u''+u = 1.5*rs*u^2;
 * b_crit = 3*sqrt(3)*M_geom (photon sphere r = 1.5*rs). The geodesic
 * leg and the fate probes run through the shipped stc_geodesic_* API;
 * oracle_photon_* re-implementations give second opinions. */
static int check_photon_bending(void)
{
    const double rs = stc_schwarzschild_radius(STC_M_SUN);
    const double arcsec = stc_arcsec();
    const double form = 4.0 * STC_G * STC_M_SUN / (STC_C * STC_C * STC_R_SUN);
    const double lib = stc_lensing_deflection(STC_M_SUN, STC_R_SUN);
    /* Shipped geodesic legs: init + integrate + status, inbound from
     * r_start = 1e12 m; also the shipped one-shot deflection wrapper. */
    const double r_start = 1.0e12;
    const double geo = shipped_photon_bend(rs, STC_R_SUN, r_start);
    const double geo1 = stc_geodesic_photon_deflection(
        rs, STC_R_SUN, r_start, STC_VERIFY_BEND_DPHI,
        STC_VERIFY_BEND_STEPS);
    /* Photon-sphere probe around an abstract M_geom = 1000 m hole. */
    const double mm = 1000.0;
    const double rs_t = 2.0 * mm;
    const double bcrit = 3.0 * sqrt(3.0) * mm;
    const double fate_start = 1.0e6; /* > max(rs_t, b) for init */
    double rmin_sub = 0.0;
    double rmin_sup = 0.0;
    double rmin_crit = 0.0;
    double ormin_sub = 0.0;
    double ormin_sup = 0.0;
    double ormin_crit = 0.0;
    const int fate_sub =
        shipped_photon_fate(rs_t, 0.9 * bcrit, fate_start, &rmin_sub);
    const int fate_sup =
        shipped_photon_fate(rs_t, 1.1 * bcrit, fate_start, &rmin_sup);
    const int fate_crit =
        shipped_photon_fate(rs_t, bcrit, fate_start, &rmin_crit);
    /* Oracle second opinions (independent re-implementations). */
    const double ora = oracle_photon_deflection(rs, STC_R_SUN);
    const int ofate_sub = oracle_photon_fate(rs_t, 0.9 * bcrit, &ormin_sub);
    const int ofate_sup = oracle_photon_fate(rs_t, 1.1 * bcrit, &ormin_sup);
    const int ofate_crit = oracle_photon_fate(rs_t, bcrit, &ormin_crit);
    const double ora_err = stc_verify_reldiff(ora, geo);
    int ok = 1;

    (void)rmin_sub;
    (void)rmin_sup;
    (void)ormin_sub;
    (void)ormin_sup;
    ok = ok && isfinite(lib) && isfinite(geo) && isfinite(geo1);
    ok = ok && (stc_verify_reldiff(lib, form) < 1.0e-12);
    ok = ok && (stc_verify_reldiff(lib / arcsec, 1.75) < 0.01);
    ok = ok && (stc_verify_reldiff(geo / arcsec, 1.75) < 0.01);
    ok = ok && (stc_verify_reldiff(geo, form) < 0.01);
    ok = ok && (stc_verify_reldiff(geo1, geo) < 1.0e-9); /* shipped paths */
    ok = ok && (fate_sub == 2);
    ok = ok && (fate_sup == 1);
    ok = ok && (fate_crit == 0);
    ok = ok && (stc_verify_reldiff(rmin_crit, 1.5 * rs_t) < 0.01);
    /* Oracle agreement (stated tolerances): deflection 1e-3 rel, fate
     * codes exact, critical rmin 5% rel. */
    ok = ok && isfinite(ora) && (ora_err < 1.0e-3);
    ok = ok && (ofate_sub == fate_sub) && (ofate_sup == fate_sup) &&
         (ofate_crit == fate_crit);
    ok = ok && (stc_verify_reldiff(ormin_crit, rmin_crit) < 0.05);
    printf("%s photon_bending: computed=lens %.6f\" geo %.6f\" geo1 %.6f\" "
           "sphere %d/%d/%d oracle %.6f\" %d/%d/%d oracle-relerr %.3e "
           "expected=1.75\" 1.75\" 2/1/0 tol=1%% 1%% exact\n",
           ok ? "PASS" : "FAIL", lib / arcsec, geo / arcsec, geo1 / arcsec,
           fate_sub, fate_sup, fate_crit, ora / arcsec, ofate_sub,
           ofate_sup, ofate_crit, ora_err);
    return ok;
}

/* Check 3. Equations: dphi = 6*pi*GM/(c^2 a (1-e^2)); Kepler-III period;
 * timelike u''+u = GM/L^2 + 1.5*rs*u^2. The pass/fail consumes the
 * shipped clocks API plus a shipped timelike-geodesic leg (orbit_init +
 * integrate); oracle_mercury_dphi is the second-opinion oracle. */
static int check_mercury_precession(void)
{
    const double gm = STC_G * STC_M_SUN;
    const double ma = STC_CLOCKS_MERCURY_A;
    const double me = STC_CLOCKS_MERCURY_ECC;
    const double rs = stc_schwarzschild_radius(STC_M_SUN);
    const double arcsec = stc_arcsec();
    const double per_orbit = 6.0 * STC_PI * gm /
        (STC_C * STC_C * ma * (1.0 - me * me));
    const double period = 2.0 * STC_PI * sqrt(ma * ma * ma / gm);
    const double orbits = STC_CLOCKS_JULIAN_CENTURY_D * STC_DAY / period;
    const double form_cy = per_orbit * orbits / arcsec;
    const double lib_orbit = stc_clocks_perihelion_shift(gm, ma, me);
    const double lib_cy = stc_clocks_mercury_arcsec_per_century();
    const double pp = ma * (1.0 - me * me);
    /* Shipped geodesic leg: perihelion advance of the timelike orbit. */
    const double geo_orbit = shipped_mercury_dphi(gm, rs, ma, me);
    const double geo_cy = geo_orbit * orbits / arcsec;
    /* Oracle second opinion (independent x = u*p integrator). */
    const double ora_orbit = oracle_mercury_dphi(1.5 * rs / pp, me);
    const double ora_err = stc_verify_reldiff(ora_orbit, geo_orbit);
    int ok = 1;

    ok = ok && isfinite(lib_cy) && isfinite(geo_cy) && isfinite(ora_orbit);
    ok = ok && (stc_verify_reldiff(lib_orbit, per_orbit) < 1.0e-9);
    ok = ok && (stc_verify_reldiff(lib_cy, form_cy) < 1.0e-9);
    ok = ok && (stc_verify_reldiff(form_cy, 43.0) < 0.05);
    ok = ok && (stc_verify_reldiff(lib_cy, 43.0) < 0.05);
    ok = ok && (stc_verify_reldiff(geo_cy, 43.0) < 0.05);
    ok = ok && (stc_verify_reldiff(geo_orbit, per_orbit) < 1.0e-3);
    ok = ok && (ora_err < 1.0e-3); /* oracle agrees with shipped geo */
    printf("%s mercury_precession: computed=clocks %.4f\"/cy geo %.4f\"/cy "
           "oracle-relerr %.3e expected=43\"/cy tol=5%% 1e-3\n",
           ok ? "PASS" : "FAIL", lib_cy, geo_cy, ora_err);
    return ok;
}

/* Check 4. Equations: v = sqrt(GM/r); SR dF = (vg^2-vs^2)/(2c^2);
 * GR dF = sqrt(1-rs/rsat)/sqrt(1-rs/rgnd) - 1. All offsets COMPUTED. */
static int check_gps_offset(void)
{
    const double gm = STC_G * STC_M_EARTH;
    const double r_sat = STC_CLOCKS_GPS_RADIUS;
    const double r_gnd = STC_CLOCKS_R_EARTH;
    const double oracle_v_sat = sqrt(gm / r_sat);
    const double oracle_v_gnd = STC_CLOCKS_OMEGA_EARTH * r_gnd;
    const double c2 = STC_C * STC_C;
    const double oracle_sr = (oracle_v_gnd * oracle_v_gnd -
        oracle_v_sat * oracle_v_sat) / (2.0 * c2) * STC_DAY * 1.0e6;
    const double oracle_rs = 2.0 * gm / c2;
    const double oracle_gr =
        (sqrt(1.0 - oracle_rs / r_sat) / sqrt(1.0 - oracle_rs / r_gnd) -
         1.0) * STC_DAY * 1.0e6;
    const double oracle_net = oracle_sr + oracle_gr;
    const double lib_sr = stc_clocks_gps_sr_us_per_day();
    const double lib_gr = stc_clocks_gps_gr_us_per_day();
    const double lib_net = stc_clocks_gps_net_us_per_day();
    int ok = 1;

    ok = ok && (stc_verify_reldiff(lib_sr, oracle_sr) < 1.0e-9);
    ok = ok && (stc_verify_reldiff(lib_gr, oracle_gr) < 1.0e-9);
    ok = ok && (stc_verify_reldiff(lib_net, oracle_net) < 1.0e-9);
    ok = ok && (stc_verify_reldiff(oracle_sr, -7.2) < 0.05);
    ok = ok && (stc_verify_reldiff(oracle_gr, 45.7) < 0.05);
    ok = ok && (stc_verify_reldiff(oracle_net, 38.0) < 0.05);
    ok = ok && (fabs(oracle_net - 38.0) < 1.0);
    printf("%s gps_offset: computed=SR %.4f GR %.4f net %.4f us/day "
           "(lib SR %.4f GR %.4f net %.4f) expected=-7.2 +45.7 +38 us/day "
           "tol=5%% + net 1.0us\n",
           ok ? "PASS" : "FAIL", oracle_sr, oracle_gr, oracle_net,
           lib_sr, lib_gr, lib_net);
    return ok;
}

/* Check 5. Equations: grazing dt = (K GM/c^3) ln(4 r1 r2/b^2), K = 2 one-way
 * and 4 round-trip; PPN triangle dt = (2GM/c^3) ln((s+d)/(s-d)).
 * NOTE: ~120 us is the one-way Earth-Venus value; round-trip is twice that. */
static int check_shapiro_delay(void)
{
    const double gm = STC_G * STC_M_SUN;
    const double c3 = STC_C * STC_C * STC_C;
    const double r1 = STC_AU;
    const double r2 = STC_CLOCKS_VENUS_R;
    const double bb = STC_R_SUN;
    const double oracle_one_way =
        (2.0 * gm / c3) * log(4.0 * r1 * r2 / (bb * bb));
    const double lib_one = stc_clocks_shapiro_grazing(gm, r1, r2, bb, 0);
    const double lib_rt = stc_clocks_shapiro_grazing(gm, r1, r2, bb, 1);
    const double lib_sg1 = stc_clocks_shapiro_sun_grazing_us(0);
    const double lib_sg2 = stc_clocks_shapiro_sun_grazing_us(1);
    const double oracle_sep =
        sqrt(r1 * r1 - bb * bb) + sqrt(r2 * r2 - bb * bb);
    const double oracle_tri =
        stc_clocks_shapiro_one_way(gm, r1, r2, oracle_sep);
    const double rm = STC_VERIFY_MARS_AU * STC_AU;
    const double oracle_mars_rt =
        (4.0 * gm / c3) * log(4.0 * r1 * rm / (bb * bb));
    int ok = 1;

    ok = ok && isfinite(lib_one) && isfinite(lib_rt) && isfinite(oracle_tri);
    ok = ok && (stc_verify_reldiff(lib_one, oracle_one_way) < 1.0e-9);
    ok = ok && (stc_verify_reldiff(lib_rt, 2.0 * oracle_one_way) < 1.0e-12);
    ok = ok && (stc_verify_reldiff(lib_sg1, oracle_one_way * 1.0e6) < 1.0e-9);
    ok = ok &&
         (stc_verify_reldiff(lib_sg2, 2.0 * oracle_one_way * 1.0e6) < 1.0e-9);
    ok = ok && (stc_verify_reldiff(oracle_tri, oracle_one_way) < 1.0e-3);
    ok = ok && (stc_verify_reldiff(oracle_one_way * 1.0e6, 120.0) < 0.10);
    ok = ok && (stc_verify_reldiff(oracle_mars_rt * 1.0e6, 250.0) < 0.50);
    printf("%s shapiro_delay: computed=EV 1-way %.3fus 2-way %.3fus Mars "
           "2-way %.3fus expected=120us ~233us 250us tol=10%% order 50%%\n",
           ok ? "PASS" : "FAIL", oracle_one_way * 1.0e6,
           2.0 * oracle_one_way * 1.0e6, oracle_mars_rt * 1.0e6);
    return ok;
}

/* Check 6. Equations: h = 4G^2 m1 m2/(c^4 a R); f(tau) chirp law and its
 * inverse tau(f); seeded phase determinism. */
static int check_gw_chirp(void)
{
    const double m1 = 30.0 * STC_M_SUN;
    const double m2 = 30.0 * STC_M_SUN;
    const double mc = stc_waves_chirp_mass(m1, m2);
    const double dist = 1.0e25;
    const double freq = 50.0;
    const double h1 = stc_waves_strain_amplitude(mc, freq, dist);
    const double h2 = stc_waves_strain_amplitude(mc, freq, 2.0 * dist);
    const double sep = 1.0e9;
    const double f_kep = sqrt(STC_G * (m1 + m2) / (sep * sep * sep)) / STC_PI;
    const double h_sep = stc_waves_quadrupole_strain(m1, m2, sep, dist);
    const double h_amp = stc_waves_strain_amplitude(mc, f_kep, dist);
    const double tau0 = 100.0;
    const double fwd = stc_waves_chirp_freq(mc, tau0);
    const double rt = stc_waves_time_to_coalescence(mc, fwd);
    /* Seeded determinism probe: fixed seed -> identical waveform. */
    StcRng64 rng;
    double hp_a = 0.0;
    double hx_a = 0.0;
    double hp_b = 0.0;
    double hx_b = 0.0;
    double incl;
    double phic;
    int ok = 1;

    stc_rng64_seed(&rng, 0x123456789ABCDEFULL);
    incl = stc_rng64_uniform_range(&rng, 0.0, STC_PI);
    phic = stc_waves_phase0_from_seed(0x123456789ABCDEFULL);
    stc_waves_chirp_waveform(mc, dist, tau0, phic, incl, &hp_a, &hx_a);
    stc_waves_chirp_waveform(mc, dist, tau0, phic, incl, &hp_b, &hx_b);
    ok = ok && isfinite(h1) && isfinite(h2) && isfinite(h_sep);
    ok = ok && isfinite(h_amp) && isfinite(fwd) && isfinite(rt);
    ok = ok && (stc_verify_reldiff(h2 / h1, 0.5) < 1.0e-12);
    ok = ok && (stc_verify_reldiff(h_sep, h_amp) < 1.0e-6);
    ok = ok && (stc_verify_reldiff(rt, tau0) < 1.0e-9);
    ok = ok && (hp_a == hp_b) && (hx_a == hx_b);
    ok = ok && (stc_waves_phase0_from_seed(7ULL) ==
        stc_waves_phase0_from_seed(7ULL));
    printf("%s gw_chirp: computed=scale %.12f ident-relerr %.3e roundtrip "
           "%.9f expected=0.5 0.0 100.0 tol=1e-12 1e-6 1e-9\n",
           ok ? "PASS" : "FAIL", h2 / h1, stc_verify_reldiff(h_sep, h_amp),
           rt);
    return ok;
}

/* Check 7. Equation: Bardeen-Press-Teukolsky ISCO, 6M/1M/9M limits. */
static int check_kerr_isco(void)
{
    const double mm = 1000.0;
    const double isco0 = stc_kerr_isco(mm, 0.0);
    const double isco_p = stc_kerr_isco(mm, mm);
    const double isco_m = stc_kerr_isco(mm, -mm);
    const double pro_p = stc_kerr_isco_prograde(mm, mm);
    const double ret_m = stc_kerr_isco_retrograde(mm, mm);
    const double pro_0 = stc_kerr_isco_prograde(mm, 0.0);
    const double ret_0 = stc_kerr_isco_retrograde(mm, 0.0);
    int ok = 1;

    ok = ok && isfinite(isco0) && isfinite(isco_p) && isfinite(isco_m);
    ok = ok && (stc_verify_reldiff(isco0 / mm, 6.0) < 1.0e-6);
    ok = ok && (stc_verify_reldiff(isco_p / mm, 1.0) < 1.0e-6);
    ok = ok && (stc_verify_reldiff(isco_m / mm, 9.0) < 1.0e-6);
    ok = ok && (stc_verify_reldiff(pro_p / mm, 1.0) < 1.0e-6);
    ok = ok && (stc_verify_reldiff(ret_m / mm, 9.0) < 1.0e-6);
    ok = ok && (stc_verify_reldiff(pro_0 / mm, 6.0) < 1.0e-6);
    ok = ok && (stc_verify_reldiff(ret_0 / mm, 6.0) < 1.0e-6);
    printf("%s kerr_isco: computed=%.9fM %.9fM %.9fM expected=6M 1M 9M "
           "tol=1e-6\n",
           ok ? "PASS" : "FAIL", isco0 / mm, isco_p / mm, isco_m / mm);
    return ok;
}

/* Check 8. Equations: circular velocity v = sqrt(GM/r); Kepler-III period
 * T = 2*pi*sqrt(r^3/GM); velocity Verlet orbit + total energy conservation.
 * Earth-mass planet, circular r=1e7 about 5.972e24 fixed center, dt=T/2000,
 * one period: |r-r0|/r0 < 1e-3 AND |E-E0|/|E0| < 1e-3. */
static int check_nbody_orbit(void)
{
    const double m_center = 5.972e24;
    const double m_planet = 5.972e24;
    const double r0 = 1.0e7;
    const double mu = STC_G * m_center;
    const double vcirc = sqrt(mu / r0);
    const double period = 2.0 * STC_PI * sqrt(r0 * r0 * r0 / mu);
    const double dt = period / 2000.0;
    const int steps = 2000;
    StcSim s;
    StcVec3 pc = { 0.0, 0.0, 0.0 };
    StcVec3 vc = { 0.0, 0.0, 0.0 };
    StcVec3 pp = { 1.0e7, 0.0, 0.0 };
    StcVec3 vp = { 0.0, 0.0, 0.0 };
    double e0 = 0.0;
    double e1 = 0.0;
    double r1 = 0.0;
    double rerr = 1.0;
    double eerr = 1.0;
    int ok = 1;
    int i;

    vp.y = vcirc;
    stc_sim_init(&s, STC_G, 0.0, 0.0, 0);
    ok = ok && (stc_sim_add(&s, "center", m_center, 0.0, pc, vc, 1) == 0);
    ok = ok && (stc_sim_add(&s, "planet", m_planet, 0.0, pp, vp, 0) == 1);
    if (ok) {
        e0 = stc_sim_energy(&s);
        for (i = 0; i < steps; i++) {
            stc_sim_step(&s, dt);
        }
        ok = ok && (s.n == 2);
        ok = ok && isfinite(e0);
        if (s.n == 2) {
            r1 = stc_vec3_dist(s.b[1].p, s.b[0].p);
            e1 = stc_sim_energy(&s);
            rerr = fabs(r1 - r0) / r0;
            if (isfinite(e1) && e0 != 0.0) {
                eerr = fabs(e1 - e0) / fabs(e0);
            } else {
                eerr = 1.0;
            }
            ok = ok && isfinite(r1) && isfinite(e1);
            ok = ok && (rerr < 1.0e-3);
            ok = ok && (eerr < 1.0e-3);
        }
    }
    printf("%s nbody_orbit: computed=r %.6e rerr %.3e eerr %.3e "
           "expected=r0 %.6e 0.0 0.0 tol=1e-3 1e-3\n",
           ok ? "PASS" : "FAIL", r1, rerr, eerr, r0);
    return ok;
}

/* Check 9. Equation: momentum-conserving merge v = (m1*v1+m2*v2)/(m1+m2).
 * Two equal free masses head-on with merge on: |v_merge| < 1e-9 and mass
 * doubled (n == 1). */
static int check_momentum_merge(void)
{
    const double m0 = 1.0e12;
    const double rr = 10.0;
    StcSim s;
    StcVec3 p0 = { -5.0, 0.0, 0.0 };
    StcVec3 v0 = { 100.0, 0.0, 0.0 };
    StcVec3 p1 = { 5.0, 0.0, 0.0 };
    StcVec3 v1 = { -100.0, 0.0, 0.0 };
    double vlen = HUGE_VAL;
    double merr = 1.0;
    int ok = 1;

    stc_sim_init(&s, STC_G, 0.0, 0.0, 1);
    ok = ok && (stc_sim_add(&s, "a", m0, rr, p0, v0, 0) == 0);
    ok = ok && (stc_sim_add(&s, "b", m0, rr, p1, v1, 0) == 1);
    if (ok) {
        stc_sim_step(&s, 0.01);
        ok = ok && (s.n == 1);
        if (s.n == 1) {
            vlen = stc_vec3_len(s.b[0].v);
            merr = fabs(s.b[0].m - 2.0 * m0) / (2.0 * m0);
            ok = ok && isfinite(vlen) && isfinite(merr);
            ok = ok && (vlen < 1.0e-9);
            ok = ok && (merr < 1.0e-12);
        }
    }
    printf("%s momentum_merge: computed=|v| %.3e merr %.3e n %d "
           "expected=0.0 0.0 1 tol=1e-9 1e-12\n",
           ok ? "PASS" : "FAIL", vlen, merr, s.n);
    return ok;
}

/* Check 10. Same physics and shape as check 9, but the integrator is
 * stc_sim_step_rk4 with the collision/cap pass applied explicitly via
 * stc_sim_poststep after every step (the CONTRACT Gate-0 seam: rk4 has
 * no collisions inside; poststep is caller-owned). Two equal free
 * masses head-on with merge on: |v_merge| < 1e-9, mass doubled, n == 1. */
static int check_rk4_poststep(void)
{
    const double m0 = 1.0e12;
    const double rr = 10.0;
    const int steps = 5;
    StcSim s;
    StcVec3 p0 = { -5.0, 0.0, 0.0 };
    StcVec3 v0 = { 100.0, 0.0, 0.0 };
    StcVec3 p1 = { 5.0, 0.0, 0.0 };
    StcVec3 v1 = { -100.0, 0.0, 0.0 };
    double vlen = HUGE_VAL;
    double merr = 1.0;
    int ok = 1;
    int i;

    stc_sim_init(&s, STC_G, 0.0, 0.0, 1);
    ok = ok && (stc_sim_add(&s, "a", m0, rr, p0, v0, 0) == 0);
    ok = ok && (stc_sim_add(&s, "b", m0, rr, p1, v1, 0) == 1);
    if (ok) {
        for (i = 0; i < steps; i++) {
            stc_sim_step_rk4(&s, 0.01);
            stc_sim_poststep(&s);
        }
        ok = ok && (s.n == 1);
        if (s.n == 1) {
            vlen = stc_vec3_len(s.b[0].v);
            merr = fabs(s.b[0].m - 2.0 * m0) / (2.0 * m0);
            ok = ok && isfinite(vlen) && isfinite(merr);
            ok = ok && (vlen < 1.0e-9);
            ok = ok && (merr < 1.0e-12);
        }
    }
    printf("%s rk4_poststep: computed=|v| %.3e merr %.3e n %d "
           "expected=0.0 0.0 1 tol=1e-9 1e-12\n",
           ok ? "PASS" : "FAIL", vlen, merr, s.n);
    return ok;
}

/* Check 11. Equations: TaylorT4 3.5PN inspiral (arXiv:0907.0700 Eq. 6,
 * tabulated by RK4 in stc_bbh_init): internal consistency only --
 * dphi/dt == omega = v^3 by central difference at 5 interior inspiral
 * times, a(t_m) = A_ISCO = 6 and omega(t_m) = 6^(-3/2) (ISCO boundary),
 * phase AND frequency continuity across the inspiral->plunge (t_m) and
 * plunge->ringdown (t_peak) seams, remnant mass/spin phenomenological
 * fits, determinism of stc_bbh_eval. Every measured quantity comes from
 * the shipped stc_bbh_* API (src/merger/merger.c), q = 1, a0 = 11 M.
 * Tolerances (chosen before the run): omega-relerr < 1e-6, seam jumps
 * < 1e-6, ISCO boundary values < 1e-9 relative. */
static int check_bbh_timeline(void)
{
    StcBbh b;
    StcBbhState s;
    StcBbhState sp;
    StcBbhState sm;
    StcBbhState sa;
    StcBbhState sc;
    StcBbhState s1;
    StcBbhState s2;
    const double h = 1.0e-3;
    const double eps = 1.0e-9;
    static const double frac[5] = {0.06, 0.25, 0.50, 0.75, 0.95};
    double t;
    double a_isco;
    double dphi_dt;
    double omega_err = 0.0;
    double om_err_m;
    double jump_m;
    double jump_peak;
    double om_jump_m;
    int ok = 1;
    int j;

    if (stc_bbh_init(&b, 1.0, 11.0, 3.4) != 0) {
        printf("FAIL bbh_timeline: stc_bbh_init rejected q=1 a0=11\n");
        return 0;
    }

    /* dphi/dt vs the state's omega = v^3 at 5 interior inspiral times. */
    for (j = 0; j < 5; j++) {
        double e;

        t = frac[j] * b.t_m;
        stc_bbh_eval(&b, t, &s);
        stc_bbh_eval(&b, t + h, &sp);
        stc_bbh_eval(&b, t - h, &sm);
        dphi_dt = (sp.phi - sm.phi) / (2.0 * h);
        e = stc_verify_reldiff(dphi_dt, s.omega);
        if (e > omega_err) {
            omega_err = e;
        }
    }

    /* Phase and frequency continuity across the stage boundaries
     * (inspiral->plunge at t_m, plunge->ringdown at t_peak): the
     * piecewise model is built to match at the seams, so a t +/- eps
     * probe should show no jump. */
    stc_bbh_eval(&b, b.t_m - eps, &sa);
    a_isco = sa.a;
    stc_bbh_eval(&b, b.t_m + eps, &sc);
    jump_m = fabs(sc.phi - sa.phi);
    om_err_m = stc_verify_reldiff(sc.omega, sa.omega);
    stc_bbh_eval(&b, b.t_peak - eps, &sa);
    stc_bbh_eval(&b, b.t_peak + eps, &sc);
    jump_peak = fabs(sc.phi - sa.phi);
    om_jump_m = om_err_m;

    /* Determinism: identical (b, t) -> byte-identical state. */
    stc_bbh_eval(&b, 0.5 * b.t_m, &s1);
    stc_bbh_eval(&b, 0.5 * b.t_m, &s2);

    ok = ok && isfinite(a_isco) && isfinite(omega_err) &&
         isfinite(b.mf) && isfinite(b.af);
    ok = ok && (omega_err < 1.0e-6);
    ok = ok && (stc_verify_reldiff(a_isco, STC_BBH_A_ISCO) < 1.0e-9);
    ok = ok && (stc_verify_reldiff(b.om_m, pow(STC_BBH_V_ISCO, 3.0)) <
                1.0e-12);
    ok = ok && (b.mf > 0.94) && (b.mf < 0.96);
    ok = ok && (b.af > 0.66) && (b.af < 0.71);
    ok = ok && (jump_m < 1.0e-6) && (jump_peak < 1.0e-6);
    ok = ok && (om_jump_m < 1.0e-6);
    ok = ok && (memcmp(&s1, &s2, sizeof s1) == 0);
    printf("%s bbh_timeline: computed=omega-relerr %.3e a(t_m) %.6f "
           "om_m %.6f mf %.4f af %.4f jump_m %.3e jump_peak %.3e "
           "omjump_m %.3e "
           "expected=0 6 0.0680 (0.94,0.96) (0.66,0.71) 0 0 0 "
           "tol=1e-6 1e-9 1e-12 range range 1e-6 1e-6 1e-6\n",
           ok ? "PASS" : "FAIL", omega_err, a_isco, b.om_m, b.mf, b.af,
           jump_m, jump_peak, om_jump_m);
    stc_bbh_free(&b);
    return ok;
}

/* TaylorT2 3.5PN closed-form t(v) bracket B(v), t(v) = t_ref -
 * (5M/256 nu) v^(-8) B(v): Buonanno et al. 2009, arXiv:0907.0700,
 * Eq. (8b). Independent oracle for the shipped TaylorT4 table -- the
 * two approximants differ only at >3.5PN order. */
static double verify_t2_bracket(double eta, double v)
{
    const double gamma = 0.57721566490153286061; /* Euler-Mascheroni */
    double v2 = v * v;
    double d2;
    double d3;
    double d4;
    double d5;
    double d6;
    double d7;

    d2 = 743.0 / 252.0 + 11.0 / 3.0 * eta;
    d3 = -32.0 / 5.0 * M_PI;
    d4 = 3058673.0 / 508032.0 + 5429.0 / 504.0 * eta +
         617.0 / 72.0 * eta * eta;
    d5 = -(7729.0 / 252.0 - 13.0 / 3.0 * eta) * M_PI;
    d6 = -10052469856691.0 / 23471078400.0 + 128.0 / 3.0 * M_PI * M_PI +
         6848.0 / 105.0 * gamma +
         (3147553127.0 / 3048192.0 - 451.0 / 12.0 * M_PI * M_PI) * eta -
         15211.0 / 1728.0 * eta * eta + 25565.0 / 1296.0 * eta * eta *
         eta + 3424.0 / 105.0 * log(16.0 * v2);
    d7 = (-15419335.0 / 127008.0 - 75703.0 / 756.0 * eta +
          14809.0 / 378.0 * eta * eta) * M_PI;
    return 1.0 + d2 * v2 + d3 * v2 * v + d4 * v2 * v2 +
           d5 * v2 * v2 * v + d6 * v2 * v2 * v2 +
           d7 * v2 * v2 * v2 * v;
}

/* Check 12. Cross-verification of the production TaylorT4 inspiral
 * against the independent TaylorT2 closed form: for a GW150914-style
 * binary (q = 1.24, M = 65 Msun), time from f_GW = 35 Hz
 * (v0 = (pi M f)^(1/3), a0 = v0^(-2), leading-order Kepler) to the ISCO
 * v = 6^(-1/2) from the shipped table (b.t_m) vs the T2 oracle Eq. (8b).
 * Tolerance chosen before the run: 10% -- the approximants differ at
 * >3.5PN order, and the report records the actual difference. */
static int check_bbh_t2_crosscheck(void)
{
    StcBbh b;
    const double q = 1.24;
    const double f_gw = 35.0;
    const double v_isco = STC_BBH_V_ISCO;
    double eta = q / ((1.0 + q) * (1.0 + q));
    double m_sec = 65.0 * STC_GM_SUN / (STC_C * STC_C * STC_C);
    double v0;
    double a0;
    double t2_m;
    double t4_s;
    double t2_s;
    double rel;
    int ok;

    v0 = pow(STC_PI * m_sec * f_gw, 1.0 / 3.0);
    a0 = 1.0 / (v0 * v0);
    if (stc_bbh_init(&b, q, a0, 3.4) != 0) {
        printf("FAIL bbh_t2_crosscheck: stc_bbh_init rejected "
               "q=1.24 a0=%.4f\n", a0);
        return 0;
    }
    t4_s = b.t_m * m_sec;
    t2_m = 5.0 / (256.0 * eta) *
           (verify_t2_bracket(eta, v0) / pow(v0, 8.0) -
            verify_t2_bracket(eta, v_isco) / pow(v_isco, 8.0));
    t2_s = t2_m * m_sec;
    rel = fabs(t4_s - t2_s) / t2_s;
    ok = isfinite(rel) && (rel < 0.10);
    printf("%s bbh_t2_crosscheck: computed=T4 %.4f s vs T2 %.4f s "
           "(%.2f%%) expected=agree tol=10%%\n",
           ok ? "PASS" : "FAIL", t4_s, t2_s, 100.0 * rel);
    stc_bbh_free(&b);
    return ok;
}

/* INFO (not pass/fail). GW150914 context: for the same q = 1.24,
 * 65 Msun binary, the production model's 35 Hz -> t_peak duration and
 * accumulated GW cycle count (Delta phi_GW / 2 pi = phi_peak / pi since
 * phi(0) = 0 at 35 Hz). LIGO (Abbott et al. 2016, arXiv:1602.03837)
 * quotes ">0.2 s ... from 35 to 150 Hz" and ~8 cycles -- that interval
 * includes the post-ISCO plunge/merger a PN inspiral does not describe,
 * so this is context for the eyeball, not a gate. */
static void check_bbh_gw150914_context(void)
{
    StcBbh b;
    StcBbhState s;
    const double q = 1.24;
    const double f_gw = 35.0;
    double m_sec = 65.0 * STC_GM_SUN / (STC_C * STC_C * STC_C);
    double v0;
    double a0;
    double cycles;

    v0 = pow(STC_PI * m_sec * f_gw, 1.0 / 3.0);
    a0 = 1.0 / (v0 * v0);
    if (stc_bbh_init(&b, q, a0, 3.4) != 0) {
        printf("INFO bbh_gw150914_context: init failed (q=1.24, "
               "a0=%.4f)\n", a0);
        return;
    }
    stc_bbh_eval(&b, b.t_peak, &s);
    cycles = s.phi / STC_PI;
    printf("INFO bbh_gw150914_context: 35 Hz -> ISCO %.4f s, "
           "35 Hz -> t_peak %.4f s, %.1f GW cycles "
           "(LIGO context: ~0.2 s and ~8 cycles span 35-150 Hz "
           "including plunge/merger)\n",
           b.t_m * m_sec, b.t_peak * m_sec, cycles);
    stc_bbh_free(&b);
}

int main(void)
{
    int passed = 0;

    passed += check_sr_invariance();
    passed += check_photon_bending();
    passed += check_mercury_precession();
    passed += check_gps_offset();
    passed += check_shapiro_delay();
    passed += check_gw_chirp();
    passed += check_kerr_isco();
    passed += check_nbody_orbit();
    passed += check_momentum_merge();
    passed += check_rk4_poststep();
    passed += check_bbh_timeline();
    passed += check_bbh_t2_crosscheck();
    check_bbh_gw150914_context();
    printf("stverify: %d/12 PASS\n", passed);
    return (passed == 12) ? 0 : 1;
}
