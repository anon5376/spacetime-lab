/* SPACETIME LAB — merger: tabulated TaylorT4 3.5PN BBH
 * inspiral-merger-ringdown timeline (see stc/merger.h for units,
 * the model equations, and claim classes).
 *
 * The inspiral table is built once in stc_bbh_init by fixed-step
 * classic RK4 on (v, phi) with v' = F_T4(v) (arXiv:0907.0700 Eq. 6)
 * and phi' = v^3, then evaluated by cubic Hermite. All step counts and
 * bisection iterations are fixed: same (q, a0) -> byte-identical table.
 */

#include "stc/merger.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

/* Table budget: backward fill to t_lo = -400 M for retarded lookups,
 * then forward to the ISCO. <= 200k nodes total (4 doubles each). */
#define BBH_T_LO (-400.0)
#define BBH_N_BACK 4096                  /* fixed backward steps */
#define BBH_H_BACK (400.0 / BBH_N_BACK)  /* 0.09765625 M */
#define BBH_N_FWD_CAP (200000 - BBH_N_BACK - 8)
#define BBH_H_FWD_MIN 0.05               /* M */
#define BBH_H_FWD_MAX 0.5                /* M: << v/(dv/dt) at the ISCO */
#define BBH_H_FWD_SPAN 60000.0           /* target nodes = t_est/this */
#define BBH_ISCO_BISECT 60               /* fixed root-solve iterations */

/* TaylorT4 3.5PN dv/dt = (32/5) nu v^9 B(v), M = 1
 * (arXiv:0907.0700 Eq. 6). */
static double bbh_dvdt(double eta, double v)
{
    const double gamma = 0.57721566490153286061; /* Euler-Mascheroni */
    double v2 = v * v;
    double b2;
    double b4;
    double b5;
    double b6;
    double b7;
    double br;

    b2 = -(743.0 / 336.0 + 11.0 / 4.0 * eta);
    b4 = 34103.0 / 18144.0 + 13661.0 / 2016.0 * eta +
         59.0 / 18.0 * eta * eta;
    b5 = -(4159.0 / 672.0 + 189.0 / 8.0 * eta) * M_PI;
    b6 = 16447322263.0 / 139708800.0 + 16.0 / 3.0 * M_PI * M_PI -
         1712.0 / 105.0 * gamma +
         (451.0 / 48.0 * M_PI * M_PI - 56198689.0 / 217728.0) * eta +
         541.0 / 896.0 * eta * eta - 5605.0 / 2592.0 * eta * eta * eta -
         856.0 / 105.0 * log(16.0 * v2);
    b7 = -(4415.0 / 4032.0 - 358675.0 / 6048.0 * eta -
           91495.0 / 1512.0 * eta * eta) * M_PI;
    br = 1.0 + b2 * v2 + 4.0 * M_PI * v2 * v + b4 * v2 * v2 +
         b5 * v2 * v2 * v + b6 * v2 * v2 * v2 +
         b7 * v2 * v2 * v2 * v;
    return 32.0 / 5.0 * eta * v2 * v2 * v2 * v2 * v * br;
}

/* One fixed RK4 step on y = (v, phi): v' = F(v), phi' = v^3. */
static void bbh_rk4(double eta, double h, double *v, double *phi)
{
    double k1v = bbh_dvdt(eta, *v);
    double k1p = *v * *v * *v;
    double v2 = *v + 0.5 * h * k1v;
    double k2v = bbh_dvdt(eta, v2);
    double k2p = v2 * v2 * v2;
    double v3 = *v + 0.5 * h * k2v;
    double k3v = bbh_dvdt(eta, v3);
    double k3p = v3 * v3 * v3;
    double v4 = *v + h * k3v;
    double k4v = bbh_dvdt(eta, v4);
    double k4p = v4 * v4 * v4;

    *v += h / 6.0 * (k1v + 2.0 * k2v + 2.0 * k3v + k4v);
    *phi += h / 6.0 * (k1p + 2.0 * k2p + 2.0 * k3p + k4p);
}

/* Cubic Hermite on interval i (node fields: t, v, dv/dt, phi; the phi
 * derivative is v^3). s in [0,1] normally, extrapolates outside. */
static void bbh_hermite(const StcBbh *b, int i, double t, double *v,
                        double *phi)
{
    const double *n0 = b->tab + (size_t)i * 4u;
    const double *n1 = n0 + 4;
    double h = n1[0] - n0[0];
    double s = (t - n0[0]) / h;
    double s2 = s * s;
    double s3 = s2 * s;
    double h00 = 2.0 * s3 - 3.0 * s2 + 1.0;
    double h10 = s3 - 2.0 * s2 + s;
    double h01 = -2.0 * s3 + 3.0 * s2;
    double h11 = s3 - s2;
    double m0;
    double m1;

    *v = h00 * n0[1] + h10 * h * n0[2] + h01 * n1[1] + h11 * h * n1[2];
    m0 = n0[1] * n0[1] * n0[1];
    m1 = n1[1] * n1[1] * n1[1];
    *phi = h00 * n0[3] + h10 * h * m0 + h01 * n1[3] + h11 * h * m1;
}

/* Inspiral (v, phi) at t: binary search + Hermite over the tabulated
 * TaylorT4 trajectory. Requires b->tab_n >= 2. */
static void bbh_inspiral(const StcBbh *b, double t, double *v, double *phi)
{
    int lo = 0;
    int hi = b->tab_n - 1;

    while (hi - lo > 1) {
        int mid = lo + (hi - lo) / 2;

        if (b->tab[(size_t)mid * 4u] <= t) {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    if (t >= b->tab[(size_t)(b->tab_n - 1) * 4u]) {
        lo = b->tab_n - 2;
    }
    bbh_hermite(b, lo, t, v, phi);
}

static double bbh_smooth(double x)
{
    if (x <= 0.0) {
        return 0.0;
    }
    if (x >= 1.0) {
        return 1.0;
    }
    return x * x * (3.0 - 2.0 * x);
}

int stc_bbh_init(StcBbh *b, double q, double a0, double a_merge)
{
    double eta;
    double om_r;
    double tp = STC_BBH_PLUNGE_T;
    double erad;
    double one_minus_a;
    double v0;
    double v;
    double phi;
    double t_est;
    double h_fwd;
    double t;
    double lo;
    double hi;
    int n_fwd_alloc;
    int n;
    int i;
    int steps;

    if (b == NULL) {
        return -1;
    }
    memset(b, 0, sizeof(*b));
    if (!(q >= 1.0) || !(q <= 1.0e6) || !(a_merge > 0.0) ||
        !(a0 > a_merge) || !(a0 > STC_BBH_A_ISCO) || !isfinite(a0)) {
        return -1;
    }
    b->m1 = q / (1.0 + q);
    b->m2 = 1.0 / (1.0 + q);
    eta = b->m1 * b->m2;
    b->eta = eta;
    b->a0 = a0;
    b->a_merge = STC_BBH_A_ISCO;
    b->t_c = 5.0 * pow(a0, 4.0) / (256.0 * eta); /* Peters reference */
    b->t_lo = BBH_T_LO;

    /* Fixed forward step: sized so the 0PN estimate of the inspiral
     * span lands near BBH_H_FWD_SPAN nodes, clamped to
     * [BBH_H_FWD_MIN, BBH_H_FWD_MAX] (the ISCO dynamical time
     * v/(dv/dt) ~ 140 M is always >> h). The 3.5PN inspiral is faster
     * than 0PN, so t_est overestimates the span. */
    t_est = 5.0 * (pow(a0, 4.0) - 1296.0) / (256.0 * eta);
    h_fwd = t_est / BBH_H_FWD_SPAN;
    if (h_fwd < BBH_H_FWD_MIN) {
        h_fwd = BBH_H_FWD_MIN;
    }
    if (h_fwd > BBH_H_FWD_MAX) {
        h_fwd = BBH_H_FWD_MAX;
    }
    n_fwd_alloc = (int)(t_est / h_fwd) + 4096;
    if (n_fwd_alloc > BBH_N_FWD_CAP) {
        n_fwd_alloc = BBH_N_FWD_CAP;
    }
    if (t_est / BBH_H_FWD_MIN > (double)BBH_N_FWD_CAP) {
        return -1; /* span exceeds the 200k-node table budget */
    }

    b->tab = malloc((size_t)(BBH_N_BACK + n_fwd_alloc) * 4u *
                    sizeof(double));
    if (b->tab == NULL) {
        return -1;
    }

    /* Node BBH_N_BACK is t = 0, v = v0 = a0^(-1/2) (leading-order
     * Kepler omega = a^(-3/2)), phi = 0. */
    v0 = pow(a0, -0.5);
    b->tab[(size_t)BBH_N_BACK * 4u + 0u] = 0.0;
    b->tab[(size_t)BBH_N_BACK * 4u + 1u] = v0;
    b->tab[(size_t)BBH_N_BACK * 4u + 2u] = bbh_dvdt(eta, v0);
    b->tab[(size_t)BBH_N_BACK * 4u + 3u] = 0.0;
    if (!(b->tab[(size_t)BBH_N_BACK * 4u + 2u] > 0.0)) {
        free(b->tab);
        memset(b, 0, sizeof(*b));
        return -1;
    }

    /* Backward fill to t_lo = -400 M (retarded-time lookups): fixed
     * step -BBH_H_BACK, nodes stored in reverse so t stays increasing
     * with index. */
    v = v0;
    phi = 0.0;
    for (i = 1; i <= BBH_N_BACK; i++) {
        bbh_rk4(eta, -BBH_H_BACK, &v, &phi);
        if (!(v > 0.0) || !isfinite(v) || !isfinite(phi)) {
            free(b->tab);
            memset(b, 0, sizeof(*b));
            return -1;
        }
        n = BBH_N_BACK - i;
        b->tab[(size_t)n * 4u + 0u] = -(double)i * BBH_H_BACK;
        b->tab[(size_t)n * 4u + 1u] = v;
        b->tab[(size_t)n * 4u + 2u] = bbh_dvdt(eta, v);
        b->tab[(size_t)n * 4u + 3u] = phi;
    }

    /* Forward fill to just past the ISCO v = 6^(-1/2). */
    v = v0;
    phi = 0.0;
    t = 0.0;
    steps = 0;
    n = BBH_N_BACK;
    while (v < STC_BBH_V_ISCO) {
        if (n >= BBH_N_BACK + n_fwd_alloc - 1) {
            free(b->tab);
            memset(b, 0, sizeof(*b));
            return -1; /* 0PN estimate undershot; span exceeds budget */
        }
        bbh_rk4(eta, h_fwd, &v, &phi);
        t += h_fwd;
        steps++;
        if (!isfinite(v) || !isfinite(phi) || !(v > 0.0)) {
            free(b->tab);
            memset(b, 0, sizeof(*b));
            return -1;
        }
        n++;
        b->tab[(size_t)n * 4u + 0u] = t;
        b->tab[(size_t)n * 4u + 1u] = v;
        b->tab[(size_t)n * 4u + 2u] = bbh_dvdt(eta, v);
        b->tab[(size_t)n * 4u + 3u] = phi;
    }
    b->tab_n = n + 1;
    (void)steps;

    /* t_m: bisection for v(t) = V_ISCO on the last (straddling)
     * interval -- fixed iteration count, deterministic. */
    lo = b->tab[(size_t)(b->tab_n - 2) * 4u];
    hi = b->tab[(size_t)(b->tab_n - 1) * 4u];
    {
        double vm;
        double pm;
        double mid;
        int it;
        int idx = b->tab_n - 2;

        for (it = 0; it < BBH_ISCO_BISECT; it++) {
            mid = 0.5 * (lo + hi);
            bbh_hermite(b, idx, mid, &vm, &pm);
            if (vm < STC_BBH_V_ISCO) {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        b->t_m = 0.5 * (lo + hi);
        bbh_hermite(b, idx, b->t_m, &vm, &pm);
        b->phi_m = pm;
    }
    b->om_m = STC_BBH_V_ISCO * STC_BBH_V_ISCO * STC_BBH_V_ISCO;
    b->amp_m = 4.0 * eta * STC_BBH_V_ISCO * STC_BBH_V_ISCO;
    b->amp_peak = 1.4 * b->amp_m;
    b->t_peak = b->t_m + tp;

    /* Remnant (ILLUSTRATIVE phenomenological fits, non-spinning
     * progenitors): radiated energy and final spin as polynomials in eta
     * tuned to ~4.5% and ~0.69 at equal mass. */
    erad = eta * (1.0 - sqrt(8.0 / 9.0)) + 0.498 * eta * eta;
    b->mf = 1.0 - erad;
    b->af = 2.0 * sqrt(3.0) * eta - 3.5171 * eta * eta +
            2.5763 * eta * eta * eta;
    /* l=m=2 n=0 Kerr QNM (phenomenological fit in spin). */
    one_minus_a = 1.0 - b->af;
    b->om_qnm = (1.5251 - 1.1568 * pow(one_minus_a, 0.1292)) / b->mf;
    b->tau_qnm = 2.0 * (0.7 + 1.4187 * pow(one_minus_a, -0.4990)) /
                 b->om_qnm;

    om_r = 0.5 * b->om_qnm;
    b->phi_peak = b->phi_m + tp * (b->om_m + 0.5 * (om_r - b->om_m));
    return 0;
}

void stc_bbh_free(StcBbh *b)
{
    if (b == NULL) {
        return;
    }
    free(b->tab);
    b->tab = NULL;
    b->tab_n = 0;
}

/* Phase, orbital frequency, amplitude only (no positions): hot path. */
static void bbh_wave(const StcBbh *b, double t, double *phi, double *amp)
{
    double om_r = 0.5 * b->om_qnm;

    if (t < b->t_m) {
        double v;

        bbh_inspiral(b, t, &v, phi);
        *amp = 4.0 * b->eta * v * v;
    } else if (t < b->t_peak) {
        double s = (t - b->t_m) / STC_BBH_PLUNGE_T;

        *phi = b->phi_m + STC_BBH_PLUNGE_T *
               (b->om_m * s + 0.5 * (om_r - b->om_m) * s * s);
        *amp = b->amp_m + (b->amp_peak - b->amp_m) * bbh_smooth(s);
    } else {
        double u = t - b->t_peak;

        *phi = b->phi_peak + om_r * u;
        *amp = b->amp_peak * exp(-u / b->tau_qnm);
    }
}

void stc_bbh_eval(const StcBbh *b, double t, StcBbhState *s)
{
    double c;
    double sn;

    if (s == NULL) {
        return;
    }
    memset(s, 0, sizeof(*s));
    if (b == NULL || !(b->eta > 0.0) || b->tab == NULL ||
        b->tab_n < 2) {
        return;
    }
    bbh_wave(b, t, &s->phi, &s->amp);
    if (t < b->t_m) {
        double v;

        s->stage = STC_BBH_INSPIRAL;
        bbh_inspiral(b, t, &v, &sn);
        s->a = 1.0 / (v * v);
        s->omega = v * v * v;
        s->mass = 1.0;
        s->fuse = 0.0;
    } else if (t < b->t_peak) {
        double u = (t - b->t_m) / STC_BBH_PLUNGE_T;

        s->stage = STC_BBH_PLUNGE;
        s->a = b->a_merge * (1.0 - u * u);
        s->omega = b->om_m + (0.5 * b->om_qnm - b->om_m) * u;
        s->fuse = bbh_smooth(u);
        s->mass = 1.0 + (b->mf - 1.0) * s->fuse;
    } else {
        s->stage = STC_BBH_RINGDOWN;
        s->a = 0.0;
        s->omega = 0.5 * b->om_qnm;
        s->fuse = 1.0;
        s->mass = b->mf;
    }
    c = cos(s->phi);
    sn = sin(s->phi);
    s->x1[0] = b->m2 * s->a * c;
    s->x1[1] = b->m2 * s->a * sn;
    s->x2[0] = -b->m1 * s->a * c;
    s->x2[1] = -b->m1 * s->a * sn;
}

double stc_bbh_rh_plus(const StcBbh *b, double t_ret, double az)
{
    double phi;
    double amp;

    if (b == NULL || !(b->eta > 0.0) || b->tab == NULL) {
        return 0.0;
    }
    bbh_wave(b, t_ret, &phi, &amp);
    return 0.5 * amp * cos(2.0 * (phi - az));
}

double stc_bbh_rh_faceon(const StcBbh *b, double t_ret)
{
    double phi;
    double amp;

    if (b == NULL || !(b->eta > 0.0) || b->tab == NULL) {
        return 0.0;
    }
    bbh_wave(b, t_ret, &phi, &amp);
    return amp * cos(2.0 * phi);
}
