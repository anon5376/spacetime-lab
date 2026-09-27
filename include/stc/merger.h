#ifndef STC_MERGER_H
#define STC_MERGER_H

/* SPACETIME LAB — merger: binary black hole inspiral-merger-ringdown
 * timeline for the cinematic merger view (--merger).
 *
 * Geometric units G = c = 1, total initial mass M = m1 + m2 = 1, so
 * lengths and times are in M (1 M = G M / c^3 seconds, G M / c^2 meters).
 * The inspiral is a tabulated solution of the TaylorT4 3.5PN equations,
 * evaluated by cubic Hermite interpolation; that makes retarded-time
 * lookups (t - r) deterministic and cheap.
 *
 * Inspiral model (MEASURED): v(t) = (M omega)^(1/3) and phi(t) are
 * integrated with fixed-step classic RK4 from
 *
 *   dv/dt = (32/5) nu v^9 [ 1 - (743/336 + 11 nu/4) v^2 + 4 pi v^3
 *           + (34103/18144 + 13661 nu/2016 + 59 nu^2/18) v^4
 *           - (4159/672 + 189 nu/8) pi v^5
 *           + (16447322263/139708800 + 16 pi^2/3 - 1712 gamma/105
 *              + (451 pi^2/48 - 56198689/217728) nu + 541 nu^2/896
 *              - 5605 nu^3/2592 - (856/105) ln(16 v^2)) v^6
 *           - (4415/4032 - 358675 nu/6048 - 91495 nu^2/1512) pi v^7 ],
 *   dphi/dt = v^3,
 *
 * Buonanno, Iyer, Ochsner, Pan & Sathyaprakash 2009, arXiv:0907.0700,
 * Eq. (6) (M = 1, nu = eta). The integration starts at t = 0 with
 * v0 = a0^(-1/2) (leading-order Kepler omega = a^(-3/2), phi(0) = 0),
 * runs forward to the ISCO v = 6^(-1/2) and backward to t_lo = -400 M
 * for retarded-time lookups, and is tabulated on <= 200k nodes.
 * Evaluation is cubic Hermite on (v, dv/dt) and (phi, v^3) node data.
 * Separation a = v^(-2) and strain scale r h = 4 eta v^2 are
 * leading-order (Newtonian) relations carried on the PN trajectory.
 *
 * Claim classes:
 *   MEASURED      inspiral v(t), phi(t): RK4 evolution of the stated
 *                 3.5PN TaylorT4 equations, Hermite-evaluated.
 *   ANALYTIC      a = v^(-2), omega = v^3, r h = 4 eta v^2 on that
 *                 trajectory; remnant mass/spin/QNM phenomenological fits.
 *   ILLUSTRATIVE  plunge (smooth blend from the ISCO to coalescence),
 *                 ringdown (single l=m=2 damped sinusoid). Shapes and
 *                 scales are right; this is not a numerical-relativity
 *                 waveform.
 */

/* Stage codes. */
#define STC_BBH_INSPIRAL 0
#define STC_BBH_PLUNGE 1
#define STC_BBH_RINGDOWN 2

/* Plunge duration from the ISCO to the common-horizon peak, M. */
#define STC_BBH_PLUNGE_T 14.0

/* ISCO velocity and separation (Schwarzschild, M = 1): the PN inspiral
 * ends where PN itself breaks down; the plunge blend takes over there. */
#define STC_BBH_V_ISCO 0.40824829046386301637 /* 6^(-1/2) */
#define STC_BBH_A_ISCO 6.0

typedef struct {
    double m1;       /* primary mass, M (m1 >= m2, m1 + m2 = 1) */
    double m2;       /* secondary mass, M */
    double eta;      /* symmetric mass ratio m1 m2 / M^2 */
    double a0;       /* separation at t = 0, M (leading-order Kepler) */
    double a_merge;  /* separation where the plunge starts = A_ISCO = 6 M */
    double t_c;      /* 0PN (Peters) coalescence time from t = 0, M --
                        retained as a reference, not used in the timeline */
    double t_m;      /* plunge start: v reaches V_ISCO, M */
    double t_peak;   /* strain peak / common horizon, M */
    double phi_m;    /* orbital phase at t_m, rad */
    double om_m;     /* orbital angular frequency at t_m = 6^(-3/2), rad/M */
    double phi_peak; /* orbital phase at t_peak, rad */
    double amp_m;    /* r h at t_m = 4 eta v^2 = 2 eta / 3 */
    double amp_peak; /* r h at t_peak */
    double mf;       /* remnant mass, M */
    double af;       /* remnant dimensionless spin */
    double om_qnm;   /* l=m=2 n=0 QNM angular frequency, rad/M */
    double tau_qnm;  /* QNM e-folding time, M */

    /* Inspiral table, heap-allocated by stc_bbh_init, owned by this
     * struct; release with stc_bbh_free. tab[i*4 + {0,1,2,3}] =
     * {t, v, dv/dt, phi} at node i, t strictly increasing over
     * [t_lo, last node at/just past t_m]. NULL on a zeroed/failed b. */
    double *tab;
    int tab_n;
    double t_lo;     /* earliest tabulated time (~-400 M) */
} StcBbh;

/* Instantaneous state at coordinate time t (all ILLUSTRATIVE after t_m). */
typedef struct {
    int stage;       /* STC_BBH_* */
    double a;        /* separation, M (0 once coalesced) */
    double phi;      /* orbital phase, rad (GW phase = 2 phi) */
    double omega;    /* orbital angular frequency, rad/M */
    double amp;      /* strain scale r h (dimensionless, times M/r) */
    double mass;     /* horizon mass budget, M (1 -> mf across plunge) */
    double fuse;     /* 0 two holes .. 1 single remnant (plunge blend) */
    double x1[2];    /* hole 1 position in the orbital plane, M */
    double x2[2];    /* hole 2 position in the orbital plane, M */
} StcBbhState;

/* Configure a system. 1 <= q = m1/m2 <= 1e6, a_merge > 0, a0 > a_merge,
 * a0 > 6 (the inspiral must start outside the ISCO). The a_merge argument
 * is retained (and validated) for signature compatibility, but the PN
 * inspiral now ends at the ISCO, so the plunge always starts at
 * a = 6 M; b->a_merge is set to STC_BBH_A_ISCO. Returns 0 on success,
 * -1 on bad arguments or if the inspiral span exceeds the table cap
 * (b zeroed). On success b owns a heap table: pair with stc_bbh_free,
 * and free before re-initing the same struct. */
int stc_bbh_init(StcBbh *b, double q, double a0, double a_merge);

/* Release the inspiral table owned by b. NULL- and zero-safe. */
void stc_bbh_free(StcBbh *b);

/* Evaluate the timeline at t (any real t; t < 0 extrapolates the
 * inspiral backwards via the Hermite table, which is what retarded-time
 * lookups need; below t_lo the first-interval cubic continues).
 * NULL-safe (zeroed state). */
void stc_bbh_eval(const StcBbh *b, double t, StcBbhState *s);

/* Plus-polarization strain times r, as seen at azimuth az in the orbital
 * plane (inclination 90 deg), at retarded time t_ret = t - r:
 *   r h+ = amp(t_ret) cos(2 phi(t_ret) - 2 az) / 2.
 * Cheap path used per grid vertex. NULL-safe (0). */
double stc_bbh_rh_plus(const StcBbh *b, double t_ret, double az);

/* Face-on (inclination 0) plus-polarization r h+ = amp cos(2 phi), the
 * strain an observer on the orbital axis records. NULL-safe (0). */
double stc_bbh_rh_faceon(const StcBbh *b, double t_ret);

#endif /* STC_MERGER_H */
