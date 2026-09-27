/* SPACETIME LAB — scenario presets: deterministic demo table implementation.
 * Self-contained analytic steppers (no sibling physics includes).
 */

#include "stc/scenario.h"

#include <math.h>
#include <stdio.h>
#include <string.h>

#include "stc/units.h"

/* --- Preset table (CONTRACT Sec 5 names; fixed seeds; fixed frames) --- */

static const StcScenarioInfo k_scenarios[STC_SCENARIO_COUNT] = {
    {"twin", NULL, "Twin paradox out-and-back at 0.8c",
     "Time dilation d_tau = dt / gamma", 20260901ULL, 100,
     0.1 * STC_YEAR},
    {"pole", "ladder", "Barn-pole paradox at gamma = 2",
     "Length contraction L = L0 / gamma", 20260902ULL, 60, 2.0e-9},
    {"doppler", "minkowski", "Receding source at 0.5c, flat spacetime",
     "Relativistic Doppler 1+z = sqrt((1+b)/(1-b))", 20260903ULL, 80,
     1.0e-9},
    {"mercury", NULL, "Mercury perihelion advance, 415 orbits",
     "Einstein precession dphi = 6*pi*G*M/(a*(1-e^2)*c^2)", 20260904ULL,
     415, 87.969 * STC_DAY},
    {"gps", NULL, "GPS clock SR+GR rate offsets over one day",
     "Proper rate sqrt(1-rs/r) + SR gamma", 20260905ULL, 100,
     STC_DAY / 100.0},
    {"bend", "photon", "Sun-grazing photon deflection",
     "Weak deflection alpha = 4*G*M/(c^2*b)", 20260906ULL, 64,
     16.0 * STC_R_SUN / STC_C / 64.0},
    {"ring", "lensing", "Point-lens images sweeping past alignment",
     "Thin-lens beta = theta - theta_E^2/theta", 20260907ULL, 40, 1.0},
    {"flamm", NULL, "Flamm paraboloid radial scan, 10 Msun",
     "Flamm w = 2*sqrt(rs*(r-rs))", 20260908ULL, 32, 1.0},
    {"chirp", "waves", "1.4+1.4 Msun inspiral chirp + ring response",
     "Quadrupole chirp f(tau) ~ tau^(-3/8)", 20260909ULL, 200,
     19.5 / 200.0},
    {"shadow", "kerr", "Kerr spin sweep: ISCO, horizon, shadow",
     "Bardeen ISCO + horizon r+ = M + sqrt(M^2-a^2)", 20260910ULL, 41,
     1.0},
};

/* --- Small local helpers (closed forms; no iteration) --- */

/* Lorentz factor gamma = 1/sqrt(1-b^2); |b|>=1 clamps to gamma 1. */
static double s_gamma(double beta) {
    double b2;
    if (beta >= 1.0 || beta <= -1.0) {
        return 1.0;
    }
    b2 = beta * beta;
    return 1.0 / sqrt(1.0 - b2);
}

/* Arcseconds per radian, from pi (matches stc_arcsec reciprocal). */
static double s_arcsec_per_rad(void) {
    return 180.0 * 3600.0 / STC_PI;
}

/* Bardeen prograde (coordinate Lz>0) ISCO radius in meters.
 * M_geom in meters, chi = a/M in [-1,1]. Hits 6M/1M/9M at 0/+1/-1. */
static double s_isco_prograde(double m_geom, double chi) {
    double c, t0, t1, t2, z1, z2, root, sgn;
    if (!(m_geom > 0.0)) {
        return 0.0;
    }
    c = stc_clamp(chi, -1.0, 1.0);
    t0 = cbrt(1.0 - c * c);
    t1 = cbrt(1.0 + c);
    t2 = cbrt(1.0 - c);
    z1 = 1.0 + t0 * (t1 + t2);
    z2 = sqrt(3.0 * c * c + z1 * z1);
    root = sqrt((3.0 - z1) * (3.0 + z1 + 2.0 * z2));
    sgn = (c < 0.0) ? 1.0 : -1.0; /* coordinate-prograde branch */
    return m_geom * (3.0 + z2 + sgn * root);
}

/* --- Per-scenario init/step pairs; state is a pure function of frame --- */

static void s_init_twin(StcScenarioCtx *c) {
    double v = 0.8 * STC_C;
    c->obs[2] = s_gamma(0.8);
    c->obs[3] = v;
}

/* Time dilation d_tau = dt/gamma on each leg; turnaround instantaneous. */
static void s_step_twin(StcScenarioCtx *c) {
    double leg = 4.0 * STC_C * STC_YEAR; /* 4 light-years, m */
    double half_t = 5.0 * STC_YEAR;      /* 5 yr per leg, s */
    double v = c->obs[3];
    double g = c->obs[2];
    double x;
    c->obs[0] = c->t;
    c->obs[1] = c->t / g;
    if (c->t <= half_t) {
        x = v * c->t;
        c->vel[0].x = v;
    } else {
        x = leg - v * (c->t - half_t);
        c->vel[0].x = -v;
    }
    c->pos[0].x = x;
    c->pos[0].y = 0.0;
    c->pos[0].z = 0.0;
    c->pos[1].x = 0.0;
    c->pos[1].y = 0.0;
    c->pos[1].z = 0.0;
}

static void s_init_pole(StcScenarioCtx *c) {
    double beta = 0.8660254037844386; /* sqrt(3)/2, gamma = 2 */
    c->obs[1] = 20.0;                 /* ladder rest length, m */
    c->obs[2] = s_gamma(beta);
    c->obs[3] = 10.0; /* barn length, m */
    c->obs[0] = c->obs[1] / c->obs[2];
    c->vel[0].x = beta * STC_C;
    c->vel[0].y = 0.0;
    c->vel[0].z = 0.0;
}

/* Length contraction L = L0/gamma; ladder center sweeps past the barn. */
static void s_step_pole(StcScenarioCtx *c) {
    double v = c->vel[0].x;
    double xc = -16.0 + v * c->t;
    double half = 0.5 * c->obs[0];
    c->pos[0].x = xc + half;
    c->pos[0].y = 0.0;
    c->pos[0].z = 0.0;
    c->pos[1].x = xc - half;
    c->pos[1].y = 0.0;
    c->pos[1].z = 0.0;
}

static void s_init_doppler(StcScenarioCtx *c) {
    double beta = 0.5;
    double z = sqrt((1.0 + beta) / (1.0 - beta)) - 1.0;
    c->obs[3] = beta;
    c->obs[0] = z;
    c->obs[2] = 5.0e14; /* emitted frequency, Hz */
    c->obs[1] = c->obs[2] / (1.0 + z);
    c->vel[0].x = beta * STC_C;
    c->vel[0].y = 0.0;
    c->vel[0].z = 0.0;
}

/* Relativistic Doppler: received pulse rate = emitted/(1+z). */
static void s_step_doppler(StcScenarioCtx *c) {
    double d0 = 1.0e6;
    c->pos[0].x = d0 + c->vel[0].x * c->t;
    c->pos[0].y = 0.0;
    c->pos[0].z = 0.0;
    c->obs[4] = c->obs[2] * c->t;
    c->obs[5] = c->obs[1] * c->t;
    c->phase = fmod(2.0 * STC_PI * c->obs[1] * c->t, 2.0 * STC_PI);
}

static void s_init_mercury(StcScenarioCtx *c) {
    c->obs[2] = 5.7909e10; /* semi-major axis, m */
    c->obs[3] = 0.2056;    /* eccentricity */
}

/* Einstein perihelion precession accumulated orbit by orbit. */
static void s_step_mercury(StcScenarioCtx *c) {
    double a = c->obs[2];
    double e = c->obs[3];
    double gm = STC_G * STC_M_SUN;
    double c2 = STC_C * STC_C;
    double per_orbit = 6.0 * STC_PI * gm / (a * (1.0 - e * e) * c2);
    double acc = per_orbit * (double)c->frame;
    double rp = a * (1.0 - e);
    c->obs[1] = (double)c->frame;
    c->obs[0] = acc * s_arcsec_per_rad();
    c->phase = fmod(2.0 * STC_PI * c->t / (87.969 * STC_DAY), 2.0 * STC_PI);
    c->pos[0].x = rp * cos(acc);
    c->pos[0].y = rp * sin(acc);
    c->pos[0].z = 0.0;
}

static void s_init_gps(StcScenarioCtx *c) {
    double r = 2.656e7; /* GPS orbital radius, m */
    double re = STC_R_EARTH; /* Earth mean radius, m */
    double gm = STC_GM_EARTH; /* GM_earth, m^3/s^2 (units.h; was WGS84 literal) */
    double c2 = STC_C * STC_C;
    double v2 = gm / r;
    double sr = -v2 / (2.0 * c2);
    double gr = (gm / c2) * (1.0 / re - 1.0 / r);
    c->obs[0] = sr * STC_DAY * 1.0e6; /* SR offset, us/day */
    c->obs[1] = gr * STC_DAY * 1.0e6; /* GR offset, us/day */
    c->obs[2] = (sr + gr) * STC_DAY * 1.0e6;
    c->pos[1].x = r;
    c->pos[1].y = 0.0;
    c->pos[1].z = 0.0;
}

/* Proper rate sqrt(1-rs/r) vs SR gamma: net offset integrated over a day. */
static void s_step_gps(StcScenarioCtx *c) {
    double r = 2.656e7;
    double gm = STC_GM_EARTH; /* GM_earth, m^3/s^2 */
    double period = 2.0 * STC_PI * sqrt(r * r * r / gm);
    double th = 2.0 * STC_PI * c->t / period;
    double net_frac = c->obs[2] / (STC_DAY * 1.0e6);
    c->phase = fmod(th, 2.0 * STC_PI);
    c->pos[0].x = r * cos(th);
    c->pos[0].y = r * sin(th);
    c->pos[0].z = 0.0;
    c->obs[3] = net_frac * c->t * 1.0e9; /* accumulated, ns */
}

static void s_init_bend(StcScenarioCtx *c) {
    double u = stc_rng64_uniform(&c->rng);
    double b = STC_R_SUN * (1.0 + 0.05 * u);
    double gm = STC_GM_SUN;
    double alpha = 4.0 * gm / (STC_C * STC_C * b);
    c->obs[1] = b;
    c->obs[0] = alpha * s_arcsec_per_rad();
}

/* Weak deflection alpha = 4GM/(c^2 b) as a small-angle kink at the limb. */
static void s_step_bend(StcScenarioCtx *c) {
    double b = c->obs[1];
    double alpha = c->obs[0] / s_arcsec_per_rad();
    double span = 16.0 * STC_R_SUN;
    double x = -0.5 * span + STC_C * c->t;
    double y = (x > 0.0) ? (b - alpha * x) : b;
    c->obs[2] = x;
    c->obs[3] = y;
    c->pos[0].x = x;
    c->pos[0].y = y;
    c->pos[0].z = 0.0;
    c->vel[0].x = STC_C;
    c->vel[0].y = 0.0;
    c->vel[0].z = 0.0;
}

static void s_init_ring(StcScenarioCtx *c) {
    double dl = 3.0857e25; /* lens distance, m (~1 Gpc) */
    double ds = 6.1714e25; /* source distance, m (~2 Gpc) */
    double gm = STC_G * 1.0e12 * STC_M_SUN;
    double c2 = STC_C * STC_C;
    double dls = ds - dl;
    double te = sqrt(4.0 * gm / c2 * dls / (dl * ds));
    double u = stc_rng64_uniform(&c->rng);
    c->obs[0] = te * s_arcsec_per_rad();
    c->obs[4] = (u - 0.5) * 0.2 * c->obs[0]; /* seeded offset, arcsec */
    c->pos[1].x = dl;
    c->pos[1].y = 0.0;
    c->pos[1].z = 0.0;
}

/* Thin-lens equation: theta_pm = (beta +- sqrt(beta^2+4 te^2))/2. */
static void s_step_ring(StcScenarioCtx *c) {
    double te = c->obs[0];
    double frac = (c->frames > 1)
                      ? ((double)(c->frame - 1) / (double)(c->frames - 1))
                      : 0.0;
    double beta = -2.0 * te + 4.0 * te * frac + c->obs[4];
    double disc = sqrt(beta * beta + 4.0 * te * te);
    c->obs[3] = beta;
    c->obs[1] = 0.5 * (beta + disc);
    c->obs[2] = 0.5 * (beta - disc);
    c->pos[0].x = c->obs[1];
    c->pos[0].y = c->obs[2];
    c->pos[0].z = 0.0;
}

static void s_init_flamm(StcScenarioCtx *c) {
    int i;
    double rs = stc_schwarzschild_radius(10.0 * STC_M_SUN);
    double rmax = 8.0 * rs;
    c->obs[0] = rs;
    c->obs[3] = 2.0 * sqrt(rs * (rmax - rs));
    for (i = 0; i < 4; i++) {
        c->obs[4 + i] = 2.0 * STC_PI * stc_rng64_uniform(&c->rng);
    }
}

/* Flamm paraboloid w = 2*sqrt(rs*(r-rs)); radius sweeps throat to 8 rs. */
static void s_step_flamm(StcScenarioCtx *c) {
    int i;
    double rs = c->obs[0];
    double frac = (c->frames > 1)
                      ? ((double)(c->frame - 1) / (double)(c->frames - 1))
                      : 0.0;
    double r = rs * (1.05 + (8.0 - 1.05) * frac);
    double w = 2.0 * sqrt(rs * (r - rs));
    c->obs[2] = r;
    c->obs[1] = w;
    for (i = 0; i < 4; i++) {
        double ph = c->obs[4 + i];
        c->pos[i].x = r * cos(ph);
        c->pos[i].y = w;
        c->pos[i].z = r * sin(ph);
    }
}

static void s_init_chirp(StcScenarioCtx *c) {
    double m1 = 1.4 * STC_M_SUN;
    double m2 = 1.4 * STC_M_SUN;
    double mc = pow(m1 * m2, 0.6) / pow(m1 + m2, 0.2);
    c->obs[4] = mc;
    c->phase = 2.0 * STC_PI * stc_rng64_uniform(&c->rng);
}

/* Quadrupole chirp f(tau)=(1/pi)(5/256/tau)^(3/8)(G Mc/c^3)^(-5/8)
 * with TT plus/cross projection; phase integrated in fixed steps. */
static void s_step_chirp(StcScenarioCtx *c) {
    double dist = 3.0857e24; /* 100 Mpc, m */
    double mc = c->obs[4];
    double tcoal = 20.0 - c->t; /* time to coalescence, s */
    double tc;
    double f, h, ci = 0.5; /* cos(inclination), i = 60 deg */
    double hp, hx;
    if (tcoal < 0.25) {
        tcoal = 0.25;
    }
    tc = STC_G * mc / (STC_C * STC_C * STC_C);
    f = (1.0 / STC_PI) * pow(5.0 / 256.0 / tcoal, 0.375) * pow(tc, -0.625);
    h = (4.0 / dist) * pow(STC_G * mc / (STC_C * STC_C), 5.0 / 3.0) *
        pow(STC_PI * f / STC_C, 2.0 / 3.0);
    c->phase = fmod(c->phase + 2.0 * STC_PI * f * c->dt, 2.0 * STC_PI);
    hp = h * 0.5 * (1.0 + ci * ci) * cos(c->phase);
    hx = h * ci * sin(c->phase);
    c->obs[0] = f;
    c->obs[1] = h;
    c->obs[2] = hp;
    c->obs[3] = hx;
    c->pos[0].x = hp;
    c->pos[0].y = hx;
    c->pos[0].z = 0.0;
}

static void s_init_shadow(StcScenarioCtx *c) {
    (void)c;
}

/* Bardeen ISCO sweep over chi in [-1,1]; shadow estimate 3*sqrt(3)*M. */
static void s_step_shadow(StcScenarioCtx *c) {
    double mg = stc_geometric_M(10.0 * STC_M_SUN);
    double frac = (c->frames > 1)
                      ? ((double)(c->frame - 1) / (double)(c->frames - 1))
                      : 0.0;
    double chi = -1.0 + 2.0 * frac;
    double a = chi * mg;
    double rplus = mg + sqrt(mg * mg - a * a);
    double risco = s_isco_prograde(mg, chi);
    c->obs[0] = chi;
    c->obs[1] = risco;
    c->obs[2] = rplus;
    c->obs[3] = 3.0 * sqrt(3.0) * mg;
    c->pos[0].x = risco;
    c->pos[0].y = 0.0;
    c->pos[0].z = 0.0;
    c->pos[1].x = rplus;
    c->pos[1].y = 0.0;
    c->pos[1].z = 0.0;
}

/* --- Dispatch --- */

static void s_init_by_index(StcScenarioCtx *c) {
    switch (c->index) {
    case 0: s_init_twin(c); break;
    case 1: s_init_pole(c); break;
    case 2: s_init_doppler(c); break;
    case 3: s_init_mercury(c); break;
    case 4: s_init_gps(c); break;
    case 5: s_init_bend(c); break;
    case 6: s_init_ring(c); break;
    case 7: s_init_flamm(c); break;
    case 8: s_init_chirp(c); break;
    case 9: s_init_shadow(c); break;
    default: break;
    }
}

static void s_step_by_index(StcScenarioCtx *c) {
    switch (c->index) {
    case 0: s_step_twin(c); break;
    case 1: s_step_pole(c); break;
    case 2: s_step_doppler(c); break;
    case 3: s_step_mercury(c); break;
    case 4: s_step_gps(c); break;
    case 5: s_step_bend(c); break;
    case 6: s_step_ring(c); break;
    case 7: s_step_flamm(c); break;
    case 8: s_step_chirp(c); break;
    case 9: s_step_shadow(c); break;
    default: break;
    }
}

/* --- Public API --- */

size_t stc_scenario_count(void) {
    return STC_SCENARIO_COUNT;
}

const StcScenarioInfo *stc_scenario_list(void) {
    return k_scenarios;
}

const StcScenarioInfo *stc_scenario_find(const char *name) {
    size_t i;
    if (name == NULL) {
        return NULL;
    }
    for (i = 0; i < STC_SCENARIO_COUNT; i++) {
        if (strcmp(name, k_scenarios[i].name) == 0) {
            return &k_scenarios[i];
        }
        if (k_scenarios[i].alias != NULL &&
            strcmp(name, k_scenarios[i].alias) == 0) {
            return &k_scenarios[i];
        }
    }
    return NULL;
}

int stc_scenario_setup_with_seed(StcScenarioCtx *ctx, const char *name,
                                 uint64_t seed) {
    const StcScenarioInfo *info;
    size_t i;
    if (ctx == NULL) {
        return -1;
    }
    info = stc_scenario_find(name);
    if (info == NULL) {
        ctx->index = -1;
        return -1;
    }
    memset(ctx, 0, sizeof(*ctx));
    (void)snprintf(ctx->name, sizeof(ctx->name), "%s", info->name);
    ctx->seed = seed;
    ctx->frames = info->frames;
    ctx->dt = info->dt;
    ctx->index = -1;
    for (i = 0; i < STC_SCENARIO_COUNT; i++) {
        if (info == &k_scenarios[i]) {
            ctx->index = (int)i;
        }
    }
    stc_rng64_seed(&ctx->rng, seed);
    s_init_by_index(ctx);
    if (ctx->frames <= 0) {
        ctx->done = 1;
    }
    return 0;
}

int stc_scenario_setup(StcScenarioCtx *ctx, const char *name) {
    const StcScenarioInfo *info = stc_scenario_find(name);
    if (ctx == NULL || info == NULL) {
        if (ctx != NULL) {
            ctx->index = -1;
        }
        return -1;
    }
    return stc_scenario_setup_with_seed(ctx, name, info->seed);
}

int stc_scenario_step(StcScenarioCtx *ctx) {
    if (ctx == NULL || ctx->index < 0 ||
        ctx->index >= (int)STC_SCENARIO_COUNT) {
        return -1;
    }
    if (ctx->done || ctx->frame >= ctx->frames) {
        ctx->done = 1;
        return 1;
    }
    ctx->frame += 1;
    ctx->t = (double)ctx->frame * ctx->dt;
    s_step_by_index(ctx);
    if (ctx->frame >= ctx->frames) {
        ctx->done = 1;
        return 1;
    }
    return 0;
}

int stc_scenario_run(StcScenarioCtx *ctx, int max_frames) {
    int n = 0;
    int cap;
    if (ctx == NULL || ctx->index < 0 ||
        ctx->index >= (int)STC_SCENARIO_COUNT) {
        return -1;
    }
    cap = (max_frames < 0) ? ctx->frames : max_frames;
    while (n < cap && !ctx->done) {
        int rc = stc_scenario_step(ctx);
        if (rc < 0) {
            return -1;
        }
        n++;
        if (rc > 0) {
            break;
        }
    }
    return n;
}
