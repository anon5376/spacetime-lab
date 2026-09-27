/* SPACETIME LAB — body-core: N-body sandbox physics core.
 * Method: velocity-Verlet gravity with Plummer softening, merge/bounce
 * collisions, vis-viva Kepler elements, seeded presets. SI doubles,
 * deterministic, no I/O. Depends only on stc/core (vec/units/rng) + stdlib.
 */

#include "stc/body.h"
#include "stc/rng.h"
#include "stc/units.h"

#include <math.h>
#include <stdio.h>
#include <string.h>

/* Method: sim init — zero state, store G/softening/restitution/merge. */
void stc_sim_init(StcSim *s, double G, double soft, double restitution,
                  int merge)
{
    if (s == NULL) {
        return;
    }
    memset(s, 0, sizeof *s);
    s->G = G;
    s->soft2 = (soft > 0.0 && isfinite(soft)) ? soft * soft : 0.0;
    if (!isfinite(restitution)) {
        restitution = 0.0;
    }
    s->restitution = stc_clamp(restitution, 0.0, 1.0);
    s->merge = merge ? 1 : 0;
    s->vmax = 0.0;
    s->t = 0.0;
    s->n = 0;
}

/* Method: spawn — append a body; index or -1 (full/invalid). */
int stc_sim_add(StcSim *s, const char *name, double m, double r, StcVec3 p,
                StcVec3 v, int fixed)
{
    StcBody *b;

    if (s == NULL) {
        return -1;
    }
    if (!(m > 0.0) || !isfinite(m)) {
        return -1;
    }
    if (!(r >= 0.0) || !isfinite(r)) {
        return -1;
    }
    if (!isfinite(p.x) || !isfinite(p.y) || !isfinite(p.z)) {
        return -1;
    }
    if (!isfinite(v.x) || !isfinite(v.y) || !isfinite(v.z)) {
        return -1;
    }
    if (s->n < 0) {
        s->n = 0;
    }
    if (s->n >= STC_BODY_MAX) {
        return -1;
    }
    if (name == NULL) {
        name = "body";
    }
    b = &s->b[s->n];
    snprintf(b->name, sizeof b->name, "%s", name);
    b->m = m;
    b->r = r;
    b->p = p;
    b->v = v;
    b->fixed = fixed ? 1 : 0;
    b->alive = 1;
    return s->n++;
}

/* Method: swap-with-last removal; compacts n. */
int stc_sim_remove(StcSim *s, int idx)
{
    if (s == NULL) {
        return -1;
    }
    if (idx < 0 || idx >= s->n) {
        return -1;
    }
    if (idx != s->n - 1) {
        s->b[idx] = s->b[s->n - 1];
    }
    s->n--;
    memset(&s->b[s->n], 0, sizeof s->b[s->n]);
    return 0;
}

/* Method: Plummer-softened N-body acceleration of body idx.
 * Equation: a_i = G * sum_j m_j (p_j - p_i) / (|p_j - p_i|^2 + soft2)^(3/2).
 * Fixed bodies exert but are never integrated by the stepper. */
void stc_sim_accel(const StcSim *s, int idx, StcVec3 *a)
{
    StcVec3 acc = { 0.0, 0.0, 0.0 };
    const StcBody *bi;
    int j;

    if (a == NULL) {
        return;
    }
    *a = acc;
    if (s == NULL || idx < 0 || idx >= s->n) {
        return;
    }
    bi = &s->b[idx];
    if (!bi->alive) {
        return;
    }
    for (j = 0; j < s->n; j++) {
        const StcBody *bj = &s->b[j];
        double dx, dy, dz, d2, inv, inv3, f;

        if (j == idx || !bj->alive) {
            continue;
        }
        dx = bj->p.x - bi->p.x;
        dy = bj->p.y - bi->p.y;
        dz = bj->p.z - bi->p.z;
        d2 = dx * dx + dy * dy + dz * dz + s->soft2;
        if (!(d2 > 0.0) || !isfinite(d2)) {
            continue; /* coincident with zero softening: skip singular pair */
        }
        inv = 1.0 / sqrt(d2);
        inv3 = inv * inv * inv;
        f = s->G * bj->m * inv3;
        acc.x += f * dx;
        acc.y += f * dy;
        acc.z += f * dz;
    }
    *a = acc;
}

/* Method: pairwise collision pass (dist < r1 + r2).
 * Merge mode: momentum-conserving merger, lower index survives, radius =
 * cbrt of volume sum; if either body is fixed the survivor is fixed at the
 * fixed position with v = 0. Bounce mode: impulse along the line of centers
 * with coefficient of restitution (fixed = infinite mass). */
static void stc_body_collide(StcSim *s)
{
    int i, j;

    if (s->merge) {
        for (i = 0; i < s->n; i++) {
            if (!s->b[i].alive) {
                continue;
            }
            for (j = i + 1; j < s->n; j++) {
                double m1, m2, M, r3, rad;
                double rr, d;
                StcVec3 p, v;
                int fix;

                if (!s->b[j].alive) {
                    continue;
                }
                rr = s->b[i].r + s->b[j].r;
                d = stc_vec3_dist(s->b[i].p, s->b[j].p);
                if (!(d < rr)) {
                    continue;
                }
                m1 = s->b[i].m;
                m2 = s->b[j].m;
                M = m1 + m2;
                p.x = (m1 * s->b[i].p.x + m2 * s->b[j].p.x) / M;
                p.y = (m1 * s->b[i].p.y + m2 * s->b[j].p.y) / M;
                p.z = (m1 * s->b[i].p.z + m2 * s->b[j].p.z) / M;
                v.x = (m1 * s->b[i].v.x + m2 * s->b[j].v.x) / M;
                v.y = (m1 * s->b[i].v.y + m2 * s->b[j].v.y) / M;
                v.z = (m1 * s->b[i].v.z + m2 * s->b[j].v.z) / M;
                r3 = s->b[i].r * s->b[i].r * s->b[i].r +
                     s->b[j].r * s->b[j].r * s->b[j].r;
                rad = cbrt(r3);
                fix = s->b[i].fixed || s->b[j].fixed;
                if (fix) {
                    p = s->b[i].fixed ? s->b[i].p : s->b[j].p;
                    v.x = 0.0;
                    v.y = 0.0;
                    v.z = 0.0;
                }
                s->b[i].m = M;
                s->b[i].p = p;
                s->b[i].v = v;
                s->b[i].r = rad;
                s->b[i].fixed = fix;
                /* Lower index survives; swap-with-last removal of j > i
                 * keeps slot i stable. Recheck the swapped-in slot j. */
                stc_sim_remove(s, j);
                j--;
            }
        }
        return;
    }

    {
        double e = stc_clamp(s->restitution, 0.0, 1.0);

        for (i = 0; i < s->n; i++) {
            if (!s->b[i].alive) {
                continue;
            }
            for (j = i + 1; j < s->n; j++) {
                StcVec3 n, dv;
                double rr, d, vn, m1, m2;

                if (!s->b[j].alive) {
                    continue;
                }
                rr = s->b[i].r + s->b[j].r;
                d = stc_vec3_dist(s->b[i].p, s->b[j].p);
                if (!(d < rr)) {
                    continue;
                }
                if (d > 0.0) {
                    n = stc_vec3_sub(s->b[j].p, s->b[i].p);
                    n = stc_vec3_scale(n, 1.0 / d);
                } else {
                    /* Coincident centers: arbitrary collision normal. */
                    n.x = 1.0;
                    n.y = 0.0;
                    n.z = 0.0;
                }
                dv = stc_vec3_sub(s->b[j].v, s->b[i].v);
                vn = stc_vec3_dot(dv, n);
                if (!(vn < 0.0)) {
                    continue; /* separating or resting: no impulse */
                }
                if (s->b[i].fixed && s->b[j].fixed) {
                    continue;
                } else if (s->b[i].fixed) {
                    /* m_i = infinity: reflect j off a static wall. */
                    s->b[j].v = stc_vec3_sub(
                        s->b[j].v, stc_vec3_scale(n, (1.0 + e) * vn));
                } else if (s->b[j].fixed) {
                    s->b[i].v = stc_vec3_add(
                        s->b[i].v, stc_vec3_scale(n, (1.0 + e) * vn));
                } else {
                    /* 1-D restitution impulse along line of centers. */
                    double imp;
                    m1 = s->b[i].m;
                    m2 = s->b[j].m;
                    imp = -(1.0 + e) * vn / (1.0 / m1 + 1.0 / m2);
                    s->b[i].v = stc_vec3_sub(
                        s->b[i].v, stc_vec3_scale(n, imp / m1));
                    s->b[j].v = stc_vec3_add(
                        s->b[j].v, stc_vec3_scale(n, imp / m2));
                }
            }
        }
    }
}

/* Method: velocity Verlet (kick-drift-kick) + collisions + speed cap.
 * Fixed bodies never move; dead slots skipped; t advances by dt. */
void stc_sim_step(StcSim *s, double dt)
{
    StcVec3 a0[STC_BODY_MAX];
    StcVec3 a1[STC_BODY_MAX];
    StcVec3 zero = { 0.0, 0.0, 0.0 };
    double h;
    int i;

    if (s == NULL || !isfinite(dt)) {
        return;
    }
    if (s->n < 0) {
        s->n = 0;
    }
    if (s->n > STC_BODY_MAX) {
        s->n = STC_BODY_MAX;
    }
    h = 0.5 * dt;

    /* Kick (first half) with accelerations at current positions. */
    for (i = 0; i < s->n; i++) {
        if (s->b[i].alive && !s->b[i].fixed) {
            stc_sim_accel(s, i, &a0[i]);
        } else {
            a0[i] = zero;
        }
    }
    for (i = 0; i < s->n; i++) {
        if (s->b[i].alive && !s->b[i].fixed) {
            s->b[i].v = stc_vec3_add(s->b[i].v, stc_vec3_scale(a0[i], h));
        }
    }
    /* Drift to new positions. */
    for (i = 0; i < s->n; i++) {
        if (s->b[i].alive && !s->b[i].fixed) {
            s->b[i].p = stc_vec3_add(s->b[i].p, stc_vec3_scale(s->b[i].v, dt));
        }
    }
    /* Kick (second half) with accelerations at new positions. */
    for (i = 0; i < s->n; i++) {
        if (s->b[i].alive && !s->b[i].fixed) {
            stc_sim_accel(s, i, &a1[i]);
        } else {
            a1[i] = zero;
        }
    }
    for (i = 0; i < s->n; i++) {
        if (s->b[i].alive && !s->b[i].fixed) {
            s->b[i].v = stc_vec3_add(s->b[i].v, stc_vec3_scale(a1[i], h));
        }
    }

    stc_sim_poststep(s);

    s->t += dt;
}

/* Method: post-integration pass shared by all integrators — the pairwise
 * collision pass (merge or restitution bounce; merge mode compacts via
 * swap-with-last removal) followed by the vmax speed cap. This is exactly
 * the tail stc_sim_step runs after its velocity-Verlet update, exported so
 * stc_sim_step_rk4 callers can run the identical collision/cap policy.
 * Advances no time: the owning stepper owns s->t. */
void stc_sim_poststep(StcSim *s)
{
    int i;

    if (s == NULL) {
        return;
    }
    if (s->n < 0) {
        s->n = 0;
    }
    if (s->n > STC_BODY_MAX) {
        s->n = STC_BODY_MAX;
    }

    stc_body_collide(s);

    /* Speed cap (vmax > 0): clamp every live body post-step. */
    if (s->vmax > 0.0 && isfinite(s->vmax)) {
        for (i = 0; i < s->n; i++) {
            double sp;

            if (!s->b[i].alive) {
                continue;
            }
            sp = stc_vec3_len(s->b[i].v);
            if (sp > s->vmax && sp > 0.0) {
                s->b[i].v = stc_vec3_scale(s->b[i].v, s->vmax / sp);
            }
        }
    }
}

/* Method: total energy — kinetic (free bodies) + pairwise Plummer potential.
 * Equation: E = sum 1/2 m v^2 + sum_{i<j} -G m_i m_j / sqrt(r^2 + soft2). */
double stc_sim_energy(const StcSim *s)
{
    double ke = 0.0, pe = 0.0;
    int i, j;

    if (s == NULL) {
        return 0.0;
    }
    for (i = 0; i < s->n; i++) {
        if (s->b[i].alive && !s->b[i].fixed) {
            double v2 = stc_vec3_dot(s->b[i].v, s->b[i].v);
            ke += 0.5 * s->b[i].m * v2;
        }
    }
    for (i = 0; i < s->n; i++) {
        for (j = i + 1; j < s->n; j++) {
            double d2;

            if (!s->b[i].alive || !s->b[j].alive) {
                continue;
            }
            d2 = stc_vec3_dist(s->b[i].p, s->b[j].p);
            d2 = d2 * d2 + s->soft2;
            if (d2 > 0.0 && isfinite(d2)) {
                pe += -s->G * s->b[i].m * s->b[j].m / sqrt(d2);
            }
        }
    }
    return ke + pe;
}

/* Method: vis-viva Kepler elements of idx about center.
 * Equations: vis-viva 1/a = 2/r - v^2/mu, eccentricity vector
 * e_vec = ((v^2 - mu/r) r_vec - (r_vec . v_vec) v_vec)/mu, Kepler period
 * T = 2 pi sqrt(a^3/mu). Period is 0 when unbound or invalid. */
void stc_sim_elements(const StcSim *s, int idx, int center, double *a,
                      double *e, double *period)
{
    double la = 0.0, le = 0.0, lp = 0.0;

    if (s != NULL && idx >= 0 && idx < s->n && center >= 0 && center < s->n &&
        idx != center && s->b[idx].alive && s->b[center].alive) {
        double mu = s->G * (s->b[center].m + s->b[idx].m);
        StcVec3 rp = stc_vec3_sub(s->b[idx].p, s->b[center].p);
        StcVec3 vv = stc_vec3_sub(s->b[idx].v, s->b[center].v);
        double r = stc_vec3_len(rp);
        double v2 = stc_vec3_dot(vv, vv);

        if (mu > 0.0 && isfinite(mu) && r > 0.0 && isfinite(r) &&
            isfinite(v2)) {
            double inv_a = 2.0 / r - v2 / mu;

            if (isfinite(inv_a)) {
                /* Eccentricity vector (Laplace-Runge-Lenz / mu). */
                double rdotv = stc_vec3_dot(rp, vv);
                double c1 = v2 - mu / r;
                StcVec3 ev;

                ev.x = (c1 * rp.x - rdotv * vv.x) / mu;
                ev.y = (c1 * rp.y - rdotv * vv.y) / mu;
                ev.z = (c1 * rp.z - rdotv * vv.z) / mu;
                le = stc_vec3_len(ev);
                if (!isfinite(le)) {
                    le = 0.0;
                }
                if (inv_a > 0.0) {
                    la = 1.0 / inv_a;
                    lp = 2.0 * STC_PI * sqrt(la * la * la / mu);
                    if (!isfinite(lp)) {
                        lp = 0.0;
                    }
                } else if (inv_a < 0.0) {
                    la = 1.0 / inv_a; /* hyperbolic: negative a, no period */
                } else {
                    la = INFINITY; /* parabolic: unbounded, no period */
                }
            }
        }
    }
    if (a != NULL) {
        *a = la;
    }
    if (e != NULL) {
        *e = le;
    }
    if (period != NULL) {
        *period = lp;
    }
}

/* Seeded presets in the XZ plane (y-up), prograde toward +Z. */

int stc_preset_count(void)
{
    return 4;
}

const char *stc_preset_name(int i)
{
    static const char *names[4] = { "solar", "binary", "hole_disk",
                                    "slingshot" };

    if (i < 0 || i >= 4) {
        return NULL;
    }
    return names[i];
}

/* Real-planet table (SI): name, mass kg, radius m, semi-major axis m. */
struct stc_planet_row {
    const char *name;
    double m;
    double r;
    double a;
};

static const struct stc_planet_row stc_solar_planets[8] = {
    { "Mercury", 3.3011e23, 2.4397e6, 5.7909e10 },
    { "Venus", 4.8675e24, 6.0518e6, 1.0821e11 },
    { "Earth", STC_M_EARTH, STC_R_EARTH, STC_AU },
    { "Mars", 6.4171e23, 3.3895e6, 2.2794382e11 },
    { "Jupiter", 1.8982e27, 6.9911e7, 7.785472e11 },
    { "Saturn", 5.6834e26, 5.8232e7, 1.43344937e12 },
    { "Uranus", 8.6810e25, 2.5362e7, 2.872463e12 },
    { "Neptune", 1.02413e26, 2.4622e7, 4.49506e12 },
};

/* Preset 0: solar — free Sun + 8 planets on circular orbits; the Sun's
 * velocity momentum-compensates the planets so total momentum is ~0. */
static int stc_preset_solar(StcSim *s)
{
    const double msun = STC_M_SUN;
    StcVec3 origin = { 0.0, 0.0, 0.0 };
    double px = 0.0, py = 0.0, pz = 0.0;
    int i;

    stc_sim_init(s, STC_G, 0.0, 0.5, 1);
    if (stc_sim_add(s, "Sun", msun, STC_R_SUN, origin, origin, 0) < 0) {
        return -1;
    }
    for (i = 0; i < 8; i++) {
        StcVec3 p = { stc_solar_planets[i].a, 0.0, 0.0 };
        StcVec3 v;
        double vc;

        /* Circular speed about the Sun. */
        vc = sqrt(STC_G * (msun + stc_solar_planets[i].m) /
                  stc_solar_planets[i].a);
        v.x = 0.0;
        v.y = 0.0;
        v.z = vc;
        px += stc_solar_planets[i].m * v.x;
        py += stc_solar_planets[i].m * v.y;
        pz += stc_solar_planets[i].m * v.z;
        if (stc_sim_add(s, stc_solar_planets[i].name, stc_solar_planets[i].m,
                        stc_solar_planets[i].r, p, v, 0) < 0) {
            return -1;
        }
    }
    s->b[0].v.x = -px / msun;
    s->b[0].v.y = -py / msun;
    s->b[0].v.z = -pz / msun;
    return 0;
}

/* Preset 1: binary — two 1e30 stars at +-1e10 with v = sqrt(G*M/(4*d)) plus
 * a 6e24 planet at 1e11 circular about the total mass. */
static int stc_preset_binary(StcSim *s)
{
    const double mstar = 1.0e30;
    const double d = 1.0e10;
    const double mplanet = 6.0e24;
    const double ap = 1.0e11;
    double vstar, vplanet;
    StcVec3 pa = { d, 0.0, 0.0 }, va;
    StcVec3 pb = { -d, 0.0, 0.0 }, vb;
    StcVec3 pp = { ap, 0.0, 0.0 }, vp;

    stc_sim_init(s, STC_G, 0.0, 0.5, 1);
    vstar = sqrt(STC_G * mstar / (4.0 * d));
    va.x = 0.0;
    va.y = 0.0;
    va.z = vstar;
    vb.x = 0.0;
    vb.y = 0.0;
    vb.z = -vstar;
    if (stc_sim_add(s, "starA", mstar, STC_R_SUN, pa, va, 0) < 0) {
        return -1;
    }
    if (stc_sim_add(s, "starB", mstar, STC_R_SUN, pb, vb, 0) < 0) {
        return -1;
    }
    vplanet = sqrt(STC_G * (2.0 * mstar + mplanet) / ap);
    vp.x = 0.0;
    vp.y = 0.0;
    vp.z = vplanet;
    if (stc_sim_add(s, "planet", mplanet, STC_R_EARTH, pp, vp, 0) < 0) {
        return -1;
    }
    return 0;
}

/* Preset 2: hole_disk — FIXED 10 Msun hole + 24 seeded ring particles on
 * circular orbits (seeded radii 1e6..3e6 m) + one eccentric probe. */
static int stc_preset_holedisk(StcSim *s, uint64_t seed)
{
    const double mhole = 10.0 * STC_M_SUN;
    StcVec3 origin = { 0.0, 0.0, 0.0 };
    StcRng64 rng;
    char name[32];
    int i;

    stc_sim_init(s, STC_G, 0.0, 0.5, 1);
    /* Collision radius ~ Schwarzschild radius of a 10 Msun hole. */
    if (stc_sim_add(s, "hole", mhole, 3.0e4, origin, origin, 1) < 0) {
        return -1;
    }
    stc_rng64_seed(&rng, seed);
    for (i = 0; i < 24; i++) {
        double rad = stc_rng64_uniform_range(&rng, 1.0e6, 3.0e6);
        double th = stc_rng64_uniform(&rng) * 2.0 * STC_PI;
        double vc = sqrt(STC_G * mhole / rad);
        StcVec3 p = { rad * cos(th), 0.0, rad * sin(th) };
        /* Prograde circular tangent (-sin th, 0, cos th) * vc. */
        StcVec3 v = { -sin(th) * vc, 0.0, cos(th) * vc };

        snprintf(name, sizeof name, "d%02d", i);
        if (stc_sim_add(s, name, 1.0e15, 1.0e3, p, v, 0) < 0) {
            return -1;
        }
    }
    {
        /* Eccentric probe: starts at 2.5e6 m with 55% of circular speed. */
        double rp = 2.5e6;
        double vc = sqrt(STC_G * mhole / rp);
        StcVec3 p = { rp, 0.0, 0.0 };
        StcVec3 v = { 0.0, 0.0, 0.55 * vc };

        if (stc_sim_add(s, "probe", 1.0e3, 1.0e3, p, v, 0) < 0) {
            return -1;
        }
    }
    return 0;
}

/* Preset 3: slingshot — FIXED Jupiter-mass planet at origin + probe.
 * Restricted model: the planet has infinite inertia (fixed), so the probe
 * feels its gravity but does not perturb it; total momentum of the free
 * bodies alone is not conserved. */
static int stc_preset_slingshot(StcSim *s)
{
    StcVec3 origin = { 0.0, 0.0, 0.0 };
    StcVec3 p = { -6.0e9, 0.0, 3.0e9 };
    StcVec3 v = { 1.0e4, 0.0, 0.0 };

    stc_sim_init(s, STC_G, 0.0, 0.5, 1);
    if (stc_sim_add(s, "jupiter", 1.8982e27, 6.9911e7, origin, origin, 1) <
        0) {
        return -1;
    }
    if (stc_sim_add(s, "probe", 1.0e3, 1.0e1, p, v, 0) < 0) {
        return -1;
    }
    return 0;
}

int stc_preset_load(StcSim *s, int i, uint64_t seed)
{
    if (s == NULL) {
        return -1;
    }
    switch (i) {
    case 0:
        return stc_preset_solar(s);
    case 1:
        return stc_preset_binary(s);
    case 2:
        return stc_preset_holedisk(s, seed);
    case 3:
        return stc_preset_slingshot(s);
    default:
        return -1;
    }
}

/* Method: gravitational redshift — proper rate sqrt(1 - rs/r), clamped
 * to [0,1]; rs = 2GM/c^2. Returns 0 at or inside the horizon. */
double stc_body_proper_rate(double r, double M)
{
    double rs, q;

    if (!(r > 0.0) || !isfinite(r)) {
        return 0.0;
    }
    if (!(M > 0.0) || !isfinite(M)) {
        return 1.0;
    }
    rs = 2.0 * STC_G * M / (STC_C * STC_C);
    if (!(rs > 0.0) || !isfinite(rs)) {
        return 1.0;
    }
    if (r <= rs) {
        return 0.0;
    }
    q = 1.0 - rs / r;
    if (q <= 0.0) {
        return 0.0;
    }
    if (q >= 1.0) {
        return 1.0;
    }
    return sqrt(q);
}

/* Method: Lorentz factor gamma = 1/sqrt(1-(v/c)^2); -1 when v >= c. */
double stc_gamma(double v)
{
    double b;

    if (!isfinite(v)) {
        return -1.0;
    }
    b = fabs(v) / STC_C;
    if (b >= 1.0) {
        return -1.0;
    }
    return 1.0 / sqrt(1.0 - b * b);
}
