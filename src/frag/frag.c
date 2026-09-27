#include "stc/frag.h"

#include <math.h>
#include <stdint.h>
#include <stdio.h>

#include "stc/rng.h"
#include "stc/units.h"
#include "stc/vec.h"

/* Finite-vector check for input validation. */
static int stc_frag_finite_vec(StcVec3 v)
{
    return isfinite(v.x) && isfinite(v.y) && isfinite(v.z);
}

int stc_frag_impact(StcSim *s, int a, int b, uint64_t seed)
{
    /* Method: center-of-mass impact energy + COM-frame zero-net-kick
     * debris dispersal (see stc/frag.h for the equation index). */
    StcBody pa;
    StcBody pb;
    StcRng64 rng;
    StcVec3 pcm;
    StcVec3 vcm;
    StcVec3 vrel;
    StcVec3 mk;
    StcVec3 mo;
    double w[STC_BODY_MAX];
    StcVec3 kick[STC_BODY_MAX];
    StcVec3 off[STC_BODY_MAX];
    double m_tot;
    double mu;
    double vrel_mag;
    double e_imp;
    double rc;
    double gg;
    double v_esc;
    double sigma;
    double per;
    double scatter;
    double w_sum;
    double v_sum;
    char name[32];
    int count;
    int free_slots;
    int hi;
    int lo;
    int i;
    int spawned;

    /* Validate without mutating. */
    if (s == NULL) {
        return -1;
    }
    if (s->n < 2 || s->n > STC_BODY_MAX) {
        return -1;
    }
    if (a < 0 || a >= s->n || b < 0 || b >= s->n || a == b) {
        return -1;
    }
    pa = s->b[a];
    pb = s->b[b];
    if (!isfinite(pa.m) || !isfinite(pb.m) || pa.m <= 0.0 || pb.m <= 0.0) {
        return -1;
    }
    if (!isfinite(pa.r) || !isfinite(pb.r) || pa.r < 0.0 || pb.r < 0.0) {
        return -1;
    }
    if (!stc_frag_finite_vec(pa.p) || !stc_frag_finite_vec(pb.p) ||
        !stc_frag_finite_vec(pa.v) || !stc_frag_finite_vec(pb.v)) {
        return -1;
    }
    m_tot = pa.m + pb.m;
    if (!isfinite(m_tot) || m_tot <= 0.0) {
        return -1;
    }

    /* Equation: center of mass (position + velocity). Fixed parents
     * contribute their stored velocity (documented approximation). */
    pcm = stc_vec3_scale(stc_vec3_add(stc_vec3_scale(pa.p, pa.m),
                                      stc_vec3_scale(pb.p, pb.m)),
                         1.0 / m_tot);
    vcm = stc_vec3_scale(stc_vec3_add(stc_vec3_scale(pa.v, pa.m),
                                      stc_vec3_scale(pb.v, pb.m)),
                         1.0 / m_tot);
    if (!stc_frag_finite_vec(pcm) || !stc_frag_finite_vec(vcm)) {
        return -1;
    }

    /* Equation: center-of-mass impact energy E = 0.5*mu*vrel^2. */
    vrel = stc_vec3_sub(pa.v, pb.v);
    vrel_mag = stc_vec3_len(vrel);
    mu = pa.m * pb.m / m_tot;
    e_imp = 0.5 * mu * vrel_mag * vrel_mag;
    if (!isfinite(e_imp) || e_imp < 0.0) {
        return -1;
    }

    /* Debris count: min 2, log-scaled by impact energy, capped at the
     * free slots left after removing both parents. */
    count = STC_FRAG_MIN_DEBRIS +
        (int)(3.0 * log10(1.0 + e_imp / STC_FRAG_E_REF));
    free_slots = STC_BODY_MAX - (s->n - 2);
    if (free_slots < STC_FRAG_MIN_DEBRIS) {
        return -1;
    }
    if (count > free_slots) {
        count = free_slots;
    }

    /* Equation: mutual escape velocity at contact, sqrt(2*G*M/(ra+rb));
     * sets the debris dispersal scale together with |vrel|. */
    rc = pa.r + pb.r;
    if (!isfinite(rc) || rc < 0.0) {
        rc = 0.0;
    }
    gg = (s->G > 0.0 && isfinite(s->G)) ? s->G : STC_G;
    v_esc = 0.0;
    if (rc > 0.0) {
        v_esc = sqrt(2.0 * gg * m_tot / rc);
        if (!isfinite(v_esc) || v_esc < 0.0) {
            v_esc = 0.0;
        }
    }
    sigma = 0.10 * vrel_mag + 0.10 * v_esc;
    if (!isfinite(sigma) || sigma < 0.0) {
        sigma = 0.0;
    }
    per = sigma / sqrt(3.0);
    scatter = 0.5 * rc;

    /* Remove parents (larger index first: swap-with-last removal). */
    hi = (a > b) ? a : b;
    lo = (a > b) ? b : a;
    if (stc_sim_remove(s, hi) != 0) {
        return -1;
    }
    if (stc_sim_remove(s, lo) != 0) {
        return -1;
    }

    /* Deterministic stream: mass weights, then per-debris Gaussian
     * velocity kicks + uniform-in-sphere offsets, in fixed draw order. */
    stc_rng64_seed(&rng, seed);
    w_sum = 0.0;
    for (i = 0; i < count; i++) {
        w[i] = 0.5 + stc_rng64_uniform(&rng);
        w_sum += w[i];
    }
    for (i = 0; i < count; i++) {
        StcVec3 k;
        StcVec3 d;
        double dl;
        double rr;

        k.x = stc_rng64_normal(&rng) * per;
        k.y = stc_rng64_normal(&rng) * per;
        k.z = stc_rng64_normal(&rng) * per;
        kick[i] = k;
        d.x = stc_rng64_normal(&rng);
        d.y = stc_rng64_normal(&rng);
        d.z = stc_rng64_normal(&rng);
        dl = stc_vec3_len(d);
        if (dl > 0.0) {
            d = stc_vec3_scale(d, 1.0 / dl);
        } else {
            d.x = 1.0;
            d.y = 0.0;
            d.z = 0.0;
        }
        rr = scatter * cbrt(stc_rng64_uniform(&rng));
        off[i] = stc_vec3_scale(d, rr);
    }

    /* Equation: COM-frame zero-net-kick correction — subtract the
     * mass-weighted mean kick/offset so debris carry exactly M*vcm
     * (to roundoff) and sit exactly on pcm (to roundoff). */
    mk.x = 0.0;
    mk.y = 0.0;
    mk.z = 0.0;
    mo.x = 0.0;
    mo.y = 0.0;
    mo.z = 0.0;
    for (i = 0; i < count; i++) {
        double mi = m_tot * w[i] / w_sum;

        mk = stc_vec3_add(mk, stc_vec3_scale(kick[i], mi));
        mo = stc_vec3_add(mo, stc_vec3_scale(off[i], mi));
    }
    mk = stc_vec3_scale(mk, 1.0 / m_tot);
    mo = stc_vec3_scale(mo, 1.0 / m_tot);

    /* Equation: volume-conserving radii, r_i = cbrt(V_tot * m_i/M). */
    v_sum = pa.r * pa.r * pa.r + pb.r * pb.r * pb.r;
    if (!isfinite(v_sum) || v_sum < 0.0) {
        v_sum = 0.0;
    }
    spawned = 0;
    for (i = 0; i < count; i++) {
        double mi = m_tot * w[i] / w_sum;
        double ri = (v_sum > 0.0) ? cbrt(v_sum * (mi / m_tot)) : 0.0;
        StcVec3 pi = stc_vec3_add(pcm, stc_vec3_sub(off[i], mo));
        StcVec3 vi = stc_vec3_add(vcm, stc_vec3_sub(kick[i], mk));

        if (!isfinite(mi) || mi <= 0.0 || !isfinite(ri) || ri < 0.0 ||
            !stc_frag_finite_vec(pi) || !stc_frag_finite_vec(vi)) {
            break;
        }
        snprintf(name, sizeof name, "frag%d", i);
        if (stc_sim_add(s, name, mi, ri, pi, vi, 0) < 0) {
            break;
        }
        spawned++;
    }
    return spawned;
}
