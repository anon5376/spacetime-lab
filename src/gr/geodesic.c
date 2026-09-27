#include "stc/geodesic.h"

#include "stc/units.h"

#include <math.h>
#include <stddef.h>

static int stc_geo_bad(double x) {
    return !(x == x && x != 1.0 / 0.0 && x != -1.0 / 0.0);
}

static void stc_geo_reset(StcGeodesic *g, StcGeoKind kind) {
    double nan = 0.0 / 0.0;
    g->kind = kind;
    g->rs = nan;
    g->gm = nan;
    g->angmom = nan;
    g->impact = nan;
    g->c0 = nan;
    g->u = nan;
    g->up = nan;
    g->phi = nan;
}

int stc_geodesic_photon_init(StcGeodesic *g, double rs, double b,
                             double r_start, double phi_start) {
    if (g == NULL) {
        return -1;
    }
    stc_geo_reset(g, STC_GEO_PHOTON);
    if (stc_geo_bad(rs) || stc_geo_bad(b) || stc_geo_bad(r_start) ||
        stc_geo_bad(phi_start)) {
        return -1;
    }
    if (rs <= 0.0 || b <= 0.0 || r_start <= rs || r_start <= b) {
        return -1;
    }
    /* Null energy equation: up^2 = 1/b^2 - u^2 + rs*u^3. */
    {
        double u0 = 1.0 / r_start;
        double inv_b = 1.0 / b;
        double disc = inv_b * inv_b - u0 * u0 + rs * u0 * u0 * u0;
        double c = STC_C;
        if (!(disc >= 0.0) || disc == 1.0 / 0.0) {
            return -1;
        }
        g->rs = rs;
        g->gm = 0.5 * rs * c * c; /* rs = 2GM/c^2 */
        g->angmom = 0.0 / 0.0;
        g->impact = b;
        g->c0 = 0.0;
        g->u = u0;
        g->up = sqrt(disc); /* inbound: u grows with phi */
        g->phi = phi_start;
    }
    return 0;
}

int stc_geodesic_orbit_init(StcGeodesic *g, double gm, double rs,
                            double rp, double ecc, double phi_start) {
    if (g == NULL) {
        return -1;
    }
    stc_geo_reset(g, STC_GEO_TIMELIKE);
    if (stc_geo_bad(gm) || stc_geo_bad(rs) || stc_geo_bad(rp) ||
        stc_geo_bad(ecc) || stc_geo_bad(phi_start)) {
        return -1;
    }
    if (gm <= 0.0 || rs < 0.0 || rp <= 0.0 || rp <= rs || ecc < 0.0) {
        return -1;
    }
    {
        double p = rp * (1.0 + ecc); /* semi-latus rectum, m */
        double L;
        if (!(p > 0.0) || p == 1.0 / 0.0) {
            return -1;
        }
        L = sqrt(gm * p);
        if (stc_geo_bad(L) || L <= 0.0) {
            return -1;
        }
        g->rs = rs;
        g->gm = gm;
        g->angmom = L;
        g->impact = 0.0 / 0.0;
        g->c0 = gm / (L * L); /* = 1/p */
        g->u = 1.0 / rp;      /* periapsis passage */
        g->up = 0.0;
        g->phi = phi_start;
    }
    return 0;
}

double stc_geodesic_rhs(const StcGeodesic *g, double u) {
    if (g == NULL || stc_geo_bad(u)) {
        return 0.0 / 0.0;
    }
    if (stc_geo_bad(g->c0) || stc_geo_bad(g->rs)) {
        return 0.0 / 0.0;
    }
    /* Relativistic Binet equation: f(u) = c0 - u + 1.5*rs*u^2. */
    return g->c0 - u + 1.5 * g->rs * u * u;
}

void stc_geodesic_step(StcGeodesic *g, double dphi) {
    if (g == NULL || stc_geo_bad(dphi)) {
        return;
    }
    if (stc_geo_bad(g->u) || stc_geo_bad(g->up) || stc_geo_bad(g->phi)) {
        return;
    }
    {
        double h = dphi;
        double u = g->u;
        double w = g->up;
        double c0 = g->c0;
        double rs15 = 1.5 * g->rs;
        double k1u, k1w, k2u, k2w, k3u, k3w, k4u, k4w, uu;
        if (stc_geo_bad(c0) || stc_geo_bad(rs15)) {
            return;
        }
        /* Classical RK4 on (u, w) with w' = c0 - u + 1.5*rs*u^2. */
        k1u = w;
        k1w = c0 - u + rs15 * u * u;
        uu = u + 0.5 * h * k1u;
        k2u = w + 0.5 * h * k1w;
        k2w = c0 - uu + rs15 * uu * uu;
        uu = u + 0.5 * h * k2u;
        k3u = w + 0.5 * h * k2w;
        k3w = c0 - uu + rs15 * uu * uu;
        uu = u + h * k3u;
        k4u = w + h * k3w;
        k4w = c0 - uu + rs15 * uu * uu;
        g->u = u + (h / 6.0) * (k1u + 2.0 * k2u + 2.0 * k3u + k4u);
        g->up = w + (h / 6.0) * (k1w + 2.0 * k2w + 2.0 * k3w + k4w);
        g->phi = g->phi + h;
    }
}

double stc_geodesic_radius(const StcGeodesic *g) {
    if (g == NULL || stc_geo_bad(g->u)) {
        return 0.0 / 0.0;
    }
    if (g->u <= 0.0) {
        return 1.0 / 0.0;
    }
    return 1.0 / g->u;
}

StcGeoStatus stc_geodesic_status(const StcGeodesic *g, double r_escape) {
    double r;
    if (g == NULL) {
        return STC_GEO_ORBITING;
    }
    if (stc_geo_bad(g->u)) {
        /* Non-finite u only arises from a horizon plunge overflow. */
        return STC_GEO_CAPTURED;
    }
    if (g->u <= 0.0) {
        return STC_GEO_ESCAPED;
    }
    r = 1.0 / g->u;
    if (!stc_geo_bad(g->rs) && g->rs >= 0.0 && r <= g->rs) {
        return STC_GEO_CAPTURED;
    }
    if (!stc_geo_bad(r_escape) && r_escape > 0.0 && r >= r_escape) {
        return STC_GEO_ESCAPED;
    }
    return STC_GEO_ORBITING;
}

StcGeoStatus stc_geodesic_integrate(StcGeodesic *g, double dphi, int steps,
                                    double r_escape,
                                    double *r_hist, double *phi_hist) {
    int i;
    if (g == NULL) {
        return STC_GEO_ORBITING;
    }
    if (steps <= 0 || stc_geo_bad(dphi) || dphi <= 0.0) {
        return stc_geodesic_status(g, r_escape);
    }
    for (i = 0; i < steps; i++) {
        StcGeoStatus s;
        stc_geodesic_step(g, dphi);
        if (r_hist != NULL) {
            r_hist[i] = stc_geodesic_radius(g);
        }
        if (phi_hist != NULL) {
            phi_hist[i] = g->phi;
        }
        s = stc_geodesic_status(g, r_escape);
        if (s != STC_GEO_ORBITING) {
            return s;
        }
    }
    return stc_geodesic_status(g, r_escape);
}

double stc_geodesic_energy(const StcGeodesic *g) {
    if (g == NULL || g->kind != STC_GEO_TIMELIKE) {
        return 0.0 / 0.0;
    }
    if (stc_geo_bad(g->u) || stc_geo_bad(g->up) || stc_geo_bad(g->rs) ||
        stc_geo_bad(g->angmom) || g->angmom <= 0.0) {
        return 0.0 / 0.0;
    }
    {
        double c = STC_C;
        double L = g->angmom;
        double c_over_L = c / L;
        double rad = g->up * g->up +
                     (1.0 - g->rs * g->u) * (c_over_L * c_over_L + g->u * g->u);
        if (!(rad >= 0.0) || rad == 1.0 / 0.0) {
            return 0.0 / 0.0;
        }
        return c * L * sqrt(rad);
    }
}

double stc_geodesic_angmom(const StcGeodesic *g) {
    if (g == NULL || g->kind != STC_GEO_TIMELIKE) {
        return 0.0 / 0.0;
    }
    return g->angmom;
}

double stc_geodesic_impact(const StcGeodesic *g) {
    if (g == NULL || g->kind != STC_GEO_PHOTON) {
        return 0.0 / 0.0;
    }
    return g->impact;
}

double stc_geodesic_photon_deflection(double rs, double b, double r_start,
                                      double dphi, int max_steps) {
    StcGeodesic g;
    double phi_start = 0.0;
    double swept;
    StcGeoStatus s;
    if (stc_geo_bad(rs) || stc_geo_bad(b) || stc_geo_bad(r_start) ||
        stc_geo_bad(dphi)) {
        return 0.0 / 0.0;
    }
    if (rs <= 0.0 || b <= 0.0 || r_start <= rs || r_start <= b ||
        dphi <= 0.0 || max_steps <= 0) {
        return 0.0 / 0.0;
    }
    if (stc_geodesic_photon_init(&g, rs, b, r_start, phi_start) != 0) {
        return 0.0 / 0.0;
    }
    /* Own loop (not integrate) so the escape crossing u = 1/r_start
     * can be linearly interpolated: stopping on the coarse post-step
     * sample would quantize phi to dphi, dwarfing the deflection. */
    {
        double u_esc = 1.0 / r_start;
        double u_prev = g.u;
        double phi_prev = g.phi;
        int i;
        for (i = 0; i < max_steps; i++) {
            double u_now, phi_now, denom, frac;
            stc_geodesic_step(&g, dphi);
            s = stc_geodesic_status(&g, r_start);
            if (s == STC_GEO_CAPTURED) {
                return 0.0 / 0.0;
            }
            if (s != STC_GEO_ESCAPED) {
                u_prev = g.u;
                phi_prev = g.phi;
                continue;
            }
            u_now = g.u;
            phi_now = g.phi;
            denom = u_prev - u_now;
            frac = (denom != 0.0) ? (u_prev - u_esc) / denom : 1.0;
            if (!(frac >= 0.0)) {
                frac = 0.0;
            }
            if (frac > 1.0) {
                frac = 1.0;
            }
            swept = (phi_prev + frac * (phi_now - phi_prev)) - phi_start;
            /* Undeflected straight-line sweep from r_start to r_start. */
            return swept - 2.0 * acos(b / r_start);
        }
    }
    return 0.0 / 0.0;
}

double stc_geodesic_cone_tilt(double r, double rs) {
    double tilt;
    if (stc_geo_bad(r) || stc_geo_bad(rs)) {
        return 0.0 / 0.0;
    }
    if (r <= 0.0 || rs < 0.0) {
        return 0.0 / 0.0;
    }
    if (r <= rs) {
        return 0.25 * STC_PI;
    }
    tilt = 0.25 * STC_PI * (rs / r);
    if (tilt < 0.0) {
        tilt = 0.0;
    }
    if (tilt > 0.25 * STC_PI) {
        tilt = 0.25 * STC_PI;
    }
    return tilt;
}
