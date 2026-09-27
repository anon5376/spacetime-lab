/* SPACETIME LAB — Kerr (rotating black hole) quantities.
 * Geometric meters; frame-drag rates in SI rad/s. See stc/kerr.h.
 * Pure C11, libm only, deterministic (closed forms, no iteration).
 */

#include "stc/kerr.h"
#include "stc/units.h"

#include <math.h>

/* sqrt(27) = 3*sqrt(3): Schwarzschild shadow radius in units of M. */
#define STC_KERR_SHADOW_PARABOLIC_COEF 0.06

double stc_kerr_clamp_spin(double m_m, double a_m)
{
    if (!isfinite(m_m) || !isfinite(a_m) || m_m <= 0.0) {
        return 0.0;
    }
    if (a_m > m_m) {
        return m_m;
    }
    if (a_m < -m_m) {
        return -m_m;
    }
    return a_m;
}

double stc_kerr_horizon(double m_m, double a_m)
{
    if (!isfinite(m_m) || !isfinite(a_m) || m_m <= 0.0) {
        return NAN;
    }
    /* Outer horizon: r+ = M + sqrt(M^2 - a^2). */
    double a = stc_kerr_clamp_spin(m_m, a_m);
    double disc = m_m * m_m - a * a;
    if (disc < 0.0) {
        disc = 0.0;
    }
    return m_m + sqrt(disc);
}

double stc_kerr_horizon_inner(double m_m, double a_m)
{
    if (!isfinite(m_m) || !isfinite(a_m) || m_m <= 0.0) {
        return NAN;
    }
    /* Inner (Cauchy) horizon: r- = M - sqrt(M^2 - a^2). */
    double a = stc_kerr_clamp_spin(m_m, a_m);
    double disc = m_m * m_m - a * a;
    if (disc < 0.0) {
        disc = 0.0;
    }
    return m_m - sqrt(disc);
}

double stc_kerr_ergosphere(double m_m, double a_m, double theta_rad)
{
    if (!isfinite(m_m) || !isfinite(a_m) || !isfinite(theta_rad) ||
        m_m <= 0.0) {
        return NAN;
    }
    /* Outer ergosurface: r_e = M + sqrt(M^2 - a^2*cos^2(theta)). */
    double a = stc_kerr_clamp_spin(m_m, a_m);
    double c = cos(theta_rad);
    double disc = m_m * m_m - a * a * c * c;
    if (disc < 0.0) {
        disc = 0.0;
    }
    return m_m + sqrt(disc);
}

/* Bardeen-Press-Teukolsky (1972) ISCO core. s = a/M in [-1, 1];
 * prograde selects the minus branch, retrograde the plus branch. */
static double stc_kerr_isco_core(double m_m, double s, int prograde)
{
    double t_plus = cbrt(1.0 + s);
    double t_minus = cbrt(1.0 - s);
    double t_prod = cbrt(1.0 - s * s);
    double z1 = 1.0 + t_prod * (t_plus + t_minus);
    double z2 = sqrt(3.0 * s * s + z1 * z1);
    double inner = (3.0 - z1) * (3.0 + z1 + 2.0 * z2);
    if (inner < 0.0) {
        inner = 0.0;
    }
    double root = sqrt(inner);
    if (prograde) {
        return m_m * (3.0 + z2 - root);
    }
    return m_m * (3.0 + z2 + root);
}

double stc_kerr_isco(double m_m, double a_m)
{
    if (!isfinite(m_m) || !isfinite(a_m) || m_m <= 0.0) {
        return NAN;
    }
    /* Signed Bardeen ISCO: prograde branch for a >= 0. */
    double a = stc_kerr_clamp_spin(m_m, a_m);
    double s = a / m_m;
    return stc_kerr_isco_core(m_m, s, s >= 0.0);
}

double stc_kerr_isco_prograde(double m_m, double a_abs_m)
{
    if (!isfinite(m_m) || !isfinite(a_abs_m) || m_m <= 0.0) {
        return NAN;
    }
    /* Bardeen prograde branch with the spin magnitude. */
    double a = stc_kerr_clamp_spin(m_m, fabs(a_abs_m));
    return stc_kerr_isco_core(m_m, a / m_m, 1);
}

double stc_kerr_isco_retrograde(double m_m, double a_abs_m)
{
    if (!isfinite(m_m) || !isfinite(a_abs_m) || m_m <= 0.0) {
        return NAN;
    }
    /* Bardeen retrograde branch with the spin magnitude. */
    double a = stc_kerr_clamp_spin(m_m, fabs(a_abs_m));
    return stc_kerr_isco_core(m_m, -(a / m_m), 0);
}

double stc_kerr_frame_drag_theta(double m_m, double a_m, double r_m,
                                 double theta_rad)
{
    if (!isfinite(m_m) || !isfinite(a_m) || !isfinite(r_m) ||
        !isfinite(theta_rad) || m_m <= 0.0 || r_m <= 0.0) {
        return NAN;
    }
    double a = stc_kerr_clamp_spin(m_m, a_m);
    if (a == 0.0) {
        return 0.0;
    }
    /* ZAMO rate w = 2*M*a*r / ((r^2+a^2)^2 - a^2*Δ*sin^2(theta)). */
    double r2 = r_m * r_m;
    double a2 = a * a;
    double delta = r2 - 2.0 * m_m * r_m + a2;
    double sum = r2 + a2;
    double sn = sin(theta_rad);
    double denom = sum * sum - a2 * delta * sn * sn;
    if (!isfinite(denom) || denom == 0.0) {
        return NAN;
    }
    double w_geom = 2.0 * m_m * a * r_m / denom; /* 1/m */
    return STC_C * w_geom; /* SI rad/s */
}

double stc_kerr_frame_drag(double m_m, double a_m, double r_m)
{
    if (!isfinite(m_m) || !isfinite(a_m) || !isfinite(r_m) ||
        m_m <= 0.0 || r_m <= 0.0) {
        return NAN;
    }
    double a = stc_kerr_clamp_spin(m_m, a_m);
    if (a == 0.0) {
        return 0.0;
    }
    /* Equatorial ZAMO rate: sin^2(theta) = 1. */
    double r2 = r_m * r_m;
    double a2 = a * a;
    double delta = r2 - 2.0 * m_m * r_m + a2;
    double sum = r2 + a2;
    double denom = sum * sum - a2 * delta;
    if (!isfinite(denom) || denom == 0.0) {
        return NAN;
    }
    double w_geom = 2.0 * m_m * a * r_m / denom; /* 1/m */
    return STC_C * w_geom; /* SI rad/s */
}

double stc_kerr_shadow_radius(double m_m, double a_m)
{
    if (!isfinite(m_m) || !isfinite(a_m) || m_m <= 0.0) {
        return NAN;
    }
    /* Mean-shadow estimator: 3*sqrt(3)*M*(1 - 0.06*(a/M)^2). */
    double a = stc_kerr_clamp_spin(m_m, a_m);
    double s = a / m_m;
    return 3.0 * sqrt(3.0) * m_m *
           (1.0 - STC_KERR_SHADOW_PARABOLIC_COEF * s * s);
}
