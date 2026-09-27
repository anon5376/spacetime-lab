/* SPACETIME LAB — lensing: weak-field deflection + point-mass thin lens.
 * SI doubles; closed-form only (no iteration, no RNG). See lensing.h.
 */

#include "stc/lensing.h"
#include "stc/units.h"

#include <math.h>
#include <stddef.h>

double stc_lensing_deflection(double m_kg, double b_m)
{
    if (!isfinite(m_kg) || !isfinite(b_m) || m_kg < 0.0 || b_m <= 0.0) {
        return NAN;
    }
    if (m_kg == 0.0) {
        return 0.0;
    }
    /* alpha = 4GM / (c^2 b) */
    double c = STC_C;
    return 4.0 * STC_G * m_kg / (c * c * b_m);
}

double stc_lensing_einstein_angle(double m_kg, double dl_m,
                                  double ds_m, double dls_m)
{
    if (!isfinite(m_kg) || !isfinite(dl_m) || !isfinite(ds_m) ||
        !isfinite(dls_m) || m_kg < 0.0 ||
        dl_m <= 0.0 || ds_m <= 0.0 || dls_m <= 0.0) {
        return NAN;
    }
    if (m_kg == 0.0) {
        return 0.0;
    }
    /* thetaE^2 = (4GM/c^2) * Dls / (Dl * Ds) */
    double c = STC_C;
    double rs = 2.0 * STC_G * m_kg / (c * c); /* Schwarzschild radius */
    double t2 = 2.0 * rs * dls_m / (dl_m * ds_m);
    return sqrt(t2);
}

double stc_lensing_einstein_radius(double m_kg, double dl_m,
                                   double ds_m, double dls_m)
{
    double th = stc_lensing_einstein_angle(m_kg, dl_m, ds_m, dls_m);
    if (!isfinite(th)) {
        return NAN;
    }
    return th * dl_m;
}

double stc_lensing_beta_of_theta(double theta_rad, double thetaE_rad)
{
    if (!isfinite(theta_rad) || !isfinite(thetaE_rad) ||
        theta_rad == 0.0 || thetaE_rad <= 0.0) {
        return NAN;
    }
    /* beta = theta - thetaE^2 / theta */
    double t2 = thetaE_rad * thetaE_rad;
    return theta_rad - t2 / theta_rad;
}

StcLensImages stc_lensing_images(double beta_rad, double thetaE_rad)
{
    StcLensImages out;
    if (!isfinite(beta_rad) || !isfinite(thetaE_rad) || thetaE_rad <= 0.0) {
        out.theta_plus = NAN;
        out.theta_minus = NAN;
        out.n_images = 0;
        return out;
    }
    if (beta_rad == 0.0) {
        out.theta_plus = thetaE_rad;
        out.theta_minus = -thetaE_rad;
        out.n_images = 1;
        return out;
    }
    /* theta_pm = (beta +- sqrt(beta^2 + 4 thetaE^2)) / 2 */
    double disc = sqrt(beta_rad * beta_rad + 4.0 * thetaE_rad * thetaE_rad);
    out.theta_plus = 0.5 * (beta_rad + disc);
    out.theta_minus = 0.5 * (beta_rad - disc);
    out.n_images = 2;
    return out;
}

StcLensMagnification stc_lensing_magnification(double beta_rad,
                                               double thetaE_rad)
{
    StcLensMagnification out;
    if (!isfinite(beta_rad) || !isfinite(thetaE_rad) || thetaE_rad <= 0.0) {
        out.mu_plus = NAN;
        out.mu_minus = NAN;
        out.mu_total = NAN;
        return out;
    }
    if (beta_rad == 0.0) {
        out.mu_plus = INFINITY;
        out.mu_minus = -INFINITY;
        out.mu_total = INFINITY;
        return out;
    }
    /* u = beta/thetaE; mu_pm = 1/2 +- (u^2+2)/(2 u sqrt(u^2+4)) */
    double u = beta_rad / thetaE_rad;
    double au = fabs(u);
    double s = sqrt(au * au + 4.0);
    double k = (au * au + 2.0) / (2.0 * au * s);
    if (u > 0.0) {
        out.mu_plus = 0.5 + k;
        out.mu_minus = 0.5 - k;
    } else {
        /* Signed-u form: mu_plus = 1/2 - k, mu_minus = 1/2 + k.
         * For beta < 0 the minus-branch image sits outside the
         * Einstein ring (positive parity), the plus-branch inside. */
        out.mu_plus = 0.5 - k;
        out.mu_minus = 0.5 + k;
    }
    out.mu_total = fabs(out.mu_plus) + fabs(out.mu_minus);
    return out;
}

double stc_lensing_total_magnification(double beta_rad, double thetaE_rad)
{
    return stc_lensing_magnification(beta_rad, thetaE_rad).mu_total;
}

int stc_lensing_ring_points(double thetaE_rad, double dl_m,
                            int n, StcVec3 *out)
{
    if (n <= 0 || out == NULL || !isfinite(thetaE_rad) ||
        !isfinite(dl_m) || thetaE_rad <= 0.0 || dl_m <= 0.0) {
        return 0;
    }
    double r = thetaE_rad * dl_m;
    if (!isfinite(r)) {
        return 0;
    }
    for (int i = 0; i < n; i++) {
        double phi = 2.0 * STC_PI * (double)i / (double)n;
        out[i].x = r * cos(phi);
        out[i].y = r * sin(phi);
        out[i].z = dl_m;
    }
    return n;
}
