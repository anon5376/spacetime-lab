#include "stc/tides.h"

#include <math.h>
#include <stddef.h>

#include "stc/units.h"

void stc_tidal_accel(double M_central, const double r_vec[3], double R_body,
                     double out[3])
{
    /* Tidal acceleration: near-side-minus-center Newtonian gravity,
     * first-order expansion in R/r: a = 2 G M R / r^3 along +r_hat. */
    double r2, r, mag;

    if (out == NULL) {
        return;
    }
    if (r_vec == NULL || !(M_central >= 0.0) || !isfinite(M_central) ||
        !(R_body >= 0.0) || !isfinite(R_body) ||
        !isfinite(r_vec[0]) || !isfinite(r_vec[1]) || !isfinite(r_vec[2])) {
        out[0] = NAN;
        out[1] = NAN;
        out[2] = NAN;
        return;
    }
    r2 = r_vec[0] * r_vec[0] + r_vec[1] * r_vec[1] + r_vec[2] * r_vec[2];
    r = sqrt(r2);
    if (!(r > 0.0) || !isfinite(r)) {
        out[0] = NAN;
        out[1] = NAN;
        out[2] = NAN;
        return;
    }
    mag = 2.0 * STC_G * M_central * R_body / (r2 * r);
    out[0] = mag * r_vec[0] / r;
    out[1] = mag * r_vec[1] / r;
    out[2] = mag * r_vec[2] / r;
}

double stc_roche_limit_fluid(double R_primary, double rho_primary,
                             double rho_satellite)
{
    /* Roche limit, fluid satellite: d = 2.44 R (rho_p / rho_s)^(1/3). */
    if (!(R_primary > 0.0) || !isfinite(R_primary) ||
        !(rho_primary > 0.0) || !isfinite(rho_primary) ||
        !(rho_satellite > 0.0) || !isfinite(rho_satellite)) {
        return NAN;
    }
    return 2.44 * R_primary * cbrt(rho_primary / rho_satellite);
}

double stc_roche_limit_rigid(double R_primary, double rho_primary,
                             double rho_satellite)
{
    /* Roche limit, rigid satellite: d = 1.26 R (rho_p / rho_s)^(1/3). */
    if (!(R_primary > 0.0) || !isfinite(R_primary) ||
        !(rho_primary > 0.0) || !isfinite(rho_primary) ||
        !(rho_satellite > 0.0) || !isfinite(rho_satellite)) {
        return NAN;
    }
    return 1.26 * R_primary * cbrt(rho_primary / rho_satellite);
}

double stc_roche_limit(double R_primary, double rho_primary,
                       double rho_satellite)
{
    /* Roche limit, default: the fluid value is the conventional quote. */
    return stc_roche_limit_fluid(R_primary, rho_primary, rho_satellite);
}

double stc_tidal_heating(double M_central, double R_body, double a, double e,
                         double k2_over_Q)
{
    /* Eccentricity tidal heating: dE/dt = (21/2)(k2/Q) G M^2 R^5 n e^2/a^6
     * with n = sqrt(G M/a^3); evaluated as (21/2)(k2/Q) G M^2 q^5 n e^2/a
     * with q = R/a for numerical stability. */
    double q, n, rate;

    if (!(M_central > 0.0) || !isfinite(M_central) ||
        !(R_body > 0.0) || !isfinite(R_body) ||
        !(a > 0.0) || !isfinite(a) ||
        !(e >= 0.0) || !(e < 1.0) || !isfinite(e) ||
        !(k2_over_Q >= 0.0) || !isfinite(k2_over_Q)) {
        return NAN;
    }
    if (e == 0.0 || k2_over_Q == 0.0) {
        return 0.0;
    }
    q = R_body / a;
    n = sqrt(STC_G * M_central / (a * a * a));
    rate = 10.5 * k2_over_Q * STC_G * M_central * M_central *
        q * q * q * q * q * n * e * e / a;
    return rate;
}
