#include "stc/drag.h"

#include <math.h>
#include <stddef.h>

#include "stc/units.h"

double stc_atmo_density(double h, double rho0, double H)
{
    /* Barometric formula (isothermal exponential atmosphere):
     * rho(h) = rho0 * exp(-h / H). */
    if (!(rho0 >= 0.0) || !isfinite(rho0)) {
        return NAN;
    }
    if (!(H > 0.0) || !isfinite(H)) {
        return NAN;
    }
    return rho0 * exp(-h / H);
}

void stc_drag_accel(const double v_vec[3], double rho, double Cd, double A,
                    double m, double out[3])
{
    /* Drag equation: F = -(1/2) rho |v|^2 Cd A v_hat, so
     * a = F/m = -(1/2) rho |v| Cd A / m * v. */
    double speed, k;

    if (out == NULL) {
        return;
    }
    if (v_vec == NULL || !(m > 0.0) || !isfinite(m) ||
        !(rho >= 0.0) || !isfinite(rho) ||
        !(Cd >= 0.0) || !isfinite(Cd) ||
        !(A >= 0.0) || !isfinite(A) ||
        !isfinite(v_vec[0]) || !isfinite(v_vec[1]) || !isfinite(v_vec[2])) {
        out[0] = NAN;
        out[1] = NAN;
        out[2] = NAN;
        return;
    }
    speed = sqrt(v_vec[0] * v_vec[0] + v_vec[1] * v_vec[1] +
                 v_vec[2] * v_vec[2]);
    if (!(speed > 0.0) || !isfinite(speed)) {
        out[0] = 0.0;
        out[1] = 0.0;
        out[2] = 0.0;
        return;
    }
    k = -0.5 * rho * speed * Cd * A / m;
    out[0] = k * v_vec[0] + 0.0; /* +0.0 normalizes -0.0 */
    out[1] = k * v_vec[1] + 0.0;
    out[2] = k * v_vec[2] + 0.0;
}

void stc_radiation_pressure_accel(const double r_vec_sun[3],
                                  double luminosity, double area, double mass,
                                  double Q, double out[3])
{
    /* Radiation pressure: flux L/(4 pi r^2), force Q A flux/c radially
     * outward, so a = Q A L / (4 pi r^2 c m) * r_hat. */
    double r2, r, k;

    if (out == NULL) {
        return;
    }
    if (r_vec_sun == NULL || !(mass > 0.0) || !isfinite(mass) ||
        !(luminosity >= 0.0) || !isfinite(luminosity) ||
        !(area >= 0.0) || !isfinite(area) ||
        !(Q >= 0.0) || !isfinite(Q) ||
        !isfinite(r_vec_sun[0]) || !isfinite(r_vec_sun[1]) ||
        !isfinite(r_vec_sun[2])) {
        out[0] = NAN;
        out[1] = NAN;
        out[2] = NAN;
        return;
    }
    r2 = r_vec_sun[0] * r_vec_sun[0] + r_vec_sun[1] * r_vec_sun[1] +
         r_vec_sun[2] * r_vec_sun[2];
    r = sqrt(r2);
    if (!(r > 0.0) || !isfinite(r)) {
        out[0] = NAN;
        out[1] = NAN;
        out[2] = NAN;
        return;
    }
    k = Q * area * luminosity / (4.0 * STC_PI * r2 * STC_C * mass);
    out[0] = k * r_vec_sun[0] / r;
    out[1] = k * r_vec_sun[1] / r;
    out[2] = k * r_vec_sun[2] / r;
}
