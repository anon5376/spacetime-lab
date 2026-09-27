/* SPACETIME LAB — elements: Hill/SOI radii, escape/circular velocity,
 * Lagrange points. See stc/elements.h for equations and units.
 * Pure C11, libm only, deterministic (closed forms, no iteration).
 */

#include "stc/elements.h"
#include "stc/units.h"

#include <math.h>
#include <stddef.h>

double stc_hill_radius(double a_m, double m_kg, double M_kg)
{
    if (!isfinite(a_m) || !isfinite(m_kg) || !isfinite(M_kg) ||
        a_m <= 0.0 || m_kg < 0.0 || M_kg <= 0.0) {
        return NAN;
    }
    if (m_kg == 0.0) {
        return 0.0;
    }
    /* Hill radius: r_H = a * (m / (3*M))^(1/3). Ratio first so a huge
     * M_kg cannot overflow 3*M before the division. */
    return a_m * cbrt(m_kg / M_kg / 3.0);
}

double stc_soi_radius(double a_m, double m_kg, double M_kg)
{
    if (!isfinite(a_m) || !isfinite(m_kg) || !isfinite(M_kg) ||
        a_m <= 0.0 || m_kg < 0.0 || M_kg <= 0.0) {
        return NAN;
    }
    if (m_kg == 0.0) {
        return 0.0;
    }
    /* Laplace sphere of influence: r_SOI = a * (m / M)^(2/5). */
    return a_m * pow(m_kg / M_kg, 0.4);
}

double stc_escape_velocity(double M_kg, double r_m)
{
    if (!isfinite(M_kg) || !isfinite(r_m) || M_kg < 0.0 || r_m <= 0.0) {
        return NAN;
    }
    if (M_kg == 0.0) {
        return 0.0;
    }
    /* Escape velocity: v_esc = sqrt(2*G*M / r). */
    return sqrt(2.0 * STC_G * M_kg / r_m);
}

double stc_circular_velocity(double mu_m3s2, double r_m)
{
    if (!isfinite(mu_m3s2) || !isfinite(r_m) || mu_m3s2 < 0.0 ||
        r_m <= 0.0) {
        return NAN;
    }
    if (mu_m3s2 == 0.0) {
        return 0.0;
    }
    /* Circular velocity: v_circ = sqrt(mu / r). */
    return sqrt(mu_m3s2 / r_m);
}

void stc_lagrange_points(double mu_ratio, double dist_m,
                         double L1[3], double L2[3], double L3[3],
                         double L4[3], double L5[3])
{
    double alpha, a2, a3, g1, g2, l3x, h;
    double *outs[5];
    int i;

    outs[0] = L1;
    outs[1] = L2;
    outs[2] = L3;
    outs[3] = L4;
    outs[4] = L5;
    if (!isfinite(mu_ratio) || !isfinite(dist_m) || mu_ratio <= 0.0 ||
        mu_ratio > 0.5 || dist_m <= 0.0) {
        /* Invalid scalars: NaN-fill every requested point. */
        for (i = 0; i < 5; i++) {
            if (outs[i] != NULL) {
                outs[i][0] = NAN;
                outs[i][1] = NAN;
                outs[i][2] = NAN;
            }
        }
        return;
    }

    /* Collinear series in alpha = (mu/3)^(1/3) through O(alpha^3),
     * L3 through O(mu); triangular points exact equilateral. */
    alpha = cbrt(mu_ratio / 3.0);
    a2 = alpha * alpha;
    a3 = a2 * alpha;
    g1 = alpha - a2 / 3.0 - a3 / 9.0;
    g2 = alpha + a2 / 3.0 - a3 / 9.0;
    l3x = -1.0 + 7.0 * mu_ratio / 12.0;
    h = dist_m * sqrt(3.0) / 2.0;

    if (L1 != NULL) {
        L1[0] = dist_m * (1.0 - g1);
        L1[1] = 0.0;
        L1[2] = 0.0;
    }
    if (L2 != NULL) {
        L2[0] = dist_m * (1.0 + g2);
        L2[1] = 0.0;
        L2[2] = 0.0;
    }
    if (L3 != NULL) {
        L3[0] = dist_m * l3x;
        L3[1] = 0.0;
        L3[2] = 0.0;
    }
    if (L4 != NULL) {
        L4[0] = dist_m * 0.5;
        L4[1] = h;
        L4[2] = 0.0;
    }
    if (L5 != NULL) {
        L5[0] = dist_m * 0.5;
        L5[1] = -h;
        L5[2] = 0.0;
    }
}
