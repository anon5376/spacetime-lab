#ifndef STC_LENSING_H
#define STC_LENSING_H

/* SPACETIME LAB — gravitational lensing: weak-field deflection + thin lens.
 * Pure C11, SI doubles (kg, m, s, rad). All formulas closed-form and
 * deterministic (no RNG, no iteration). Angles in radians.
 */

#include "stc/vec.h"

/* Image positions of a point-mass lens along the lens axis (radians).
 * n_images: 2 = two distinct images, 1 = aligned Einstein ring
 * (theta_plus = +thetaE, theta_minus = -thetaE diametral markers),
 * 0 = invalid input (both fields NaN).
 */
typedef struct {
    double theta_plus;
    double theta_minus;
    int n_images;
} StcLensImages;

/* Point-mass lens magnifications (dimensionless, signed per image).
 * mu_total = |mu_plus| + |mu_minus| >= 1. Perfect alignment
 * (beta = 0) gives mu_total = +inf (point source on caustic).
 */
typedef struct {
    double mu_plus;
    double mu_minus;
    double mu_total;
} StcLensMagnification;

/* Equation: Einstein weak-field deflection.
 * alpha = 4*G*M / (c^2 * b).
 * m_kg: lens mass in kg; b_m: impact parameter in meters.
 * Returns deflection angle in radians. m_kg = 0 returns 0.
 * Returns NaN if m_kg < 0, b_m <= 0, or any input non-finite.
 */
double stc_lensing_deflection(double m_kg, double b_m);

/* Equation: Einstein angle (point-mass thin lens).
 * thetaE = sqrt(4*G*M/c^2 * Dls / (Dl * Ds)).
 * m_kg: lens mass in kg; dl_m: observer-lens distance (m);
 * ds_m: observer-source distance (m); dls_m: lens-source distance (m).
 * Returns Einstein angle in radians. m_kg = 0 returns 0.
 * Returns NaN if m_kg < 0, any distance <= 0, or any input non-finite.
 */
double stc_lensing_einstein_angle(double m_kg, double dl_m,
                                  double ds_m, double dls_m);

/* Equation: Einstein radius at the lens plane.
 * RE = thetaE * Dl (meters), with thetaE the Einstein angle above.
 * Same parameters/units as stc_lensing_einstein_angle; returns meters.
 * Returns NaN on the same invalid inputs.
 */
double stc_lensing_einstein_radius(double m_kg, double dl_m,
                                   double ds_m, double dls_m);

/* Equation: point-mass thin-lens mapping (lens equation).
 * beta = theta - thetaE^2 / theta.
 * theta_rad: image angle (rad); thetaE_rad: Einstein angle (rad).
 * Returns source angle beta in radians.
 * Returns NaN if theta_rad = 0, thetaE_rad <= 0, or input non-finite.
 */
double stc_lensing_beta_of_theta(double theta_rad, double thetaE_rad);

/* Equation: point-mass image positions (quadratic lens-equation roots).
 * theta_pm = (beta +- sqrt(beta^2 + 4*thetaE^2)) / 2.
 * beta_rad: source angle (rad); thetaE_rad: Einstein angle (rad).
 * beta = 0 yields n_images = 1 (aligned ring; markers +-thetaE).
 * Invalid (thetaE <= 0 or non-finite input) yields n_images = 0, NaNs.
 */
StcLensImages stc_lensing_images(double beta_rad, double thetaE_rad);

/* Equation: point-mass magnifications.
 * mu_pm = 1/2 +- (u^2 + 2) / (2*u*sqrt(u^2 + 4)), u = beta / thetaE.
 * beta_rad: source angle (rad); thetaE_rad: Einstein angle (rad).
 * beta = 0 gives mu_plus = +inf, mu_minus = -inf, mu_total = +inf.
 * Invalid (thetaE <= 0 or non-finite input) gives NaN fields.
 */
StcLensMagnification stc_lensing_magnification(double beta_rad,
                                               double thetaE_rad);

/* Equation: point-mass total magnification (magnification invariant).
 * mu_total = |mu_plus| + |mu_minus| = (u^2 + 2) / (u*sqrt(u^2 + 4)).
 * Same inputs as stc_lensing_magnification; returns mu_total only.
 */
double stc_lensing_total_magnification(double beta_rad, double thetaE_rad);

/* Equation: aligned Einstein ring (beta = 0 circle).
 * Writes n evenly spaced points on the ring of angular radius thetaE_rad
 * at lens-plane distance dl_m: out[i] = (R*cos(phi), R*sin(phi), dl_m),
 * R = thetaE_rad * dl_m, phi = 2*pi*i/n, i = 0..n-1.
 * Returns n on success, 0 on error (n <= 0, out NULL, thetaE <= 0,
 * dl <= 0, or non-finite input; buffer untouched on error).
 */
int stc_lensing_ring_points(double thetaE_rad, double dl_m,
                            int n, StcVec3 *out);

#endif /* STC_LENSING_H */
