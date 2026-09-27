#ifndef STC_TIDES_H
#define STC_TIDES_H

/* SPACETIME LAB — tides: tidal acceleration, Roche limits, tidal heating.
 * Pure C11, SI doubles, no state, no sim mutation, no I/O.
 * Every function below names its equation.
 */

/* Equation: tidal acceleration (differential Newtonian gravity, first-order
 * expansion in R/r) —
 * a_tidal = 2 G M R / r^3, directed along +r_hat (center-to-center).
 * M_central: gravitating mass, kg (must be finite and >= 0; 0 yields the
 * zero vector); r_vec: separation vector from the central mass to the body
 * center, m (3-vector, must be finite and nonzero); R_body: radius of the
 * tidally distorted body, m (must be finite and >= 0; 0 yields the zero
 * vector); out: characteristic tidal acceleration vector, m/s^2
 * (3-vector): the near-side-minus-center gravity difference to leading
 * order. NULL out is a no-op; NULL r_vec, r == 0, or invalid /
 * non-finite scalars fill out with NaN.
 * Approximation: valid for R << r (terms of order (R/r)^2 dropped); the
 * body is a test mass (M_central >> m_body).
 * Cite: tidal acceleration (Newtonian differential gravity). */
void stc_tidal_accel(double M_central, const double r_vec[3], double R_body,
                     double out[3]);

/* Equation: Roche limit, fluid satellite (strengthless body in synchronous
 * rotation) —
 * d = 2.44 R_primary (rho_primary / rho_satellite)^(1/3).
 * R_primary: primary radius, m (must be finite and > 0); rho_primary:
 * primary bulk density, kg/m^3 (must be finite and > 0); rho_satellite:
 * satellite bulk density, kg/m^3 (must be finite and > 0).
 * Returns the disruption distance in m. NaN if any input is invalid.
 * This is the FLUID variant (larger distance): a fluid satellite held
 * together by self-gravity alone disrupts inside d.
 * Cite: Roche limit (E. Roche, 1848), fluid case. */
double stc_roche_limit_fluid(double R_primary, double rho_primary,
                             double rho_satellite);

/* Equation: Roche limit, rigid satellite (strengthless rigid body) —
 * d = 1.26 R_primary (rho_primary / rho_satellite)^(1/3).
 * Parameters, units, and NaN rule are identical to the fluid variant.
 * This is the RIGID variant (smaller distance): a rigid satellite held
 * together by self-gravity alone disrupts inside d.
 * Cite: Roche limit (E. Roche, 1848), rigid case. */
double stc_roche_limit_rigid(double R_primary, double rho_primary,
                             double rho_satellite);

/* Equation: Roche limit, default (fluid) —
 * identical to stc_roche_limit_fluid; the fluid value is the conventional
 * quoted Roche limit. Parameters, units, and NaN rule are identical. */
double stc_roche_limit(double R_primary, double rho_primary,
                       double rho_satellite);

/* Equation: eccentricity tidal heating (synchronous satellite, leading
 * order in e) —
 * dE/dt = (21/2) (k2/Q) G M^2 R^5 n e^2 / a^6,
 * with mean motion n = sqrt(G M / a^3).
 * M_central: gravitating mass, kg (must be finite and > 0); R_body:
 * heated-body radius, m (must be finite and > 0); a: orbital
 * semi-major axis, m (must be finite and > 0); e: orbital eccentricity,
 * dimensionless (must satisfy 0 <= e < 1; 0 yields 0, circular orbits
 * have no eccentricity heating); k2_over_Q: tidal dissipation factor,
 * dimensionless Love number over quality factor (must be finite and
 * >= 0; 0 yields 0; rocky satellites ~1e-3..1e-2, gas giants ~1e-5..1e-4).
 * Returns the dissipated power in W (J/s). NaN if any input is invalid.
 * Approximation: ORDER OF MAGNITUDE only. Assumes M_central >> m_body
 * (mean motion from the central mass alone), synchronous rotation, zero
 * obliquity, small e (leading e^2 term only), and frequency-independent
 * k2/Q. Factor-level differences between published forms are expected.
 * Cite: eccentricity tidal dissipation (Peale-Cassen-Segatz mechanism;
 * Io/Jupiter is the type example). */
double stc_tidal_heating(double M_central, double R_body, double a, double e,
                         double k2_over_Q);

#endif /* STC_TIDES_H */
