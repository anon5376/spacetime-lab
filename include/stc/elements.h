#ifndef STC_ELEMENTS_H
#define STC_ELEMENTS_H

/* SPACETIME LAB — elements: Hill/SOI radii, escape/circular velocity,
 * Lagrange points (CR3BP).
 * Pure C11, SI doubles (kg, m, s, m/s), no state, no sim mutation, no I/O.
 * Closed forms + truncated series only (no iteration, no RNG);
 * deterministic. Every function below names its equation.
 */

/* Equation: Hill radius (Hill 1878) —
 * r_H = a * (m / (3*M))^(1/3).
 * a_m: orbiter semi-major axis about the central mass, m (must be finite
 * and > 0); m_kg: orbiter mass, kg (must be finite and >= 0; 0 yields 0);
 * M_kg: central mass, kg (must be finite and > 0).
 * Returns the Hill radius in m. NaN if any input is invalid.
 * Approximation: circular orbiter orbit, m << M (test-particle limit);
 * the true stability boundary is fuzzy (orbits beyond ~r_H/2..r_H/3 are
 * typically unstable). Cite: Hill sphere. */
double stc_hill_radius(double a_m, double m_kg, double M_kg);

/* Equation: Laplace sphere of influence (SOI) —
 * r_SOI = a * (m / M)^(2/5).
 * Parameters, units, and NaN rule are identical to stc_hill_radius:
 * a_m finite and > 0; m_kg finite and >= 0 (0 yields 0); M_kg finite
 * and > 0. Returns the SOI radius in m.
 * Approximation: circular orbiter orbit; boundary where the central
 * body's tidal perturbation equals the orbiter's own gravity (ratio of
 * perturbing to primary acceleration). Patched-conic heuristic, not a
 * sharp dynamical boundary. Cite: Laplace sphere of influence. */
double stc_soi_radius(double a_m, double m_kg, double M_kg);

/* Equation: escape velocity (parabolic orbit / energy conservation) —
 * v_esc = sqrt(2*G*M / r), G = STC_G.
 * M_kg: central mass, kg (must be finite and >= 0; 0 yields 0);
 * r_m: distance from the central mass center, m (must be finite
 * and > 0). Returns the escape speed in m/s. NaN if any input invalid.
 * Newtonian point mass from rest at infinity; ignores atmosphere,
 * rotation, and other bodies. */
double stc_escape_velocity(double M_kg, double r_m);

/* Equation: circular orbital velocity (Kepler / vis-viva at e = 0) —
 * v_circ = sqrt(mu / r).
 * mu_m3s2: gravitational parameter G*M, m^3/s^2 (must be finite and
 * >= 0; 0 yields 0); r_m: orbital radius, m (must be finite and > 0).
 * Returns the circular speed in m/s. NaN if any input is invalid.
 * Newtonian two-body circular orbit. */
double stc_circular_velocity(double mu_m3s2, double r_m);

/* Equation: Lagrange points of the circular restricted three-body
 * problem (Lagrange 1772).
 * mu_ratio: secondary mass fraction m2/(m1+m2) with m2 the SMALLER
 * mass, dimensionless (must be finite with 0 < mu_ratio <= 0.5);
 * dist_m: separation of the two primaries, m (must be finite
 * and > 0); L1..L5: outputs, one 3-vector each (m).
 * Frame: primary m1 at the origin (0,0,0), secondary m2 at
 * (dist_m,0,0); orbit plane is z = 0. L4 (+y) leads and L5 (-y)
 * trails for counterclockwise motion seen from +z.
 * Method: L1/L2/L3 (collinear) use a truncated power series in
 * alpha = (mu/3)^(1/3):
 *   gamma1 = alpha - alpha^2/3 - alpha^3/9,  L1 = ((1-gamma1)*d,0,0);
 *   gamma2 = alpha + alpha^2/3 - alpha^3/9,  L2 = ((1+gamma2)*d,0,0);
 *   L3 = ((-1 + 7*mu/12)*d,0,0);
 * L4/L5 (triangular) are the exact equilateral points
 * (d/2,+/-d*sqrt(3)/2,0).
 * Series order: L1/L2 through O(alpha^3) = O(mu) (cubic in alpha,
 * error O(alpha^4)); L3 through O(mu) (linear, error O(mu^2)).
 * Accuracy: collinear error < 0.1% for mu <~ 1e-3 (Sun-planet
 * scale), ~0.3% at the Earth-Moon mu = 0.012, a few percent at
 * mu = 0.1; exact only in the mu -> 0 limit. L4/L5 exact for all mu.
 * Each out pointer may be NULL (that point is skipped); all-NULL is
 * a no-op. Invalid / non-finite mu_ratio or dist_m fills every
 * non-NULL out with NaN.
 * Cite: Lagrange points (J.-L. Lagrange, 1772). */
void stc_lagrange_points(double mu_ratio, double dist_m,
                         double L1[3], double L2[3], double L3[3],
                         double L4[3], double L5[3]);

#endif /* STC_ELEMENTS_H */
