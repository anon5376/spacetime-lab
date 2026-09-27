#ifndef STC_KERR_H
#define STC_KERR_H

/* SPACETIME LAB — Kerr (rotating black hole) quantities.
 *
 * Units: geometric meters for every length (M = G*M_kg/c^2, spin a,
 * radii). Use stc_geometric_M() from stc/units.h to convert kg -> m.
 * Frame-drag rates are SI rad/s (geometric 1/m times STC_C).
 * Spin sign: a > 0 prograde, a < 0 retrograde, |a| <= M enforced by
 * clamping (see stc_kerr_clamp_spin).
 *
 * Pure C11, doubles only, no state, deterministic (no iteration).
 */

 /* Spin clamp: Kerr bound |a| <= M.
  * Returns a clamped to [-M, +M]. M <= 0 or non-finite input -> 0.0. */
double stc_kerr_clamp_spin(double m_m, double a_m);

/* Outer horizon: r+ = M + sqrt(M^2 - a^2) (Boyer-Lindquist).
 * Spin clamped to |a| <= M (superextremal -> extremal).
 * M <= 0 or non-finite input -> NaN. */
double stc_kerr_horizon(double m_m, double a_m);

/* Inner (Cauchy) horizon: r- = M - sqrt(M^2 - a^2).
 * Same clamping/NaN rules as stc_kerr_horizon. */
double stc_kerr_horizon_inner(double m_m, double a_m);

/* Outer ergosurface: r_e(theta) = M + sqrt(M^2 - a^2*cos^2(theta)).
 * theta_rad is the polar angle from the spin axis (0 = pole, pi/2 =
 * equator). Spin clamped to |a| <= M. M <= 0 or non-finite input
 * (including theta) -> NaN. */
double stc_kerr_ergosphere(double m_m, double a_m, double theta_rad);

/* ISCO radius: Bardeen-Press-Teukolsky (1972) formula with signed spin.
 * a >= 0 takes the prograde branch (6M -> 1M), a < 0 the retrograde
 * branch (6M -> 9M). Spin clamped to |a| <= M.
 * M <= 0 or non-finite input -> NaN. */
double stc_kerr_isco(double m_m, double a_m);

/* ISCO radius: Bardeen prograde branch (a >= 0 corotating).
 * a_abs_m is the spin magnitude; |a_abs_m| is used, clamped to M.
 * Limits: a = 0 -> 6M, a = +M -> 1M.
 * M <= 0 or non-finite input -> NaN. */
double stc_kerr_isco_prograde(double m_m, double a_abs_m);

/* ISCO radius: Bardeen retrograde branch (a <= 0 counterrotating).
 * a_abs_m is the spin magnitude; |a_abs_m| is used, clamped to M.
 * Limits: a = 0 -> 6M, a = -M -> 9M.
 * M <= 0 or non-finite input -> NaN. */
double stc_kerr_isco_retrograde(double m_m, double a_abs_m);

/* ZAMO frame-dragging rate on the equator:
 *   w = -g_tφ/g_φφ = 2*M*a*r / ((r^2+a^2)^2 - a^2*Δ),
 *   Δ = r^2 - 2*M*r + a^2,
 * converted to SI rad/s (times STC_C). a = 0 -> 0 exactly.
 * Spin clamped to |a| <= M. M <= 0, r <= 0, or non-finite input -> NaN.
 * Negative for a < 0 (dragging in the -φ sense). */
double stc_kerr_frame_drag(double m_m, double a_m, double r_m);

/* ZAMO frame-dragging rate at polar angle theta:
 *   w = 2*M*a*r / ((r^2+a^2)^2 - a^2*Δ*sin^2(theta)),
 * in SI rad/s. theta_rad = pi/2 matches stc_kerr_frame_drag.
 * Same clamping/NaN rules, plus non-finite theta -> NaN. */
double stc_kerr_frame_drag_theta(double m_m, double a_m, double r_m,
                                 double theta_rad);

/* Shadow (photon-ring) radius estimator, geometric meters.
 * APPROXIMATION (documented): r_sh = 3*sqrt(3)*M*(1 - 0.06*(a/M)^2),
 * a parabolic fit to the inclination-averaged mean shadow radius. Exact
 * 3*sqrt(3)*M at a = 0; ~4.85M at |a| = M, within a few percent of
 * ray-traced Kerr shadows. Ignores the prograde-side asymmetry (dent)
 * and inclination dependence; NOT a ray-traced shadow.
 * Spin clamped to |a| <= M. M <= 0 or non-finite input -> NaN. */
double stc_kerr_shadow_radius(double m_m, double a_m);

#endif /* STC_KERR_H */
