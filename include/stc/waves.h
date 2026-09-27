#ifndef STC_WAVES_H
#define STC_WAVES_H

/* SPACETIME LAB — gravitational waves: quadrupole strain, TT-gauge ring
 * response, Newtonian binary chirp.
 * Pure C11, SI doubles (kg, m, s, Hz, rad; strain dimensionless).
 *
 * APPROXIMATION (applies to every function below): Newtonian quadrupole,
 * leading-order (0PN), circular orbit of point masses, adiabatic inspiral.
 * Valid for v << c and separation a >> rs; breaks down near merger (cut
 * off at ISCO in the scenario layer). No PN corrections, spin,
 * eccentricity, or merger/ringdown. f(tau) and Phi(tau) diverge as
 * tau -> 0, which is the approximation failing, not a physical prediction.
 */

#include "stc/vec.h"

#include <stdint.h>

/* Chirp mass: Mc = (m1*m2)^(3/5) / (m1+m2)^(1/5).
 * Inputs in kg, returns kg. Non-finite or non-positive input -> NaN. */
double stc_waves_chirp_mass(double m1_kg, double m2_kg);

/* Circular-binary quadrupole strain amplitude (Newtonian, leading-order):
 *   h0 = (4/d) * (G*Mc/c^2)^(5/3) * (pi*f/c)^(2/3).
 * Mc in kg, f in Hz, d in m; returns dimensionless strain.
 * Non-finite or non-positive input -> NaN. */
double stc_waves_strain_amplitude(double chirp_mass_kg, double freq_hz,
                                  double dist_m);

/* Newtonian quadrupole strain, separation form:
 *   h = 4*G^2*m1*m2 / (c^4*a*R).
 * Masses in kg, orbital separation a and distance R in m; dimensionless.
 * Algebraically identical to stc_waves_strain_amplitude via Kepler's law
 * f = (1/pi)*sqrt(G*(m1+m2)/a^3). Non-finite/non-positive input -> NaN. */
double stc_waves_quadrupole_strain(double m1_kg, double m2_kg,
                                   double separation_m, double dist_m);

/* Newtonian chirp frequency evolution (quadrupole energy loss):
 *   f(tau) = (1/pi) * (5/(256*tau))^(3/8) * (G*Mc/c^3)^(-5/8),
 * where tau = t_c - t is the time to coalescence.
 * Mc in kg, tau in s; returns Hz. Non-finite/non-positive input -> NaN. */
double stc_waves_chirp_freq(double chirp_mass_kg, double t_to_coal_s);

/* Time to coalescence (inverse chirp relation):
 *   tau(f) = (5/256) * (G*Mc/c^3)^(-5/3) * (pi*f)^(-8/3).
 * Mc in kg, f in Hz; returns s. Non-finite/non-positive input -> NaN. */
double stc_waves_time_to_coalescence(double chirp_mass_kg, double freq_hz);

/* Newtonian inspiral GW phase:
 *   Phi(tau) = Phi_c - 2 * (tau / (5*G*Mc/c^3))^(5/8).
 * Mc in kg, tau in s, coalescence phase Phi_c in rad; returns rad.
 * Non-finite Mc/tau/Phi_c or non-positive Mc/tau -> NaN. */
double stc_waves_phase(double chirp_mass_kg, double t_to_coal_s,
                       double phase_c_rad);

/* Deterministic coalescence phase from seed: splitmix64 hash of the seed
 * mapped to [0, 2*pi). Same seed = same phase on any platform. */
double stc_waves_phase0_from_seed(uint64_t seed);

/* TT-gauge geodesic-deviation ring response (plus/cross), wave
 * propagating along +z acting on a particle at p0 (meters):
 *   x = x0 + (h+*x0 + hx*y0)/2
 *   y = y0 + (hx*x0 - h+*y0)/2
 *   z = z0 (transverse: no response along propagation).
 * NULL out -> no-op. Non-finite h+/hx -> identity copy of p0. */
void stc_waves_ring_response(double h_plus, double h_cross, StcVec3 p0,
                             StcVec3 *out);

/* Newtonian chirp waveform (restricted quadrupole, inclination-projected):
 *   f  = chirp_freq(Mc, tau), Phi = phase(Mc, tau, Phi_c),
 *   h0 = strain_amplitude(Mc, f, d),
 *   h+ = h0*(1+cos^2 i)/2 * cos(Phi),  hx = h0*cos i * sin(Phi).
 * Mc in kg, d in m, tau in s, Phi_c and inclination i in rad.
 * NULL h_plus/h_cross tolerated individually; invalid physics input
 * writes NaN into each non-NULL output. */
void stc_waves_chirp_waveform(double chirp_mass_kg, double dist_m,
                              double t_to_coal_s, double phase_c_rad,
                              double inclination_rad, double *h_plus,
                              double *h_cross);

#endif /* STC_WAVES_H */
