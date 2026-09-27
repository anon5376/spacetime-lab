#ifndef STC_SR_H
#define STC_SR_H

/* SPACETIME LAB — special relativity: flat-spacetime kinematics.
 * Pure C11, SI doubles, deterministic (closed-form only, plus one
 * caller-bounded worldline sum; no RNG, no uncapped iteration).
 *
 * Units on every boundary: velocities m/s, times s, lengths m,
 * masses kg, energies J, momenta kg m/s, frequencies Hz, angles rad,
 * rapidity/beta/gamma/z dimensionless.
 *
 * Speed-limit guard: every velocity input is clamped internally to
 * |v| <= STC_SR_BETA_MAX * c (or mapped to a documented limit value),
 * and NaN inputs map to documented defaults, so these functions never
 * return NaN. Superluminal gamma/rapidity saturate at +/-HUGE_VAL.
 */

/* Internal |beta| clamp: 1 - 1e-12 (gamma ~ 707107, safely finite). */
#define STC_SR_BETA_MAX 0.999999999999

/* 1+1-D spacetime event: coordinate time in s, position in m. */
typedef struct {
    double t_s;
    double x_m;
} StcSrEvent;

/* Twin-paradox result: out-and-back at cruise speed, instant turnaround.
 * earth_time_s: coordinate time elapsed on Earth (s).
 * traveler_time_s: proper time on the traveler's clock (s).
 * gamma: Lorentz factor at cruise speed (dimensionless).
 * v_used_ms: cruise speed actually used after guarding (m/s). */
typedef struct {
    double earth_time_s;
    double traveler_time_s;
    double gamma;
    double v_used_ms;
} StcSrTwin;

/* Barn-pole (ladder) paradox result at speed |v|.
 * gamma: Lorentz factor (dimensionless).
 * pole_in_barn_m: pole length measured in the barn frame (m).
 * barn_in_pole_m: barn length measured in the pole frame (m).
 * simultaneity_gap_s: barn-door simultaneity offset in the pole frame (s).
 * fits_barn_frame: 1 if the contracted pole fits the barn (0/1).
 * fits_pole_frame: 1 if the pole fits the contracted barn (0/1).
 * The fits flags use a 1e-9 relative tolerance so textbook exact-fit
 * boundaries (e.g. 20 m pole at gamma = 2 vs 10 m barn) count as fits. */
typedef struct {
    double gamma;
    double pole_in_barn_m;
    double barn_in_pole_m;
    double simultaneity_gap_s;
    int fits_barn_frame;
    int fits_pole_frame;
} StcSrPole;

/* Lorentz factor: gamma = 1/sqrt(1-beta^2), beta = v/c.
 * v_ms in m/s. |v| >= c -> +HUGE_VAL; NaN -> 1.0. */
double stc_sr_gamma(double v_ms);

/* Lorentz factor from beta: gamma = 1/sqrt(1-beta^2).
 * |beta| >= 1 -> +HUGE_VAL; NaN -> 1.0. */
double stc_sr_gamma_from_beta(double beta);

/* Lorentz factor from rapidity: gamma = cosh(phi). NaN -> 1.0. */
double stc_sr_gamma_from_rapidity(double phi);

/* Beta factor: beta = v/c, clamped to [-1, 1].
 * v_ms in m/s. NaN -> 0.0; +/-inf -> +/-1.0. */
double stc_sr_beta(double v_ms);

/* Velocity from gamma: v = c*sqrt(1-1/gamma^2), in m/s.
 * gamma <= 1 or NaN -> 0.0; +inf -> c. */
double stc_sr_velocity_from_gamma(double gamma);

/* Lorentz transformation (x-boost of a 1+1-D event):
 * t' = gamma*(t - v*x/c^2), x' = gamma*(x - v*t).
 * e: event (s, m); v_ms: boost velocity in m/s (guarded to < c).
 * Non-finite event coordinates are treated as 0. Never NaN. */
StcSrEvent stc_sr_boost(StcSrEvent e, double v_ms);

/* Lorentz transformation (x-boost), pointer-output form.
 * Same equations as stc_sr_boost. NULL outputs are skipped. */
void stc_sr_transform(double t_s, double x_m, double v_ms,
                      double *out_t, double *out_x);

/* Lorentz transformation (x-boost), alias of stc_sr_boost. */
StcSrEvent stc_sr_lorentz(StcSrEvent e, double v_ms);

/* Spacetime interval: s^2 = (c*dt)^2 - dx^2, in m^2 (timelike positive).
 * dt_s in s, dx_m in m. NaN inputs -> 0; astronomic inputs saturate
 * at signed HUGE_VAL instead of overflowing to NaN. */
double stc_sr_interval(double dt_s, double dx_m);

/* Spacetime interval between two events: s^2 = (c*dt)^2 - dx^2, in m^2. */
double stc_sr_interval_events(StcSrEvent a, StcSrEvent b);

/* Einstein velocity addition (collinear): w = (u+v)/(1+u*v/c^2), in m/s.
 * Inputs guarded to < c; result is always subluminal (up to rounding). */
double stc_sr_add_velocity(double u_ms, double v_ms);

/* Einstein velocity addition (collinear), alias of stc_sr_add_velocity. */
double stc_sr_velocity_add(double u_ms, double v_ms);

/* Relativistic longitudinal Doppler factor:
 * f_obs/f_emit = sqrt((1-beta)/(1+beta)), v > 0 receding.
 * v_ms in m/s (guarded). Dimensionless, always finite and positive. */
double stc_sr_doppler_factor(double v_ms);

/* Relativistic longitudinal Doppler shift: f_obs = f0*factor(v), in Hz.
 * f0_hz <= 0 or NaN -> 0.0. */
double stc_sr_doppler(double f0_hz, double v_ms);

/* Relativistic transverse Doppler shift: f_obs = f0/gamma, in Hz.
 * f0_hz <= 0 or NaN -> 0.0. */
double stc_sr_doppler_transverse(double f0_hz, double v_ms);

/* Cosmological-style redshift from recession speed:
 * z = sqrt((1+beta)/(1-beta)) - 1, v > 0 receding.
 * v_ms in m/s (guarded). Dimensionless, always > -1. */
double stc_sr_redshift(double v_ms);

/* Time dilation: dt_lab = gamma*dt_proper.
 * dt_proper_s in s; returns lab-frame time in s. NaN -> 0.0. */
double stc_sr_time_dilation(double dt_proper_s, double v_ms);

/* Proper time from lab time: d_tau = dt_lab/gamma.
 * dt_lab_s in s; returns proper time in s. NaN -> 0.0. */
double stc_sr_proper_time(double dt_lab_s, double v_ms);

/* Length contraction: L = L0/gamma.
 * l_proper_m in m; returns contracted length in m.
 * L0 <= 0 or NaN -> 0.0. */
double stc_sr_length_contraction(double l_proper_m, double v_ms);

/* Rest energy (mass-energy equivalence): E0 = m*c^2, in J.
 * m_kg <= 0 or NaN -> 0.0. */
double stc_sr_rest_energy(double m_kg);

/* Total relativistic energy: E = gamma*m*c^2, in J.
 * m_kg <= 0 or NaN -> 0.0. */
double stc_sr_energy(double m_kg, double v_ms);

/* Relativistic kinetic energy: K = (gamma-1)*m*c^2, in J.
 * Uses the cancellation-free form (gamma-1) = beta^2*gamma^2/(1+gamma).
 * m_kg <= 0 or NaN -> 0.0. */
double stc_sr_kinetic_energy(double m_kg, double v_ms);

/* Relativistic momentum: p = gamma*m*v, in kg m/s.
 * m_kg <= 0 or NaN -> 0.0. */
double stc_sr_momentum(double m_kg, double v_ms);

/* Energy-momentum relation: E = sqrt((p*c)^2 + (m*c^2)^2), in J.
 * Computed with hypot(); m_kg <= 0 or NaN -> treated as 0. */
double stc_sr_energy_from_momentum(double m_kg, double p_si);

/* Invariant mass: m = sqrt(E^2-(p*c)^2)/c^2, in kg.
 * E <= 0, spacelike input (p*c >= E), or indeterminate inf/inf -> 0.0. */
double stc_sr_invariant_mass(double e_J, double p_si);

/* Rapidity: phi = artanh(beta) = artanh(v/c), dimensionless.
 * |v| >= c -> +/-HUGE_VAL (sign of v); NaN -> 0.0. */
double stc_sr_rapidity(double v_ms);

/* Inverse rapidity: v = c*tanh(phi), in m/s. NaN -> 0.0. */
double stc_sr_rapidity_to_velocity(double phi);

/* Twin paradox: out-and-back over one-way Earth distance at cruise speed.
 * distance_m: one-way Earth-frame distance (m); v_ms: cruise speed (m/s,
 * magnitude is used, guarded to < c). Instant turnaround; proper clock
 * on the traveler reads earth_time/gamma. distance <= 0, NaN, or
 * v == 0/NaN -> zero trip (gamma = 1). */
StcSrTwin stc_sr_twin(double distance_m, double v_ms);

/* Barn-pole paradox: pole rushing through a barn at speed |v|.
 * barn_length_m / pole_length_m: rest lengths (m).
 * simultaneity_gap = gamma*|v|*barn_length/c^2 (relativity of
 * simultaneity between the barn doors in the pole frame), in s.
 * Non-positive/NaN lengths -> zero result with gamma = 1, fits = 0. */
StcSrPole stc_sr_pole(double barn_length_m, double pole_length_m,
                      double v_ms);

/* Ladder paradox, alias of stc_sr_pole (ladder == pole). */
StcSrPole stc_sr_ladder(double barn_length_m, double ladder_length_m,
                        double v_ms);

/* Barn-pole paradox, alias of stc_sr_pole. */
StcSrPole stc_sr_barn_pole(double barn_length_m, double pole_length_m,
                           double v_ms);

/* Worldline proper time (integrated clock):
 * tau = sum over segments of sqrt(dt^2-(dx/c)^2), in s.
 * t_s/x_m: n samples (s, m); fixed loop of n-1 segments.
 * NULL buffers or n < 2 -> 0.0; spacelike/null/past-directed or
 * non-finite segments contribute 0. Never NaN. */
double stc_sr_worldline_proper_time(const double *t_s, const double *x_m,
                                    int n);

#endif /* STC_SR_H */
