#ifndef STC_DRAG_H
#define STC_DRAG_H

/* SPACETIME LAB — drag: atmosphere density, aerodynamic drag, radiation
 * pressure. Pure C11, SI doubles, no state, no sim mutation, no I/O.
 * Every function below names its equation.
 */

/* Equation: barometric formula (isothermal exponential atmosphere) —
 * rho(h) = rho0 * exp(-h / H).
 * h: altitude above the reference level, m (any finite value; negative h
 * extrapolates above rho0); rho0: reference density at h = 0, kg/m^3
 * (must be finite and >= 0); H: pressure scale height, m (must be finite
 * and > 0). Returns density in kg/m^3. NaN if rho0/H invalid.
 * Approximation: isothermal, constant-g, single-species atmosphere. */
double stc_atmo_density(double h, double rho0, double H);

/* Equation: quadratic drag equation —
 * a = -(1/2) rho |v| Cd A / m * v.
 * v_vec: velocity relative to the atmosphere, m/s (3-vector); rho:
 * atmosphere density, kg/m^3 (>= 0); Cd: drag coefficient, dimensionless
 * (>= 0); A: reference area, m^2 (>= 0); m: body mass, kg (> 0);
 * out: acceleration, m/s^2 (3-vector). Zero velocity (or rho/Cd/A == 0)
 * yields the zero vector. NULL out is a no-op; NULL v_vec or invalid /
 * non-finite scalars fill out with NaN. */
void stc_drag_accel(const double v_vec[3], double rho, double Cd, double A,
                    double m, double out[3]);

/* Equation: radiation pressure (radial inverse-square) —
 * a = Q A L / (4 pi r^2 c m) * r_hat, directed away from the source.
 * r_vec_sun: body position relative to the Sun/source, m (3-vector,
 * must be nonzero); luminosity: source luminosity, W (>= 0); area:
 * exposed area, m^2 (>= 0); mass: body mass, kg (> 0); Q: radiation
 * pressure efficiency, dimensionless (>= 0; 1 = perfect absorber,
 * 2 = perfect normal reflector); out: acceleration, m/s^2 (3-vector).
 * NULL out is a no-op; NULL r_vec_sun, zero r, or invalid / non-finite
 * scalars fill out with NaN. */
void stc_radiation_pressure_accel(const double r_vec_sun[3],
                                  double luminosity, double area, double mass,
                                  double Q, double out[3]);

#endif /* STC_DRAG_H */
