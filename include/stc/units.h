#ifndef STC_UNITS_H
#define STC_UNITS_H

/* SPACETIME LAB — physical constants (SI) + small helpers.
 * Every constant: value + unit documented below. All doubles.
 */

/* Speed of light in vacuum: 299792458 m/s (exact, SI defining). */
#define STC_C 299792458.0

/* Newtonian gravitational constant: 6.67430e-11 m^3 kg^-1 s^-2. */
#define STC_G 6.67430e-11

/* Solar mass: 1.98847e30 kg (IAU 2015 nominal). */
#define STC_M_SUN 1.98847e30

/* Earth mass: 5.97237e24 kg. */
#define STC_M_EARTH 5.97237e24

/* Solar radius: 6.957e8 m (IAU 2015 nominal). */
#define STC_R_SUN 6.957e8

/* Earth mean radius: 6.371e6 m (WGS84-derived mean (2a+b)/3 = 6371008.8 m). */
#define STC_R_EARTH 6.371e6

/* Solar standard gravitational parameter GM_sun = STC_G * STC_M_SUN =
 * 1.3271645e20 m^3/s^2 (CODATA 2018 G composed with IAU 2015 nominal mass;
 * the IAU nominal GM_sun 1.3271244e20 differs by ~3e-5 relative). */
#define STC_GM_SUN (STC_G * STC_M_SUN)

/* Earth standard gravitational parameter GM_earth = STC_G * STC_M_EARTH =
 * 3.9861389e14 m^3/s^2 (CODATA/IAU composition). WGS84's 3.986004418e14
 * differs by ~3.4e-5 relative; all in-tree code uses this composed value. */
#define STC_GM_EARTH (STC_G * STC_M_EARTH)

/* Astronomical unit: 149597870700 m (exact, IAU 2012). */
#define STC_AU 149597870700.0

/* Julian year: 31557600 s (365.25 d, exact by definition). */
#define STC_YEAR 31557600.0

/* Mean solar day: 86400 s (exact by definition). */
#define STC_DAY 86400.0

/* Pi (dimensionless). */
#define STC_PI 3.14159265358979323846

/* Schwarzschild radius rs = 2*G*M/c^2. M in kg, returns meters. */
double stc_schwarzschild_radius(double m_kg);

/* Geometric mass M_geom = G*M/c^2. M in kg, returns meters. */
double stc_geometric_M(double m_kg);

/* One arcsecond in radians: pi / (180 * 3600). */
double stc_arcsec(void);

/* One microarcsecond in radians: arcsec / 1e6. */
double stc_microarcsec(void);

#endif /* STC_UNITS_H */
