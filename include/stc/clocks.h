#ifndef STC_CLOCKS_H
#define STC_CLOCKS_H

/* SPACETIME LAB — clocks & classic tests of GR.
 *
 * Schwarzschild proper-time rate, GPS clock offsets COMPUTED from orbit
 * mechanics (circular-orbit velocity + Schwarzschild rate ratio, never
 * hardcoded), Shapiro time delay, Einstein perihelion precession, and a
 * crude Hafele-Keating model. Pure C11, SI doubles, no state.
 *
 * Every function below names its equation. Time results carry their unit
 * in the function or field name (SI seconds unless suffixed).
 */

#include "stc/units.h"

/* Mean Earth radius: 6371000 m. */
#define STC_CLOCKS_R_EARTH 6371000.0

/* Earth rotation rate: 7.2921159e-5 rad/s (IERS nominal). */
#define STC_CLOCKS_OMEGA_EARTH 7.2921159e-5

/* GPS constellation nominal orbit radius (semi-major axis): 26560000 m. */
#define STC_CLOCKS_GPS_RADIUS 26560000.0

/* Mercury semi-major axis: 0.38710 AU, in meters. */
#define STC_CLOCKS_MERCURY_A (0.38710 * STC_AU)

/* Mercury orbital eccentricity: 0.20563 (dimensionless). */
#define STC_CLOCKS_MERCURY_ECC 0.20563

/* Venus orbital radius (circular approx): 0.72333 AU, in meters. */
#define STC_CLOCKS_VENUS_R (0.72333 * STC_AU)

/* Hafele-Keating nominal cruise altitude: 10000 m. */
#define STC_CLOCKS_HK_ALTITUDE_M 10000.0

/* Hafele-Keating nominal ground speed: 250 m/s. */
#define STC_CLOCKS_HK_AIRSPEED_M_S 250.0

/* Julian century: 36525 days (exact by definition). */
#define STC_CLOCKS_JULIAN_CENTURY_D 36525.0

/* Hafele-Keating prediction. Positive = plane clock gains vs ground. */
typedef struct {
    double sr_ns;    /* kinematic (SR) part, nanoseconds */
    double gr_ns;    /* gravitational (GR) part, nanoseconds */
    double net_ns;   /* sr_ns + gr_ns, nanoseconds */
    double flight_s; /* time aloft, SI seconds */
} StcHafeleKeating;

/* Schwarzschild proper-time rate: d(tau)/dt = sqrt(1 - rs/r).
 * r_m: areal radius in m; rs_m: Schwarzschild radius in m.
 * Returns dimensionless factor in [0,1]. r<=rs -> 0.0 (static clock
 * frozen; no static observers at/below the horizon). NaN if r<=0 or
 * rs<0. */
double stc_clocks_proper_rate(double r_m, double rs_m);

/* Gravitational redshift: z = 1/sqrt(1-rs/r) - 1 (dimensionless).
 * Derived from stc_clocks_proper_rate. +INFINITY at/below horizon,
 * NaN on the same invalid inputs as proper_rate. */
double stc_clocks_redshift(double r_m, double rs_m);

/* Newtonian circular-orbit velocity: v = sqrt(GM/r).
 * gm_si: G*M in m^3/s^2; r_m: orbit radius in m. Returns m/s.
 * NaN if gm<0 or r<=0. */
double stc_clocks_circular_velocity(double gm_si, double r_m);

/* Second-order time-dilation (transverse Doppler) difference:
 * dF = (v_gnd^2 - v_sat^2) / (2 c^2), dimensionless fractional
 * frequency offset of moving clock vs ground clock (negative = moving
 * clock runs slow). Inputs in m/s. */
double stc_clocks_sr_fraction(double v_sat_m_s, double v_gnd_m_s);

/* Schwarzschild gravitational rate ratio:
 * dF = sqrt(1-rs/r_sat)/sqrt(1-rs/r_gnd) - 1 with rs = 2GM/c^2,
 * dimensionless fractional offset of high clock vs ground clock
 * (positive = high clock runs fast). Distances in m, gm in m^3/s^2.
 * NaN if gm<0, either r<=0, or either r<=rs. */
double stc_clocks_gr_fraction(double gm_si, double r_sat_m, double r_gnd_m);

/* GPS special-relativistic offset, COMPUTED from orbit mechanics:
 * circular velocity at STC_CLOCKS_GPS_RADIUS vs equatorial ground
 * speed STC_CLOCKS_OMEGA_EARTH*STC_CLOCKS_R_EARTH through
 * stc_clocks_sr_fraction. Returns microseconds per day, expect ~-7.2. */
double stc_clocks_gps_sr_us_per_day(void);

/* GPS general-relativistic offset, COMPUTED from the Schwarzschild rate
 * ratio between STC_CLOCKS_GPS_RADIUS and STC_CLOCKS_R_EARTH.
 * Returns microseconds per day, expect ~+45.7. */
double stc_clocks_gps_gr_us_per_day(void);

/* GPS net clock offset = SR + GR parts above. Returns microseconds per
 * day, expect ~+38. */
double stc_clocks_gps_net_us_per_day(void);

/* Shapiro time delay, one-way PPN form (gamma = 1):
 * dt = (2GM/c^3) ln((r_emit + r_recv + sep)/(r_emit + r_recv - sep)).
 * r_emit/r_recv: emitter/receiver distances from the mass in m; sep:
 * straight-line emitter-receiver distance in m; gm in m^3/s^2.
 * Returns SI seconds. NaN unless |r_emit-r_recv| < sep < r_emit+r_recv
 * (triangle inequality) with gm>=0 and both r>0. */
double stc_clocks_shapiro_one_way(double gm_si, double r_emit_m,
                                   double r_recv_m, double sep_m);

/* Shapiro delay, Sun-grazing logarithmic approximation:
 * dt = (K GM/c^3) ln(4 r1 r2 / b^2), K = 2 one-way, K = 4 round-trip.
 * r1/r2: endpoint distances from the mass in m; b: ray impact
 * parameter in m; round_trip: 0 = one-way, nonzero = round-trip.
 * Returns SI seconds. NaN unless 0 < b < min(r1,r2) with gm>=0. */
double stc_clocks_shapiro_grazing(double gm_si, double r1_m, double r2_m,
                                   double b_m, int round_trip);

/* Canonical Shapiro case: Sun-grazing Earth-Venus link at superior
 * conjunction (r1 = 1 AU, r2 = Venus radius, b = solar radius).
 * round_trip: 0 = one-way (~116 us), nonzero = round-trip (~232 us).
 * Returns microseconds. */
double stc_clocks_shapiro_sun_grazing_us(int round_trip);

/* Einstein perihelion advance: dphi = 6 pi GM / (c^2 a (1-e^2)).
 * gm in m^3/s^2, a in m, ecc dimensionless. Returns radians per orbit.
 * NaN if gm<0, a<=0, or ecc outside [0,1). */
double stc_clocks_perihelion_shift(double gm_si, double a_m, double ecc);

/* Mercury perihelion advance per Julian century, COMPUTED: per-orbit
 * Einstein shift at Mercury's elements times orbits-per-century from
 * the Kepler-III period T = 2 pi sqrt(a^3/GM). Returns arcseconds per
 * Julian century, expect ~43. */
double stc_clocks_mercury_arcsec_per_century(void);

/* Hafele-Keating circumnavigation, crude equatorial model.
 * Kinematic part: second-order time dilation of ECI velocities
 * (v_gnd +/- airspeed) over the flight; gravitational part: gh/c^2
 * with surface g = GM/R^2 over the flight; flight time 2 pi R/u.
 * eastward: nonzero = eastward (plane faster in ECI, net loss),
 * 0 = westward (plane slower in ECI, net gain).
 * Approximations (documented): equatorial great-circle track, constant
 * altitude and ground speed, no stops/winds/latitude variation, ECI
 * comparison (Sagnac-corrected). Expect east ~-90 ns, west ~+330 ns
 * at nominal inputs (observed -59/+273 ns: signs and order right).
 * Returns NaN-filled struct if altitude<0 or airspeed<=0. */
StcHafeleKeating stc_clocks_hafele_keating(double altitude_m,
                                           double airspeed_m_s, int eastward);

/* Hafele-Keating at STC_CLOCKS_HK_ALTITUDE_M/_AIRSPEED_M_S. Same
 * convention: eastward nonzero = eastward. */
StcHafeleKeating stc_clocks_hafele_keating_nominal(int eastward);

#endif /* STC_CLOCKS_H */
