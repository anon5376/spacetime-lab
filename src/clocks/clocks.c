/* SPACETIME LAB — clocks & classic tests of GR (implementation).
 * Own header + stc/units.h (via clocks.h) + math.h only. Pure functions,
 * no state, no I/O. Error policy: NaN (or NaN-filled struct) on invalid
 * input as documented per function; libs never exit().
 */

#include "stc/clocks.h"

#include <math.h>

double stc_clocks_proper_rate(double r_m, double rs_m)
{
    if (!(r_m > 0.0) || !(rs_m >= 0.0)) {
        return NAN;
    }
    if (r_m <= rs_m) {
        return 0.0;
    }
    return sqrt(1.0 - rs_m / r_m);
}

double stc_clocks_redshift(double r_m, double rs_m)
{
    double g = stc_clocks_proper_rate(r_m, rs_m);

    if (!(g >= 0.0)) {
        return NAN;
    }
    if (g == 0.0) {
        return INFINITY;
    }
    return 1.0 / g - 1.0;
}

double stc_clocks_circular_velocity(double gm_si, double r_m)
{
    if (!(gm_si >= 0.0) || !(r_m > 0.0)) {
        return NAN;
    }
    return sqrt(gm_si / r_m);
}

double stc_clocks_sr_fraction(double v_sat_m_s, double v_gnd_m_s)
{
    double c2 = STC_C * STC_C;

    return (v_gnd_m_s * v_gnd_m_s - v_sat_m_s * v_sat_m_s) / (2.0 * c2);
}

double stc_clocks_gr_fraction(double gm_si, double r_sat_m, double r_gnd_m)
{
    double rs, g_sat, g_gnd;

    if (!(gm_si >= 0.0) || !(r_sat_m > 0.0) || !(r_gnd_m > 0.0)) {
        return NAN;
    }
    rs = 2.0 * gm_si / (STC_C * STC_C);
    g_sat = stc_clocks_proper_rate(r_sat_m, rs);
    g_gnd = stc_clocks_proper_rate(r_gnd_m, rs);
    if (!(g_sat > 0.0) || !(g_gnd > 0.0)) {
        return NAN;
    }
    return g_sat / g_gnd - 1.0;
}

double stc_clocks_gps_sr_us_per_day(void)
{
    double gm = STC_G * STC_M_EARTH;
    double v_sat = stc_clocks_circular_velocity(gm, STC_CLOCKS_GPS_RADIUS);
    double v_gnd = STC_CLOCKS_OMEGA_EARTH * STC_CLOCKS_R_EARTH;

    return stc_clocks_sr_fraction(v_sat, v_gnd) * STC_DAY * 1.0e6;
}

double stc_clocks_gps_gr_us_per_day(void)
{
    double gm = STC_G * STC_M_EARTH;
    double frac = stc_clocks_gr_fraction(gm, STC_CLOCKS_GPS_RADIUS,
                                         STC_CLOCKS_R_EARTH);

    return frac * STC_DAY * 1.0e6;
}

double stc_clocks_gps_net_us_per_day(void)
{
    return stc_clocks_gps_sr_us_per_day() + stc_clocks_gps_gr_us_per_day();
}

double stc_clocks_shapiro_one_way(double gm_si, double r_emit_m,
                                   double r_recv_m, double sep_m)
{
    double sum, diff, c3;

    if (!(gm_si >= 0.0) || !(r_emit_m > 0.0) || !(r_recv_m > 0.0)) {
        return NAN;
    }
    sum = r_emit_m + r_recv_m;
    diff = fabs(r_emit_m - r_recv_m);
    if (!(sep_m > diff) || !(sep_m < sum)) {
        return NAN;
    }
    c3 = STC_C * STC_C * STC_C;
    return (2.0 * gm_si / c3) * log((sum + sep_m) / (sum - sep_m));
}

double stc_clocks_shapiro_grazing(double gm_si, double r1_m, double r2_m,
                                   double b_m, int round_trip)
{
    double rmin, pre, c3;

    if (!(gm_si >= 0.0) || !(r1_m > 0.0) || !(r2_m > 0.0) ||
        !(b_m > 0.0)) {
        return NAN;
    }
    rmin = r1_m < r2_m ? r1_m : r2_m;
    if (!(b_m < rmin)) {
        return NAN;
    }
    c3 = STC_C * STC_C * STC_C;
    pre = round_trip ? 4.0 : 2.0;
    return (pre * gm_si / c3) * log(4.0 * r1_m * r2_m / (b_m * b_m));
}

double stc_clocks_shapiro_sun_grazing_us(int round_trip)
{
    double gm = STC_G * STC_M_SUN;
    double s = stc_clocks_shapiro_grazing(gm, STC_AU, STC_CLOCKS_VENUS_R,
                                           STC_R_SUN, round_trip);

    return s * 1.0e6;
}

double stc_clocks_perihelion_shift(double gm_si, double a_m, double ecc)
{
    double c2;

    if (!(gm_si >= 0.0) || !(a_m > 0.0) || !(ecc >= 0.0) || !(ecc < 1.0)) {
        return NAN;
    }
    c2 = STC_C * STC_C;
    return 6.0 * STC_PI * gm_si / (c2 * a_m * (1.0 - ecc * ecc));
}

double stc_clocks_mercury_arcsec_per_century(void)
{
    double gm = STC_G * STC_M_SUN;
    double a = STC_CLOCKS_MERCURY_A;
    double per_orbit = stc_clocks_perihelion_shift(gm, a,
                                                   STC_CLOCKS_MERCURY_ECC);
    double period = 2.0 * STC_PI * sqrt(a * a * a / gm);
    double orbits = STC_CLOCKS_JULIAN_CENTURY_D * STC_DAY / period;

    return per_orbit * orbits / stc_arcsec();
}

StcHafeleKeating stc_clocks_hafele_keating(double altitude_m,
                                           double airspeed_m_s, int eastward)
{
    StcHafeleKeating out = { NAN, NAN, NAN, NAN };
    double gm, radius, v_gnd, v_plane, flight, c2, sr_frac, gr_frac;

    if (!(altitude_m >= 0.0) || !(airspeed_m_s > 0.0)) {
        return out;
    }
    gm = STC_G * STC_M_EARTH;
    radius = STC_CLOCKS_R_EARTH;
    v_gnd = STC_CLOCKS_OMEGA_EARTH * radius;
    v_plane = eastward ? v_gnd + airspeed_m_s : v_gnd - airspeed_m_s;
    flight = 2.0 * STC_PI * radius / airspeed_m_s;
    c2 = STC_C * STC_C;
    sr_frac = (v_gnd * v_gnd - v_plane * v_plane) / (2.0 * c2);
    gr_frac = (gm / (radius * radius)) * altitude_m / c2;
    out.sr_ns = sr_frac * flight * 1.0e9;
    out.gr_ns = gr_frac * flight * 1.0e9;
    out.net_ns = out.sr_ns + out.gr_ns;
    out.flight_s = flight;
    return out;
}

StcHafeleKeating stc_clocks_hafele_keating_nominal(int eastward)
{
    return stc_clocks_hafele_keating(STC_CLOCKS_HK_ALTITUDE_M,
                                     STC_CLOCKS_HK_AIRSPEED_M_S, eastward);
}
