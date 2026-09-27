#include "stc/sr.h"
#include "stc/units.h"

#include <float.h>
#include <math.h>
#include <stddef.h>

/* Map any velocity (m/s) to a guarded beta with |b| <= STC_SR_BETA_MAX. */
static double sr_beta_guarded(double v_ms)
{
    double b;
    if (isnan(v_ms)) {
        return 0.0;
    }
    if (isinf(v_ms)) {
        return (v_ms > 0.0) ? STC_SR_BETA_MAX : -STC_SR_BETA_MAX;
    }
    b = v_ms / STC_C;
    if (b > STC_SR_BETA_MAX) {
        return STC_SR_BETA_MAX;
    }
    if (b < -STC_SR_BETA_MAX) {
        return -STC_SR_BETA_MAX;
    }
    return b;
}

/* Gamma from an already-guarded beta; always finite, >= 1. */
static double sr_gamma_guarded(double b)
{
    return 1.0 / sqrt(1.0 - b * b);
}

/* Mass sanitizer: NaN/non-positive -> 0; +inf passes through. */
static double sr_mass(double m_kg)
{
    if (isnan(m_kg) || m_kg <= 0.0) {
        return 0.0;
    }
    return m_kg;
}

double stc_sr_gamma(double v_ms)
{
    double b;
    if (isnan(v_ms)) {
        return 1.0;
    }
    if (isinf(v_ms)) {
        return HUGE_VAL;
    }
    b = fabs(v_ms) / STC_C;
    if (b >= 1.0) {
        return HUGE_VAL;
    }
    if (b == 0.0) {
        return 1.0;
    }
    return 1.0 / sqrt(1.0 - b * b);
}

double stc_sr_gamma_from_beta(double beta)
{
    double b;
    if (isnan(beta)) {
        return 1.0;
    }
    if (isinf(beta)) {
        return HUGE_VAL;
    }
    b = fabs(beta);
    if (b >= 1.0) {
        return HUGE_VAL;
    }
    if (b == 0.0) {
        return 1.0;
    }
    return 1.0 / sqrt(1.0 - b * b);
}

double stc_sr_gamma_from_rapidity(double phi)
{
    if (isnan(phi)) {
        return 1.0;
    }
    return cosh(phi);
}

double stc_sr_beta(double v_ms)
{
    double b;
    if (isnan(v_ms)) {
        return 0.0;
    }
    if (isinf(v_ms)) {
        return (v_ms > 0.0) ? 1.0 : -1.0;
    }
    b = v_ms / STC_C;
    if (b > 1.0) {
        return 1.0;
    }
    if (b < -1.0) {
        return -1.0;
    }
    return b;
}

double stc_sr_velocity_from_gamma(double gamma)
{
    double inv;
    if (isnan(gamma) || gamma <= 1.0) {
        return 0.0;
    }
    if (isinf(gamma)) {
        return STC_C;
    }
    inv = 1.0 / gamma;
    return STC_C * sqrt(1.0 - inv * inv);
}

StcSrEvent stc_sr_boost(StcSrEvent e, double v_ms)
{
    StcSrEvent r;
    double t, x, b, g, vc, c2;
    t = isfinite(e.t_s) ? e.t_s : 0.0;
    x = isfinite(e.x_m) ? e.x_m : 0.0;
    b = sr_beta_guarded(v_ms);
    g = sr_gamma_guarded(b);
    vc = b * STC_C;
    c2 = STC_C * STC_C;
    r.t_s = g * (t - vc * x / c2);
    r.x_m = g * (x - vc * t);
    return r;
}

void stc_sr_transform(double t_s, double x_m, double v_ms,
                      double *out_t, double *out_x)
{
    StcSrEvent e, r;
    e.t_s = t_s;
    e.x_m = x_m;
    r = stc_sr_boost(e, v_ms);
    if (out_t != NULL) {
        *out_t = r.t_s;
    }
    if (out_x != NULL) {
        *out_x = r.x_m;
    }
}

StcSrEvent stc_sr_lorentz(StcSrEvent e, double v_ms)
{
    return stc_sr_boost(e, v_ms);
}

double stc_sr_interval(double dt_s, double dx_m)
{
    double q, adt, aq, d2;
    if (isnan(dt_s)) {
        dt_s = 0.0;
    }
    if (isnan(dx_m)) {
        dx_m = 0.0;
    }
    if (isinf(dt_s)) {
        dt_s = (dt_s > 0.0) ? DBL_MAX : -DBL_MAX;
    }
    if (isinf(dx_m)) {
        dx_m = (dx_m > 0.0) ? DBL_MAX : -DBL_MAX;
    }
    q = dx_m / STC_C;
    adt = fabs(dt_s);
    aq = fabs(q);
    if (adt > 1e140 || aq > 1e140) {
        /* Astronomic inputs: keep the sign, saturate the magnitude. */
        if (adt == aq) {
            return 0.0;
        }
        return (adt > aq) ? HUGE_VAL : -HUGE_VAL;
    }
    d2 = dt_s * dt_s - q * q;
    return d2 * (STC_C * STC_C);
}

double stc_sr_interval_events(StcSrEvent a, StcSrEvent b)
{
    double dt = a.t_s - b.t_s;
    double dx = a.x_m - b.x_m;
    return stc_sr_interval(dt, dx);
}

double stc_sr_add_velocity(double u_ms, double v_ms)
{
    double bu = sr_beta_guarded(u_ms);
    double bv = sr_beta_guarded(v_ms);
    double w = (bu + bv) / (1.0 + bu * bv);
    return w * STC_C;
}

double stc_sr_velocity_add(double u_ms, double v_ms)
{
    return stc_sr_add_velocity(u_ms, v_ms);
}

double stc_sr_doppler_factor(double v_ms)
{
    double b = sr_beta_guarded(v_ms);
    return sqrt((1.0 - b) / (1.0 + b));
}

double stc_sr_doppler(double f0_hz, double v_ms)
{
    if (isnan(f0_hz) || f0_hz <= 0.0) {
        return 0.0;
    }
    if (isinf(f0_hz)) {
        return HUGE_VAL;
    }
    return f0_hz * stc_sr_doppler_factor(v_ms);
}

double stc_sr_doppler_transverse(double f0_hz, double v_ms)
{
    double b, g;
    if (isnan(f0_hz) || f0_hz <= 0.0) {
        return 0.0;
    }
    if (isinf(f0_hz)) {
        return HUGE_VAL;
    }
    b = sr_beta_guarded(v_ms);
    g = sr_gamma_guarded(b);
    return f0_hz / g;
}

double stc_sr_redshift(double v_ms)
{
    double b = sr_beta_guarded(v_ms);
    return sqrt((1.0 + b) / (1.0 - b)) - 1.0;
}

double stc_sr_time_dilation(double dt_proper_s, double v_ms)
{
    double b, g;
    if (isnan(dt_proper_s)) {
        return 0.0;
    }
    if (isinf(dt_proper_s)) {
        return dt_proper_s;
    }
    b = sr_beta_guarded(v_ms);
    g = sr_gamma_guarded(b);
    return dt_proper_s * g;
}

double stc_sr_proper_time(double dt_lab_s, double v_ms)
{
    double b, g;
    if (isnan(dt_lab_s)) {
        return 0.0;
    }
    if (isinf(dt_lab_s)) {
        return dt_lab_s;
    }
    b = sr_beta_guarded(v_ms);
    g = sr_gamma_guarded(b);
    return dt_lab_s / g;
}

double stc_sr_length_contraction(double l_proper_m, double v_ms)
{
    double b, g;
    if (isnan(l_proper_m) || l_proper_m <= 0.0) {
        return 0.0;
    }
    if (isinf(l_proper_m)) {
        return HUGE_VAL;
    }
    b = sr_beta_guarded(v_ms);
    g = sr_gamma_guarded(b);
    return l_proper_m / g;
}

double stc_sr_rest_energy(double m_kg)
{
    double m = sr_mass(m_kg);
    if (m == 0.0) {
        return 0.0;
    }
    return m * STC_C * STC_C;
}

double stc_sr_energy(double m_kg, double v_ms)
{
    double m = sr_mass(m_kg);
    double b, g;
    if (m == 0.0) {
        return 0.0;
    }
    b = sr_beta_guarded(v_ms);
    g = sr_gamma_guarded(b);
    return g * m * STC_C * STC_C;
}

double stc_sr_kinetic_energy(double m_kg, double v_ms)
{
    double m = sr_mass(m_kg);
    double b, g, gm1;
    if (m == 0.0) {
        return 0.0;
    }
    b = sr_beta_guarded(v_ms);
    g = sr_gamma_guarded(b);
    gm1 = b * b * g * g / (1.0 + g);
    return gm1 * m * STC_C * STC_C;
}

double stc_sr_momentum(double m_kg, double v_ms)
{
    double m = sr_mass(m_kg);
    double b = sr_beta_guarded(v_ms);
    if (m == 0.0 || b == 0.0) {
        return 0.0;
    }
    return sr_gamma_guarded(b) * m * (b * STC_C);
}

double stc_sr_energy_from_momentum(double m_kg, double p_si)
{
    double m = sr_mass(m_kg);
    double pc, mc2;
    if (isnan(p_si)) {
        p_si = 0.0;
    }
    if (isinf(p_si)) {
        return HUGE_VAL;
    }
    if (m == 0.0 && p_si == 0.0) {
        return 0.0;
    }
    pc = p_si * STC_C;
    mc2 = m * STC_C * STC_C;
    return hypot(pc, mc2);
}

double stc_sr_invariant_mass(double e_J, double p_si)
{
    double pc;
    if (isnan(e_J) || e_J <= 0.0) {
        return 0.0;
    }
    if (isnan(p_si)) {
        p_si = 0.0;
    }
    if (!isfinite(p_si)) {
        return 0.0;
    }
    if (isinf(e_J)) {
        return HUGE_VAL;
    }
    pc = fabs(p_si) * STC_C;
    if (e_J <= pc) {
        return 0.0;
    }
    return sqrt((e_J - pc) * (e_J + pc)) / (STC_C * STC_C);
}

double stc_sr_rapidity(double v_ms)
{
    double b;
    if (isnan(v_ms)) {
        return 0.0;
    }
    if (isinf(v_ms)) {
        return (v_ms > 0.0) ? HUGE_VAL : -HUGE_VAL;
    }
    b = v_ms / STC_C;
    if (b >= 1.0) {
        return HUGE_VAL;
    }
    if (b <= -1.0) {
        return -HUGE_VAL;
    }
    return atanh(b);
}

double stc_sr_rapidity_to_velocity(double phi)
{
    if (isnan(phi)) {
        return 0.0;
    }
    return STC_C * tanh(phi);
}

StcSrTwin stc_sr_twin(double distance_m, double v_ms)
{
    StcSrTwin r = { 0.0, 0.0, 1.0, 0.0 };
    double speed, b, g, vc;
    if (isnan(distance_m) || distance_m <= 0.0) {
        return r;
    }
    if (isnan(v_ms)) {
        return r;
    }
    speed = fabs(v_ms);
    if (speed <= 0.0) {
        return r;
    }
    b = speed / STC_C;
    if (!isfinite(b) || b > STC_SR_BETA_MAX) {
        b = STC_SR_BETA_MAX;
    }
    vc = b * STC_C;
    g = sr_gamma_guarded(b);
    r.v_used_ms = vc;
    r.gamma = g;
    r.earth_time_s = 2.0 * distance_m / vc;
    r.traveler_time_s = r.earth_time_s / g;
    return r;
}

StcSrPole stc_sr_pole(double barn_length_m, double pole_length_m,
                      double v_ms)
{
    StcSrPole r = { 1.0, 0.0, 0.0, 0.0, 0, 0 };
    double b, g, vc;
    if (isnan(barn_length_m) || barn_length_m <= 0.0) {
        return r;
    }
    if (isnan(pole_length_m) || pole_length_m <= 0.0) {
        return r;
    }
    b = sr_beta_guarded(v_ms);
    if (b < 0.0) {
        b = -b;
    }
    g = sr_gamma_guarded(b);
    vc = b * STC_C;
    r.gamma = g;
    r.pole_in_barn_m = pole_length_m / g;
    r.barn_in_pole_m = barn_length_m / g;
    if (vc == 0.0) {
        r.simultaneity_gap_s = 0.0;
    } else {
        r.simultaneity_gap_s = g * vc * barn_length_m / (STC_C * STC_C);
    }
    r.fits_barn_frame = (r.pole_in_barn_m <= barn_length_m * 1.000000001) ? 1 : 0;
    r.fits_pole_frame = (pole_length_m <= r.barn_in_pole_m * 1.000000001) ? 1 : 0;
    return r;
}

StcSrPole stc_sr_ladder(double barn_length_m, double ladder_length_m,
                        double v_ms)
{
    return stc_sr_pole(barn_length_m, ladder_length_m, v_ms);
}

StcSrPole stc_sr_barn_pole(double barn_length_m, double pole_length_m,
                           double v_ms)
{
    return stc_sr_pole(barn_length_m, pole_length_m, v_ms);
}

double stc_sr_worldline_proper_time(const double *t_s, const double *x_m,
                                    int n)
{
    double tau = 0.0;
    int i;
    if (t_s == NULL || x_m == NULL || n < 2) {
        return 0.0;
    }
    for (i = 0; i < n - 1; i++) {
        double dt = t_s[i + 1] - t_s[i];
        double dx = x_m[i + 1] - x_m[i];
        double q, rr, s;
        if (!isfinite(dt) || !isfinite(dx)) {
            continue;
        }
        if (dt <= 0.0) {
            continue;
        }
        q = dx / STC_C;
        if (fabs(q) >= dt) {
            continue;
        }
        rr = q / dt;
        s = 1.0 - rr * rr;
        if (!(s > 0.0)) {
            continue;
        }
        tau += dt * sqrt(s);
    }
    return tau;
}
