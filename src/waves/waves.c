#include "stc/waves.h"

#include "stc/units.h"
#include "stc/vec.h"

#include <math.h>
#include <stddef.h>
#include <stdint.h>

double stc_waves_chirp_mass(double m1_kg, double m2_kg)
{
    if (!isfinite(m1_kg) || !isfinite(m2_kg) || m1_kg <= 0.0 || m2_kg <= 0.0) {
        return NAN;
    }
    /* Mc = (m1*m2)^(3/5) / (m1+m2)^(1/5). */
    return pow(m1_kg * m2_kg, 3.0 / 5.0) / pow(m1_kg + m2_kg, 1.0 / 5.0);
}

double stc_waves_strain_amplitude(double chirp_mass_kg, double freq_hz,
                                  double dist_m)
{
    if (!isfinite(chirp_mass_kg) || !isfinite(freq_hz) || !isfinite(dist_m) ||
        chirp_mass_kg <= 0.0 || freq_hz <= 0.0 || dist_m <= 0.0) {
        return NAN;
    }
    {
        double c = STC_C;
        double gm = STC_G * chirp_mass_kg / (c * c);
        double pf = STC_PI * freq_hz / c;

        /* h0 = (4/d) * (G*Mc/c^2)^(5/3) * (pi*f/c)^(2/3). */
        return (4.0 / dist_m) * pow(gm, 5.0 / 3.0) * pow(pf, 2.0 / 3.0);
    }
}

double stc_waves_quadrupole_strain(double m1_kg, double m2_kg,
                                   double separation_m, double dist_m)
{
    if (!isfinite(m1_kg) || !isfinite(m2_kg) || !isfinite(separation_m) ||
        !isfinite(dist_m) || m1_kg <= 0.0 || m2_kg <= 0.0 ||
        separation_m <= 0.0 || dist_m <= 0.0) {
        return NAN;
    }
    {
        double c = STC_C;
        double c4 = c * c * c * c;

        /* h = 4*G^2*m1*m2 / (c^4*a*R). */
        return 4.0 * STC_G * STC_G * m1_kg * m2_kg / (c4 * separation_m * dist_m);
    }
}

double stc_waves_chirp_freq(double chirp_mass_kg, double t_to_coal_s)
{
    if (!isfinite(chirp_mass_kg) || !isfinite(t_to_coal_s) ||
        chirp_mass_kg <= 0.0 || t_to_coal_s <= 0.0) {
        return NAN;
    }
    {
        double c = STC_C;
        double tc = STC_G * chirp_mass_kg / (c * c * c);

        /* f(tau) = (1/pi) * (5/(256*tau))^(3/8) * (G*Mc/c^3)^(-5/8). */
        return pow(5.0 / (256.0 * t_to_coal_s), 3.0 / 8.0) * pow(tc, -5.0 / 8.0) /
               STC_PI;
    }
}

double stc_waves_time_to_coalescence(double chirp_mass_kg, double freq_hz)
{
    if (!isfinite(chirp_mass_kg) || !isfinite(freq_hz) || chirp_mass_kg <= 0.0 ||
        freq_hz <= 0.0) {
        return NAN;
    }
    {
        double c = STC_C;
        double tc = STC_G * chirp_mass_kg / (c * c * c);

        /* tau(f) = (5/256) * (G*Mc/c^3)^(-5/3) * (pi*f)^(-8/3). */
        return (5.0 / 256.0) * pow(tc, -5.0 / 3.0) *
               pow(STC_PI * freq_hz, -8.0 / 3.0);
    }
}

double stc_waves_phase(double chirp_mass_kg, double t_to_coal_s,
                       double phase_c_rad)
{
    if (!isfinite(chirp_mass_kg) || !isfinite(t_to_coal_s) ||
        !isfinite(phase_c_rad) || chirp_mass_kg <= 0.0 || t_to_coal_s <= 0.0) {
        return NAN;
    }
    {
        double c = STC_C;
        double tc = STC_G * chirp_mass_kg / (c * c * c);

        /* Phi(tau) = Phi_c - 2*(tau/(5*G*Mc/c^3))^(5/8). */
        return phase_c_rad - 2.0 * pow(t_to_coal_s / (5.0 * tc), 5.0 / 8.0);
    }
}

double stc_waves_phase0_from_seed(uint64_t seed)
{
    /* splitmix64 finalizer (deterministic integer hash, no state). */
    uint64_t z = seed + 0x9e3779b97f4a7c15ULL;

    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
    z = z ^ (z >> 31);
    /* Top 53 bits -> u in [0,1), then scale to [0, 2*pi). */
    return (double)(z >> 11) * (1.0 / 9007199254740992.0) * 2.0 * STC_PI;
}

void stc_waves_ring_response(double h_plus, double h_cross, StcVec3 p0,
                             StcVec3 *out)
{
    StcVec3 r;

    if (out == NULL) {
        return;
    }
    if (!isfinite(h_plus) || !isfinite(h_cross)) {
        *out = p0;
        return;
    }
    /* TT gauge, +z propagation: transverse-traceless displacement. */
    r.x = p0.x + 0.5 * (h_plus * p0.x + h_cross * p0.y);
    r.y = p0.y + 0.5 * (h_cross * p0.x - h_plus * p0.y);
    r.z = p0.z;
    *out = r;
}

void stc_waves_chirp_waveform(double chirp_mass_kg, double dist_m,
                              double t_to_coal_s, double phase_c_rad,
                              double inclination_rad, double *h_plus,
                              double *h_cross)
{
    double freq;
    double phi;
    double h0;
    double ci;

    if (h_plus == NULL && h_cross == NULL) {
        return;
    }
    freq = stc_waves_chirp_freq(chirp_mass_kg, t_to_coal_s);
    phi = stc_waves_phase(chirp_mass_kg, t_to_coal_s, phase_c_rad);
    if (!isfinite(freq) || !isfinite(phi) || !isfinite(dist_m) ||
        !isfinite(inclination_rad) || dist_m <= 0.0) {
        if (h_plus != NULL) {
            *h_plus = NAN;
        }
        if (h_cross != NULL) {
            *h_cross = NAN;
        }
        return;
    }
    h0 = stc_waves_strain_amplitude(chirp_mass_kg, freq, dist_m);
    if (!isfinite(h0)) {
        if (h_plus != NULL) {
            *h_plus = NAN;
        }
        if (h_cross != NULL) {
            *h_cross = NAN;
        }
        return;
    }
    ci = cos(inclination_rad);
    if (h_plus != NULL) {
        *h_plus = h0 * (1.0 + ci * ci) * 0.5 * cos(phi);
    }
    if (h_cross != NULL) {
        *h_cross = h0 * ci * sin(phi);
    }
}
