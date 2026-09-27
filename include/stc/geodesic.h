#ifndef STC_GEODESIC_H
#define STC_GEODESIC_H

/* SPACETIME LAB — Schwarzschild equatorial geodesics in u = 1/r.
 * Master equation (relativistic Binet / orbit equation), RK4 in phi:
 *   u'' + u = c0 + 1.5 * rs * u^2
 *   photon (null):   c0 = 0
 *   timelike:        c0 = GM / L^2   (L = specific angular momentum)
 * SI doubles, radians. phi increases along the trajectory. Pure C11.
 */

typedef enum {
    STC_GEO_PHOTON = 0,  /* null geodesic, init from impact parameter */
    STC_GEO_TIMELIKE = 1 /* massive particle, init from rp + ecc */
} StcGeoKind;

typedef enum {
    STC_GEO_ORBITING = 0, /* r inside (rs, r_escape): keep integrating */
    STC_GEO_ESCAPED = 1,  /* r >= r_escape (or u <= 0): left the system */
    STC_GEO_CAPTURED = 2  /* r <= rs (horizon crossing / plunge) */
} StcGeoStatus;

typedef struct {
    StcGeoKind kind; /* photon vs timelike */
    double rs;       /* Schwarzschild radius, meters */
    double gm;       /* GM, m^3/s^2 (timelike input; photon: rs*c^2/2) */
    double angmom;   /* timelike specific L, m^2/s; photon: NaN */
    double impact;   /* photon impact parameter b, meters; timelike: NaN */
    double c0;       /* orbit-eq constant: 0 (photon) or GM/L^2, 1/m */
    double u;        /* 1/r, 1/m */
    double up;       /* du/dphi, 1/m */
    double phi;      /* azimuth, radians */
} StcGeodesic;

/* Photon first integral (null energy equation):
 *   (du/dphi)^2 = 1/b^2 - u^2 + rs*u^3.
 * Inbound photon at radius r_start (m) with impact parameter b (m),
 * rs (m), starting azimuth phi_start (rad). Inbound => up = +sqrt(...).
 * Returns 0 on success, -1 on bad input (NULL g, non-finite/non-positive
 * rs/b, r_start <= max(rs, b), or negative radicand); on failure g is
 * reset to a safe all-NaN photon state. */
int stc_geodesic_photon_init(StcGeodesic *g, double rs, double b,
                             double r_start, double phi_start);

/* Keplerian periapsis elements (semi-latus rectum p = rp*(1+ecc),
 * L = sqrt(GM*p)) as relativistic initial data: u = 1/rp, up = 0
 * at phi_start (periapsis passage). gm (m^3/s^2), rs (m, >= 0;
 * rs = 0 gives the Newtonian orbit), rp (m, must exceed rs),
 * ecc (dimensionless, >= 0; < 1 bound, >= 1 unbound).
 * Returns 0 on success, -1 on bad input (g reset to safe state). */
int stc_geodesic_orbit_init(StcGeodesic *g, double gm, double rs,
                            double rp, double ecc, double phi_start);

/* Orbit-equation RHS (relativistic Binet equation):
 *   f(u) = c0 - u + 1.5*rs*u^2   (derivatives d/dphi).
 * NULL g or non-finite u => NaN. */
double stc_geodesic_rhs(const StcGeodesic *g, double u);

/* One classical 4th-order Runge-Kutta step in phi (fixed step dphi, rad).
 * No-op on NULL g or non-finite dphi. */
void stc_geodesic_step(StcGeodesic *g, double dphi);

/* Fixed-N RK4 integration: at most `steps` steps of `dphi` (rad),
 * stopping early on escape/capture. r_hist/phi_hist (each `steps`
 * long) receive post-step (r, phi) samples, either may be NULL.
 * r_escape (m) <= 0 disables the escape check. Deterministic:
 * identical input => identical output; the step loop is capped at
 * `steps` (no error-driven adaptivity).
 * Returns the terminal StcGeoStatus; NULL g => STC_GEO_ORBITING. */
StcGeoStatus stc_geodesic_integrate(StcGeodesic *g, double dphi, int steps,
                                    double r_escape,
                                    double *r_hist, double *phi_hist);

/* Current radius r = 1/u (m); +Inf when u <= 0 (at infinity);
 * NULL g => NaN. */
double stc_geodesic_radius(const StcGeodesic *g);

/* Termination flags: r <= rs => STC_GEO_CAPTURED, r >= r_escape
 * (when r_escape > 0) or u <= 0 => STC_GEO_ESCAPED, else
 * STC_GEO_ORBITING. NULL g => STC_GEO_ORBITING. */
StcGeoStatus stc_geodesic_status(const StcGeodesic *g, double r_escape);

/* Conserved specific energy (Schwarzschild energy equation):
 *   E^2 = c^2*L^2*(up^2 + (1-rs*u)*(c^2/L^2 + u^2)).
 * Timelike => E in m^2/s^2 (J/kg); photon or NULL or negative
 * radicand => NaN. */
double stc_geodesic_energy(const StcGeodesic *g);

/* Conserved specific angular momentum L (m^2/s). Timelike => L
 * from orbit_init; photon or NULL => NaN. */
double stc_geodesic_angmom(const StcGeodesic *g);

/* Conserved photon impact parameter b (m). Photon => b;
 * timelike or NULL => NaN. */
double stc_geodesic_impact(const StcGeodesic *g);

/* Total light deflection (geodesic photon bending):
 * integrate an inbound photon (rs, b, r_start as in photon_init)
 * and return swept_phi - 2*acos(b/r_start) (rad), the finite-start
 * corrected deflection (exact Schwarzschild, all orders).
 * Deterministic: at most max_steps of dphi (documented cap).
 * NaN on bad input, capture, or failure to escape within the cap. */
double stc_geodesic_photon_deflection(double rs, double b, double r_start,
                                      double dphi, int max_steps);

/* Light-cone tipping (rendering aid, Schwarzschild light cones):
 *   tilt = (pi/4) * (rs/r), clamped to [0, pi/4].
 * 0 far away (upright cone), pi/4 at r <= rs (fully tipped inward).
 * Returns radians; NaN on non-finite input, r <= 0, or rs < 0. */
double stc_geodesic_cone_tilt(double r, double rs);

#endif /* STC_GEODESIC_H */
