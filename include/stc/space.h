#ifndef STC_SPACE_H
#define STC_SPACE_H

/* SPACETIME LAB — space: large-N gravitational dynamics core for the
 * space simulator (--space).
 *
 * Frame and units: SI (m, m/s, s, m^3/s^2). Positions/velocities are
 * barycentric, J2000 ecliptic axes (x toward the J2000 equinox, z toward
 * the ecliptic north pole). Time t is TDB seconds since J2000.0
 * (JD 2451545.0 TDB).
 *
 * Two particle classes:
 *   massive  (nm)  mutually interacting; each has GM, radius, kind.
 *   tracer   (nt)  massless test particles; feel massive bodies only
 *                  (asteroid belts, debris, "a million Earths" clouds).
 *
 * Dynamics:
 *   Newtonian pairwise gravity with optional Plummer softening (soft = 0
 *   means exact point masses — the default for real systems). Direct
 *   O(nm^2) summation up to tree_threshold massive bodies, Barnes-Hut
 *   octree (opening angle theta) above it. Tracers: direct over massive
 *   bodies, tree above the threshold.
 *   Optional 1PN correction (pn != 0) from every massive body flagged
 *   as a relativistic source (kind STAR or BLACKHOLE), applied to every
 *   other particle in the source's frame (harmonic-gauge Schwarzschild
 *   test-particle form):
 *     a_1PN = GM/(c^2 r^3) [ (4 GM/r - v^2) r_vec + 4 (r_vec . v_vec) v_vec ]
 *   with r_vec, v_vec relative to the source. Reaction on the source is
 *   neglected (documented). This reproduces Mercury's ~43"/century
 *   perihelion advance.
 *
 *   Optional zonal harmonics (J2, J4) of any massive body with
 *   oblateness set: the standard axisymmetric-potential acceleration
 *   about its spin pole, felt by every other particle (massive and
 *   tracer), with the equal and opposite reaction applied to the
 *   oblate body for massive targets (momentum conserving). Required for
 *   giant-planet moons: without Jupiter's J2, Io's phase is lost within
 *   a few years.
 *   1PN kicks evaluate the velocity-dependent term at the kick's
 *   implicit midpoint velocity (fixed-point iteration), keeping the
 *   composition's accuracy when 1PN is on.
 *
 * Integrator: fixed-substep 6th-order Yoshida composition of the
 * kick-drift-kick leapfrog (Yoshida 1990, Phys. Lett. A 150, 262).
 * With pn = 0 and no events it is symplectic and time-reversible:
 * advancing by +T then -T returns the initial state to roundoff.
 * Negative dt integrates backwards in time (physically valid for these
 * conservative dynamics). The velocity-dependent 1PN term breaks exact
 * reversibility at the 1e-8-relative level per orbit (measured, see
 * stspace); exact rewind is provided by checkpoints instead.
 *
 * Events (irreversible, logged): massive-massive contact -> perfectly
 * inelastic merge (mass, momentum, centre of mass conserved; radius from
 * summed volume, or the Schwarzschild radius if either is a black hole);
 * tracer hitting a massive body -> absorbed; any particle closer than
 * STC_SPACE_CAPTURE_RS Schwarzschild radii to a black hole -> captured;
 * massive body inside the tidal radius r_t = R (M_bh/m)^(1/3) of a black
 * hole (and outside capture) -> tidally disrupted: removed and replaced
 * by `disrupt_tracers` tracers spread over its volume with its velocity.
 *
 * Exact time travel: checkpoints every ckpt_interval seconds of
 * integration (ring of ckpt_max snapshots, both time directions).
 * stc_space_seek(s, t) restores the nearest checkpoint on the correct
 * side of t and re-integrates deterministically to t, so scrubbing back
 * across an event reproduces the pre-event state bitwise.
 *
 * Determinism: same inputs = bitwise identical state, independent of
 * thread count (force loops write per-index slots only; fixed
 * summation order).
 */

#include <stddef.h>

/* Particle kinds (massive). */
#define STC_SK_STAR 0
#define STC_SK_PLANET 1
#define STC_SK_DWARF 2
#define STC_SK_MOON 3
#define STC_SK_ASTEROID 4
#define STC_SK_BLACKHOLE 5
#define STC_SK_ROCK 6 /* generic user-spawned massive body */

/* Capture radius for black holes, in Schwarzschild radii (ILLUSTRATIVE
 * threshold: 2 r_s = 4 GM/c^2, inside the photon sphere). */
#define STC_SPACE_CAPTURE_RS 2.0

#define STC_SPACE_NAME 24

/* Event kinds. */
#define STC_SEV_MERGE 1    /* a absorbed b (both massive) */
#define STC_SEV_ABSORB 2   /* tracer b hit massive a */
#define STC_SEV_CAPTURE 3  /* black hole a captured massive b (b >= 0) or
                              tracer (-1 - tracer index) */
#define STC_SEV_DISRUPT 4  /* black hole a tidally disrupted massive b */

typedef struct {
    int kind;
    int a;
    int b;
    double t;
} StcSpaceEvent;

#define STC_SPACE_EVLOG 256

typedef struct StcSpaceCkpt StcSpaceCkpt; /* opaque, owned by space.c */

typedef struct {
    /* --- massive bodies (structure of arrays, capacity cap_m) --- */
    int nm;
    int cap_m;
    double *x, *y, *z;       /* m */
    double *vx, *vy, *vz;    /* m/s */
    double *gm;              /* m^3/s^2 */
    double *radius;          /* m (black hole: Schwarzschild radius) */
    int *kind;               /* STC_SK_* */
    int *parent;             /* display primary index or -1 */
    unsigned char *alive;    /* 0 after merge/capture/disrupt */
    char (*name)[STC_SPACE_NAME];
    /* Zonal gravity (oblateness) of body i, all 0 = point mass:
     * J2, J4 (unnormalised), reference (equatorial) radius rref in m,
     * spin-pole unit vector (polx, poly, polz) in the J2000 ecliptic
     * frame. Set with stc_space_set_oblateness. */
    double *j2, *j4, *rref;
    double *polx, *poly, *polz;

    /* --- tracers (capacity cap_t) --- */
    int nt;
    int cap_t;
    double *tx, *ty, *tz;
    double *tvx, *tvy, *tvz;
    unsigned char *talive;

    /* --- clock + integrator settings --- */
    double t;               /* TDB s since J2000.0 */
    double h;               /* max |substep|, s (> 0) */
    double soft;            /* Plummer softening, m (0 = none) */
    int pn;                 /* 1PN correction on/off */
    double theta;           /* Barnes-Hut opening angle */
    int tree_threshold;     /* use the tree when live nm exceeds this */
    int disrupt_tracers;    /* tracers spawned per tidal disruption */
    int collisions;         /* merge/absorb/capture/disrupt on/off */

    /* --- diagnostics (updated by stc_space_advance) --- */
    long steps;             /* substeps taken since init */
    double e0;              /* total energy at the last reset_diag */

    /* --- event log (ring) --- */
    StcSpaceEvent ev[STC_SPACE_EVLOG];
    int nev;                /* total events ever (ring index nev % LOG) */

    /* --- checkpoints (opaque) --- */
    StcSpaceCkpt *ck;
    int ckpt_max;
    double ckpt_interval;
} StcSpace;

/* Allocate capacities (cap_m >= 1, cap_t >= 0). Defaults: h = 3600 s,
 * soft = 0, pn = 1, theta = 0.5, tree_threshold = 1024,
 * disrupt_tracers = 400, collisions = 1, checkpoints 64 x 30 days.
 * Returns 0, -1 bad args, -2 allocation failure (s zeroed). */
int stc_space_init(StcSpace *s, int cap_m, int cap_t);

/* Free everything. NULL-safe. */
void stc_space_free(StcSpace *s);

/* Clear all particles, events, checkpoints; set t. Capacities kept. */
void stc_space_clear(StcSpace *s, double t);

/* Grow capacities (keeps contents). 0 ok, -2 allocation failure. */
int stc_space_reserve(StcSpace *s, int cap_m, int cap_t);

/* Add one massive body. pos/vel SI barycentric. gm > 0, radius >= 0.
 * Returns its index, or -1 on bad args / full (call reserve first). */
int stc_space_add(StcSpace *s, const char *name, int kind, int parent,
                  double gm, double radius, const double pos[3],
                  const double vel[3]);

/* Add one tracer. Returns its index or -1. */
int stc_space_add_tracer(StcSpace *s, const double pos[3],
                         const double vel[3]);

/* Turn massive body i into a black hole of the same GM (radius becomes
 * r_s = 2 GM/c^2, kind BLACKHOLE, name kept unless new_name != NULL).
 * gm_new > 0 also changes its GM (momentum kept: velocity rescaled by
 * gm_old/gm_new). Returns 0 or -1. */
int stc_space_make_blackhole(StcSpace *s, int i, double gm_new,
                             const char *new_name);

/* Set body i's zonal harmonics: j2, j4 (>= 0 typical), rref > 0 (m),
 * pole = spin axis (normalised internally; ecliptic J2000). j2 = j4 = 0
 * restores a point mass. Returns 0 or -1 on bad args. */
int stc_space_set_oblateness(StcSpace *s, int i, double j2, double j4,
                             double rref, const double pole[3]);

/* Remove the net barycentric momentum and centre-of-mass offset by a
 * uniform shift of all massive bodies and tracers. */
void stc_space_to_barycentre(StcSpace *s);

/* Advance by dt seconds (either sign) in ceil(|dt|/h) equal substeps,
 * but at most max_substeps (then stops early; returns the number taken
 * as a negative value so the caller can show "warp limited"). Processes
 * events after every substep and records checkpoints. Returns substeps
 * taken (>= 0) when dt was fully covered. */
long stc_space_advance(StcSpace *s, double dt, long max_substeps);

/* Exact time travel: restore the best checkpoint and re-integrate to t.
 * Returns 0 ok, -1 when t lies outside the checkpointed span (then the
 * state is unchanged; the caller may integrate there with advance). */
int stc_space_seek(StcSpace *s, double t);

/* Accelerations of all live massive bodies into ax/ay/az (length nm),
 * same force model as the integrator (Newton [+ tree] [+ 1PN]). Dead
 * slots get 0. Used by verification. */
void stc_space_accel(const StcSpace *s, double *ax, double *ay,
                     double *az);

/* Total Newtonian energy of the live massive bodies, joules:
 * sum m_i v_i^2 / 2 - sum_{i<j} G m_i m_j / r_ij, with m_i = gm_i / G
 * (Plummer-softened when soft > 0). The 1PN term is not included. */
double stc_space_energy(const StcSpace *s);

/* Total linear momentum (kg m/s) of massive bodies into p[3]. */
void stc_space_momentum(const StcSpace *s, double p[3]);

/* Number of live massive bodies / live tracers. */
int stc_space_live_massive(const StcSpace *s);
int stc_space_live_tracers(const StcSpace *s);

/* Osculating two-body elements of body i about body j (j = -1: about
 * its parent, or the most massive live body if none). Angles rad.
 * Returns 0, -1 if undefined (unbound -> a < 0, e >= 1 still filled). */
typedef struct {
    double a;       /* m (negative when hyperbolic) */
    double e;
    double inc;     /* to the J2000 ecliptic */
    double node;    /* longitude of ascending node */
    double argp;    /* argument of periapsis */
    double M;       /* mean anomaly */
    double period;  /* s (0 when unbound) */
    double r;       /* current distance, m */
    double v;       /* current relative speed, m/s */
} StcOrbit;
int stc_space_orbit(const StcSpace *s, int i, int j, StcOrbit *o);

#endif /* STC_SPACE_H */
