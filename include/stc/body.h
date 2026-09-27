#ifndef STC_BODY_H
#define STC_BODY_H

/* SPACETIME LAB — body: N-body sandbox physics core.
 * Velocity-Verlet gravity with Plummer softening, merge/bounce collisions,
 * Kepler elements, seeded presets. SI doubles, deterministic, no I/O.
 * Depends only on stc/core (vec/units/rng) + C stdlib.
 */

#include <stddef.h>
#include <stdint.h>

#include "stc/vec.h"

/* Maximum live bodies in one sim (fixed stack array, no malloc). */
#define STC_BODY_MAX 64

/* One gravitating body. `fixed` bodies exert gravity but never move;
 * `alive == 0` slots are skipped (merge/remove compaction). */
typedef struct {
    char name[32];
    double m;      /* mass, kg (> 0) */
    double r;      /* collision radius, m (>= 0) */
    StcVec3 p;     /* position, m */
    StcVec3 v;     /* velocity, m/s */
    int fixed;     /* nonzero: infinite inertia, never integrated */
    int alive;     /* nonzero: participates */
} StcBody;

/* Sandbox universe state. */
typedef struct {
    StcBody b[STC_BODY_MAX];
    int n;             /* live count (compacted, 0..STC_BODY_MAX) */
    double t;          /* sim time, s */
    double G;          /* gravitational constant, m^3/kg/s^2 */
    double soft2;      /* Plummer softening squared, m^2 */
    double restitution;/* bounce e in [0,1] when merge == 0 */
    int merge;         /* nonzero: collisions merge, else bounce */
    double vmax;       /* speed cap m/s after each step (0 = none) */
} StcSim;

/* Equation: sim init — zero state, store G/softening/restitution/merge. */
void stc_sim_init(StcSim *s, double G, double soft, double restitution,
                 int merge);

/* Equation: spawn — append a body; returns index or -1 (full/invalid).
 * Requires m > 0, r >= 0; name copied (truncated, NULL -> "body"). */
int stc_sim_add(StcSim *s, const char *name, double m, double r, StcVec3 p,
               StcVec3 v, int fixed);

/* Equation: swap-with-last removal; compacts n. Returns 0 / -1. */
int stc_sim_remove(StcSim *s, int idx);

/* Equation: velocity Verlet (kick-drift-kick) + collisions + speed cap.
 * Advances t by dt. Fixed bodies never move; dead slots skipped.
 * Equals the bare Verlet update followed by stc_sim_poststep. */
void stc_sim_step(StcSim *s, double dt);

/* Post-integration pass shared by all integrators.
 * Runs the merge/bounce collision pass and vmax cap exactly as
 * stc_sim_step does internally; stc_sim_step = verlet + poststep.
 * Callers pairing stc_sim_step_rk4 must call this after each step.
 * Advances no sim time: the owning stepper owns s->t. */
void stc_sim_poststep(StcSim *s);

/* Equation: Plummer-softened N-body acceleration of body idx. */
void stc_sim_accel(const StcSim *s, int idx, StcVec3 *a);

/* Equation: total energy — kinetic (free bodies) + pairwise potential. */
double stc_sim_energy(const StcSim *s);

/* Equation: vis-viva Kepler elements of idx about center (a, e, period;
 * period 0 when unbound/invalid). Any out-pointer may be NULL. */
void stc_sim_elements(const StcSim *s, int idx, int center, double *a,
                     double *e, double *period);

/* Seeded presets: 0 solar, 1 binary, 2 hole_disk, 3 slingshot. */
int stc_preset_count(void);
const char *stc_preset_name(int i);
int stc_preset_load(StcSim *s, int i, uint64_t seed);

/* Equation: gravitational redshift — proper rate sqrt(1-rs/r), clamped. */
double stc_body_proper_rate(double r, double M);

/* Equation: Lorentz factor 1/sqrt(1-(v/c)^2); returns -1 for v >= c. */
double stc_gamma(double v);

#endif /* STC_BODY_H */
