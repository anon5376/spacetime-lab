#ifndef STC_FRAG_H
#define STC_FRAG_H

/* SPACETIME LAB — frag: collisional fragmentation debris model.
 * Destroys two colliding bodies and spawns deterministic debris around
 * their center of mass. SI doubles, seeded RNG only, no I/O.
 * Depends only on stc/core (vec/units/rng) + stc/body.h + C stdlib.
 */

#include <stdint.h>

#include "stc/body.h"

/* Minimum debris spawned per impact (even at zero impact energy). */
#define STC_FRAG_MIN_DEBRIS 2

/* Impact-energy scale, J: debris count grows as 3*log10(1+E/E_REF). */
#define STC_FRAG_E_REF 1.0e20

/* Equation: center-of-mass impact energy + COM-frame zero-net-kick debris
 * dispersal (isotropic Gaussian velocity kicks, uniform-in-sphere
 * position offsets, volume-conserving radii).
 *
 * Collisional fragmentation of bodies a and b in s:
 *  - Parents are removed via stc_sim_remove (larger index first); their
 *    total mass M = m_a + m_b is partitioned over N debris with random
 *    weights in [0.5, 1.5), normalized to sum exactly to M.
 *  - N = STC_FRAG_MIN_DEBRIS + floor(3*log10(1+E_imp/STC_FRAG_E_REF)),
 *    where E_imp = 0.5*mu*|vrel|^2 is the center-of-mass impact energy
 *    (mu = reduced mass); N is at least 2 and is capped at the free
 *    slots left after removing the parents.
 *  - Debris velocities are v_cm + Gaussian kicks with RMS |kick| sigma =
 *    0.1*|vrel| + 0.1*v_esc (v_esc = mutual escape velocity at contact,
 *    sqrt(2*G*M/(r_a+r_b))), corrected so the mass-weighted mean kick is
 *    zero: total momentum equals the parents' M*v_cm to floating-point
 *    roundoff ("approximately": fixed parents contribute their stored
 *    velocity although fixed means infinite inertia, plus roundoff).
 *  - Debris positions are p_cm + uniform-in-sphere offsets of radius
 *    (r_a+r_b)/2, likewise mean-corrected to preserve the center of mass.
 *  - Debris radii conserve volume: r_i = cbrt((r_a^3+r_b^3) * m_i/M).
 *  - All debris are free (fixed = 0), named "frag0", "frag1", ...;
 *    fixed parents are consumed like any other (stored velocity used).
 *  - Deterministic: the seed drives one StcRng64 (splitmix64) stream in a
 *    fixed draw order, so the same sim + seed gives identical debris.
 *
 * s in/out: sim mutated via stc_sim_add/remove only (t, G, softening,
 *   collision mode untouched). a, b: parent indices in [0, s->n), a != b.
 * seed: determinism seed. No I/O.
 *
 * Returns the debris count (>= STC_FRAG_MIN_DEBRIS), or -1 on invalid
 * input (NULL sim, bad indices, non-positive/non-finite masses,
 * negative/non-finite radii or coordinates); -1 leaves s untouched.
 */
int stc_frag_impact(StcSim *s, int a, int b, uint64_t seed);

#endif /* STC_FRAG_H */
