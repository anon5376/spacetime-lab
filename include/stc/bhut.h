#ifndef STC_BHUT_H
#define STC_BHUT_H

/* SPACETIME LAB — bhut: Barnes-Hut octree approximate N-body gravity.
 * Method: Barnes-Hut hierarchical tree approximation (Barnes & Hut,
 * "A hierarchical O(N log N) force-calculation algorithm", Nature 324,
 * pp. 446-449, 1986). Far cells are approximated by their total mass at
 * the center of mass; near cells are opened recursively down to direct
 * Plummer-softened pairwise summation, matching stc_sim_accel.
 * SI doubles, deterministic (fixed insertion/traversal order), no I/O,
 * no globals. Memory: internal malloc/free per call (node pool + link
 * array), both freed before return; on allocation failure the call falls
 * back to exact direct O(N^2) summation.
 */

#include "stc/body.h"
#include "stc/vec.h"

/* Default opening-angle guidance: theta ~ 0.7. Smaller theta is more
 * accurate (theta <= 0 opens every cell: exact, direct-equivalent);
 * larger theta is faster and cruder. */
#define STC_BHUT_THETA_DEFAULT 0.7

/* Method: Barnes-Hut 1986 octree gravity with Plummer softening.
 * Approximates the acceleration of every body in s (fixed bodies exert
 * and receive values exactly like stc_sim_accel; dead slots get zero).
 * A cell of full size L is accepted as a monopole when
 * L < (theta/8) * (dist - rho), with dist the target-to-cell-center-
 * of-mass distance and rho the bounding radius of the cell's mass about
 * its center of mass; otherwise the cell is opened (conservative
 * distance-to-boundary variant of the Barnes-Hut opening criterion:
 * the documented 1/8 factor holds each accepted cell's quadrupole-level
 * error near 1% of its own contribution at the default theta, trading
 * pruning for accuracy in the N <= 64 regime; theta <= 0 opens every
 * cell and is exactly direct-equivalent). out must hold s->n entries;
 * out[i] is the acceleration of body i in m/s^2. Every return path
 * leaves out[0..s->n) fully written: the range is zero-filled up
 * front, so bodies that get no kick (dead slots, or a sim with no
 * live bodies) read as {0,0,0}, never uninitialized memory. NULL s
 * or NULL out is a no-op (nothing can be written); n <= 0 is a no-op
 * (the range is empty). */
void stc_bhut_accel(const StcSim *s, double theta, StcVec3 *out);

#endif /* STC_BHUT_H */
