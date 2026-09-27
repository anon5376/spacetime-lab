#ifndef STC_RK4_H
#define STC_RK4_H

/* SPACETIME LAB — rk4: classic 4th-order Runge-Kutta N-body stepper.
 * Same StcSim state as the velocity-Verlet core (stc/body.h). Pure C11,
 * SI doubles, deterministic, no I/O, no allocation, no globals.
 * Every function below names its equation.
 */

#include "stc/body.h"

/* Equation: classic 4th-order Runge-Kutta (RK4) step of the full N-body
 * state y = (p, v) with dy/dt = (v, a(p)) —
 *   k1 = f(y0), k2 = f(y0 + dt/2 k1), k3 = f(y0 + dt/2 k2),
 *   k4 = f(y0 + dt k3),
 *   y1 = y0 + dt/6 (k1 + 2 k2 + 2 k3 + k4),
 * where the acceleration a(p) at each stage is the same Plummer-softened
 * N-body gravity used by the Verlet core (queried via stc_sim_accel on
 * the stage-shifted positions, so the force model matches exactly).
 * s: sim state (must be non-NULL; NULL is a no-op); dt: step size, s
 * (added to s->t unconditionally, even when n == 0).
 * Respects fixed/alive: fixed bodies exert gravity but are never moved,
 * dead (alive == 0) slots are skipped. Out-of-range n is clamped to
 * [0, STC_BODY_MAX] for the step.
 * NO collisions inside: unlike stc_sim_step, this performs no merge /
 * bounce pass and no vmax cap — the caller pairs it with a collision
 * pass (and any speed cap) when the scenario needs one.
 * Deterministic: fixed evaluation order, no RNG, no allocation.
 * Cite: classic 4th-order Runge-Kutta (Kutta, 1901). */
void stc_sim_step_rk4(StcSim *s, double dt);

#endif /* STC_RK4_H */
