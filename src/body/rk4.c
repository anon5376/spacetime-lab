/* SPACETIME LAB — rk4: classic 4th-order Runge-Kutta N-body stepper.
 * Integrates the full state y = (p, v) with dy/dt = (v, a(p)); the stage
 * accelerations are the Verlet core's Plummer-softened N-body gravity via
 * stc_sim_accel. No collisions, no speed cap (caller-owned post-passes).
 */

#include "stc/rk4.h"

#include "stc/body.h"
#include "stc/vec.h"

/* Method: integrated-body predicate — alive and free to move. */
static int rk4_moves(const StcSim *s, int i)
{
    return s->b[i].alive && !s->b[i].fixed;
}

void stc_sim_step_rk4(StcSim *s, double dt)
{
    /* Method: classic 4th-order Runge-Kutta (Kutta, 1901) on the coupled
     * first-order system dy/dt = (v, a(p)) for all free bodies at once.
     * Four acceleration evaluations per step (one per stage); the k_p
     * slopes are rebuilt from the stored velocities and accelerations. */
    static const StcVec3 ZERO = { 0.0, 0.0, 0.0 };
    StcVec3 p0[STC_BODY_MAX];
    StcVec3 v0[STC_BODY_MAX];
    StcVec3 a1[STC_BODY_MAX];
    StcVec3 a2[STC_BODY_MAX];
    StcVec3 a3[STC_BODY_MAX];
    StcVec3 a4[STC_BODY_MAX];
    StcSim tmp;
    int n;
    int i;

    if (!s) {
        return;
    }
    n = s->n;
    if (n < 0) {
        n = 0;
    }
    if (n > STC_BODY_MAX) {
        n = STC_BODY_MAX;
    }
    if (n == 0) {
        s->t += dt;
        return;
    }

    for (i = 0; i < n; i++) {
        p0[i] = s->b[i].p;
        v0[i] = s->b[i].v;
        a1[i] = ZERO;
        a2[i] = ZERO;
        a3[i] = ZERO;
        a4[i] = ZERO;
    }

    /* Stage 1: k1p = v0, k1v = a(p0). */
    for (i = 0; i < n; i++) {
        if (rk4_moves(s, i)) {
            stc_sim_accel(s, i, &a1[i]);
        }
    }

    /* Stage 2: k2v = a(p0 + dt/2 k1p). */
    tmp = *s;
    for (i = 0; i < n; i++) {
        if (rk4_moves(s, i)) {
            tmp.b[i].p =
                stc_vec3_add(p0[i], stc_vec3_scale(v0[i], 0.5 * dt));
        }
    }
    for (i = 0; i < n; i++) {
        if (rk4_moves(s, i)) {
            stc_sim_accel(&tmp, i, &a2[i]);
        }
    }

    /* Stage 3: k3v = a(p0 + dt/2 k2p), k2p = v0 + dt/2 k1v. */
    tmp = *s;
    for (i = 0; i < n; i++) {
        if (rk4_moves(s, i)) {
            StcVec3 k2p =
                stc_vec3_add(v0[i], stc_vec3_scale(a1[i], 0.5 * dt));
            tmp.b[i].p =
                stc_vec3_add(p0[i], stc_vec3_scale(k2p, 0.5 * dt));
        }
    }
    for (i = 0; i < n; i++) {
        if (rk4_moves(s, i)) {
            stc_sim_accel(&tmp, i, &a3[i]);
        }
    }

    /* Stage 4: k4v = a(p0 + dt k3p), k3p = v0 + dt/2 k2v. */
    tmp = *s;
    for (i = 0; i < n; i++) {
        if (rk4_moves(s, i)) {
            StcVec3 k3p =
                stc_vec3_add(v0[i], stc_vec3_scale(a2[i], 0.5 * dt));
            tmp.b[i].p = stc_vec3_add(p0[i], stc_vec3_scale(k3p, dt));
        }
    }
    for (i = 0; i < n; i++) {
        if (rk4_moves(s, i)) {
            stc_sim_accel(&tmp, i, &a4[i]);
        }
    }

    /* Combine: y1 = y0 + dt/6 (k1 + 2 k2 + 2 k3 + k4), both halves. */
    for (i = 0; i < n; i++) {
        if (rk4_moves(s, i)) {
            StcVec3 k2p =
                stc_vec3_add(v0[i], stc_vec3_scale(a1[i], 0.5 * dt));
            StcVec3 k3p =
                stc_vec3_add(v0[i], stc_vec3_scale(a2[i], 0.5 * dt));
            StcVec3 k4p = stc_vec3_add(v0[i], stc_vec3_scale(a3[i], dt));
            StcVec3 dp = stc_vec3_add(
                stc_vec3_add(v0[i], stc_vec3_scale(k2p, 2.0)),
                stc_vec3_add(stc_vec3_scale(k3p, 2.0), k4p));
            StcVec3 dv = stc_vec3_add(
                stc_vec3_add(a1[i], stc_vec3_scale(a2[i], 2.0)),
                stc_vec3_add(stc_vec3_scale(a3[i], 2.0), a4[i]));
            s->b[i].p = stc_vec3_add(p0[i], stc_vec3_scale(dp, dt / 6.0));
            s->b[i].v = stc_vec3_add(v0[i], stc_vec3_scale(dv, dt / 6.0));
        }
    }

    s->t += dt;
}
