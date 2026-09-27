#ifndef STC_SANDBOX_H
#define STC_SANDBOX_H

/* SPACETIME LAB — sandbox: interactive N-body gravity playground.
 * Spawns bodies with arbitrary physical characteristics, integrates with
 * the stc/body.h core (velocity Verlet), renders a character-grid viewport,
 * exports CSV/PPM. Pure C11 over core + body + render + ui.
 */

#include "stc/app.h"

/* Equation: sandbox run — interactive TTY screen, or headless batch when
 * stdin is not a TTY or --smoke is given. Returns the process exit code
 * (0 ok, 1 runtime failure, 2 bad options). NULL-safe (uses defaults). */
int stc_sandbox_run(const StcAppOpts *o);

#endif /* STC_SANDBOX_H */
