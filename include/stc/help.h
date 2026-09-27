#ifndef STC_HELP_H
#define STC_HELP_H

/* SPACETIME LAB — help: static sandbox key guide + splash title block.
 * Pure C11, no state, no sim mutation, no I/O (caller prints).
 * SI doubles N/A (text only); deterministic (fixed string literals).
 * The integrator wires '?' to toggle the key guide and shows a one-line
 * startup hint pointing at it.
 */

/* Method: static string table — sandbox key guide.
 * Returns a static multi-line NUL-terminated guide covering the
 * interactive sandbox keys (camera, follow, spawn/remove, preset, dt,
 * warp, trails, pause, export, help, quit, plus spawn-form editing).
 * Never returns NULL; the pointer stays valid for the program lifetime. */
const char *stc_help_sandbox_keys(void);

/* Method: static string table — splash title block.
 * Returns a static multi-line NUL-terminated ASCII title block for
 * SPACETIME LAB / sandbox startup. Never returns NULL; the pointer
 * stays valid for the program lifetime. */
const char *stc_splash(void);

#endif /* STC_HELP_H */
