#ifndef STC_FMT_H
#define STC_FMT_H

/* SPACETIME LAB — fmt: human-readable SI formatting helpers.
 * Pure C11, no state, no sim mutation, no I/O (snprintf into caller buffer).
 */

#include <stddef.h>

/* Method: auto-unit time formatting (same bands as the sandbox HUD) —
 * |t| < 120s -> "%.1fs"; < 7200s -> "%.1fmin" (t/60);
 * < 172800s -> "%.2fhr" (t/3600); < 63072000s -> "%.2fday" (t/86400);
 * else -> "%.3fyr" (t/31557600, Julian year). Band selection uses |t|;
 * the sign of t is preserved in the printed value.
 * t: seconds, SI double. buf: caller-owned output buffer (may be NULL);
 * cap: buffer capacity in bytes. NULL buf or cap == 0 is a no-op. */
void stc_fmt_time(double t, char *buf, size_t cap);

/* Method: "%.3g %s" scientific formatting with NULL-unit safety —
 * unit != NULL and non-empty -> "<%.3g> <unit>"; otherwise just "<%.3g>".
 * x: value to format; unit: unit suffix (may be NULL); buf/cap as above.
 * NULL buf or cap == 0 is a no-op. */
void stc_fmt_sci(double x, const char *unit, char *buf, size_t cap);

#endif /* STC_FMT_H */
