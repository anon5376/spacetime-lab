#include "stc/fmt.h"

#include <math.h>
#include <stdio.h>

/* Method: auto-unit time formatting (same bands as the sandbox HUD). */
void stc_fmt_time(double t, char *buf, size_t cap)
{
    double a;

    if (buf == NULL || cap == 0) {
        return;
    }
    a = fabs(t);
    if (a < 120.0) {
        snprintf(buf, cap, "%.1fs", t);
    } else if (a < 7200.0) {
        snprintf(buf, cap, "%.1fmin", t / 60.0);
    } else if (a < 172800.0) {
        snprintf(buf, cap, "%.2fhr", t / 3600.0);
    } else if (a < 63072000.0) {
        snprintf(buf, cap, "%.2fday", t / 86400.0);
    } else {
        snprintf(buf, cap, "%.3fyr", t / 31557600.0);
    }
}

/* Method: "%.3g %s" scientific formatting with NULL-unit safety. */
void stc_fmt_sci(double x, const char *unit, char *buf, size_t cap)
{
    if (buf == NULL || cap == 0) {
        return;
    }
    if (unit != NULL && unit[0] != '\0') {
        snprintf(buf, cap, "%.3g %s", x, unit);
    } else {
        snprintf(buf, cap, "%.3g", x);
    }
}
