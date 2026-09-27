/* SPACETIME LAB — help: static sandbox key guide + splash title block.
 * Pure strings, no I/O, no sim mutation, no RNG/clock (deterministic).
 */

#include "stc/help.h"

/* Method: static string table — sandbox key guide (see header). */
const char *stc_help_sandbox_keys(void)
{
    static const char kHelp[] =
        "sandbox keys:\n"
        "  W/S or up/down     orbit camera pitch up/down\n"
        "  A or left/right    orbit camera yaw (D = grid in visual)\n"
        "  +/= or R, -/_ or F zoom camera in/out\n"
        "  tab                cycle follow: centroid -> body 0..n-1\n"
        "  N                  spawn body form (enter on SPAWN> spawns)\n"
        "  K                  remove followed body\n"
        "  G                  next preset (solar/binary/hole_disk/slingshot)\n"
        "  Z/X                halve/double dt (1e-6..1e9 s)\n"
        "  C                  cycle warp x1/x10/x100 steps per frame\n"
        "  T                  toggle trails\n"
        "  D                  grid sheet + axis aids (visual only;\n"
        "                     illustrative sheet, not the metric)\n"
        "  P                  pause/resume\n"
        "  E                  export bodies CSV + PPM view\n"
        "  ?                  toggle this help\n"
        "  Q or esc or ctrl-c quit\n"
        "spawn form: up/down field, type to edit, backspace deletes,\n"
        "  enter next field, enter on SPAWN> spawns, esc cancels\n";

    return kHelp;
}

/* Method: static string table — splash title block (see header). */
const char *stc_splash(void)
{
    static const char kSplash[] =
        "+================================+\n"
        "|          SPACETIME LAB         |\n"
        "|    N-body gravity sandbox      |\n"
        "|  SI units * seeded * 60 fps    |\n"
        "|     press ? for key help       |\n"
        "+================================+\n";

    return kSplash;
}
