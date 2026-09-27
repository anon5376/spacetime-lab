#ifndef STC_APP_H
#define STC_APP_H

/* SPACETIME LAB — app: CLI options + frame loop + deterministic export.
 * Pure C11 app layer over core + own (ui, csv) + scenario + render.
 * SI doubles on every boundary; seeded RNG only; same argv = byte-identical
 * CSV/PPM. Frame pacing: dt = 1/60 s unless --no-adapt is given.
 *
 * NOTE (integration): main.c drives the stc/scenario.h layer directly
 * (setup/step per frame, looping the preset when --frames exceeds it) and
 * renders through stc/render.h (ASCII live view, PPM export). No CLI change.
 */

#include "stc/vec.h"
#include "stc/units.h"
#include "stc/rng.h"
#include "stc/ui.h"
#include "stc/csv.h"

/* Equation: frame pacing — fixed step dt = 1/60 s (60 fps). */
#define STC_APP_FRAME_DT (1.0 / 60.0)

/* Default frames per scenario for batch runs (10 s at 60 fps). */
#define STC_APP_DEFAULT_FRAMES 600L

/* Frames per scenario for --smoke without an explicit --frames. */
#define STC_APP_SMOKE_FRAMES 1L

/* Upper bound accepted by --frames (determinism cap, documented). */
#define STC_APP_MAX_FRAMES 1000000L

/* Contract Sec 5 preset count (twin..shadow). */
#define STC_APP_NSCEN 10

/* PPM export raster size (P6 binary via stc_render_write_ppm). */
#define STC_APP_PPM_W 320
#define STC_APP_PPM_H 200

/* Interactive ASCII view size (columns x rows). */
#define STC_APP_ASCII_W 76
#define STC_APP_ASCII_H 22

/* Sandbox stepper select (--integrator verlet|rk4). Verlet is the default:
 * velocity Verlet; RK4 is the classic 4th-order Runge-Kutta step (see
 * stc/rk4.h). The app runs stc_sim_poststep (collisions + vmax cap) after
 * every step under either choice. Scenarios ignore it. */
#define STC_APP_INTEG_VERLET 0
#define STC_APP_INTEG_RK4 1

/* Parsed CLI state (contract Sec 6 flags plus --no-adapt, sandbox flags). */
typedef struct {
    const char *scenario; /* preset name or "all"; default "twin" */
    int scenario_given;   /* nonzero iff --scenario was seen */
    long frames;          /* frames per scenario; default 600 */
    int frames_given;     /* nonzero iff --frames was seen */
    const char *ppm;      /* --ppm path (file, or dir when scenario=all) */
    const char *csv;      /* --csv path (file, or dir when scenario=all) */
    int smoke;            /* headless run, exit 0 */
    int list;             /* print preset names and exit 0 */
    int help;             /* print help and exit 0 */
    int no_adapt;         /* disable 60fps pacing */
    int sandbox;          /* N-body sandbox mode instead of scenarios */
    long preset;          /* sandbox preset index; default 0 */
    int preset_given;     /* nonzero iff --preset was seen */
    unsigned long long seed; /* sandbox RNG seed; default 20260927 */
    int seed_given;       /* nonzero iff --seed was seen */
    int integrator;       /* sandbox stepper; default STC_APP_INTEG_VERLET */
    int integrator_given; /* nonzero iff --integrator was seen */
    int visual;           /* native-window visual sandbox (macOS) */
    const char *visual_shot; /* --visual-shot path (headless render) */
    int merger;           /* cinematic binary black hole merger view */
    const char *movie;    /* --movie dir (headless merger/space frames) */
    long movie_every;     /* --movie-every N; default 2 */
    int space;            /* full space simulator (real Solar System) */
} StcAppOpts;

/* Equation: CLI defaults — scenario "twin", 600 frames, no outputs. NULL-safe. */
void stc_app_opts_init(StcAppOpts *o);

/* Equation: argv parse — exact contract flags, --flag value and --flag=value.
 * Returns 0 on success (then call stc_app_run), 2 on usage error after
 * printing diagnostics to stderr. NULL o/argv or argc <= 0 -> 2. */
int stc_app_parse_args(StcAppOpts *o, int argc, char *argv[]);

/* Equation: help text — documents every flag (contract Sec 6). */
void stc_app_print_help(const char *prog);

/* Equation: scenario list — prints the 10 contract Sec 5 preset names. */
void stc_app_print_list(void);

/* Equation: frame loop — N frames at dt = 1/60 s per selected scenario.
 * Dispatches help/list/smoke/batch/interactive; interactive iff stdin is a
 * TTY and --smoke is absent (else headless batch). Returns the process
 * exit code: 0 ok, 1 runtime failure (e.g. output unwritable), 2 bad options. */
int stc_app_run(const StcAppOpts *o);

/* Equation: preset table — contract Sec 5 count (always 10). */
int stc_app_scenario_count(void);

/* Equation: preset table — name of index 0..9, NULL when out of range. */
const char *stc_app_scenario_name(int idx);

/* Equation: preset lookup — index of name, -2 for "all", -1 if unknown/NULL. */
int stc_app_scenario_index(const char *name);

#endif /* STC_APP_H */
