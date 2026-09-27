#ifndef STC_VIS_SPACE_H
#define STC_VIS_SPACE_H

/* SPACETIME LAB — vis_space: the space simulator front end (--space).
 * JPL mission-operations style display of an stc/space.h simulation:
 * black field, thin grey osculating-orbit ellipses, monospace labels,
 * tick-marked ecliptic reference grid, selected-body readout panel
 * (distance, speed, osculating a/e/i/period), calendar date (TDB),
 * time warp both directions, exact rewind, and a command console for
 * sandbox edits (replace the Sun with a black hole, spawn N Earths,
 * jump to a date, ...). Claim classes are shown on screen: initial
 * conditions EXTERNAL (JPL Horizons), dynamics MEASURED (integrated),
 * sampled belts/swarms ILLUSTRATIVE.
 */

#include "stc/app.h"
#include "stc/render.h"
#include "stc/space.h"

#define STC_SVIS_NSCEN 6
#define STC_SVIS_CMD 96

typedef struct {
    StcSpace sim;
    int scenario;
    long frames;
    double warp;        /* simulated seconds per displayed second (signed) */
    int paused;
    int focus;          /* massive index the camera tracks, -1 = barycentre */
    int select;         /* massive index shown in the readout, -1 none */
    double yaw;         /* camera orbit, rad */
    double pitch;
    double dist;        /* camera distance, m */
    int show_orbits;
    int show_labels;
    int show_grid;
    int show_tracers;
    int show_hud;
    int console;        /* command console open */
    char cmd[STC_SVIS_CMD];
    int cmd_len;
    char msg[160];      /* last console/status message */
    long limited;       /* frames where warp was substep-limited */
    int w;
    int h;
    StcRgbBuf fb;
    float *acc;         /* w*h*3 linear accumulation buffer */
} StcSVis;

/* Init scenario 0..STC_SVIS_NSCEN-1 at w*h. 0 ok, -1 bad args, -2 alloc. */
int stc_svis_init(StcSVis *v, int scenario, int w, int h);
void stc_svis_free(StcSVis *v);
int stc_svis_resize(StcSVis *v, int w, int h);

/* Advance one displayed frame (1/60 s of warp). */
void stc_svis_step(StcSVis *v);
void stc_svis_render(StcSVis *v);

/* Keys (stc/ui.h codes). Returns 1 to quit. */
int stc_svis_key(StcSVis *v, int key);

/* Mouse: orbit by radians, zoom by factor, pick the body nearest the
 * framebuffer pixel (x, y) (origin top-left). */
void stc_svis_orbit(StcSVis *v, double dyaw, double dpitch);
void stc_svis_zoom(StcSVis *v, double factor);
void stc_svis_pick(StcSVis *v, int x, int y);

/* Scenario name for help text; NULL out of range. */
const char *stc_svis_scenario_name(int k);

/* Windowed session (macOS). */
int stc_svis_run(const StcAppOpts *o, int scenario);

/* Headless: like stc_mvis_shot (last frame to ppm_path and/or every
 * dir_every-th frame to dir/space_NNNNN.ppm). */
int stc_svis_shot(int scenario, long nframes, int w, int h,
                  const char *ppm_path, const char *dir, long dir_every);

#endif /* STC_VIS_SPACE_H */
