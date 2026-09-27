#ifndef STC_VIS_MERGER_H
#define STC_VIS_MERGER_H

/* SPACETIME LAB — vis_merger: cinematic binary black hole merger view.
 *
 * Renders the stc/merger.h timeline as an HDR scene:
 *   - a spacetime sheet whose depth is the (softened) Newtonian potential
 *     of the holes and whose ripples are the retarded-time quadrupole
 *     strain r h+(t - r, azimuth) — the two-armed gravitational-wave
 *     spiral, exaggerated so it is visible;
 *   - black hole shadows as a screen-space metaball field, so the two
 *     shadows fuse through a peanut-shaped common horizon at merger;
 *   - thin-lens deflection of everything behind the holes (sky + far
 *     sheet), photon-ring glow, orbital trails;
 *   - bloom, ACES tone map, vignette; HUD with a live chirp trace.
 *
 * Claim class: ILLUSTRATIVE throughout (see stc/merger.h for the
 * ANALYTIC inspiral pieces). Sheet depth, ripple amplitude and falloff,
 * lens strength, and shadow sizes are visual choices, not a metric.
 *
 * Deterministic: frame-count driven, no wall clock; row passes run on a
 * worker pool but every pixel is computed independently, so the output
 * is identical for any thread count.
 */

#include "stc/app.h"
#include "stc/merger.h"
#include "stc/render.h"

#define STC_MVIS_NPRESET 4

typedef struct {
    StcBbh bbh;
    int preset;
    double msun;     /* total mass, solar masses (HUD unit readouts) */
    double t;        /* coordinate time, M */
    long frames;
    double speed;    /* playback multiplier */
    int paused;
    int grid_on;
    int lens_on;
    int hud_on;
    int trails_on;
    int autocam;
    double wave_gain;
    double yaw;      /* camera orbit, rad */
    double pitch;
    double dist;     /* camera distance, M */
    int w;
    int h;
    StcRgbBuf fb;    /* final 8-bit output */
    float *far;      /* w*h*3 HDR: behind the holes (lensed) */
    float *nearl;    /* w*h*3 HDR: in front (not lensed) */
    float *hdr;      /* w*h*3 composite */
    float *blm[7];   /* bloom pyramid */
    int blw[7];
    int blh[7];
    int nbl;
    char status[160];
} StcMVis;

/* Init preset 0..STC_MVIS_NPRESET-1 at w*h. Returns 0 ok, -1 bad preset
 * or size, -2 allocation failure (state freed). */
int stc_mvis_init(StcMVis *v, int preset, int w, int h);

/* Free buffers. NULL-safe. */
void stc_mvis_free(StcMVis *v);

/* Reallocate the framebuffers at w*h. 0 ok, -1 failure. */
int stc_mvis_resize(StcMVis *v, int w, int h);

/* Advance one frame (playback pacing slows through merger; the timeline
 * loops after the ringdown waves leave the sheet). */
void stc_mvis_step(StcMVis *v);

/* Paint the current state into v->fb. */
void stc_mvis_render(StcMVis *v);

/* Keys (stc/ui.h codes). Returns 1 to quit. */
int stc_mvis_key(StcMVis *v, int key);

/* Mouse orbit / zoom helpers for the window loop. */
void stc_mvis_orbit(StcMVis *v, double dyaw, double dpitch);
void stc_mvis_zoom(StcMVis *v, double factor);

/* Windowed session (macOS; elsewhere prints unsupported, returns 1). */
int stc_mvis_run(const StcAppOpts *o, int preset);

/* Headless: render nframes (frame index 0 = first) and write the last
 * frame as P6 PPM. When dir_every > 0 and dir != NULL, additionally
 * writes every dir_every-th frame to dir/merger_NNNNN.ppm (movie
 * frames). Returns 0 ok, -1 bad args, -2 alloc, -3 write failure. */
int stc_mvis_shot(int preset, long nframes, int w, int h,
                  const char *ppm_path, const char *dir, long dir_every);

#endif /* STC_VIS_MERGER_H */
