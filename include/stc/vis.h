#ifndef STC_VIS_H
#define STC_VIS_H

/* SPACETIME LAB — vis: graphical N-body sandbox frontend.
 * Same sim, same keys, same spawn form as the terminal sandbox, rendered
 * into an RGB framebuffer: shaded body discs, fading trail polylines,
 * starfield, and a bitmap-font HUD. Windowed display is macOS-native
 * (AppKit via the Objective-C runtime from pure C, system frameworks
 * only, no third-party dependencies); the rasterizer itself is portable
 * C11 and headless-testable through stc_vis_shot.
 */

#include "stc/app.h"
#include "stc/body.h"
#include "stc/render.h"
#include "stc/spawn.h"
#include "stc/trails.h"
#include "stc/ui.h"

/* Default window/viewport size, pixels. */
#define STC_VIS_DEF_W 1280
#define STC_VIS_DEF_H 800

/* Visual sandbox state: one window, one sim, one framebuffer. */
typedef struct {
    StcSim sim;
    StcUiCamera cam;
    StcTrailBank bank;
    int bank_ok;
    StcRgbBuf fb;
    SbxForm form;
    double dt;
    int warp_i;
    int follow;
    int paused;
    int forming;
    int help_on;
    int trail_on;
    int grid_on;   /* spacetime-grid sheet + axis aids (ILLUSTRATIVE) */
    int integrator;
    int preset;
    unsigned long long seed;
    long frames;
    const char *csv_path; /* borrowed export paths (NULL = defaults) */
    const char *ppm_path;
    char status[256];
} StcVis;

/* Method: visual state init — sim + preset + camera + trails + w*h
 * framebuffer. dt = 1000 s, warp x1, Verlet unless integrator says RK4.
 * Returns 0 on success, -1 on bad preset, -2 on framebuffer/trails
 * allocation failure (state left freed/empty). */
int stc_vis_init(StcVis *v, int preset, unsigned long long seed,
                 int integrator, int w, int h);

/* Method: visual state free — framebuffer + trails. NULL-safe. */
void stc_vis_free(StcVis *v);

/* Method: one frame — warp integrator substeps + trail push + frame
 * count. No-op while paused or while the spawn form is open. */
void stc_vis_step(StcVis *v);

/* Method: full repaint — background, trails, bodies, HUD, overlays —
 * into v->fb. NULL-safe no-op. */
void stc_vis_render(StcVis *v);

/* Method: key machine — full sandbox keymap (spawn/remove/preset/dt/
 * warp/trails/pause/export/help/quit/camera/form editing) on stc/ui.h
 * key codes. Returns 1 when the session should quit, 0 to continue. */
int stc_vis_key(StcVis *v, int key);

/* Method: framebuffer resize — reallocates v->fb at w*h (contents lost).
 * Returns 0 on success, -1 on failure (fb left empty). */
int stc_vis_resize(StcVis *v, int w, int h);

/* Method: windowed session — opens the native window and runs the
 * frame loop until quit (Q, window close, or --frames cap). Returns the
 * process exit code (0 ok, 1 runtime failure). Outside macOS prints
 * "visual unsupported" and returns 1. NULL-safe (uses defaults). */
int stc_vis_run(const StcAppOpts *o, int preset, unsigned long long seed);

/* Method: headless visual render — N frames of preset + seed at w*h,
 * then write the final framebuffer to ppm_path (P6). No window, no
 * pacing, deterministic: same args = byte-identical file. Returns 0 on
 * success, -1 bad args, -2 render failure, -3 PPM write failure. */
int stc_vis_shot(int preset, unsigned long long seed, int integrator,
                 long nframes, int w, int h, const char *ppm_path);

/* 5x7 glyph rows for ASCII c (row-major, bit4 = left pixel); unknown ->
 * space. Never NULL. */
const unsigned char *stc_vis_glyph(char c);

#endif /* STC_VIS_H */
