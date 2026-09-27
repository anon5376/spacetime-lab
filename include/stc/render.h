#ifndef STC_RENDER_H
#define STC_RENDER_H

/* SPACETIME LAB — render: pinhole camera + ASCII/ANSI canvas + RGB rasterizer.
 *
 * Pure C11, zero dependencies (C stdlib + libm only). Camera geometry uses
 * SI doubles (meters, radians); pixels are integers with origin at top-left
 * (x right, y down). All functions are deterministic: fixed loop bounds, no
 * RNG, no wall clock, no global state. Terminal output uses ANSI SGR color
 * escapes only (no cursor motion, no screen clearing).
 *
 * Equation index (each public function also names its equation):
 *   camera   = lookAt view x perspective projection (core StcMat4)
 *   project  = clip = VP * world, NDC = clip / w, viewport map
 *   line     = Bresenham integer line
 *   triangle = edge-function flat fill
 *   shade    = flat Lambert factor: C = base * clamp(shade, 0, 1)
 *   ppm      = PPM P6 binary raster
 *   ansi     = ANSI SGR color escapes (3-bit + bright)
 *   palette  = phosphor-green / amber linear ramps
 */

#include "stc/vec.h"

#include <stdint.h>
#include <stdio.h>

/* ---------- ANSI color indices (equation: ANSI SGR 3-bit + bright) ---------- */
/* Values 0..7 map to SGR 30..37; 8..15 map to SGR 90..97. Any other value
 * passed as a color means "default terminal color" (no escape emitted). */
typedef enum {
    STC_RENDER_BLACK = 0,
    STC_RENDER_RED = 1,
    STC_RENDER_GREEN = 2,
    STC_RENDER_YELLOW = 3,
    STC_RENDER_BLUE = 4,
    STC_RENDER_MAGENTA = 5,
    STC_RENDER_CYAN = 6,
    STC_RENDER_WHITE = 7,
    STC_RENDER_BRIGHT_BLACK = 8,
    STC_RENDER_BRIGHT_RED = 9,
    STC_RENDER_BRIGHT_GREEN = 10,
    STC_RENDER_BRIGHT_YELLOW = 11,
    STC_RENDER_BRIGHT_BLUE = 12,
    STC_RENDER_BRIGHT_MAGENTA = 13,
    STC_RENDER_BRIGHT_CYAN = 14,
    STC_RENDER_BRIGHT_WHITE = 15
} StcRenderColor;

/* Out-of-range alias meaning "default terminal color". */
#define STC_RENDER_DEFAULT (-1)

/* ---------- Camera (equation: VP = perspective * lookAt, core StcMat4) ---------- */
typedef struct {
    StcVec3 eye;      /* camera position, meters (world) */
    StcVec3 center;   /* look-at target, meters (world) */
    StcVec3 up;       /* view-up direction (need not be unit) */
    double fovy_rad;  /* vertical field of view, radians */
    double aspect;    /* viewport width / height, dimensionless */
    double znear;     /* near plane distance, meters */
    double zfar;      /* far plane distance, meters */
    int vp_w;         /* viewport width, pixels */
    int vp_h;         /* viewport height, pixels */
    StcMat4 view;     /* lookAt(eye, center, up) */
    StcMat4 proj;     /* perspective(fovy, aspect, znear, zfar) */
    StcMat4 viewproj; /* proj * view */
} StcCamera;

/* Build a pinhole camera. Equation: view = lookAt, proj = perspective,
 * viewproj = proj * view (all core StcMat4). Degenerate inputs fall back to
 * core semantics (identity matrices); aspect <= 0 (or NaN) is derived from
 * vp_w / vp_h when valid, else 1.0; non-positive vp dims become 0. */
StcCamera stc_render_camera(StcVec3 eye, StcVec3 center, StcVec3 up,
                           double fovy_rad, double aspect,
                           double znear, double zfar,
                           int vp_w, int vp_h);

/* Project a world point to viewport pixels. Equation: clip = VP * [world, 1],
 * ndc = clip.xyz / clip.w, sx = (ndc.x/2+1/2)*w, sy = (1-(ndc.y/2+1/2))*h.
 * world in meters; sx/sy in pixels (origin top-left, unclamped); depth_ndc
 * dimensionless in [-1, 1] on success. Any out-pointer may be NULL.
 * Returns 1 if in front of camera (clip w > 0) and inside NDC + viewport,
 * else 0 (outputs still written when non-NULL: 0/0/0.0 on NULL cam). */
int stc_render_project_point(const StcCamera *cam, StcVec3 world,
                             int *sx, int *sy, double *depth_ndc);

/* ---------- ASCII canvas (row-major, origin top-left) ---------- */
typedef struct {
    int w;        /* width, pixels (columns) */
    int h;        /* height, pixels (rows) */
    char *ch;     /* w*h glyphs, (x,y) -> y*w+x */
    uint8_t *fg;  /* w*h ANSI colors 0..15, 255 = default */
} StcAsciiBuf;

/* Allocate a w*h ASCII buffer (documented malloc/free pair with
 * stc_render_ascii_free). Returns 0 on success, -1 on NULL buf, non-positive
 * dims, size overflow, or allocation failure (buf zeroed on failure). */
int stc_render_ascii_init(StcAsciiBuf *b, int w, int h);

/* Release an ASCII buffer; NULL buf or already-freed buf is a no-op. */
void stc_render_ascii_free(StcAsciiBuf *b);

/* Fill every cell. Equation: constant clear. Out-of-range color -> default. */
void stc_render_ascii_clear(StcAsciiBuf *b, char fill, int color);

/* Set one cell; out-of-bounds plots are silently clipped (no-op). */
void stc_render_ascii_plot(StcAsciiBuf *b, int x, int y, char ch, int color);

/* Draw a glyph run. Equation: Bresenham integer line, bounded by
 * |dx| + |dy| + 2 plotted points; endpoints clipped by plot. */
void stc_render_ascii_line(StcAsciiBuf *b, int x0, int y0, int x1, int y1,
                           char ch, int color);

/* Write the buffer to a stream. Equation: ANSI SGR color runs + '\n' per row;
 * only SGR escapes ("ESC[<n>m", "ESC[0m") are emitted, no cursor control.
 * Non-printable cells render as spaces. Returns 0 on success, -1 on NULL
 * args, empty/uninitialized buffer, or stream error. */
int stc_render_ascii_present(const StcAsciiBuf *b, FILE *out);

/* ---------- RGB rasterizer (row-major, origin top-left, row 0 = top) ---------- */
typedef struct {
    uint8_t r;
    uint8_t g;
    uint8_t b;
} StcRgb;

typedef struct {
    int w;      /* width, pixels */
    int h;      /* height, pixels */
    StcRgb *px; /* w*h pixels, (x,y) -> y*w+x */
} StcRgbBuf;

/* Pack three channels. Equation: RGB identity constructor. */
StcRgb stc_render_rgb(uint8_t r, uint8_t g, uint8_t b);

/* Allocate a w*h RGB buffer (documented malloc/free pair with
 * stc_render_rgb_free). Returns 0 on success, -1 on NULL buf, non-positive
 * dims, size overflow, or allocation failure (buf zeroed on failure). */
int stc_render_rgb_init(StcRgbBuf *b, int w, int h);

/* Release an RGB buffer; NULL buf or already-freed buf is a no-op. */
void stc_render_rgb_free(StcRgbBuf *b);

/* Fill every pixel. Equation: constant clear. */
void stc_render_rgb_clear(StcRgbBuf *b, StcRgb c);

/* Set one pixel; out-of-bounds plots are silently clipped (no-op). */
void stc_render_rgb_plot(StcRgbBuf *b, int x, int y, StcRgb c);

/* Draw a pixel run. Equation: Bresenham integer line, bounded by
 * |dx| + |dy| + 2 plotted points; endpoints clipped by plot. */
void stc_render_rgb_line(StcRgbBuf *b, int x0, int y0, int x1, int y1,
                         StcRgb c);

/* Fill a triangle with one flat color. Equation: edge-function half-space
 * fill (sign-consistent winding, boundary-inclusive). Degenerate
 * (zero-area) triangles fall back to stroking the three edges. */
void stc_render_rgb_triangle(StcRgbBuf *b,
                             int x0, int y0, int x1, int y1, int x2, int y2,
                             StcRgb c);

/* Fill a triangle with a uniformly shaded color. Equation: flat shading
 * C = base * clamp(shade, 0, 1) (NaN shade -> 0), then edge-function fill. */
void stc_render_rgb_triangle_shaded(StcRgbBuf *b,
                                    int x0, int y0,
                                    int x1, int y1,
                                    int x2, int y2,
                                    StcRgb base, double shade);

/* Write the buffer to a file. Equation: PPM P6 binary raster
 * ("P6\\n<w> <h>\\n255\\n" + w*h raw RGB triplets, top row first).
 * Returns 0 on success, -1 on NULL args, empty buffer, or I/O failure. */
int stc_render_write_ppm(const StcRgbBuf *b, const char *path);

/* ---------- ANSI + palette helpers ---------- */

/* SGR escape for a color index. Equation: ANSI SGR table lookup
 * (0..7 -> "ESC[3<k>m", 8..15 -> "ESC[9<k>m"); out-of-range -> "". */
const char *stc_render_ansi_escape(int color);

/* Phosphor-green ramp black -> (51, 255, 51). Equation: C = full * clamp01(t)
 * with round-half-up; NaN t -> black. t dimensionless in [0, 1]. */
StcRgb stc_render_phosphor(double t);

/* Amber ramp black -> (255, 176, 0). Equation: C = full * clamp01(t) with
 * round-half-up; NaN t -> black. t dimensionless in [0, 1]. */
StcRgb stc_render_amber(double t);

#endif /* STC_RENDER_H */
