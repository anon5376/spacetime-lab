/* SPACETIME LAB — render implementation: camera + ASCII/ANSI + RGB + PPM.
 * Pure C11 (stdlib + libm only). Deterministic: fixed loop bounds, no RNG,
 * no wall clock, no global state. Includes only own + stc/core + stdlib.
 */

#include "stc/render.h"
#include "stc/vec.h"

#include <math.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* PPM writer blits pixels as raw triplets; StcRgb must have no padding. */
_Static_assert(sizeof(StcRgb) == 3, "StcRgb must be 3 packed bytes");

/* ---------- internal helpers ---------- */

static int render_color_valid(int color)
{
    return color >= 0 && color <= 15;
}

static uint8_t render_color_norm(int color)
{
    return render_color_valid(color) ? (uint8_t)color : 255u;
}

/* NaN-safe clamp to [0, 1]: NaN -> 0. */
static double render_clamp01(double t)
{
    if (!(t >= 0.0)) {
        return 0.0;
    }
    if (t > 1.0) {
        return 1.0;
    }
    return t;
}

static uint8_t render_ramp_channel(double full, double t)
{
    double v = full * render_clamp01(t) + 0.5;

    if (v <= 0.0) {
        return 0;
    }
    if (v >= 255.0) {
        return 255;
    }
    return (uint8_t)v;
}

/* Multiply-accumulate-free size check: 0 on success with *n = w*h. */
static int render_pixel_count(int w, int h, size_t *n)
{
    if (w <= 0 || h <= 0) {
        return -1;
    }
    *n = (size_t)w * (size_t)h;
    if ((size_t)w != 0 && *n / (size_t)w != (size_t)h) {
        return -1;
    }
    return 0;
}

/* ---------- camera ---------- */

StcCamera stc_render_camera(StcVec3 eye, StcVec3 center, StcVec3 up,
                           double fovy_rad, double aspect,
                           double znear, double zfar,
                           int vp_w, int vp_h)
{
    StcCamera cam;

    cam.eye = eye;
    cam.center = center;
    cam.up = up;
    cam.fovy_rad = fovy_rad;
    if (!(aspect > 0.0)) {
        if (vp_w > 0 && vp_h > 0) {
            aspect = (double)vp_w / (double)vp_h;
        } else {
            aspect = 1.0;
        }
    }
    cam.aspect = aspect;
    cam.znear = znear;
    cam.zfar = zfar;
    cam.vp_w = vp_w > 0 ? vp_w : 0;
    cam.vp_h = vp_h > 0 ? vp_h : 0;
    /* Equation: view = lookAt, proj = perspective, viewproj = proj * view. */
    cam.view = stc_mat4_lookat(eye, center, up);
    cam.proj = stc_mat4_perspective(fovy_rad, aspect, znear, zfar);
    cam.viewproj = stc_mat4_mul(cam.proj, cam.view);
    return cam;
}

int stc_render_project_point(const StcCamera *cam, StcVec3 world,
                             int *sx, int *sy, double *depth_ndc)
{
    const double *m;
    double cx, cy, cz, cw;
    double nx, ny, nz, fx, fy;
    int px, py;

    if (sx != NULL) {
        *sx = 0;
    }
    if (sy != NULL) {
        *sy = 0;
    }
    if (depth_ndc != NULL) {
        *depth_ndc = 0.0;
    }
    if (cam == NULL || cam->vp_w <= 0 || cam->vp_h <= 0) {
        return 0;
    }
    /* Equation: clip = VP * [world, 1] (column-major, column vector). */
    m = cam->viewproj.m;
    cx = m[0] * world.x + m[4] * world.y + m[8] * world.z + m[12];
    cy = m[1] * world.x + m[5] * world.y + m[9] * world.z + m[13];
    cz = m[2] * world.x + m[6] * world.y + m[10] * world.z + m[14];
    cw = m[3] * world.x + m[7] * world.y + m[11] * world.z + m[15];
    if (!(cw > 0.0)) {
        return 0; /* behind camera, degenerate, or NaN */
    }
    /* Equation: ndc = clip.xyz / clip.w, then viewport map (y down). */
    nx = cx / cw;
    ny = cy / cw;
    nz = cz / cw;
    fx = (nx * 0.5 + 0.5) * (double)cam->vp_w;
    fy = (1.0 - (ny * 0.5 + 0.5)) * (double)cam->vp_h;
    px = (int)floor(fx);
    py = (int)floor(fy);
    if (sx != NULL) {
        *sx = px;
    }
    if (sy != NULL) {
        *sy = py;
    }
    if (depth_ndc != NULL) {
        *depth_ndc = nz;
    }
    if (nx < -1.0 || nx > 1.0 || ny < -1.0 || ny > 1.0 || nz < -1.0 ||
        nz > 1.0) {
        return 0;
    }
    if (px < 0 || px >= cam->vp_w || py < 0 || py >= cam->vp_h) {
        return 0;
    }
    return 1;
}

/* ---------- ASCII canvas ---------- */

int stc_render_ascii_init(StcAsciiBuf *b, int w, int h)
{
    size_t n = 0;

    if (b == NULL) {
        return -1;
    }
    b->w = 0;
    b->h = 0;
    b->ch = NULL;
    b->fg = NULL;
    if (render_pixel_count(w, h, &n) != 0) {
        return -1;
    }
    b->ch = (char *)malloc(n * sizeof(char));
    b->fg = (uint8_t *)malloc(n * sizeof(uint8_t));
    if (b->ch == NULL || b->fg == NULL) {
        free(b->ch);
        free(b->fg);
        b->ch = NULL;
        b->fg = NULL;
        return -1;
    }
    b->w = w;
    b->h = h;
    memset(b->ch, ' ', n);
    memset(b->fg, 255, n);
    return 0;
}

void stc_render_ascii_free(StcAsciiBuf *b)
{
    if (b == NULL) {
        return;
    }
    free(b->ch);
    free(b->fg);
    b->w = 0;
    b->h = 0;
    b->ch = NULL;
    b->fg = NULL;
}

static int render_ascii_ready(const StcAsciiBuf *b)
{
    return b != NULL && b->w > 0 && b->h > 0 && b->ch != NULL && b->fg != NULL;
}

void stc_render_ascii_clear(StcAsciiBuf *b, char fill, int color)
{
    size_t n;

    if (!render_ascii_ready(b)) {
        return;
    }
    n = (size_t)b->w * (size_t)b->h;
    memset(b->ch, fill, n);
    memset(b->fg, render_color_norm(color), n);
}

void stc_render_ascii_plot(StcAsciiBuf *b, int x, int y, char ch, int color)
{
    size_t i;

    if (!render_ascii_ready(b)) {
        return;
    }
    if (x < 0 || x >= b->w || y < 0 || y >= b->h) {
        return;
    }
    i = (size_t)y * (size_t)b->w + (size_t)x;
    b->ch[i] = ch;
    b->fg[i] = render_color_norm(color);
}

void stc_render_ascii_line(StcAsciiBuf *b, int x0, int y0, int x1, int y1,
                           char ch, int color)
{
    /* Equation: Bresenham integer line in 64-bit deltas (no int overflow). */
    long long dx = llabs((long long)x1 - (long long)x0);
    long long dy = -llabs((long long)y1 - (long long)y0);
    long long err = dx + dy;
    long long cap = dx - dy + 2; /* max plotted points: |dx|+|dy|+2 */
    long long i;
    int sx = x0 < x1 ? 1 : -1;
    int sy = y0 < y1 ? 1 : -1;
    int x = x0;
    int y = y0;

    if (!render_ascii_ready(b)) {
        return;
    }
    for (i = 0; i < cap; i++) {
        long long e2;

        stc_render_ascii_plot(b, x, y, ch, color);
        if (x == x1 && y == y1) {
            break;
        }
        e2 = 2 * err;
        if (e2 >= dy) {
            err += dy;
            x += sx;
        }
        if (e2 <= dx) {
            err += dx;
            y += sy;
        }
    }
}

int stc_render_ascii_present(const StcAsciiBuf *b, FILE *out)
{
    int y, x;
    unsigned cur;

    if (!render_ascii_ready(b) || out == NULL) {
        return -1;
    }
    for (y = 0; y < b->h; y++) {
        cur = 255u; /* default: no escape in effect */
        for (x = 0; x < b->w; x++) {
            size_t i = (size_t)y * (size_t)b->w + (size_t)x;
            unsigned want = b->fg[i];
            char c = b->ch[i];

            if (want != cur) {
                /* Equation: ANSI SGR select on color change only. */
                fputs(stc_render_ansi_escape((int)want), out);
                cur = want;
            }
            if (c < 32 || c > 126) {
                c = ' ';
            }
            fputc(c, out);
        }
        if (cur != 255u) {
            fputs("\x1b[0m", out);
        }
        fputc('\n', out);
    }
    return ferror(out) ? -1 : 0;
}

/* ---------- RGB rasterizer ---------- */

StcRgb stc_render_rgb(uint8_t r, uint8_t g, uint8_t b)
{
    StcRgb c;

    c.r = r;
    c.g = g;
    c.b = b;
    return c;
}

int stc_render_rgb_init(StcRgbBuf *b, int w, int h)
{
    size_t n = 0;

    if (b == NULL) {
        return -1;
    }
    b->w = 0;
    b->h = 0;
    b->px = NULL;
    if (render_pixel_count(w, h, &n) != 0) {
        return -1;
    }
    b->px = (StcRgb *)malloc(n * sizeof(StcRgb));
    if (b->px == NULL) {
        return -1;
    }
    b->w = w;
    b->h = h;
    memset(b->px, 0, n * sizeof(StcRgb));
    return 0;
}

void stc_render_rgb_free(StcRgbBuf *b)
{
    if (b == NULL) {
        return;
    }
    free(b->px);
    b->w = 0;
    b->h = 0;
    b->px = NULL;
}

static int render_rgb_ready(const StcRgbBuf *b)
{
    return b != NULL && b->w > 0 && b->h > 0 && b->px != NULL;
}

void stc_render_rgb_clear(StcRgbBuf *b, StcRgb c)
{
    size_t n, i;

    if (!render_rgb_ready(b)) {
        return;
    }
    n = (size_t)b->w * (size_t)b->h;
    for (i = 0; i < n; i++) {
        b->px[i] = c;
    }
}

void stc_render_rgb_plot(StcRgbBuf *b, int x, int y, StcRgb c)
{
    size_t i;

    if (!render_rgb_ready(b)) {
        return;
    }
    if (x < 0 || x >= b->w || y < 0 || y >= b->h) {
        return;
    }
    i = (size_t)y * (size_t)b->w + (size_t)x;
    b->px[i] = c;
}

void stc_render_rgb_line(StcRgbBuf *b, int x0, int y0, int x1, int y1,
                         StcRgb c)
{
    /* Equation: Bresenham integer line in 64-bit deltas (no int overflow). */
    long long dx = llabs((long long)x1 - (long long)x0);
    long long dy = -llabs((long long)y1 - (long long)y0);
    long long err = dx + dy;
    long long cap = dx - dy + 2; /* max plotted points: |dx|+|dy|+2 */
    long long i;
    int sx = x0 < x1 ? 1 : -1;
    int sy = y0 < y1 ? 1 : -1;
    int x = x0;
    int y = y0;

    if (!render_rgb_ready(b)) {
        return;
    }
    for (i = 0; i < cap; i++) {
        long long e2;

        stc_render_rgb_plot(b, x, y, c);
        if (x == x1 && y == y1) {
            break;
        }
        e2 = 2 * err;
        if (e2 >= dy) {
            err += dy;
            x += sx;
        }
        if (e2 <= dx) {
            err += dx;
            y += sy;
        }
    }

}

static double render_edge(double ax, double ay, double bx, double by,
                          double px, double py)
{
    /* Equation: 2D edge function E = (b-a) x (p-a) (z-component). */
    return (bx - ax) * (py - ay) - (by - ay) * (px - ax);
}

void stc_render_rgb_triangle(StcRgbBuf *b,
                             int x0, int y0, int x1, int y1, int x2, int y2,
                             StcRgb c)
{
    double area2;
    int lo_x, hi_x, lo_y, hi_y;
    int y, x;

    if (!render_rgb_ready(b)) {
        return;
    }
    area2 = render_edge((double)x0, (double)y0, (double)x1, (double)y1,
                        (double)x2, (double)y2);
    if (area2 == 0.0) {
        /* Degenerate: stroke the edges so the simplex stays visible. */
        stc_render_rgb_line(b, x0, y0, x1, y1, c);
        stc_render_rgb_line(b, x1, y1, x2, y2, c);
        stc_render_rgb_line(b, x2, y2, x0, y0, c);
        return;
    }
    /* Equation: edge-function half-space fill over a clipped bbox. */
    lo_x = x0 < x1 ? (x0 < x2 ? x0 : x2) : (x1 < x2 ? x1 : x2);
    hi_x = x0 > x1 ? (x0 > x2 ? x0 : x2) : (x1 > x2 ? x1 : x2);
    lo_y = y0 < y1 ? (y0 < y2 ? y0 : y2) : (y1 < y2 ? y1 : y2);
    hi_y = y0 > y1 ? (y0 > y2 ? y0 : y2) : (y1 > y2 ? y1 : y2);
    if (lo_x < 0) {
        lo_x = 0;
    }
    if (lo_y < 0) {
        lo_y = 0;
    }
    if (hi_x >= b->w) {
        hi_x = b->w - 1;
    }
    if (hi_y >= b->h) {
        hi_y = b->h - 1;
    }
    for (y = lo_y; y <= hi_y; y++) {
        for (x = lo_x; x <= hi_x; x++) {
            double px = (double)x;
            double py = (double)y;
            double e0 = render_edge((double)x0, (double)y0,
                                    (double)x1, (double)y1, px, py);
            double e1 = render_edge((double)x1, (double)y1,
                                    (double)x2, (double)y2, px, py);
            double e2 = render_edge((double)x2, (double)y2,
                                    (double)x0, (double)y0, px, py);

            if (area2 > 0.0) {
                if (e0 >= 0.0 && e1 >= 0.0 && e2 >= 0.0) {
                    stc_render_rgb_plot(b, x, y, c);
                }
            } else {
                if (e0 <= 0.0 && e1 <= 0.0 && e2 <= 0.0) {
                    stc_render_rgb_plot(b, x, y, c);
                }
            }
        }
    }
}

void stc_render_rgb_triangle_shaded(StcRgbBuf *b,
                                    int x0, int y0,
                                    int x1, int y1,
                                    int x2, int y2,
                                    StcRgb base, double shade)
{
    StcRgb c;
    double s = render_clamp01(shade);

    /* Equation: flat shading C = base * clamp(shade, 0, 1). */
    c.r = render_ramp_channel((double)base.r, s);
    c.g = render_ramp_channel((double)base.g, s);
    c.b = render_ramp_channel((double)base.b, s);
    stc_render_rgb_triangle(b, x0, y0, x1, y1, x2, y2, c);
}

int stc_render_write_ppm(const StcRgbBuf *b, const char *path)
{
    FILE *f;
    size_t n;

    if (!render_rgb_ready(b) || path == NULL) {
        return -1;
    }
    n = (size_t)b->w * (size_t)b->h;
    f = fopen(path, "wb");
    if (f == NULL) {
        return -1;
    }
    /* Equation: PPM P6 binary raster, top row first. */
    if (fprintf(f, "P6\n%d %d\n255\n", b->w, b->h) < 0) {
        fclose(f);
        return -1;
    }
    if (n > 0 && fwrite(b->px, sizeof(StcRgb), n, f) != n) {
        fclose(f);
        return -1;
    }
    if (fclose(f) != 0) {
        return -1;
    }
    return 0;
}

/* ---------- ANSI + palette ---------- */

const char *stc_render_ansi_escape(int color)
{
    /* Equation: ANSI SGR table lookup. */
    static const char *table[16] = {
        "\x1b[30m", "\x1b[31m", "\x1b[32m", "\x1b[33m",
        "\x1b[34m", "\x1b[35m", "\x1b[36m", "\x1b[37m",
        "\x1b[90m", "\x1b[91m", "\x1b[92m", "\x1b[93m",
        "\x1b[94m", "\x1b[95m", "\x1b[96m", "\x1b[97m"
    };

    if (!render_color_valid(color)) {
        return "";
    }
    return table[color];
}

StcRgb stc_render_phosphor(double t)
{
    /* Equation: phosphor P1 green ramp, C = (51, 255, 51) * clamp01(t). */
    return stc_render_rgb(render_ramp_channel(51.0, t),
                          render_ramp_channel(255.0, t),
                          render_ramp_channel(51.0, t));
}

StcRgb stc_render_amber(double t)
{
    /* Equation: amber ramp, C = (255, 176, 0) * clamp01(t). */
    return stc_render_rgb(render_ramp_channel(255.0, t),
                          render_ramp_channel(176.0, t),
                          render_ramp_channel(0.0, t));
}
