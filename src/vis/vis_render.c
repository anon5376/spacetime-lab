/* SPACETIME LAB — vis_render: portable visual-sandbox rasterizer.
 * Starfield, shaded body discs, fading trails, bitmap-font HUD, spawn-form
 * and help overlays. Pure C11, deterministic, no I/O except PPM export.
 */

#include "stc/vis.h"

#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "stc/body.h"
#include "stc/elements.h"
#include "stc/fmt.h"
#include "stc/help.h"
#include "stc/rk4.h"

/* 5x7 bitmap font, row-major, bit4 = left pixel, ASCII 0x20..0x7E. */
static const uint8_t kFont[95][7] = {
    {0x00,0x00,0x00,0x00,0x00,0x00,0x00}, /* space */
    {0x04,0x04,0x04,0x04,0x04,0x00,0x04}, /* ! */
    {0x0A,0x0A,0x0A,0x00,0x00,0x00,0x00}, /* " */
    {0x0A,0x0A,0x1F,0x0A,0x1F,0x0A,0x0A}, /* # */
    {0x04,0x0F,0x14,0x0E,0x05,0x1E,0x04}, /* $ */
    {0x18,0x19,0x02,0x04,0x08,0x13,0x03}, /* % */
    {0x0C,0x12,0x14,0x08,0x15,0x12,0x0D}, /* & */
    {0x04,0x04,0x08,0x00,0x00,0x00,0x00}, /* ' */
    {0x02,0x04,0x08,0x08,0x08,0x04,0x02}, /* ( */
    {0x08,0x04,0x02,0x02,0x02,0x04,0x08}, /* ) */
    {0x00,0x04,0x15,0x0E,0x15,0x04,0x00}, /* * */
    {0x00,0x04,0x04,0x1F,0x04,0x04,0x00}, /* + */
    {0x00,0x00,0x00,0x00,0x04,0x04,0x08}, /* , */
    {0x00,0x00,0x00,0x1F,0x00,0x00,0x00}, /* - */
    {0x00,0x00,0x00,0x00,0x00,0x00,0x04}, /* . */
    {0x01,0x02,0x02,0x04,0x08,0x08,0x10}, /* / */
    {0x0E,0x11,0x13,0x15,0x19,0x11,0x0E}, /* 0 */
    {0x04,0x0C,0x04,0x04,0x04,0x04,0x0E}, /* 1 */
    {0x0E,0x11,0x01,0x02,0x04,0x08,0x1F}, /* 2 */
    {0x1F,0x02,0x04,0x02,0x01,0x11,0x0E}, /* 3 */
    {0x02,0x06,0x0A,0x12,0x1F,0x02,0x02}, /* 4 */
    {0x1F,0x10,0x1E,0x01,0x01,0x11,0x0E}, /* 5 */
    {0x06,0x08,0x10,0x1E,0x11,0x11,0x0E}, /* 6 */
    {0x1F,0x01,0x02,0x04,0x08,0x08,0x08}, /* 7 */
    {0x0E,0x11,0x11,0x0E,0x11,0x11,0x0E}, /* 8 */
    {0x0E,0x11,0x11,0x0F,0x01,0x02,0x0C}, /* 9 */
    {0x00,0x00,0x04,0x00,0x00,0x04,0x00}, /* : */
    {0x00,0x00,0x04,0x00,0x00,0x04,0x08}, /* ; */
    {0x02,0x04,0x08,0x10,0x08,0x04,0x02}, /* < */
    {0x00,0x00,0x1F,0x00,0x1F,0x00,0x00}, /* = */
    {0x08,0x04,0x02,0x01,0x02,0x04,0x08}, /* > */
    {0x0E,0x11,0x01,0x02,0x04,0x00,0x04}, /* ? */
    {0x0E,0x11,0x01,0x0D,0x15,0x15,0x0E}, /* @ */
    {0x0E,0x11,0x11,0x1F,0x11,0x11,0x11}, /* A */
    {0x1E,0x11,0x11,0x1E,0x11,0x11,0x1E}, /* B */
    {0x0E,0x11,0x10,0x10,0x10,0x11,0x0E}, /* C */
    {0x1E,0x11,0x11,0x11,0x11,0x11,0x1E}, /* D */
    {0x1F,0x10,0x10,0x1E,0x10,0x10,0x1F}, /* E */
    {0x1F,0x10,0x10,0x1E,0x10,0x10,0x10}, /* F */
    {0x0E,0x11,0x10,0x17,0x11,0x11,0x0F}, /* G */
    {0x11,0x11,0x11,0x1F,0x11,0x11,0x11}, /* H */
    {0x0E,0x04,0x04,0x04,0x04,0x04,0x0E}, /* I */
    {0x07,0x02,0x02,0x02,0x02,0x12,0x0C}, /* J */
    {0x11,0x12,0x14,0x18,0x14,0x12,0x11}, /* K */
    {0x10,0x10,0x10,0x10,0x10,0x10,0x1F}, /* L */
    {0x11,0x1B,0x15,0x15,0x11,0x11,0x11}, /* M */
    {0x11,0x19,0x19,0x15,0x13,0x13,0x11}, /* N */
    {0x0E,0x11,0x11,0x11,0x11,0x11,0x0E}, /* O */
    {0x1E,0x11,0x11,0x1E,0x10,0x10,0x10}, /* P */
    {0x0E,0x11,0x11,0x11,0x15,0x12,0x0D}, /* Q */
    {0x1E,0x11,0x11,0x1E,0x14,0x12,0x11}, /* R */
    {0x0E,0x11,0x10,0x0E,0x01,0x11,0x0E}, /* S */
    {0x1F,0x04,0x04,0x04,0x04,0x04,0x04}, /* T */
    {0x11,0x11,0x11,0x11,0x11,0x11,0x0E}, /* U */
    {0x11,0x11,0x11,0x11,0x11,0x0A,0x04}, /* V */
    {0x11,0x11,0x11,0x15,0x15,0x1B,0x11}, /* W */
    {0x11,0x11,0x0A,0x04,0x0A,0x11,0x11}, /* X */
    {0x11,0x11,0x0A,0x04,0x04,0x04,0x04}, /* Y */
    {0x1F,0x01,0x02,0x04,0x08,0x10,0x1F}, /* Z */
    {0x0E,0x08,0x08,0x08,0x08,0x08,0x0E}, /* [ */
    {0x10,0x08,0x08,0x04,0x02,0x02,0x01}, /* \ */
    {0x0E,0x02,0x02,0x02,0x02,0x02,0x0E}, /* ] */
    {0x04,0x0A,0x11,0x00,0x00,0x00,0x00}, /* ^ */
    {0x00,0x00,0x00,0x00,0x00,0x00,0x1F}, /* _ */
    {0x04,0x08,0x00,0x00,0x00,0x00,0x00}, /* ` */
    {0x00,0x00,0x0E,0x01,0x0F,0x11,0x0F}, /* a */
    {0x10,0x10,0x1E,0x11,0x11,0x11,0x1E}, /* b */
    {0x00,0x00,0x0E,0x10,0x10,0x11,0x0E}, /* c */
    {0x01,0x01,0x0F,0x11,0x11,0x11,0x0F}, /* d */
    {0x00,0x00,0x0E,0x11,0x1F,0x10,0x0E}, /* e */
    {0x06,0x09,0x08,0x1C,0x08,0x08,0x08}, /* f */
    {0x00,0x0F,0x11,0x11,0x0F,0x01,0x0E}, /* g */
    {0x10,0x10,0x1E,0x11,0x11,0x11,0x11}, /* h */
    {0x04,0x00,0x0C,0x04,0x04,0x04,0x0E}, /* i */
    {0x02,0x00,0x06,0x02,0x02,0x12,0x0C}, /* j */
    {0x10,0x10,0x12,0x14,0x18,0x14,0x12}, /* k */
    {0x0C,0x04,0x04,0x04,0x04,0x04,0x0E}, /* l */
    {0x00,0x00,0x1A,0x15,0x15,0x15,0x15}, /* m */
    {0x00,0x00,0x1E,0x11,0x11,0x11,0x11}, /* n */
    {0x00,0x00,0x0E,0x11,0x11,0x11,0x0E}, /* o */
    {0x00,0x1E,0x11,0x11,0x1E,0x10,0x10}, /* p */
    {0x00,0x0F,0x11,0x11,0x0F,0x01,0x01}, /* q */
    {0x00,0x00,0x16,0x19,0x10,0x10,0x10}, /* r */
    {0x00,0x00,0x0F,0x10,0x0E,0x01,0x1E}, /* s */
    {0x08,0x08,0x1C,0x08,0x08,0x09,0x06}, /* t */
    {0x00,0x00,0x11,0x11,0x11,0x11,0x0F}, /* u */
    {0x00,0x00,0x11,0x11,0x11,0x0A,0x04}, /* v */
    {0x00,0x00,0x11,0x11,0x15,0x15,0x0A}, /* w */
    {0x00,0x00,0x11,0x0A,0x04,0x0A,0x11}, /* x */
    {0x00,0x11,0x11,0x11,0x0F,0x01,0x0E}, /* y */
    {0x00,0x00,0x1F,0x02,0x04,0x08,0x1F}, /* z */
    {0x02,0x04,0x04,0x08,0x04,0x04,0x02}, /* { */
    {0x04,0x04,0x04,0x04,0x04,0x04,0x04}, /* | */
    {0x08,0x04,0x04,0x02,0x04,0x04,0x08}, /* } */
    {0x00,0x09,0x16,0x00,0x00,0x00,0x00}  /* ~ */
};

/* --- 2D primitives -------------------------------------------------- */

static void vis_glyph(StcRgbBuf *fb, int x, int y, char c, StcRgb col,
                      int scale)
{
    const uint8_t *g;
    int r;
    int cc;

    if (fb == NULL || fb->px == NULL || scale < 1) {
        return;
    }
    g = kFont[(c >= 32 && c < 127) ? (size_t)(c - 32) : 0];
    for (r = 0; r < 7; r++) {
        for (cc = 0; cc < 5; cc++) {
            int dx;
            int dy;

            if (!((g[r] >> (4 - cc)) & 1u)) {
                continue;
            }
            for (dy = 0; dy < scale; dy++) {
                for (dx = 0; dx < scale; dx++) {
                    stc_render_rgb_plot(fb, x + cc * scale + dx,
                                        y + r * scale + dy, col);
                }
            }
        }
    }
}

/* Text run with '\n' support. Returns the y of the next line start. */
static int vis_text(StcRgbBuf *fb, int x, int y, const char *s, StcRgb col,
                    int scale)
{
    int cx = x;
    int step = 6 * (scale < 1 ? 1 : scale);
    int lh = 8 * (scale < 1 ? 1 : scale);

    if (fb == NULL || s == NULL) {
        return y;
    }
    if (scale < 1) {
        scale = 1;
    }
    for (; *s != '\0'; s++) {
        if (*s == '\n') {
            cx = x;
            y += lh;
            continue;
        }
        vis_glyph(fb, cx, y, *s, col, scale);
        cx += step;
    }
    return y + lh;
}

static void vis_rect(StcRgbBuf *fb, int x, int y, int w, int h, StcRgb col)
{
    int i;
    int j;

    if (fb == NULL || w <= 0 || h <= 0) {
        return;
    }
    for (j = 0; j < h; j++) {
        for (i = 0; i < w; i++) {
            stc_render_rgb_plot(fb, x + i, y + j, col);
        }
    }
}

static void vis_frame(StcRgbBuf *fb, int x, int y, int w, int h, StcRgb col)
{
    if (fb == NULL || w <= 0 || h <= 0) {
        return;
    }
    stc_render_rgb_line(fb, x, y, x + w - 1, y, col);
    stc_render_rgb_line(fb, x, y + h - 1, x + w - 1, y + h - 1, col);
    stc_render_rgb_line(fb, x, y, x, y + h - 1, col);
    stc_render_rgb_line(fb, x + w - 1, y, x + w - 1, y + h - 1, col);
}

/* Shaded disc, light from top: shade 1.0 top -> 0.35 bottom. */
static void vis_disc(StcRgbBuf *fb, int cx, int cy, int r, StcRgb base)
{
    int dy;

    if (fb == NULL || r <= 0) {
        return;
    }
    for (dy = -r; dy <= r; dy++) {
        double t = (double)(dy + r) / (double)(2 * r);
        double k = 1.0 - 0.65 * t;
        int half = (int)sqrt((double)(r * r - dy * dy));
        int dx;
        StcRgb c;

        c.r = (uint8_t)(base.r * k + 0.5);
        c.g = (uint8_t)(base.g * k + 0.5);
        c.b = (uint8_t)(base.b * k + 0.5);
        for (dx = -half; dx <= half; dx++) {
            stc_render_rgb_plot(fb, cx + dx, cy + dy, c);
        }
    }
}

/* 1px ring (midpoint circle outline). */
static void vis_ring(StcRgbBuf *fb, int cx, int cy, int r, StcRgb col)
{
    int x = r;
    int y = 0;
    int err = 0;

    if (fb == NULL || r <= 0) {
        return;
    }
    while (x >= y) {
        stc_render_rgb_plot(fb, cx + x, cy + y, col);
        stc_render_rgb_plot(fb, cx + y, cy + x, col);
        stc_render_rgb_plot(fb, cx - y, cy + x, col);
        stc_render_rgb_plot(fb, cx - x, cy + y, col);
        stc_render_rgb_plot(fb, cx - x, cy - y, col);
        stc_render_rgb_plot(fb, cx - y, cy - x, col);
        stc_render_rgb_plot(fb, cx + y, cy - x, col);
        stc_render_rgb_plot(fb, cx + x, cy - y, col);
        y++;
        err += 1 + 2 * y;
        if (2 * (err - x) + 1 > 0) {
            x--;
            err += 1 - 2 * x;
        }
    }
}

/* --- View fit + camera (same model as the terminal sandbox) ---------- */

typedef struct {
    double cx;
    double cy;
    double cz;
    double spread;
} VisFit;

static VisFit vis_fit(const StcSim *s, int follow)
{
    VisFit f;
    int i;
    int n = 0;

    f.cx = 0.0;
    f.cy = 0.0;
    f.cz = 0.0;
    if (s == NULL || s->n <= 0) {
        f.spread = 1.0;
        return f;
    }
    if (follow >= 0 && follow < s->n) {
        f.cx = s->b[follow].p.x;
        f.cy = s->b[follow].p.y;
        f.cz = s->b[follow].p.z;
    } else {
        for (i = 0; i < s->n; i++) {
            f.cx += s->b[i].p.x;
            f.cy += s->b[i].p.y;
            f.cz += s->b[i].p.z;
        }
        n = s->n > 0 ? s->n : 1;
        f.cx /= (double)n;
        f.cy /= (double)n;
        f.cz /= (double)n;
    }
    f.spread = 0.0;
    for (i = 0; i < s->n; i++) {
        double dx = s->b[i].p.x - f.cx;
        double dy = s->b[i].p.y - f.cy;
        double dz = s->b[i].p.z - f.cz;
        double d = sqrt(dx * dx + dy * dy + dz * dz);

        if (d > f.spread) {
            f.spread = d;
        }
    }
    if (!(f.spread > 0.0)) {
        f.spread = 1.0;
    }
    return f;
}

static StcCamera vis_camera(double yaw, double pitch, double dist, int w,
                            int h)
{
    StcVec3 eye;
    StcVec3 center;
    StcVec3 up;

    eye.x = dist * cos(pitch) * cos(yaw);
    eye.y = dist * sin(pitch);
    eye.z = dist * cos(pitch) * sin(yaw);
    center.x = 0.0;
    center.y = 0.0;
    center.z = 0.0;
    up.x = 0.0;
    up.y = 1.0;
    up.z = 0.0;
    return stc_render_camera(eye, center, up, 60.0 * STC_PI / 180.0,
                            (double)w / (double)h, 0.01, 100.0, w, h);
}

/* Same dispatch as the terminal sandbox: stc_sim_step_rk4 holds no
 * collision/cap pass itself, so stc_sim_poststep runs after each rk4
 * step — identical post-step physics under both integrators (R3). */
static void vis_step_sim(StcSim *s, double dt, int integrator)
{
    if (integrator == STC_APP_INTEG_RK4) {
        stc_sim_step_rk4(s, dt);
        stc_sim_poststep(s);
    } else {
        stc_sim_step(s, dt);
    }
}

/* Deterministic starfield: index-hashed positions, no RNG state. */
static void vis_stars(StcRgbBuf *fb)
{
    static const int kStars = 160;
    int i;

    if (fb == NULL || fb->w <= 0 || fb->h <= 0) {
        return;
    }
    for (i = 0; i < kStars; i++) {
        uint32_t h1 = (uint32_t)(i + 1) * 2654435761u;
        uint32_t h2 = h1 * 2246822519u + 3266489917u;
        uint32_t h3 = h2 * 3266489917u + 668265263u;
        int x = (int)(h1 % (uint32_t)fb->w);
        int y = (int)(h2 % (uint32_t)fb->h);
        int tier = (int)(h3 % 4u);
        StcRgb c;

        if (tier == 0) {
            c = stc_render_rgb(150, 170, 200);
        } else if (tier == 1) {
            c = stc_render_rgb(110, 130, 160);
        } else if (tier == 2) {
            c = stc_render_rgb(200, 190, 170);
        } else {
            c = stc_render_rgb(70, 85, 110);
        }
        stc_render_rgb_plot(fb, x, y, c);
    }
}

/* --- State ---------------------------------------------------------- */

int stc_vis_init(StcVis *v, int preset, unsigned long long seed,
                 int integrator, int w, int h)
{
    if (v == NULL) {
        return -1;
    }
    memset(v, 0, sizeof *v);
    v->bank.trails = NULL;
    v->bank.nbodies = 0;
    v->bank.depth = 0;
    v->bank_ok = 0;
    v->fb.px = NULL;
    v->fb.w = 0;
    v->fb.h = 0;
    stc_sim_init(&v->sim, STC_G, 1000.0, 0.4, 1);
    if (stc_preset_load(&v->sim, preset, (uint64_t)seed) != 0) {
        return -1;
    }
    stc_ui_camera_init(&v->cam);
    v->cam.yaw_rad = 0.6;
    v->cam.pitch_rad = 0.5;
    v->cam.dist_m = 3.2;
    stc_ui_camera_update(&v->cam);
    v->bank_ok = (stc_trail_bank_init(&v->bank, STC_BODY_MAX,
                                     STC_TRAIL_DEPTH_DEFAULT) == 0);
    if (stc_render_rgb_init(&v->fb, w, h) != 0) {
        stc_trail_bank_free(&v->bank);
        v->bank_ok = 0;
        return -2;
    }
    v->dt = 1000.0;
    v->warp_i = 0;
    v->follow = -1;
    v->paused = 0;
    v->forming = 0;
    v->help_on = 0;
    v->trail_on = 1;
    v->grid_on = 1;
    v->integrator = integrator;
    v->preset = preset;
    v->seed = seed;
    v->frames = 0;
    v->status[0] = '\0';
    snprintf(v->status, sizeof v->status, "press ? for key help");
    return 0;
}

void stc_vis_free(StcVis *v)
{
    if (v == NULL) {
        return;
    }
    stc_render_rgb_free(&v->fb);
    stc_trail_bank_free(&v->bank);
    v->bank_ok = 0;
}

int stc_vis_resize(StcVis *v, int w, int h)
{
    if (v == NULL) {
        return -1;
    }
    stc_render_rgb_free(&v->fb);
    if (stc_render_rgb_init(&v->fb, w, h) != 0) {
        return -1;
    }
    return 0;
}

void stc_vis_step(StcVis *v)
{
    static const int kWarps[3] = { 1, 10, 100 };
    int w;
    int i;

    if (v == NULL || v->paused || v->forming || v->help_on) {
        return;
    }
    for (w = 0; w < kWarps[v->warp_i]; w++) {
        vis_step_sim(&v->sim, v->dt, v->integrator);
    }
    if (v->trail_on && v->bank_ok) {
        for (i = 0; i < v->sim.n && i < STC_BODY_MAX; i++) {
            stc_trail_push(&v->bank, (size_t)i, v->sim.b[i].p.x,
                           v->sim.b[i].p.y, v->sim.b[i].p.z);
        }
    }
    v->frames++;
}

/* Body snapshot CSV: t_s + one row per body (same shape as terminal E). */
static int vis_write_csv(const char *path, const StcSim *s)
{
    FILE *fp;
    int i;

    if (path == NULL || s == NULL) {
        return -1;
    }
    fp = fopen(path, "w");
    if (fp == NULL) {
        return -1;
    }
    if (fprintf(fp, "t_s=%g\nname,m_kg,r_m,px_m,py_m,pz_m,vx_m_s,vy_m_s,"
               "vz_m_s,fixed\n", s->t) < 0) {
        fclose(fp);
        return -1;
    }
    for (i = 0; i < s->n; i++) {
        const StcBody *b = &s->b[i];

        if (fprintf(fp, "%s,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,"
                   "%.17g,%d\n", b->name, b->m, b->r, b->p.x, b->p.y,
                   b->p.z, b->v.x, b->v.y, b->v.z, b->fixed ? 1 : 0) < 0) {
            fclose(fp);
            return -1;
        }
    }
    return fclose(fp) == 0 ? 0 : -1;
}

/* --- Painters ------------------------------------------------------- */

static int vis_body_radius(double mass)
{
    double r;

    if (!(mass > 0.0)) {
        return 2;
    }
    r = 1.5 * cbrt(mass / 1.0e24);
    if (r < 2.0) {
        r = 2.0;
    }
    if (r > 14.0) {
        r = 14.0;
    }
    return (int)(r + 0.5);
}

typedef struct {
    int idx;
    int sx;
    int sy;
    double depth;
} VisProj;

static void vis_paint_trails(StcVis *v, const VisFit *fit,
                             const StcCamera *cam)
{
    int i;

    if (!v->trail_on || !v->bank_ok) {
        return;
    }
    for (i = 0; i < v->sim.n && i < STC_BODY_MAX; i++) {
        size_t count = stc_trail_count(&v->bank, (size_t)i);
        size_t pos;
        int px = -1000000;
        int py = -1000000;

        for (pos = 0; pos < count; pos++) {
            StcVec3 raw;
            StcVec3 p;
            int qx;
            int qy;
            double t;
            StcRgb c;

            if (stc_trail_get(&v->bank, (size_t)i, pos, &raw) != 0) {
                px = -1000000;
                continue;
            }
            p.x = (raw.x - fit->cx) / fit->spread;
            p.y = (raw.y - fit->cy) / fit->spread;
            p.z = (raw.z - fit->cz) / fit->spread;
            if (!stc_render_project_point(cam, p, &qx, &qy, NULL)) {
                px = -1000000;
                continue;
            }
            if (px > -1000000 && (px != qx || py != qy)) {
                t = count > 1 ? (double)pos / (double)(count - 1) : 1.0;
                c = stc_render_rgb((uint8_t)(10 + 30 * t),
                                   (uint8_t)(40 + 90 * t),
                                   (uint8_t)(20 + 40 * t));
                stc_render_rgb_line(&v->fb, px, py, qx, qy, c);
            }
            px = qx;
            py = qy;
        }
    }
}

/* --- Spacetime grid + orientation aids ------------------------------
 * ILLUSTRATIVE rubber-sheet metaphor, flat sharp chrome style: a square
 * lattice in the orbital plane (y = 0, centered on the view fit), each
 * node pushed along +y by a potential-shaped dent of the live bodies,
 *     h(x,z) = -K * sum_i  w_i * sigma^2 / (dx_i^2 + dz_i^2 + sigma^2),
 * a Lorentzian well: it has the monotone potential-dent shape but its
 * tail dies at ~5 sigma, so the funnel is LOCAL -- a plain
 * 1/sqrt(d^2+sigma^2) kernel at any view-scale sigma dishes the whole
 * sheet instead. Kernel width sigma is set in VIEW units (a couple of
 * lattice spacings, so the funnel spans ~5 nodes at any zoom) and mass
 * weights w_i = (m_i / M_max)^VIS_GRID_MPOW -- a compressed
 * quarter-power fraction: the heaviest body digs the deep funnel while
 * lighter bodies still leave visible dimples (linear m_i would put
 * planets ~4 decades under the Sun, i.e. invisible). K is then
 * normalized directly on the deepest SAMPLED node: |h|max over the
 * lattice == VIS_GRID_DIP of the view half-span. Arbitrary view-scale
 * depth -- shaped like the Newtonian potential, not equal to it; it is
 * NOT the spacetime metric, claims nothing about GR, and feeds nothing
 * back into the physics.
 */
#define VIS_GRID_N 33      /* uniform core lines each way */
#define VIS_GRID_EXT 14    /* doubled-coord far lines per side (to 16384x) */
#define VIS_GRID_M (VIS_GRID_N + 2 * VIS_GRID_EXT)
#define VIS_GRID_SPAN 1.5  /* core half-extent, x fit->spread */
#define VIS_GRID_DIP 0.18  /* deepest sampled dent, x fit->spread */
#define VIS_GRID_SOFT 0.08 /* kernel width, x full core span */
#define VIS_GRID_MPOW 0.25 /* mass-weight exponent (dynamic range) */

/* Mass weights: (m_i / M_max)^VIS_GRID_MPOW. */
static void vis_grid_weights(const StcVis *v, double *w)
{
    double mmax = 0.0;
    int i;

    for (i = 0; i < v->sim.n && i < STC_BODY_MAX; i++) {
        if (!v->sim.b[i].alive) {
            continue;
        }
        if (v->sim.b[i].m > mmax) {
            mmax = v->sim.b[i].m;
        }
    }
    for (i = 0; i < v->sim.n && i < STC_BODY_MAX; i++) {
        w[i] = (mmax > 0.0 && v->sim.b[i].alive)
                 ? pow(v->sim.b[i].m / mmax, VIS_GRID_MPOW)
                 : 0.0;
    }
}

/* Flat grey band for one grid segment: 4 discrete steps across this
 * frame's observed depth range, near -> far (chrome greys 150..96;
 * banded, not a smooth gradient). */
static StcRgb vis_grid_grey(double depth_ndc, double dmin, double dmax)
{
    double t = dmax > dmin ? (depth_ndc - dmin) / (dmax - dmin) : 0.5;
    int band = (int)(stc_clamp(t, 0.0, 1.0) * 4.0);
    int g;

    if (band > 3) {
        band = 3;
    }
    g = 150 - 18 * band;
    return stc_render_rgb((uint8_t)g, (uint8_t)g, (uint8_t)g);
}

/* Clip-space w for a point (view space depth): cw > 0 == in front of
 * the camera plane. Same matrix math as stc_render_project_point. */
static double vis_grid_clipw(const StcCamera *cam, StcVec3 p)
{
    const double *m = cam->viewproj.m;

    return m[3] * p.x + m[7] * p.y + m[11] * p.z + m[15];
}

/* Project like stc_render_project_point but meant for the open-ended
 * lattice: returns 1 iff in front of the camera (cw > 0), and the NDC
 * is clamped to [-4, 4] so far-off-screen nodes get bounded finite
 * pixel coords (Bresenham then clips the actual pixels). This is what
 * lets far lines run to the screen edge / horizon instead of stopping
 * at the first off-viewport node. */
static int vis_grid_project(const StcCamera *cam, StcVec3 w, int *sx,
                            int *sy, double *dep)
{
    const double *m = cam->viewproj.m;
    double cx = m[0] * w.x + m[4] * w.y + m[8] * w.z + m[12];
    double cy = m[1] * w.x + m[5] * w.y + m[9] * w.z + m[13];
    double cz = m[2] * w.x + m[6] * w.y + m[10] * w.z + m[14];
    double cw = m[3] * w.x + m[7] * w.y + m[11] * w.z + m[15];
    double nx;
    double ny;

    *sx = 0;
    *sy = 0;
    *dep = 0.0;
    if (!(cw > 0.0)) {
        return 0;
    }
    nx = stc_clamp(cx / cw, -4.0, 4.0);
    ny = stc_clamp(cy / cw, -4.0, 4.0);
    /* ndc z unclamped for banding only; nodes just in front of the
     * camera plane (cw -> 0+) would sink to -inf and crush the range,
     * so keep it in the valid [-1, 1]. */
    *dep = stc_clamp(cz / cw, -1.0, 1.0);
    *sx = (int)floor((nx * 0.5 + 0.5) * (double)cam->vp_w);
    *sy = (int)floor((1.0 - (ny * 0.5 + 0.5)) * (double)cam->vp_h);
    /* 2 = in front AND inside the viewport; only those nodes set the
     * depth-band range, so the extremes of the far sheet and the
     * under-camera region don't collapse the greys to one band. */
    return fabs(cx / cw) <= 1.0 && fabs(cy / cw) <= 1.0 ? 2 : 1;
}

/* One lattice segment: draws when both endpoints are in front of the
 * camera; when exactly one is, the segment is clipped at the camera
 * plane (interpolated world point at cw -> 0+, which projects beyond
 * the viewport in the correct direction and clamps off-screen). */
static void vis_grid_seg(StcVis *v, const StcCamera *cam, double dmin,
                         double dmax, StcVec3 wa, StcVec3 wb, int fa,
                         int fb, int ax, int ay, int bx, int by,
                         double dep)
{
    StcRgb c = vis_grid_grey(dep, dmin, dmax);

    if (fa && fb) {
        stc_render_rgb_line(&v->fb, ax, ay, bx, by, c);
    } else if (fa || fb) {
        StcVec3 pf = fa ? wa : wb;
        StcVec3 pb = fa ? wb : wa;
        int pfx = fa ? ax : bx;
        int pfy = fa ? ay : by;
        double cwf = vis_grid_clipw(cam, pf);
        double cwb = vis_grid_clipw(cam, pb);
        double t = cwf / (cwf - cwb);
        StcVec3 pc;
        int cx;
        int cy;
        double dep2;

        if (t < 0.0) {
            t = 0.0;
        }
        if (t > 1.0) {
            t = 1.0;
        }
        t *= 0.995; /* stay epsilon in front of the camera plane */
        pc.x = pf.x + (pb.x - pf.x) * t;
        pc.y = pf.y + (pb.y - pf.y) * t;
        pc.z = pf.z + (pb.z - pf.z) * t;
        (void)vis_grid_project(cam, pc, &cx, &cy, &dep2);
        stc_render_rgb_line(&v->fb, pfx, pfy, cx, cy, c);
    }
}

/* Kernel sum at (lx, lz), offsets from the fit center in meters. */
static double vis_grid_well(const StcVis *v, const VisFit *fit,
                            const double *w, double sigma,
                            double lx, double lz)
{
    double s = 0.0;
    double s2 = sigma * sigma;
    int i;

    for (i = 0; i < v->sim.n && i < STC_BODY_MAX; i++) {
        double dx;
        double dz;

        if (!v->sim.b[i].alive) {
            continue;
        }
        dx = lx - (v->sim.b[i].p.x - fit->cx);
        dz = lz - (v->sim.b[i].p.z - fit->cz);
        s += w[i] * s2 / (dx * dx + dz * dz + s2);
    }
    return s;
}

/* Orientation aids, tied to the same grid_on switch: a world-axis triad
 * top-right (project unit +X/+Y/+Z directions through the camera, x/y/z
 * labels in the bitmap font) and a bottom-right scale bar giving the
 * real-world span of exactly 100 px at the view-center depth. */
static void vis_paint_aids(StcVis *v, const VisFit *fit,
                           const StcCamera *cam)
{
    static const char kLbl[3] = { 'x', 'y', 'z' };
    StcRgb ink = stc_render_rgb(170, 170, 170);
    StcRgb lbl = stc_render_rgb(210, 210, 210);
    StcVec3 org;
    int c0x;
    int c0y;
    int a;

    org.x = 0.0;
    org.y = 0.0;
    org.z = 0.0;
    if (!stc_render_project_point(cam, org, &c0x, &c0y, NULL)) {
        return;
    }
    {
        int ax = v->fb.w - 58;
        int ay = 52;

        for (a = 0; a < 3; a++) {
            StcVec3 tip;
            int qx;
            int qy;
            double dx;
            double dy;
            double l;
            int ex;
            int ey;
            char lbuf[2];

            tip.x = a == 0 ? 0.3 : 0.0;
            tip.y = a == 1 ? 0.3 : 0.0;
            tip.z = a == 2 ? 0.3 : 0.0;
            if (!stc_render_project_point(cam, tip, &qx, &qy, NULL)) {
                continue;
            }
            dx = (double)(qx - c0x);
            dy = (double)(qy - c0y);
            l = sqrt(dx * dx + dy * dy);
            if (l < 1.5) {
                /* Axis points at the camera: draw a dot, still label it. */
                vis_rect(&v->fb, ax - 1, ay - 1, 3, 3, ink);
                lbuf[0] = kLbl[a];
                lbuf[1] = '\0';
                vis_text(&v->fb, ax + 6, ay - 4, lbuf, lbl, 1);
                continue;
            }
            ex = ax + (int)(dx / l * 22.0);
            ey = ay + (int)(dy / l * 22.0);
            stc_render_rgb_line(&v->fb, ax, ay, ex, ey, ink);
            lbuf[0] = kLbl[a];
            lbuf[1] = '\0';
            vis_text(&v->fb, ax + (int)(dx / l * 30.0) - 2,
                     ay + (int)(dy / l * 30.0) - 4, lbuf, lbl, 1);
        }
    }

    /* Scale bar: meters per pixel along the camera-right direction at
     * the view-center depth, times a fixed 100 px bar. */
    {
        StcVec3 fwd = stc_vec3_sub(cam->center, cam->eye);
        StcVec3 rgt;
        StcVec3 p1;
        int qx;
        int qy;
        double dpx;
        double raw;
        double unit;
        const char *uname;
        char buf[32];
        int y0;
        int x1;
        int tw;

        fwd = stc_vec3_norm(fwd);
        rgt = stc_vec3_norm(stc_vec3_cross(fwd, cam->up));
        p1 = stc_vec3_scale(rgt, 0.02);
        if (stc_render_project_point(cam, p1, &qx, &qy, NULL)) {
            dpx = hypot((double)(qx - c0x), (double)(qy - c0y));
            if (dpx > 0.5 && fit->spread > 0.0) {
                raw = 100.0 * 0.02 * fit->spread / dpx;
                if (raw >= 0.1 * STC_AU) {
                    unit = STC_AU;
                    uname = "AU";
                } else if (raw >= 1.0e6) {
                    unit = 1.0e3;
                    uname = "km";
                } else {
                    unit = 1.0;
                    uname = "m";
                }
                stc_fmt_sci(raw / unit, uname, buf, sizeof buf);
                y0 = v->fb.h - 24;
                x1 = v->fb.w - 40;
                stc_render_rgb_line(&v->fb, x1 - 100, y0, x1, y0, ink);
                stc_render_rgb_line(&v->fb, x1 - 100, y0 - 3, x1 - 100,
                                    y0 + 3, ink);
                stc_render_rgb_line(&v->fb, x1, y0 - 3, x1, y0 + 3, ink);
                tw = 6 * (int)strlen(buf);
                vis_text(&v->fb, x1 - 50 - tw / 2, y0 - 12, buf, lbl, 1);
            }
        }
    }
}

static void vis_paint_grid(StcVis *v, const VisFit *fit,
                           const StcCamera *cam)
{
    int px[VIS_GRID_M][VIS_GRID_M];
    int py[VIS_GRID_M][VIS_GRID_M];
    int pk[VIS_GRID_M][VIS_GRID_M];
    double pd[VIS_GRID_M][VIS_GRID_M];
    double ws[VIS_GRID_M][VIS_GRID_M]; /* kernel sums per node */
    double gw[STC_BODY_MAX];
    double ax[VIS_GRID_M]; /* node offsets from the fit center, meters */
    double half;
    double step;
    double sigma;
    double k;
    double sp;
    double smax;
    double dmin;
    double dmax;
    int na;
    int ix;
    int iz;
    int j;

    if (!v->grid_on || !(fit->spread > 0.0)) {
        return;
    }
    half = VIS_GRID_SPAN * fit->spread;
    step = 2.0 * half / (double)(VIS_GRID_N - 1);
    sigma = VIS_GRID_SOFT * 2.0 * half; /* ~2.6 core spacings wide */
    sp = fit->spread;
    /* Node coords along one axis: the uniform core, then the sheet
     * keeps going with geometric spacing -- node coordinate doubles
     * each step to 2^VIS_GRID_EXT * half (16384x the half-span). At
     * that radius the finite square edge sits within ~0.1 deg of the
     * horizon line from any sane camera and the far lines pack to
     * sub-pixel density there, so the sheet reads as unbounded rather
     * than stopping at a visible boundary. The dent tails
     * ~sigma^2/d^2 to zero out there, so the far sheet flattens to
     * the y=0 plane on its own. */
    na = 0;
    for (j = VIS_GRID_EXT; j >= 1; j--) {
        ax[na++] = -half * (double)(1 << j);
    }
    for (j = 0; j < VIS_GRID_N; j++) {
        ax[na++] = -half + step * (double)j;
    }
    for (j = 1; j <= VIS_GRID_EXT; j++) {
        ax[na++] = half * (double)(1 << j);
    }
    vis_grid_weights(v, gw);
    smax = 0.0;
    for (iz = 0; iz < na; iz++) {
        for (ix = 0; ix < na; ix++) {
            ws[iz][ix] = vis_grid_well(v, fit, gw, sigma, ax[ix],
                                       ax[iz]);
            if (ws[iz][ix] > smax) {
                smax = ws[iz][ix];
            }
        }
    }
    /* K from the deepest SAMPLED node: |h|max == DIP * view half-span. */
    k = smax > 0.0 ? VIS_GRID_DIP * sp / smax : 0.0;
    for (iz = 0; iz < na; iz++) {
        for (ix = 0; ix < na; ix++) {
            StcVec3 w;

            w.x = ax[ix] / sp;
            w.y = -k * ws[iz][ix] / sp;
            w.z = ax[iz] / sp;
            pk[iz][ix] = vis_grid_project(cam, w, &px[iz][ix],
                                          &py[iz][ix], &pd[iz][ix]);
        }
    }
    /* Band greys across the depths actually projected this frame. */
    dmin = 1.0;
    dmax = -1.0;
    for (iz = 0; iz < na; iz++) {
        for (ix = 0; ix < na; ix++) {
            if (pk[iz][ix] == 2) {
                if (pd[iz][ix] < dmin) {
                    dmin = pd[iz][ix];
                }
                if (pd[iz][ix] > dmax) {
                    dmax = pd[iz][ix];
                }
            }
        }
    }
    if (dmin > dmax) {
        /* Nothing inside the viewport this frame: fixed fallback. */
        dmin = 0.0;
        dmax = 1.0;
    }
    /* Lines of constant z (x runs) then constant x (z runs): per-segment
     * flat grey from the mean endpoint depth. */
    for (iz = 0; iz < na; iz++) {
        for (ix = 1; ix < na; ix++) {
            StcVec3 wa;
            StcVec3 wb;

            wa.x = ax[ix - 1] / sp;
            wa.y = -k * ws[iz][ix - 1] / sp;
            wa.z = ax[iz] / sp;
            wb.x = ax[ix] / sp;
            wb.y = -k * ws[iz][ix] / sp;
            wb.z = ax[iz] / sp;
            vis_grid_seg(v, cam, dmin, dmax, wa, wb, pk[iz][ix - 1],
                         pk[iz][ix], px[iz][ix - 1], py[iz][ix - 1],
                         px[iz][ix], py[iz][ix],
                         0.5 * (pd[iz][ix - 1] + pd[iz][ix]));
        }
    }
    for (ix = 0; ix < na; ix++) {
        for (iz = 1; iz < na; iz++) {
            StcVec3 wa;
            StcVec3 wb;

            wa.x = ax[ix] / sp;
            wa.y = -k * ws[iz - 1][ix] / sp;
            wa.z = ax[iz - 1] / sp;
            wb.x = ax[ix] / sp;
            wb.y = -k * ws[iz][ix] / sp;
            wb.z = ax[iz] / sp;
            vis_grid_seg(v, cam, dmin, dmax, wa, wb, pk[iz - 1][ix],
                         pk[iz][ix], px[iz - 1][ix], py[iz - 1][ix],
                         px[iz][ix], py[iz][ix],
                         0.5 * (pd[iz - 1][ix] + pd[iz][ix]));
        }
    }
    vis_paint_aids(v, fit, cam);
}

static void vis_paint_bodies(StcVis *v, const VisFit *fit,
                             const StcCamera *cam)
{
    VisProj projs[STC_BODY_MAX];
    int n = 0;
    int i;
    int j;
    StcRgb amber = stc_render_amber(0.9);
    StcRgb phos = stc_render_phosphor(1.0);
    StcRgb white = stc_render_rgb(220, 228, 240);

    for (i = 0; i < v->sim.n && i < STC_BODY_MAX; i++) {
        StcVec3 p;
        int sx;
        int sy;
        double depth;

        p.x = (v->sim.b[i].p.x - fit->cx) / fit->spread;
        p.y = (v->sim.b[i].p.y - fit->cy) / fit->spread;
        p.z = (v->sim.b[i].p.z - fit->cz) / fit->spread;
        if (!stc_render_project_point(cam, p, &sx, &sy, &depth)) {
            continue;
        }
        projs[n].idx = i;
        projs[n].sx = sx;
        projs[n].sy = sy;
        projs[n].depth = depth;
        n++;
    }
    /* Far-to-near (insertion sort, n <= 64). */
    for (i = 1; i < n; i++) {
        VisProj tmp = projs[i];

        j = i - 1;
        while (j >= 0 && projs[j].depth < tmp.depth) {
            projs[j + 1] = projs[j];
            j--;
        }
        projs[j + 1] = tmp;
    }
    for (i = 0; i < n; i++) {
        const StcBody *b = &v->sim.b[projs[i].idx];
        int r = vis_body_radius(b->m);
        StcRgb base = b->fixed ? amber : white;

        vis_disc(&v->fb, projs[i].sx, projs[i].sy, r, base);
        if (projs[i].idx == v->follow) {
            vis_ring(&v->fb, projs[i].sx, projs[i].sy, r + 3, phos);
            vis_text(&v->fb, projs[i].sx + r + 5, projs[i].sy - 4,
                     b->name, phos, 1);
        } else if (b->fixed) {
            vis_ring(&v->fb, projs[i].sx, projs[i].sy, r + 2, amber);
        }
    }
}

static void vis_paint_hud(StcVis *v)
{
    static const int kWarps[3] = { 1, 10, 100 };
    char line[256];
    char tbuf[32];
    char dtbuf[32];
    char ebuf[32];
    int y = 10;
    int i;
    int shown = 0;
    int nrows;
    int fe;
    /* Default-software chrome: flat 0xC0C0C0 face + 0x808080 1px edge,
     * near-black ink, dark-grey selection strip with light ink. */
    StcRgb panel = stc_render_rgb(192, 192, 192);
    StcRgb edge = stc_render_rgb(128, 128, 128);
    StcRgb txt = stc_render_rgb(16, 16, 16);
    StcRgb dim = stc_render_rgb(80, 80, 80);
    StcRgb sel = stc_render_rgb(64, 64, 64);
    StcRgb selink = stc_render_rgb(245, 245, 245);
    const char *iname = v->integrator == STC_APP_INTEG_RK4 ? "rk4"
                                                          : "verlet";

    vis_rect(&v->fb, 6, 6, 452, 30, panel);
    vis_frame(&v->fb, 6, 6, 452, 30, edge);
    stc_fmt_time(v->sim.t, tbuf, sizeof tbuf);
    stc_fmt_time(v->dt, dtbuf, sizeof dtbuf);
    snprintf(ebuf, sizeof ebuf, "%.3g", stc_sim_energy(&v->sim));
    snprintf(line, sizeof line, "%s t=%s dt=%s x%d %s n=%d E=%s%s",
             stc_preset_name(v->preset), tbuf, dtbuf, kWarps[v->warp_i],
             iname, v->sim.n, ebuf, v->paused ? " PAUSED" : "");
    vis_text(&v->fb, 12, 12, line, txt, 1);

    /* Body list on its own grey block (height tracks rows shown). */
    nrows = v->sim.n < 8 ? v->sim.n : 8;
    fe = (v->follow >= 0 && v->follow < v->sim.n && v->sim.n > 1) ? 12 : 0;
    vis_rect(&v->fb, 6, 40, 452, 10 * nrows + 6 + fe, panel);
    vis_frame(&v->fb, 6, 40, 452, 10 * nrows + 6 + fe, edge);
    y = 44;
    for (i = 0; i < v->sim.n && shown < 8; i++) {
        const StcBody *b = &v->sim.b[i];
        double spd = sqrt(b->v.x * b->v.x + b->v.y * b->v.y +
                          b->v.z * b->v.z);
        double pr = 1.0;

        if (v->sim.n > 0 && v->sim.b[0].m > 0.0) {
            double dx = b->p.x - v->sim.b[0].p.x;
            double dy = b->p.y - v->sim.b[0].p.y;
            double dz = b->p.z - v->sim.b[0].p.z;

            pr = stc_body_proper_rate(sqrt(dx * dx + dy * dy + dz * dz),
                                     v->sim.b[0].m);
        }
        snprintf(line, sizeof line, "%s%2d %-10.10s m=%.2g v=%.3g taurate=%.5f%s",
                 i == v->follow ? ">" : " ", i, b->name, b->m, spd, pr,
                 b->fixed ? " FIX" : "");
        if (i == v->follow) {
            vis_rect(&v->fb, 8, y - 1, 448, 9, sel);
            vis_text(&v->fb, 12, y, line, selink, 1);
        } else {
            vis_text(&v->fb, 12, y, line, dim, 1);
        }
        y += 10;
        shown++;
    }

    /* Follow elements. */
    if (v->follow >= 0 && v->follow < v->sim.n && v->sim.n > 1) {
        int center = (v->follow == 0) ? 1 : 0;
        double a = 0.0;
        double e = 0.0;
        double period = 0.0;

        stc_sim_elements(&v->sim, v->follow, center, &a, &e, &period);
        if (period > 0.0) {
            stc_fmt_time(period, tbuf, sizeof tbuf);
            snprintf(line, sizeof line, "follow #%d: a=%.4g e=%.4g T=%s",
                     v->follow, a, e, tbuf);
        } else {
            snprintf(line, sizeof line, "follow #%d: a=%.4g e=%.4g unbound",
                     v->follow, a, e);
        }
        vis_rect(&v->fb, 8, y + 1, 448, 9, sel);
        vis_text(&v->fb, 12, y + 2, line, selink, 1);
        y += 12;
    }

    /* Status + key hints on a bottom grey block. */
    vis_rect(&v->fb, 6, v->fb.h - 42, 452, 40, panel);
    vis_frame(&v->fb, 6, v->fb.h - 42, 452, 40, edge);
    if (v->status[0] != '\0') {
        snprintf(line, sizeof line, "last: %s", v->status);
        vis_text(&v->fb, 12, v->fb.h - 34, line, dim, 1);
    }
    vis_text(&v->fb, 12, v->fb.h - 20,
             "tab follow N spawn K remove G preset D grid T trails",
             dim, 1);
    vis_text(&v->fb, 12, v->fb.h - 10,
             "Z/X dt C warp P pause E export ? help Q quit", dim, 1);
}

static void vis_paint_form(StcVis *v)
{
    SbxForm *f = &v->form;
    int i;
    int bw = 460;
    int bh = 15 * 12 + 56;
    int bx = (v->fb.w - bw) / 2;
    int by = (v->fb.h - bh) / 2;
    int y;
    char line[160];
    StcRgb panel = stc_render_rgb(192, 192, 192);
    StcRgb edge = stc_render_rgb(128, 128, 128);
    StcRgb txt = stc_render_rgb(16, 16, 16);
    StcRgb sel = stc_render_rgb(64, 64, 64);
    StcRgb selink = stc_render_rgb(245, 245, 245);

    if (bx < 0) {
        bx = 0;
    }
    if (by < 0) {
        by = 0;
    }
    vis_rect(&v->fb, bx, by, bw, bh, panel);
    vis_frame(&v->fb, bx, by, bw, bh, edge);
    vis_rect(&v->fb, bx + 2, by + 4, bw - 4, 14, sel);
    y = by + 10;
    y = vis_text(&v->fb, bx + 12, y, "SPAWN BODY", selink, 1) + 2;
    for (i = 0; i < SBX_NFIELDS; i++) {
        if (i == SBX_NFIELDS - 1) {
            snprintf(line, sizeof line, "%s [SPAWN>]  enter here to spawn",
                     f->field == i ? ">" : " ");
        } else {
            snprintf(line, sizeof line, "%s %-8s : %s%s",
                     f->field == i ? ">" : " ", sbx_field_name(i),
                     f->buf[i], f->field == i ? "_" : "");
        }
        if (f->field == i) {
            vis_rect(&v->fb, bx + 8, y - 1, bw - 16, 9, sel);
            y = vis_text(&v->fb, bx + 12, y, line, selink, 1);
        } else {
            y = vis_text(&v->fb, bx + 12, y, line, txt, 1);
        }
    }
    if (f->err[0] != '\0') {
        snprintf(line, sizeof line, "! %s", f->err);
        vis_text(&v->fb, bx + 12, y + 2, line,
                 stc_render_rgb(150, 20, 20), 1);
    }
}

static void vis_paint_help(StcVis *v)
{
    const char *guide = stc_help_sandbox_keys();
    const char *nl = strchr(guide, '\n');
    char title[64];
    size_t tl;
    int bw = 560;
    int bh = 19 * 10 + 20;
    int bx = (v->fb.w - bw) / 2;
    int by = (v->fb.h - bh) / 2;
    StcRgb panel = stc_render_rgb(192, 192, 192);
    StcRgb edge = stc_render_rgb(128, 128, 128);
    StcRgb txt = stc_render_rgb(16, 16, 16);
    StcRgb sel = stc_render_rgb(64, 64, 64);
    StcRgb selink = stc_render_rgb(245, 245, 245);

    if (bx < 0) {
        bx = 0;
    }
    if (by < 0) {
        by = 0;
    }
    vis_rect(&v->fb, bx, by, bw, bh, panel);
    vis_frame(&v->fb, bx, by, bw, bh, edge);
    /* First guide line becomes the dark title strip; rest unchanged. */
    tl = nl != NULL ? (size_t)(nl - guide) : strlen(guide);
    if (tl >= sizeof title) {
        tl = sizeof title - 1;
    }
    memcpy(title, guide, tl);
    title[tl] = '\0';
    vis_rect(&v->fb, bx + 2, by + 4, bw - 4, 14, sel);
    vis_text(&v->fb, bx + 12, by + 10, title, selink, 1);
    vis_text(&v->fb, bx + 12, by + 18, nl != NULL ? nl + 1 : "", txt, 1);
}

static void vis_render_scene(StcVis *v);

void stc_vis_render(StcVis *v)
{
    const char *calib;

    if (v == NULL || v->fb.px == NULL) {
        return;
    }
    calib = getenv("STC_VIS_CALIB");
    if (calib != NULL && calib[0] != '\0') {
        /* TEMPORARY display-path calibration: solid thirds + crosshair. */
        int y;
        int x;

        for (y = 0; y < v->fb.h; y++) {
            StcRgb c;

            if (y < v->fb.h / 3) {
                c = stc_render_rgb(255, 0, 0);
            } else if (y < 2 * v->fb.h / 3) {
                c = stc_render_rgb(0, 255, 0);
            } else {
                c = stc_render_rgb(0, 0, 255);
            }
            for (x = 0; x < v->fb.w; x++) {
                stc_render_rgb_plot(&v->fb, x, y, c);
            }
        }
        for (x = 0; x < v->fb.w; x++) {
            stc_render_rgb_plot(&v->fb, x, v->fb.h / 2,
                                stc_render_rgb(255, 255, 255));
        }
        for (y = 0; y < v->fb.h; y++) {
            stc_render_rgb_plot(&v->fb, v->fb.w / 2, y,
                                stc_render_rgb(255, 255, 255));
        }
        return;
    }
    vis_render_scene(v);
}

/* Full scene painter. File-producing paths (stc_vis_shot) call this
 * directly so env state can never reach exported bytes (R4). */
static void vis_render_scene(StcVis *v)
{
    VisFit fit;
    StcCamera cam;

    stc_render_rgb_clear(&v->fb, stc_render_rgb(3, 6, 11));
    vis_stars(&v->fb);
    fit = vis_fit(&v->sim, v->follow);
    cam = vis_camera(v->cam.yaw_rad, v->cam.pitch_rad, v->cam.dist_m,
                     v->fb.w, v->fb.h);
    vis_paint_trails(v, &fit, &cam);
    vis_paint_grid(v, &fit, &cam);
    vis_paint_bodies(v, &fit, &cam);
    vis_paint_hud(v);
    if (v->forming) {
        vis_paint_form(v);
    } else if (v->help_on) {
        vis_paint_help(v);
    }
}

/* --- Key machine ---------------------------------------------------- */

static void vis_trail_clear(StcVis *v)
{
    if (v->bank_ok) {
        stc_trail_clear(&v->bank);
    }
}

static void vis_preset_load(StcVis *v, int preset)
{
    stc_sim_init(&v->sim, STC_G, 1000.0, 0.4, 1);
    if (stc_preset_load(&v->sim, preset, (uint64_t)v->seed) == 0) {
        v->preset = preset;
    }
    v->follow = -1;
    vis_trail_clear(v);
    snprintf(v->status, sizeof v->status, "preset=%s",
             stc_preset_name(v->preset));
}

int stc_vis_key(StcVis *v, int key)
{
    if (v == NULL) {
        return 0;
    }
    if (key == STC_UI_KEY_ERROR) {
        snprintf(v->status, sizeof v->status, "input error");
        return 0;
    }
    if (key == 'Q' || key == 'q' || key == STC_UI_KEY_CTRL_C) {
        return 1;
    }
    if (v->forming) {
        int rc = sbx_form_key(&v->sim, &v->form, key);

        if (rc > 0) {
            v->forming = 0;
            v->follow = v->sim.n - 1;
            vis_trail_clear(v);
            snprintf(v->status, sizeof v->status, "spawned %s (#%d)",
                     v->sim.b[v->follow].name, v->follow);
        } else if (rc < 0) {
            v->forming = 0;
            snprintf(v->status, sizeof v->status, "spawn cancelled");
        }
        return 0;
    }
    if (key == '?' && !v->forming) {
        v->help_on = !v->help_on;
        v->status[0] = '\0';
        return 0;
    }
    if (key == STC_UI_KEY_ESC) {
        return 1;
    }
    if (v->help_on) {
        v->help_on = 0;
        v->status[0] = '\0';
        return 0;
    }
    if (key == 'P' || key == 'p') {
        v->paused = !v->paused;
    } else if (key == 'N' || key == 'n') {
        VisFit fit = vis_fit(&v->sim, v->follow);

        v->forming = 1;
        sbx_form_prefill(&v->form, &v->sim, fit.cx, fit.cy, fit.cz,
                         fit.spread);
    } else if (key == 'K' || key == 'k') {
        if (v->follow >= 0 && v->follow < v->sim.n) {
            snprintf(v->status, sizeof v->status, "removed %s",
                     v->sim.b[v->follow].name);
            stc_sim_remove(&v->sim, v->follow);
            vis_trail_clear(v);
            if (v->follow >= v->sim.n) {
                v->follow = v->sim.n - 1;
            }
        }
    } else if (key == 'G' || key == 'g') {
        vis_preset_load(v, (v->preset + 1) % stc_preset_count());
    } else if (key == '\t') {
        v->follow = v->sim.n > 0 ? (v->follow + 2) % (v->sim.n + 1) - 1
                                : -1;
    } else if (key == 'Z' || key == 'z') {
        v->dt = v->dt / 2.0 < 1e-6 ? 1e-6 : v->dt / 2.0;
    } else if (key == 'X' || key == 'x') {
        v->dt = v->dt * 2.0 > 1e9 ? 1e9 : v->dt * 2.0;
    } else if (key == 'C' || key == 'c') {
        v->warp_i = (v->warp_i + 1) % 3;
    } else if (key == 'T' || key == 't') {
        v->trail_on = !v->trail_on;
        if (!v->trail_on) {
            vis_trail_clear(v);
        }
    } else if (key == 'D' || key == 'd') {
        /* Steals the WASD-D yaw-right binding in the visual sandbox only
         * (terminal untouched; Right arrow still yaws). */
        v->grid_on = !v->grid_on;
        snprintf(v->status, sizeof v->status, "grid sheet %s",
                 v->grid_on ? "on" : "off");
    } else if (key == 'E' || key == 'e') {
        const char *cpath = v->csv_path != NULL ? v->csv_path
                                                : "sandbox.csv";
        const char *ppath = v->ppm_path != NULL ? v->ppm_path
                                                : "sandbox.ppm";
        int okc;
        int okp;

        vis_render_scene(v);
        okc = vis_write_csv(cpath, &v->sim) == 0;
        okp = stc_render_write_ppm(&v->fb, ppath) == 0;
        snprintf(v->status, sizeof v->status, "export csv=%s ppm=%s",
                 okc ? cpath : "FAIL", okp ? ppath : "FAIL");
    } else if (key != STC_UI_KEY_NONE) {
        stc_ui_camera_key(&v->cam, key);
    }
    return 0;
}

/* --- Headless shot -------------------------------------------------- */

int stc_vis_shot(int preset, unsigned long long seed, int integrator,
                 long nframes, int w, int h, const char *ppm_path)
{
    StcVis v;
    long f;

    if (ppm_path == NULL || nframes < 0) {
        return -1;
    }
    if (stc_vis_init(&v, preset, seed, integrator, w, h) != 0) {
        return -2;
    }
    for (f = 0; f < nframes; f++) {
        stc_vis_step(&v);
    }
    vis_render_scene(&v);
    if (stc_render_write_ppm(&v.fb, ppm_path) != 0) {
        stc_vis_free(&v);
        return -3;
    }
    stc_vis_free(&v);
    return 0;
}

const unsigned char *stc_vis_glyph(char c)
{
    return kFont[(c >= 32 && c < 127) ? (size_t)(c - 32) : 0];
}
