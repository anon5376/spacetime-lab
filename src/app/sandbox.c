/* SPACETIME LAB — sandbox: interactive N-body gravity playground.
 * Includes core + body + render + ui + own only. SI doubles; rendering uses
 * a normalized unitless view fit (labels stay SI). Deterministic: presets
 * seed from --seed, spawn state is user-driven, headless runs are repeatable.
 */

#include "stc/sandbox.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "stc/body.h"
#include "stc/fmt.h"
#include "stc/help.h"
#include "stc/render.h"
#include "stc/spawn.h"
#include "stc/vis.h"
#include "stc/rk4.h"
#include "stc/trails.h"
#include "stc/ui.h"
#include "stc/units.h"
#include "stc/vec.h"

/* Trails: per-body ring buffers (stc/trails.h, depth 128), cleared on
 * structural change. Initialized per run; if allocation fails the bank
 * stays unusable and push/paint degrade to trails-off (no crash, no leak). */
static StcTrailBank g_bank;
static int g_bank_ok = 0;
static int g_trail_on = 1;

static void sbx_trail_init(void)
{
    g_bank.trails = NULL;
    g_bank.nbodies = 0;
    g_bank.depth = 0;
    g_bank_ok = (stc_trail_bank_init(&g_bank, STC_BODY_MAX,
                                    STC_TRAIL_DEPTH_DEFAULT) == 0);
}

static void sbx_trail_free(void)
{
    stc_trail_bank_free(&g_bank);
    g_bank_ok = 0;
}

static void sbx_trail_clear(void)
{
    if (g_bank_ok) {
        stc_trail_clear(&g_bank);
    }
}

static void sbx_trail_push(const StcSim *s)
{
    int i;

    if (!g_trail_on || !g_bank_ok || s == NULL) {
        return;
    }
    for (i = 0; i < s->n && i < STC_BODY_MAX; i++) {
        stc_trail_push(&g_bank, (size_t)i, s->b[i].p.x, s->b[i].p.y,
                       s->b[i].p.z);
    }
}

/* Equation: sandbox step dispatch — velocity Verlet (default) or classic
 * RK4 per --integrator. stc_sim_step_rk4 performs no collision pass and
 * no vmax cap itself (stc/rk4.h), so the shared post-step runs after it:
 * stc_sim_poststep applies the identical merge/bounce + speed-cap policy
 * stc_sim_step does internally — both integrators differ only in the
 * stepping scheme, never in collision/cap physics (R3). */
static void sbx_step(StcSim *s, double dt, int integrator)
{
    if (integrator == STC_APP_INTEG_RK4) {
        stc_sim_step_rk4(s, dt);
        stc_sim_poststep(s);
    } else {
        stc_sim_step(s, dt);
    }
}

static const char *sbx_integ_name(int integrator)
{
    return (integrator == STC_APP_INTEG_RK4) ? "rk4" : "verlet";
}

/* Normalized view fit around a center (followed body or centroid). */
typedef struct {
    double cx;
    double cy;
    double cz;
    double spread;
} SbxFit;

static SbxFit sbx_fit(const StcSim *s, int follow)
{
    SbxFit f;
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

static StcVec3 sbx_norm(double x, double y, double z, const SbxFit *f)
{
    StcVec3 p;

    p.x = (x - f->cx) / f->spread;
    p.y = (y - f->cy) / f->spread;
    p.z = (z - f->cz) / f->spread;
    return p;
}

static StcCamera sbx_camera(double yaw, double pitch, double dist, int w,
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

/* Human time rendering uses stc_fmt_time (stc/fmt.h, own core helper). */

/* Equation: frame pacing — busy-wait on clock() until dt = 1/60 s elapses. */
static void sbx_pace(clock_t start, int no_adapt)
{
    double elapsed;

    if (no_adapt || start == (clock_t)-1) {
        return;
    }
    do {
        clock_t now = clock();

        if (now == (clock_t)-1) {
            return;
        }
        elapsed = (double)(now - start) / (double)CLOCKS_PER_SEC;
    } while (elapsed < STC_APP_FRAME_DT);
}

/* ASCII viewport: trails dim, bodies bright, followed body '@'. */
static int sbx_present_ascii(const StcSim *s, const SbxFit *fit,
                            const StcCamera *cam, int follow, FILE *out)
{
    StcAsciiBuf fb;
    int i;
    int rc;

    if (s == NULL || fit == NULL || cam == NULL || out == NULL) {
        return -1;
    }
    if (stc_render_ascii_init(&fb, STC_APP_ASCII_W, STC_APP_ASCII_H) != 0) {
        return -2;
    }
    stc_render_ascii_clear(&fb, ' ', STC_RENDER_DEFAULT);
    if (g_trail_on && g_bank_ok) {
        for (i = 0; i < s->n && i < STC_BODY_MAX; i++) {
            size_t count = stc_trail_count(&g_bank, (size_t)i);
            size_t pos;
            int px = -1000000;
            int py = -1000000;

            for (pos = 0; pos < count; pos++) {
                StcVec3 raw;
                StcVec3 p;
                int qx;
                int qy;

                if (stc_trail_get(&g_bank, (size_t)i, pos, &raw) != 0) {
                    px = -1000000;
                    continue;
                }
                p = sbx_norm(raw.x, raw.y, raw.z, fit);
                if (!stc_render_project_point(cam, p, &qx, &qy, NULL)) {
                    px = -1000000;
                    continue;
                }
                if (px > -1000000 && (px != qx || py != qy)) {
                    stc_render_ascii_line(&fb, px, py, qx, qy, '.',
                                         STC_RENDER_GREEN);
                }
                px = qx;
                py = qy;
            }
        }
    }
    for (i = 0; i < s->n; i++) {
        StcVec3 p = sbx_norm(s->b[i].p.x, s->b[i].p.y, s->b[i].p.z, fit);
        int sx;
        int sy;

        if (!stc_render_project_point(cam, p, &sx, &sy, NULL)) {
            continue;
        }
        stc_render_ascii_plot(&fb, sx, sy, i == follow ? '@' : 'o',
                             i == follow ? STC_RENDER_BRIGHT_GREEN
                             : (s->b[i].fixed ? STC_RENDER_BRIGHT_YELLOW
                                              : STC_RENDER_BRIGHT_WHITE));
    }
    rc = stc_render_ascii_present(&fb, out);
    stc_render_ascii_free(&fb);
    return rc;
}

/* PPM export painter: same projection, RGB points + trail polylines. */
static int sbx_write_ppm(const char *path, const StcSim *s, const SbxFit *fit,
                        const StcCamera *cam, int follow)
{
    StcRgbBuf fb;
    StcRgb amber;
    StcRgb phos;
    StcRgb white;
    StcRgb dim;
    int i;
    int rc;

    if (path == NULL || s == NULL || fit == NULL || cam == NULL) {
        return -1;
    }
    if (stc_render_rgb_init(&fb, STC_APP_PPM_W, STC_APP_PPM_H) != 0) {
        return -2;
    }
    stc_render_rgb_clear(&fb, stc_render_rgb(4, 8, 12));
    amber = stc_render_amber(0.9);
    phos = stc_render_phosphor(1.0);
    white = stc_render_rgb(220, 228, 240);
    dim = stc_render_rgb(20, 60, 30);
    if (g_trail_on && g_bank_ok) {
        for (i = 0; i < s->n && i < STC_BODY_MAX; i++) {
            size_t count = stc_trail_count(&g_bank, (size_t)i);
            size_t pos;
            int px = -1000000;
            int py = -1000000;

            for (pos = 0; pos < count; pos++) {
                StcVec3 raw;
                StcVec3 p;
                int qx;
                int qy;

                if (stc_trail_get(&g_bank, (size_t)i, pos, &raw) != 0) {
                    px = -1000000;
                    continue;
                }
                p = sbx_norm(raw.x, raw.y, raw.z, fit);
                if (!stc_render_project_point(cam, p, &qx, &qy, NULL)) {
                    px = -1000000;
                    continue;
                }
                if (px > -1000000 && (px != qx || py != qy)) {
                    stc_render_rgb_line(&fb, px, py, qx, qy, dim);
                }
                px = qx;
                py = qy;
            }
        }
    }
    for (i = 0; i < s->n; i++) {
        StcVec3 p = sbx_norm(s->b[i].p.x, s->b[i].p.y, s->b[i].p.z, fit);
        StcRgb c = (i == follow) ? phos : (s->b[i].fixed ? amber : white);
        int sx;
        int sy;
        int dx;
        int dy;

        if (!stc_render_project_point(cam, p, &sx, &sy, NULL)) {
            continue;
        }
        for (dy = -1; dy <= 1; dy++) {
            for (dx = -1; dx <= 1; dx++) {
                stc_render_rgb_plot(&fb, sx + dx, sy + dy, c);
            }
        }
    }
    rc = stc_render_write_ppm(&fb, path);
    stc_render_rgb_free(&fb);
    return rc == 0 ? 0 : -2;
}

/* Equation: body snapshot CSV — one row per body (names quoted raw; names
 * with commas/quotes are sanitized to spaces on spawn). */
static int sbx_write_bodies_csv(const char *path, const StcSim *s)
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

/* --- Interactive sandbox screen --- */

/* Help overlay: splash title + key guide ('?' toggles, any key returns). */
static void sbx_draw_help(void)
{
    printf("\x1b[2J\x1b[H");
    printf("%s\n%s", stc_splash(), stc_help_sandbox_keys());
    printf("\npress ? or any key to return\n");
    fflush(stdout);
}

static void sbx_draw(const StcSim *s, const StcUiCamera *cam, int follow,
                    double dt, int warp, int paused, int integrator,
                    const char *status)
{
    char tbuf[32];
    char dtbuf[32];
    SbxFit fit;
    StcCamera rcam;
    int i;
    int shown;

    stc_fmt_time(s->t, tbuf, sizeof tbuf);
    stc_fmt_time(dt, dtbuf, sizeof dtbuf);
    printf("\x1b[2J\x1b[H");
    printf("sandbox  t=%s  dt=%s x%d  integ=%s  bodies=%d  E=%.6g J%s\n",
           tbuf, dtbuf, warp, sbx_integ_name(integrator), s->n,
           stc_sim_energy(s), paused ? "  [PAUSED]" : "");
    fit = sbx_fit(s, follow);
    rcam = sbx_camera(cam->yaw_rad, cam->pitch_rad, cam->dist_m,
                     STC_APP_ASCII_W, STC_APP_ASCII_H);
    (void)sbx_present_ascii(s, &fit, &rcam, follow, stdout);
    shown = 0;
    for (i = 0; i < s->n && shown < 8; i++) {
        const StcBody *b = &s->b[i];
        double spd = sqrt(b->v.x * b->v.x + b->v.y * b->v.y +
                         b->v.z * b->v.z);
        double pr = 1.0;

        if (s->n > 0 && s->b[0].m > 0.0) {
            double dx = b->p.x - s->b[0].p.x;
            double dy = b->p.y - s->b[0].p.y;
            double dz = b->p.z - s->b[0].p.z;

            pr = stc_body_proper_rate(sqrt(dx * dx + dy * dy + dz * dz),
                                     s->b[0].m);
        }
        printf("%s%2d %-10.10s m=%.3g r=%.3g |v|=%.4g taurate=%.6f%s\n",
               i == follow ? ">" : " ", i, b->name, b->m, b->r, spd, pr,
               b->fixed ? " FIX" : "");
        shown++;
    }
    if (follow >= 0 && follow < s->n && s->n > 1) {
        int center = (follow == 0) ? 1 : 0;
        double a = 0.0;
        double e = 0.0;
        double period = 0.0;

        stc_sim_elements(s, follow, center, &a, &e, &period);
        printf("follow #%d about #%d: a=%.4g m e=%.4g T=%s\n", follow,
               center, a, e,
               period > 0.0 ? (stc_fmt_time(period, tbuf, sizeof tbuf), tbuf)
                            : "unbound");
    }
    if (status[0] != '\0') {
        printf("last: %s\n", status);
    }
    printf("WASD/arrows+zoom  tab follow  N spawn  K remove  G preset  Z/X dt"
           "  C warp  T trails  P pause  E export  ? help  Q quit\n");
    fflush(stdout);
}

static int sbx_run_interactive(const StcAppOpts *o, int preset,
                              unsigned long long seed)
{
    static const int kWarps[3] = { 1, 10, 100 };
    StcSim sim;
    StcUiCamera cam;
    char status[256];
    double dt = 1000.0;
    int warp_i = 0;
    int follow = -1;
    int paused = 0;
    int forming = 0;
    int help_on = 0;
    int integrator = (o != NULL) ? o->integrator : STC_APP_INTEG_VERLET;
    SbxForm form;
    long cap = o->frames_given ? o->frames : -1L;
    long f = 0;

    stc_sim_init(&sim, STC_G, 1000.0, 0.4, 1);
    if (stc_preset_load(&sim, preset, (uint64_t)seed) != 0) {
        fprintf(stderr, "stlab: bad preset %d (see --help)\n", preset);
        return 2;
    }
    if (stc_ui_enter_raw() != 0) {
        return -99; /* caller falls back to headless */
    }
    atexit(stc_ui_exit_raw);
    stc_ui_camera_init(&cam);
    cam.yaw_rad = 0.6;
    cam.pitch_rad = 0.5;
    cam.dist_m = 3.2;
    stc_ui_camera_update(&cam);
    sbx_trail_init();
    snprintf(status, sizeof status, "press ? for key help");

    for (;;) {
        clock_t t0 = clock();
        int key;

        if (help_on && !forming) {
            sbx_draw_help();
        } else if (forming) {
            sbx_form_draw(&form);
        } else {
            sbx_draw(&sim, &cam, follow, dt, kWarps[warp_i], paused,
                     integrator, status);
        }
        key = stc_ui_poll_key();
        if (key == STC_UI_KEY_ERROR) {
            snprintf(status, sizeof status, "input error");
        } else if (key == 'Q' || key == 'q' ||
                   key == STC_UI_KEY_CTRL_C) {
            break;
        } else if (forming) {
            int rc = sbx_form_key(&sim, &form, key);

            if (rc > 0) {
                forming = 0;
                follow = sim.n - 1;
                sbx_trail_clear();
                snprintf(status, sizeof status, "spawned %s (#%d)",
                         sim.b[follow].name, follow);
            } else if (rc < 0) {
                forming = 0;
                snprintf(status, sizeof status, "spawn cancelled");
            }
        } else if (key == '?' && !forming) {
            help_on = !help_on;
            status[0] = '\0';
        } else if (key == STC_UI_KEY_ESC) {
            break;
        } else if (help_on && !forming) {
            /* Any waiting key returns from help; idle frames keep it. */
            if (key != STC_UI_KEY_NONE) {
                help_on = 0;
                status[0] = '\0';
            }
        } else if (key == 'P' || key == 'p') {
            paused = !paused;
        } else if (key == 'N' || key == 'n') {
            SbxFit fit = sbx_fit(&sim, follow);

            forming = 1;
            sbx_form_prefill(&form, &sim, fit.cx, fit.cy, fit.cz,
                             fit.spread);
        } else if (key == 'K' || key == 'k') {
            if (follow >= 0 && follow < sim.n) {
                snprintf(status, sizeof status, "removed %s",
                         sim.b[follow].name);
                stc_sim_remove(&sim, follow);
                sbx_trail_clear();
                if (follow >= sim.n) {
                    follow = sim.n - 1;
                }
            }
        } else if (key == 'G' || key == 'g') {
            preset = (preset + 1) % stc_preset_count();
            stc_sim_init(&sim, STC_G, 1000.0, 0.4, 1);
            stc_preset_load(&sim, preset, (uint64_t)seed);
            follow = -1;
            sbx_trail_clear();
            snprintf(status, sizeof status, "preset=%s",
                     stc_preset_name(preset));
        } else if (key == '\t') {
            follow = sim.n > 0 ? (follow + 2) % (sim.n + 1) - 1 : -1;
        } else if (key == 'Z' || key == 'z') {
            dt = dt / 2.0 < 1e-6 ? 1e-6 : dt / 2.0;
        } else if (key == 'X' || key == 'x') {
            dt = dt * 2.0 > 1e9 ? 1e9 : dt * 2.0;
        } else if (key == 'C' || key == 'c') {
            warp_i = (warp_i + 1) % 3;
        } else if (key == 'T' || key == 't') {
            g_trail_on = !g_trail_on;
            if (!g_trail_on) {
                sbx_trail_clear();
            }
        } else if (key == 'E' || key == 'e') {
            char cpath[1024];
            char ppath[1024];
            int okc = 0;
            int okp = 0;
            SbxFit fit = sbx_fit(&sim, follow);
            StcCamera rcam = sbx_camera(cam.yaw_rad, cam.pitch_rad,
                                       cam.dist_m, STC_APP_PPM_W,
                                       STC_APP_PPM_H);

            snprintf(cpath, sizeof cpath, "%s",
                     o->csv != NULL ? o->csv : "sandbox.csv");
            snprintf(ppath, sizeof ppath, "%s",
                     o->ppm != NULL ? o->ppm : "sandbox.ppm");
            okc = sbx_write_bodies_csv(cpath, &sim) == 0;
            okp = sbx_write_ppm(ppath, &sim, &fit, &rcam, follow) == 0;
            snprintf(status, sizeof status, "export csv=%s ppm=%s",
                     okc ? cpath : "FAIL", okp ? ppath : "FAIL");
        } else if (key != STC_UI_KEY_NONE) {
            stc_ui_camera_key(&cam, key);
        }

        if (!paused && !forming && !help_on) {
            int w;

            for (w = 0; w < kWarps[warp_i]; w++) {
                sbx_step(&sim, dt, integrator);
            }
            sbx_trail_push(&sim);
            f++;
            if (cap >= 0 && f >= cap) {
                break;
            }
        }
        sbx_pace(t0, o->no_adapt);
    }
    sbx_trail_free();
    stc_ui_exit_raw();
    printf("\nbye: t=%g s steps=%ld bodies=%d\n", sim.t, f, sim.n);
    return 0;
}

/* Headless: preset + N frames of single steps + time-series CSV + PPM. */
static int sbx_run_headless(const StcAppOpts *o, int preset,
                           unsigned long long seed, long nframes)
{
    StcSim sim;
    FILE *fp = NULL;
    const char *csv_use = o->csv != NULL ? o->csv : NULL;
    const char *ppm_use = o->ppm != NULL ? o->ppm : NULL;
    int integrator = (o != NULL) ? o->integrator : STC_APP_INTEG_VERLET;
    double dt = 1000.0;
    long f;

    stc_sim_init(&sim, STC_G, 1000.0, 0.4, 1);
    if (stc_preset_load(&sim, preset, (uint64_t)seed) != 0) {
        fprintf(stderr, "stlab: bad preset %d (see --help)\n", preset);
        return 2;
    }
    if (csv_use != NULL) {
        fp = fopen(csv_use, "w");
        if (fp == NULL) {
            fprintf(stderr, "stlab: cannot write csv '%s'\n", csv_use);
            return 1;
        }
        if (fprintf(fp, "frame,t_s,E_J,n\n") < 0) {
            fclose(fp);
            return 1;
        }
    }
    sbx_trail_init();
    for (f = 0; f < nframes; f++) {
        clock_t t0 = clock();

        sbx_step(&sim, dt, integrator);
        sbx_trail_push(&sim);
        if (fp != NULL) {
            if (fprintf(fp, "%ld,%.17g,%.17g,%d\n", f, sim.t,
                       stc_sim_energy(&sim), sim.n) < 0) {
                fclose(fp);
                sbx_trail_free();
                return 1;
            }
        }
        sbx_pace(t0, o->no_adapt);
    }
    if (fp != NULL && fclose(fp) != 0) {
        sbx_trail_free();
        return 1;
    }
    if (ppm_use != NULL) {
        SbxFit fit = sbx_fit(&sim, -1);
        StcCamera rcam = sbx_camera(0.6, 0.5, 3.2, STC_APP_PPM_W,
                                   STC_APP_PPM_H);

        if (sbx_write_ppm(ppm_use, &sim, &fit, &rcam, -1) != 0) {
            fprintf(stderr, "stlab: cannot write ppm '%s'\n", ppm_use);
            sbx_trail_free();
            return 1;
        }
    }
    printf("sandbox preset=%s integ=%s frames=%ld csv=%s ppm=%s E=%.6g J\n",
           stc_preset_name(preset), sbx_integ_name(integrator), nframes,
           csv_use != NULL ? csv_use : "-", ppm_use != NULL ? ppm_use : "-",
           stc_sim_energy(&sim));
    sbx_trail_free();
    return 0;
}

int stc_sandbox_run(const StcAppOpts *o)
{
    StcAppOpts def;
    int preset;
    unsigned long long seed;
    long nframes;

    if (o == NULL) {
        stc_app_opts_init(&def);
        o = &def;
    }
    preset = o->preset_given ? (int)o->preset : 0;
    seed = o->seed_given ? o->seed : 20260927ULL;
    nframes = o->frames_given ? o->frames : STC_APP_DEFAULT_FRAMES;
    if (o->smoke && !o->frames_given) {
        nframes = 60L;
    }
    if (preset < 0 || preset >= stc_preset_count()) {
        fprintf(stderr, "stlab: bad --preset %d (0..%d)\n", preset,
                stc_preset_count() - 1);
        return 2;
    }
    if (o->visual_shot != NULL) {
        long shot_frames = nframes;

        if (!o->frames_given) {
            shot_frames = 60L;
        }
        if (stc_vis_shot(preset, seed, o->integrator, shot_frames,
                         STC_VIS_DEF_W, STC_VIS_DEF_H,
                         o->visual_shot) != 0) {
            fprintf(stderr, "stlab: cannot write visual shot '%s'\n",
                    o->visual_shot);
            return 1;
        }
        printf("visual preset=%s frames=%ld ppm=%s\n",
               stc_preset_name(preset), shot_frames, o->visual_shot);
        return 0;
    }
    if (o->visual) {
        return stc_vis_run(o, preset, seed);
    }
    if (!o->smoke) {
        int rc = sbx_run_interactive(o, preset, seed);

        if (rc != -99) {
            return rc;
        }
    }
    return sbx_run_headless(o, preset, seed, nframes);
}
