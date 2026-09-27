/* SPACETIME LAB — app: CLI + frame loop + deterministic export.
 * Includes core + own + scenario + render only; physics is reached through
 * the scenario layer (setup/step per frame). Pure C11 here (POSIX termios
 * lives in ui.c only). SI doubles, seeded RNG.
 */

#include <errno.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "stc/vec.h"
#include "stc/units.h"
#include "stc/rng.h"
#include "stc/app.h"
#include "stc/sandbox.h"
#include "stc/scenario.h"
#include "stc/render.h"
#include "stc/vis_merger.h"
#include "stc/vis_space.h"

/* Equation: preset table — contract Sec 5 names in CLI/CSV/verify order. */
static const char *const kScenarioNames[STC_APP_NSCEN] = {
    "twin", "pole", "doppler", "mercury", "gps",
    "bend", "ring", "flamm", "chirp", "shadow"
};

/* Live simulation state: one scenario context, stepped per frame. When a
 * run outlasts the preset frame count the preset loops (re-setup), so any
 * --frames value stays deterministic. */
typedef struct {
    StcScenarioCtx sc;
    const char *name;
    int ok;
} StcAppSim;

/* Equation: sim setup — scenario preset initial state (fixed seed). */
static int app_sim_setup(StcAppSim *s, const char *name)
{
    if (s == NULL || name == NULL) {
        return -1;
    }
    s->name = name;
    s->ok = (stc_scenario_setup(&s->sc, name) == 0);
    return s->ok ? 0 : -1;
}

/* Equation: sim advance — one preset step; loops the preset at completion. */
static void app_sim_step(StcAppSim *s)
{
    int rc;

    if (s == NULL || !s->ok) {
        return;
    }
    rc = stc_scenario_step(&s->sc);
    if (rc > 0) {
        /* Preset complete: loop deterministically. */
        if (stc_scenario_setup(&s->sc, s->name) == 0) {
            (void)stc_scenario_step(&s->sc);
        } else {
            s->ok = 0;
        }
    } else if (rc < 0) {
        s->ok = 0;
    }
}

/* Normalized view fit: centroid + spread over pos[0..3], used to scale every
 * scenario (meters to light-years to strain) into a unit rendering box.
 * Labels and CSV stay in raw SI; only the projection is normalized. */
typedef struct {
    double cx;
    double cy;
    double cz;
    double spread;
} StcAppFit;

static StcAppFit app_fit(const StcScenarioCtx *c)
{
    StcAppFit f;
    int i;

    f.cx = 0.0;
    f.cy = 0.0;
    f.cz = 0.0;
    for (i = 0; i < STC_SCENARIO_NPOS; i++) {
        f.cx += c->pos[i].x;
        f.cy += c->pos[i].y;
        f.cz += c->pos[i].z;
    }
    f.cx /= (double)STC_SCENARIO_NPOS;
    f.cy /= (double)STC_SCENARIO_NPOS;
    f.cz /= (double)STC_SCENARIO_NPOS;
    f.spread = 0.0;
    for (i = 0; i < STC_SCENARIO_NPOS; i++) {
        double dx = c->pos[i].x - f.cx;
        double dy = c->pos[i].y - f.cy;
        double dz = c->pos[i].z - f.cz;
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

/* Equation: orbit eye — eye = dist*(cos p*cos y, sin p, cos p*sin y), target
 * origin, up +Y (same model as stc_ui_camera_update, unitless space). */
static StcCamera app_camera(const StcAppFit *fit, double yaw, double pitch,
                           double dist, int w, int h)
{
    StcVec3 eye;
    StcVec3 center;
    StcVec3 up;

    (void)fit;
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

/* Project the scenario points into pixels. Returns plotted count. */
static int app_project(const StcScenarioCtx *c, const StcAppFit *fit,
                      const StcCamera *cam, int *sx, int *sy, int n)
{
    int i;
    int plotted = 0;

    for (i = 0; i < n && i < STC_SCENARIO_NPOS; i++) {
        StcVec3 p;

        p.x = (c->pos[i].x - fit->cx) / fit->spread;
        p.y = (c->pos[i].y - fit->cy) / fit->spread;
        p.z = (c->pos[i].z - fit->cz) / fit->spread;
        if (stc_render_project_point(cam, p, &sx[i], &sy[i], NULL)) {
            plotted++;
        } else {
            sx[i] = -1000000;
            sy[i] = -1000000;
        }
    }
    return plotted;
}

/* Observable display labels per preset (first 4 slots; see scenario.h). */
static const char *const kObsLabels[STC_APP_NSCEN][4] = {
    { "earth_t_s", "ship_tau_s", "gamma", "v_m_s" },
    { "L_contr_m", "L0_m", "gamma", "barn_m" },
    { "z", "f_obs_Hz", "f_emit_Hz", "beta" },
    { "prec_arcsec", "orbits", "a_m", "e" },
    { "SR_us_d", "GR_us_d", "net_us_d", "accum_ns" },
    { "defl_arcsec", "b_m", "photon_x_m", "photon_y_m" },
    { "thetaE_as", "theta+_as", "theta-_as", "beta_as" },
    { "rs_m", "w_m", "r_m", "wmax_m" },
    { "f_Hz", "h", "h_plus", "h_cross" },
    { "chi", "isco_m", "rplus_m", "shadow_m" },
};

/* Equation: frame pacing — busy-wait on clock() until dt = 1/60 s elapses. */
static void app_pace(clock_t start, int no_adapt)
{
    double elapsed;

    if (no_adapt) {
        return;
    }
    if (start == (clock_t)-1) {
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

/* Equation: output path — single scenario uses `given` verbatim (default
 * "<scen>.<ext>"); "all" appends "<scen>.<ext>" under `given` as a directory
 * (parent dirs must already exist). Returns 0 ok, -1 on bad args/truncation. */
static int app_out_path(char *out, size_t cap, const char *given,
                        const char *scen, const char *ext, int multi)
{
    int r;

    if (out == NULL || cap == 0 || scen == NULL || ext == NULL) {
        return -1;
    }
    out[0] = '\0';
    if (!multi) {
        if (given != NULL && given[0] != '\0') {
            r = snprintf(out, cap, "%s", given);
        } else {
            r = snprintf(out, cap, "%s.%s", scen, ext);
        }
        if (r < 0 || (size_t)r >= cap) {
            out[0] = '\0';
            return -1;
        }
        return 0;
    }
    if (given != NULL && given[0] != '\0') {
        size_t len = strlen(given);

        if (len > 0 && given[len - 1] == '/') {
            r = snprintf(out, cap, "%s%s.%s", given, scen, ext);
        } else {
            r = snprintf(out, cap, "%s/%s.%s", given, scen, ext);
        }
    } else {
        r = snprintf(out, cap, "%s.%s", scen, ext);
    }
    if (r < 0 || (size_t)r >= cap) {
        out[0] = '\0';
        return -1;
    }
    return 0;
}

/* Equation: scenario PPM — project pos[0..3] through the orbit camera,
 * paint phosphor points + amber path + dim axes, write P6 via the renderer.
 * Byte-deterministic for the same scenario + frame. */
static int app_write_ppm(const char *path, const StcScenarioCtx *c,
                        double yaw, double pitch, double dist)
{
    StcRgbBuf fb;
    StcAppFit fit;
    StcCamera cam;
    StcRgb bg;
    StcRgb amber;
    StcRgb phos;
    StcRgb dim;
    int sx[STC_SCENARIO_NPOS];
    int sy[STC_SCENARIO_NPOS];
    int i;
    int rc;

    if (path == NULL || c == NULL) {
        return -1;
    }
    if (stc_render_rgb_init(&fb, STC_APP_PPM_W, STC_APP_PPM_H) != 0) {
        return -2;
    }
    bg = stc_render_rgb(4, 8, 12);
    amber = stc_render_amber(0.9);
    phos = stc_render_phosphor(1.0);
    dim = stc_render_rgb(28, 40, 52);
    stc_render_rgb_clear(&fb, bg);
    fit = app_fit(c);
    cam = app_camera(&fit, yaw, pitch, dist, fb.w, fb.h);
    app_project(c, &fit, &cam, sx, sy, STC_SCENARIO_NPOS);
    /* Axes cross at the projected centroid (centroid maps near center). */
    stc_render_rgb_line(&fb, fb.w / 2 - 40, fb.h / 2, fb.w / 2 + 40, fb.h / 2,
                       dim);
    stc_render_rgb_line(&fb, fb.w / 2, fb.h / 2 - 25, fb.w / 2, fb.h / 2 + 25,
                       dim);
    for (i = 0; i + 1 < STC_SCENARIO_NPOS; i++) {
        if (sx[i] > -1000000 && sx[i + 1] > -1000000 &&
            (sx[i] != sx[i + 1] || sy[i] != sy[i + 1])) {
            stc_render_rgb_line(&fb, sx[i], sy[i], sx[i + 1], sy[i + 1],
                               amber);
        }
    }
    for (i = 0; i < STC_SCENARIO_NPOS; i++) {
        int dx;
        int dy;

        if (sx[i] <= -1000000) {
            continue;
        }
        for (dy = -1; dy <= 1; dy++) {
            for (dx = -1; dx <= 1; dx++) {
                stc_render_rgb_plot(&fb, sx[i] + dx, sy[i] + dy,
                                   i == 0 ? phos : amber);
            }
        }
    }
    rc = stc_render_write_ppm(&fb, path);
    stc_render_rgb_free(&fb);
    return rc == 0 ? 0 : -2;
}

/* Equation: scenario ASCII — same projection as the PPM painter, glyphs +
 * ANSI color, presented to `out`. Returns 0 on success. */
static int app_present_ascii(const StcScenarioCtx *c, double yaw, double pitch,
                            double dist, FILE *out)
{
    StcAsciiBuf fb;
    StcAppFit fit;
    StcCamera cam;
    int sx[STC_SCENARIO_NPOS];
    int sy[STC_SCENARIO_NPOS];
    int i;
    int rc;

    if (c == NULL || out == NULL) {
        return -1;
    }
    if (stc_render_ascii_init(&fb, STC_APP_ASCII_W, STC_APP_ASCII_H) != 0) {
        return -2;
    }
    stc_render_ascii_clear(&fb, ' ', STC_RENDER_DEFAULT);
    fit = app_fit(c);
    cam = app_camera(&fit, yaw, pitch, dist, fb.w, fb.h);
    app_project(c, &fit, &cam, sx, sy, STC_SCENARIO_NPOS);
    stc_render_ascii_line(&fb, fb.w / 2 - 20, fb.h / 2, fb.w / 2 + 20,
                         fb.h / 2, '-', STC_RENDER_BRIGHT_BLACK);
    stc_render_ascii_line(&fb, fb.w / 2, 1, fb.w / 2, fb.h - 2, '|',
                         STC_RENDER_BRIGHT_BLACK);
    for (i = 0; i + 1 < STC_SCENARIO_NPOS; i++) {
        if (sx[i] > -1000000 && sx[i + 1] > -1000000 &&
            (sx[i] != sx[i + 1] || sy[i] != sy[i + 1])) {
            stc_render_ascii_line(&fb, sx[i], sy[i], sx[i + 1], sy[i + 1],
                                 '.', STC_RENDER_YELLOW);
        }
    }
    for (i = 0; i < STC_SCENARIO_NPOS; i++) {
        if (sx[i] <= -1000000) {
            continue;
        }
        stc_render_ascii_plot(&fb, sx[i], sy[i], i == 0 ? '@' : 'o',
                             i == 0 ? STC_RENDER_BRIGHT_GREEN
                                    : STC_RENDER_BRIGHT_YELLOW);
    }
    rc = stc_render_ascii_present(&fb, out);
    stc_render_ascii_free(&fb);
    return rc;
}

/* Equation: CSV row — frame, sim time, head position, head observable,
 * then all 8 observable slots (SI; see scenario.h per-preset table). */
static int app_csv_row(StcCsv *csv, long f, const StcScenarioCtx *c)
{
    double row[6 + STC_SCENARIO_NOBS];
    int i;

    row[0] = (double)f;
    row[1] = c->t;
    row[2] = c->pos[0].x;
    row[3] = c->pos[0].y;
    row[4] = c->pos[0].z;
    row[5] = c->obs[0];
    for (i = 0; i < STC_SCENARIO_NOBS; i++) {
        row[6 + i] = c->obs[i];
    }
    return stc_csv_row(csv, row, 6 + STC_SCENARIO_NOBS);
}

#define STC_APP_CSV_HEADER \
    "frame,t_s,x_m,y_m,z_m,q,obs0,obs1,obs2,obs3,obs4,obs5,obs6,obs7"

/* Default orbit-camera angle for headless rendering (radians, unitless). */
#define STC_APP_DEF_YAW 0.6
#define STC_APP_DEF_PITCH 0.5
#define STC_APP_DEF_DIST 3.2

/* Equation: headless scenario run — N paced frames + optional CSV/PPM export. */
static int app_run_headless_one(const StcAppOpts *o, int idx, long nframes,
                                int multi)
{
    const char *scen = stc_app_scenario_name(idx);
    StcAppSim sim;
    char csv_path[1024];
    char ppm_path[1024];
    const char *csv_use = NULL;
    const char *ppm_use = NULL;
    StcCsv csv;
    long f;

    sim.ok = 0;
    sim.name = NULL;
    csv.fp = NULL;
    csv.ncols = -1;
    csv.nrows = 0;
    csv.open = 0;
    csv_path[0] = '\0';
    ppm_path[0] = '\0';

    if (app_sim_setup(&sim, scen) != 0) {
        fprintf(stderr, "stlab: scenario setup failed '%s'\n", scen);
        return 1;
    }

    if (o->csv != NULL) {
        if (app_out_path(csv_path, sizeof csv_path, o->csv, scen, "csv",
                         multi) != 0) {
            fprintf(stderr, "stlab: csv path too long\n");
            return 1;
        }
        csv_use = csv_path;
        if (stc_csv_open(&csv, csv_use, STC_APP_CSV_HEADER) != 0) {
            fprintf(stderr, "stlab: cannot write csv '%s'\n", csv_use);
            return 1;
        }
    }
    if (o->ppm != NULL) {
        if (app_out_path(ppm_path, sizeof ppm_path, o->ppm, scen, "ppm",
                         multi) != 0) {
            fprintf(stderr, "stlab: ppm path too long\n");
            stc_csv_close(&csv);
            return 1;
        }
        ppm_use = ppm_path;
    }
    for (f = 0; f < nframes; f++) {
        clock_t t0 = clock();

        app_sim_step(&sim);
        if (!sim.ok) {
            fprintf(stderr, "stlab: scenario stepped into failure '%s'\n",
                    scen);
            stc_csv_close(&csv);
            return 1;
        }
        if (csv_use != NULL) {
            if (app_csv_row(&csv, f, &sim.sc) != 0) {
                fprintf(stderr, "stlab: csv write failed '%s'\n", csv_use);
                stc_csv_close(&csv);
                return 1;
            }
        }
        app_pace(t0, o->no_adapt);
    }
    stc_csv_close(&csv);
    if (ppm_use != NULL) {
        if (app_write_ppm(ppm_use, &sim.sc, STC_APP_DEF_YAW,
                         STC_APP_DEF_PITCH, STC_APP_DEF_DIST) != 0) {
            fprintf(stderr, "stlab: cannot write ppm '%s'\n", ppm_use);
            return 1;
        }
    }
    printf("scenario=%s frames=%ld csv=%s ppm=%s\n", scen, nframes,
           csv_use != NULL ? csv_use : "-",
           ppm_use != NULL ? ppm_use : "-");
    return 0;
}

/* Equation: batch dispatch — selected set headless (bare --smoke runs all). */
static int app_run_batch(const StcAppOpts *o)
{
    long nframes;
    int scen_idx;
    int i;

    if (o->smoke && !o->scenario_given) {
        nframes = o->frames_given ? o->frames : STC_APP_SMOKE_FRAMES;
        for (i = 0; i < STC_APP_NSCEN; i++) {
            int rc = app_run_headless_one(o, i, nframes, 1);

            if (rc != 0) {
                return rc;
            }
        }
        return 0;
    }
    scen_idx = stc_app_scenario_index(o->scenario);
    if (scen_idx == -2) {
        nframes = o->frames_given ? o->frames : STC_APP_DEFAULT_FRAMES;
        if (o->smoke && !o->frames_given) {
            nframes = STC_APP_SMOKE_FRAMES;
        }
        for (i = 0; i < STC_APP_NSCEN; i++) {
            int rc = app_run_headless_one(o, i, nframes, 1);

            if (rc != 0) {
                return rc;
            }
        }
        return 0;
    }
    if (scen_idx < 0) {
        fprintf(stderr, "stlab: unknown scenario '%s' (see --list)\n",
                o->scenario != NULL ? o->scenario : "(null)");
        return 2;
    }
    nframes = o->frames_given ? o->frames : STC_APP_DEFAULT_FRAMES;
    if (o->smoke && !o->frames_given) {
        nframes = STC_APP_SMOKE_FRAMES;
    }
    return app_run_headless_one(o, scen_idx, nframes, 0);
}

/* Equation: export current frame — one CSV row + one PPM via E key. */
static void app_export_current(const StcAppOpts *o, int cur, long f,
                               const StcScenarioCtx *c,
                               const StcUiCamera *cam, char *status,
                               size_t cap)
{
    const char *scen = stc_app_scenario_name(cur);
    char cpath[1024];
    char ppath[1024];
    int ok_csv = 0;
    int ok_ppm = 0;

    cpath[0] = '\0';
    ppath[0] = '\0';
    if (app_out_path(cpath, sizeof cpath, o->csv, scen, "csv", 0) == 0) {
        StcCsv ex;

        ex.fp = NULL;
        ex.ncols = -1;
        ex.nrows = 0;
        ex.open = 0;
        if (stc_csv_open(&ex, cpath, STC_APP_CSV_HEADER) == 0 &&
            app_csv_row(&ex, f, c) == 0) {
            ok_csv = 1;
        }
        stc_csv_close(&ex);
    }
    if (app_out_path(ppath, sizeof ppath, o->ppm, scen, "ppm", 0) == 0 &&
        app_write_ppm(ppath, c, cam->yaw_rad, cam->pitch_rad,
                     cam->dist_m) == 0) {
        ok_ppm = 1;
    }
    snprintf(status, cap, "export frame %ld: csv=%s ppm=%s", f,
             ok_csv ? cpath : "FAIL", ok_ppm ? ppath : "FAIL");
}

/* Equation: interactive loop — raw keys drive camera/scenario/pause/export. */
static int app_run_interactive(const StcAppOpts *o)
{
    StcUiCamera cam;
    StcAppSim sim;
    long cap = o->frames_given ? o->frames : -1L;
    char status[256];
    int cur;
    int cur_idx;
    int paused = 0;
    long f = 0;

    sim.ok = 0;
    sim.name = NULL;
    if (o->scenario != NULL && stc_app_scenario_index(o->scenario) == -1) {
        fprintf(stderr, "stlab: unknown scenario '%s' (see --list)\n",
                o->scenario);
        return 2;
    }
    cur_idx = stc_app_scenario_index(o->scenario);
    cur = cur_idx >= 0 ? cur_idx : 0; /* "all" starts on the first preset */
    if (stc_ui_enter_raw() != 0) {
        /* Not a TTY (or termios failed): fall back to headless batch. */
        return app_run_batch(o);
    }
    atexit(stc_ui_exit_raw);
    stc_ui_camera_init(&cam);
    cam.yaw_rad = STC_APP_DEF_YAW;
    cam.pitch_rad = STC_APP_DEF_PITCH;
    cam.dist_m = STC_APP_DEF_DIST;
    stc_ui_camera_update(&cam);
    if (app_sim_setup(&sim, stc_app_scenario_name(cur)) != 0) {
        fprintf(stderr, "stlab: scenario setup failed\n");
        stc_ui_exit_raw();
        return 1;
    }
    status[0] = '\0';

    for (;;) {
        clock_t t0 = clock();
        int key;

        printf("\x1b[2J\x1b[H");
        printf("stlab  scenario=%s (%d/%d)  frame=%ld  t=%.6gs%s\n",
               stc_app_scenario_name(cur), cur + 1, STC_APP_NSCEN, f,
               sim.sc.t, paused ? "  [PAUSED]" : "");
        printf("cam yaw=%.2f pitch=%.2f dist=%.2f (unitless view)\n",
               cam.yaw_rad, cam.pitch_rad, cam.dist_m);
        printf("obs %s=%.6g %s=%.6g | %s=%.6g %s=%.6g (SI)\n",
               kObsLabels[cur][0], sim.sc.obs[0],
               kObsLabels[cur][1], sim.sc.obs[1],
               kObsLabels[cur][2], sim.sc.obs[2],
               kObsLabels[cur][3], sim.sc.obs[3]);
        (void)app_present_ascii(&sim.sc, cam.yaw_rad, cam.pitch_rad,
                               cam.dist_m, stdout);
        if (status[0] != '\0') {
            printf("last: %s\n", status);
        }
        printf("keys: WASD/arrows camera  1-9,0 scenario  P pause  E export  Q quit\n");
        fflush(stdout);

        key = stc_ui_poll_key();
        if (key == STC_UI_KEY_ERROR) {
            snprintf(status, sizeof status, "input error");
        } else if (key == 'Q' || key == 'q' ||
                   key == STC_UI_KEY_CTRL_C || key == STC_UI_KEY_ESC) {
            break;
        } else if (key == 'P' || key == 'p') {
            paused = !paused;
        } else if (key == 'E' || key == 'e') {
            app_export_current(o, cur, f, &sim.sc, &cam, status,
                               sizeof status);
        } else if (key >= '1' && key <= '9') {
            cur = key - '1';
            f = 0;
            if (app_sim_setup(&sim, stc_app_scenario_name(cur)) != 0) {
                snprintf(status, sizeof status, "setup failed");
            } else {
                snprintf(status, sizeof status, "scenario=%s",
                         stc_app_scenario_name(cur));
            }
        } else if (key == '0') {
            cur = 9;
            f = 0;
            if (app_sim_setup(&sim, stc_app_scenario_name(cur)) != 0) {
                snprintf(status, sizeof status, "setup failed");
            } else {
                snprintf(status, sizeof status, "scenario=%s",
                         stc_app_scenario_name(cur));
            }
        } else if (key != STC_UI_KEY_NONE) {
            stc_ui_camera_key(&cam, key);
            stc_ui_camera_update(&cam);
        }

        if (!paused) {
            app_sim_step(&sim);
            f++;
            if (cap >= 0 && f >= cap) {
                break;
            }
        }
        app_pace(t0, o->no_adapt);
    }
    stc_ui_exit_raw();
    printf("\nbye: scenario=%s frames=%ld\n", stc_app_scenario_name(cur), f);
    return 0;
}

/* Equation: --movie-every bound — strtol with full-consumption, v >= 0. */
static int app_parse_nonneg(const char *val, long *out)
{
    char *end = NULL;
    long v;

    if (val == NULL || *val == '\0' || out == NULL) {
        return -1;
    }
    errno = 0;
    v = strtol(val, &end, 10);
    if (errno != 0 || end == val || *end != '\0' || v < 0) {
        return -1;
    }
    *out = v;
    return 0;
}

/* Equation: --frames bound — strtol with full-consumption + 0..MAX check. */
static int app_parse_frames(const char *val, long *out)
{
    char *end = NULL;
    long v;

    if (val == NULL || *val == '\0' || out == NULL) {
        return -1;
    }
    errno = 0;
    v = strtol(val, &end, 10);
    if (errno != 0 || end == val || *end != '\0') {
        return -1;
    }
    if (v < 0 || v > STC_APP_MAX_FRAMES) {
        return -1;
    }
    *out = v;
    return 0;
}

/* Equation: value-flag match — "--flag v" or "--flag=v"; 1 matched (sets *val),
 * 0 no match, -1 matched but value missing. */
static int app_take_value(const char *arg, const char *flag, const char **val,
                          int argc, char *argv[], int *idx)
{
    size_t n = strlen(flag);

    if (strncmp(arg, flag, n) != 0) {
        return 0;
    }
    if (arg[n] == '=') {
        *val = arg + n + 1;
        return 1;
    }
    if (arg[n] == '\0') {
        if (*idx + 1 >= argc || argv[*idx + 1] == NULL) {
            return -1;
        }
        (*idx)++;
        *val = argv[*idx];
        return 1;
    }
    return 0;
}

void stc_app_opts_init(StcAppOpts *o)
{
    if (o == NULL) {
        return;
    }
    o->scenario = "twin";
    o->scenario_given = 0;
    o->frames = STC_APP_DEFAULT_FRAMES;
    o->frames_given = 0;
    o->ppm = NULL;
    o->csv = NULL;
    o->smoke = 0;
    o->list = 0;
    o->help = 0;
    o->no_adapt = 0;
    o->sandbox = 0;
    o->preset = 0;
    o->preset_given = 0;
    o->seed = 20260927ULL;
    o->seed_given = 0;
    o->integrator = STC_APP_INTEG_VERLET;
    o->integrator_given = 0;
    o->visual = 0;
    o->visual_shot = NULL;
    o->merger = 0;
    o->movie = NULL;
    o->movie_every = 2L;
    o->space = 0;
}

int stc_app_parse_args(StcAppOpts *o, int argc, char *argv[])
{
    int i;

    if (o == NULL || argv == NULL || argc <= 0) {
        return 2;
    }
    for (i = 1; i < argc; i++) {
        const char *a = argv[i];

        if (a == NULL) {
            fprintf(stderr, "stlab: null argument\n");
            return 2;
        }
        if (strcmp(a, "--smoke") == 0) {
            o->smoke = 1;
        } else if (strcmp(a, "--sandbox") == 0) {
            o->sandbox = 1;
        } else if (strcmp(a, "--list") == 0) {
            o->list = 1;
        } else if (strcmp(a, "--help") == 0 || strcmp(a, "-h") == 0) {
            o->help = 1;
        } else if (strcmp(a, "--no-adapt") == 0) {
            o->no_adapt = 1;
        } else if (strcmp(a, "--visual") == 0) {
            o->visual = 1;
        } else if (strcmp(a, "--merger") == 0) {
            o->merger = 1;
        } else if (strcmp(a, "--space") == 0) {
            o->space = 1;
        } else {
            const char *val = NULL;
            int m;

            m = app_take_value(a, "--scenario", &val, argc, argv, &i);
            if (m < 0) {
                fprintf(stderr, "stlab: --scenario needs a value\n");
                return 2;
            }
            if (m > 0) {
                o->scenario = val;
                o->scenario_given = 1;
                continue;
            }
            m = app_take_value(a, "--frames", &val, argc, argv, &i);
            if (m < 0) {
                fprintf(stderr, "stlab: --frames needs a value\n");
                return 2;
            }
            if (m > 0) {
                if (app_parse_frames(val, &o->frames) != 0) {
                    fprintf(stderr, "stlab: bad --frames '%s' (0..%ld)\n",
                            val, STC_APP_MAX_FRAMES);
                    return 2;
                }
                o->frames_given = 1;
                continue;
            }
            m = app_take_value(a, "--preset", &val, argc, argv, &i);
            if (m < 0) {
                fprintf(stderr, "stlab: --preset needs a value\n");
                return 2;
            }
            if (m > 0) {
                if (app_parse_frames(val, &o->preset) != 0) {
                    fprintf(stderr, "stlab: bad --preset '%s'\n", val);
                    return 2;
                }
                o->preset_given = 1;
                continue;
            }
            m = app_take_value(a, "--seed", &val, argc, argv, &i);
            if (m < 0) {
                fprintf(stderr, "stlab: --seed needs a value\n");
                return 2;
            }
            if (m > 0) {
                char *end = NULL;
                unsigned long long v;

                errno = 0;
                v = strtoull(val, &end, 10);
                if (errno != 0 || end == val || *end != '\0') {
                    fprintf(stderr, "stlab: bad --seed '%s'\n", val);
                    return 2;
                }
                o->seed = v;
                o->seed_given = 1;
                continue;
            }
            m = app_take_value(a, "--integrator", &val, argc, argv, &i);
            if (m < 0) {
                fprintf(stderr, "stlab: --integrator needs a value\n");
                return 2;
            }
            if (m > 0) {
                if (strcmp(val, "verlet") == 0) {
                    o->integrator = STC_APP_INTEG_VERLET;
                } else if (strcmp(val, "rk4") == 0) {
                    o->integrator = STC_APP_INTEG_RK4;
                } else {
                    fprintf(stderr,
                            "stlab: bad --integrator '%s' (verlet|rk4)\n",
                            val);
                    return 2;
                }
                o->integrator_given = 1;
                continue;
            }
            m = app_take_value(a, "--ppm", &val, argc, argv, &i);
            if (m < 0) {
                fprintf(stderr, "stlab: --ppm needs a value\n");
                return 2;
            }
            if (m > 0) {
                o->ppm = val;
                continue;
            }
            m = app_take_value(a, "--csv", &val, argc, argv, &i);
            if (m < 0) {
                fprintf(stderr, "stlab: --csv needs a value\n");
                return 2;
            }
            if (m > 0) {
                o->csv = val;
                continue;
            }
            m = app_take_value(a, "--visual-shot", &val, argc, argv, &i);
            if (m < 0) {
                fprintf(stderr, "stlab: --visual-shot needs a value\n");
                return 2;
            }
            if (m > 0) {
                o->visual_shot = val;
                continue;
            }
            m = app_take_value(a, "--movie-every", &val, argc, argv, &i);
            if (m < 0) {
                fprintf(stderr, "stlab: --movie-every needs a value\n");
                return 2;
            }
            if (m > 0) {
                if (app_parse_nonneg(val, &o->movie_every) != 0 ||
                    o->movie_every <= 0) {
                    fprintf(stderr,
                            "stlab: bad --movie-every '%s' (must be >= 1)\n",
                            val);
                    return 2;
                }
                continue;
            }
            m = app_take_value(a, "--movie", &val, argc, argv, &i);
            if (m < 0) {
                fprintf(stderr, "stlab: --movie needs a value\n");
                return 2;
            }
            if (m > 0) {
                o->movie = val;
                continue;
            }
            fprintf(stderr, "stlab: unknown flag '%s' (see --help)\n", a);
            return 2;
        }
    }
    return 0;
}

void stc_app_print_help(const char *prog)
{
    const char *p = (prog != NULL && prog[0] != '\0') ? prog : "stlab";
    int k;

    printf("usage: %s [options]\n", p);
    printf("\nSPACETIME LAB relativity demos (10 scenarios, 60fps frame loop).\n\n");
    printf("options:\n");
    printf("  --scenario <name|all>  preset: twin pole doppler mercury gps bend ring\n");
    printf("                         flamm chirp shadow, or all (default: twin;\n");
    printf("                         bare --smoke runs all)\n");
    printf("  --frames N             frames per scenario, 0..1000000 (default: 600;\n");
    printf("                         smoke default: 1)\n");
    printf("  --ppm <path>           PPM output: file for one scenario, directory\n");
    printf("                         when scenario=all (parent dirs must exist)\n");
    printf("  --csv <path>           CSV output: same file-or-directory rule as --ppm\n");
    printf("  --sandbox              N-body sandbox mode (spawn bodies, time warp)\n");
    printf("  --preset N             sandbox preset 0..3 (solar, binary, hole_disk,\n");
    printf("                         slingshot; default 0)\n");
    printf("  --seed N               sandbox RNG seed (default 20260927)\n");
    printf("  --integrator verlet|rk4\n");
    printf("                         sandbox stepper (default verlet; both run\n");
    printf("                         the same collision/cap poststep)\n");
    printf("  --visual               native-window visual sandbox (macOS only)\n");
    printf("  --visual-shot <ppm>    headless visual render to PPM, no window\n");
    printf("  --merger               cinematic binary black hole merger view,\n");
    printf("                         presets 0..3: GW150914-like, equal mass,\n");
    printf("                         q=3, q=8; opens a window by default\n");
    printf("                         (--visual is accepted but not required);\n");
    printf("                         --visual-shot/--movie run headless instead\n");
    printf("  --space                full space simulator: real Solar System from\n");
    printf("                         JPL Horizons (J2000 / today), N-body + 1PN,\n");
    printf("                         time warp both ways, black-hole and swarm\n");
    printf("                         sandbox; presets 0..5 (--preset selects the\n");
    printf("                         scenario); opens a window by default\n");
    printf("                         (--visual is accepted but not required);\n");
    printf("                         --visual-shot/--movie run headless instead:\n");
    for (k = 0; k < STC_SVIS_NSCEN; k++) {
        printf("                           %d: %s\n", k,
               stc_svis_scenario_name(k));
    }
    printf("  --movie <dir>          --merger/--space headless: write a PPM frame\n");
    printf("                         sequence to dir (parent dirs must exist)\n");
    printf("  --movie-every N        --movie frame stride (default 2)\n");
    printf("  --smoke                headless run, no UI, exit 0 on success\n");
    printf("  --list                 print the 10 preset names and exit 0\n");
    printf("  --help                 print this help and exit 0\n");
    printf("  --no-adapt             disable 60fps pacing (run as fast as possible)\n");
    printf("\ninteractive keys (TTY): WASD/arrows camera, +/- zoom, 1-9,0 scenario,\n");
    printf("  P pause, E export current frame, Q quit (Ctrl-C/Esc also quit).\n");
    printf("outputs: CSV header frame,t_s,x_m,y_m,z_m,q,obs0..obs7 (SI);\n");
    printf("  PPM is P6 binary %dx%d. Same argv = byte-identical files.\n",
           STC_APP_PPM_W, STC_APP_PPM_H);
}

void stc_app_print_list(void)
{
    int i;

    for (i = 0; i < STC_APP_NSCEN; i++) {
        printf("%s\n", kScenarioNames[i]);
    }
}

int stc_app_run(const StcAppOpts *o)
{
    StcAppOpts def;

    if (o == NULL) {
        stc_app_opts_init(&def);
        o = &def;
    }
    if (o->help) {
        stc_app_print_help("stlab");
        return 0;
    }
    if (o->list) {
        stc_app_print_list();
        return 0;
    }
    if (o->merger) {
        int preset = o->preset_given ? (int)o->preset : 0;

        if (preset < 0 || preset >= STC_MVIS_NPRESET) {
            fprintf(stderr, "stlab: bad --preset %d (0..%d)\n", preset,
                    STC_MVIS_NPRESET - 1);
            return 2;
        }
        if (o->visual_shot != NULL || o->movie != NULL) {
            long frames = o->frames_given ? o->frames : 600L;
            int rc = stc_mvis_shot(preset, frames, 1280, 800,
                                  o->visual_shot, o->movie,
                                  o->movie != NULL ? o->movie_every : 0L);

            return rc != 0 ? 1 : 0;
        }
        return stc_mvis_run(o, preset);
    }
    if (o->space) {
        int scenario = o->preset_given ? (int)o->preset : 0;

        if (scenario < 0 || scenario >= STC_SVIS_NSCEN) {
            fprintf(stderr, "stlab: bad --preset %d (0..%d)\n", scenario,
                    STC_SVIS_NSCEN - 1);
            return 2;
        }
        if (o->visual_shot != NULL || o->movie != NULL) {
            long frames = o->frames_given ? o->frames : 600L;
            int rc = stc_svis_shot(scenario, frames, 1280, 800,
                                   o->visual_shot, o->movie,
                                   o->movie != NULL ? o->movie_every : 0L);

            return rc != 0 ? 1 : 0;
        }
        return stc_svis_run(o, scenario);
    }
    if (o->sandbox) {
        return stc_sandbox_run(o);
    }
    if (o->smoke) {
        return app_run_batch(o);
    }
    return app_run_interactive(o);
}

int stc_app_scenario_count(void)
{
    return STC_APP_NSCEN;
}

const char *stc_app_scenario_name(int idx)
{
    if (idx < 0 || idx >= STC_APP_NSCEN) {
        return NULL;
    }
    return kScenarioNames[idx];
}

int stc_app_scenario_index(const char *name)
{
    int i;

    if (name == NULL) {
        return -1;
    }
    if (strcmp(name, "all") == 0) {
        return -2;
    }
    for (i = 0; i < STC_APP_NSCEN; i++) {
        if (strcmp(name, kScenarioNames[i]) == 0) {
            return i;
        }
    }
    return -1;
}

int main(int argc, char *argv[])
{
    StcAppOpts o;

    stc_app_opts_init(&o);
    if (stc_app_parse_args(&o, argc, argv) != 0) {
        return 2;
    }
    return stc_app_run(&o);
}
