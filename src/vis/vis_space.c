/* SPACETIME LAB — vis_space: space simulator front end (see
 * stc/vis_space.h). JPL mission-operations display conventions: black
 * field, thin grey orbit traces, monospace labels, tick-marked reference
 * rings, framed readout panels, no glow. Rendering reads the simulation
 * state only (R6); the command console edits the simulation through the
 * stc/space.h API.
 */

#include "stc/vis_space.h"

#include <ctype.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "stc/ephem.h"
#include "stc/par.h"
#include "stc/ui.h"
#include "stc/vis.h"

#define SV_AU 149597870700.0
#define SV_DAY 86400.0
#define SV_YEAR (365.25 * 86400.0)
#define SV_C 299792458.0
#define SV_G 6.67430e-11
#define SV_GM_SUN 1.32712440041e20
#define SV_GM_EARTH 3.986004418e14
#define SV_R_EARTH 6.371e6
#define SV_FOV 0.70
#define SV_PI 3.14159265358979323846

/* --- scenarios -------------------------------------------------------- */

static const char *kScen[STC_SVIS_NSCEN] = {
    "SOLAR SYSTEM  J2000.0",
    "SOLAR SYSTEM  2026-09-27",
    "SUN REPLACED BY 1 MSUN BLACK HOLE",
    "ROGUE 10 MSUN BLACK HOLE FLYBY",
    "EARTH SWARM  2000 EARTHS AT 1 AU",
    "BLACK HOLE + 1,000,000 TEST EARTHS",
};

const char *stc_svis_scenario_name(int k)
{
    return k >= 0 && k < STC_SVIS_NSCEN ? kScen[k] : NULL;
}

/* splitmix64 for sandbox spawns (seeded, deterministic). */
static uint64_t sv_mix(uint64_t *s)
{
    uint64_t z = (*s += 0x9e3779b97f4a7c15ULL);

    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
    return z ^ (z >> 31);
}

static double sv_u01(uint64_t *s)
{
    return ((double)(sv_mix(s) >> 11) + 0.5) / 9007199254740992.0;
}

static int sv_find(const StcSpace *s, const char *name)
{
    int i;

    for (i = 0; i < s->nm; i++) {
        const char *a = s->name[i];
        const char *b = name;

        if (!s->alive[i]) {
            continue;
        }
        while (*a != '\0' && *b != '\0' &&
               tolower((unsigned char)*a) == tolower((unsigned char)*b)) {
            a++;
            b++;
        }
        if (*a == '\0' && *b == '\0') {
            return i;
        }
    }
    return -1;
}

static int sv_heaviest(const StcSpace *s)
{
    int i;
    int best = -1;

    for (i = 0; i < s->nm; i++) {
        if (s->alive[i] && (best < 0 || s->gm[i] > s->gm[best])) {
            best = i;
        }
    }
    return best;
}

/* Kepler elements about a primary -> state (m, m/s). */
static void sv_kep(double mu, double a, double e, double inc, double node,
                   double argp, double M, double *p, double *v)
{
    double E = M;
    double cn = cos(node);
    double sn = sin(node);
    double ci = cos(inc);
    double si = sin(inc);
    double cw = cos(argp);
    double sw = sin(argp);
    double xo;
    double yo;
    double vxo;
    double vyo;
    double n;
    double r;
    int k;

    for (k = 0; k < 40; k++) {
        double d = (E - e * sin(E) - M) / (1.0 - e * cos(E));

        E -= d;
        if (fabs(d) < 1e-14) {
            break;
        }
    }
    xo = a * (cos(E) - e);
    yo = a * sqrt(1.0 - e * e) * sin(E);
    n = sqrt(mu / (a * a * a));
    r = a * (1.0 - e * cos(E));
    vxo = -a * n * sin(E) / (r / a);
    vyo = a * n * sqrt(1.0 - e * e) * cos(E) / (r / a);
    {
        double P11 = cn * cw - sn * sw * ci;
        double P12 = -cn * sw - sn * cw * ci;
        double P21 = sn * cw + cn * sw * ci;
        double P22 = -sn * sw + cn * cw * ci;
        double P31 = sw * si;
        double P32 = cw * si;

        p[0] = P11 * xo + P12 * yo;
        p[1] = P21 * xo + P22 * yo;
        p[2] = P31 * xo + P32 * yo;
        v[0] = P11 * vxo + P12 * vyo;
        v[1] = P21 * vxo + P22 * vyo;
        v[2] = P31 * vxo + P32 * vyo;
    }
}

/* Per-frame substep budget from the force cost (deterministic: depends
 * only on particle counts, never on wall-clock time). */
static long sv_budget(const StcSpace *s)
{
    double nm = stc_space_live_massive(s);
    double nt = stc_space_live_tracers(s);
    double cost;

    if (nm > s->tree_threshold) {
        cost = 7.0 * (nm + nt) * 40.0 * log2(nm + 2.0);
    } else {
        cost = 7.0 * (nm * nm + nt * nm) + 1.0;
    }
    /* calibrated on the M5 Pro build: ~1.1e-6 ms per cost unit, so
     * 1.0e7 units ~ 11 ms of integration per displayed frame */
    cost = 1.0e7 / cost;
    if (cost < 1.0) {
        cost = 1.0;
    }
    if (cost > 200000.0) {
        cost = 200000.0;
    }
    return (long)cost;
}

/* --- sandbox operations (console + scenarios) ------------------------- */

typedef enum { SV_RING, SV_DISK, SV_CLOUD, SV_SHELL } SvShape;

/* Spawn n Earth-like bodies about primary c on near-circular orbits
 * (ring: r0 +- 10%, i < 2 deg; disk: 0.3..1.5 r0 flat; cloud: isotropic
 * inclinations; shell: r0 isotropic). massive != 0 -> real Earths
 * (GM_Earth, R_Earth), else massless tracers. ILLUSTRATIVE population. */
static int sv_spawn(StcSpace *s, int c, long n, SvShape shape, double r0,
                    int massive, uint64_t seed, char *msg, int msgn)
{
    long k;
    long made = 0;
    double mu;
    uint64_t st = seed;

    if (c < 0 || n <= 0 || !(r0 > 0.0)) {
        snprintf(msg, (size_t)msgn, "SPAWN: BAD ARGUMENTS");
        return -1;
    }
    if (massive) {
        if (stc_space_reserve(s, s->nm + (int)n, s->cap_t) != 0) {
            snprintf(msg, (size_t)msgn, "SPAWN: OUT OF MEMORY");
            return -1;
        }
    } else if (stc_space_reserve(s, s->cap_m, s->nt + (int)n) != 0) {
        snprintf(msg, (size_t)msgn, "SPAWN: OUT OF MEMORY");
        return -1;
    }
    mu = s->gm[c];
    for (k = 0; k < n; k++) {
        double a;
        double inc;
        double e = 0.02 * sv_u01(&st);
        double node = 2.0 * SV_PI * sv_u01(&st);
        double argp = 2.0 * SV_PI * sv_u01(&st);
        double M = 2.0 * SV_PI * sv_u01(&st);
        double p[3];
        double v[3];
        int j;

        switch (shape) {
        case SV_DISK:
            a = r0 * (0.3 + 1.2 * sqrt(sv_u01(&st)));
            inc = 0.01 * sv_u01(&st);
            break;
        case SV_CLOUD:
            a = r0 * (0.5 + sv_u01(&st));
            inc = acos(1.0 - 2.0 * sv_u01(&st));
            break;
        case SV_SHELL:
            a = r0;
            inc = acos(1.0 - 2.0 * sv_u01(&st));
            break;
        default:
            a = r0 * (0.9 + 0.2 * sv_u01(&st));
            inc = 0.035 * sv_u01(&st);
            break;
        }
        sv_kep(mu + (massive ? SV_GM_EARTH : 0.0), a, e, inc, node, argp, M,
               p, v);
        for (j = 0; j < 3; j++) {
            p[j] += (j == 0 ? s->x[c] : j == 1 ? s->y[c] : s->z[c]);
            v[j] += (j == 0 ? s->vx[c] : j == 1 ? s->vy[c] : s->vz[c]);
        }
        if (massive) {
            char nm[STC_SPACE_NAME];

            snprintf(nm, sizeof nm, "EARTH-%ld", k + 1);
            if (stc_space_add(s, nm, STC_SK_ROCK, c, SV_GM_EARTH, SV_R_EARTH,
                              p, v) < 0) {
                break;
            }
        } else if (stc_space_add_tracer(s, p, v) < 0) {
            break;
        }
        made++;
    }
    snprintf(msg, (size_t)msgn, "SPAWNED %ld %s ABOUT %s AT %.3g AU", made,
             massive ? "MASSIVE EARTHS" : "TEST EARTHS (MASSLESS)",
             s->name[c], r0 / SV_AU);
    return made == n ? 0 : -1;
}

static void sv_default_view(StcSVis *v)
{
    v->yaw = -1.9;
    v->pitch = 0.55;
    v->dist = 12.0 * SV_AU;
    v->focus = -1;
    v->select = -1;
}

static int sv_load(StcSVis *v, int k)
{
    StcSpace *s = &v->sim;
    int rc = 0;
    int sun;

    v->scenario = k;
    sv_default_view(v);
    v->warp = 10.0 * SV_DAY;
    v->msg[0] = '\0';
    switch (k) {
    case 0:
    case 1:
    case 2:
        rc = stc_ephem_load(s, k == 0 ? 0 : 1, STC_EPH_ALL, 5000, 20260927ULL);
        s->h = 1800.0;
        v->dist = 7.0 * SV_AU;
        if (k == 2 && rc > 0) {
            sun = sv_find(s, "Sun");
            (void)stc_space_make_blackhole(s, sun, 0.0, "BLACK HOLE");
            v->focus = -1;
            snprintf(v->msg, sizeof v->msg,
                     "SUN -> BLACK HOLE, SAME MASS: R_S = %.2f KM. ORBITS "
                     "UNCHANGED (BIRKHOFF).", 2.0 * s->gm[sun] / (SV_C * SV_C)
                     / 1000.0);
        }
        break;
    case 3:
        rc = stc_ephem_load(s, 1, STC_EPH_CORE | STC_EPH_DWARFS, 3000,
                            20260927ULL);
        s->h = 6.0 * 3600.0;
        if (rc > 0) {
            double gm = 10.0 * SV_GM_SUN;
            double p[3] = {-60.0 * SV_AU, 15.0 * SV_AU, 22.0 * SV_AU};
            double vv[3] = {45000.0, 0.0, -14000.0};

            if (stc_space_add(s, "ROGUE BH", STC_SK_BLACKHOLE, -1, gm,
                              2.0 * gm / (SV_C * SV_C), p, vv) < 0) {
                rc = -2;
            }
        }
        v->dist = 90.0 * SV_AU;
        v->warp = 0.5 * SV_YEAR;
        break;
    case 4:
        rc = stc_ephem_load(s, 1, STC_EPH_CORE, 0, 20260927ULL);
        s->h = 4.0 * 3600.0;
        if (rc > 0) {
            sun = sv_find(s, "Sun");
            rc = sv_spawn(s, sun, 2000, SV_RING, SV_AU, 1, 7ULL, v->msg,
                          (int)sizeof v->msg) == 0 ? rc : -2;
        }
        v->dist = 4.0 * SV_AU;
        v->warp = 5.0 * SV_DAY;
        break;
    default: {
        double p[3] = {0.0, 0.0, 0.0};
        double vv[3] = {0.0, 0.0, 0.0};
        double gm = 10.0 * SV_GM_SUN;
        int bh;

        stc_space_clear(s, 0.0);
        s->h = 3600.0;
        bh = stc_space_add(s, "BLACK HOLE", STC_SK_BLACKHOLE, -1, gm,
                           2.0 * gm / (SV_C * SV_C), p, vv);
        rc = bh < 0 ? -2 : 1;
        if (bh >= 0) {
            char dummy[160];

            if (sv_spawn(s, bh, 1000000, SV_DISK, 1.0 * SV_AU, 0, 11ULL,
                         v->msg, (int)sizeof v->msg) != 0) {
                rc = -2;
            }
            (void)sv_spawn(s, bh, 12, SV_RING, 0.6 * SV_AU, 1, 13ULL, dummy,
                           (int)sizeof dummy);
        }
        v->dist = 3.2 * SV_AU;
        v->warp = 2.0 * SV_DAY;
        break;
    }
    }
    s->e0 = stc_space_energy(s);
    if (rc <= 0) {
        snprintf(v->msg, sizeof v->msg, "SCENARIO %d FAILED TO LOAD (%d)", k,
                 rc);
        return -1;
    }
    return 0;
}

/* --- state ------------------------------------------------------------ */

int stc_svis_resize(StcSVis *v, int w, int h)
{
    if (v == NULL || w < 64 || h < 64) {
        return -1;
    }
    stc_render_rgb_free(&v->fb);
    free(v->acc);
    v->acc = NULL;
    v->w = w;
    v->h = h;
    if (stc_render_rgb_init(&v->fb, w, h) != 0) {
        return -1;
    }
    v->acc = (float *)calloc((size_t)w * (size_t)h * 3u, sizeof(float));
    if (v->acc == NULL) {
        stc_render_rgb_free(&v->fb);
        return -1;
    }
    return 0;
}

int stc_svis_init(StcSVis *v, int scenario, int w, int h)
{
    if (v == NULL) {
        return -1;
    }
    memset(v, 0, sizeof(*v));
    if (scenario < 0 || scenario >= STC_SVIS_NSCEN) {
        return -1;
    }
    if (stc_space_init(&v->sim, 256, 8192) != 0) {
        return -2;
    }
    v->show_orbits = 1;
    v->show_labels = 1;
    v->show_grid = 1;
    v->show_tracers = 1;
    v->show_hud = 1;
    if (stc_svis_resize(v, w, h) != 0) {
        stc_space_free(&v->sim);
        return -2;
    }
    (void)sv_load(v, scenario);
    return 0;
}

void stc_svis_free(StcSVis *v)
{
    if (v == NULL) {
        return;
    }
    stc_space_free(&v->sim);
    stc_render_rgb_free(&v->fb);
    free(v->acc);
    v->acc = NULL;
}

void stc_svis_step(StcSVis *v)
{
    long r;

    if (v == NULL) {
        return;
    }
    v->frames++;
    if (v->paused || v->console) {
        return;
    }
    r = stc_space_advance(&v->sim, v->warp / 60.0, sv_budget(&v->sim));
    if (r < 0) {
        v->limited++;
    } else {
        v->limited = 0;
    }
    if (v->focus >= 0 && !v->sim.alive[v->focus]) {
        v->focus = -1;
    }
    if (v->select >= 0 && !v->sim.alive[v->select]) {
        v->select = -1;
    }
}

/* --- camera ------------------------------------------------------------- */

typedef struct {
    double c[3];   /* camera target (world, m) */
    double e[3];   /* eye offset from target (m) */
    double f[3];   /* forward */
    double r[3];   /* right */
    double u[3];   /* up */
    double foc;    /* focal length, px */
    double cx;
    double cy;
    double ws;     /* ui scale */
} SvCam;

typedef struct {
    double x;
    double y;
    double z;      /* view depth, m */
    int ok;
} SvP;

static void sv_center(const StcSVis *v, double *c)
{
    const StcSpace *s = &v->sim;

    if (v->focus >= 0 && v->focus < s->nm && s->alive[v->focus]) {
        c[0] = s->x[v->focus];
        c[1] = s->y[v->focus];
        c[2] = s->z[v->focus];
    } else {
        c[0] = c[1] = c[2] = 0.0;
    }
}

static void sv_camera(const StcSVis *v, SvCam *C)
{
    double cp = cos(v->pitch);
    double l;

    sv_center(v, C->c);
    C->e[0] = v->dist * cp * cos(v->yaw);
    C->e[1] = v->dist * cp * sin(v->yaw);
    C->e[2] = v->dist * sin(v->pitch);
    C->f[0] = -C->e[0] / v->dist;
    C->f[1] = -C->e[1] / v->dist;
    C->f[2] = -C->e[2] / v->dist;
    /* right = f x z */
    C->r[0] = C->f[1];
    C->r[1] = -C->f[0];
    C->r[2] = 0.0;
    l = sqrt(C->r[0] * C->r[0] + C->r[1] * C->r[1]);
    C->r[0] /= l;
    C->r[1] /= l;
    /* up = r x f */
    C->u[0] = C->r[1] * C->f[2] - C->r[2] * C->f[1];
    C->u[1] = C->r[2] * C->f[0] - C->r[0] * C->f[2];
    C->u[2] = C->r[0] * C->f[1] - C->r[1] * C->f[0];
    C->foc = 0.5 * v->h / tan(0.5 * SV_FOV);
    C->cx = 0.5 * v->w;
    C->cy = 0.5 * v->h;
    C->ws = v->h / 800.0;
    if (C->ws < 1.0) {
        C->ws = 1.0;
    }
}

/* Project a world point given relative to the camera target (keeps
 * precision: callers subtract the target in double first). */
static SvP sv_proj_rel(const SvCam *C, double dx, double dy, double dz)
{
    SvP o;
    double px = dx - C->e[0];
    double py = dy - C->e[1];
    double pz = dz - C->e[2];

    o.z = px * C->f[0] + py * C->f[1] + pz * C->f[2];
    o.ok = o.z > 1e-6 * (fabs(C->e[0]) + fabs(C->e[1]) + fabs(C->e[2]));
    if (!o.ok) {
        o.x = o.y = 0.0;
        return o;
    }
    o.x = C->cx + C->foc * (px * C->r[0] + py * C->r[1] + pz * C->r[2]) / o.z;
    o.y = C->cy - C->foc * (px * C->u[0] + py * C->u[1] + pz * C->u[2]) / o.z;
    return o;
}

static SvP sv_proj(const SvCam *C, double x, double y, double z)
{
    return sv_proj_rel(C, x - C->c[0], y - C->c[1], z - C->c[2]);
}

/* --- linear accumulation raster -------------------------------------- */

static void sv_px(StcSVis *v, int x, int y, float r, float g, float b)
{
    float *p;

    if (x < 0 || y < 0 || x >= v->w || y >= v->h) {
        return;
    }
    p = v->acc + ((size_t)y * (size_t)v->w + (size_t)x) * 3u;
    p[0] += r;
    p[1] += g;
    p[2] += b;
}

/* Anti-aliased 1-px line (Wu), additive, colour in linear units. */
static void sv_line(StcSVis *v, double x0, double y0, double x1, double y1,
                    float r, float g, float b)
{
    int steep = fabs(y1 - y0) > fabs(x1 - x0);
    double t;
    double grad;
    int i0;
    int i1;
    int i;
    int lim;
    double k;

    if (!(isfinite(x0) && isfinite(y0) && isfinite(x1) && isfinite(y1))) {
        return;
    }
    if (steep) {
        t = x0; x0 = y0; y0 = t;
        t = x1; x1 = y1; y1 = t;
    }
    if (x0 > x1) {
        t = x0; x0 = x1; x1 = t;
        t = y0; y0 = y1; y1 = t;
    }
    lim = steep ? v->h : v->w;
    if (x1 < 0.0 || x0 > lim) {
        return;
    }
    if (x1 - x0 > 8.0 * (v->w + v->h)) {
        return;
    }
    grad = x1 - x0 > 1e-9 ? (y1 - y0) / (x1 - x0) : 0.0;
    k = sqrt(1.0 + grad * grad);
    i0 = (int)ceil(x0 - 0.5);
    i1 = (int)floor(x1 - 0.5);
    if (i0 < 0) {
        i0 = 0;
    }
    if (i1 > lim - 1) {
        i1 = lim - 1;
    }
    for (i = i0; i <= i1; i++) {
        double yc = y0 + grad * (i + 0.5 - x0) - 0.5;
        int yi = (int)floor(yc);
        float fr = (float)(yc - yi);
        float a0 = (float)((1.0 - fr) * k);
        float a1 = (float)(fr * k);

        if (steep) {
            sv_px(v, yi, i, r * a0, g * a0, b * a0);
            sv_px(v, yi + 1, i, r * a1, g * a1, b * a1);
        } else {
            sv_px(v, i, yi, r * a0, g * a0, b * a0);
            sv_px(v, i, yi + 1, r * a1, g * a1, b * a1);
        }
    }
}

static void sv_dot(StcSVis *v, double x, double y, float r, float g, float b)
{
    int xi = (int)floor(x - 0.5);
    int yi = (int)floor(y - 0.5);
    float fx = (float)(x - 0.5 - xi);
    float fy = (float)(y - 0.5 - yi);

    sv_px(v, xi, yi, r * (1 - fx) * (1 - fy), g * (1 - fx) * (1 - fy),
          b * (1 - fx) * (1 - fy));
    sv_px(v, xi + 1, yi, r * fx * (1 - fy), g * fx * (1 - fy),
          b * fx * (1 - fy));
    sv_px(v, xi, yi + 1, r * (1 - fx) * fy, g * (1 - fx) * fy,
          b * (1 - fx) * fy);
    sv_px(v, xi + 1, yi + 1, r * fx * fy, g * fx * fy, b * fx * fy);
}

static void sv_circle(StcSVis *v, double cx, double cy, double rad, float r,
                      float g, float b)
{
    int n = (int)(rad * 0.8) + 16;
    int k;
    double px = cx + rad;
    double py = cy;

    if (n > 720) {
        n = 720;
    }
    for (k = 1; k <= n; k++) {
        double a = 2.0 * SV_PI * k / n;
        double qx = cx + rad * cos(a);
        double qy = cy + rad * sin(a);

        sv_line(v, px, py, qx, qy, r, g, b);
        px = qx;
        py = qy;
    }
}

/* --- colours -------------------------------------------------------------- */

static void sv_body_col(const StcSpace *s, int i, float *c)
{
    const char *n = s->name[i];

    c[0] = c[1] = c[2] = 0.78f;
    switch (s->kind[i]) {
    case STC_SK_STAR:
        c[0] = 1.00f; c[1] = 0.96f; c[2] = 0.86f;
        return;
    case STC_SK_BLACKHOLE:
        c[0] = c[1] = c[2] = 0.0f;
        return;
    case STC_SK_ROCK:
        c[0] = 0.42f; c[1] = 0.60f; c[2] = 0.82f;
        return;
    default:
        break;
    }
    if (strcmp(n, "Earth") == 0) {
        c[0] = 0.36f; c[1] = 0.55f; c[2] = 0.86f;
    } else if (strcmp(n, "Mars") == 0) {
        c[0] = 0.80f; c[1] = 0.46f; c[2] = 0.30f;
    } else if (strcmp(n, "Jupiter") == 0) {
        c[0] = 0.82f; c[1] = 0.70f; c[2] = 0.56f;
    } else if (strcmp(n, "Saturn") == 0) {
        c[0] = 0.86f; c[1] = 0.78f; c[2] = 0.58f;
    } else if (strcmp(n, "Uranus") == 0) {
        c[0] = 0.62f; c[1] = 0.82f; c[2] = 0.86f;
    } else if (strcmp(n, "Neptune") == 0) {
        c[0] = 0.40f; c[1] = 0.55f; c[2] = 0.90f;
    } else if (strcmp(n, "Venus") == 0) {
        c[0] = 0.88f; c[1] = 0.82f; c[2] = 0.66f;
    } else if (strcmp(n, "Io") == 0) {
        c[0] = 0.88f; c[1] = 0.80f; c[2] = 0.45f;
    }
}

/* --- scene painters ----------------------------------------------------- */

typedef struct {
    StcSVis *v;
    SvCam C;
} SvCtx;

/* Nice ring spacing for a view scale, in metres. */
static double sv_nice(double x)
{
    double unit = x >= 0.05 * SV_AU ? SV_AU : 1000.0;
    double d;

    x /= unit;
    d = pow(10.0, floor(log10(x))) * unit;
    x *= unit;
    double m = x / d;

    if (m < 2.0) {
        return d;
    }
    if (m < 5.0) {
        return 2.0 * d;
    }
    return 5.0 * d;
}

static void sv_fmt_len(double m, char *buf, int n)
{
    double a = fabs(m);

    if (a >= 0.05 * SV_AU) {
        snprintf(buf, (size_t)n, "%.6f AU", m / SV_AU);
    } else if (a >= 1000.0) {
        snprintf(buf, (size_t)n, "%.1f KM", m / 1000.0);
    } else {
        snprintf(buf, (size_t)n, "%.1f M", m);
    }
}

/* Reference grid: rings in the ecliptic plane about the view centre with
 * tick marks and range labels, 12 longitude spokes. */
static void sv_paint_grid(SvCtx *X, double *ring_step, int *ring_n)
{
    StcSVis *v = X->v;
    const SvCam *C = &X->C;
    double step = sv_nice(v->dist * 0.22);
    int nr = 6;
    int k;
    float g = 0.022f;

    *ring_step = step;
    *ring_n = nr;
    for (k = 1; k <= nr; k++) {
        double r = step * k;
        int n = 240;
        int j;
        SvP prev;
        int have = 0;

        for (j = 0; j <= n; j++) {
            double a = 2.0 * SV_PI * j / n;
            SvP p = sv_proj_rel(C, r * cos(a), r * sin(a), 0.0);

            if (p.ok && have && (j % 2 == 0)) {
                sv_line(v, prev.x, prev.y, p.x, p.y, g, g, g * 1.15f);
            }
            prev = p;
            have = p.ok;
        }
    }
    for (k = 0; k < 12; k++) {
        double a = 2.0 * SV_PI * k / 12.0;
        SvP p0 = sv_proj_rel(C, step * 0.25 * cos(a), step * 0.25 * sin(a),
                             0.0);
        SvP p1 = sv_proj_rel(C, step * (nr + 0.4) * cos(a),
                             step * (nr + 0.4) * sin(a), 0.0);

        if (p0.ok && p1.ok) {
            sv_line(v, p0.x, p0.y, p1.x, p1.y, g * 0.7f, g * 0.7f,
                    g * 0.8f);
        }
        /* ticks every half ring on the spoke */
        {
            int t;

            for (t = 1; t <= 2 * nr; t++) {
                double rr = step * 0.5 * t;
                double tl = step * (t % 2 == 0 ? 0.06 : 0.03);
                double ca = cos(a);
                double sa = sin(a);
                SvP q0 = sv_proj_rel(C, rr * ca - tl * sa, rr * sa + tl * ca,
                                     0.0);
                SvP q1 = sv_proj_rel(C, rr * ca + tl * sa, rr * sa - tl * ca,
                                     0.0);

                if (q0.ok && q1.ok) {
                    sv_line(v, q0.x, q0.y, q1.x, q1.y, 0.05f, 0.05f, 0.055f);
                }
            }
        }
    }
}

/* Osculating orbit trace of body i about its primary. */
static void sv_paint_orbit(SvCtx *X, int i, float br)
{
    StcSVis *v = X->v;
    const StcSpace *s = &v->sim;
    const SvCam *C = &X->C;
    StcOrbit o;
    int j;
    int n = 256;
    double cO;
    double sO;
    double ci;
    double si;
    double cw;
    double sw;
    double pc[3];
    SvP prev;
    int have = 0;
    double numax;
    int prim;
    float col[3];

    if (stc_space_orbit(s, i, -1, &o) != 0 || !(o.e >= 0.0)) {
        return;
    }
    prim = s->parent[i] >= 0 && s->alive[s->parent[i]] ? s->parent[i]
                                                        : sv_heaviest(s);
    if (prim < 0 || prim == i) {
        return;
    }
    pc[0] = s->x[prim] - C->c[0];
    pc[1] = s->y[prim] - C->c[1];
    pc[2] = s->z[prim] - C->c[2];
    /* skip orbits smaller than ~2 px on screen */
    {
        SvP a = sv_proj_rel(C, pc[0], pc[1], pc[2]);
        double rr = o.e < 1.0 ? o.a * (1.0 + o.e) : o.r;

        if (!a.ok || C->foc * rr / a.z < 2.5) {
            return;
        }
    }
    cO = cos(o.node);
    sO = sin(o.node);
    ci = cos(o.inc);
    si = sin(o.inc);
    cw = cos(o.argp);
    sw = sin(o.argp);
    numax = o.e < 1.0 ? SV_PI : acos(-1.0 / o.e) * 0.97;
    sv_body_col(s, i, col);
    col[0] = br * (0.55f + 0.45f * col[0]);
    col[1] = br * (0.55f + 0.45f * col[1]);
    col[2] = br * (0.55f + 0.45f * col[2]);
    for (j = 0; j <= n; j++) {
        double nu = -numax + 2.0 * numax * j / n;
        double p = o.a * (1.0 - o.e * o.e);
        double r;
        double xo;
        double yo;
        double w[3];
        SvP q;

        if (o.e >= 1.0) {
            p = fabs(o.a) * (o.e * o.e - 1.0);
        }
        r = p / (1.0 + o.e * cos(nu));
        if (!(r > 0.0) || r > 1e4 * SV_AU) {
            have = 0;
            continue;
        }
        xo = r * cos(nu);
        yo = r * sin(nu);
        w[0] = (cO * cw - sO * sw * ci) * xo + (-cO * sw - sO * cw * ci) * yo;
        w[1] = (sO * cw + cO * sw * ci) * xo + (-sO * sw + cO * cw * ci) * yo;
        w[2] = (sw * si) * xo + (cw * si) * yo;
        q = sv_proj_rel(C, pc[0] + w[0], pc[1] + w[1], pc[2] + w[2]);
        if (q.ok && have) {
            sv_line(v, prev.x, prev.y, q.x, q.y, col[0], col[1], col[2]);
        }
        prev = q;
        have = q.ok;
    }
}

/* Tracer projection runs on the pool (pure per-index writes into a
 * screen-coordinate scratch array); the additive splat stays serial so
 * the image is identical for any thread count. */
typedef struct {
    const StcSpace *s;
    const SvCam *C;
    float *xy;
} SvTr;

static float *g_sv_trxy = NULL;
static int g_sv_trcap = 0;

static void sv_tr_rows(void *ctx, int i0, int i1)
{
    const SvTr *T = (const SvTr *)ctx;
    int i;

    for (i = i0; i < i1; i++) {
        float *o = T->xy + 2 * (size_t)i;

        o[0] = -1.0f;
        o[1] = -1.0f;
        if (T->s->talive[i]) {
            SvP p = sv_proj(T->C, T->s->tx[i], T->s->ty[i], T->s->tz[i]);

            if (p.ok && p.x > -2.0 && p.y > -2.0 && p.x < 1e6 &&
                p.y < 1e6) {
                o[0] = (float)p.x;
                o[1] = (float)p.y;
            }
        }
    }
}

static void sv_paint_tracers(SvCtx *X)
{
    StcSVis *v = X->v;
    const StcSpace *s = &v->sim;
    SvTr T;
    int i;
    int big = s->nt > 200000;
    /* per-particle brightness falls with population so dense swarms
     * read as structure, not a saturated blob */
    float k = (float)(0.22 * (s->nt > 20000 ? pow(20000.0 / s->nt, 0.6)
                                            : 1.0));

    if (s->nt <= 0) {
        return;
    }
    if (s->nt > g_sv_trcap) {
        float *n = (float *)realloc(g_sv_trxy, 2u * (size_t)s->nt *
                                                   sizeof(float));

        if (n == NULL) {
            return;
        }
        g_sv_trxy = n;
        g_sv_trcap = s->nt;
    }
    T.s = s;
    T.C = &X->C;
    T.xy = g_sv_trxy;
    stc_par_for(sv_tr_rows, &T, s->nt, 0);
    for (i = 0; i < s->nt; i++) {
        float x = g_sv_trxy[2 * (size_t)i];
        float y = g_sv_trxy[2 * (size_t)i + 1];

        if (x < 0.0f) {
            continue;
        }
        if (big) {
            sv_px(v, (int)x, (int)y, k * 0.85f, k * 0.88f, k);
        } else {
            sv_dot(v, x, y, k * 0.85f, k * 0.88f, k);
        }
    }
}

/* Filled, Sun-lit sphere (Lambert, sharp terminator) with AA edge. */
static void sv_sphere(StcSVis *v, double cx, double cy, double rad,
                      const float *col, const double *ldir, int emissive)
{
    int x0 = (int)floor(cx - rad - 1.0);
    int x1 = (int)ceil(cx + rad + 1.0);
    int y0 = (int)floor(cy - rad - 1.0);
    int y1 = (int)ceil(cy + rad + 1.0);
    int x;
    int y;

    if (x1 < 0 || y1 < 0 || x0 >= v->w || y0 >= v->h) {
        return;
    }
    if (x0 < 0) {
        x0 = 0;
    }
    if (y0 < 0) {
        y0 = 0;
    }
    if (x1 >= v->w) {
        x1 = v->w - 1;
    }
    if (y1 >= v->h) {
        y1 = v->h - 1;
    }
    for (y = y0; y <= y1; y++) {
        for (x = x0; x <= x1; x++) {
            double dx = (x + 0.5 - cx) / rad;
            double dy = (y + 0.5 - cy) / rad;
            double d2 = dx * dx + dy * dy;
            double cov = (1.0 - sqrt(d2)) * rad + 0.5;
            double sh;
            float *p;

            if (cov <= 0.0) {
                continue;
            }
            if (cov > 1.0) {
                cov = 1.0;
            }
            if (emissive) {
                sh = 1.0 - 0.25 * d2;
            } else {
                double nz = sqrt(d2 < 1.0 ? 1.0 - d2 : 0.0);
                double l = dx * ldir[0] - dy * ldir[1] + nz * ldir[2];

                sh = 0.03 + 0.97 * (l > 0.0 ? l : 0.0);
            }
            p = v->acc + ((size_t)y * (size_t)v->w + (size_t)x) * 3u;
            p[0] = (float)(p[0] * (1.0 - cov) + col[0] * sh * cov);
            p[1] = (float)(p[1] * (1.0 - cov) + col[1] * sh * cov);
            p[2] = (float)(p[2] * (1.0 - cov) + col[2] * sh * cov);
        }
    }
}

/* Screen record of a painted body, for labels and picking. */
typedef struct {
    int i;
    double x;
    double y;
    double rad;
    double z;
} SvMark;

static int sv_mark_cmp(const void *a, const void *b)
{
    const SvMark *p = (const SvMark *)a;
    const SvMark *q = (const SvMark *)b;

    if (p->z > q->z) {
        return -1;
    }
    if (p->z < q->z) {
        return 1;
    }
    return p->i - q->i;
}

static int sv_collect(const StcSVis *v, const SvCam *C, SvMark *m, int cap)
{
    const StcSpace *s = &v->sim;
    int i;
    int n = 0;

    for (i = 0; i < s->nm && n < cap; i++) {
        SvP p;

        if (!s->alive[i]) {
            continue;
        }
        p = sv_proj(C, s->x[i], s->y[i], s->z[i]);
        if (!p.ok) {
            continue;
        }
        m[n].i = i;
        m[n].x = p.x;
        m[n].y = p.y;
        m[n].z = p.z;
        m[n].rad = C->foc * s->radius[i] / p.z;
        n++;
    }
    return n;
}

static void sv_paint_bodies(SvCtx *X, SvMark *m, int n)
{
    StcSVis *v = X->v;
    const StcSpace *s = &v->sim;
    int star = -1;
    int k;

    for (k = 0; k < s->nm; k++) {
        if (s->alive[k] && s->kind[k] == STC_SK_STAR &&
            (star < 0 || s->gm[k] > s->gm[star])) {
            star = k;
        }
    }
    qsort(m, (size_t)n, sizeof *m, sv_mark_cmp);
    for (k = 0; k < n; k++) {
        int i = m[k].i;
        float col[3];
        double l[3] = {0.0, 0.0, 1.0};

        sv_body_col(s, i, col);
        if (star >= 0 && star != i) {
            /* light direction in screen basis (x right, y up, z to viewer) */
            double d[3];
            double len;
            const SvCam *C = &X->C;

            d[0] = s->x[star] - s->x[i];
            d[1] = s->y[star] - s->y[i];
            d[2] = s->z[star] - s->z[i];
            len = sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
            if (len > 0.0) {
                l[0] = (d[0] * C->r[0] + d[1] * C->r[1] + d[2] * C->r[2]) / len;
                l[1] = (d[0] * C->u[0] + d[1] * C->u[1] + d[2] * C->u[2]) / len;
                l[2] = -(d[0] * C->f[0] + d[1] * C->f[1] + d[2] * C->f[2]) /
                       len;
            }
        }
        if (s->kind[i] == STC_SK_BLACKHOLE) {
            double rsh = 2.598 * m[k].rad; /* shadow: sqrt(27)/2 r_s */
            float z0[3] = {0.0f, 0.0f, 0.0f};

            if (rsh >= 1.5) {
                sv_sphere(v, m[k].x, m[k].y, rsh, z0, l, 1);
                sv_circle(v, m[k].x, m[k].y, rsh + 0.6, 0.55f, 0.55f, 0.58f);
            } else {
                sv_circle(v, m[k].x, m[k].y, 4.0 * X->C.ws, 0.70f, 0.70f,
                          0.72f);
                sv_line(v, m[k].x - 7 * X->C.ws, m[k].y, m[k].x - 4.5 *
                        X->C.ws, m[k].y, 0.7f, 0.7f, 0.72f);
                sv_line(v, m[k].x + 4.5 * X->C.ws, m[k].y, m[k].x + 7 *
                        X->C.ws, m[k].y, 0.7f, 0.7f, 0.72f);
            }
            continue;
        }
        if (m[k].rad >= 1.2) {
            sv_sphere(v, m[k].x, m[k].y, m[k].rad, col, l,
                      s->kind[i] == STC_SK_STAR);
        } else if (s->kind[i] == STC_SK_ROCK) {
            sv_dot(v, m[k].x, m[k].y, col[0] * 0.9f, col[1] * 0.9f,
                   col[2] * 0.9f);
        } else {
            double r = (s->kind[i] == STC_SK_STAR ? 3.5 :
                        s->kind[i] == STC_SK_PLANET ? 2.6 : 1.8) * X->C.ws;

            sv_circle(v, m[k].x, m[k].y, r, col[0] * 0.9f, col[1] * 0.9f,
                      col[2] * 0.9f);
            sv_dot(v, m[k].x, m[k].y, col[0], col[1], col[2]);
        }
    }
}

/* --- final conversion ------------------------------------------------- */

static uint8_t g_sv_lut[4097];
static int g_sv_lut_ok = 0;

static void sv_lut(void)
{
    int i;

    if (g_sv_lut_ok) {
        return;
    }
    for (i = 0; i <= 4096; i++) {
        double a = i / 4096.0;
        double c = a <= 0.0031308 ? 12.92 * a
                                  : 1.055 * pow(a, 1.0 / 2.4) - 0.055;

        g_sv_lut[i] = (uint8_t)(c * 255.0 + 0.5);
    }
    g_sv_lut_ok = 1;
}

static void sv_final_rows(void *ctx, int y0, int y1)
{
    StcSVis *v = (StcSVis *)ctx;
    int y;

    for (y = y0; y < y1; y++) {
        int x;

        for (x = 0; x < v->w; x++) {
            size_t o = (size_t)y * (size_t)v->w + (size_t)x;
            const float *p = v->acc + o * 3u;
            int c[3];
            int k;

            for (k = 0; k < 3; k++) {
                float a = p[k];

                if (!(a > 0.0f)) {
                    a = 0.0f;
                }
                if (a > 1.0f) {
                    a = 1.0f;
                }
                c[k] = g_sv_lut[(int)(a * 4096.0f + 0.5f)];
            }
            v->fb.px[o].r = (uint8_t)c[0];
            v->fb.px[o].g = (uint8_t)c[1];
            v->fb.px[o].b = (uint8_t)c[2];
        }
    }
}

/* --- HUD (8-bit, after conversion) ---------------------------------- */

static void hud_blend(StcRgbBuf *fb, int x, int y, StcRgb c, double a)
{
    StcRgb *p;

    if (x < 0 || y < 0 || x >= fb->w || y >= fb->h || a <= 0.0) {
        return;
    }
    if (a > 1.0) {
        a = 1.0;
    }
    p = &fb->px[(size_t)y * (size_t)fb->w + (size_t)x];
    p->r = (uint8_t)(p->r + (c.r - p->r) * a + 0.5);
    p->g = (uint8_t)(p->g + (c.g - p->g) * a + 0.5);
    p->b = (uint8_t)(p->b + (c.b - p->b) * a + 0.5);
}

/* Glyphs at fractional scale with box-filter coverage (crisp mono text
 * at any UI scale). Returns the pen x after the string. */
static double hud_text(StcRgbBuf *fb, double x, double y, const char *s,
                       StcRgb c, double sc)
{
    double adv = 6.0 * sc;

    for (; *s != '\0'; s++) {
        const unsigned char *g = stc_vis_glyph(*s);
        int x0 = (int)floor(x);
        int x1 = (int)ceil(x + 5.0 * sc);
        int y0 = (int)floor(y);
        int y1 = (int)ceil(y + 7.0 * sc);
        int px;
        int py;

        for (py = y0; py < y1; py++) {
            for (px = x0; px < x1; px++) {
                double cov = 0.0;
                int r;

                for (r = 0; r < 7; r++) {
                    double gy0 = y + r * sc;
                    double oy = (py + 1 < gy0 + sc ? py + 1 : gy0 + sc) -
                                (py > gy0 ? py : gy0);
                    int q;

                    if (oy <= 0.0) {
                        continue;
                    }
                    for (q = 0; q < 5; q++) {
                        double gx0;
                        double ox;

                        if (!((g[r] >> (4 - q)) & 1u)) {
                            continue;
                        }
                        gx0 = x + q * sc;
                        ox = (px + 1 < gx0 + sc ? px + 1 : gx0 + sc) -
                             (px > gx0 ? px : gx0);
                        if (ox > 0.0) {
                            cov += ox * oy;
                        }
                    }
                }
                if (cov > 0.0) {
                    hud_blend(fb, px, py, c, cov);
                }
            }
        }
        x += adv;
    }
    return x;
}

static double hud_tw(const char *s, double sc)
{
    return (double)strlen(s) * 6.0 * sc;
}

static void hud_hline(StcRgbBuf *fb, int x0, int x1, int y, StcRgb c,
                      double a)
{
    int x;

    for (x = x0; x <= x1; x++) {
        hud_blend(fb, x, y, c, a);
    }
}

static void hud_vline(StcRgbBuf *fb, int x, int y0, int y1, StcRgb c,
                      double a)
{
    int y;

    for (y = y0; y <= y1; y++) {
        hud_blend(fb, x, y, c, a);
    }
}

/* Framed panel: darkened interior, 1 px rule, corner ticks. */
static void hud_panel(StcRgbBuf *fb, int x, int y, int w, int h)
{
    StcRgb blk = stc_render_rgb(0, 0, 0);
    StcRgb rule = stc_render_rgb(120, 124, 132);
    int i;
    int j;

    for (j = y; j < y + h; j++) {
        for (i = x; i < x + w; i++) {
            hud_blend(fb, i, j, blk, 0.72);
        }
    }
    hud_hline(fb, x, x + w - 1, y, rule, 0.55);
    hud_hline(fb, x, x + w - 1, y + h - 1, rule, 0.55);
    hud_vline(fb, x, y, y + h - 1, rule, 0.55);
    hud_vline(fb, x + w - 1, y, y + h - 1, rule, 0.55);
    for (i = 0; i < 6; i++) {
        hud_blend(fb, x + i, y + 1, rule, 1.0);
        hud_blend(fb, x + 1, y + i, rule, 1.0);
        hud_blend(fb, x + w - 1 - i, y + h - 2, rule, 1.0);
        hud_blend(fb, x + w - 2, y + h - 1 - i, rule, 1.0);
    }
}

static const char *sv_kind_name(int k)
{
    switch (k) {
    case STC_SK_STAR: return "STAR";
    case STC_SK_PLANET: return "PLANET";
    case STC_SK_DWARF: return "DWARF PLANET";
    case STC_SK_MOON: return "MOON";
    case STC_SK_ASTEROID: return "ASTEROID";
    case STC_SK_BLACKHOLE: return "BLACK HOLE";
    default: return "BODY";
    }
}

static void sv_warp_str(double w, char *buf, int n)
{
    double a = fabs(w);
    const char *sg = w < 0.0 ? "-" : "+";

    if (a >= SV_YEAR) {
        snprintf(buf, (size_t)n, "%s%.3g YR/S", sg, a / SV_YEAR);
    } else if (a >= SV_DAY) {
        snprintf(buf, (size_t)n, "%s%.3g D/S", sg, a / SV_DAY);
    } else if (a >= 3600.0) {
        snprintf(buf, (size_t)n, "%s%.3g H/S", sg, a / 3600.0);
    } else if (a >= 60.0) {
        snprintf(buf, (size_t)n, "%s%.3g MIN/S", sg, a / 60.0);
    } else {
        snprintf(buf, (size_t)n, "%s%.3g S/S", sg, a);
    }
}

static void sv_upper(char *s)
{
    for (; *s != '\0'; s++) {
        *s = (char)toupper((unsigned char)*s);
    }
}

static void sv_paint_labels(SvCtx *X, const SvMark *m, int n)
{
    StcSVis *v = X->v;
    const StcSpace *s = &v->sim;
    double sc = 1.25 * X->C.ws;
    StcRgb lab = stc_render_rgb(196, 200, 208);
    StcRgb dim = stc_render_rgb(120, 124, 132);
    StcRgb lead = stc_render_rgb(90, 94, 100);
    int k;

    for (k = 0; k < n; k++) {
        int i = m[k].i;
        char nm[STC_SPACE_NAME + 8];
        int major;
        double off;
        int t;

        major = s->kind[i] == STC_SK_STAR || s->kind[i] == STC_SK_PLANET ||
                s->kind[i] == STC_SK_BLACKHOLE;
        if (s->kind[i] == STC_SK_ROCK && i != v->select) {
            continue;
        }
        if (!major && i != v->select && i != v->focus) {
            /* minor bodies: label only when clear of their primary */
            int p = s->parent[i];

            if (p >= 0 && s->alive[p]) {
                SvP q = sv_proj(&X->C, s->x[p], s->y[p], s->z[p]);

                if (!q.ok || hypot(q.x - m[k].x, q.y - m[k].y) < 28.0 *
                    X->C.ws) {
                    continue;
                }
            } else if (s->kind[i] == STC_SK_ASTEROID) {
                continue;
            }
        }
        if (m[k].x < -50 || m[k].y < -50 || m[k].x > v->w + 50 ||
            m[k].y > v->h + 50) {
            continue;
        }
        snprintf(nm, sizeof nm, "%s", s->name[i]);
        sv_upper(nm);
        off = (m[k].rad > 3.0 ? m[k].rad : 3.0) + 5.0 * X->C.ws;
        for (t = 0; t < (int)(6 * X->C.ws); t++) {
            hud_blend(&v->fb, (int)(m[k].x + off * 0.7 + t),
                      (int)(m[k].y - off * 0.7 - t), lead, 0.8);
        }
        hud_text(&v->fb, m[k].x + off * 0.7 + 7.0 * X->C.ws,
                 m[k].y - off * 0.7 - 6.0 * X->C.ws - 7.0 * sc * 0.5, nm,
                 major ? lab : dim, sc);
    }
    /* selection reticle */
    if (v->select >= 0 && s->alive[v->select]) {
        SvP p = sv_proj(&X->C, s->x[v->select], s->y[v->select],
                        s->z[v->select]);

        if (p.ok) {
            double r = fmax(X->C.foc * s->radius[v->select] / p.z, 4.0) +
                       7.0 * X->C.ws;
            StcRgb ret = stc_render_rgb(230, 232, 236);
            int t;
            int L = (int)(6 * X->C.ws);

            for (t = 0; t < L; t++) {
                hud_blend(&v->fb, (int)(p.x - r + t), (int)(p.y - r), ret, 1);
                hud_blend(&v->fb, (int)(p.x - r), (int)(p.y - r + t), ret, 1);
                hud_blend(&v->fb, (int)(p.x + r - t), (int)(p.y - r), ret, 1);
                hud_blend(&v->fb, (int)(p.x + r), (int)(p.y - r + t), ret, 1);
                hud_blend(&v->fb, (int)(p.x - r + t), (int)(p.y + r), ret, 1);
                hud_blend(&v->fb, (int)(p.x - r), (int)(p.y + r - t), ret, 1);
                hud_blend(&v->fb, (int)(p.x + r - t), (int)(p.y + r), ret, 1);
                hud_blend(&v->fb, (int)(p.x + r), (int)(p.y + r - t), ret, 1);
            }
        }
    }
}

static void sv_paint_hud(SvCtx *X, double ring_step, int ring_n)
{
    StcSVis *v = X->v;
    const StcSpace *s = &v->sim;
    StcRgbBuf *fb = &v->fb;
    double ws = X->C.ws;
    double sc = 1.35 * ws;
    double sm = 1.15 * ws;
    double lh = 10.5 * sc;
    StcRgb hi = stc_render_rgb(226, 228, 232);
    StcRgb md = stc_render_rgb(168, 172, 180);
    StcRgb lo = stc_render_rgb(112, 116, 124);
    StcRgb amber = stc_render_rgb(236, 176, 72);
    char buf[200];
    char tb[64];
    char wb[32];
    int pad = (int)(10 * ws);
    double y;
    double jd = stc_ephem_jd(s->t);

    if (!v->show_hud) {
        return;
    }
    /* top bar */
    hud_panel(fb, pad, pad, v->w - 2 * pad, (int)(2 * lh + 10 * ws));
    hud_text(fb, pad + 8 * ws, pad + 6 * ws,
             "SPACETIME LAB  //  N-BODY DYNAMICS", hi, sc);
    stc_ephem_format(s->t, tb, (int)sizeof tb);
    snprintf(buf, sizeof buf, "%s   JD %.5f", tb, jd);
    hud_text(fb, v->w - pad - 8 * ws - hud_tw(buf, sc), pad + 6 * ws, buf,
             hi, sc);
    snprintf(buf, sizeof buf, "%s", kScen[v->scenario]);
    hud_text(fb, pad + 8 * ws, pad + 6 * ws + lh, buf, md, sm);
    sv_warp_str(v->warp, wb, (int)sizeof wb);
    {
        double e = stc_space_energy(s);
        double de = s->e0 != 0.0 ? (e - s->e0) / fabs(s->e0) : 0.0;

        snprintf(buf, sizeof buf,
                 "WARP %s%s  STEP %.0f S  N %d + %d TR  1PN %s  DE/E %+.1E",
                 wb, v->paused ? " PAUSED" : v->limited > 2 ? " LIMITED" : "",
                 s->h, stc_space_live_massive(s), stc_space_live_tracers(s),
                 s->pn ? "ON" : "OFF", de);
    }
    hud_text(fb, v->w - pad - 8 * ws - hud_tw(buf, sm), pad + 6 * ws + lh,
             buf, v->limited > 2 ? amber : md, sm);

    /* selected body panel */
    if (v->select >= 0 && s->alive[v->select]) {
        int i = v->select;
        StcOrbit o;
        int ok = stc_space_orbit(s, i, -1, &o) == 0;
        int pw = (int)(300 * ws);
        int ph = (int)(11.5 * lh);
        int px0 = v->w - pad - pw;
        int py0 = (int)(pad + 2 * lh + 18 * ws);
        char nm[STC_SPACE_NAME + 8];
        int prim = s->parent[i] >= 0 && s->alive[s->parent[i]]
                       ? s->parent[i] : sv_heaviest(s);

        hud_panel(fb, px0, py0, pw, ph);
        snprintf(nm, sizeof nm, "%s", s->name[i]);
        sv_upper(nm);
        y = py0 + 7 * ws;
        hud_text(fb, px0 + 8 * ws, y, nm, hi, sc);
        hud_text(fb, px0 + pw - 8 * ws - hud_tw(sv_kind_name(s->kind[i]), sm),
                 y, sv_kind_name(s->kind[i]), lo, sm);
        y += lh * 1.3;
        snprintf(buf, sizeof buf, "GM      %.6E M3/S2", s->gm[i]);
        hud_text(fb, px0 + 8 * ws, y, buf, md, sm);
        y += lh;
        snprintf(buf, sizeof buf, "MASS    %.4E KG", s->gm[i] / SV_G);
        hud_text(fb, px0 + 8 * ws, y, buf, md, sm);
        y += lh;
        snprintf(buf, sizeof buf, "RADIUS  %.1f KM%s", s->radius[i] / 1000.0,
                 s->kind[i] == STC_SK_BLACKHOLE ? " (R_S)" : "");
        hud_text(fb, px0 + 8 * ws, y, buf, md, sm);
        y += lh;
        if (ok && prim >= 0 && prim != i) {
            char pn[STC_SPACE_NAME + 8];
            char lb[32];

            snprintf(pn, sizeof pn, "%s", s->name[prim]);
            sv_upper(pn);
            snprintf(buf, sizeof buf, "ORBITS  %s", pn);
            hud_text(fb, px0 + 8 * ws, y, buf, md, sm);
            y += lh;
            sv_fmt_len(o.r, lb, (int)sizeof lb);
            snprintf(buf, sizeof buf, "R       %s", lb);
            hud_text(fb, px0 + 8 * ws, y, buf, hi, sm);
            y += lh;
            snprintf(buf, sizeof buf, "V       %.4f KM/S", o.v / 1000.0);
            hud_text(fb, px0 + 8 * ws, y, buf, hi, sm);
            y += lh;
            sv_fmt_len(o.a, lb, (int)sizeof lb);
            snprintf(buf, sizeof buf, "A %s  E %.6f", lb, o.e);
            hud_text(fb, px0 + 8 * ws, y, buf, md, sm);
            y += lh;
            snprintf(buf, sizeof buf, "I %.4f DEG  NODE %.2f DEG",
                     o.inc * 180.0 / SV_PI, o.node * 180.0 / SV_PI);
            hud_text(fb, px0 + 8 * ws, y, buf, md, sm);
            y += lh;
            if (o.period > 0.0) {
                if (o.period > 2.0 * SV_YEAR) {
                    snprintf(buf, sizeof buf, "PERIOD  %.4f YR",
                             o.period / SV_YEAR);
                } else {
                    snprintf(buf, sizeof buf, "PERIOD  %.4f D",
                             o.period / SV_DAY);
                }
            } else {
                snprintf(buf, sizeof buf, "PERIOD  UNBOUND");
            }
            hud_text(fb, px0 + 8 * ws, y, buf, md, sm);
            y += lh;
            snprintf(buf, sizeof buf, "LIGHT TIME TO PRIMARY %.1f S",
                     o.r / SV_C);
            hud_text(fb, px0 + 8 * ws, y, buf, lo, sm);
        }
    }

    /* event log */
    {
        int shown = s->nev < 4 ? s->nev : 4;
        int k;
        double yy = v->h - pad - (shown + 1) * lh - 4 * lh;

        for (k = 0; k < shown; k++) {
            const StcSpaceEvent *e = &s->ev[(s->nev - 1 - k) % STC_SPACE_EVLOG];
            const char *what = e->kind == STC_SEV_MERGE ? "MERGE" :
                               e->kind == STC_SEV_ABSORB ? "IMPACT" :
                               e->kind == STC_SEV_CAPTURE ? "CAPTURE" :
                               "TIDAL DISRUPTION";
            char an[STC_SPACE_NAME + 8];
            char bn[STC_SPACE_NAME + 8];

            stc_ephem_format(e->t, tb, (int)sizeof tb);
            tb[10] = '\0';
            snprintf(an, sizeof an, "%s", e->a >= 0 && e->a < s->nm
                     ? s->name[e->a] : "?");
            if (e->b >= 0 && e->b < s->nm && e->kind != STC_SEV_ABSORB) {
                snprintf(bn, sizeof bn, "%s", s->name[e->b]);
            } else {
                snprintf(bn, sizeof bn, "TEST PARTICLE");
            }
            sv_upper(an);
            sv_upper(bn);
            snprintf(buf, sizeof buf, "%s  %-16s %s <- %s", tb, what, an, bn);
            hud_text(fb, pad + 8 * ws, yy, buf, k == 0 ? amber : lo, sm);
            yy += lh;
        }
    }

    /* scale bar (grid ring spacing) */
    {
        SvP a = sv_proj_rel(&X->C, 0.0, 0.0, 0.0);
        char lb[32];

        if (a.ok) {
            double px = X->C.foc * ring_step / a.z;
            int x1 = (int)(v->w - pad - 12 * ws);
            int x0 = (int)(x1 - px);
            int yb = (int)(v->h - pad - 3.2 * lh);
            int t;

            if (px > 12 && px < v->w * 0.4) {
                hud_hline(fb, x0, x1, yb, md, 1.0);
                for (t = 0; t <= 4; t++) {
                    int xt = x0 + (int)(px * t / 4.0);

                    hud_vline(fb, xt, yb - (int)((t % 2 ? 3 : 6) * ws), yb, md,
                              1.0);
                }
                if (ring_step >= 0.05 * SV_AU) {
                    snprintf(lb, sizeof lb, "%g AU", ring_step / SV_AU);
                } else {
                    snprintf(lb, sizeof lb, "%g KM", ring_step / 1000.0);
                }
                snprintf(buf, sizeof buf, "RING %s", lb);
                hud_text(fb, x1 - hud_tw(buf, sm), yb - 8 * ws - 7 * sm, buf,
                         md, sm);
            }
        }
        (void)ring_n;
    }

    /* footer: console or key strip, provenance */
    {
        int fh = (int)(2 * lh + 10 * ws);
        int fy = v->h - pad - fh;

        hud_panel(fb, pad, fy, v->w - 2 * pad, fh);
        if (v->console) {
            snprintf(buf, sizeof buf, "> %s_", v->cmd);
            hud_text(fb, pad + 8 * ws, fy + 6 * ws, buf, amber, sc);
        } else {
            hud_text(fb, pad + 8 * ws, fy + 6 * ws,
                     "/ CONSOLE  SPACE PAUSE  , . WARP  B REVERSE  [ ] SEEK"
                     "  TAB FOCUS  CLICK SELECT  O L G T H  Q QUIT", md, sm);
        }
        hud_text(fb, pad + 8 * ws, fy + 6 * ws + lh,
                 v->msg[0] != '\0' ? v->msg :
                 "IC: JPL HORIZONS (EXTERNAL)  DYNAMICS: N-BODY + 1PN "
                 "(MEASURED)  BELTS/SWARMS: SAMPLED (ILLUSTRATIVE)",
                 v->msg[0] != '\0' ? hi : lo, sm);
    }
}

/* --- render ------------------------------------------------------------- */

#define SV_MAXMARK 65536

static SvMark g_marks[SV_MAXMARK];

void stc_svis_render(StcSVis *v)
{
    SvCtx X;
    const StcSpace *s;
    int i;
    int n;
    double ring_step = 0.0;
    int ring_n = 0;

    if (v == NULL || v->acc == NULL || v->fb.px == NULL) {
        return;
    }
    sv_lut();
    s = &v->sim;
    X.v = v;
    sv_camera(v, &X.C);
    memset(v->acc, 0, (size_t)v->w * (size_t)v->h * 3u * sizeof(float));

    /* faint fixed star field (direction-hashed, deterministic) */
    {
        uint64_t st = 99ULL;
        int k;

        for (k = 0; k < 2400; k++) {
            double zz = 2.0 * sv_u01(&st) - 1.0;
            double ph = 2.0 * SV_PI * sv_u01(&st);
            double rr = sqrt(1.0 - zz * zz);
            double d[3];
            double b = 0.02 + 0.10 * pow(sv_u01(&st), 6.0);
            double vz;

            d[0] = rr * cos(ph);
            d[1] = rr * sin(ph);
            d[2] = zz;
            vz = d[0] * X.C.f[0] + d[1] * X.C.f[1] + d[2] * X.C.f[2];
            if (vz > 0.05) {
                double sx = X.C.cx + X.C.foc * (d[0] * X.C.r[0] + d[1] *
                            X.C.r[1] + d[2] * X.C.r[2]) / vz;
                double sy = X.C.cy - X.C.foc * (d[0] * X.C.u[0] + d[1] *
                            X.C.u[1] + d[2] * X.C.u[2]) / vz;

                sv_dot(v, sx, sy, (float)b, (float)b, (float)(b * 1.05));
            }
        }
    }
    if (v->show_grid) {
        sv_paint_grid(&X, &ring_step, &ring_n);
    } else {
        ring_step = sv_nice(v->dist * 0.22);
    }
    if (v->show_orbits) {
        for (i = 0; i < s->nm; i++) {
            if (!s->alive[i]) {
                continue;
            }
            if (s->kind[i] == STC_SK_ROCK && i != v->select &&
                stc_space_live_massive(s) > 40) {
                continue;
            }
            if (s->kind[i] == STC_SK_STAR && i == sv_heaviest(s)) {
                continue;
            }
            sv_paint_orbit(&X, i, i == v->select ? 0.55f :
                           s->kind[i] == STC_SK_PLANET ? 0.20f :
                           s->kind[i] == STC_SK_BLACKHOLE ? 0.30f : 0.11f);
        }
    }
    if (v->show_tracers) {
        sv_paint_tracers(&X);
    }
    n = sv_collect(v, &X.C, g_marks, SV_MAXMARK);
    sv_paint_bodies(&X, g_marks, n);
    stc_par_for(sv_final_rows, v, v->h, 16);
    if (v->show_labels) {
        sv_paint_labels(&X, g_marks, n);
    }
    sv_paint_hud(&X, ring_step, ring_n);
}

/* --- picking ---------------------------------------------------------- */

void stc_svis_pick(StcSVis *v, int x, int y)
{
    SvCam C;
    const StcSpace *s;
    int i;
    int best = -1;
    double bd = 0.0;

    if (v == NULL) {
        return;
    }
    s = &v->sim;
    sv_camera(v, &C);
    for (i = 0; i < s->nm; i++) {
        SvP p;
        double d;
        double rad;

        if (!s->alive[i]) {
            continue;
        }
        p = sv_proj(&C, s->x[i], s->y[i], s->z[i]);
        if (!p.ok) {
            continue;
        }
        rad = C.foc * s->radius[i] / p.z;
        d = hypot(p.x - x, p.y - y) - rad;
        if (d < 0.0) {
            d = 0.0;
        }
        /* prefer major bodies on near ties */
        d += s->kind[i] == STC_SK_ROCK ? 3.0 : s->kind[i] == STC_SK_MOON
             ? 1.5 : 0.0;
        if (best < 0 || d < bd) {
            best = i;
            bd = d;
        }
    }
    if (best >= 0 && bd < 22.0 * C.ws) {
        v->select = best;
    } else {
        v->select = -1;
    }
}

void stc_svis_orbit(StcSVis *v, double dyaw, double dpitch)
{
    if (v == NULL || !isfinite(dyaw) || !isfinite(dpitch)) {
        return;
    }
    v->yaw -= dyaw;
    v->pitch += dpitch;
    if (v->pitch > 1.5607) {
        v->pitch = 1.5607;
    }
    if (v->pitch < -1.5607) {
        v->pitch = -1.5607;
    }
}

void stc_svis_zoom(StcSVis *v, double factor)
{
    double minr = 2.0e6;

    if (v == NULL || !(factor > 0.0) || !isfinite(factor)) {
        return;
    }
    if (v->focus >= 0 && v->sim.alive[v->focus]) {
        minr = 3.0 * v->sim.radius[v->focus] + 1.0e5;
    }
    v->dist *= factor;
    if (v->dist < minr) {
        v->dist = minr;
    }
    if (v->dist > 2.0e4 * SV_AU) {
        v->dist = 2.0e4 * SV_AU;
    }
}

/* --- console ---------------------------------------------------------- */

static double sv_parse_dur(const char *s, int *ok)
{
    char *end;
    double x = strtod(s, &end);

    *ok = end != s;
    if (!*ok) {
        return 0.0;
    }
    switch (tolower((unsigned char)*end)) {
    case 'y': return x * SV_YEAR;
    case 'd': return x * SV_DAY;
    case 'h': return x * 3600.0;
    case 'm': return x * 60.0;
    default: return x;
    }
}

static void sv_travel(StcSVis *v, double t_target)
{
    StcSpace *s = &v->sim;
    double dt = t_target - s->t;
    long budget = 20000000L;

    if (stc_space_seek(s, t_target) == 0) {
        snprintf(v->msg, sizeof v->msg, "SEEK: RESTORED FROM CHECKPOINT");
        return;
    }
    if (stc_space_advance(s, dt, budget) < 0) {
        snprintf(v->msg, sizeof v->msg,
                 "TRAVEL STOPPED EARLY (STEP BUDGET); RUN AGAIN TO CONTINUE");
    } else {
        char tb[64];

        stc_ephem_format(s->t, tb, (int)sizeof tb);
        snprintf(v->msg, sizeof v->msg, "INTEGRATED TO %s", tb);
    }
}

static void sv_exec(StcSVis *v, const char *line)
{
    StcSpace *s = &v->sim;
    char a0[32] = "";
    char a1[48] = "";
    char a2[32] = "";
    char a3[32] = "";
    char a4[32] = "";
    int na = sscanf(line, "%31s %47s %31s %31s %31s", a0, a1, a2, a3, a4);
    int ok;

    if (na <= 0) {
        return;
    }
    sv_upper(a0);
    if (strcmp(a0, "HELP") == 0 || strcmp(a0, "?") == 0) {
        snprintf(v->msg, sizeof v->msg,
                 "PRESET N | DATE Y-M-D [H:M] | WARP -1Y | BH BODY [MSUN] | "
                 "MASS BODY X | SPAWN EARTH N RING|DISK|CLOUD|SHELL AU | "
                 "ROGUE MSUN KMS AU | FOCUS BODY | DEL BODY | PN ON|OFF | STEP S");
    } else if (strcmp(a0, "PRESET") == 0 || strcmp(a0, "SCENARIO") == 0) {
        int k = atoi(a1);

        if (k >= 0 && k < STC_SVIS_NSCEN) {
            (void)sv_load(v, k);
        } else {
            snprintf(v->msg, sizeof v->msg, "PRESET 0..%d", STC_SVIS_NSCEN - 1);
        }
    } else if (strcmp(a0, "DATE") == 0 || strcmp(a0, "GOTO") == 0) {
        int y = 0;
        int mo = 1;
        int d = 1;
        int hh = 0;
        int mi = 0;

        if (sscanf(a1, "%d-%d-%d", &y, &mo, &d) >= 1) {
            if (na >= 3) {
                (void)sscanf(a2, "%d:%d", &hh, &mi);
            }
            sv_travel(v, stc_ephem_t(stc_ephem_jd_from_cal(y, mo, d, hh, mi,
                                                           0.0)));
        } else {
            snprintf(v->msg, sizeof v->msg, "DATE YYYY-MM-DD [HH:MM]");
        }
    } else if (strcmp(a0, "WARP") == 0) {
        double w = sv_parse_dur(a1, &ok);

        if (ok) {
            v->warp = w;
        }
    } else if (strcmp(a0, "BH") == 0 || strcmp(a0, "MASS") == 0) {
        int i = sv_find(s, a1);
        double x = na >= 3 ? atof(a2) : 0.0;

        if (i < 0) {
            snprintf(v->msg, sizeof v->msg, "NO BODY NAMED %s", a1);
        } else if (strcmp(a0, "BH") == 0) {
            (void)stc_space_make_blackhole(s, i, x > 0.0 ? x * SV_GM_SUN : 0.0,
                                           NULL);
            snprintf(v->msg, sizeof v->msg,
                     "%s IS NOW A %.3g MSUN BLACK HOLE, R_S %.3g KM", a1,
                     s->gm[i] / SV_GM_SUN, s->radius[i] / 1000.0);
            s->e0 = stc_space_energy(s);
        } else if (x > 0.0) {
            s->gm[i] *= x;
            s->e0 = stc_space_energy(s);
            snprintf(v->msg, sizeof v->msg, "%s MASS X %.3g", a1, x);
        }
    } else if (strcmp(a0, "SPAWN") == 0) {
        long n = atol(a2);
        double r = na >= 5 ? atof(a4) : 1.0;
        SvShape sh = SV_RING;
        int c = v->focus >= 0 ? v->focus : sv_heaviest(s);

        sv_upper(a3);
        if (strcmp(a3, "DISK") == 0) {
            sh = SV_DISK;
        } else if (strcmp(a3, "CLOUD") == 0) {
            sh = SV_CLOUD;
        } else if (strcmp(a3, "SHELL") == 0) {
            sh = SV_SHELL;
        }
        if (n <= 0 || n > 5000000L) {
            snprintf(v->msg, sizeof v->msg, "SPAWN EARTH N (1..5000000)");
        } else {
            /* Up to 3000 Earths interact fully; beyond that they are test
             * particles (massless) so the step stays interactive. */
            (void)sv_spawn(s, c, n, sh, r * SV_AU, n <= 3000,
                           (uint64_t)(s->nm * 7919 + s->nt + 1), v->msg,
                           (int)sizeof v->msg);
            s->e0 = stc_space_energy(s);
        }
    } else if (strcmp(a0, "ROGUE") == 0) {
        double m = na >= 2 ? atof(a1) : 1.0;
        double kms = na >= 3 ? atof(a2) : 30.0;
        double b = na >= 4 ? atof(a3) : 10.0;
        double gm = m * SV_GM_SUN;
        double p[3] = {-150.0 * SV_AU, b * SV_AU, 5.0 * SV_AU};
        double vv[3] = {kms * 1000.0, 0.0, 0.0};
        int bh = m >= 3.0;

        if (stc_space_reserve(s, s->nm + 1, s->cap_t) == 0 &&
            stc_space_add(s, bh ? "ROGUE BH" : "ROGUE STAR",
                          bh ? STC_SK_BLACKHOLE : STC_SK_STAR, -1, gm,
                          bh ? 2.0 * gm / (SV_C * SV_C)
                             : 6.957e8 * pow(m, 0.8),
                          p, vv) >= 0) {
            snprintf(v->msg, sizeof v->msg,
                     "ROGUE %.3g MSUN AT 150 AU, %.3g KM/S, AIMED %.3g AU",
                     m, kms, b);
            s->e0 = stc_space_energy(s);
        }
    } else if (strcmp(a0, "FOCUS") == 0) {
        int i = sv_find(s, a1);

        v->focus = i;
        v->select = i;
    } else if (strcmp(a0, "DEL") == 0 || strcmp(a0, "DELETE") == 0) {
        int i = sv_find(s, a1);

        if (i >= 0) {
            s->alive[i] = 0;
            s->e0 = stc_space_energy(s);
        }
    } else if (strcmp(a0, "PN") == 0) {
        sv_upper(a1);
        s->pn = strcmp(a1, "OFF") != 0;
    } else if (strcmp(a0, "STEP") == 0) {
        double h = sv_parse_dur(a1, &ok);

        if (ok && h > 0.0) {
            s->h = h;
        }
    } else {
        snprintf(v->msg, sizeof v->msg, "UNKNOWN COMMAND %s (TRY HELP)", a0);
    }
}

static const double kWarps[] = {1.0, 60.0, 3600.0, 21600.0, 86400.0,
                                7.0 * 86400.0, 30.0 * 86400.0, SV_YEAR,
                                10.0 * SV_YEAR, 100.0 * SV_YEAR};

static void sv_warp_step(StcSVis *v, int dir)
{
    int n = (int)(sizeof kWarps / sizeof kWarps[0]);
    double a = fabs(v->warp);
    double sg = v->warp < 0.0 ? -1.0 : 1.0;
    int k;
    int best = 0;

    for (k = 0; k < n; k++) {
        if (fabs(kWarps[k] - a) < fabs(kWarps[best] - a)) {
            best = k;
        }
    }
    best += dir;
    if (best < 0) {
        best = 0;
    }
    if (best >= n) {
        best = n - 1;
    }
    v->warp = sg * kWarps[best];
}

static void sv_cycle_focus(StcSVis *v)
{
    const StcSpace *s = &v->sim;
    int k;
    int i = v->focus;

    for (k = 0; k < s->nm; k++) {
        i = (i + 1) % s->nm;
        if (s->alive[i] && (s->kind[i] == STC_SK_STAR ||
                            s->kind[i] == STC_SK_PLANET ||
                            s->kind[i] == STC_SK_BLACKHOLE ||
                            s->kind[i] == STC_SK_DWARF)) {
            v->focus = i;
            v->select = i;
            if (v->dist > 2.0 * SV_AU && s->kind[i] != STC_SK_STAR) {
                v->dist = 400.0 * s->radius[i] + 2.0e8;
            }
            return;
        }
    }
    v->focus = -1;
}

int stc_svis_key(StcSVis *v, int key)
{
    if (v == NULL) {
        return 0;
    }
    if (v->console) {
        if (key == STC_UI_KEY_ESC) {
            v->console = 0;
        } else if (key == '\r' || key == '\n') {
            v->console = 0;
            sv_exec(v, v->cmd);
            v->cmd_len = 0;
            v->cmd[0] = '\0';
        } else if (key == 127 || key == 8) {
            if (v->cmd_len > 0) {
                v->cmd[--v->cmd_len] = '\0';
            }
        } else if (key >= 32 && key < 127 && v->cmd_len < STC_SVIS_CMD - 1) {
            v->cmd[v->cmd_len++] = (char)key;
            v->cmd[v->cmd_len] = '\0';
        }
        return 0;
    }
    switch (key) {
    case 'q':
    case 'Q':
    case STC_UI_KEY_ESC:
    case STC_UI_KEY_CTRL_C:
        return 1;
    case '/':
    case ':':
        v->console = 1;
        v->cmd_len = 0;
        v->cmd[0] = '\0';
        break;
    case ' ':
        v->paused = !v->paused;
        break;
    case '.':
        sv_warp_step(v, 1);
        break;
    case ',':
        sv_warp_step(v, -1);
        break;
    case 'b':
    case 'B':
        v->warp = -v->warp;
        break;
    case '[':
    case ']': {
        double dt = fabs(v->warp) * (key == '[' ? -1.0 : 1.0);

        sv_travel(v, v->sim.t + dt);
        break;
    }
    case '\t':
        sv_cycle_focus(v);
        break;
    case 'f':
    case 'F':
        v->focus = v->select;
        break;
    case '0':
        v->focus = -1;
        break;
    case 'o':
    case 'O':
        v->show_orbits = !v->show_orbits;
        break;
    case 'l':
    case 'L':
        v->show_labels = !v->show_labels;
        break;
    case 'g':
    case 'G':
        v->show_grid = !v->show_grid;
        break;
    case 't':
    case 'T':
        v->show_tracers = !v->show_tracers;
        break;
    case 'h':
    case 'H':
        v->show_hud = !v->show_hud;
        break;
    case '+':
    case '=':
        stc_svis_zoom(v, 0.8);
        break;
    case '-':
    case '_':
        stc_svis_zoom(v, 1.25);
        break;
    case 'a':
    case 'A':
    case STC_UI_KEY_LEFT:
        stc_svis_orbit(v, 0.06, 0.0);
        break;
    case 'd':
    case 'D':
    case STC_UI_KEY_RIGHT:
        stc_svis_orbit(v, -0.06, 0.0);
        break;
    case 'w':
    case 'W':
    case STC_UI_KEY_UP:
        stc_svis_orbit(v, 0.0, 0.05);
        break;
    case 's':
    case 'S':
    case STC_UI_KEY_DOWN:
        stc_svis_orbit(v, 0.0, -0.05);
        break;
    default:
        if (key >= '1' && key <= '6') {
            (void)sv_load(v, key - '1');
        }
        break;
    }
    return 0;
}

/* --- headless ------------------------------------------------------------ */

int stc_svis_shot(int scenario, long nframes, int w, int h,
                  const char *ppm_path, const char *dir, long dir_every)
{
    StcSVis *v;
    long f;
    int rc;

    if (nframes < 1 || (ppm_path == NULL && (dir == NULL ||
                                               dir_every <= 0))) {
        return -1;
    }
    v = (StcSVis *)malloc(sizeof *v);
    if (v == NULL) {
        return -2;
    }
    rc = stc_svis_init(v, scenario, w, h);
    if (rc != 0) {
        free(v);
        return rc == -1 ? -1 : -2;
    }
    v->select = sv_find(&v->sim, "Earth");
    for (f = 0; f < nframes; f++) {
        int last = f == nframes - 1;
        int dump = dir != NULL && dir_every > 0 && f % dir_every == 0;

        if (last || dump) {
            stc_svis_render(v);
        }
        if (dump) {
            char path[1024];

            snprintf(path, sizeof path, "%s/space_%05ld.ppm", dir,
                     f / dir_every);
            if (stc_render_write_ppm(&v->fb, path) != 0) {
                stc_svis_free(v);
                free(v);
                return -3;
            }
        }
        if (last && ppm_path != NULL &&
            stc_render_write_ppm(&v->fb, ppm_path) != 0) {
            stc_svis_free(v);
            free(v);
            return -3;
        }
        stc_svis_step(v);
    }
    stc_svis_free(v);
    free(v);
    return 0;
}
