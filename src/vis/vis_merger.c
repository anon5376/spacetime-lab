/* SPACETIME LAB — vis_merger: cinematic binary black hole merger view.
 * HDR software renderer (see stc/vis_merger.h). Pipeline per frame:
 *   sheet + trails -> far/near HDR layers (depth split at the holes)
 *   lens pass: metaball shadows, thin-lens deflection of sky + far layer,
 *              photon ring, near layer on top
 *   bloom pyramid -> ACES tone map + vignette + dither -> 8-bit
 *   HUD: readouts, stage badge, scrolling chirp trace.
 * Row passes run on a small pthread pool; each pixel is independent, so
 * the image does not depend on the thread count.
 */

#include "stc/vis_merger.h"

#include <math.h>
#include <pthread.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "stc/ui.h"
#include "stc/vis.h"

/* --- presets --------------------------------------------------------- */

typedef struct {
    const char *name;
    double q;
    double a0;
    double a_merge;
    double msun;
} MvPreset;

static const MvPreset kPre[STC_MVIS_NPRESET] = {
    {"GW150914-like  36 + 29 Msun", 1.24, 11.0, 3.4, 65.0},
    {"equal mass  30 + 30 Msun", 1.0, 11.0, 3.4, 60.0},
    {"unequal  q = 3  30 + 10 Msun", 3.0, 10.0, 3.4, 40.0},
    {"extreme  q = 8  24 + 3 Msun", 8.0, 8.5, 3.4, 27.0},
};

/* Sheet + look constants (ILLUSTRATIVE visual choices, units of M). */
#define MV_R 170.0        /* sheet half-extent */
#define MV_SP 4.0         /* grid line spacing */
#define MV_DS 1.0         /* sample step along a line */
#define MV_DEPTH 30.0     /* potential-well depth scale */
#define MV_EPS 2.6        /* well softening */
#define MV_RT 0.25        /* retarded-time table step */
#define MV_NRT 1032 /* >= 1.5 * MV_R / MV_RT + 8 */
#define MV_SHADOW (5.196 * 0.62) /* 3 sqrt(3) m, scaled down to read */
#define MV_LENS2 2.1      /* (Einstein radius / shadow radius)^2 */
#define MV_FOV 0.80       /* vertical field of view, rad */
#define MV_GSAT 4.5       /* ripple soft-saturation height */

/* --- worker pool ------------------------------------------------------ */

typedef void (*MvRowFn)(void *ctx, int y0, int y1);

#define MV_MAXT 16
#define MV_CHUNK 8

typedef struct {
    pthread_t th[MV_MAXT];
    int nt;
    int ready;
    pthread_mutex_t mu;
    pthread_cond_t go;
    pthread_cond_t done;
    unsigned long gen;
    int pending;
    MvRowFn fn;
    void *ctx;
    int n;
    int next;
} MvPool;

static MvPool g_pool;

static int mv_take(MvPool *p)
{
    int r;

    pthread_mutex_lock(&p->mu);
    r = p->next;
    p->next += MV_CHUNK;
    pthread_mutex_unlock(&p->mu);
    return r;
}

static void mv_drain(MvPool *p, MvRowFn fn, void *ctx, int n)
{
    for (;;) {
        int y0 = mv_take(p);
        int y1;

        if (y0 >= n) {
            return;
        }
        y1 = y0 + MV_CHUNK < n ? y0 + MV_CHUNK : n;
        fn(ctx, y0, y1);
    }
}

static void *mv_worker(void *arg)
{
    MvPool *p = (MvPool *)arg;
    unsigned long seen = 0;

    for (;;) {
        MvRowFn fn;
        void *ctx;
        int n;

        pthread_mutex_lock(&p->mu);
        while (p->gen == seen) {
            pthread_cond_wait(&p->go, &p->mu);
        }
        seen = p->gen;
        fn = p->fn;
        ctx = p->ctx;
        n = p->n;
        pthread_mutex_unlock(&p->mu);
        mv_drain(p, fn, ctx, n);
        pthread_mutex_lock(&p->mu);
        p->pending--;
        if (p->pending == 0) {
            pthread_cond_signal(&p->done);
        }
        pthread_mutex_unlock(&p->mu);
    }
    return NULL;
}

static void mv_pool_start(void)
{
    long nc;
    int i;

    if (g_pool.ready) {
        return;
    }
    g_pool.ready = 1;
    pthread_mutex_init(&g_pool.mu, NULL);
    pthread_cond_init(&g_pool.go, NULL);
    pthread_cond_init(&g_pool.done, NULL);
    nc = sysconf(_SC_NPROCESSORS_ONLN);
    if (nc < 1) {
        nc = 1;
    }
    g_pool.nt = (int)(nc - 1 < MV_MAXT ? nc - 1 : MV_MAXT);
    for (i = 0; i < g_pool.nt; i++) {
        if (pthread_create(&g_pool.th[i], NULL, mv_worker, &g_pool) != 0) {
            g_pool.nt = i;
            break;
        }
    }
}

/* Run fn over rows [0, n) on the pool (caller participates). */
static void mv_par(MvRowFn fn, void *ctx, int n)
{
    MvPool *p = &g_pool;

    mv_pool_start();
    if (p->nt == 0 || n < 2 * MV_CHUNK) {
        fn(ctx, 0, n);
        return;
    }
    pthread_mutex_lock(&p->mu);
    p->fn = fn;
    p->ctx = ctx;
    p->n = n;
    p->next = 0;
    p->pending = p->nt;
    p->gen++;
    pthread_cond_broadcast(&p->go);
    pthread_mutex_unlock(&p->mu);
    mv_drain(p, fn, ctx, n);
    pthread_mutex_lock(&p->mu);
    while (p->pending > 0) {
        pthread_cond_wait(&p->done, &p->mu);
    }
    pthread_mutex_unlock(&p->mu);
}

/* --- small math -------------------------------------------------------- */

typedef struct {
    double x;
    double y;
    double z;
} V3;

static V3 v3(double x, double y, double z)
{
    V3 r;

    r.x = x;
    r.y = y;
    r.z = z;
    return r;
}

static V3 v3sub(V3 a, V3 b)
{
    return v3(a.x - b.x, a.y - b.y, a.z - b.z);
}

static double v3dot(V3 a, V3 b)
{
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

static V3 v3cross(V3 a, V3 b)
{
    return v3(a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z,
              a.x * b.y - a.y * b.x);
}

static V3 v3norm(V3 a)
{
    double l = sqrt(v3dot(a, a));

    return l > 0.0 ? v3(a.x / l, a.y / l, a.z / l) : a;
}

static double mv_smooth(double e0, double e1, double x)
{
    double t = (x - e0) / (e1 - e0);

    if (t <= 0.0) {
        return 0.0;
    }
    if (t >= 1.0) {
        return 1.0;
    }
    return t * t * (3.0 - 2.0 * t);
}

static uint32_t mv_hash(uint32_t x)
{
    x ^= x >> 16;
    x *= 0x7feb352dU;
    x ^= x >> 15;
    x *= 0x846ca68bU;
    x ^= x >> 16;
    return x;
}

static double mv_u01(uint32_t *s)
{
    *s = mv_hash(*s + 0x9e3779b9U);
    return ((double)(*s >> 8) + 0.5) / 16777216.0;
}

/* --- sky cubemap (generated once, shared) ----------------------------- */

#define MV_SKY_N 1024
#define MV_SKY_SCALE 4096.0f

static uint16_t *g_sky = NULL; /* 6 faces * N * N * 3 */

static double mv_vnoise(double x, double y, double z)
{
    int ix = (int)floor(x);
    int iy = (int)floor(y);
    int iz = (int)floor(z);
    double fx = x - ix;
    double fy = y - iy;
    double fz = z - iz;
    double v[8];
    int k;

    fx = fx * fx * (3.0 - 2.0 * fx);
    fy = fy * fy * (3.0 - 2.0 * fy);
    fz = fz * fz * (3.0 - 2.0 * fz);
    for (k = 0; k < 8; k++) {
        uint32_t h = mv_hash((uint32_t)(ix + (k & 1)) * 73856093U ^
                             (uint32_t)(iy + ((k >> 1) & 1)) * 19349663U ^
                             (uint32_t)(iz + ((k >> 2) & 1)) * 83492791U);
        v[k] = (double)(h & 0xffffU) / 65535.0;
    }
    {
        double a = v[0] + (v[1] - v[0]) * fx;
        double b = v[2] + (v[3] - v[2]) * fx;
        double c = v[4] + (v[5] - v[4]) * fx;
        double d = v[6] + (v[7] - v[6]) * fx;
        double e = a + (b - a) * fy;
        double f = c + (d - c) * fy;

        return e + (f - e) * fz;
    }
}

static double mv_fbm(double x, double y, double z, int oct)
{
    double s = 0.0;
    double a = 0.5;
    double n = 0.0;
    int i;

    for (i = 0; i < oct; i++) {
        s += a * mv_vnoise(x, y, z);
        n += a;
        a *= 0.5;
        x = x * 2.03 + 11.7;
        y = y * 2.03 - 5.3;
        z = z * 2.03 + 3.1;
    }
    return s / n;
}

static V3 mv_face_dir(int f, double u, double v)
{
    switch (f) {
    case 0: return v3(1.0, v, u);
    case 1: return v3(-1.0, v, u);
    case 2: return v3(u, 1.0, v);
    case 3: return v3(u, -1.0, v);
    case 4: return v3(u, v, 1.0);
    default: return v3(u, v, -1.0);
    }
}

/* Direction -> face, texel coords (continuous, texel centers at +0.5). */
static int mv_face_uv(double x, double y, double z, double *tu, double *tv)
{
    double ax = fabs(x);
    double ay = fabs(y);
    double az = fabs(z);
    double u;
    double v;
    int f;

    if (ax >= ay && ax >= az) {
        f = x > 0.0 ? 0 : 1;
        u = z / ax;
        v = y / ax;
    } else if (ay >= az) {
        f = y > 0.0 ? 2 : 3;
        u = x / ay;
        v = z / ay;
    } else {
        f = z > 0.0 ? 4 : 5;
        u = x / az;
        v = y / az;
    }
    *tu = (u + 1.0) * 0.5 * MV_SKY_N;
    *tv = (v + 1.0) * 0.5 * MV_SKY_N;
    return f;
}

static void mv_sky_rows(void *ctx, int r0, int r1)
{
    int row;

    (void)ctx;
    for (row = r0; row < r1; row++) {
        int f = row / MV_SKY_N;
        int j = row % MV_SKY_N;
        int i;

        for (i = 0; i < MV_SKY_N; i++) {
            double u = (i + 0.5) / MV_SKY_N * 2.0 - 1.0;
            double v = (j + 0.5) / MV_SKY_N * 2.0 - 1.0;
            V3 n = v3norm(mv_face_dir(f, u, v));
            V3 g = v3norm(v3(0.32, 0.78, -0.54));
            double d = v3dot(n, g);
            double band = exp(-d * d / (2.0 * 0.15 * 0.15));
            double f1 = mv_fbm(n.x * 2.1 + 5.0, n.y * 2.1, n.z * 2.1, 5);
            double f2 = mv_fbm(n.x * 5.3, n.y * 5.3 + 9.0, n.z * 5.3, 5);
            double f3 = mv_fbm(n.x * 1.3 - 4.0, n.y * 1.3, n.z * 1.3 + 2.0, 4);
            double dust = 1.0 - 0.85 * band * mv_smooth(0.48, 0.66, f2);
            double neb = pow(mv_smooth(0.35, 0.85, f3), 2.2);
            double r = 0.0010;
            double gg = 0.0010;
            double b = 0.0012;
            uint16_t *px;

            r += band * (0.006 + 0.010 * f1) * dust;
            gg += band * (0.006 + 0.010 * f1) * dust;
            b += band * (0.007 + 0.011 * f1) * dust;
            (void)neb;
            px = g_sky + ((size_t)row * MV_SKY_N + (size_t)i) * 3u;
            px[0] = (uint16_t)(r * MV_SKY_SCALE);
            px[1] = (uint16_t)(gg * MV_SKY_SCALE);
            px[2] = (uint16_t)(b * MV_SKY_SCALE);
        }
    }
}

static void mv_sky_add(int f, int x, int y, double r, double g, double b)
{
    uint16_t *px;
    double v[3];
    int k;

    if (x < 0 || y < 0 || x >= MV_SKY_N || y >= MV_SKY_N) {
        return;
    }
    px = g_sky + (((size_t)f * MV_SKY_N + (size_t)y) * MV_SKY_N +
                  (size_t)x) * 3u;
    v[0] = r;
    v[1] = g;
    v[2] = b;
    for (k = 0; k < 3; k++) {
        double o = px[k] + v[k] * MV_SKY_SCALE;

        px[k] = (uint16_t)(o > 65535.0 ? 65535.0 : o);
    }
}

static int mv_sky_build(void)
{
    uint32_t s = 20260927U;
    int i;

    if (g_sky != NULL) {
        return 0;
    }
    g_sky = (uint16_t *)calloc((size_t)6 * MV_SKY_N * MV_SKY_N * 3u,
                               sizeof(uint16_t));
    if (g_sky == NULL) {
        return -1;
    }
    mv_par(mv_sky_rows, NULL, 6 * MV_SKY_N);
    /* Stars: power-law brightness, crowded toward the band. */
    for (i = 0; i < 60000; i++) {
        double z = mv_u01(&s) * 2.0 - 1.0;
        double ph = mv_u01(&s) * 6.283185307179586;
        double rr = sqrt(1.0 - z * z);
        V3 n = v3(rr * cos(ph), rr * sin(ph), z);
        V3 g = v3norm(v3(0.32, 0.78, -0.54));
        double d = v3dot(n, g);
        double band = exp(-d * d / (2.0 * 0.2 * 0.2));
        double u = mv_u01(&s);
        double br;
        double t = mv_u01(&s);
        double cr;
        double cg;
        double cb;
        double tu;
        double tv;
        double sig;
        int f;
        int x0;
        int y0;
        int dx;
        int dy;
        int rad;
        double wsum = 0.0;

        if (mv_u01(&s) > 0.30 + 0.70 * band) {
            continue;
        }
        br = 0.010 * pow(u, -0.80);
        if (br > 30.0) {
            br = 30.0;
        }
        if (t < 0.22) {
            cr = 0.62; cg = 0.74; cb = 1.00;
        } else if (t < 0.60) {
            cr = 0.95; cg = 0.96; cb = 1.00;
        } else if (t < 0.85) {
            cr = 1.00; cg = 0.88; cb = 0.68;
        } else {
            cr = 1.00; cg = 0.66; cb = 0.42;
        }
        f = mv_face_uv(n.x, n.y, n.z, &tu, &tv);
        sig = br > 1.0 ? 0.75 : 0.55;
        rad = br > 1.0 ? 2 : 1;
        x0 = (int)floor(tu);
        y0 = (int)floor(tv);
        for (dy = -rad; dy <= rad; dy++) {
            for (dx = -rad; dx <= rad; dx++) {
                double ex = x0 + dx + 0.5 - tu;
                double ey = y0 + dy + 0.5 - tv;

                wsum += exp(-(ex * ex + ey * ey) / (2.0 * sig * sig));
            }
        }
        for (dy = -rad; dy <= rad; dy++) {
            for (dx = -rad; dx <= rad; dx++) {
                double ex = x0 + dx + 0.5 - tu;
                double ey = y0 + dy + 0.5 - tv;
                double w = exp(-(ex * ex + ey * ey) / (2.0 * sig * sig)) /
                           wsum * br;

                mv_sky_add(f, x0 + dx, y0 + dy, cr * w, cg * w, cb * w);
            }
        }
    }
    return 0;
}

static void mv_sky_sample(double x, double y, double z, float *out)
{
    double tu;
    double tv;
    int f = mv_face_uv(x, y, z, &tu, &tv);
    int x0;
    int y0;
    double fx;
    double fy;
    const uint16_t *base;
    const uint16_t *p00;
    const uint16_t *p10;
    const uint16_t *p01;
    const uint16_t *p11;
    int k;
    int x1;
    int y1;

    tu -= 0.5;
    tv -= 0.5;
    if (tu < 0.0) {
        tu = 0.0;
    }
    if (tv < 0.0) {
        tv = 0.0;
    }
    if (tu > MV_SKY_N - 1.001) {
        tu = MV_SKY_N - 1.001;
    }
    if (tv > MV_SKY_N - 1.001) {
        tv = MV_SKY_N - 1.001;
    }
    x0 = (int)tu;
    y0 = (int)tv;
    x1 = x0 + 1;
    y1 = y0 + 1;
    fx = tu - x0;
    fy = tv - y0;
    base = g_sky + (size_t)f * MV_SKY_N * MV_SKY_N * 3u;
    p00 = base + ((size_t)y0 * MV_SKY_N + (size_t)x0) * 3u;
    p10 = base + ((size_t)y0 * MV_SKY_N + (size_t)x1) * 3u;
    p01 = base + ((size_t)y1 * MV_SKY_N + (size_t)x0) * 3u;
    p11 = base + ((size_t)y1 * MV_SKY_N + (size_t)x1) * 3u;
    for (k = 0; k < 3; k++) {
        double a = p00[k] + (p10[k] - (double)p00[k]) * fx;
        double b = p01[k] + (p11[k] - (double)p01[k]) * fx;

        out[k] = (float)((a + (b - a) * fy) / MV_SKY_SCALE);
    }
}

/* --- tone map LUT ------------------------------------------------------ */

#define MV_LUT_N 4096
#define MV_LUT_MAX 24.0

static uint8_t g_lut[MV_LUT_N + 1];
static int g_lut_ok = 0;

static void mv_lut_build(void)
{
    int i;

    if (g_lut_ok) {
        return;
    }
    for (i = 0; i <= MV_LUT_N; i++) {
        double s = (double)i / MV_LUT_N;
        double x = s * s * MV_LUT_MAX;
        double a = x * (2.51 * x + 0.03) / (x * (2.43 * x + 0.59) + 0.14);
        double c;

        if (a < 0.0) {
            a = 0.0;
        }
        if (a > 1.0) {
            a = 1.0;
        }
        c = a <= 0.0031308 ? 12.92 * a : 1.055 * pow(a, 1.0 / 2.4) - 0.055;
        g_lut[i] = (uint8_t)(c * 255.0 + 0.5);
    }
    g_lut_ok = 1;
}

static int mv_tm(float x)
{
    double s;

    if (!(x > 0.0f)) {
        return 0;
    }
    s = sqrt((double)x / MV_LUT_MAX);
    if (s >= 1.0) {
        return g_lut[MV_LUT_N];
    }
    return g_lut[(int)(s * MV_LUT_N + 0.5)];
}

/* --- per-frame scene description ---------------------------------------- */

typedef struct {
    StcMVis *v;
    StcBbhState st;
    /* camera */
    V3 eye;
    V3 fwd;
    V3 right;
    V3 up;
    double foc;       /* focal length, px */
    double cx;
    double cy;
    double zcom;      /* view depth of the system center */
    /* holes */
    int nh;
    double hm[2];     /* well mass */
    double hx[2];     /* orbital-plane positions */
    double hz[2];
    V3 hw[2];         /* world centers */
    double sx[2];     /* screen centers */
    double sy[2];
    double r2[2];     /* shadow radius^2, px^2 */
    double e2[2];     /* Einstein radius^2, px^2 */
    double ring_eps;  /* ringdown shadow wobble amplitude */
    double ring_ph;   /* ringdown wobble phase */
    double glow;      /* photon ring brightness */
    /* retarded strain table */
    float tc[MV_NRT];
    float ts[MV_NRT];
    double fade;
} MvFrame;

static double mv_rate(const StcMVis *v)
{
    double d = fabs(v->t - v->bbh.t_peak);

    return 1.5 * (0.24 + 0.76 * mv_smooth(12.0, 170.0, d));
}

static double mv_t_end(const StcMVis *v)
{
    return v->bbh.t_peak + MV_R + 60.0;
}

static void mv_frame_setup(StcMVis *v, MvFrame *F)
{
    const StcBbh *b = &v->bbh;
    double bump;
    double dist;
    double pitch;
    V3 tgt;
    int i;
    double wg;

    memset(F, 0, sizeof(*F));
    F->v = v;
    stc_bbh_eval(b, v->t, &F->st);

    /* Camera: slow dolly-in around the merger. */
    bump = exp(-pow((v->t - b->t_peak) / 160.0, 2.0));
    dist = v->dist * (1.0 - 0.22 * bump);
    pitch = v->pitch + 0.06 * bump;
    tgt = v3(0.0, -4.0, 0.0);
    F->eye = v3(tgt.x + dist * cos(pitch) * sin(v->yaw),
                tgt.y + dist * sin(pitch),
                tgt.z + dist * cos(pitch) * cos(v->yaw));
    F->fwd = v3norm(v3sub(tgt, F->eye));
    F->right = v3norm(v3cross(F->fwd, v3(0.0, 1.0, 0.0)));
    F->up = v3cross(F->right, F->fwd);
    F->foc = 0.5 * v->h / tan(0.5 * MV_FOV);
    F->cx = 0.5 * v->w;
    F->cy = 0.5 * v->h;
    F->zcom = v3dot(v3sub(v3(0.0, 0.0, 0.0), F->eye), F->fwd);

    /* Retarded strain table: r h+(t - r) split as C cos2az + S sin2az,
     * premultiplied by the visual gain, near-zone taper and falloff. */
    wg = 75.0 * v->wave_gain;
    for (i = 0; i < MV_NRT; i++) {
        double r = i * MV_RT;
        StcBbhState s;
        double k;

        stc_bbh_eval(b, v->t - r, &s);
        k = wg * (1.0 - exp(-(r * r) / (14.0 * 14.0))) /
            pow(r + 12.0, 0.45);
        F->tc[i] = (float)(k * 0.5 * s.amp * cos(2.0 * s.phi));
        F->ts[i] = (float)(k * 0.5 * s.amp * sin(2.0 * s.phi));
    }

    /* Holes. Mass budget m_i * mass; shadows grow into the remnant's. */
    F->nh = 2;
    {
        double m[2];
        double rf = MV_SHADOW * b->mf * (1.0 - 0.06 * b->af);
        double fu = F->st.fuse;

        m[0] = b->m1;
        m[1] = b->m2;
        F->hx[0] = F->st.x1[0];
        F->hz[0] = F->st.x1[1];
        F->hx[1] = F->st.x2[0];
        F->hz[1] = F->st.x2[1];
        for (i = 0; i < 2; i++) {
            double r0 = MV_SHADOW * m[i];
            double rw2 = (1.0 - fu) * r0 * r0 + fu * m[i] * rf * rf;
            double rw = sqrt(rw2);
            double y;
            V3 d;
            double zv;

            F->hm[i] = m[i] * F->st.mass;
            y = 0.0;
            {
                int j;

                for (j = 0; j < 2; j++) {
                    double ddx = F->hx[i] - F->st.x1[0] * (j == 0) -
                                 F->st.x2[0] * (j == 1);
                    double ddz = F->hz[i] - F->st.x1[1] * (j == 0) -
                                 F->st.x2[1] * (j == 1);

                    y -= MV_DEPTH * m[j] * F->st.mass /
                         sqrt(ddx * ddx + ddz * ddz + MV_EPS * MV_EPS);
                }
            }
            F->hw[i] = v3(F->hx[i], y + 0.55 * rw + 0.5, F->hz[i]);
            d = v3sub(F->hw[i], F->eye);
            zv = v3dot(d, F->fwd);
            if (zv < 1.0) {
                zv = 1.0;
            }
            F->sx[i] = F->cx + F->foc * v3dot(d, F->right) / zv;
            F->sy[i] = F->cy - F->foc * v3dot(d, F->up) / zv;
            F->r2[i] = pow(F->foc * rw / zv, 2.0);
            if (i == 0 || zv - 1.2 * rw < F->zcom) {
                F->zcom = zv - 1.2 * rw;
            }
            F->e2[i] = v->lens_on ? MV_LENS2 * F->r2[i] : 0.0;
        }
    }
    if (F->st.stage == STC_BBH_RINGDOWN) {
        double u = v->t - b->t_peak;

        F->ring_eps = 0.16 * exp(-u / b->tau_qnm);
        F->ring_ph = 2.0 * F->st.phi;
    } else if (F->st.stage == STC_BBH_PLUNGE) {
        F->ring_eps = 0.10 * F->st.fuse;
        F->ring_ph = 2.0 * F->st.phi;
    }
    F->glow = 1.0 + 0.8 * exp(-pow((v->t - b->t_peak) / 25.0, 2.0));

    /* Loop fade in/out. */
    F->fade = mv_smooth(0.0, 45.0, v->t) *
              (1.0 - mv_smooth(mv_t_end(v) - 70.0, mv_t_end(v), v->t));
}

/* Sheet height (well + retarded ripple) and ripple value at (x, z). */
static double mv_height(const MvFrame *F, double x, double z, double *gw,
                        double *well)
{
    double w = 0.0;
    double r;
    double g = 0.0;
    int i;

    for (i = 0; i < F->nh; i++) {
        double dx = x - F->hx[i];
        double dz = z - F->hz[i];

        w += F->hm[i] / sqrt(dx * dx + dz * dz + MV_EPS * MV_EPS);
    }
    r = sqrt(x * x + z * z);
    if (r > 1e-6) {
        double fi = r / MV_RT;
        int k = (int)fi;
        double c2 = (x * x - z * z) / (r * r);
        double s2 = 2.0 * x * z / (r * r);

        if (k < MV_NRT - 1) {
            double fr = fi - k;
            double tc = F->tc[k] + (F->tc[k + 1] - F->tc[k]) * fr;
            double ts = F->ts[k] + (F->ts[k + 1] - F->ts[k]) * fr;

            g = tc * c2 + ts * s2;
            /* soft-saturate so the ringdown crest stays a ripple */
            g = MV_GSAT * tanh(g / MV_GSAT);
        }
    }
    *gw = g;
    *well = w;
    return -MV_DEPTH * w + g;
}

/* --- HDR line splat ------------------------------------------------------- */

static void mv_add(float *buf, int w, int h, int x, int y, float r, float g,
                   float b)
{
    float *p;

    if (x < 0 || y < 0 || x >= w || y >= h) {
        return;
    }
    p = buf + ((size_t)y * (size_t)w + (size_t)x) * 3u;
    p[0] += r;
    p[1] += g;
    p[2] += b;
}

/* Anti-aliased line of half-width hw px (continuous pixel coords, pixel i
 * spans [i, i+1)), constant color, additive. */
static void mv_line(float *buf, int w, int h, double x0, double y0,
                    double x1, double y1, double hw, float r, float g,
                    float b)
{
    int steep = fabs(y1 - y0) > fabs(x1 - x0);
    double t;
    double grad;
    double lenk;
    int i0;
    int i1;
    int i;
    int lim_major;
    int lim_minor;

    if (steep) {
        t = x0; x0 = y0; y0 = t;
        t = x1; x1 = y1; y1 = t;
    }
    if (x0 > x1) {
        t = x0; x0 = x1; x1 = t;
        t = y0; y0 = y1; y1 = t;
    }
    lim_major = steep ? h : w;
    lim_minor = steep ? w : h;
    if (x1 < 0.0 || x0 > lim_major) {
        return;
    }
    if ((y0 < -4.0 && y1 < -4.0) ||
        (y0 > lim_minor + 4.0 && y1 > lim_minor + 4.0)) {
        return;
    }
    grad = x1 - x0 > 1e-9 ? (y1 - y0) / (x1 - x0) : 0.0;
    lenk = sqrt(1.0 + grad * grad);
    r *= (float)lenk;
    g *= (float)lenk;
    b *= (float)lenk;
    i0 = (int)ceil(x0 - 0.5);
    i1 = (int)floor(x1 - 0.5);
    if (i0 < 0) {
        i0 = 0;
    }
    if (i1 > lim_major - 1) {
        i1 = lim_major - 1;
    }
    if (i1 < i0) {
        /* Sub-pixel segment: deposit its length at the midpoint. */
        double xm = 0.5 * (x0 + x1);
        double ym = 0.5 * (y0 + y1);
        float k = (float)(x1 - x0);
        int xi = (int)floor(xm);
        int yi = (int)floor(ym);

        if (steep) {
            mv_add(buf, w, h, yi, xi, r * k, g * k, b * k);
        } else {
            mv_add(buf, w, h, xi, yi, r * k, g * k, b * k);
        }
        return;
    }
    for (i = i0; i <= i1; i++) {
        double yc = y0 + grad * (i + 0.5 - x0);
        double lo = yc - hw;
        double hi = yc + hw;
        int j0 = (int)floor(lo);
        int j1 = (int)floor(hi);
        int j;
        double inv = 1.0 / (2.0 * hw);

        for (j = j0; j <= j1; j++) {
            double a = (j + 1 < hi ? j + 1 : hi) - (j > lo ? j : lo);
            float c;

            if (a <= 0.0) {
                continue;
            }
            c = (float)(a * inv);
            if (steep) {
                mv_add(buf, w, h, j, i, r * c, g * c, b * c);
            } else {
                mv_add(buf, w, h, i, j, r * c, g * c, b * c);
            }
        }
    }
}

typedef struct {
    double sx;
    double sy;
    double zv;
} MvP;

static MvP mv_proj(const MvFrame *F, V3 p)
{
    V3 d = v3sub(p, F->eye);
    MvP o;

    o.zv = v3dot(d, F->fwd);
    o.sx = v3dot(d, F->right);
    o.sy = v3dot(d, F->up);
    return o;
}

/* 3D segment -> near-clip -> screen -> layer by depth. a/b hold view
 * coords (sx, sy = view x, y before division). */
static void mv_seg3(const MvFrame *F, MvP a, MvP b, double hw, float r,
                    float g, float bl)
{
    const double zn = 2.0;
    StcMVis *v = F->v;
    float *buf;
    double ax;
    double ay;
    double bx;
    double by;

    if (a.zv < zn && b.zv < zn) {
        return;
    }
    if (a.zv < zn || b.zv < zn) {
        double t = (zn - a.zv) / (b.zv - a.zv);
        MvP c;

        c.sx = a.sx + (b.sx - a.sx) * t;
        c.sy = a.sy + (b.sy - a.sy) * t;
        c.zv = zn;
        if (a.zv < zn) {
            a = c;
        } else {
            b = c;
        }
    }
    ax = F->cx + F->foc * a.sx / a.zv;
    ay = F->cy - F->foc * a.sy / a.zv;
    bx = F->cx + F->foc * b.sx / b.zv;
    by = F->cy - F->foc * b.sy / b.zv;
    if (!(fabs(ax) < 1e7 && fabs(ay) < 1e7 && fabs(bx) < 1e7 &&
          fabs(by) < 1e7)) {
        return; /* near-plane blow-up; mv_line clips the rest */
    }
    buf = 0.5 * (a.zv + b.zv) > F->zcom ? v->far : v->nearl;
    mv_line(buf, v->w, v->h, ax, ay, bx, by, hw, r, g, bl);
}

/* Grid color from well depth + ripple (ILLUSTRATIVE palette). */
static void mv_grid_col(double x, double z, double gw,
                        double well, double zv, float *c)
{
    double rr = sqrt(x * x + z * z);
    double edge = 1.0 - mv_smooth(0.55 * MV_R, MV_R, rr);
    double wd = mv_smooth(0.02, 0.55, well);
    double crest = tanh(gw / 1.1);
    double ib;
    double cr;
    double cg;
    double cb;
    double depth_fade = 1.0 / (1.0 + pow(zv / 260.0, 2.0));

    /* mission-ops monochrome: neutral grey, brighter on ripple crests and
     * down the wells */
    cr = 0.86;
    cg = 0.88;
    cb = 0.92;
    ib = 0.16 * (1.0 + 1.2 * (crest > 0.0 ? crest : 0.6 * crest)) *
         (1.0 + 1.2 * wd);
    if (ib < 0.05) {
        ib = 0.05;
    }
    ib *= edge * depth_fade;
    c[0] = (float)(cr * ib);
    c[1] = (float)(cg * ib);
    c[2] = (float)(cb * ib);
}

static void mv_paint_grid(const MvFrame *F)
{
    int nl = (int)(2.0 * MV_R / MV_SP + 0.5) + 1;
    int ns = (int)(2.0 * MV_R / MV_DS + 0.5) + 1;
    double hw = F->v->h > 1000 ? 0.85 : 0.6;
    int dir;
    int l;
    int s;

    for (dir = 0; dir < 2; dir++) {
        for (l = 0; l < nl; l++) {
            double c = -MV_R + l * MV_SP;
            MvP prev;
            float pc[3];
            int have = 0;

            for (s = 0; s < ns; s++) {
                double a = -MV_R + s * MV_DS;
                double x = dir == 0 ? a : c;
                double z = dir == 0 ? c : a;
                double gw;
                double well;
                double y;
                MvP p;
                float col[3];

                if (x * x + z * z > MV_R * MV_R * 1.02) {
                    have = 0;
                    continue;
                }
                y = mv_height(F, x, z, &gw, &well);
                p = mv_proj(F, v3(x, y, z));
                mv_grid_col(x, z, gw, well, p.zv, col);
                if (have) {
                    mv_seg3(F, prev, p, hw, 0.5f * (pc[0] + col[0]),
                            0.5f * (pc[1] + col[1]),
                            0.5f * (pc[2] + col[2]));
                }
                prev = p;
                pc[0] = col[0];
                pc[1] = col[1];
                pc[2] = col[2];
                have = 1;
            }
        }
    }
}

/* Orbital trails: last ~1.3 orbits, fading, per hole. */
static void mv_paint_trails(const MvFrame *F)
{
    const StcMVis *v = F->v;
    const StcBbh *b = &v->bbh;
    double span;
    int n = 160;
    int i;
    int h;
    static const float kc[2][3] = {{0.80f, 0.80f, 0.82f},
                                   {0.55f, 0.56f, 0.60f}};

    if (F->st.stage == STC_BBH_RINGDOWN || F->st.omega <= 0.0) {
        return;
    }
    span = 1.3 * 6.283185307179586 / F->st.omega;
    if (span > v->t + 400.0) {
        span = v->t + 400.0;
    }
    for (h = 0; h < 2; h++) {
        MvP prev;
        int have = 0;

        for (i = 0; i <= n; i++) {
            double tt = v->t - span * (1.0 - (double)i / n);
            StcBbhState s;
            double x;
            double z;
            double gw;
            double well;
            double y;
            MvP p;
            double a = (double)i / n;
            float k;

            stc_bbh_eval(b, tt, &s);
            x = h == 0 ? s.x1[0] : s.x2[0];
            z = h == 0 ? s.x1[1] : s.x2[1];
            y = mv_height(F, x, z, &gw, &well) + 0.6;
            p = mv_proj(F, v3(x, y, z));
            k = (float)(0.7 * a * a * (1.0 - F->st.fuse));
            if (have) {
                mv_seg3(F, prev, p, v->h > 1000 ? 1.2 : 0.9, kc[h][0] * k,
                        kc[h][1] * k, kc[h][2] * k);
            }
            prev = p;
            have = 1;
        }
        (void)b;
    }
}

/* --- lens + composite pass ----------------------------------------------- */

static void mv_far_sample(const StcMVis *v, double x, double y, float *o)
{
    int x0;
    int y0;
    double fx;
    double fy;
    int k;
    const float *p00;
    const float *p10;
    const float *p01;
    const float *p11;

    x -= 0.5;
    y -= 0.5;
    if (x < 0.0 || y < 0.0 || x >= v->w - 1 || y >= v->h - 1) {
        o[0] = o[1] = o[2] = 0.0f;
        return;
    }
    x0 = (int)x;
    y0 = (int)y;
    fx = x - x0;
    fy = y - y0;
    p00 = v->far + ((size_t)y0 * (size_t)v->w + (size_t)x0) * 3u;
    p10 = p00 + 3;
    p01 = p00 + (size_t)v->w * 3u;
    p11 = p01 + 3;
    for (k = 0; k < 3; k++) {
        double a = p00[k] + (p10[k] - p00[k]) * fx;
        double b = p01[k] + (p11[k] - p01[k]) * fx;

        o[k] = (float)(a + (b - a) * fy);
    }
}

static void mv_lens_rows(void *ctx, int y0, int y1)
{
    const MvFrame *F = (const MvFrame *)ctx;
    const StcMVis *v = F->v;
    int y;

    for (y = y0; y < y1; y++) {
        int x;

        for (x = 0; x < v->w; x++) {
            double px = x + 0.5;
            double py = y + 0.5;
            double S = 0.0;
            double bx = px;
            double by = py;
            size_t o = ((size_t)y * (size_t)v->w + (size_t)x) * 3u;
            float *dst = v->hdr + o;
            const float *nl = v->nearl + o;
            int i;

            for (i = 0; i < F->nh; i++) {
                double dx = px - F->sx[i];
                double dy = py - F->sy[i];
                double d2 = dx * dx + dy * dy + 1e-6;
                double r2 = F->r2[i];

                if (F->ring_eps != 0.0) {
                    double c2 = (dx * dx - dy * dy) / d2;
                    double s2 = 2.0 * dx * dy / d2;

                    r2 *= 1.0 + F->ring_eps *
                                (c2 * cos(F->ring_ph) + s2 * sin(F->ring_ph));
                }
                S += r2 / d2;
                {
                    /* thin-lens 1/d deflection, tapered beyond ~4 shadow
                     * radii so distant lines stay put (ILLUSTRATIVE) */
                    double k = F->e2[i] / d2 /
                               (1.0 + d2 / (16.0 * F->r2[i] + 1e-6));

                    bx -= k * dx;
                    by -= k * dy;
                }
            }
            if (S >= 1.0) {
                dst[0] = dst[1] = dst[2] = 0.0f;
            } else {
                float sky[3];
                float fr[3];
                V3 d;
                double q;
                double ring;

                d.x = F->fwd.x * F->foc + (bx - F->cx) * F->right.x -
                      (by - F->cy) * F->up.x;
                d.y = F->fwd.y * F->foc + (bx - F->cx) * F->right.y -
                      (by - F->cy) * F->up.y;
                d.z = F->fwd.z * F->foc + (bx - F->cx) * F->right.z -
                      (by - F->cy) * F->up.z;
                mv_sky_sample(d.x, d.y, d.z, sky);
                mv_far_sample(v, bx, by, fr);
                q = 1.0 / sqrt(S) - 1.0;
                ring = F->glow * (0.9 * exp(-q * q / (0.012 * 0.012)) +
                                  0.06 * exp(-q / 0.08));
                dst[0] = sky[0] + fr[0] + (float)(ring * 0.95);
                dst[1] = sky[1] + fr[1] + (float)(ring * 0.95);
                dst[2] = sky[2] + fr[2] + (float)(ring * 0.98);
            }
            dst[0] += nl[0];
            dst[1] += nl[1];
            dst[2] += nl[2];
        }
    }
}

/* --- bloom ------------------------------------------------------------------ */

typedef struct {
    const float *src;
    int sw;
    int sh;
    float *dst;
    int dw;
    int dh;
    int thresh;
} MvDown;

static const float kTap[4] = {1.0f / 8.0f, 3.0f / 8.0f, 3.0f / 8.0f,
                              1.0f / 8.0f};

static void mv_down_rows(void *ctx, int y0, int y1)
{
    const MvDown *D = (const MvDown *)ctx;
    int y;

    for (y = y0; y < y1; y++) {
        int x;

        for (x = 0; x < D->dw; x++) {
            float acc[3] = {0.0f, 0.0f, 0.0f};
            int j;
            int i;
            float *o = D->dst + ((size_t)y * (size_t)D->dw + (size_t)x) * 3u;

            for (j = 0; j < 4; j++) {
                int sy = 2 * y - 1 + j;

                if (sy < 0) {
                    sy = 0;
                }
                if (sy >= D->sh) {
                    sy = D->sh - 1;
                }
                for (i = 0; i < 4; i++) {
                    int sx = 2 * x - 1 + i;
                    const float *p;
                    float wgt = kTap[i] * kTap[j];

                    if (sx < 0) {
                        sx = 0;
                    }
                    if (sx >= D->sw) {
                        sx = D->sw - 1;
                    }
                    p = D->src + ((size_t)sy * (size_t)D->sw + (size_t)sx) * 3u;
                    acc[0] += p[0] * wgt;
                    acc[1] += p[1] * wgt;
                    acc[2] += p[2] * wgt;
                }
            }
            if (D->thresh) {
                float l = 0.2126f * acc[0] + 0.7152f * acc[1] +
                          0.0722f * acc[2];
                float k = l > 1e-6f ? (l - 0.35f) / l : 0.0f;

                if (k < 0.0f) {
                    k = 0.0f;
                }
                acc[0] *= k;
                acc[1] *= k;
                acc[2] *= k;
            }
            o[0] = acc[0];
            o[1] = acc[1];
            o[2] = acc[2];
        }
    }
}

static void mv_bilin(const float *s, int w, int h, double x, double y,
                     float *o)
{
    int x0;
    int y0;
    int x1;
    int y1;
    double fx;
    double fy;
    int k;

    x -= 0.5;
    y -= 0.5;
    if (x < 0.0) {
        x = 0.0;
    }
    if (y < 0.0) {
        y = 0.0;
    }
    if (x > w - 1) {
        x = w - 1;
    }
    if (y > h - 1) {
        y = h - 1;
    }
    x0 = (int)x;
    y0 = (int)y;
    x1 = x0 + 1 < w ? x0 + 1 : x0;
    y1 = y0 + 1 < h ? y0 + 1 : y0;
    fx = x - x0;
    fy = y - y0;
    for (k = 0; k < 3; k++) {
        double a = s[((size_t)y0 * w + x0) * 3 + k] +
                   (s[((size_t)y0 * w + x1) * 3 + k] -
                    s[((size_t)y0 * w + x0) * 3 + k]) * fx;
        double b = s[((size_t)y1 * w + x0) * 3 + k] +
                   (s[((size_t)y1 * w + x1) * 3 + k] -
                    s[((size_t)y1 * w + x0) * 3 + k]) * fx;

        o[k] = (float)(a + (b - a) * fy);
    }
}

typedef struct {
    const float *src;
    int sw;
    int sh;
    float *dst;
    int dw;
    int dh;
} MvUp;

static void mv_up_rows(void *ctx, int y0, int y1)
{
    const MvUp *U = (const MvUp *)ctx;
    int y;
    double kx = (double)U->sw / U->dw;
    double ky = (double)U->sh / U->dh;

    for (y = y0; y < y1; y++) {
        int x;

        for (x = 0; x < U->dw; x++) {
            float acc[3] = {0.0f, 0.0f, 0.0f};
            float t[3];
            int j;
            int i;
            float *o = U->dst + ((size_t)y * (size_t)U->dw + (size_t)x) * 3u;
            double sx = (x + 0.5) * kx;
            double sy = (y + 0.5) * ky;

            for (j = -1; j <= 1; j++) {
                for (i = -1; i <= 1; i++) {
                    float wgt = (float)((2 - abs(i)) * (2 - abs(j))) / 16.0f;

                    mv_bilin(U->src, U->sw, U->sh, sx + i, sy + j, t);
                    acc[0] += t[0] * wgt;
                    acc[1] += t[1] * wgt;
                    acc[2] += t[2] * wgt;
                }
            }
            o[0] += acc[0];
            o[1] += acc[1];
            o[2] += acc[2];
        }
    }
}

typedef struct {
    const MvFrame *F;
    StcMVis *v;
} MvFin;

static void mv_final_rows(void *ctx, int y0, int y1)
{
    const MvFin *G = (const MvFin *)ctx;
    StcMVis *v = G->v;
    const MvFrame *F = G->F;
    double kx = (double)v->blw[0] / v->w;
    double ky = (double)v->blh[0] / v->h;
    int y;

    for (y = y0; y < y1; y++) {
        int x;

        for (x = 0; x < v->w; x++) {
            size_t o = (size_t)y * (size_t)v->w + (size_t)x;
            const float *p = v->hdr + o * 3u;
            float bl[3];
            double nx = (x + 0.5) / v->w * 2.0 - 1.0;
            double ny = (y + 0.5) / v->h * 2.0 - 1.0;
            double vig = 1.0 - 0.12 * (nx * nx * 0.8 + ny * ny);
            double e = 1.15 * F->fade * vig;
            uint32_t hs = mv_hash((uint32_t)o * 2654435761U ^
                                  (uint32_t)v->frames * 97U);
            int dith = (int)(hs & 3u) - 1;
            int c[3];
            int k;

            mv_bilin(v->blm[0], v->blw[0], v->blh[0], (x + 0.5) * kx,
                     (y + 0.5) * ky, bl);
            for (k = 0; k < 3; k++) {
                float val = (float)(e * (p[k] + 0.04f * bl[k]));

                c[k] = mv_tm(val) + (dith > 0 && k == (int)((hs >> 4) % 3u)
                                         ? 1 : 0);
                if (c[k] > 255) {
                    c[k] = 255;
                }
            }
            v->fb.px[o].r = (uint8_t)c[0];
            v->fb.px[o].g = (uint8_t)c[1];
            v->fb.px[o].b = (uint8_t)c[2];
        }
    }
}

/* --- HUD (8-bit) -------------------------------------------------------------- */

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
    p->r = (uint8_t)(p->r + (c.r - p->r) * a);
    p->g = (uint8_t)(p->g + (c.g - p->g) * a);
    p->b = (uint8_t)(p->b + (c.b - p->b) * a);
}

static void hud_dim(StcRgbBuf *fb, int x0, int y0, int w, int h, double k)
{
    int x;
    int y;

    for (y = y0; y < y0 + h; y++) {
        for (x = x0; x < x0 + w; x++) {
            StcRgb *p;

            if (x < 0 || y < 0 || x >= fb->w || y >= fb->h) {
                continue;
            }
            p = &fb->px[(size_t)y * (size_t)fb->w + (size_t)x];
            p->r = (uint8_t)(p->r * k);
            p->g = (uint8_t)(p->g * k);
            p->b = (uint8_t)(p->b * k);
        }
    }
}

static int hud_text(StcRgbBuf *fb, int x, int y, const char *s, StcRgb c,
                    int sc, double a)
{
    StcRgb sh = stc_render_rgb(0, 0, 0);
    int cx = x;

    for (; *s != '\0'; s++) {
        const unsigned char *g = stc_vis_glyph(*s);
        int r;
        int q;

        for (r = 0; r < 7; r++) {
            for (q = 0; q < 5; q++) {
                int dx;
                int dy;

                if (!((g[r] >> (4 - q)) & 1u)) {
                    continue;
                }
                for (dy = 0; dy < sc; dy++) {
                    for (dx = 0; dx < sc; dx++) {
                        int px = cx + q * sc + dx;
                        int py = y + r * sc + dy;

                        hud_blend(fb, px + 1, py + 1, sh, 0.55 * a);
                        hud_blend(fb, px, py, c, a);
                    }
                }
            }
        }
        cx += 6 * sc;
    }
    return cx;
}

static int hud_tw(const char *s, int sc)
{
    return (int)strlen(s) * 6 * sc;
}

/* Ops panel: darkened interior, 1 px rule, corner ticks. */
static void hud_panel(StcRgbBuf *fb, int x, int y, int w, int h, double a)
{
    StcRgb rule = stc_render_rgb(120, 124, 132);
    int i;

    hud_dim(fb, x, y, w, h, 1.0 - 0.7 * a);
    for (i = x; i < x + w; i++) {
        hud_blend(fb, i, y, rule, 0.55 * a);
        hud_blend(fb, i, y + h - 1, rule, 0.55 * a);
    }
    for (i = y; i < y + h; i++) {
        hud_blend(fb, x, i, rule, 0.55 * a);
        hud_blend(fb, x + w - 1, i, rule, 0.55 * a);
    }
    for (i = 0; i < 6; i++) {
        hud_blend(fb, x + i, y + 1, rule, a);
        hud_blend(fb, x + 1, y + i, rule, a);
        hud_blend(fb, x + w - 1 - i, y + h - 2, rule, a);
        hud_blend(fb, x + w - 2, y + h - 1 - i, rule, a);
    }
}

/* Strain plot: face-on r h+(t) normalised to the peak, scrolling window
 * of the last `win` M ending now, with a time axis in ms relative to
 * the strain peak, tick marks every 20 ms, and the peak marker. */
static void hud_chirp(StcMVis *v, int x0, int y0, int w, int h, int sc,
                      double fade)
{
    const StcBbh *b = &v->bbh;
    double win = 520.0;
    double tms = v->msun * 4.925491e-6 * 1000.0; /* ms per M */
    int s1 = sc > 1 ? sc - 1 : 1;
    int px0 = x0 + 10 * s1 * 4;
    int px1 = x0 + w - 10 * s1;
    int py0 = y0 + 18 * s1;
    int py1 = y0 + h - 12 * s1;
    int pw = px1 - px0;
    int mid = (py0 + py1) / 2;
    double amp = 0.5 * (py1 - py0) * 0.9;
    StcRgb trace = stc_render_rgb(232, 234, 238);
    StcRgb axis = stc_render_rgb(120, 124, 132);
    StcRgb lab = stc_render_rgb(150, 154, 162);
    StcRgb amber = stc_render_rgb(236, 176, 72);
    char buf[64];
    int prev = -1;
    int x;
    double t_lo = v->t - win;
    double first;
    double tk;

    hud_panel(&v->fb, x0, y0, w, h, fade);
    /* axes */
    for (x = px0; x <= px1; x++) {
        hud_blend(&v->fb, x, mid, axis, 0.35 * fade);
        hud_blend(&v->fb, x, py1, axis, 0.8 * fade);
    }
    for (x = py0; x <= py1; x++) {
        hud_blend(&v->fb, px0, x, axis, 0.8 * fade);
    }
    hud_text(&v->fb, px0 - 6 * s1 * 4 + 2, py0 - 3, "+1", lab, s1, fade);
    hud_text(&v->fb, px0 - 6 * s1 * 4 + 2, mid - 3, " 0", lab, s1, fade);
    hud_text(&v->fb, px0 - 6 * s1 * 4 + 2, py1 - 7 * s1, "-1", lab, s1,
             fade);
    /* time ticks every 20 ms relative to the peak */
    first = ceil(((t_lo - b->t_peak) * tms) / 20.0) * 20.0;
    for (tk = first; tk <= (v->t - b->t_peak) * tms; tk += 20.0) {
        double tt = b->t_peak + tk / tms;
        int xt = px0 + (int)((tt - t_lo) / win * pw);
        int k;
        int major = fmod(fabs(tk), 100.0) < 1e-6;

        for (k = 0; k < (major ? 5 : 3) * s1; k++) {
            hud_blend(&v->fb, xt, py1 + k, axis, 0.9 * fade);
        }
        if (major && xt > px0 + 20 && xt < px1 - 20) {
            snprintf(buf, sizeof buf, "%+.0f", tk);
            hud_text(&v->fb, xt - hud_tw(buf, s1) / 2, py1 + 6 * s1, buf, lab,
                     s1, fade);
        }
    }
    /* trace */
    for (x = 0; x <= pw; x++) {
        double tt = t_lo + win * x / pw;
        double hv = stc_bbh_rh_faceon(b, tt) / b->amp_peak;
        int yy = mid - (int)floor(hv * amp + 0.5);
        int lo;
        int hi;
        int yv;

        if (tt < 0.0) {
            prev = -1;
            continue;
        }
        lo = prev < 0 ? yy : (prev < yy ? prev : yy);
        hi = prev < 0 ? yy : (prev > yy ? prev : yy);
        for (yv = lo; yv <= hi; yv++) {
            hud_blend(&v->fb, px0 + x, yv, trace, 0.9 * fade);
        }
        prev = yy;
    }
    if (b->t_peak <= v->t && b->t_peak >= t_lo) {
        int mx = px0 + (int)((b->t_peak - t_lo) / win * pw);
        int yv;

        for (yv = py0; yv < py1; yv += 3) {
            hud_blend(&v->fb, mx, yv, amber, 0.8 * fade);
        }
        hud_text(&v->fb, mx + 4, py0, "PEAK", amber, s1, fade);
    }
    hud_text(&v->fb, x0 + 8, y0 + 4, "STRAIN H+ / PEAK  (FACE-ON)", lab, s1,
             fade);
    snprintf(buf, sizeof buf, "T - T_PEAK (MS)");
    hud_text(&v->fb, x0 + w - hud_tw(buf, s1) - 8, y0 + 4, buf, lab, s1,
             fade);
}

static void mv_paint_hud(StcMVis *v, const MvFrame *F)
{
    const StcBbh *b = &v->bbh;
    int sc = v->h >= 1100 ? 3 : 2;
    int s1 = sc > 2 ? 2 : 1;
    int lh = 10 * s1 + 3;
    int pad = 10 * s1;
    int x = pad + 8;
    int y;
    StcRgb hi = stc_render_rgb(226, 228, 232);
    StcRgb md = stc_render_rgb(168, 172, 180);
    StcRgb lo = stc_render_rgb(112, 116, 124);
    StcRgb amber = stc_render_rgb(236, 176, 72);
    StcRgb badge = md;
    const char *stage;
    char buf[160];
    double tsec = v->msun * 4.925491e-6;
    double km = v->msun * 1.4766250;
    double a = F->fade;
    int barh = 2 * lh + 10;

    if (!v->hud_on) {
        return;
    }
    if (F->st.stage == STC_BBH_INSPIRAL) {
        stage = "INSPIRAL";
    } else if (F->st.stage == STC_BBH_PLUNGE) {
        stage = "PLUNGE / COMMON HORIZON";
        badge = amber;
    } else if (v->t - b->t_peak > 6.0 * b->tau_qnm) {
        stage = "REMNANT KERR BLACK HOLE";
    } else {
        stage = "RINGDOWN";
        badge = amber;
    }
    /* top bar */
    hud_panel(&v->fb, pad, pad, v->w - 2 * pad, barh, a);
    hud_text(&v->fb, x, pad + 5, "SPACETIME LAB  //  BINARY BLACK HOLE "
             "COALESCENCE", hi, s1 + 1, a);
    snprintf(buf, sizeof buf, "%s", kPre[v->preset].name);
    {
        char *c;

        for (c = buf; *c; c++) {
            if (*c >= 'a' && *c <= 'z') {
                *c = (char)(*c - 32);
            }
        }
    }
    hud_text(&v->fb, x, pad + 5 + lh + 2, buf, md, s1, a);
    hud_text(&v->fb, v->w - pad - 8 - hud_tw(stage, s1 + 1), pad + 5, stage,
             badge, s1 + 1, a);
    snprintf(buf, sizeof buf, "SPEED X%.3g%s   Q %.3g   ETA %.4f",
             v->speed, v->paused ? "  PAUSED" : "", b->m1 / b->m2, b->eta);
    hud_text(&v->fb, v->w - pad - 8 - hud_tw(buf, s1), pad + 5 + lh + 2, buf,
             md, s1, a);

    /* readout panel */
    {
        int pw = 36 * 6 * s1 + 16;
        int ph = 7 * lh + 12;
        int py = pad + barh + 8;

        hud_panel(&v->fb, pad, py, pw, ph, a);
        y = py + 6;
        snprintf(buf, sizeof buf, "M1 / M2      %.1f / %.1f MSUN",
                 b->m1 * v->msun, b->m2 * v->msun);
        hud_text(&v->fb, x, y, buf, md, s1, a);
        y += lh;
        if (F->st.stage != STC_BBH_RINGDOWN) {
            snprintf(buf, sizeof buf, "SEPARATION   %7.2f M  %7.0f KM",
                     F->st.a, F->st.a * km);
        } else {
            snprintf(buf, sizeof buf, "SEPARATION   --");
        }
        hud_text(&v->fb, x, y, buf, hi, s1, a);
        y += lh;
        snprintf(buf, sizeof buf, "F GW         %7.1f HZ",
                 F->st.omega / 3.141592653589793 / tsec);
        hud_text(&v->fb, x, y, buf, hi, s1, a);
        y += lh;
        snprintf(buf, sizeof buf, "T - T_PEAK   %+8.1f MS",
                 (v->t - b->t_peak) * tsec * 1000.0);
        hud_text(&v->fb, x, y, buf, hi, s1, a);
        y += lh;
        snprintf(buf, sizeof buf, "REMNANT      %.1f MSUN  A %.3f",
                 b->mf * v->msun, b->af);
        hud_text(&v->fb, x, y, buf,
                 F->st.stage == STC_BBH_RINGDOWN ? hi : lo, s1, a);
        y += lh;
        snprintf(buf, sizeof buf, "RADIATED     %.2f MSUN C^2",
                 (1.0 - b->mf) * v->msun);
        hud_text(&v->fb, x, y, buf,
                 F->st.stage == STC_BBH_RINGDOWN ? hi : lo, s1, a);
        y += lh;
        snprintf(buf, sizeof buf, "QNM 220      %.1f HZ  TAU %.2f MS",
                 b->om_qnm / (2.0 * 3.141592653589793) / tsec,
                 b->tau_qnm * tsec * 1000.0);
        hud_text(&v->fb, x, y, buf,
                 F->st.stage == STC_BBH_RINGDOWN ? hi : lo, s1, a);
    }
    /* strain plot + footer */
    {
        int fh = 2 * lh + 10;
        int cw = v->w * 50 / 100;
        int ch = v->h * 17 / 100;
        int cx0 = (v->w - cw) / 2;
        int cy0 = v->h - pad - fh - 8 - ch;

        hud_chirp(v, cx0, cy0, cw, ch, sc, a);
        hud_panel(&v->fb, pad, v->h - pad - fh, v->w - 2 * pad, fh, a);
        hud_text(&v->fb, x, v->h - pad - fh + 5,
                 "SPACE PAUSE  R RESTART  1-4 SYSTEM  +/- SPEED  G SHEET  "
                 "L LENS  W WAVES  T TRAILS  C CAM  H HUD  DRAG/SCROLL  Q QUIT",
                 md, s1, a);
        hud_text(&v->fb, x, v->h - pad - fh + 5 + lh,
                 v->status[0] != '\0' ? v->status :
                 "INSPIRAL: 3.5PN TAYLORT4 (ANALYTIC)  MERGER/RINGDOWN/REMNANT: "
                 "ILLUSTRATIVE  SHEET DEPTH + RIPPLES EXAGGERATED",
                 v->status[0] != '\0' ? hi : lo, s1, a);
    }
}

/* --- state -------------------------------------------------------------------- */

static void mv_free_bufs(StcMVis *v)
{
    int i;

    stc_render_rgb_free(&v->fb);
    free(v->far);
    free(v->nearl);
    free(v->hdr);
    v->far = v->nearl = v->hdr = NULL;
    for (i = 0; i < 7; i++) {
        free(v->blm[i]);
        v->blm[i] = NULL;
    }
    v->nbl = 0;
}

int stc_mvis_resize(StcMVis *v, int w, int h)
{
    size_t n;
    int i;
    int bw;
    int bh;

    if (v == NULL || w < 16 || h < 16) {
        return -1;
    }
    mv_free_bufs(v);
    v->w = w;
    v->h = h;
    n = (size_t)w * (size_t)h * 3u;
    if (stc_render_rgb_init(&v->fb, w, h) != 0) {
        return -1;
    }
    v->far = (float *)calloc(n, sizeof(float));
    v->nearl = (float *)calloc(n, sizeof(float));
    v->hdr = (float *)calloc(n, sizeof(float));
    if (v->far == NULL || v->nearl == NULL || v->hdr == NULL) {
        mv_free_bufs(v);
        return -1;
    }
    bw = w;
    bh = h;
    for (i = 0; i < 7; i++) {
        bw = (bw + 1) / 2;
        bh = (bh + 1) / 2;
        if (bw < 4 || bh < 4) {
            break;
        }
        v->blw[i] = bw;
        v->blh[i] = bh;
        v->blm[i] = (float *)calloc((size_t)bw * (size_t)bh * 3u,
                                    sizeof(float));
        if (v->blm[i] == NULL) {
            mv_free_bufs(v);
            return -1;
        }
        v->nbl = i + 1;
    }
    return 0;
}

static void mv_load(StcMVis *v, int preset)
{
    const MvPreset *p = &kPre[preset];

    v->preset = preset;
    v->msun = p->msun;
    stc_bbh_free(&v->bbh); /* preset switch re-inits the same struct */
    (void)stc_bbh_init(&v->bbh, p->q, p->a0, p->a_merge);
    v->t = 0.0;
}

int stc_mvis_init(StcMVis *v, int preset, int w, int h)
{
    if (v == NULL) {
        return -1;
    }
    memset(v, 0, sizeof(*v));
    if (preset < 0 || preset >= STC_MVIS_NPRESET || w < 16 || h < 16) {
        return -1;
    }
    mv_lut_build();
    if (mv_sky_build() != 0) {
        return -2;
    }
    mv_load(v, preset);
    v->speed = 1.0;
    v->grid_on = 1;
    v->lens_on = 1;
    v->hud_on = 1;
    v->trails_on = 1;
    v->autocam = 1;
    v->wave_gain = 1.0;
    v->yaw = 0.55;
    v->pitch = 0.58;
    v->dist = 66.0;
    if (stc_mvis_resize(v, w, h) != 0) {
        return -2;
    }
    return 0;
}

void stc_mvis_free(StcMVis *v)
{
    if (v == NULL) {
        return;
    }
    stc_bbh_free(&v->bbh);
    mv_free_bufs(v);
}

void stc_mvis_step(StcMVis *v)
{
    if (v == NULL) {
        return;
    }
    v->frames++;
    if (v->autocam) {
        v->yaw += 0.0011;
    }
    if (v->paused) {
        return;
    }
    v->t += mv_rate(v) * v->speed;
    if (v->t > mv_t_end(v)) {
        v->t = 0.0;
    }
}

void stc_mvis_orbit(StcMVis *v, double dyaw, double dpitch)
{
    if (v == NULL || !isfinite(dyaw) || !isfinite(dpitch)) {
        return;
    }
    v->yaw -= dyaw;
    v->pitch += dpitch;
    if (v->pitch < 0.05) {
        v->pitch = 0.05;
    }
    if (v->pitch > 1.45) {
        v->pitch = 1.45;
    }
}

void stc_mvis_zoom(StcMVis *v, double factor)
{
    if (v == NULL || !(factor > 0.0) || !isfinite(factor)) {
        return;
    }
    v->dist *= factor;
    if (v->dist < 25.0) {
        v->dist = 25.0;
    }
    if (v->dist > 320.0) {
        v->dist = 320.0;
    }
}

void stc_mvis_render(StcMVis *v)
{
    MvFrame F;
    MvFin G;
    size_t n;
    int i;

    if (v == NULL || v->fb.px == NULL || v->hdr == NULL) {
        return;
    }
    mv_frame_setup(v, &F);
    n = (size_t)v->w * (size_t)v->h * 3u;
    memset(v->far, 0, n * sizeof(float));
    memset(v->nearl, 0, n * sizeof(float));
    if (v->grid_on) {
        mv_paint_grid(&F);
    }
    if (v->trails_on) {
        mv_paint_trails(&F);
    }
    mv_par(mv_lens_rows, &F, v->h);

    /* bloom: threshold-downsample, chain down, tent-upsample back up */
    if (v->nbl > 0) {
        MvDown D;
        MvUp U;

        D.src = v->hdr;
        D.sw = v->w;
        D.sh = v->h;
        D.dst = v->blm[0];
        D.dw = v->blw[0];
        D.dh = v->blh[0];
        D.thresh = 1;
        mv_par(mv_down_rows, &D, D.dh);
        for (i = 1; i < v->nbl; i++) {
            D.src = v->blm[i - 1];
            D.sw = v->blw[i - 1];
            D.sh = v->blh[i - 1];
            D.dst = v->blm[i];
            D.dw = v->blw[i];
            D.dh = v->blh[i];
            D.thresh = 0;
            mv_par(mv_down_rows, &D, D.dh);
        }
        for (i = v->nbl - 1; i >= 1; i--) {
            U.src = v->blm[i];
            U.sw = v->blw[i];
            U.sh = v->blh[i];
            U.dst = v->blm[i - 1];
            U.dw = v->blw[i - 1];
            U.dh = v->blh[i - 1];
            mv_par(mv_up_rows, &U, U.dh);
        }
    }
    G.F = &F;
    G.v = v;
    mv_par(mv_final_rows, &G, v->h);
    mv_paint_hud(v, &F);
}

int stc_mvis_key(StcMVis *v, int key)
{
    if (v == NULL) {
        return 0;
    }
    v->status[0] = '\0';
    switch (key) {
    case 'q':
    case 'Q':
    case STC_UI_KEY_ESC:
    case STC_UI_KEY_CTRL_C:
        return 1;
    case ' ':
    case 'p':
    case 'P':
        v->paused = !v->paused;
        break;
    case 'r':
    case 'R':
        v->t = 0.0;
        break;
    case '1':
    case '2':
    case '3':
    case '4':
        mv_load(v, key - '1');
        break;
    case '+':
    case '=':
        v->speed = v->speed * 2.0 > 16.0 ? 16.0 : v->speed * 2.0;
        break;
    case '-':
    case '_':
        v->speed = v->speed * 0.5 < 0.0625 ? 0.0625 : v->speed * 0.5;
        break;
    case 'g':
    case 'G':
        v->grid_on = !v->grid_on;
        break;
    case 'l':
    case 'L':
        v->lens_on = !v->lens_on;
        snprintf(v->status, sizeof v->status, "lensing %s",
                 v->lens_on ? "on" : "off");
        break;
    case 'w':
    case 'W':
        v->wave_gain = v->wave_gain >= 4.0 ? 0.5 : v->wave_gain * 2.0;
        snprintf(v->status, sizeof v->status, "ripple gain x%.1f",
                 v->wave_gain);
        break;
    case 't':
    case 'T':
        v->trails_on = !v->trails_on;
        break;
    case 'h':
    case 'H':
        v->hud_on = !v->hud_on;
        break;
    case 'c':
    case 'C':
        v->autocam = !v->autocam;
        break;
    case 'a':
    case 'A':
    case STC_UI_KEY_LEFT:
        stc_mvis_orbit(v, 0.05, 0.0);
        break;
    case 'd':
    case 'D':
    case STC_UI_KEY_RIGHT:
        stc_mvis_orbit(v, -0.05, 0.0);
        break;
    case 's':
    case 'S':
    case STC_UI_KEY_DOWN:
        stc_mvis_orbit(v, 0.0, -0.04);
        break;
    case STC_UI_KEY_UP:
        stc_mvis_orbit(v, 0.0, 0.04);
        break;
    default:
        break;
    }
    return 0;
}

int stc_mvis_shot(int preset, long nframes, int w, int h,
                  const char *ppm_path, const char *dir, long dir_every)
{
    StcMVis v;
    long f;
    int rc;

    if (nframes < 1 || (ppm_path == NULL && (dir == NULL ||
                                               dir_every <= 0))) {
        return -1;
    }
    rc = stc_mvis_init(&v, preset, w, h);
    if (rc != 0) {
        return rc == -1 ? -1 : -2;
    }
    for (f = 0; f < nframes; f++) {
        int last = f == nframes - 1;
        int dump = dir != NULL && dir_every > 0 && f % dir_every == 0;

        if (last || dump) {
            stc_mvis_render(&v);
        }
        if (dump) {
            char path[1024];

            snprintf(path, sizeof path, "%s/merger_%05ld.ppm", dir,
                     f / dir_every);
            if (stc_render_write_ppm(&v.fb, path) != 0) {
                stc_mvis_free(&v);
                return -3;
            }
        }
        if (last && ppm_path != NULL &&
            stc_render_write_ppm(&v.fb, ppm_path) != 0) {
            stc_mvis_free(&v);
            return -3;
        }
        stc_mvis_step(&v);
    }
    stc_mvis_free(&v);
    return 0;
}
