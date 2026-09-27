/* SPACETIME LAB — space: large-N gravitational dynamics core.
 * Implements include/stc/space.h exactly (SI units, SoA, stable
 * indices, no compaction).
 *
 * Integrator: 6th-order Yoshida "solution A" (Yoshida 1990,
 * Phys. Lett. A 150, 262): a symmetric 7-stage product of KDK
 * leapfrogs S2(w*h) with
 *   w = {w3, w2, w1, w0, w1, w2, w3},
 *   w1 = -1.17767998417887, w2 = 0.235573213359357,
 *   w3 =  0.784513610477560, w0 = 1 - 2(w1+w2+w3).
 * Adjacent half-kicks merge, giving 8 kicks + 7 drifts per substep.
 * The acceleration from the last kick of a substep is cached and
 * reused as the first kick of the next while the force is
 * position-only (pn == 0, or no live 1PN source), so a steady
 * Newtonian run costs 7 force evaluations per substep. With pn != 0
 * the acceleration is velocity-dependent and a kick that froze v at
 * its start value would drop the composition's effective order
 * (measured defect: Mercury's 1PN error at h = 6 h was ~50x the
 * Newton-only error). Every PN kick therefore solves the implicit
 * midpoint v_mid = v + (tau/2) a(x, v_mid), then applies
 * v += tau a(x, v_mid) = 2 v_mid - v. Only the 1PN term depends on
 * velocity, so the kick splits a(x,v) = a_static(x) + a_pn(x,v):
 * a_static (Newton + zonal) is evaluated once per kick and cached
 * across substep boundaries (avalid == 2), while a_pn involves only
 * the few live PN sources and reads source velocities from the
 * frozen state — so each particle's fixed point is independent and
 * iterated locally in registers (bounded by SP_MID_ITERS, exits when
 * the relative velocity change is < 1e-15 — a deterministic stopping
 * rule). Exact rewind is still provided by checkpoints.
 *
 * Forces: direct O(nm^2) per-target sums while the live massive count
 * is <= tree_threshold, else a monopole Barnes-Hut octree rebuilt at
 * every force evaluation. Tree cells are accepted when
 *   size < theta * (d - rho),   rho = sqrt(3)*half + |com - centre|,
 * i.e. the opening angle is compared against the distance to the
 * NEAREST point mass could occupy — strictly more conservative than
 * the classic size/d < theta (measured force errors are dominated by
 * the strict form; see stspace "tree" check). The cell containing the
 * target is always opened so self-mass never contributes.
 *
 * Optional 1PN (header formula, test-particle Schwarzschild,
 * harmonic gauge): every live STAR/BLACKHOLE sources
 *   a += GM/(c^2 r^3) [ (4GM/r - v^2) r_vec + 4 (r_vec . v_vec) v_vec ]
 * on every other live particle, r_vec/v_vec measured relative to the
 * source. Reaction on the source is neglected (documented in header).
 *
 * Optional zonal harmonics: every live massive body with j2 or j4
 * nonzero (stc_space_set_oblateness) adds the axisymmetric multipole
 * acceleration about its unit spin pole (derivation at sp_zonal) to
 * every other particle — massive and tracer, in both the direct and
 * tree paths (the tree only ever carries the monopole). The
 * equal-and-opposite reaction is applied to the oblate body for
 * massive targets (m_target a_target = -m_planet a_planet), so total
 * momentum is conserved.
 *
 * Events (collisions != 0, applied after every substep piece):
 *   massive pair: BH capture first (dist < 2 r_s), then tidal
 *   disruption (r_t = R_b (M_bh/m_b)^(1/3), only when r_t exceeds the
 *   capture radius), then contact merge (dist < R_i + R_j). Two black
 *   holes merge inside the larger capture radius; the heavier one is
 *   the survivor. Merges/captures conserve mass, momentum and the
 *   centre of mass; radius = summed-volume sphere, or r_s = 2GM/c^2
 *   when the survivor is (or becomes) a black hole.
 *   tracers: inside a non-BH radius -> ABSORB; inside 2 r_s of a BH
 *   -> CAPTURE. When several bodies claim one tracer the lowest
 *   massive index wins (fixed rule -> deterministic).
 *   Detection: massive pairs via an x-sorted sweep + per-pair reach
 *   window (reach_i = R_i, or max(2 r_s, R_max*(gm_i/gm_min)^(1/3))
 *   for black holes), candidates sorted by (i,j) then applied with a
 *   fresh distance test. Tracers are tested in parallel against the
 *   x-sorted radius-bearing bodies (binary-search window + lowest
 *   index win), applied sequentially in tracer index order.
 *   DISRUPT spawns min(disrupt_tracers, free tracer slots) tracers
 *   uniformly in the victim's volume carrying the victim's velocity;
 *   the stream is seeded by (event index, body index) — no rand().
 *
 * Checkpoints: full-state snapshots (massive + tracers + t + steps +
 * event log + e0) on an absolute time grid k*ckpt_interval (both
 * directions). A substep that would straddle a grid time is SPLIT at
 * the grid line, so every crossed grid time gets an exact checkpoint
 * at a real substep boundary; an initial snapshot is stored at the
 * first advance's start time, and one at every advance() end. The
 * ring keeps the ckpt_max snapshots nearest the current time.
 *
 *   SEEK RULE: seek(t) picks the NEWEST stored checkpoint with
 *   t_ck <= t and re-integrates forward to t. Re-integration replays
 *   the deterministic dynamics — including events — so a forward
 *   seek reproduces the original timeline bitwise when t is a stored
 *   checkpoint time (trivially) or a state on the same substep grid.
 *   Checkpoints are never taken ahead of t and used backward (a
 *   merge cannot be un-done going backward), which is what "the
 *   correct side" means here. If every stored checkpoint is later
 *   than t the seek returns -1 (outside the checkpointed span) and
 *   the state is untouched.
 *
 * Determinism (R4): every parallel loop writes only its own index
 * slots and sums in a fixed index order; reductions are ordered.
 * STC_PAR_THREADS=1 gives bitwise-identical state to any other count.
 */

#include "stc/space.h"

#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include "stc/par.h"
#include "stc/rng.h"
#include "stc/units.h"

/* ---------------------------------------------------------------- */
/* Yoshida solution-A coefficients                                   */

#define SP_W1 (-1.17767998417887)
#define SP_W2 0.235573213359357
#define SP_W3 0.784513610477560
#define SP_W0 (1.0 - 2.0 * (SP_W1 + SP_W2 + SP_W3))

/* drift coefficients, in application order */
static const double sp_yw[7] = {
    SP_W3, SP_W2, SP_W1, SP_W0, SP_W1, SP_W2, SP_W3
};
/* merged kick coefficients (first/last are half kicks) */
static const double sp_yc[8] = {
    SP_W3 * 0.5,
    (SP_W3 + SP_W2) * 0.5,
    (SP_W2 + SP_W1) * 0.5,
    (SP_W1 + SP_W0) * 0.5,
    (SP_W0 + SP_W1) * 0.5,
    (SP_W1 + SP_W2) * 0.5,
    (SP_W2 + SP_W3) * 0.5,
    SP_W3 * 0.5
};

#define SP_MAX_DEPTH 40
#define SP_SQRT3 1.7320508075688772
#define SP_MASS_G (1.0 / STC_G) /* kg per unit of GM */

/* ---------------------------------------------------------------- */
/* Octree node + snapshot + private scratch                          */

typedef struct {
    double cx, cy, cz;  /* cell centre */
    double half;        /* half edge length */
    double gm;          /* summed GM in the cell */
    double mx, my, mz;  /* centre of mass */
    int child[8];       /* -1: leaf or empty */
    int head;           /* leaf: head of body chain in tnext[], -1 */
} SpNode;

typedef struct {
    double t;
    long steps;
    int nm, nt, nev;
    double e0;
    StcSpaceEvent ev[STC_SPACE_EVLOG];
    int mcap, tcap;
    double *mv;          /* 14*mcap: x y z vx vy vz gm radius
                            j2 j4 rref polx poly polz */
    int *mi;             /* 2*mcap: kind parent */
    unsigned char *ma;   /* mcap alive */
    char *mn;            /* STC_SPACE_NAME*mcap */
    double *tv;          /* 6*tcap: tx ty tz tvx tvy tvz */
    unsigned char *ta;   /* tcap */
} SpSnap;

struct StcSpaceCkpt {
    /* snapshot ring, sorted by t ascending */
    SpSnap *sn;
    int count;
    int sncap;
    /* integrator scratch (sized to capacities, grown lazily) */
    int acap_m, acap_t;
    int avalid;              /* accel cache: 0 stale, 1 full a(x,v),
                                2 static a_static(x) (Newton+zonal) */
    double *ax, *ay, *az;    /* massive accelerations */
    double *tax, *tay, *taz; /* tracer accelerations */
    /* tree arena */
    SpNode *tn;
    int tncap;
    int *tnext;              /* per-body leaf chain, cap_m */
    int tnxc;
    /* event scratch */
    int *ilist;              /* live indices, cap_m */
    double *reach;           /* per-body event reach, cap_m */
    int *src;                /* 1PN source indices, cap_m */
    int *obl;                /* oblate (j2/j4) body indices, cap_m */
    int *thit;               /* tracer -> massive hit, cap_t */
    uint64_t *pairs;         /* massive candidate pairs */
    long paircap;
    double *part;            /* per-index reduction partials, cap_m */
    double *part3;           /* momentum partials, 3*cap_m */
    double *tpart;           /* tracer reduction partials, cap_t */
    double *mvm;             /* midpoint kick: staged source v, 3*cap_m */
    int ilcap, srcflag, oblflag;
};

typedef struct {
    double key;
    int idx;
} SpEnt;

/* ---------------------------------------------------------------- */
/* small utilities                                                   */

static int sp_ent_cmp(const void *pa, const void *pb)
{
    const SpEnt *a = (const SpEnt *)pa;
    const SpEnt *b = (const SpEnt *)pb;
    if (a->key < b->key) {
        return -1;
    }
    if (a->key > b->key) {
        return 1;
    }
    return a->idx - b->idx;
}

static int sp_u64_cmp(const void *pa, const void *pb)
{
    uint64_t a = *(const uint64_t *)pa;
    uint64_t b = *(const uint64_t *)pb;
    return (a > b) - (a < b);
}

static void *sp_xrealloc(void *p, size_t n)
{
    return realloc(p, n ? n : 1);
}

static int sp_is_bh_kind(int k)
{
    return k == STC_SK_BLACKHOLE;
}

static int sp_is_pn_source(int k)
{
    return k == STC_SK_STAR || k == STC_SK_BLACKHOLE;
}

/* ---------------------------------------------------------------- */
/* capacity / scratch                                                */

/* Grow the private scratch to the current capacities. Returns 0 ok,
 * -2 on allocation failure (partial growth is harmless). */
static int sp_fit_scratch(const StcSpace *s)
{
    StcSpaceCkpt *ck = s->ck;
    int m = s->cap_m;
    int t = s->cap_t;
    if (ck->acap_m != m) {
        double *ax = (double *)sp_xrealloc(ck->ax, sizeof(double) * (size_t)m);
        double *ay = (double *)sp_xrealloc(ck->ay, sizeof(double) * (size_t)m);
        double *az = (double *)sp_xrealloc(ck->az, sizeof(double) * (size_t)m);
        int *tnext = (int *)sp_xrealloc(ck->tnext, sizeof(int) * (size_t)m);
        int *ilist = (int *)sp_xrealloc(ck->ilist, sizeof(int) * (size_t)m);
        int *src = (int *)sp_xrealloc(ck->src, sizeof(int) * (size_t)m);
        int *obl = (int *)sp_xrealloc(ck->obl, sizeof(int) * (size_t)m);
        double *reach = (double *)sp_xrealloc(ck->reach, sizeof(double) * (size_t)m);
        double *part = (double *)sp_xrealloc(ck->part, sizeof(double) * (size_t)m);
        double *part3 = (double *)sp_xrealloc(ck->part3, sizeof(double) * (size_t)m * 3u);
        double *mvm = (double *)sp_xrealloc(ck->mvm, sizeof(double) * (size_t)m * 3u);
        if (!ax || !ay || !az || !tnext || !ilist || !src || !obl ||
            !reach || !part || !part3 || !mvm) {
            free(ax); free(ay); free(az); free(tnext); free(ilist);
            free(src); free(obl); free(reach); free(part); free(part3);
            free(mvm);
            return -2;
        }
        ck->ax = ax; ck->ay = ay; ck->az = az;
        ck->tnext = tnext; ck->ilist = ilist; ck->src = src; ck->obl = obl;
        ck->reach = reach; ck->part = part; ck->part3 = part3;
        ck->mvm = mvm;
        ck->acap_m = m;
    }
    if (ck->acap_t != t) {
        double *tax = (double *)sp_xrealloc(ck->tax, sizeof(double) * (size_t)(t ? t : 1));
        double *tay = (double *)sp_xrealloc(ck->tay, sizeof(double) * (size_t)(t ? t : 1));
        double *taz = (double *)sp_xrealloc(ck->taz, sizeof(double) * (size_t)(t ? t : 1));
        int *thit = (int *)sp_xrealloc(ck->thit, sizeof(int) * (size_t)(t ? t : 1));
        double *tpart = (double *)sp_xrealloc(ck->tpart, sizeof(double) * (size_t)(t ? t : 1));
        if (!tax || !tay || !taz || !thit || !tpart) {
            free(tax); free(tay); free(taz); free(thit); free(tpart);
            return -2;
        }
        ck->tax = tax; ck->tay = tay; ck->taz = taz; ck->thit = thit;
        ck->tpart = tpart;
        ck->acap_t = t;
    }
    return 0;
}

/* ---------------------------------------------------------------- */
/* lifecycle                                                         */

int stc_space_init(StcSpace *s, int cap_m, int cap_t)
{
    int ok;
    if (s == NULL || cap_m < 1 || cap_t < 0) {
        return -1;
    }
    memset(s, 0, sizeof(*s));
    s->cap_m = cap_m;
    s->cap_t = cap_t;
    s->x = (double *)calloc((size_t)cap_m, sizeof(double));
    s->y = (double *)calloc((size_t)cap_m, sizeof(double));
    s->z = (double *)calloc((size_t)cap_m, sizeof(double));
    s->vx = (double *)calloc((size_t)cap_m, sizeof(double));
    s->vy = (double *)calloc((size_t)cap_m, sizeof(double));
    s->vz = (double *)calloc((size_t)cap_m, sizeof(double));
    s->gm = (double *)calloc((size_t)cap_m, sizeof(double));
    s->radius = (double *)calloc((size_t)cap_m, sizeof(double));
    s->kind = (int *)calloc((size_t)cap_m, sizeof(int));
    s->parent = (int *)calloc((size_t)cap_m, sizeof(int));
    s->alive = (unsigned char *)calloc((size_t)cap_m, 1);
    s->name = (char (*)[STC_SPACE_NAME])calloc((size_t)cap_m, STC_SPACE_NAME);
    s->j2 = (double *)calloc((size_t)cap_m, sizeof(double));
    s->j4 = (double *)calloc((size_t)cap_m, sizeof(double));
    s->rref = (double *)calloc((size_t)cap_m, sizeof(double));
    s->polx = (double *)calloc((size_t)cap_m, sizeof(double));
    s->poly = (double *)calloc((size_t)cap_m, sizeof(double));
    s->polz = (double *)calloc((size_t)cap_m, sizeof(double));
    s->tx = (double *)calloc((size_t)(cap_t ? cap_t : 1), sizeof(double));
    s->ty = (double *)calloc((size_t)(cap_t ? cap_t : 1), sizeof(double));
    s->tz = (double *)calloc((size_t)(cap_t ? cap_t : 1), sizeof(double));
    s->tvx = (double *)calloc((size_t)(cap_t ? cap_t : 1), sizeof(double));
    s->tvy = (double *)calloc((size_t)(cap_t ? cap_t : 1), sizeof(double));
    s->tvz = (double *)calloc((size_t)(cap_t ? cap_t : 1), sizeof(double));
    s->talive = (unsigned char *)calloc((size_t)(cap_t ? cap_t : 1), 1);
    s->ck = (StcSpaceCkpt *)calloc(1, sizeof(StcSpaceCkpt));
    ok = s->x && s->y && s->z && s->vx && s->vy && s->vz && s->gm &&
         s->radius && s->kind && s->parent && s->alive && s->name &&
         s->j2 && s->j4 && s->rref && s->polx && s->poly && s->polz &&
         s->tx && s->ty && s->tz && s->tvx && s->tvy && s->tvz &&
         s->talive && s->ck;
    if (!ok) {
        stc_space_free(s);
        return -2;
    }
    s->h = 3600.0;
    s->soft = 0.0;
    s->pn = 1;
    s->theta = 0.5;
    s->tree_threshold = 1024;
    s->disrupt_tracers = 400;
    s->collisions = 1;
    s->ckpt_max = 64;
    s->ckpt_interval = 30.0 * STC_DAY;
    s->t = 0.0;
    return 0;
}

static void sp_snap_free(SpSnap *sn)
{
    free(sn->mv);
    free(sn->mi);
    free(sn->ma);
    free(sn->mn);
    free(sn->tv);
    free(sn->ta);
    memset(sn, 0, sizeof(*sn));
}

void stc_space_free(StcSpace *s)
{
    int k;
    if (s == NULL) {
        return;
    }
    free(s->x); free(s->y); free(s->z);
    free(s->vx); free(s->vy); free(s->vz);
    free(s->gm); free(s->radius); free(s->kind); free(s->parent);
    free(s->alive); free(s->name);
    free(s->j2); free(s->j4); free(s->rref);
    free(s->polx); free(s->poly); free(s->polz);
    free(s->tx); free(s->ty); free(s->tz);
    free(s->tvx); free(s->tvy); free(s->tvz); free(s->talive);
    if (s->ck) {
        StcSpaceCkpt *ck = s->ck;
        for (k = 0; k < ck->sncap; k++) {
            sp_snap_free(&ck->sn[k]);
        }
        free(ck->sn);
        free(ck->ax); free(ck->ay); free(ck->az);
        free(ck->tax); free(ck->tay); free(ck->taz);
        free(ck->tn); free(ck->tnext); free(ck->ilist); free(ck->reach);
        free(ck->src); free(ck->obl); free(ck->thit); free(ck->pairs);
        free(ck->part); free(ck->part3); free(ck->tpart);
        free(ck->mvm);
        free(ck);
    }
    memset(s, 0, sizeof(*s));
}

void stc_space_clear(StcSpace *s, double t)
{
    if (s == NULL) {
        return;
    }
    s->nm = 0;
    s->nt = 0;
    s->t = t;
    s->steps = 0;
    s->e0 = 0.0;
    s->nev = 0;
    memset(s->ev, 0, sizeof(s->ev));
    if (s->ck) {
        s->ck->count = 0;
        s->ck->avalid = 0;
    }
}

int stc_space_reserve(StcSpace *s, int cap_m, int cap_t)
{
    double *x, *y, *z, *vx, *vy, *vz, *gm, *radius;
    double *j2, *j4, *rref, *polx, *poly, *polz;
    double *tx, *ty, *tz, *tvx, *tvy, *tvz;
    int *kind, *parent;
    unsigned char *alive, *talive;
    char (*name)[STC_SPACE_NAME];

    if (s == NULL || cap_m < 1 || cap_t < 0) {
        return -1;
    }
    if (cap_m <= s->cap_m && cap_t <= s->cap_t) {
        return 0;
    }
    if (cap_m > s->cap_m) {
        x = (double *)sp_xrealloc(s->x, sizeof(double) * (size_t)cap_m);
        y = (double *)sp_xrealloc(s->y, sizeof(double) * (size_t)cap_m);
        z = (double *)sp_xrealloc(s->z, sizeof(double) * (size_t)cap_m);
        vx = (double *)sp_xrealloc(s->vx, sizeof(double) * (size_t)cap_m);
        vy = (double *)sp_xrealloc(s->vy, sizeof(double) * (size_t)cap_m);
        vz = (double *)sp_xrealloc(s->vz, sizeof(double) * (size_t)cap_m);
        gm = (double *)sp_xrealloc(s->gm, sizeof(double) * (size_t)cap_m);
        radius = (double *)sp_xrealloc(s->radius, sizeof(double) * (size_t)cap_m);
        kind = (int *)sp_xrealloc(s->kind, sizeof(int) * (size_t)cap_m);
        parent = (int *)sp_xrealloc(s->parent, sizeof(int) * (size_t)cap_m);
        alive = (unsigned char *)sp_xrealloc(s->alive, (size_t)cap_m);
        name = (char (*)[STC_SPACE_NAME])sp_xrealloc(
            s->name, sizeof(*name) * (size_t)cap_m);
        j2 = (double *)sp_xrealloc(s->j2, sizeof(double) * (size_t)cap_m);
        j4 = (double *)sp_xrealloc(s->j4, sizeof(double) * (size_t)cap_m);
        rref = (double *)sp_xrealloc(s->rref, sizeof(double) * (size_t)cap_m);
        polx = (double *)sp_xrealloc(s->polx, sizeof(double) * (size_t)cap_m);
        poly = (double *)sp_xrealloc(s->poly, sizeof(double) * (size_t)cap_m);
        polz = (double *)sp_xrealloc(s->polz, sizeof(double) * (size_t)cap_m);
        if (!x || !y || !z || !vx || !vy || !vz || !gm || !radius ||
            !kind || !parent || !alive || !name ||
            !j2 || !j4 || !rref || !polx || !poly || !polz) {
            return -2;
        }
        s->x = x; s->y = y; s->z = z;
        s->vx = vx; s->vy = vy; s->vz = vz;
        s->gm = gm; s->radius = radius; s->kind = kind;
        s->parent = parent; s->alive = alive; s->name = name;
        s->j2 = j2; s->j4 = j4; s->rref = rref;
        s->polx = polx; s->poly = poly; s->polz = polz;
        memset(s->x + s->cap_m, 0, sizeof(double) * (size_t)(cap_m - s->cap_m));
        memset(s->kind + s->cap_m, 0, sizeof(int) * (size_t)(cap_m - s->cap_m));
        memset(s->alive + s->cap_m, 0, (size_t)(cap_m - s->cap_m));
        memset(s->name + s->cap_m, 0,
               sizeof(*s->name) * (size_t)(cap_m - s->cap_m));
        memset(s->j2 + s->cap_m, 0,
               sizeof(double) * (size_t)(cap_m - s->cap_m));
        memset(s->j4 + s->cap_m, 0,
               sizeof(double) * (size_t)(cap_m - s->cap_m));
        memset(s->rref + s->cap_m, 0,
               sizeof(double) * (size_t)(cap_m - s->cap_m));
        memset(s->polx + s->cap_m, 0,
               sizeof(double) * (size_t)(cap_m - s->cap_m));
        memset(s->poly + s->cap_m, 0,
               sizeof(double) * (size_t)(cap_m - s->cap_m));
        memset(s->polz + s->cap_m, 0,
               sizeof(double) * (size_t)(cap_m - s->cap_m));
        s->cap_m = cap_m;
    }
    if (cap_t > s->cap_t) {
        tx = (double *)sp_xrealloc(s->tx, sizeof(double) * (size_t)cap_t);
        ty = (double *)sp_xrealloc(s->ty, sizeof(double) * (size_t)cap_t);
        tz = (double *)sp_xrealloc(s->tz, sizeof(double) * (size_t)cap_t);
        tvx = (double *)sp_xrealloc(s->tvx, sizeof(double) * (size_t)cap_t);
        tvy = (double *)sp_xrealloc(s->tvy, sizeof(double) * (size_t)cap_t);
        tvz = (double *)sp_xrealloc(s->tvz, sizeof(double) * (size_t)cap_t);
        talive = (unsigned char *)sp_xrealloc(s->talive, (size_t)cap_t);
        if (!tx || !ty || !tz || !tvx || !tvy || !tvz || !talive) {
            return -2;
        }
        s->tx = tx; s->ty = ty; s->tz = tz;
        s->tvx = tvx; s->tvy = tvy; s->tvz = tvz; s->talive = talive;
        memset(s->talive + s->cap_t, 0, (size_t)(cap_t - s->cap_t));
        s->cap_t = cap_t;
    }
    if (s->ck) {
        s->ck->avalid = 0; /* scratch sizes changed */
    }
    return 0;
}

/* ---------------------------------------------------------------- */
/* add / mutate                                                      */

int stc_space_add(StcSpace *s, const char *name, int kind, int parent,
                  double gm, double radius, const double pos[3],
                  const double vel[3])
{
    int i;
    if (s == NULL || pos == NULL || vel == NULL || !(gm > 0.0) ||
        !(radius >= 0.0) || s->nm >= s->cap_m) {
        return -1;
    }
    i = s->nm++;
    s->x[i] = pos[0];
    s->y[i] = pos[1];
    s->z[i] = pos[2];
    s->vx[i] = vel[0];
    s->vy[i] = vel[1];
    s->vz[i] = vel[2];
    s->gm[i] = gm;
    s->radius[i] = radius;
    s->kind[i] = kind;
    s->parent[i] = parent;
    s->alive[i] = 1;
    s->j2[i] = 0.0;
    s->j4[i] = 0.0;
    s->rref[i] = 0.0;
    s->polx[i] = 0.0;
    s->poly[i] = 0.0;
    s->polz[i] = 0.0;
    memset(s->name[i], 0, STC_SPACE_NAME);
    if (name != NULL) {
        strncpy(s->name[i], name, STC_SPACE_NAME - 1);
    }
    if (s->ck) {
        s->ck->avalid = 0;
    }
    return i;
}

int stc_space_add_tracer(StcSpace *s, const double pos[3],
                         const double vel[3])
{
    int i;
    if (s == NULL || pos == NULL || vel == NULL || s->nt >= s->cap_t) {
        return -1;
    }
    i = s->nt++;
    s->tx[i] = pos[0];
    s->ty[i] = pos[1];
    s->tz[i] = pos[2];
    s->tvx[i] = vel[0];
    s->tvy[i] = vel[1];
    s->tvz[i] = vel[2];
    s->talive[i] = 1;
    if (s->ck) {
        s->ck->avalid = 0;
    }
    return i;
}

int stc_space_make_blackhole(StcSpace *s, int i, double gm_new,
                             const char *new_name)
{
    if (s == NULL || i < 0 || i >= s->nm || !s->alive[i]) {
        return -1;
    }
    if (gm_new > 0.0 && gm_new != s->gm[i]) {
        double f = s->gm[i] / gm_new; /* keep m*v: v *= m_old/m_new */
        s->vx[i] *= f;
        s->vy[i] *= f;
        s->vz[i] *= f;
        s->gm[i] = gm_new;
    }
    s->radius[i] = 2.0 * s->gm[i] / (STC_C * STC_C);
    s->kind[i] = STC_SK_BLACKHOLE;
    if (new_name != NULL) {
        memset(s->name[i], 0, STC_SPACE_NAME);
        strncpy(s->name[i], new_name, STC_SPACE_NAME - 1);
    }
    if (s->ck) {
        s->ck->avalid = 0;
    }
    return 0;
}

int stc_space_set_oblateness(StcSpace *s, int i, double j2, double j4,
                             double rref, const double pole[3])
{
    double n2;
    if (s == NULL || i < 0 || i >= s->nm || !s->alive[i] ||
        pole == NULL || !isfinite(j2) || !isfinite(j4) ||
        !isfinite(rref) || !(rref >= 0.0)) {
        return -1;
    }
    n2 = pole[0] * pole[0] + pole[1] * pole[1] + pole[2] * pole[2];
    if (!isfinite(n2)) {
        return -1;
    }
    if ((j2 != 0.0 || j4 != 0.0) && (!(rref > 0.0) || !(n2 > 0.0))) {
        return -1; /* an oblate body needs a reference radius and pole */
    }
    s->j2[i] = j2;
    s->j4[i] = j4;
    s->rref[i] = rref;
    if (n2 > 0.0) {
        double n = sqrt(n2);
        s->polx[i] = pole[0] / n;
        s->poly[i] = pole[1] / n;
        s->polz[i] = pole[2] / n;
    } else {
        s->polx[i] = s->poly[i] = s->polz[i] = 0.0;
    }
    if (s->ck) {
        s->ck->avalid = 0;
    }
    return 0;
}

void stc_space_to_barycentre(StcSpace *s)
{
    double m = 0.0;
    double cx = 0.0, cy = 0.0, cz = 0.0;
    double px = 0.0, py = 0.0, pz = 0.0;
    int i;
    if (s == NULL) {
        return;
    }
    for (i = 0; i < s->nm; i++) {
        if (s->alive[i]) {
            double g = s->gm[i];
            m += g;
            cx += g * s->x[i];
            cy += g * s->y[i];
            cz += g * s->z[i];
            px += g * s->vx[i];
            py += g * s->vy[i];
            pz += g * s->vz[i];
        }
    }
    if (!(m > 0.0)) {
        return;
    }
    cx /= m; cy /= m; cz /= m;
    px /= m; py /= m; pz /= m;
    for (i = 0; i < s->nm; i++) {
        s->x[i] -= cx; s->y[i] -= cy; s->z[i] -= cz;
        s->vx[i] -= px; s->vy[i] -= py; s->vz[i] -= pz;
    }
    for (i = 0; i < s->nt; i++) {
        s->tx[i] -= cx; s->ty[i] -= cy; s->tz[i] -= cz;
        s->tvx[i] -= px; s->tvy[i] -= py; s->tvz[i] -= pz;
    }
    if (s->ck) {
        s->ck->avalid = 0;
    }
}

/* ---------------------------------------------------------------- */
/* octree                                                            */

static void sp_node_init(SpNode *nd, double cx, double cy, double cz,
                         double half)
{
    int k;
    nd->cx = cx; nd->cy = cy; nd->cz = cz; nd->half = half;
    nd->gm = 0.0;
    nd->mx = nd->my = nd->mz = 0.0;
    for (k = 0; k < 8; k++) {
        nd->child[k] = -1;
    }
    nd->head = -1;
}

static int sp_node_leaf(const SpNode *nd)
{
    return nd->child[0] == -1;
}

static int sp_octant(const SpNode *nd, double x, double y, double z)
{
    int o = 0;
    if (x >= nd->cx) {
        o |= 1;
    }
    if (y >= nd->cy) {
        o |= 2;
    }
    if (z >= nd->cz) {
        o |= 4;
    }
    return o;
}

static int sp_subdivide(SpNode *nodes, int *count, int cap, int idx)
{
    double h;
    int k;
    SpNode *nd = &nodes[idx];
    if (*count < 0 || *count + 8 > cap) {
        return 0;
    }
    h = nd->half * 0.5;
    if (!(h > 0.0)) {
        return 0;
    }
    for (k = 0; k < 8; k++) {
        double ox = (k & 1) ? h : -h;
        double oy = (k & 2) ? h : -h;
        double oz = (k & 4) ? h : -h;
        int c = (*count)++;
        sp_node_init(&nodes[c], nd->cx + ox, nd->cy + oy, nd->cz + oz, h);
        nd->child[k] = c;
    }
    return 1;
}

/* Bodies are inserted in ascending index order into a fixed
 * centre/half root cell, so the tree (and therefore every traversal
 * order) is a pure function of the state — deterministic build. */
static void sp_insert(SpNode *nodes, int *count, int cap, int *next,
                      const StcSpace *s, int idx, int b, int depth)
{
    SpNode *nd = &nodes[idx];
    if (sp_node_leaf(nd)) {
        if (nd->head == -1) {
            nd->head = b;
            next[b] = -1;
            return;
        }
        if (depth >= SP_MAX_DEPTH || !sp_subdivide(nodes, count, cap, idx)) {
            next[b] = nd->head;
            nd->head = b;
            return;
        }
        {
            int cur = nd->head;
            nd->head = -1;
            while (cur != -1) {
                int nxt = next[cur];
                int o;
                next[cur] = -1;
                o = sp_octant(nd, s->x[cur], s->y[cur], s->z[cur]);
                sp_insert(nodes, count, cap, next, s, nd->child[o],
                          cur, depth + 1);
                cur = nxt;
            }
        }
        {
            int o = sp_octant(nd, s->x[b], s->y[b], s->z[b]);
            sp_insert(nodes, count, cap, next, s, nd->child[o], b,
                      depth + 1);
        }
        return;
    }
    {
        int o = sp_octant(nd, s->x[b], s->y[b], s->z[b]);
        sp_insert(nodes, count, cap, next, s, nd->child[o], b, depth + 1);
    }
}

static void sp_summarize(SpNode *nodes, const StcSpace *s,
                         const int *next, int idx)
{
    SpNode *nd = &nodes[idx];
    double m = 0.0, x = 0.0, y = 0.0, z = 0.0;
    if (sp_node_leaf(nd)) {
        int cur = nd->head;
        while (cur != -1) {
            double bm = s->gm[cur];
            m += bm;
            x += bm * s->x[cur];
            y += bm * s->y[cur];
            z += bm * s->z[cur];
            cur = next[cur];
        }
    } else {
        int k;
        for (k = 0; k < 8; k++) {
            int c = nd->child[k];
            if (c != -1) {
                sp_summarize(nodes, s, next, c);
                m += nodes[c].gm;
                x += nodes[c].gm * nodes[c].mx;
                y += nodes[c].gm * nodes[c].my;
                z += nodes[c].gm * nodes[c].mz;
            }
        }
    }
    nd->gm = m;
    if (m > 0.0) {
        nd->mx = x / m;
        nd->my = y / m;
        nd->mz = z / m;
    } else {
        nd->mx = nd->my = nd->mz = 0.0;
    }
}

static void sp_kick_one(double dx, double dy, double dz, double gm,
                        double soft2, double *ax, double *ay, double *az)
{
    double r2 = dx * dx + dy * dy + dz * dz + soft2;
    if (r2 > 0.0) {
        double inv = 1.0 / sqrt(r2);
        double q = gm * inv * inv * inv;
        *ax += q * dx;
        *ay += q * dy;
        *az += q * dz;
    }
}

static int sp_cell_contains(const SpNode *nd, double x, double y,
                            double z)
{
    double h = nd->half * (1.0 + 1e-12) + 1e-30;
    return fabs(x - nd->cx) <= h && fabs(y - nd->cy) <= h &&
           fabs(z - nd->cz) <= h;
}

/* Accept cell when size < theta * (d - rho): rho bounds how far any
 * mass in the cell can sit from the centre of mass on the near side,
 * so the comparison is against the worst-case source distance. Leaves
 * sum their bodies directly; the target's own index is skipped. */
static void sp_tree_eval(const SpNode *nodes, const StcSpace *s,
                         const int *next, int idx, double px, double py,
                         double pz, int self, double theta, double soft2,
                         double *ax, double *ay, double *az)
{
    const SpNode *nd = &nodes[idx];
    if (!(nd->gm > 0.0)) {
        return;
    }
    if (sp_node_leaf(nd)) {
        int cur = nd->head;
        while (cur != -1) {
            if (cur != self) {
                sp_kick_one(s->x[cur] - px, s->y[cur] - py,
                            s->z[cur] - pz, s->gm[cur], soft2,
                            ax, ay, az);
            }
            cur = next[cur];
        }
        return;
    }
    {
        double dx = nd->mx - px;
        double dy = nd->my - py;
        double dz = nd->mz - pz;
        double d2 = dx * dx + dy * dy + dz * dz;
        double size = nd->half * 2.0;
        double ox = nd->mx - nd->cx;
        double oy = nd->my - nd->cy;
        double oz = nd->mz - nd->cz;
        double rho = nd->half * SP_SQRT3 + sqrt(ox * ox + oy * oy + oz * oz);
        double near;
        if (theta > 0.0 && d2 > 0.0 && !sp_cell_contains(nd, px, py, pz)) {
            near = sqrt(d2) - rho;
            if (near > 0.0 && size < theta * near) {
                sp_kick_one(dx, dy, dz, nd->gm, soft2, ax, ay, az);
                return;
            }
        }
    }
    {
        int k;
        for (k = 0; k < 8; k++) {
            if (nd->child[k] != -1) {
                sp_tree_eval(nodes, s, next, nd->child[k], px, py, pz,
                             self, theta, soft2, ax, ay, az);
            }
        }
    }
}

/* Build the octree over the live massive bodies into the scratch
 * arena. Returns node count, 0 when there is nothing to build. */
static int sp_tree_build(const StcSpace *s, int nlive)
{
    StcSpaceCkpt *ck = s->ck;
    int i, count, cap;
    double minx, miny, minz, maxx, maxy, maxz, cx, cy, cz, ext, ac, half;
    int any = 0;

    minx = miny = minz = maxx = maxy = maxz = 0.0;
    for (i = 0; i < s->nm; i++) {
        if (!s->alive[i]) {
            continue;
        }
        if (!any) {
            minx = maxx = s->x[i];
            miny = maxy = s->y[i];
            minz = maxz = s->z[i];
            any = 1;
        } else {
            if (s->x[i] < minx) minx = s->x[i];
            if (s->x[i] > maxx) maxx = s->x[i];
            if (s->y[i] < miny) miny = s->y[i];
            if (s->y[i] > maxy) maxy = s->y[i];
            if (s->z[i] < minz) minz = s->z[i];
            if (s->z[i] > maxz) maxz = s->z[i];
        }
    }
    if (!any) {
        return 0;
    }
    cx = (minx + maxx) * 0.5;
    cy = (miny + maxy) * 0.5;
    cz = (minz + maxz) * 0.5;
    ext = maxx - minx;
    if (maxy - miny > ext) {
        ext = maxy - miny;
    }
    if (maxz - minz > ext) {
        ext = maxz - minz;
    }
    ac = fabs(cx);
    if (fabs(cy) > ac) ac = fabs(cy);
    if (fabs(cz) > ac) ac = fabs(cz);
    half = ext * 0.5 + 1e-9 * (1.0 + ac);
    if (!(half > 0.0)) {
        half = 1.0;
    }
    cap = nlive * 8 + 64;
    if (cap > ck->tncap) {
        SpNode *nn = (SpNode *)sp_xrealloc(ck->tn, sizeof(SpNode) * (size_t)cap);
        if (nn == NULL) {
            return 0;
        }
        ck->tn = nn;
        ck->tncap = cap;
    }
    count = 1;
    sp_node_init(&ck->tn[0], cx, cy, cz, half);
    for (i = 0; i < s->nm; i++) {
        ck->tnext[i] = -1;
        if (s->alive[i]) {
            sp_insert(ck->tn, &count, ck->tncap, ck->tnext, s, 0, i, 0);
        }
    }
    sp_summarize(ck->tn, s, ck->tnext, 0);
    return count;
}

/* ---------------------------------------------------------------- */
/* force evaluation                                                  */

typedef struct {
    const StcSpace *s;
    const SpNode *nodes; /* NULL -> direct summation */
    int nsrc;
    int nobl;
    double soft2;
    int kick;            /* also midpoint-kick each tracer (PN kick) */
    double half;         /* kick half-step tau/2 when kick != 0 */
    double *ax, *ay, *az;
    double *tax, *tay, *taz;
} SpForceCtx;

/* Newtonian pair + 1PN source terms for one target; the 1PN part is
 * identical for massive and tracer targets (self is only skipped for
 * massive targets, -1 for tracers). */
static void sp_pn_add(const StcSpace *s, const int *src, int nsrc,
                      int self, double px, double py, double pz,
                      double vx, double vy, double vz,
                      double *ax, double *ay, double *az)
{
    const double c2 = STC_C * STC_C;
    int k;
    for (k = 0; k < nsrc; k++) {
        int j = src[k];
        double rx, ry, rz, rvx, rvy, rvz, r2, r, g, v2, rv, f, b1;
        if (j == self) {
            continue;
        }
        rx = px - s->x[j];
        ry = py - s->y[j];
        rz = pz - s->z[j];
        rvx = vx - s->vx[j];
        rvy = vy - s->vy[j];
        rvz = vz - s->vz[j];
        r2 = rx * rx + ry * ry + rz * rz;
        if (!(r2 > 0.0)) {
            continue;
        }
        r = sqrt(r2);
        g = s->gm[j];
        v2 = rvx * rvx + rvy * rvy + rvz * rvz;
        rv = rx * rvx + ry * rvy + rz * rvz;
        f = g / (c2 * r2 * r);
        b1 = 4.0 * g / r - v2;
        *ax += f * (b1 * rx + 4.0 * rv * rvx);
        *ay += f * (b1 * ry + 4.0 * rv * rvy);
        *az += f * (b1 * rz + 4.0 * rv * rvz);
    }
}

/* Cap on the per-particle fixed-point iteration solving v_mid. The
 * contraction per iteration is ~ (tau/2) |da_pn/dv| ~ 1e-8 for
 * Solar-System PN kicks, so the 1e-15-relative exit fires after 1-2
 * iterations; the bound only guards pathological cases (determinism:
 * a pure function of state). */
#define SP_MID_ITERS 8

/* One particle's implicit midpoint: v_mid solves
 *   v_mid = v0 + half * ( a_static + a_pn(x, v_mid) )
 * where a_static = (asx,asy,asz) (computed once per kick) and a_pn
 * comes from the few live 1PN sources in ck->src. The 1PN term reads
 * source velocities from the STATE (frozen at v_start for the whole
 * kick), so each particle's fixed point is independent — iterated
 * locally in registers, with the same stopping rule as the old joint
 * iteration (which converged every particle to the same limit).
 * Writes v_new = 2 v_mid - v0. */
static void sp_mid_one(const StcSpace *s, double half, int self,
                       double px, double py, double pz,
                       double asx, double asy, double asz,
                       double v0x, double v0y, double v0z,
                       double *outx, double *outy, double *outz)
{
    double mx, my, mz, ax, ay, az;
    double sc = sqrt(v0x * v0x + v0y * v0y + v0z * v0z);
    int it;
    ax = ay = az = 0.0;
    sp_pn_add(s, s->ck->src, s->ck->srcflag, self, px, py, pz,
              v0x, v0y, v0z, &ax, &ay, &az);
    mx = v0x + half * (asx + ax);
    my = v0y + half * (asy + ay);
    mz = v0z + half * (asz + az);
    for (it = 0; it < SP_MID_ITERS; it++) {
        double nx, ny, nz, dx, dy, dz;
        ax = ay = az = 0.0;
        sp_pn_add(s, s->ck->src, s->ck->srcflag, self, px, py, pz,
                  mx, my, mz, &ax, &ay, &az);
        nx = v0x + half * (asx + ax);
        ny = v0y + half * (asy + ay);
        nz = v0z + half * (asz + az);
        dx = nx - mx;
        dy = ny - my;
        dz = nz - mz;
        mx = nx;
        my = ny;
        mz = nz;
        if (sqrt(dx * dx + dy * dy + dz * dz) / (sc + 1e-30) <= 1e-15) {
            break;
        }
    }
    *outx = 2.0 * mx - v0x;
    *outy = 2.0 * my - v0y;
    *outz = 2.0 * mz - v0z;
}

/* Zonal-harmonic (J2/J4) acceleration on a target at displacement
 * d = (dx,dy,dz) = x_target - x_body from oblate body j.
 *   U = -(GM/r) [1 - J2 (R/r)^2 P2(u) - J4 (R/r)^4 P4(u)],
 *   u = sin(phi) = (d . p) / r,   p = unit spin pole,
 *   P2(u) = (3u^2-1)/2,  P4(u) = (35u^4-30u^2+3)/8.
 * With du/dd = (p - u rhat)/r and rhat = d/r, the gradient of one
 * zonal term U_n = GM Jn R^n r^{-(n+1)} P_n(u) is
 *   a_n = -grad U_n
 *       = GM Jn R^n r^{-(n+2)} [ ((n+1) P_n + u P_n') rhat - P_n' p ]
 * so, with P2' = 3u and P4' = (5/2) u (7u^2 - 3),
 *   a_J2 = GM J2 R^2 r^-4 [ (15u^2 - 3)/2 rhat - 3u p ],
 *   a_J4 = GM J4 R^4 r^-6 [ (15/8)(21u^4-14u^2+1) rhat
 *                          - (5/2) u (7u^2-3) p ].
 * (Pole check: at the equator u=0 -> a_J2 = -(3/2)k rhat, extra pull
 * toward the bulge; at the pole u=1 -> a_J2 = +3k rhat, reduced pull.
 * Both match the standard zonal formula with p = z-hat.) */
static void sp_zonal(const StcSpace *s, int j, double dx, double dy,
                     double dz, double *ax, double *ay, double *az)
{
    double r2 = dx * dx + dy * dy + dz * dz;
    double r, u, rx, ry, rz, px, py, pz, g, R;
    if (!(r2 > 0.0)) {
        return;
    }
    r = sqrt(r2);
    px = s->polx[j];
    py = s->poly[j];
    pz = s->polz[j];
    u = (dx * px + dy * py + dz * pz) / r;
    rx = dx / r;
    ry = dy / r;
    rz = dz / r;
    g = s->gm[j];
    R = s->rref[j];
    if (s->j2[j] != 0.0) {
        double k = g * s->j2[j] * R * R / (r2 * r2);
        double br = (15.0 * u * u - 3.0) * 0.5;
        double bp = -3.0 * u;
        *ax += k * (br * rx + bp * px);
        *ay += k * (br * ry + bp * py);
        *az += k * (br * rz + bp * pz);
    }
    if (s->j4[j] != 0.0) {
        double u2 = u * u;
        double R2 = R * R;
        double k = g * s->j4[j] * R2 * R2 / (r2 * r2 * r2);
        double br = (15.0 / 8.0) * (21.0 * u2 * u2 - 14.0 * u2 + 1.0);
        double bp = -2.5 * u * (7.0 * u2 - 3.0);
        *ax += k * (br * rx + bp * px);
        *ay += k * (br * ry + bp * py);
        *az += k * (br * rz + bp * pz);
    }
}

static void sp_fmassive(void *vctx, int i0, int i1)
{
    const SpForceCtx *c = (const SpForceCtx *)vctx;
    const StcSpace *s = c->s;
    int i, j;
    for (i = i0; i < i1; i++) {
        double axx = 0.0, ayy = 0.0, azz = 0.0;
        if (s->alive[i]) {
            double xi = s->x[i], yi = s->y[i], zi = s->z[i];
            if (c->nodes != NULL) {
                sp_tree_eval(c->nodes, s, s->ck->tnext, 0,
                             xi, yi, zi, i, s->theta, c->soft2,
                             &axx, &ayy, &azz);
            } else {
                for (j = 0; j < s->nm; j++) {
                    if (j != i && s->alive[j]) {
                        double dx = s->x[j] - xi;
                        double dy = s->y[j] - yi;
                        double dz = s->z[j] - zi;
                        double r2 = dx * dx + dy * dy + dz * dz + c->soft2;
                        if (r2 > 0.0) {
                            double inv = 1.0 / sqrt(r2);
                            double q = s->gm[j] * inv * inv * inv;
                            axx += q * dx;
                            ayy += q * dy;
                            azz += q * dz;
                        }
                    }
                }
            }
            if (c->nsrc > 0) {
                sp_pn_add(s, s->ck->src, c->nsrc, i, xi, yi, zi,
                          s->vx[i], s->vy[i], s->vz[i],
                          &axx, &ayy, &azz);
            }
            if (c->nobl > 0) {
                int k;
                for (k = 0; k < c->nobl; k++) {
                    int j = s->ck->obl[k];
                    if (j != i) {
                        sp_zonal(s, j, xi - s->x[j], yi - s->y[j],
                                 zi - s->z[j], &axx, &ayy, &azz);
                    }
                }
            }
        }
        c->ax[i] = axx;
        c->ay[i] = ayy;
        c->az[i] = azz;
    }
}

static void sp_ftracer(void *vctx, int i0, int i1)
{
    const SpForceCtx *c = (const SpForceCtx *)vctx;
    const StcSpace *s = c->s;
    int i, j;
    for (i = i0; i < i1; i++) {
        double axx = 0.0, ayy = 0.0, azz = 0.0;
        if (s->talive[i]) {
            double xi = s->tx[i], yi = s->ty[i], zi = s->tz[i];
            if (c->nodes != NULL) {
                sp_tree_eval(c->nodes, s, s->ck->tnext, 0,
                             xi, yi, zi, -1, s->theta, c->soft2,
                             &axx, &ayy, &azz);
            } else {
                for (j = 0; j < s->nm; j++) {
                    if (s->alive[j]) {
                        double dx = s->x[j] - xi;
                        double dy = s->y[j] - yi;
                        double dz = s->z[j] - zi;
                        double r2 = dx * dx + dy * dy + dz * dz + c->soft2;
                        if (r2 > 0.0) {
                            double inv = 1.0 / sqrt(r2);
                            double q = s->gm[j] * inv * inv * inv;
                            axx += q * dx;
                            ayy += q * dy;
                            azz += q * dz;
                        }
                    }
                }
            }
            if (c->nsrc > 0) {
                sp_pn_add(s, s->ck->src, c->nsrc, -1, xi, yi, zi,
                          s->tvx[i], s->tvy[i], s->tvz[i],
                          &axx, &ayy, &azz);
            }
            if (c->nobl > 0) {
                int k;
                for (k = 0; k < c->nobl; k++) {
                    int j = s->ck->obl[k];
                    sp_zonal(s, j, xi - s->x[j], yi - s->y[j],
                             zi - s->z[j], &axx, &ayy, &azz);
                }
            }
            /* PN kick: the tracer's static accel is complete here, so
             * its implicit-midpoint update runs in place (tracers are
             * massless — no reaction to wait for). */
            if (c->kick) {
                sp_mid_one(s, c->half, -1, xi, yi, zi, axx, ayy, azz,
                           s->tvx[i], s->tvy[i], s->tvz[i],
                           s->tvx + i, s->tvy + i, s->tvz + i);
            }
        }
        c->tax[i] = axx;
        c->tay[i] = ayy;
        c->taz[i] = azz;
    }
}

typedef struct {
    const StcSpace *s;
    int j;             /* oblate body receiving the reaction */
    double *part3;     /* out: per-index gm_i * a_i^zonal (3*i) */
} SpReactCtx;

static void sp_react_fn(void *vctx, int i0, int i1)
{
    const SpReactCtx *c = (const SpReactCtx *)vctx;
    const StcSpace *s = c->s;
    int i;
    for (i = i0; i < i1; i++) {
        double w = 0.0, axx = 0.0, ayy = 0.0, azz = 0.0;
        if (s->alive[i] && i != c->j) {
            sp_zonal(s, c->j, s->x[i] - s->x[c->j],
                     s->y[i] - s->y[c->j], s->z[i] - s->z[c->j],
                     &axx, &ayy, &azz);
            w = s->gm[i];
        }
        c->part3[3 * i] = w * axx;
        c->part3[3 * i + 1] = w * ayy;
        c->part3[3 * i + 2] = w * azz;
    }
}

/* Equal and opposite zonal reaction on each oblate body j:
 * m_j a_j += -sum_i m_i a_i^zonal(i,j) over MASSIVE targets (tracers
 * are massless, no reaction). Parallel fill of per-index partials,
 * then a fixed-order serial sum -> deterministic (R4). */
static void sp_zonal_react(const StcSpace *s)
{
    StcSpaceCkpt *ck = s->ck;
    SpReactCtx c;
    int k, i;
    c.s = s;
    c.part3 = ck->part3;
    for (k = 0; k < ck->oblflag; k++) {
        double sx = 0.0, sy = 0.0, sz = 0.0;
        int j = ck->obl[k];
        c.j = j;
        stc_par_for(sp_react_fn, &c, s->nm, 0);
        for (i = 0; i < s->nm; i++) {
            sx += ck->part3[3 * i];
            sy += ck->part3[3 * i + 1];
            sz += ck->part3[3 * i + 2];
        }
        ck->ax[j] -= sx / s->gm[j];
        ck->ay[j] -= sy / s->gm[j];
        ck->az[j] -= sz / s->gm[j];
    }
}

/* Refresh the live-1PN-source and oblate-body index lists (ck->src,
 * ck->obl) and return the live massive count. O(nm), no force work.
 * Caller must have fitted the scratch first. */
static int sp_scan(const StcSpace *s)
{
    StcSpaceCkpt *ck = s->ck;
    int i, nlive = 0;
    ck->srcflag = 0;
    ck->oblflag = 0;
    for (i = 0; i < s->nm; i++) {
        if (s->alive[i]) {
            nlive++;
            if (s->pn && sp_is_pn_source(s->kind[i])) {
                ck->src[ck->srcflag++] = i;
            }
            if ((s->j2[i] != 0.0 || s->j4[i] != 0.0) &&
                s->rref[i] > 0.0 &&
                s->polx[i] * s->polx[i] + s->poly[i] * s->poly[i] +
                    s->polz[i] * s->polz[i] > 0.0) {
                ck->obl[ck->oblflag++] = i;
            }
        }
    }
    return nlive;
}

/* Evaluate the force model into the scratch acceleration arrays
 * (dead slots get 0). with_pn selects the full model (Newton + zonal
 * + 1PN) or the position-only static part (Newton + zonal). The
 * velocity-dependent 1PN term always reads the state velocities here;
 * the implicit-midpoint kick evaluates a_pn at trial velocities inside
 * sp_mid_one instead. kick != 0 additionally applies the midpoint
 * update to each tracer in place (valid only with with_pn == 0:
 * tracers have no zonal reaction, so their a_static is final). */
static void sp_forces_v(const StcSpace *s, int with_pn, int kick,
                        double half)
{
    StcSpaceCkpt *ck = s->ck;
    SpForceCtx c;
    int nlive;

    if (ck == NULL || sp_fit_scratch(s) != 0) {
        return; /* allocation failure: keep previous/zero accels */
    }
    memset(ck->ax, 0, sizeof(double) * (size_t)s->cap_m);
    memset(ck->ay, 0, sizeof(double) * (size_t)s->cap_m);
    memset(ck->az, 0, sizeof(double) * (size_t)s->cap_m);
    memset(ck->tax, 0, sizeof(double) * (size_t)s->cap_t);
    memset(ck->tay, 0, sizeof(double) * (size_t)s->cap_t);
    memset(ck->taz, 0, sizeof(double) * (size_t)s->cap_t);

    nlive = sp_scan(s);

    c.s = s;
    c.nodes = NULL;
    c.nsrc = with_pn ? ck->srcflag : 0;
    c.nobl = ck->oblflag;
    c.soft2 = s->soft * s->soft;
    c.kick = kick;
    c.half = half;
    c.ax = ck->ax; c.ay = ck->ay; c.az = ck->az;
    c.tax = ck->tax; c.tay = ck->tay; c.taz = ck->taz;

    if (nlive > s->tree_threshold && sp_tree_build(s, nlive) > 0) {
        c.nodes = ck->tn;
    }
    stc_par_for(sp_fmassive, &c, s->nm, 0);
    stc_par_for(sp_ftracer, &c, s->nt, 0);
    if (ck->oblflag > 0) {
        sp_zonal_react(s);
    }
    ck->avalid = with_pn ? 1 : 2;
}

static void sp_forces(const StcSpace *s)
{
    sp_forces_v(s, 1, 0, 0.0);
}

/* ---------------------------------------------------------------- */
/* kick / drift / substep                                            */

typedef struct {
    StcSpace *s;
    double tau;
} SpKdCtx;

static void sp_kick_m(void *vctx, int i0, int i1)
{
    const SpKdCtx *c = (const SpKdCtx *)vctx;
    StcSpace *s = c->s;
    const StcSpaceCkpt *ck = s->ck;
    int i;
    for (i = i0; i < i1; i++) {
        if (s->alive[i]) {
            s->vx[i] += ck->ax[i] * c->tau;
            s->vy[i] += ck->ay[i] * c->tau;
            s->vz[i] += ck->az[i] * c->tau;
        }
    }
}

static void sp_kick_t(void *vctx, int i0, int i1)
{
    const SpKdCtx *c = (const SpKdCtx *)vctx;
    StcSpace *s = c->s;
    const StcSpaceCkpt *ck = s->ck;
    int i;
    for (i = i0; i < i1; i++) {
        if (s->talive[i]) {
            s->tvx[i] += ck->tax[i] * c->tau;
            s->tvy[i] += ck->tay[i] * c->tau;
            s->tvz[i] += ck->taz[i] * c->tau;
        }
    }
}

static void sp_drift_m(void *vctx, int i0, int i1)
{
    const SpKdCtx *c = (const SpKdCtx *)vctx;
    StcSpace *s = c->s;
    int i;
    for (i = i0; i < i1; i++) {
        if (s->alive[i]) {
            s->x[i] += s->vx[i] * c->tau;
            s->y[i] += s->vy[i] * c->tau;
            s->z[i] += s->vz[i] * c->tau;
        }
    }
}

static void sp_drift_t(void *vctx, int i0, int i1)
{
    const SpKdCtx *c = (const SpKdCtx *)vctx;
    StcSpace *s = c->s;
    int i;
    for (i = i0; i < i1; i++) {
        if (s->talive[i]) {
            s->tx[i] += s->tvx[i] * c->tau;
            s->ty[i] += s->tvy[i] * c->tau;
            s->tz[i] += s->tvz[i] * c->tau;
        }
    }
}

/* --- implicit-midpoint kick (velocity-dependent forces: 1PN) ------ */

typedef struct {
    StcSpace *s;
    double half;    /* tau/2 */
} SpMidCtx;

/* Massive pass: writes the new velocity directly for non-sources.
 * A live 1PN source's update is staged in ck->mvm and committed only
 * after every particle has iterated, so source velocities read as
 * v_start throughout (thread-count independent, R4). */
static void sp_mid_m(void *vctx, int i0, int i1)
{
    const SpMidCtx *c = (const SpMidCtx *)vctx;
    StcSpace *s = c->s;
    StcSpaceCkpt *ck = s->ck;
    int i;
    for (i = i0; i < i1; i++) {
        if (s->alive[i]) {
            if (sp_is_pn_source(s->kind[i])) {
                sp_mid_one(s, c->half, i, s->x[i], s->y[i], s->z[i],
                           ck->ax[i], ck->ay[i], ck->az[i],
                           s->vx[i], s->vy[i], s->vz[i],
                           ck->mvm + 3 * i, ck->mvm + 3 * i + 1,
                           ck->mvm + 3 * i + 2);
            } else {
                sp_mid_one(s, c->half, i, s->x[i], s->y[i], s->z[i],
                           ck->ax[i], ck->ay[i], ck->az[i],
                           s->vx[i], s->vy[i], s->vz[i],
                           s->vx + i, s->vy + i, s->vz + i);
            }
        }
    }
}

static void sp_mid_t(void *vctx, int i0, int i1)
{
    const SpMidCtx *c = (const SpMidCtx *)vctx;
    StcSpace *s = c->s;
    StcSpaceCkpt *ck = s->ck;
    int i;
    for (i = i0; i < i1; i++) {
        if (s->talive[i]) {
            sp_mid_one(s, c->half, -1, s->tx[i], s->ty[i], s->tz[i],
                       ck->tax[i], ck->tay[i], ck->taz[i],
                       s->tvx[i], s->tvy[i], s->tvz[i],
                       s->tvx + i, s->tvy + i, s->tvz + i);
        }
    }
}

/* Implicit-midpoint kick for the velocity-dependent 1PN force.
 * The position-only part a_static(x) (Newton + zonal) is evaluated at
 * most once per kick — reused outright when avalid == 2 — and each
 * particle's fixed point iterates in registers over just the PN
 * sources, so a PN kick costs one static force evaluation plus a
 * few-flop scalar iteration per particle. a_static does not depend on
 * velocity, so the cache survives the kick (avalid stays 2) and can
 * serve the first kick of the next substep when nothing moved.
 * ck->srcflag must be current (> 0 here). */
static void sp_kick_mid(StcSpace *s, double tau)
{
    StcSpaceCkpt *ck = s->ck;
    SpMidCtx c;
    int k;
    c.s = s;
    c.half = 0.5 * tau;
    if (ck->avalid != 2) {
        /* fresh static evaluation; the tracer midpoint updates fuse
         * into the same pass (their accel needs no later fix-up).
         * Massive bodies wait for the zonal reaction below. */
        sp_forces_v(s, 0, 1, c.half);
    } else {
        stc_par_for(sp_mid_t, &c, s->nt, 0);
    }
    stc_par_for(sp_mid_m, &c, s->nm, 0);
    /* commit the staged source velocities (few, serial) */
    for (k = 0; k < ck->srcflag; k++) {
        int j = ck->src[k];
        s->vx[j] = ck->mvm[3 * j];
        s->vy[j] = ck->mvm[3 * j + 1];
        s->vz[j] = ck->mvm[3 * j + 2];
    }
}

static void sp_kick(StcSpace *s, double tau)
{
    StcSpaceCkpt *ck = s->ck;
    SpKdCtx c;
    if (ck == NULL) {
        return;
    }
    /* Velocity-dependent 1PN kick: the source scan is O(nm) and needs
     * no force evaluation; the midpoint solver needs only the cached
     * static acceleration at the current positions. */
    if (s->pn && sp_fit_scratch(s) == 0 && ck->mvm != NULL) {
        sp_scan(s);
        if (ck->srcflag > 0) {
            sp_kick_mid(s, tau);
            return;
        }
    }
    if (!ck->avalid) {
        sp_forces(s);
    }
    c.s = s;
    c.tau = tau;
    stc_par_for(sp_kick_m, &c, s->nm, 0);
    stc_par_for(sp_kick_t, &c, s->nt, 0);
    /* with a live 1PN source the acceleration is velocity-dependent:
     * the cache is stale the moment any velocity changes. */
    if (s->pn && ck->srcflag > 0) {
        ck->avalid = 0;
    }
}

static void sp_drift(StcSpace *s, double tau)
{
    SpKdCtx c;
    c.s = s;
    c.tau = tau;
    stc_par_for(sp_drift_m, &c, s->nm, 0);
    stc_par_for(sp_drift_t, &c, s->nt, 0);
}

/* One 6th-order Yoshida substep of signed length tau. */
static void sp_substep(StcSpace *s, double tau)
{
    int k;
    for (k = 0; k < 8; k++) {
        sp_kick(s, sp_yc[k] * tau); /* evaluates forces when needed */
        if (k < 7) {
            sp_drift(s, sp_yw[k] * tau);
            if (s->ck) {
                s->ck->avalid = 0;
            }
        }
    }
    s->steps++;
}

/* ---------------------------------------------------------------- */
/* events                                                            */

static void sp_evlog(StcSpace *s, int kind, int a, int b)
{
    StcSpaceEvent *e = &s->ev[s->nev % STC_SPACE_EVLOG];
    e->kind = kind;
    e->a = a;
    e->b = b;
    e->t = s->t;
    s->nev++;
}

/* Minimum squared distance from the origin to the segment A->B. Used
 * for swept event detection: A = relative position at substep start,
 * B = at substep end (uniform relative motion approximation), so a
 * body cannot tunnel through an interaction sphere between substeps. */
static double sp_seg_dist2(double ax, double ay, double az,
                           double bx, double by, double bz)
{
    double dx = bx - ax, dy = by - ay, dz = bz - az;
    double dd = dx * dx + dy * dy + dz * dz;
    double u = 0.0, px, py, pz;
    if (dd > 0.0) {
        u = -(ax * dx + ay * dy + az * dz) / dd;
        if (u < 0.0) u = 0.0;
        if (u > 1.0) u = 1.0;
    }
    px = ax + u * dx;
    py = ay + u * dy;
    pz = az + u * dz;
    return px * px + py * py + pz * pz;
}

/* a absorbs b: conserve mass, momentum, centre of mass. Radius is the
 * summed-volume sphere; if either is a black hole the survivor is a
 * black hole with radius r_s = 2GM/c^2. kind 2 -> CAPTURE log line. */
static void sp_merge_into(StcSpace *s, int a, int b, int logkind)
{
    double gn = s->gm[a] + s->gm[b];
    double fa = s->gm[a] / gn;
    double fb = s->gm[b] / gn;
    s->x[a] = fa * s->x[a] + fb * s->x[b];
    s->y[a] = fa * s->y[a] + fb * s->y[b];
    s->z[a] = fa * s->z[a] + fb * s->z[b];
    s->vx[a] = fa * s->vx[a] + fb * s->vx[b];
    s->vy[a] = fa * s->vy[a] + fb * s->vy[b];
    s->vz[a] = fa * s->vz[a] + fb * s->vz[b];
    s->gm[a] = gn;
    if (sp_is_bh_kind(s->kind[a]) || sp_is_bh_kind(s->kind[b])) {
        s->kind[a] = STC_SK_BLACKHOLE;
        s->radius[a] = 2.0 * gn / (STC_C * STC_C);
    } else {
        double r3 = s->radius[a] * s->radius[a] * s->radius[a] +
                    s->radius[b] * s->radius[b] * s->radius[b];
        s->radius[a] = cbrt(r3);
    }
    s->alive[b] = 0;
    sp_evlog(s, logkind, a, b);
}

/* Tidally disrupt b near black hole bh: remove b, spawn
 * disrupt_tracers tracers uniformly inside b's volume carrying b's
 * velocity; stream seeded by (event index, body index). */
static void sp_disrupt_into(StcSpace *s, int bh, int b)
{
    StcRng64 rng;
    int k;
    uint64_t seed = (uint64_t)s->nev * 0x9E3779B97F4A7C15ULL +
                    (uint64_t)(b + 1) * 0x85EBCA6BULL;
    stc_rng64_seed(&rng, seed);
    sp_evlog(s, STC_SEV_DISRUPT, bh, b);
    for (k = 0; k < s->disrupt_tracers && s->nt < s->cap_t; k++) {
        double dx = stc_rng64_normal(&rng);
        double dy = stc_rng64_normal(&rng);
        double dz = stc_rng64_normal(&rng);
        double dn = dx * dx + dy * dy + dz * dz;
        double u = stc_rng64_uniform(&rng);
        int t = s->nt++;
        if (dn > 0.0) {
            double f = cbrt(u) * s->radius[b] / sqrt(dn);
            s->tx[t] = s->x[b] + dx * f;
            s->ty[t] = s->y[b] + dy * f;
            s->tz[t] = s->z[b] + dz * f;
        } else {
            s->tx[t] = s->x[b];
            s->ty[t] = s->y[b];
            s->tz[t] = s->z[b];
        }
        s->tvx[t] = s->vx[b];
        s->tvy[t] = s->vy[b];
        s->tvz[t] = s->vz[b];
        s->talive[t] = 1;
    }
    s->alive[b] = 0;
}

/* Massive-massive contacts: sweep live bodies sorted by x and collect
 * candidate pairs inside a conservative reach window, then resolve
 * them in strict (i,j) index order against the live state. tau is the
 * substep just taken; distances are swept over it (uniform relative
 * motion) so fast bodies cannot tunnel through contacts. */
static void sp_events_massive(StcSpace *s, double tau)
{
    StcSpaceCkpt *ck = s->ck;
    SpEnt *ord;
    int i, j, n, nl = 0, any = 0;
    double gmin = 0.0, rmax = 0.0, reach_max = 0.0;
    double vmax = 0.0, vpad;
    long npair = 0;

    if (ck == NULL || sp_fit_scratch(s) != 0) {
        return;
    }
    for (i = 0; i < s->nm; i++) {
        if (s->alive[i]) {
            double v2 = s->vx[i] * s->vx[i] + s->vy[i] * s->vy[i] +
                        s->vz[i] * s->vz[i];
            if (!any || s->gm[i] < gmin) gmin = s->gm[i];
            if (!any || s->radius[i] > rmax) rmax = s->radius[i];
            if (!any || v2 > vmax) vmax = v2;
            any = 1;
            nl++;
        }
    }
    if (nl < 2) {
        return;
    }
    vmax = sqrt(vmax);
    vpad = 2.0 * vmax * fabs(tau); /* bound on relative displacement */
    for (i = 0; i < s->nm; i++) {
        double r = 0.0;
        if (s->alive[i]) {
            if (sp_is_bh_kind(s->kind[i])) {
                double cap = STC_SPACE_CAPTURE_RS * s->radius[i];
                double rt = rmax * cbrt(s->gm[i] / gmin);
                r = cap > rt ? cap : rt;
            } else {
                r = s->radius[i];
            }
        }
        ck->reach[i] = r;
        if (r > reach_max) {
            reach_max = r;
        }
    }
    if (!(reach_max > 0.0)) {
        return; /* no contacts possible */
    }
    ord = (SpEnt *)malloc(sizeof(SpEnt) * (size_t)nl);
    if (ord == NULL) {
        return;
    }
    n = 0;
    for (i = 0; i < s->nm; i++) {
        if (s->alive[i]) {
            ord[n].key = s->x[i];
            ord[n].idx = i;
            n++;
        }
    }
    qsort(ord, (size_t)n, sizeof(SpEnt), sp_ent_cmp);
    for (i = 0; i < n; i++) {
        int bi = ord[i].idx;
        for (j = i + 1; j < n; j++) {
            int bj = ord[j].idx;
            double dx = s->x[bj] - s->x[bi];
            double w = ck->reach[bi] + ck->reach[bj] + vpad;
            double dy, dz, d2;
            if (dx > ck->reach[bi] + reach_max + vpad) {
                break;
            }
            if (dx > w) {
                continue;
            }
            dy = s->y[bj] - s->y[bi];
            if (fabs(dy) > w) {
                continue;
            }
            dz = s->z[bj] - s->z[bi];
            if (fabs(dz) > w) {
                continue;
            }
            d2 = dx * dx + dy * dy + dz * dz;
            if (d2 <= w * w) {
                uint64_t key;
                int lo = bi < bj ? bi : bj;
                int hi = bi < bj ? bj : bi;
                if (npair >= ck->paircap) {
                    long ncap = ck->paircap ? ck->paircap * 2 : 256;
                    uint64_t *np = (uint64_t *)sp_xrealloc(
                        ck->pairs, sizeof(uint64_t) * (size_t)ncap);
                    if (np == NULL) {
                        break;
                    }
                    ck->pairs = np;
                    ck->paircap = ncap;
                }
                key = ((uint64_t)(uint32_t)lo << 32) | (uint32_t)hi;
                ck->pairs[npair++] = key;
            }
        }
    }
    free(ord);
    if (npair == 0) {
        return;
    }
    qsort(ck->pairs, (size_t)npair, sizeof(uint64_t), sp_u64_cmp);
    for (i = 0; i < npair; i++) {
        int a = (int)(ck->pairs[i] >> 32);
        int b = (int)(ck->pairs[i] & 0xffffffffu);
        double dx, dy, dz, d, cap;
        double dvx, dvy, dvz;
        int bi, bj, bh, o;
        if (!s->alive[a] || !s->alive[b]) {
            continue;
        }
        /* swept relative position: end minus uniform displacement */
        dx = s->x[b] - s->x[a];
        dy = s->y[b] - s->y[a];
        dz = s->z[b] - s->z[a];
        dvx = (s->vx[b] - s->vx[a]) * tau;
        dvy = (s->vy[b] - s->vy[a]) * tau;
        dvz = (s->vz[b] - s->vz[a]) * tau;
        d = sqrt(sp_seg_dist2(dx - dvx, dy - dvy, dz - dvz,
                              dx, dy, dz));
        bi = sp_is_bh_kind(s->kind[a]);
        bj = sp_is_bh_kind(s->kind[b]);
        if (bi && bj) {
            int heavy = s->gm[a] >= s->gm[b] ? a : b;
            cap = STC_SPACE_CAPTURE_RS *
                  (s->radius[a] > s->radius[b] ? s->radius[a] : s->radius[b]);
            if (d < cap) {
                sp_merge_into(s, heavy, heavy == a ? b : a,
                              STC_SEV_CAPTURE);
                if (ck) ck->avalid = 0;
            }
            continue;
        }
        if (bi || bj) {
            bh = bi ? a : b;
            o = bi ? b : a;
            cap = STC_SPACE_CAPTURE_RS * s->radius[bh];
            if (d < cap) {
                sp_merge_into(s, bh, o, STC_SEV_CAPTURE);
                ck->avalid = 0;
            } else {
                double rt = s->radius[o] * cbrt(s->gm[bh] / s->gm[o]);
                if (rt > cap && d < rt) {
                    sp_disrupt_into(s, bh, o);
                    ck->avalid = 0;
                } else if (d < s->radius[a] + s->radius[b]) {
                    int heavy = s->gm[a] >= s->gm[b] ? a : b;
                    sp_merge_into(s, heavy, heavy == a ? b : a,
                                  STC_SEV_MERGE);
                    ck->avalid = 0;
                }
            }
            continue;
        }
        if (d < s->radius[a] + s->radius[b]) {
            int heavy = s->gm[a] >= s->gm[b] ? a : b;
            sp_merge_into(s, heavy, heavy == a ? b : a, STC_SEV_MERGE);
            ck->avalid = 0;
        }
    }
}

typedef struct {
    const StcSpace *s;
    const SpEnt *ord;   /* radius-bearing bodies sorted by x */
    int n;
    double wmax;        /* max event radius */
    double mvpad;       /* max massive-body displacement over tau */
    double tau;
    int *hit;
} SpHitCtx;

static void sp_scan_tracers(void *vctx, int i0, int i1)
{
    const SpHitCtx *c = (const SpHitCtx *)vctx;
    const StcSpace *s = c->s;
    int t;
    for (t = i0; t < i1; t++) {
        int hit = -1;
        if (s->talive[t]) {
            double xt = s->tx[t], yt = s->ty[t], zt = s->tz[t];
            double v2 = s->tvx[t] * s->tvx[t] + s->tvy[t] * s->tvy[t] +
                        s->tvz[t] * s->tvz[t];
            double w = c->wmax + (sqrt(v2) * fabs(c->tau) + c->mvpad);
            double lo = xt - w, hi = xt + w;
            int a = 0, b = c->n, k;
            while (a < b) {
                int mid = (a + b) >> 1;
                if (c->ord[mid].key < lo) a = mid + 1; else b = mid;
            }
            for (k = a; k < c->n && c->ord[k].key <= hi; k++) {
                int m = c->ord[k].idx;
                double evr = sp_is_bh_kind(s->kind[m])
                                 ? STC_SPACE_CAPTURE_RS * s->radius[m]
                                 : s->radius[m];
                /* relative segment: end = body - tracer (current),
                 * start = end - (v_body - v_tracer) * tau */
                double dx = s->x[m] - xt;
                double dy = s->y[m] - yt;
                double dz = s->z[m] - zt;
                double dvx = (s->vx[m] - s->tvx[t]) * c->tau;
                double dvy = (s->vy[m] - s->tvy[t]) * c->tau;
                double dvz = (s->vz[m] - s->tvz[t]) * c->tau;
                double d2 = sp_seg_dist2(dx - dvx, dy - dvy, dz - dvz,
                                         dx, dy, dz);
                if (d2 <= evr * evr && (hit < 0 || m < hit)) {
                    hit = m;
                }
            }
        }
        c->hit[t] = hit;
    }
}

/* Tracer contacts: detected in parallel (each tracer writes only its
 * own hit slot; lowest massive index wins on overlap), applied
 * sequentially in tracer index order. */
static void sp_events_tracer(StcSpace *s, double tau)
{
    StcSpaceCkpt *ck = s->ck;
    SpEnt *ord;
    SpHitCtx c;
    int i, t, n = 0;
    double wmax = 0.0, mvmax = 0.0;

    if (ck == NULL || s->nt == 0 || sp_fit_scratch(s) != 0) {
        return;
    }
    for (i = 0; i < s->nm; i++) {
        if (s->alive[i]) {
            double evr = sp_is_bh_kind(s->kind[i])
                             ? STC_SPACE_CAPTURE_RS * s->radius[i]
                             : s->radius[i];
            double v2 = s->vx[i] * s->vx[i] + s->vy[i] * s->vy[i] +
                        s->vz[i] * s->vz[i];
            if (v2 > mvmax) {
                mvmax = v2;
            }
            if (evr > 0.0) {
                ck->ilist[n++] = i;
                if (evr > wmax) {
                    wmax = evr;
                }
            }
        }
    }
    if (n == 0) {
        return;
    }
    ord = (SpEnt *)malloc(sizeof(SpEnt) * (size_t)n);
    if (ord == NULL) {
        return;
    }
    for (i = 0; i < n; i++) {
        ord[i].key = s->x[ck->ilist[i]];
        ord[i].idx = ck->ilist[i];
    }
    qsort(ord, (size_t)n, sizeof(SpEnt), sp_ent_cmp);
    c.s = s;
    c.ord = ord;
    c.n = n;
    c.wmax = wmax;
    c.mvpad = sqrt(mvmax) * fabs(tau);
    c.tau = tau;
    c.hit = ck->thit;
    stc_par_for(sp_scan_tracers, &c, s->nt, 0);
    free(ord);
    for (t = 0; t < s->nt; t++) {
        int m = ck->thit[t];
        if (m >= 0 && s->talive[t] && s->alive[m]) {
            if (sp_is_bh_kind(s->kind[m])) {
                sp_evlog(s, STC_SEV_CAPTURE, m, -1 - t);
            } else {
                sp_evlog(s, STC_SEV_ABSORB, m, t);
            }
            s->talive[t] = 0;
            ck->avalid = 0;
        }
    }
}

static void sp_events(StcSpace *s, double tau)
{
    if (!s->collisions) {
        return;
    }
    sp_events_massive(s, tau);
    sp_events_tracer(s, tau);
}

/* ---------------------------------------------------------------- */
/* checkpoints                                                       */

static int sp_snap_fit(SpSnap *sn, int nm, int nt)
{
    if (nm > sn->mcap) {
        double *mv = (double *)sp_xrealloc(sn->mv, sizeof(double) * 14u * (size_t)nm);
        int *mi = (int *)sp_xrealloc(sn->mi, sizeof(int) * 2u * (size_t)nm);
        unsigned char *ma = (unsigned char *)sp_xrealloc(sn->ma, (size_t)nm);
        char *mn = (char *)sp_xrealloc(sn->mn, (size_t)nm * STC_SPACE_NAME);
        if (!mv || !mi || !ma || !mn) {
            free(mv); free(mi); free(ma); free(mn);
            return -2;
        }
        sn->mv = mv; sn->mi = mi; sn->ma = ma; sn->mn = mn;
        sn->mcap = nm;
    }
    if (nt > sn->tcap) {
        double *tv = (double *)sp_xrealloc(sn->tv, sizeof(double) * 6u * (size_t)nt);
        unsigned char *ta = (unsigned char *)sp_xrealloc(sn->ta, (size_t)nt);
        if (!tv || !ta) {
            free(tv); free(ta);
            return -2;
        }
        sn->tv = tv; sn->ta = ta;
        sn->tcap = nt;
    }
    return 0;
}

static void sp_snap_save(const StcSpace *s, SpSnap *sn)
{
    size_t nm = (size_t)s->nm, nt = (size_t)s->nt;
    double *mv = sn->mv;
    int *mi = sn->mi;
    memcpy(mv + 0 * nm, s->x, nm * sizeof(double));
    memcpy(mv + 1 * nm, s->y, nm * sizeof(double));
    memcpy(mv + 2 * nm, s->z, nm * sizeof(double));
    memcpy(mv + 3 * nm, s->vx, nm * sizeof(double));
    memcpy(mv + 4 * nm, s->vy, nm * sizeof(double));
    memcpy(mv + 5 * nm, s->vz, nm * sizeof(double));
    memcpy(mv + 6 * nm, s->gm, nm * sizeof(double));
    memcpy(mv + 7 * nm, s->radius, nm * sizeof(double));
    memcpy(mv + 8 * nm, s->j2, nm * sizeof(double));
    memcpy(mv + 9 * nm, s->j4, nm * sizeof(double));
    memcpy(mv + 10 * nm, s->rref, nm * sizeof(double));
    memcpy(mv + 11 * nm, s->polx, nm * sizeof(double));
    memcpy(mv + 12 * nm, s->poly, nm * sizeof(double));
    memcpy(mv + 13 * nm, s->polz, nm * sizeof(double));
    memcpy(mi + 0 * nm, s->kind, nm * sizeof(int));
    memcpy(mi + 1 * nm, s->parent, nm * sizeof(int));
    memcpy(sn->ma, s->alive, nm);
    memcpy(sn->mn, s->name, nm * STC_SPACE_NAME);
    memcpy(sn->tv + 0 * nt, s->tx, nt * sizeof(double));
    memcpy(sn->tv + 1 * nt, s->ty, nt * sizeof(double));
    memcpy(sn->tv + 2 * nt, s->tz, nt * sizeof(double));
    memcpy(sn->tv + 3 * nt, s->tvx, nt * sizeof(double));
    memcpy(sn->tv + 4 * nt, s->tvy, nt * sizeof(double));
    memcpy(sn->tv + 5 * nt, s->tvz, nt * sizeof(double));
    memcpy(sn->ta, s->talive, nt);
    sn->t = s->t;
    sn->steps = s->steps;
    sn->nm = s->nm;
    sn->nt = s->nt;
    sn->nev = s->nev;
    sn->e0 = s->e0;
    memcpy(sn->ev, s->ev, sizeof(sn->ev));
}

static int sp_snap_load(StcSpace *s, const SpSnap *sn)
{
    size_t nm = (size_t)sn->nm, nt = (size_t)sn->nt;
    const double *mv = sn->mv;
    const int *mi = sn->mi;
    if (sn->nm > s->cap_m || sn->nt > s->cap_t) {
        return -1;
    }
    memcpy(s->x, mv + 0 * nm, nm * sizeof(double));
    memcpy(s->y, mv + 1 * nm, nm * sizeof(double));
    memcpy(s->z, mv + 2 * nm, nm * sizeof(double));
    memcpy(s->vx, mv + 3 * nm, nm * sizeof(double));
    memcpy(s->vy, mv + 4 * nm, nm * sizeof(double));
    memcpy(s->vz, mv + 5 * nm, nm * sizeof(double));
    memcpy(s->gm, mv + 6 * nm, nm * sizeof(double));
    memcpy(s->radius, mv + 7 * nm, nm * sizeof(double));
    memcpy(s->j2, mv + 8 * nm, nm * sizeof(double));
    memcpy(s->j4, mv + 9 * nm, nm * sizeof(double));
    memcpy(s->rref, mv + 10 * nm, nm * sizeof(double));
    memcpy(s->polx, mv + 11 * nm, nm * sizeof(double));
    memcpy(s->poly, mv + 12 * nm, nm * sizeof(double));
    memcpy(s->polz, mv + 13 * nm, nm * sizeof(double));
    memcpy(s->kind, mi + 0 * nm, nm * sizeof(int));
    memcpy(s->parent, mi + 1 * nm, nm * sizeof(int));
    memcpy(s->alive, sn->ma, nm);
    memcpy(s->name, sn->mn, nm * STC_SPACE_NAME);
    memcpy(s->tx, sn->tv + 0 * nt, nt * sizeof(double));
    memcpy(s->ty, sn->tv + 1 * nt, nt * sizeof(double));
    memcpy(s->tz, sn->tv + 2 * nt, nt * sizeof(double));
    memcpy(s->tvx, sn->tv + 3 * nt, nt * sizeof(double));
    memcpy(s->tvy, sn->tv + 4 * nt, nt * sizeof(double));
    memcpy(s->tvz, sn->tv + 5 * nt, nt * sizeof(double));
    memcpy(s->talive, sn->ta, nt);
    s->t = sn->t;
    s->steps = sn->steps;
    s->nm = sn->nm;
    s->nt = sn->nt;
    s->nev = sn->nev;
    s->e0 = sn->e0;
    memcpy(s->ev, sn->ev, sizeof(s->ev));
    if (s->ck) {
        s->ck->avalid = 0;
    }
    return 0;
}

/* Store a snapshot at the current time. The ring holds snapshots in
 * ascending t; when full, the endpoint farthest from the current time
 * is evicted (recent coverage in whichever direction we are going). */
static void sp_ckpt_store(StcSpace *s)
{
    StcSpaceCkpt *ck = s->ck;
    int lo, hi, pos;
    if (ck == NULL || s->ckpt_max < 1 || !(s->ckpt_interval > 0.0)) {
        return;
    }
    if (ck->sncap < s->ckpt_max) {
        int ncap = s->ckpt_max < 64 ? 64 : s->ckpt_max;
        SpSnap *ns = (SpSnap *)sp_xrealloc(ck->sn, sizeof(SpSnap) * (size_t)ncap);
        int k;
        if (ns == NULL) {
            return;
        }
        for (k = ck->sncap; k < ncap; k++) {
            memset(&ns[k], 0, sizeof(SpSnap));
        }
        ck->sn = ns;
        ck->sncap = ncap;
    }
    lo = 0;
    hi = ck->count;
    while (lo < hi) {
        int mid = (lo + hi) >> 1;
        if (ck->sn[mid].t < s->t) lo = mid + 1; else hi = mid;
    }
    pos = lo;
    if (pos < ck->count && ck->sn[pos].t == s->t) {
        if (sp_snap_fit(&ck->sn[pos], s->nm, s->nt) == 0) {
            sp_snap_save(s, &ck->sn[pos]);
        }
        return;
    }
    if (ck->count < s->ckpt_max && ck->count < ck->sncap) {
        memmove(&ck->sn[pos + 1], &ck->sn[pos],
                sizeof(SpSnap) * (size_t)(ck->count - pos));
        memset(&ck->sn[pos], 0, sizeof(SpSnap));
        ck->count++;
        if (sp_snap_fit(&ck->sn[pos], s->nm, s->nt) == 0) {
            sp_snap_save(s, &ck->sn[pos]);
        }
        return;
    }
    /* full: evict the endpoint farthest from the current time */
    {
        double dl = s->t - ck->sn[0].t;
        double dr = ck->sn[ck->count - 1].t - s->t;
        int drop_first;
        if (pos == 0) {
            return; /* t is older than everything stored */
        }
        drop_first = dl >= dr;
        if (drop_first) {
            sp_snap_free(&ck->sn[0]);
            memmove(&ck->sn[0], &ck->sn[1],
                    sizeof(SpSnap) * (size_t)(ck->count - 1));
            memset(&ck->sn[ck->count - 1], 0, sizeof(SpSnap));
            pos--;
        } else {
            sp_snap_free(&ck->sn[ck->count - 1]);
            memset(&ck->sn[ck->count - 1], 0, sizeof(SpSnap));
        }
        memmove(&ck->sn[pos + 1], &ck->sn[pos],
                sizeof(SpSnap) * (size_t)(ck->count - 1 - pos));
        memset(&ck->sn[pos], 0, sizeof(SpSnap));
        if (sp_snap_fit(&ck->sn[pos], s->nm, s->nt) == 0) {
            sp_snap_save(s, &ck->sn[pos]);
        }
    }
}

/* Next checkpoint-grid time strictly in `dir` from t (grid is the
 * absolute lines k*ckpt_interval; dir > 0 forward, < 0 backward). */
static double sp_next_grid(double t, double iv, int dir)
{
    double k;
    if (dir > 0) {
        k = floor(t / iv);
        return (k + 1.0) * iv;
    }
    k = ceil(t / iv);
    return (k - 1.0) * iv;
}

static int sp_on_grid(double t, double iv)
{
    double k = round(t / iv);
    return k * iv == t;
}

static int ck_is_active(const StcSpace *s)
{
    return s->ck != NULL && s->ckpt_max >= 1 && s->ckpt_interval > 0.0;
}

/* ---------------------------------------------------------------- */
/* advance / seek                                                    */

long stc_space_advance(StcSpace *s, double dt, long max_substeps)
{
    double h, iv, n_d;
    long n, k, taken = 0;
    int dir;

    if (s == NULL || s->ck == NULL || !(s->h > 0.0) || !(dt != 0.0) ||
        !isfinite(dt)) {
        return 0;
    }
    n_d = ceil(fabs(dt) / s->h);
    if (n_d > 9.0e15) {
        n_d = 9.0e15;
    }
    n = (long)n_d;
    if (n < 1) {
        return 0;
    }
    h = dt / (double)n;
    dir = dt > 0.0 ? 1 : -1;
    iv = ck_is_active(s) ? s->ckpt_interval : 0.0;

    if (s->steps == 0) {
        s->e0 = stc_space_energy(s);
    }
    /* anchor the timeline: first advance stores its start state */
    if (ck_is_active(s) && s->ck->count == 0) {
        sp_ckpt_store(s);
    }

    for (k = 0; k < n; k++) {
        double tend = s->t + h;
        for (;;) {
            double g = (iv > 0.0) ? sp_next_grid(s->t, iv, dir) : s->t;
            int split = iv > 0.0 &&
                        (dir > 0 ? (g > s->t && g < tend)
                                 : (g < s->t && g > tend));
            double piece = (split ? g : tend) - s->t;
            sp_substep(s, piece);
            s->t += piece;
            taken++;
            sp_events(s, piece);
            if (split || (iv > 0.0 && sp_on_grid(s->t, iv))) {
                sp_ckpt_store(s);
            }
            if (max_substeps > 0 && taken >= max_substeps) {
                return -taken;
            }
            if (!split) {
                break;
            }
        }
    }
    if (taken > 0) {
        sp_ckpt_store(s); /* end-of-advance snapshot */
    }
    return taken;
}

int stc_space_seek(StcSpace *s, double t)
{
    StcSpaceCkpt *ck;
    int k, best = -1;
    if (s == NULL || s->ck == NULL || !isfinite(t)) {
        return -1;
    }
    ck = s->ck;
    for (k = 0; k < ck->count; k++) {
        if (ck->sn[k].t <= t &&
            (best < 0 || ck->sn[k].t > ck->sn[best].t)) {
            best = k;
        }
    }
    if (best < 0) {
        return -1; /* every checkpoint is later than t: outside span */
    }
    if (sp_snap_load(s, &ck->sn[best]) != 0) {
        return -1;
    }
    if (s->t == t) {
        return 0;
    }
    (void)stc_space_advance(s, t - s->t, 0);
    return 0;
}

/* ---------------------------------------------------------------- */
/* diagnostics                                                       */

void stc_space_accel(const StcSpace *s, double *ax, double *ay,
                     double *az)
{
    int i;
    if (s == NULL || ax == NULL || ay == NULL || az == NULL) {
        return;
    }
    if (s->ck == NULL) {
        for (i = 0; i < s->nm; i++) {
            ax[i] = ay[i] = az[i] = 0.0;
        }
        return;
    }
    s->ck->avalid = 0; /* force a fresh evaluation */
    sp_forces(s);
    memcpy(ax, s->ck->ax, sizeof(double) * (size_t)s->nm);
    memcpy(ay, s->ck->ay, sizeof(double) * (size_t)s->nm);
    memcpy(az, s->ck->az, sizeof(double) * (size_t)s->nm);
}

typedef struct {
    const StcSpace *s;
    double *part;
    int what; /* 0 KE, 1 PE pair partial, 2 count */
} SpSumCtx;

static void sp_sum_fn(void *vctx, int i0, int i1)
{
    const SpSumCtx *c = (const SpSumCtx *)vctx;
    const StcSpace *s = c->s;
    int i, j;
    for (i = i0; i < i1; i++) {
        double v = 0.0;
        if (s->alive[i]) {
            if (c->what == 0) {
                v = 0.5 * s->gm[i] * SP_MASS_G *
                    (s->vx[i] * s->vx[i] + s->vy[i] * s->vy[i] +
                     s->vz[i] * s->vz[i]);
            } else if (c->what == 1) {
                double xi = s->x[i], yi = s->y[i], zi = s->z[i];
                double soft2 = s->soft * s->soft;
                for (j = i + 1; j < s->nm; j++) {
                    if (s->alive[j]) {
                        double dx = s->x[j] - xi;
                        double dy = s->y[j] - yi;
                        double dz = s->z[j] - zi;
                        double r2 = dx * dx + dy * dy + dz * dz + soft2;
                        if (r2 > 0.0) {
                            v += -s->gm[i] * s->gm[j] * SP_MASS_G /
                                 sqrt(r2);
                        }
                    }
                }
            } else {
                v = 1.0;
            }
        }
        c->part[i] = v;
    }
}

/* Deterministic ordered reduction: each index writes its own slot in
 * parallel, then a fixed-order serial sum combines them. */
static double sp_reduce(const StcSpace *s, int what)
{
    StcSpaceCkpt *ck = s->ck;
    SpSumCtx c;
    double tot = 0.0;
    int i;
    if (sp_fit_scratch(s) != 0) {
        return 0.0;
    }
    c.s = s;
    c.part = ck->part;
    c.what = what;
    stc_par_for(sp_sum_fn, &c, s->nm, 0);
    for (i = 0; i < s->nm; i++) {
        tot += ck->part[i];
    }
    return tot;
}

double stc_space_energy(const StcSpace *s)
{
    if (s == NULL || s->ck == NULL) {
        return 0.0;
    }
    return sp_reduce(s, 0) + sp_reduce(s, 1);
}

typedef struct {
    const StcSpace *s;
    double *part3;
} SpMomCtx;

static void sp_mom_fn(void *vctx, int i0, int i1)
{
    const SpMomCtx *c = (const SpMomCtx *)vctx;
    const StcSpace *s = c->s;
    int i;
    for (i = i0; i < i1; i++) {
        double m = s->alive[i] ? s->gm[i] * SP_MASS_G : 0.0;
        c->part3[3 * i] = m * s->vx[i];
        c->part3[3 * i + 1] = m * s->vy[i];
        c->part3[3 * i + 2] = m * s->vz[i];
    }
}

void stc_space_momentum(const StcSpace *s, double p[3])
{
    StcSpaceCkpt *ck;
    SpMomCtx c;
    int i;
    if (p == NULL) {
        return;
    }
    p[0] = p[1] = p[2] = 0.0;
    if (s == NULL || s->ck == NULL || sp_fit_scratch(s) != 0) {
        return;
    }
    ck = s->ck;
    c.s = s;
    c.part3 = ck->part3;
    stc_par_for(sp_mom_fn, &c, s->nm, 0);
    for (i = 0; i < s->nm; i++) {
        p[0] += ck->part3[3 * i];
        p[1] += ck->part3[3 * i + 1];
        p[2] += ck->part3[3 * i + 2];
    }
}

int stc_space_live_massive(const StcSpace *s)
{
    int n = 0, i;
    if (s == NULL) {
        return 0;
    }
    for (i = 0; i < s->nm; i++) {
        n += s->alive[i] ? 1 : 0;
    }
    return n;
}

int stc_space_live_tracers(const StcSpace *s)
{
    int n = 0, i;
    if (s == NULL) {
        return 0;
    }
    for (i = 0; i < s->nt; i++) {
        n += s->talive[i] ? 1 : 0;
    }
    return n;
}

/* ---------------------------------------------------------------- */
/* osculating elements                                               */

/* Standard state-vector -> Kepler elements, mu = gm_i + gm_j. */
int stc_space_orbit(const StcSpace *s, int i, int j, StcOrbit *o)
{
    double rx, ry, rz, vx, vy, vz;
    double r, v2, mu, hx, hy, hz, h2, h;
    double nx, ny, nn, ex, ey, ez, e;
    double xi, a, inc, node, argp, nu, M, period, v;

    if (s == NULL || o == NULL || i < 0 || i >= s->nm || !s->alive[i]) {
        return -1;
    }
    if (j < 0) {
        int p = (i < s->nm) ? s->parent[i] : -1;
        if (p >= 0 && p < s->nm && p != i && s->alive[p]) {
            j = p;
        } else {
            int k, bm = -1;
            double bg = 0.0;
            for (k = 0; k < s->nm; k++) {
                if (k != i && s->alive[k] && s->gm[k] > bg) {
                    bg = s->gm[k];
                    bm = k;
                }
            }
            j = bm;
        }
    }
    if (j < 0 || j >= s->nm || j == i || !s->alive[j]) {
        return -1;
    }
    rx = s->x[i] - s->x[j];
    ry = s->y[i] - s->y[j];
    rz = s->z[i] - s->z[j];
    vx = s->vx[i] - s->vx[j];
    vy = s->vy[i] - s->vy[j];
    vz = s->vz[i] - s->vz[j];
    r = sqrt(rx * rx + ry * ry + rz * rz);
    v2 = vx * vx + vy * vy + vz * vz;
    v = sqrt(v2);
    mu = s->gm[i] + s->gm[j];
    if (!(r > 0.0) || !(mu > 0.0)) {
        return -1;
    }
    hx = ry * vz - rz * vy;
    hy = rz * vx - rx * vz;
    hz = rx * vy - ry * vx;
    h2 = hx * hx + hy * hy + hz * hz;
    h = sqrt(h2);
    if (!(h > 0.0)) {
        return -1; /* radial trajectory: elements undefined */
    }
    /* eccentricity vector e_vec = (v x h)/mu - r/r */
    {
        double cx_ = vy * hz - vz * hy;
        double cy_ = vz * hx - vx * hz;
        double cz_ = vx * hy - vy * hx;
        ex = cx_ / mu - rx / r;
        ey = cy_ / mu - ry / r;
        ez = cz_ / mu - rz / r;
    }
    e = sqrt(ex * ex + ey * ey + ez * ez);
    nx = -hy; /* z-hat x h */
    ny = hx;
    nn = sqrt(nx * nx + ny * ny);
    xi = v2 * 0.5 - mu / r;
    a = -mu / (2.0 * xi);
    inc = acos(hz / h < -1.0 ? -1.0 : (hz / h > 1.0 ? 1.0 : hz / h));
    if (nn > 0.0) {
        node = atan2(ny, nx);
        if (node < 0.0) {
            node += 2.0 * STC_PI;
        }
        /* argp: angle from node line to e_vec */
        argp = atan2(((ny * ez - 0.0 * ey) * hx +
                      (0.0 * ex - nx * ez) * hy +
                      (nx * ey - ny * ex) * hz) / (nn * e * h),
                     (nx * ex + ny * ey) / (nn * e));
        if (!(e > 1e-12)) {
            argp = 0.0;
        }
    } else {
        node = 0.0;
        argp = atan2(ey, ex);
        if (!(e > 1e-12)) {
            argp = 0.0;
        }
    }
    /* true anomaly: angle from e_vec to r, signed by h . (e x r) */
    if (e > 1e-12) {
        double sgn = (ey * rz - ez * ry) * hx +
                     (ez * rx - ex * rz) * hy +
                     (ex * ry - ey * rx) * hz;
        double cnu = (ex * rx + ey * ry + ez * rz) / (e * r);
        if (cnu > 1.0) cnu = 1.0;
        if (cnu < -1.0) cnu = -1.0;
        nu = acos(cnu);
        if (sgn < 0.0) {
            nu = -nu;
        }
    } else {
        /* circular: measure from node line (or x axis if equatorial) */
        double ax_ = nn > 0.0 ? nx : 1.0;
        double ay_ = nn > 0.0 ? ny : 0.0;
        double az_ = 0.0;
        double an = sqrt(ax_ * ax_ + ay_ * ay_ + az_ * az_);
        double cnu = (ax_ * rx + ay_ * ry + az_ * rz) / (an * r);
        double sgn = (ay_ * rz - az_ * ry) * hx +
                     (az_ * rx - ax_ * rz) * hy +
                     (ax_ * ry - ay_ * rx) * hz;
        if (cnu > 1.0) cnu = 1.0;
        if (cnu < -1.0) cnu = -1.0;
        nu = acos(cnu);
        if (sgn < 0.0) {
            nu = -nu;
        }
    }
    if (e < 1.0) {
        double E = atan2(sqrt(1.0 - e * e) * sin(nu), e + cos(nu));
        M = E - e * sin(E);
    } else {
        double den = 1.0 + e * cos(nu);
        double sh = sqrt(e * e - 1.0) * sin(nu) /
                    (den > 1e-30 ? den : 1e-30);
        double F = asinh(sh);
        M = e * sinh(F) - F;
    }
    period = (a > 0.0 && isfinite(a)) ? 2.0 * STC_PI * sqrt(a * a * a / mu)
                                    : 0.0;
    o->a = a;
    o->e = e;
    o->inc = inc;
    o->node = node;
    o->argp = argp;
    o->M = M;
    o->period = period;
    o->r = r;
    o->v = v;
    return 0;
}
