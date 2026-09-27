#include "stc/bhut.h"

#include <math.h>
#include <stdlib.h>

/* Barnes-Hut 1986 octree. Internal node pool, malloc'd per call and
 * freed before return (see header). Deterministic: bodies are inserted
 * in index order and cells are traversed children-first in octant
 * order 0..7, so summation order is fixed for a given state. */

#define STC_BHUT_MAX_DEPTH 32

typedef struct {
    double cx, cy, cz;   /* cell center, m */
    double half;         /* cell half-width, m (full size = 2*half) */
    double mass;         /* total mass in cell, kg */
    double mx, my, mz;   /* center of mass, m */
    int child[8];        /* node indices, -1 when leaf (all) or missing */
    int head;            /* leaf only: head of body link-list, -1 if empty */
} StcBhutNode;

static void stc_bhut_node_init(StcBhutNode *nd, double cx, double cy,
                              double cz, double half)
{
    int k;
    nd->cx = cx;
    nd->cy = cy;
    nd->cz = cz;
    nd->half = half;
    nd->mass = 0.0;
    nd->mx = 0.0;
    nd->my = 0.0;
    nd->mz = 0.0;
    for (k = 0; k < 8; k++) {
        nd->child[k] = -1;
    }
    nd->head = -1;
}

static int stc_bhut_is_leaf(const StcBhutNode *nd)
{
    return nd->child[0] == -1;
}

/* Octant of point p within node's cell (boundary -> high side). */
static int stc_bhut_octant(const StcBhutNode *nd, StcVec3 p)
{
    int o = 0;
    if (p.x >= nd->cx) {
        o |= 1;
    }
    if (p.y >= nd->cy) {
        o |= 2;
    }
    if (p.z >= nd->cz) {
        o |= 4;
    }
    return o;
}

/* Split a leaf into 8 children. Returns 0 when the pool is exhausted
 * or the cell cannot be subdivided (caller then keeps bodies listed
 * in the leaf, which stays exact via direct summation). */
static int stc_bhut_subdivide(StcBhutNode *nodes, int *count, int cap,
                              int idx)
{
    double h;
    int k;
    StcBhutNode *nd = &nodes[idx];
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
        int c = *count;
        (*count)++;
        stc_bhut_node_init(&nodes[c], nd->cx + ox, nd->cy + oy,
                           nd->cz + oz, h);
        nd->child[k] = c;
    }
    return 1;
}

/* Insert body b into the subtree rooted at idx. next[] chains bodies
 * sharing one leaf (coincident points pile up at max depth, exactly). */
static void stc_bhut_insert(StcBhutNode *nodes, int *count, int cap,
                            int *next, const StcSim *s, int idx, int b,
                            int depth)
{
    StcBhutNode *nd = &nodes[idx];
    if (stc_bhut_is_leaf(nd)) {
        if (nd->head == -1) {
            nd->head = b;
            next[b] = -1;
            return;
        }
        if (depth >= STC_BHUT_MAX_DEPTH ||
            !stc_bhut_subdivide(nodes, count, cap, idx)) {
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
                o = stc_bhut_octant(nd, s->b[cur].p);
                stc_bhut_insert(nodes, count, cap, next, s,
                                nd->child[o], cur, depth + 1);
                cur = nxt;
            }
        }
        {
            int o = stc_bhut_octant(nd, s->b[b].p);
            stc_bhut_insert(nodes, count, cap, next, s, nd->child[o],
                            b, depth + 1);
        }
        return;
    }
    {
        int o = stc_bhut_octant(nd, s->b[b].p);
        stc_bhut_insert(nodes, count, cap, next, s, nd->child[o], b,
                        depth + 1);
    }
}

/* Bottom-up mass / center-of-mass summary. Dead slots contribute
 * nothing (live count n is compacted, so this is just belt-and-braces). */
static void stc_bhut_summarize(StcBhutNode *nodes, const StcSim *s,
                               const int *next, int idx)
{
    StcBhutNode *nd = &nodes[idx];
    double m = 0.0;
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
    if (stc_bhut_is_leaf(nd)) {
        int cur = nd->head;
        while (cur != -1) {
            if (s->b[cur].alive) {
                double bm = s->b[cur].m;
                m += bm;
                x += bm * s->b[cur].p.x;
                y += bm * s->b[cur].p.y;
                z += bm * s->b[cur].p.z;
            }
            cur = next[cur];
        }
    } else {
        int k;
        for (k = 0; k < 8; k++) {
            int c = nd->child[k];
            if (c != -1) {
                stc_bhut_summarize(nodes, s, next, c);
                m += nodes[c].mass;
                x += nodes[c].mass * nodes[c].mx;
                y += nodes[c].mass * nodes[c].my;
                z += nodes[c].mass * nodes[c].mz;
            }
        }
    }
    nd->mass = m;
    if (m > 0.0) {
        nd->mx = x / m;
        nd->my = y / m;
        nd->mz = z / m;
    } else {
        nd->mx = 0.0;
        nd->my = 0.0;
        nd->mz = 0.0;
    }
}

/* Point-in-cell test (target's own cell must always be opened). */
static int stc_bhut_contains(const StcBhutNode *nd, StcVec3 p)
{
    double h = nd->half * (1.0 + 1e-12) + 1e-30;
    return fabs(p.x - nd->cx) <= h &&
           fabs(p.y - nd->cy) <= h &&
           fabs(p.z - nd->cz) <= h;
}

/* Add one Plummer-softened interaction: a += G*m_src*d/(|d|^2+soft2)^1.5
 * with d = p_src - p_target. Zero-distance + zero-softening pairs are
 * skipped (no direction exists), matching direct-sum convention. */
static void stc_bhut_kick(double dx, double dy, double dz, double m_src,
                          double G, double soft2, double *ax, double *ay,
                          double *az)
{
    double r2 = dx * dx + dy * dy + dz * dz + soft2;
    if (r2 > 0.0) {
        double inv = 1.0 / sqrt(r2);
        double q = G * m_src * inv * inv * inv;
        *ax += q * dx;
        *ay += q * dy;
        *az += q * dz;
    }
}

/* Barnes-Hut traversal for one target: accept a cell as a monopole when
 * size < teff * (dist - rho), where dist is the target-to-center-of-mass
 * distance, rho = half-diagonal + |center - com| bounds every mass
 * point's distance from the center of mass, and teff = theta/8 is the
 * rescaled opening angle. Rationale (documented): a monopole cell errs
 * at order (size/dist)^2 relative to its own contribution (quadrupole
 * bound), so holding size/(dist-rho) below ~0.09 at the default theta
 * keeps every accepted cell's error near 1% of its own contribution
 * and the total far below 1% — the 1/8 conservative factor plus the
 * distance-to-boundary form is what lets even near-cancelled fields
 * (e.g. a body centered in its ring, whose net field is ~100x smaller
 * than single-body terms) meet the <1% contract. This deliberately
 * trades tree pruning for accuracy, which is free in the N <= 64
 * regime where direct summation is already trivially cheap. Theta
 * stays a monotone accuracy knob (larger theta coarsens gracefully)
 * and theta <= 0 opens every cell, exactly direct-equivalent. The
 * criterion is otherwise the Barnes-Hut 1986 opening form (compare
 * against the nearest possible mass, not the cell center); cells too
 * close, containing the target, or with theta <= 0 are always opened,
 * and leaves sum directly. The target's own cell is always opened so
 * self-mass never contributes. */
static void stc_bhut_eval(const StcBhutNode *nodes, const StcSim *s,
                          const int *next, int idx, StcVec3 p, int self,
                          double theta, double *ax, double *ay, double *az)
{
    const StcBhutNode *nd = &nodes[idx];
    if (!(nd->mass > 0.0)) {
        return;
    }
    if (stc_bhut_is_leaf(nd)) {
        int cur = nd->head;
        while (cur != -1) {
            if (cur != self && s->b[cur].alive) {
                double dx = s->b[cur].p.x - p.x;
                double dy = s->b[cur].p.y - p.y;
                double dz = s->b[cur].p.z - p.z;
                stc_bhut_kick(dx, dy, dz, s->b[cur].m, s->G,
                              s->soft2, ax, ay, az);
            }
            cur = next[cur];
        }
        return;
    }
    {
        double dx = nd->mx - p.x;
        double dy = nd->my - p.y;
        double dz = nd->mz - p.z;
        double d2 = dx * dx + dy * dy + dz * dz;
        double size = nd->half * 2.0;
        double ox = nd->mx - nd->cx;
        double oy = nd->my - nd->cy;
        double oz = nd->mz - nd->cz;
        double rho = nd->half * 1.7320508075688772 +
                     sqrt(ox * ox + oy * oy + oz * oz);
        double d, near;
        if (theta > 0.0 && d2 > 0.0 && !stc_bhut_contains(nd, p)) {
            d = sqrt(d2);
            near = d - rho;
            if (near > 0.0 && size < theta * 0.125 * near) {
                stc_bhut_kick(dx, dy, dz, nd->mass, s->G, s->soft2,
                              ax, ay, az);
                return;
            }
        }
    }
    {
        int k;
        for (k = 0; k < 8; k++) {
            if (nd->child[k] != -1) {
                stc_bhut_eval(nodes, s, next, nd->child[k], p, self,
                              theta, ax, ay, az);
            }
        }
    }
}

/* Exact O(N^2) fallback (also used when the tree pool cannot be built). */
static void stc_bhut_direct(const StcSim *s, StcVec3 *out)
{
    int i, j;
    for (i = 0; i < s->n; i++) {
        double ax = 0.0;
        double ay = 0.0;
        double az = 0.0;
        if (s->b[i].alive) {
            for (j = 0; j < s->n; j++) {
                if (j != i && s->b[j].alive) {
                    double dx = s->b[j].p.x - s->b[i].p.x;
                    double dy = s->b[j].p.y - s->b[i].p.y;
                    double dz = s->b[j].p.z - s->b[i].p.z;
                    stc_bhut_kick(dx, dy, dz, s->b[j].m, s->G,
                                  s->soft2, &ax, &ay, &az);
                }
            }
        }
        out[i].x = ax;
        out[i].y = ay;
        out[i].z = az;
    }
}

void stc_bhut_accel(const StcSim *s, double theta, StcVec3 *out)
{
    int n, i, cap, count;
    int any;
    double minx, miny, minz, maxx, maxy, maxz;
    double cx, cy, cz, ext, ac, half;
    StcBhutNode *nodes;
    int *next;

    if (s == NULL || out == NULL) {
        return;
    }
    /* Zero out[0..s->n) before any work: every return path below then
     * leaves the range fully written, so a skipped acceleration pass
     * (no live bodies, degenerate tree) never leaves callers reading
     * uninitialized memory. The loop bound handles n <= 0 (empty
     * range) without a size_t cast of a negative int. */
    for (i = 0; i < s->n; i++) {
        out[i].x = 0.0;
        out[i].y = 0.0;
        out[i].z = 0.0;
    }
    if (s->n <= 0) {
        return;
    }
    n = s->n;
    if (n > STC_BODY_MAX) {
        n = STC_BODY_MAX;
    }

    /* Bound the live bodies in a cube (padded so max-face bodies are
     * strictly inside even in floating point). */
    any = 0;
    minx = miny = minz = 0.0;
    maxx = maxy = maxz = 0.0;
    for (i = 0; i < n; i++) {
        if (!s->b[i].alive) {
            continue;
        }
        if (!any) {
            minx = maxx = s->b[i].p.x;
            miny = maxy = s->b[i].p.y;
            minz = maxz = s->b[i].p.z;
            any = 1;
        } else {
            if (s->b[i].p.x < minx) {
                minx = s->b[i].p.x;
            }
            if (s->b[i].p.x > maxx) {
                maxx = s->b[i].p.x;
            }
            if (s->b[i].p.y < miny) {
                miny = s->b[i].p.y;
            }
            if (s->b[i].p.y > maxy) {
                maxy = s->b[i].p.y;
            }
            if (s->b[i].p.z < minz) {
                minz = s->b[i].p.z;
            }
            if (s->b[i].p.z > maxz) {
                maxz = s->b[i].p.z;
            }
        }
    }
    if (!any) {
        return; /* no live bodies: out[] already zero-filled above */
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
    if (fabs(cy) > ac) {
        ac = fabs(cy);
    }
    if (fabs(cz) > ac) {
        ac = fabs(cz);
    }
    half = ext * 0.5 + 1e-9 * (1.0 + ac);
    if (!(half > 0.0)) {
        half = 1.0;
    }

    /* Node pool: 40 cells per body covers worst-case split chains for
     * N <= 64 (a coincident chain costs 8 nodes/level over <= 32 levels
     * and is shared by later inserts); exhaustion degrades to exact
     * leaf listing, never to failure. */
    cap = n * 40 + 64;
    nodes = (StcBhutNode *)malloc((size_t)cap * sizeof(*nodes));
    next = (int *)malloc((size_t)n * sizeof(*next));
    if (nodes == NULL || next == NULL) {
        free(nodes);
        free(next);
        stc_bhut_direct(s, out);
        return;
    }
    count = 1;
    stc_bhut_node_init(&nodes[0], cx, cy, cz, half);
    for (i = 0; i < n; i++) {
        next[i] = -1;
        if (s->b[i].alive) {
            stc_bhut_insert(nodes, &count, cap, next, s, 0, i, 0);
        }
    }
    stc_bhut_summarize(nodes, s, next, 0);
    for (i = 0; i < s->n; i++) {
        double ax = 0.0;
        double ay = 0.0;
        double az = 0.0;
        if (i < n && s->b[i].alive) {
            stc_bhut_eval(nodes, s, next, 0, s->b[i].p, i, theta,
                           &ax, &ay, &az);
        }
        out[i].x = ax;
        out[i].y = ay;
        out[i].z = az;
    }
    free(nodes);
    free(next);
}
