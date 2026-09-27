#include "stc/vec.h"

#include <math.h>

StcVec3 stc_vec3_add(StcVec3 a, StcVec3 b)
{
    StcVec3 r;
    r.x = a.x + b.x;
    r.y = a.y + b.y;
    r.z = a.z + b.z;
    return r;
}

StcVec3 stc_vec3_sub(StcVec3 a, StcVec3 b)
{
    StcVec3 r;
    r.x = a.x - b.x;
    r.y = a.y - b.y;
    r.z = a.z - b.z;
    return r;
}

StcVec3 stc_vec3_scale(StcVec3 v, double s)
{
    StcVec3 r;
    r.x = v.x * s;
    r.y = v.y * s;
    r.z = v.z * s;
    return r;
}

double stc_vec3_dot(StcVec3 a, StcVec3 b)
{
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

StcVec3 stc_vec3_cross(StcVec3 a, StcVec3 b)
{
    StcVec3 r;
    r.x = a.y * b.z - a.z * b.y;
    r.y = a.z * b.x - a.x * b.z;
    r.z = a.x * b.y - a.y * b.x;
    return r;
}

double stc_vec3_len(StcVec3 v)
{
    return sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
}

StcVec3 stc_vec3_norm(StcVec3 v)
{
    StcVec3 r = { 0.0, 0.0, 0.0 };
    double len = stc_vec3_len(v);
    if (len > 0.0) {
        r.x = v.x / len;
        r.y = v.y / len;
        r.z = v.z / len;
    }
    return r;
}

double stc_vec3_dist(StcVec3 a, StcVec3 b)
{
    double dx = a.x - b.x;
    double dy = a.y - b.y;
    double dz = a.z - b.z;
    return sqrt(dx * dx + dy * dy + dz * dz);
}

StcMat4 stc_mat4_identity(void)
{
    StcMat4 r;
    int i;
    for (i = 0; i < 16; i++) {
        r.m[i] = 0.0;
    }
    r.m[0] = 1.0;
    r.m[5] = 1.0;
    r.m[10] = 1.0;
    r.m[15] = 1.0;
    return r;
}

StcMat4 stc_mat4_mul(StcMat4 a, StcMat4 b)
{
    /* Column-major, column vectors: r = a * b (b applied first).
     * r(col,row) = sum_k a(k,row) * b(col,k). */
    StcMat4 r;
    int row, col, k;
    for (col = 0; col < 4; col++) {
        for (row = 0; row < 4; row++) {
            double s = 0.0;
            for (k = 0; k < 4; k++) {
                s += a.m[k * 4 + row] * b.m[col * 4 + k];
            }
            r.m[col * 4 + row] = s;
        }
    }
    return r;
}

StcMat4 stc_mat4_perspective(double fovy_rad, double aspect,
                             double znear, double zfar)
{
    /* Right-handed, NDC z in [-1, 1] (OpenGL convention), column-major.
     * Degenerate inputs -> identity (deterministic, no NaN). */
    StcMat4 r;
    int i;
    double f, nf;
    for (i = 0; i < 16; i++) {
        r.m[i] = 0.0;
    }
    if (!(fovy_rad > 0.0) || !(aspect > 0.0) ||
        !(znear > 0.0) || !(zfar > znear)) {
        return stc_mat4_identity();
    }
    f = 1.0 / tan(fovy_rad * 0.5);
    nf = 1.0 / (znear - zfar);
    r.m[0] = f / aspect;
    r.m[5] = f;
    r.m[10] = (zfar + znear) * nf;
    r.m[11] = -1.0;
    r.m[14] = 2.0 * zfar * znear * nf;
    return r;
}

StcMat4 stc_mat4_lookat(StcVec3 eye, StcVec3 center, StcVec3 up)
{
    /* Right-handed look-at, column-major. Column vectors.
     * Degenerate (eye==center or up parallel) -> identity. */
    StcVec3 f, s, u;
    StcMat4 r;
    int i;
    f = stc_vec3_norm(stc_vec3_sub(center, eye));
    if (stc_vec3_len(f) == 0.0) {
        return stc_mat4_identity();
    }
    s = stc_vec3_norm(stc_vec3_cross(f, up));
    if (stc_vec3_len(s) == 0.0) {
        return stc_mat4_identity();
    }
    u = stc_vec3_cross(s, f);
    for (i = 0; i < 16; i++) {
        r.m[i] = 0.0;
    }
    r.m[0] = s.x;
    r.m[1] = u.x;
    r.m[2] = -f.x;
    r.m[4] = s.y;
    r.m[5] = u.y;
    r.m[6] = -f.y;
    r.m[8] = s.z;
    r.m[9] = u.z;
    r.m[10] = -f.z;
    r.m[12] = -stc_vec3_dot(s, eye);
    r.m[13] = -stc_vec3_dot(u, eye);
    r.m[14] = stc_vec3_dot(f, eye);
    r.m[15] = 1.0;
    return r;
}

StcVec3 stc_mat4_transform_point(StcMat4 m, StcVec3 p)
{
    /* Homogeneous w = 1, then perspective divide. w == 0 -> input. */
    double x, y, z, w;
    StcVec3 r;
    x = m.m[0] * p.x + m.m[4] * p.y + m.m[8] * p.z + m.m[12];
    y = m.m[1] * p.x + m.m[5] * p.y + m.m[9] * p.z + m.m[13];
    z = m.m[2] * p.x + m.m[6] * p.y + m.m[10] * p.z + m.m[14];
    w = m.m[3] * p.x + m.m[7] * p.y + m.m[11] * p.z + m.m[15];
    if (w == 0.0) {
        return p;
    }
    r.x = x / w;
    r.y = y / w;
    r.z = z / w;
    return r;
}

double stc_clamp(double x, double lo, double hi)
{
    if (x < lo) {
        return lo;
    }
    if (x > hi) {
        return hi;
    }
    return x;
}

double stc_lerp(double a, double b, double t)
{
    return a + (b - a) * t;
}
