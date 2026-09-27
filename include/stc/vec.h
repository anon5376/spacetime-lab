#ifndef STC_VEC_H
#define STC_VEC_H

/* SPACETIME LAB — core math: Vec3 + Mat4 (column-major) + scalar helpers.
 * Pure C11, SI doubles. Column vectors: v' = M * v.
 */

typedef struct {
    double x;
    double y;
    double z;
} StcVec3;

typedef struct {
    /* Column-major: element (row, col) is m[col * 4 + row]. */
    double m[16];
} StcMat4;

/* --- Vec3 --- */
StcVec3 stc_vec3_add(StcVec3 a, StcVec3 b);
StcVec3 stc_vec3_sub(StcVec3 a, StcVec3 b);
StcVec3 stc_vec3_scale(StcVec3 v, double s);
double stc_vec3_dot(StcVec3 a, StcVec3 b);
StcVec3 stc_vec3_cross(StcVec3 a, StcVec3 b);
double stc_vec3_len(StcVec3 v);
StcVec3 stc_vec3_norm(StcVec3 v); /* zero vector -> zero vector */
double stc_vec3_dist(StcVec3 a, StcVec3 b);

/* --- Mat4 (column-major, column vectors, result = a * b) --- */
StcMat4 stc_mat4_identity(void);
StcMat4 stc_mat4_mul(StcMat4 a, StcMat4 b);
StcMat4 stc_mat4_perspective(double fovy_rad, double aspect,
                             double znear, double zfar);
StcMat4 stc_mat4_lookat(StcVec3 eye, StcVec3 center, StcVec3 up);
StcVec3 stc_mat4_transform_point(StcMat4 m, StcVec3 p);

/* --- Scalar helpers --- */
double stc_clamp(double x, double lo, double hi);
double stc_lerp(double a, double b, double t);

#endif /* STC_VEC_H */
