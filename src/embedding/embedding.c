/* SPACETIME LAB — embedding diagrams: Flamm paraboloid meshes.
 * Implements stc/embedding.h. SI doubles (meters), no malloc, fixed loops.
 */

#include "stc/embedding.h"

#include <math.h>
#include <stddef.h>
#include <stdint.h>

#include "stc/units.h"
#include "stc/vec.h"

double stc_embedding_flamm(double r_m, double rs_m)
{
    if (!isfinite(r_m) || !isfinite(rs_m)) {
        return NAN;
    }
    if (rs_m < 0.0) {
        return NAN;
    }
    if (rs_m == 0.0) {
        return 0.0;
    }
    if (r_m <= rs_m) {
        return 0.0; /* clamped at the throat */
    }
    return 2.0 * sqrt(rs_m * (r_m - rs_m));
}

double stc_embedding_flamm_slope(double r_m, double rs_m)
{
    if (!isfinite(r_m) || !isfinite(rs_m)) {
        return NAN;
    }
    if (rs_m < 0.0) {
        return NAN;
    }
    if (rs_m == 0.0) {
        return 0.0;
    }
    if (r_m <= rs_m) {
        return INFINITY; /* vertical throat */
    }
    return sqrt(rs_m / (r_m - rs_m));
}

double stc_embedding_height(double x_m, double y_m,
                            const StcEmbeddingWell *wells, size_t n_wells)
{
    double h;
    size_t k;

    if (!isfinite(x_m) || !isfinite(y_m)) {
        return NAN;
    }
    if (n_wells == 0) {
        return 0.0;
    }
    if (wells == NULL) {
        return NAN;
    }
    h = 0.0;
    for (k = 0; k < n_wells; k++) {
        double dx, dy, r;
        double cx = wells[k].x;
        double cy = wells[k].y;
        double rs = wells[k].rs;

        if (!isfinite(cx) || !isfinite(cy) || !isfinite(rs)) {
            return NAN;
        }
        if (rs < 0.0) {
            return NAN;
        }
        if (rs == 0.0) {
            continue;
        }
        dx = x_m - cx;
        dy = y_m - cy;
        r = hypot(dx, dy);
        h -= stc_embedding_flamm(r, rs);
    }
    return h;
}

int stc_embedding_height_grad(double x_m, double y_m,
                              const StcEmbeddingWell *wells, size_t n_wells,
                              double *h_out, double *gx_out, double *gy_out)
{
    double h, gx, gy;
    size_t k;

    if (h_out == NULL || gx_out == NULL || gy_out == NULL) {
        return STC_EMBEDDING_ERR_NULL;
    }
    if (!isfinite(x_m) || !isfinite(y_m)) {
        *h_out = NAN;
        *gx_out = NAN;
        *gy_out = NAN;
        return STC_EMBEDDING_ERR_RANGE;
    }
    if (n_wells == 0) {
        *h_out = 0.0;
        *gx_out = 0.0;
        *gy_out = 0.0;
        return STC_EMBEDDING_OK;
    }
    if (wells == NULL) {
        *h_out = NAN;
        *gx_out = NAN;
        *gy_out = NAN;
        return STC_EMBEDDING_ERR_NULL;
    }
    h = 0.0;
    gx = 0.0;
    gy = 0.0;
    for (k = 0; k < n_wells; k++) {
        double dx, dy, r, slope;
        double cx = wells[k].x;
        double cy = wells[k].y;
        double rs = wells[k].rs;

        if (!isfinite(cx) || !isfinite(cy) || !isfinite(rs) || rs < 0.0) {
            *h_out = NAN;
            *gx_out = NAN;
            *gy_out = NAN;
            return STC_EMBEDDING_ERR_RANGE;
        }
        if (rs == 0.0) {
            continue;
        }
        dx = x_m - cx;
        dy = y_m - cy;
        r = hypot(dx, dy);
        h -= stc_embedding_flamm(r, rs);
        if (r > rs && r > 0.0) {
            slope = stc_embedding_flamm_slope(r, rs);
            gx -= slope * dx / r;
            gy -= slope * dy / r;
        }
        /* r <= rs: clamped plateau, zero gradient contribution. */
    }
    *h_out = h;
    *gx_out = gx;
    *gy_out = gy;
    return STC_EMBEDDING_OK;
}

size_t stc_embedding_radial_vert_count(size_t nr, size_t nphi)
{
    if (nr == 0 || nphi == 0) {
        return 0;
    }
    if (nr > SIZE_MAX / nphi) {
        return 0;
    }
    return nr * nphi;
}

size_t stc_embedding_radial_index_count(size_t nr, size_t nphi)
{
    size_t quads;

    if (nr < 2 || nphi == 0) {
        return 0;
    }
    if ((nr - 1) > SIZE_MAX / nphi) {
        return 0;
    }
    quads = (nr - 1) * nphi;
    if (quads > SIZE_MAX / 6) {
        return 0;
    }
    return quads * 6;
}

size_t stc_embedding_cartesian_vert_count(size_t nx, size_t ny)
{
    if (nx == 0 || ny == 0) {
        return 0;
    }
    if (nx > SIZE_MAX / ny) {
        return 0;
    }
    return nx * ny;
}

size_t stc_embedding_cartesian_index_count(size_t nx, size_t ny)
{
    size_t quads;

    if (nx < 2 || ny < 2) {
        return 0;
    }
    if ((nx - 1) > SIZE_MAX / (ny - 1)) {
        return 0;
    }
    quads = (nx - 1) * (ny - 1);
    if (quads > SIZE_MAX / 6) {
        return 0;
    }
    return quads * 6;
}

/* Fill one radial-grid vertex position. r/phi/w in meters/radians/meters. */
static void stc_embedding_radial_point(double r, double phi, double w,
                                       double *x, double *y, double *z)
{
    *x = r * cos(phi);
    *y = r * sin(phi);
    *z = w;
}

int stc_embedding_radial_grid(double rs_m, double r_min_m, double r_max_m,
                              size_t nr, size_t nphi,
                              StcVec3 *verts, size_t n_verts,
                              uint32_t *indices, size_t n_indices)
{
    size_t need_v, need_i, i, j, k;

    if (nr < 2 || nphi < 3) {
        return STC_EMBEDDING_ERR_DIMS;
    }
    if (!isfinite(rs_m) || !isfinite(r_min_m) || !isfinite(r_max_m)) {
        return STC_EMBEDDING_ERR_RANGE;
    }
    if (rs_m < 0.0 || r_min_m < 0.0 || r_max_m <= r_min_m) {
        return STC_EMBEDDING_ERR_RANGE;
    }
    need_v = stc_embedding_radial_vert_count(nr, nphi);
    need_i = stc_embedding_radial_index_count(nr, nphi);
    if (need_v == 0 || need_i == 0 || need_v > (size_t)UINT32_MAX) {
        return STC_EMBEDDING_ERR_DIMS;
    }
    if (verts == NULL || indices == NULL) {
        return STC_EMBEDDING_ERR_NULL;
    }
    if (n_verts < need_v || n_indices < need_i) {
        return STC_EMBEDDING_ERR_SMALL;
    }
    for (i = 0; i < nr; i++) {
        double t = (double)i / (double)(nr - 1);
        double r = r_min_m + (r_max_m - r_min_m) * t;
        double w = stc_embedding_flamm(r, rs_m);

        for (j = 0; j < nphi; j++) {
            double phi = 2.0 * STC_PI * (double)j / (double)nphi;
            double x, y, z;

            stc_embedding_radial_point(r, phi, w, &x, &y, &z);
            verts[i * nphi + j].x = x;
            verts[i * nphi + j].y = y;
            verts[i * nphi + j].z = z;
        }
    }
    k = 0;
    for (i = 0; i < nr - 1; i++) {
        for (j = 0; j < nphi; j++) {
            uint32_t a = (uint32_t)(i * nphi + j);
            uint32_t b = (uint32_t)(i * nphi + (j + 1) % nphi);
            uint32_t c = (uint32_t)((i + 1) * nphi + j);
            uint32_t d = (uint32_t)((i + 1) * nphi + (j + 1) % nphi);

            indices[k++] = a;
            indices[k++] = c;
            indices[k++] = b;
            indices[k++] = b;
            indices[k++] = c;
            indices[k++] = d;
        }
    }
    return STC_EMBEDDING_OK;
}

int stc_embedding_radial_grid_xyz(double rs_m, double r_min_m, double r_max_m,
                                  size_t nr, size_t nphi,
                                  double *xyz, size_t n_xyz,
                                  uint32_t *indices, size_t n_indices)
{
    size_t need_v, need_i, need_xyz, i, j, k;

    if (nr < 2 || nphi < 3) {
        return STC_EMBEDDING_ERR_DIMS;
    }
    if (!isfinite(rs_m) || !isfinite(r_min_m) || !isfinite(r_max_m)) {
        return STC_EMBEDDING_ERR_RANGE;
    }
    if (rs_m < 0.0 || r_min_m < 0.0 || r_max_m <= r_min_m) {
        return STC_EMBEDDING_ERR_RANGE;
    }
    need_v = stc_embedding_radial_vert_count(nr, nphi);
    need_i = stc_embedding_radial_index_count(nr, nphi);
    if (need_v == 0 || need_i == 0 || need_v > (size_t)UINT32_MAX) {
        return STC_EMBEDDING_ERR_DIMS;
    }
    if (need_v > SIZE_MAX / 3) {
        return STC_EMBEDDING_ERR_DIMS;
    }
    need_xyz = need_v * 3;
    if (xyz == NULL || indices == NULL) {
        return STC_EMBEDDING_ERR_NULL;
    }
    if (n_xyz < need_xyz || n_indices < need_i) {
        return STC_EMBEDDING_ERR_SMALL;
    }
    for (i = 0; i < nr; i++) {
        double t = (double)i / (double)(nr - 1);
        double r = r_min_m + (r_max_m - r_min_m) * t;
        double w = stc_embedding_flamm(r, rs_m);

        for (j = 0; j < nphi; j++) {
            double phi = 2.0 * STC_PI * (double)j / (double)nphi;
            double x, y, z;

            stc_embedding_radial_point(r, phi, w, &x, &y, &z);
            xyz[(i * nphi + j) * 3 + 0] = x;
            xyz[(i * nphi + j) * 3 + 1] = y;
            xyz[(i * nphi + j) * 3 + 2] = z;
        }
    }
    k = 0;
    for (i = 0; i < nr - 1; i++) {
        for (j = 0; j < nphi; j++) {
            uint32_t a = (uint32_t)(i * nphi + j);
            uint32_t b = (uint32_t)(i * nphi + (j + 1) % nphi);
            uint32_t c = (uint32_t)((i + 1) * nphi + j);
            uint32_t d = (uint32_t)((i + 1) * nphi + (j + 1) % nphi);

            indices[k++] = a;
            indices[k++] = c;
            indices[k++] = b;
            indices[k++] = b;
            indices[k++] = c;
            indices[k++] = d;
        }
    }
    return STC_EMBEDDING_OK;
}

/* Validate shared Cartesian-grid inputs; 0 if valid, else a status code. */
static int stc_embedding_cartesian_check(const StcEmbeddingWell *wells, size_t n_wells,
                                        double x0_m, double x1_m, size_t nx,
                                        double y0_m, double y1_m, size_t ny)
{
    size_t k;

    if (nx < 2 || ny < 2) {
        return STC_EMBEDDING_ERR_DIMS;
    }
    if (!isfinite(x0_m) || !isfinite(x1_m) || !isfinite(y0_m) || !isfinite(y1_m)) {
        return STC_EMBEDDING_ERR_RANGE;
    }
    if (x1_m <= x0_m || y1_m <= y0_m) {
        return STC_EMBEDDING_ERR_RANGE;
    }
    if (n_wells > 0 && wells == NULL) {
        return STC_EMBEDDING_ERR_NULL;
    }
    for (k = 0; k < n_wells; k++) {
        if (!isfinite(wells[k].x) || !isfinite(wells[k].y) ||
            !isfinite(wells[k].rs) || wells[k].rs < 0.0) {
            return STC_EMBEDDING_ERR_RANGE;
        }
    }
    return STC_EMBEDDING_OK;
}

int stc_embedding_cartesian_grid(const StcEmbeddingWell *wells, size_t n_wells,
                                 double x0_m, double x1_m, size_t nx,
                                 double y0_m, double y1_m, size_t ny,
                                 StcVec3 *verts, size_t n_verts,
                                 uint32_t *indices, size_t n_indices)
{
    int check;
    size_t need_v, need_i, i, j, k;

    check = stc_embedding_cartesian_check(wells, n_wells, x0_m, x1_m, nx,
                                          y0_m, y1_m, ny);
    if (check != STC_EMBEDDING_OK) {
        return check;
    }
    need_v = stc_embedding_cartesian_vert_count(nx, ny);
    need_i = stc_embedding_cartesian_index_count(nx, ny);
    if (need_v == 0 || need_i == 0 || need_v > (size_t)UINT32_MAX) {
        return STC_EMBEDDING_ERR_DIMS;
    }
    if (verts == NULL || indices == NULL) {
        return STC_EMBEDDING_ERR_NULL;
    }
    if (n_verts < need_v || n_indices < need_i) {
        return STC_EMBEDDING_ERR_SMALL;
    }
    for (j = 0; j < ny; j++) {
        double ty = (double)j / (double)(ny - 1);
        double y = y0_m + (y1_m - y0_m) * ty;

        for (i = 0; i < nx; i++) {
            double tx = (double)i / (double)(nx - 1);
            double x = x0_m + (x1_m - x0_m) * tx;

            verts[j * nx + i].x = x;
            verts[j * nx + i].y = y;
            verts[j * nx + i].z = stc_embedding_height(x, y, wells, n_wells);
        }
    }
    k = 0;
    for (j = 0; j < ny - 1; j++) {
        for (i = 0; i < nx - 1; i++) {
            uint32_t a = (uint32_t)(j * nx + i);
            uint32_t b = (uint32_t)(j * nx + i + 1);
            uint32_t c = (uint32_t)((j + 1) * nx + i);
            uint32_t d = (uint32_t)((j + 1) * nx + i + 1);

            indices[k++] = a;
            indices[k++] = b;
            indices[k++] = c;
            indices[k++] = b;
            indices[k++] = d;
            indices[k++] = c;
        }
    }
    return STC_EMBEDDING_OK;
}

int stc_embedding_cartesian_grid_xyz(const StcEmbeddingWell *wells, size_t n_wells,
                                     double x0_m, double x1_m, size_t nx,
                                     double y0_m, double y1_m, size_t ny,
                                     double *xyz, size_t n_xyz,
                                     uint32_t *indices, size_t n_indices)
{
    int check;
    size_t need_v, need_i, need_xyz, i, j, k;

    check = stc_embedding_cartesian_check(wells, n_wells, x0_m, x1_m, nx,
                                          y0_m, y1_m, ny);
    if (check != STC_EMBEDDING_OK) {
        return check;
    }
    need_v = stc_embedding_cartesian_vert_count(nx, ny);
    need_i = stc_embedding_cartesian_index_count(nx, ny);
    if (need_v == 0 || need_i == 0 || need_v > (size_t)UINT32_MAX) {
        return STC_EMBEDDING_ERR_DIMS;
    }
    if (need_v > SIZE_MAX / 3) {
        return STC_EMBEDDING_ERR_DIMS;
    }
    need_xyz = need_v * 3;
    if (xyz == NULL || indices == NULL) {
        return STC_EMBEDDING_ERR_NULL;
    }
    if (n_xyz < need_xyz || n_indices < need_i) {
        return STC_EMBEDDING_ERR_SMALL;
    }
    for (j = 0; j < ny; j++) {
        double ty = (double)j / (double)(ny - 1);
        double y = y0_m + (y1_m - y0_m) * ty;

        for (i = 0; i < nx; i++) {
            double tx = (double)i / (double)(nx - 1);
            double x = x0_m + (x1_m - x0_m) * tx;

            xyz[(j * nx + i) * 3 + 0] = x;
            xyz[(j * nx + i) * 3 + 1] = y;
            xyz[(j * nx + i) * 3 + 2] = stc_embedding_height(x, y, wells, n_wells);
        }
    }
    k = 0;
    for (j = 0; j < ny - 1; j++) {
        for (i = 0; i < nx - 1; i++) {
            uint32_t a = (uint32_t)(j * nx + i);
            uint32_t b = (uint32_t)(j * nx + i + 1);
            uint32_t c = (uint32_t)((j + 1) * nx + i);
            uint32_t d = (uint32_t)((j + 1) * nx + i + 1);

            indices[k++] = a;
            indices[k++] = b;
            indices[k++] = c;
            indices[k++] = b;
            indices[k++] = d;
            indices[k++] = c;
        }
    }
    return STC_EMBEDDING_OK;
}
