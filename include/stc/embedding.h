#ifndef STC_EMBEDDING_H
#define STC_EMBEDDING_H

/* SPACETIME LAB — embedding diagrams: Flamm paraboloid meshes.
 *
 * Pure C11, SI doubles (meters everywhere). Caller-owned buffers: this
 * module never calls malloc/free (no allocation in the mesh hot path).
 * Deterministic: fixed loop counts, no RNG, no I/O, no globals.
 *
 * Single-well surface (exact GR equatorial-slice embedding):
 *   Flamm paraboloid  w = 2*sqrt(rs*(r-rs)),  r > rs,
 * clamped to w = 0 at/inside the throat (r <= rs).
 *
 * Multi-well surfaces are an ILLUSTRATIVE superposition for visualization
 * only (there is no static multi-black-hole slice of this form in GR).
 */

#include <stddef.h>
#include <stdint.h>

#include "stc/vec.h"

/* Status codes returned by the mesh/gradient entry points. */
#define STC_EMBEDDING_OK 0 /* success */
#define STC_EMBEDDING_ERR_NULL -1 /* NULL pointer (buffer/out-param/wells) */
#define STC_EMBEDDING_ERR_DIMS -2 /* bad grid dims, size_t overflow, or vert count > UINT32_MAX */
#define STC_EMBEDDING_ERR_RANGE -3 /* non-finite/negative physical inputs or inverted ranges */
#define STC_EMBEDDING_ERR_SMALL -4 /* caller buffer smaller than the required count */

/* One illustrative well: center (x, y) in meters, Schwarzschild radius in meters. */
typedef struct {
    double x;
    double y;
    double rs;
} StcEmbeddingWell;

/* Equation: Flamm paraboloid w = 2*sqrt(rs*(r-rs)).
 * r_m, rs_m in meters; returns w in meters.
 * r <= rs clamps to 0 (throat); rs == 0 is flat (0); negative rs or
 * non-finite input returns NaN. */
double stc_embedding_flamm(double r_m, double rs_m);

/* Equation: Flamm slope dw/dr = sqrt(rs/(r-rs)).
 * Dimensionless. r <= rs returns +INFINITY (vertical throat);
 * rs == 0 returns 0; negative rs or non-finite input returns NaN. */
double stc_embedding_flamm_slope(double r_m, double rs_m);

/* Equation: illustrative multi-well superposition H(x,y) = -sum_i w_i(r_i),
 * with r_i the planar distance to well i and w_i its Flamm height.
 * NOT a GR solution; visualization aid only. x, y, well fields in meters;
 * returns H in meters (wells depress toward -w, far field -> 0).
 * n_wells == 0 returns 0; NULL wells with n_wells > 0, negative/non-finite
 * well data, or non-finite x/y returns NaN. Wells with rs == 0 are flat. */
double stc_embedding_height(double x_m, double y_m,
                            const StcEmbeddingWell *wells, size_t n_wells);

/* Equation: gradient of the illustrative multi-well superposition above,
 * dH/dx = -sum_i (dw_i/dr_i)(dx_i/r_i) (and likewise for y).
 * Particle-height query: h_out gets H, gx_out/gy_out its planar gradient
 * (dimensionless) for seating particles/velocities on the sheet.
 * Contributions with r_i <= rs_i sit on the clamped plateau (zero gradient).
 * Returns STC_EMBEDDING_OK; NULL out-params or NULL wells (n_wells > 0)
 * return STC_EMBEDDING_ERR_NULL; bad well data or non-finite x/y sets all
 * outputs to NaN and returns STC_EMBEDDING_ERR_RANGE. */
int stc_embedding_height_grad(double x_m, double y_m,
                              const StcEmbeddingWell *wells, size_t n_wells,
                              double *h_out, double *gx_out, double *gy_out);

/* Equation: polar-grid sampling counts for the Flamm paraboloid.
 * Vertices nr*nphi; indices (nr-1)*nphi*6. Return 0 for degenerate dims
 * or size_t overflow (treat 0 as "do not allocate / invalid"). */
size_t stc_embedding_radial_vert_count(size_t nr, size_t nphi);
size_t stc_embedding_radial_index_count(size_t nr, size_t nphi);

/* Equation: Cartesian-grid sampling counts for the multi-well sheet.
 * Vertices nx*ny; indices (nx-1)*(ny-1)*6. Return 0 for degenerate dims
 * or size_t overflow. */
size_t stc_embedding_cartesian_vert_count(size_t nx, size_t ny);
size_t stc_embedding_cartesian_index_count(size_t nx, size_t ny);

/* Equation: polar sampling of the Flamm paraboloid.
 * Ring i samples r_i = r_min + (r_max-r_min)*i/(nr-1), segment j samples
 * phi_j = 2*pi*j/nphi; vert = (r_i*cos phi_j, r_i*sin phi_j, w(r_i)).
 * All radii in meters. Triangles are CCW seen from +w. Requires nr >= 2,
 * nphi >= 3, r_max > r_min >= 0, finite rs >= 0; vert count must fit
 * uint32_t. n_verts/n_indices are buffer capacities; undersized buffers
 * return STC_EMBEDDING_ERR_SMALL and leave buffers untouched. */
int stc_embedding_radial_grid(double rs_m, double r_min_m, double r_max_m,
                              size_t nr, size_t nphi,
                              StcVec3 *verts, size_t n_verts,
                              uint32_t *indices, size_t n_indices);

/* Equation: polar sampling of the Flamm paraboloid (flat-array variant).
 * Same surface/winding as stc_embedding_radial_grid; positions land in
 * xyz[3*k+0/1/2] (plain doubles for the renderer). n_xyz counts doubles
 * and must hold 3*vert_count. */
int stc_embedding_radial_grid_xyz(double rs_m, double r_min_m, double r_max_m,
                                  size_t nr, size_t nphi,
                                  double *xyz, size_t n_xyz,
                                  uint32_t *indices, size_t n_indices);

/* Equation: Cartesian sampling of the illustrative multi-well superposition.
 * Grid point (i,j) sits at (x0+(x1-x0)*i/(nx-1), y0+(y1-y0)*j/(ny-1), H)
 * with H from stc_embedding_height. All coordinates in meters. Triangles
 * are CCW seen from +z. Requires nx >= 2, ny >= 2, x1 > x0, y1 > y0
 * (finite), valid well data (finite, rs >= 0), vert count fitting uint32_t.
 * Undersized buffers return STC_EMBEDDING_ERR_SMALL, buffers untouched. */
int stc_embedding_cartesian_grid(const StcEmbeddingWell *wells, size_t n_wells,
                                 double x0_m, double x1_m, size_t nx,
                                 double y0_m, double y1_m, size_t ny,
                                 StcVec3 *verts, size_t n_verts,
                                 uint32_t *indices, size_t n_indices);

/* Equation: Cartesian sampling of the illustrative multi-well superposition
 * (flat-array variant). Same sheet/winding as stc_embedding_cartesian_grid;
 * positions land in xyz[3*k+0/1/2] (plain doubles for the renderer). n_xyz
 * counts doubles and must hold 3*vert_count. */
int stc_embedding_cartesian_grid_xyz(const StcEmbeddingWell *wells, size_t n_wells,
                                     double x0_m, double x1_m, size_t nx,
                                     double y0_m, double y1_m, size_t ny,
                                     double *xyz, size_t n_xyz,
                                     uint32_t *indices, size_t n_indices);

#endif /* STC_EMBEDDING_H */
