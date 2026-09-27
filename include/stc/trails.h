#ifndef STC_TRAILS_H
#define STC_TRAILS_H

/* SPACETIME LAB — trails: per-body trajectory ring buffers (pure data).
 * Method: fixed-capacity ring (circular) buffer per body; each push
 * overwrites the oldest point once full. SI doubles (meters),
 * deterministic (iteration order is insertion order; no RNG, no clock).
 * No render calls and no sim mutation: the sandbox painters read stored
 * points via stc_trail_iterate / stc_trail_get + stc_trail_count.
 * Memory: stc_trail_bank_init allocates (one block per body); every
 * successful init must be paired with stc_trail_bank_free. On allocation
 * failure init frees any partial allocation, resets the bank to empty,
 * and returns -1 (no leaks). Depends only on stc/core (vec) + C stdlib.
 */

#include <stddef.h>

#include "stc/vec.h"

/* Default ring depth (points stored per body) when init gets depth == 0. */
#define STC_TRAIL_DEPTH_DEFAULT 128

/* One body's trail: depth slots written oldest-to-newest through a ring. */
typedef struct {
    StcVec3 *pts;   /* depth position slots, m (NULL when depth == 0) */
    size_t depth;   /* ring capacity in points */
    size_t head;    /* next write slot, always < depth when depth > 0 */
    size_t count;   /* stored points, 0..depth */
} StcTrail;

/* Bank of per-body trails, indexed by body slot. */
typedef struct {
    StcTrail *trails; /* nbodies ring buffers (NULL when nbodies == 0) */
    size_t nbodies;   /* number of per-body trails */
    size_t depth;     /* ring capacity per body */
} StcTrailBank;

/* Method: ring-buffer bank init — allocate nbodies empty trails of depth
 * points each (depth == 0 selects STC_TRAIL_DEPTH_DEFAULT). The bank must
 * be uninitialized memory or a freed bank (re-init without a free leaks
 * the old buffers). Returns 0 on success, -1 on NULL bank, nbodies == 0,
 * size overflow, or allocation failure (bank reset to empty, no leaks). */
int stc_trail_bank_init(StcTrailBank *bank, size_t nbodies, size_t depth);

/* Method: ring-buffer bank free — release all trail storage and reset the
 * bank to empty (nbodies/depth zeroed, trails NULL). NULL-safe; freeing
 * an empty or already-freed bank is a no-op. Pair with every init. */
void stc_trail_bank_free(StcTrailBank *bank);

/* Method: ring-buffer push — append point (x,y,z) in m to body index's
 * trail, overwriting the oldest point when full. NULL-safe; out-of-range
 * index (or an empty bank) is a no-op. */
void stc_trail_push(StcTrailBank *bank, size_t index, double x, double y,
                    double z);

/* Method: ring-buffer clear — drop all stored points in every trail
 * (counts reset, allocations kept). NULL-safe. */
void stc_trail_clear(StcTrailBank *bank);

/* Method: ring-buffer count — number of stored points for body index
 * (0..depth). Returns 0 for NULL bank or out-of-range index. */
size_t stc_trail_count(const StcTrailBank *bank, size_t index);

/* Method: ring-buffer indexed read — fetch the pos-th oldest point of
 * body index (pos 0 = oldest) into *out (m). Returns 0 on success, -1 on
 * NULL bank/out, out-of-range index, or pos >= count (*out untouched). */
int stc_trail_get(const StcTrailBank *bank, size_t index, size_t pos,
                  StcVec3 *out);

/* Method: ring-buffer oldest-to-newest traversal — call cb(pos, x, y, z,
 * ctx) in m for each stored point of body index, pos 0 = oldest.
 * NULL-safe (NULL bank/cb or out-of-range index: no-op). */
void stc_trail_iterate(const StcTrailBank *bank, size_t index,
                       void (*cb)(size_t pos, double x, double y, double z,
                                  void *ctx),
                       void *ctx);

#endif /* STC_TRAILS_H */
