#include "stc/trails.h"

#include <stdlib.h>

int stc_trail_bank_init(StcTrailBank *bank, size_t nbodies, size_t depth)
{
    /* Ring-buffer bank init: one calloc'd point block per body so a
     * mid-loop allocation failure can roll back without leaks. */
    size_t i, j;

    if (bank == NULL) {
        return -1;
    }
    bank->trails = NULL;
    bank->nbodies = 0;
    bank->depth = 0;
    if (nbodies == 0) {
        return -1;
    }
    if (depth == 0) {
        depth = (size_t)STC_TRAIL_DEPTH_DEFAULT;
    }
    if (nbodies > (size_t)-1 / depth ||
        nbodies * depth > (size_t)-1 / sizeof(StcVec3)) {
        return -1;
    }
    bank->trails = (StcTrail *)calloc(nbodies, sizeof(StcTrail));
    if (bank->trails == NULL) {
        return -1;
    }
    for (i = 0; i < nbodies; i++) {
        bank->trails[i].pts = NULL;
        bank->trails[i].depth = depth;
        bank->trails[i].head = 0;
        bank->trails[i].count = 0;
        bank->trails[i].pts = (StcVec3 *)calloc(depth, sizeof(StcVec3));
        if (bank->trails[i].pts == NULL) {
            for (j = 0; j < i; j++) {
                free(bank->trails[j].pts);
                bank->trails[j].pts = NULL;
            }
            free(bank->trails);
            bank->trails = NULL;
            return -1;
        }
    }
    bank->nbodies = nbodies;
    bank->depth = depth;
    return 0;
}

void stc_trail_bank_free(StcTrailBank *bank)
{
    /* Ring-buffer bank free: release every per-body block, then the
     * trail array, and reset to empty so double-free is a no-op. */
    size_t i;

    if (bank == NULL || bank->trails == NULL) {
        if (bank != NULL) {
            bank->trails = NULL;
            bank->nbodies = 0;
            bank->depth = 0;
        }
        return;
    }
    for (i = 0; i < bank->nbodies; i++) {
        free(bank->trails[i].pts);
        bank->trails[i].pts = NULL;
        bank->trails[i].depth = 0;
        bank->trails[i].head = 0;
        bank->trails[i].count = 0;
    }
    free(bank->trails);
    bank->trails = NULL;
    bank->nbodies = 0;
    bank->depth = 0;
}

void stc_trail_push(StcTrailBank *bank, size_t index, double x, double y,
                    double z)
{
    /* Ring-buffer push: write at head, advance head modulo depth,
     * saturate count at depth (oldest overwritten once full). */
    StcTrail *tr;

    if (bank == NULL || bank->trails == NULL || index >= bank->nbodies) {
        return;
    }
    tr = &bank->trails[index];
    if (tr->pts == NULL || tr->depth == 0) {
        return;
    }
    tr->pts[tr->head].x = x;
    tr->pts[tr->head].y = y;
    tr->pts[tr->head].z = z;
    tr->head = (tr->head + 1) % tr->depth;
    if (tr->count < tr->depth) {
        tr->count++;
    }
}

void stc_trail_clear(StcTrailBank *bank)
{
    /* Ring-buffer clear: reset every head/count pair, keep buffers. */
    size_t i;

    if (bank == NULL || bank->trails == NULL) {
        return;
    }
    for (i = 0; i < bank->nbodies; i++) {
        bank->trails[i].head = 0;
        bank->trails[i].count = 0;
    }
}

size_t stc_trail_count(const StcTrailBank *bank, size_t index)
{
    /* Ring-buffer count: stored points for one body, 0 when unknown. */
    if (bank == NULL || bank->trails == NULL || index >= bank->nbodies) {
        return 0;
    }
    return bank->trails[index].count;
}

int stc_trail_get(const StcTrailBank *bank, size_t index, size_t pos,
                  StcVec3 *out)
{
    /* Ring-buffer indexed read: oldest slot is (head - count) mod depth,
     * so pos-th oldest is (oldest + pos) mod depth. */
    const StcTrail *tr;
    size_t oldest, slot;

    if (bank == NULL || bank->trails == NULL || out == NULL ||
        index >= bank->nbodies) {
        return -1;
    }
    tr = &bank->trails[index];
    if (tr->pts == NULL || tr->depth == 0 || pos >= tr->count) {
        return -1;
    }
    oldest = (tr->head + tr->depth - tr->count) % tr->depth;
    slot = (oldest + pos) % tr->depth;
    *out = tr->pts[slot];
    return 0;
}

void stc_trail_iterate(const StcTrailBank *bank, size_t index,
                       void (*cb)(size_t pos, double x, double y, double z,
                                  void *ctx),
                       void *ctx)
{
    /* Ring-buffer oldest-to-newest traversal: same oldest-slot origin
     * as stc_trail_get, one callback per stored point in order. */
    const StcTrail *tr;
    size_t oldest, pos, slot;

    if (bank == NULL || bank->trails == NULL || cb == NULL ||
        index >= bank->nbodies) {
        return;
    }
    tr = &bank->trails[index];
    if (tr->pts == NULL || tr->depth == 0 || tr->count == 0) {
        return;
    }
    oldest = (tr->head + tr->depth - tr->count) % tr->depth;
    for (pos = 0; pos < tr->count; pos++) {
        slot = (oldest + pos) % tr->depth;
        cb(pos, tr->pts[slot].x, tr->pts[slot].y, tr->pts[slot].z, ctx);
    }
}
