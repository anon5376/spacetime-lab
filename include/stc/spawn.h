#ifndef STC_SPAWN_H
#define STC_SPAWN_H

/* SPACETIME LAB — spawn: shared body-spawn form (terminal + visual UIs).
 * Pure C11, no I/O except sbx_form_draw (caller-printed terminal form).
 * SI doubles via strtod validation; deterministic (no RNG, no clock).
 * Both sandbox frontends share one field model, one validator, one key
 * machine so typing behaves identically everywhere.
 */

#include "stc/body.h"

/* Eleven fields: name, mass_kg, radius_m, px/py/pz_m, vx/vy/vz_m_s,
 * fixed, SPAWN> (commit line, index SBX_NFIELDS - 1). */
#define SBX_NFIELDS 11

typedef struct {
    char buf[SBX_NFIELDS][64];
    int field;
    char err[128];
} SbxForm;

/* Method: static string table — field display name.
 * Returns the name for field i, or "?" when out of range. */
const char *sbx_field_name(int i);

/* Prefill: probe on a circular orbit about body 0 through the view
 * center (cx, cy, cz, m) with view spread (m). Clears all fields. */
void sbx_form_prefill(SbxForm *f, const StcSim *s, double cx, double cy,
                      double cz, double spread);

/* Strict number parse (strtod, full-string, finite). Returns 0 ok, -1 bad. */
int sbx_parse_num(const char *txt, double *out);

/* Validate + spawn. Returns 0 spawned, -1 with f->err set. */
int sbx_form_spawn(StcSim *s, SbxForm *f);

/* Terminal form painter (stdout, ANSI clear). The visual UI reads
 * f->buf/f->field/f->err directly and paints its own overlay. */
void sbx_form_draw(const SbxForm *f);

/* Feed one key to the form (stc/ui.h codes; printable ASCII appends,
 * 127/8 backspace, up/down/enter navigate, SPAWN>+enter commits).
 * Returns 1 spawned+done, 0 still editing, -1 cancelled. */
int sbx_form_key(StcSim *s, SbxForm *f, int key);

#endif /* STC_SPAWN_H */
