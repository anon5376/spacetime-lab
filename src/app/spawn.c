/* SPACETIME LAB — spawn: shared body-spawn form (terminal + visual UIs).
 * One field model, one validator, one key machine for both frontends.
 */

#include "stc/spawn.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "stc/elements.h"
#include "stc/ui.h"

static const char *const kFieldNames[SBX_NFIELDS] = {
    "name", "mass_kg", "radius_m", "px_m", "py_m", "pz_m",
    "vx_m_s", "vy_m_s", "vz_m_s", "fixed", "SPAWN>"
};

/* Method: static string table — field display name (see header). */
const char *sbx_field_name(int i)
{
    if (i < 0 || i >= SBX_NFIELDS) {
        return "?";
    }
    return kFieldNames[i];
}

/* Prefill: probe on a circular orbit about body 0 through the view center. */
void sbx_form_prefill(SbxForm *f, const StcSim *s, double cx, double cy,
                      double cz, double spread)
{
    double px = cx + spread * 0.5;
    double py = cy;
    double pz = cz;
    double vx = 0.0;
    double vy = 0.0;
    double vz = 0.0;
    int i;

    for (i = 0; i < SBX_NFIELDS; i++) {
        f->buf[i][0] = '\0';
    }
    f->field = 0;
    f->err[0] = '\0';
    if (s != NULL && s->n > 0 && s->b[0].m > 0.0) {
        double rx = px - s->b[0].p.x;
        double ry = py - s->b[0].p.y;
        double rz = pz - s->b[0].p.z;
        double r = sqrt(rx * rx + ry * ry + rz * rz);

        if (r > 0.0) {
            double vc = stc_circular_velocity(s->G * s->b[0].m, r);
            double ux = -rz / r;
            double uy = 0.0;
            double uz = rx / r;
            double ul = sqrt(ux * ux + uz * uz);

            if (ul > 1e-9) {
                vx = vc * ux / ul;
                vz = vc * uz / ul;
            } else {
                vx = vc;
            }
            (void)uy;
        }
    }
    snprintf(f->buf[0], sizeof f->buf[0], "probe");
    snprintf(f->buf[1], sizeof f->buf[1], "1000");
    snprintf(f->buf[2], sizeof f->buf[2], "10");
    snprintf(f->buf[3], sizeof f->buf[3], "%.6g", px);
    snprintf(f->buf[4], sizeof f->buf[4], "%.6g", py);
    snprintf(f->buf[5], sizeof f->buf[5], "%.6g", pz);
    snprintf(f->buf[6], sizeof f->buf[6], "%.6g", vx);
    snprintf(f->buf[7], sizeof f->buf[7], "%.6g", vy);
    snprintf(f->buf[8], sizeof f->buf[8], "%.6g", vz);
    snprintf(f->buf[9], sizeof f->buf[9], "n");
}

int sbx_parse_num(const char *txt, double *out)
{
    char *end = NULL;
    double v;

    if (txt == NULL || *txt == '\0' || out == NULL) {
        return -1;
    }
    v = strtod(txt, &end);
    if (end == txt || *end != '\0' || !(v > -1e308 && v < 1e308)) {
        return -1;
    }
    *out = v;
    return 0;
}

/* Validate + spawn. Returns 0 spawned, -1 with f->err set. */
int sbx_form_spawn(StcSim *s, SbxForm *f)
{
    char name[32];
    double m;
    double r;
    double v[6];
    StcVec3 p;
    StcVec3 vv;
    int fixed;
    size_t i;

    if (s == NULL || f == NULL) {
        return -1;
    }
    if (sbx_parse_num(f->buf[1], &m) != 0 || !(m > 0.0)) {
        snprintf(f->err, sizeof f->err, "mass must be a positive number");
        return -1;
    }
    if (sbx_parse_num(f->buf[2], &r) != 0 || !(r >= 0.0)) {
        snprintf(f->err, sizeof f->err, "radius must be >= 0");
        return -1;
    }
    for (i = 0; i < 6; i++) {
        if (sbx_parse_num(f->buf[3 + i], &v[i]) != 0) {
            snprintf(f->err, sizeof f->err, "bad number in %s",
                     kFieldNames[3 + i]);
            return -1;
        }
    }
    if (strcmp(f->buf[9], "y") == 0 || strcmp(f->buf[9], "Y") == 0 ||
        strcmp(f->buf[9], "1") == 0) {
        fixed = 1;
    } else if (strcmp(f->buf[9], "n") == 0 ||
               strcmp(f->buf[9], "N") == 0 ||
               strcmp(f->buf[9], "0") == 0 || f->buf[9][0] == '\0') {
        fixed = 0;
    } else {
        snprintf(f->err, sizeof f->err, "fixed must be y/n");
        return -1;
    }
    snprintf(name, sizeof name, "%s", f->buf[0][0] != '\0' ? f->buf[0]
                                                            : "body");
    for (i = 0; name[i] != '\0'; i++) {
        if (name[i] == ',' || name[i] == '"' || name[i] == '\n') {
            name[i] = ' ';
        }
    }
    p.x = v[0];
    p.y = v[1];
    p.z = v[2];
    vv.x = v[3];
    vv.y = v[4];
    vv.z = v[5];
    if (stc_sim_add(s, name, m, r, p, vv, fixed) < 0) {
        snprintf(f->err, sizeof f->err, "sim full (%d bodies max)",
                 STC_BODY_MAX);
        return -1;
    }
    return 0;
}

void sbx_form_draw(const SbxForm *f)
{
    int i;

    printf("\x1b[2J\x1b[H");
    printf("SPAWN BODY  (up/down field, type to edit, enter next, esc cancel)\n\n");
    for (i = 0; i < SBX_NFIELDS; i++) {
        if (i == SBX_NFIELDS - 1) {
            printf("%s [%s]   <- enter here to spawn\n",
                   f->field == i ? ">" : " ", kFieldNames[i]);
        } else {
            printf("%s %-8s : %s%s\n", f->field == i ? ">" : " ",
                   kFieldNames[i], f->buf[i],
                   f->field == i ? "_" : "");
        }
    }
    if (f->err[0] != '\0') {
        printf("\n! %s\n", f->err);
    }
    printf("\nnumbers accept 1.989e30 form; fixed is y/n\n");
    fflush(stdout);
}

/* Feed one key to the form. Returns 1 spawned+done, 0 still editing,
 * -1 cancelled. */
int sbx_form_key(StcSim *s, SbxForm *f, int key)
{
    size_t len;

    if (key == STC_UI_KEY_ESC) {
        return -1;
    }
    if (key == STC_UI_KEY_UP) {
        f->field = (f->field + SBX_NFIELDS - 1) % SBX_NFIELDS;
        return 0;
    }
    if (key == STC_UI_KEY_DOWN || key == '\r' || key == '\n') {
        if (key != STC_UI_KEY_DOWN && f->field == SBX_NFIELDS - 1) {
            return sbx_form_spawn(s, f) == 0 ? 1 : 0;
        }
        f->field = (f->field + 1) % SBX_NFIELDS;
        return 0;
    }
    if (f->field >= SBX_NFIELDS - 1) {
        return 0;
    }
    len = strlen(f->buf[f->field]);
    if (key == 127 || key == 8) {
        if (len > 0) {
            f->buf[f->field][len - 1] = '\0';
        }
        return 0;
    }
    if (key >= 32 && key < 127 && len + 1 < sizeof f->buf[0]) {
        f->buf[f->field][len] = (char)key;
        f->buf[f->field][len + 1] = '\0';
    }
    return 0;
}
