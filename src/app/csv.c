/* SPACETIME LAB — app: deterministic CSV writer implementation.
 * Pure C11, stdio only. Fixed %.17g format => byte-identical output.
 */

#include "stc/csv.h"

int stc_csv_open(StcCsv *csv, const char *path, const char *header)
{
    if (csv == NULL || path == NULL) {
        return -1;
    }
    if (csv->open) {
        return -1;
    }
    csv->fp = NULL;
    csv->ncols = -1;
    csv->nrows = 0;
    csv->open = 0;

    csv->fp = fopen(path, "w");
    if (csv->fp == NULL) {
        return -2;
    }
    if (header != NULL) {
        if (fprintf(csv->fp, "%s\n", header) < 0) {
            fclose(csv->fp);
            csv->fp = NULL;
            return -2;
        }
    }
    csv->open = 1;
    return 0;
}

int stc_csv_row(StcCsv *csv, const double *vals, int n)
{
    int i;

    if (csv == NULL || vals == NULL || !csv->open || csv->fp == NULL) {
        return -1;
    }
    if (n <= 0) {
        return -2;
    }
    if (csv->ncols < 0) {
        csv->ncols = n;
    } else if (csv->ncols != n) {
        return -3;
    }
    for (i = 0; i < n; i++) {
        if (i > 0) {
            if (fputc(',', csv->fp) == EOF) {
                return -4;
            }
        }
        if (fprintf(csv->fp, "%.17g", vals[i]) < 0) {
            return -4;
        }
    }
    if (fputc('\n', csv->fp) == EOF) {
        return -4;
    }
    csv->nrows++;
    return 0;
}

void stc_csv_close(StcCsv *csv)
{
    if (csv == NULL || !csv->open) {
        return;
    }
    if (csv->fp != NULL) {
        fclose(csv->fp);
        csv->fp = NULL;
    }
    csv->open = 0;
}

int stc_csv_is_open(const StcCsv *csv)
{
    if (csv == NULL || !csv->open || csv->fp == NULL) {
        return 0;
    }
    return 1;
}
