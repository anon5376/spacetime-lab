#ifndef STC_CSV_H
#define STC_CSV_H

/* SPACETIME LAB — app: deterministic CSV writer.
 * Pure C11, SI doubles. One optional header line, then comma-joined
 * %.17g rows with '\n' endings. Same doubles in = byte-identical file.
 */

#include <stdio.h>

/* Open CSV handle: live FILE plus row/column counters. */
typedef struct {
    FILE *fp;
    int ncols;  /* columns per row, set by first row; -1 until then */
    long nrows; /* data rows written */
    int open;   /* nonzero while fp is live */
} StcCsv;

/* Equation: CSV open — fopen(path, "w") + optional header line.
 * csv and path are required; header may be NULL (no header line).
 * Returns 0 on success, -1 on NULL arg or already-open handle,
 * -2 if fopen or the header write fails (handle left closed). */
int stc_csv_open(StcCsv *csv, const char *path, const char *header);

/* Equation: CSV row — comma-joined %.17g doubles + '\n'.
 * The first row sets ncols; later rows must match it. vals holds n doubles.
 * Returns 0 on success, -1 if csv/vals NULL or handle closed, -2 if n <= 0,
 * -3 on column-count mismatch, -4 on write error. */
int stc_csv_row(StcCsv *csv, const double *vals, int n);

/* Equation: CSV close — fflush + fclose + state reset. Idempotent, NULL-safe. */
void stc_csv_close(StcCsv *csv);

/* Equation: CSV state query — nonzero iff open with a live fp. NULL -> 0. */
int stc_csv_is_open(const StcCsv *csv);

#endif /* STC_CSV_H */
