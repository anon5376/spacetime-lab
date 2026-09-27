#ifndef STC_PAR_H
#define STC_PAR_H

/* SPACETIME LAB — par: tiny persistent worker pool (POSIX threads).
 *
 * stc_par_for(fn, ctx, n) runs fn(ctx, i0, i1) over disjoint chunks that
 * cover [0, n) exactly once, using every online CPU (the caller thread
 * participates), and returns when all chunks are done. Chunks are handed
 * out dynamically, so callers MUST make each index's work independent of
 * which thread runs it and of chunk order (write only to slots owned by
 * the index range; no shared accumulators). Under that rule results are
 * bitwise identical for any thread count (determinism contract R4).
 *
 * Not reentrant: do not call stc_par_for from inside a running fn.
 * STC_PAR_THREADS=1 in the environment forces serial execution.
 */

typedef void (*StcParFn)(void *ctx, int i0, int i1);

/* Run fn over [0, n) in chunks of `chunk` indices (chunk < 1 -> auto:
 * >= 256 indices per chunk, so small loops run serially on the caller).
 * n <= 0 is a no-op. */
void stc_par_for(StcParFn fn, void *ctx, int n, int chunk);

/* Number of threads the pool uses (>= 1), including the caller. */
int stc_par_threads(void);

#endif /* STC_PAR_H */
