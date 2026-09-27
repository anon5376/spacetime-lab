/* SPACETIME LAB — par: persistent pthread worker pool (see stc/par.h). */

#include "stc/par.h"

#include <pthread.h>
#include <stdlib.h>
#include <unistd.h>

#define PAR_MAXT 32
#define PAR_MIN_CHUNK 256

typedef struct {
    pthread_t th[PAR_MAXT];
    int nt; /* worker threads (excluding caller) */
    int ready;
    pthread_mutex_t mu;
    pthread_cond_t go;
    pthread_cond_t done;
    unsigned long gen;
    int pending;
    StcParFn fn;
    void *ctx;
    int n;
    int chunk;
    int next;
} ParPool;

static ParPool g_par;
static pthread_once_t g_par_once = PTHREAD_ONCE_INIT;

static int par_take(ParPool *p)
{
    int r;

    pthread_mutex_lock(&p->mu);
    r = p->next;
    p->next += p->chunk;
    pthread_mutex_unlock(&p->mu);
    return r;
}

static void par_drain(ParPool *p, StcParFn fn, void *ctx, int n)
{
    for (;;) {
        int i0 = par_take(p);
        int i1;

        if (i0 >= n) {
            return;
        }
        i1 = i0 + p->chunk < n ? i0 + p->chunk : n;
        fn(ctx, i0, i1);
    }
}

static void *par_worker(void *arg)
{
    ParPool *p = (ParPool *)arg;
    unsigned long seen = 0;

    for (;;) {
        StcParFn fn;
        void *ctx;
        int n;

        pthread_mutex_lock(&p->mu);
        while (p->gen == seen) {
            pthread_cond_wait(&p->go, &p->mu);
        }
        seen = p->gen;
        fn = p->fn;
        ctx = p->ctx;
        n = p->n;
        pthread_mutex_unlock(&p->mu);
        par_drain(p, fn, ctx, n);
        pthread_mutex_lock(&p->mu);
        p->pending--;
        if (p->pending == 0) {
            pthread_cond_signal(&p->done);
        }
        pthread_mutex_unlock(&p->mu);
    }
    return NULL;
}

static void par_init(void)
{
    long nc = sysconf(_SC_NPROCESSORS_ONLN);
    const char *env = getenv("STC_PAR_THREADS");
    int i;

    pthread_mutex_init(&g_par.mu, NULL);
    pthread_cond_init(&g_par.go, NULL);
    pthread_cond_init(&g_par.done, NULL);
    if (env != NULL && atoi(env) >= 1) {
        nc = atoi(env);
    }
    if (nc < 1) {
        nc = 1;
    }
    g_par.nt = (int)(nc - 1 < PAR_MAXT ? nc - 1 : PAR_MAXT);
    for (i = 0; i < g_par.nt; i++) {
        if (pthread_create(&g_par.th[i], NULL, par_worker, &g_par) != 0) {
            g_par.nt = i;
            break;
        }
    }
    g_par.ready = 1;
}

int stc_par_threads(void)
{
    pthread_once(&g_par_once, par_init);
    return g_par.nt + 1;
}

void stc_par_for(StcParFn fn, void *ctx, int n, int chunk)
{
    ParPool *p = &g_par;

    if (fn == NULL || n <= 0) {
        return;
    }
    pthread_once(&g_par_once, par_init);
    if (chunk < 1) {
        /* Auto: at least PAR_MIN_CHUNK indices per hand-off. Dispatch
         * costs tens of microseconds; loops over a few dozen bodies are
         * far cheaper than that, so they run serially on the caller. */
        chunk = n / (8 * (p->nt + 1));
        if (chunk < PAR_MIN_CHUNK) {
            chunk = PAR_MIN_CHUNK;
        }
    }
    if (p->nt == 0 || n <= chunk) {
        fn(ctx, 0, n);
        return;
    }
    pthread_mutex_lock(&p->mu);
    p->fn = fn;
    p->ctx = ctx;
    p->n = n;
    p->chunk = chunk;
    p->next = 0;
    p->pending = p->nt;
    p->gen++;
    pthread_cond_broadcast(&p->go);
    pthread_mutex_unlock(&p->mu);
    par_drain(p, fn, ctx, n);
    pthread_mutex_lock(&p->mu);
    while (p->pending > 0) {
        pthread_cond_wait(&p->done, &p->mu);
    }
    pthread_mutex_unlock(&p->mu);
}
