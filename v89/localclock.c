/* localclock.c — v89 — per-cell clocks: READINESS SUITE
 *
 * Standalone. Does not touch cellfab.c. Build:
 *   gcc -O2 -Wall -Wextra -o localclock localclock.c -lm
 *
 * WHAT THIS IS
 * ------------
 * Design (2b) of the 2026-08-02 thread: replace the one global dt with a
 * local time and tick counter per cell. The first pass established the
 * mechanism and found three things that are now specifications:
 *
 *   - bound skew in LOCAL TIME, not tick count: cells take
 *     different-sized steps, so equal tick counts are not equal time and
 *     a tick bound distorts the ordering it is meant to protect
 *     (agreement 1.3e-2 bounded vs 1.1e-4 unbounded);
 *   - a byte suffices up to skew 127 (K < M/2 exactly), but the wrap
 *     failure is SILENT — the comparison stops ordering and the bound
 *     evaporates with no symptom, so it must be asserted;
 *   - determinism survives only under a TOTAL canonical order; arrival
 *     order diverges up to 3.7 rad, and breaking ties by scan position
 *     instead of index leaks back in at 3e-6.
 *
 * This revision applies all three and asks whether the result is fit to
 * carry the battery. Four questions, in the order that can kill it:
 *
 *   R1 CONSERVATION. The phase model tested nothing of the kind. Energy
 *      moves on channels as paired antisymmetric transfers; does the
 *      ledger close under async as exactly as under sync? This is the
 *      One Law and it is not negotiable.
 *   R2 CONVERGENCE. Does |async - sync| fall with dt? If it plateaus the
 *      local scheme is a different model, not a different integrator.
 *   R3 DETERMINISM. Bit-identical under scan reversal AND rotated scan
 *      origin (stand-ins for "which core got there first").
 *   R4 LOOKAHEAD. The bound should be PHYSICAL: a neighbour cannot
 *      influence a cell faster than the rate limit, so that limit is the
 *      lookahead and no arbitrary K is needed. Plus the payoff — how
 *      many cells are eligible at once.
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>

#define TWOPI 6.283185307179586
#define MAXDEG 16

typedef struct {
    int N, deg;
    int *nb, *nd;
    double *th, *pub, *w, *E, *t, *h;
    long *tau;
    /* CHANNELS are first-class event holders with their own clock. A
     * transfer fired once per endpoint-advance double-counts: each edge
     * would move gE*(h_i+h_j) per round instead of gE per unit time, an
     * O(1) rate error that does NOT vanish as dt->0. Channel ownership
     * makes the rate correct for any step size. */
    int ne, *ei, *ej;
    int *ce, *cen;   /* per-cell incident EDGE ids: O(deg) conflict lookup */
    double *et, *eh;
    double Kc, gE, dt;
    long cap_hits;
} Fab;

static unsigned long long RS;
static double rnd(void)
{ RS ^= RS << 13; RS ^= RS >> 7; RS ^= RS << 17;
  return (double)((RS >> 11) & 0xFFFFFFFFFFFFFULL) / 9007199254740992.0; }

static Fab *fab_new(int N, int deg, double Kc, double gE, double dt)
{
    Fab *f = calloc(1, sizeof(Fab));
    f->N = N; f->deg = deg; f->Kc = Kc; f->gE = gE; f->dt = dt;
    f->nb = calloc((size_t)N * MAXDEG, sizeof(int));
    f->nd = calloc(N, sizeof(int));
    f->th = calloc(N, sizeof(double));
    f->pub = calloc(N, sizeof(double));
    f->w = calloc(N, sizeof(double));
    f->E = calloc(N, sizeof(double));
    f->t = calloc(N, sizeof(double));
    f->h = calloc(N, sizeof(double));
    f->tau = calloc(N, sizeof(long));
    f->ei = calloc((size_t)N * MAXDEG, sizeof(int));
    f->ej = calloc((size_t)N * MAXDEG, sizeof(int));
    f->et = calloc((size_t)N * MAXDEG, sizeof(double));
    f->eh = calloc((size_t)N * MAXDEG, sizeof(double));
    f->ce = calloc((size_t)N * MAXDEG, sizeof(int));
    f->cen = calloc(N, sizeof(int));
    return f;
}
static void fab_free(Fab *f)
{ free(f->nb); free(f->nd); free(f->th); free(f->pub); free(f->w);
  free(f->E); free(f->t); free(f->h); free(f->tau);
  free(f->ei); free(f->ej); free(f->et); free(f->eh);
  free(f->ce); free(f->cen); free(f); }

static void fab_build(Fab *f, unsigned long long seed)
{
    int N = f->N;
    RS = seed;
    for (int i = 0; i < N; i++) {
        f->th[i] = TWOPI * rnd();
        f->pub[i] = f->th[i];
        f->w[i] = 1.0 + 0.30 * (rnd() - 0.5);
        f->E[i] = 1.0 + rnd();
        f->h[i] = f->dt * (0.6 + 0.8 * rnd());   /* the dilation */
    }
    for (int i = 0; i < N; i++) {               /* ring backbone */
        int j = (i + 1) % N;
        f->nb[i * MAXDEG + f->nd[i]++] = j;
        f->nb[j * MAXDEG + f->nd[j]++] = i;
    }
    for (int e = 0; e < N * (f->deg - 2) / 2; e++) {
        int i = (int)(rnd() * N), j = (int)(rnd() * N);
        if (i == j || f->nd[i] >= MAXDEG || f->nd[j] >= MAXDEG) continue;
        int dup = 0;
        for (int k = 0; k < f->nd[i]; k++) if (f->nb[i * MAXDEG + k] == j) dup = 1;
        if (dup) continue;
        f->nb[i * MAXDEG + f->nd[i]++] = j;
        f->nb[j * MAXDEG + f->nd[j]++] = i;
    }
    f->ne = 0;
    for (int i = 0; i < N; i++)
        for (int k = 0; k < f->nd[i]; k++) {
            int j = f->nb[i * MAXDEG + k];
            if (j < i) continue;
            f->ei[f->ne] = i; f->ej[f->ne] = j;
            f->eh[f->ne] = 0.5 * (f->h[i] + f->h[j]);
            f->ne++;
        }
    memset(f->cen, 0, N * sizeof(int));
    for (int e = 0; e < f->ne; e++) {
        int a = f->ei[e], b = f->ej[e];
        if (f->cen[a] < MAXDEG) f->ce[a * MAXDEG + f->cen[a]++] = e;
        if (f->cen[b] < MAXDEG) f->ce[b * MAXDEG + f->cen[b]++] = e;
    }
}

static void fab_copy_ic(Fab *d, const Fab *s)
{
    int N = s->N;
    memcpy(d->nb, s->nb, (size_t)N * MAXDEG * sizeof(int));
    memcpy(d->nd, s->nd, N * sizeof(int));
    memcpy(d->th, s->th, N * sizeof(double));
    memcpy(d->pub, s->th, N * sizeof(double));
    memcpy(d->w, s->w, N * sizeof(double));
    memcpy(d->E, s->E, N * sizeof(double));
    memcpy(d->h, s->h, N * sizeof(double));
    memset(d->t, 0, N * sizeof(double));
    memset(d->tau, 0, N * sizeof(long));
    d->ne = s->ne;
    memcpy(d->ei, s->ei, s->ne * sizeof(int));
    memcpy(d->ej, s->ej, s->ne * sizeof(int));
    memcpy(d->eh, s->eh, s->ne * sizeof(double));
    memcpy(d->ce, s->ce, (size_t)N * MAXDEG * sizeof(int));
    memcpy(d->cen, s->cen, N * sizeof(int));
    memset(d->et, 0, s->ne * sizeof(double));
    d->cap_hits = 0;
}

static double total_E(Fab *f)
{ double s = 0; for (int i = 0; i < f->N; i++) s += f->E[i]; return s; }

static double dtheta(Fab *f, int i)
{
    double s = 0;
    for (int k = 0; k < f->nd[i]; k++)
        s += sin(f->pub[f->nb[i * MAXDEG + k]] - f->th[i]);
    return f->w[i] + f->Kc * s;
}

/* Advance ONE cell by its own step. Energy moves as PAIRED ANTISYMMETRIC
 * TRANSFERS: whatever leaves i arrives at j in the same statement, so
 * conservation is exact by construction and independent of execution
 * order. That is the whole reason to structure the update as transfers
 * rather than as independently recomputed cell states. */
static void advance(Fab *f, int i)
{
    double hi = f->h[i];
    f->th[i] += hi * dtheta(f, i);
    f->pub[i] = f->th[i];
    f->t[i] += hi;
    f->tau[i]++;
}

/* the channel's own event: one paired antisymmetric transfer, fired
 * once per h_e of the CHANNEL's time. Conserves exactly, and the rate
 * is gE per unit time whatever the endpoints' steps are. */
static void advance_edge(Fab *f, int e)
{
    int i = f->ei[e], j = f->ej[e];
    double q = f->eh[e] * f->gE * (f->E[i] - f->E[j]);
    f->E[i] -= q;
    f->E[j] += q;
    f->et[e] += f->eh[e];
}

static void run_sync(Fab *f, double T)
{
    int N = f->N;
    double *nth = malloc(N * sizeof(double));
    double *dE = malloc(N * sizeof(double));
    int steps = (int)(T / f->dt + 0.5);
    for (int s = 0; s < steps; s++) {
        memset(dE, 0, N * sizeof(double));
        for (int i = 0; i < N; i++) nth[i] = f->th[i] + f->dt * dtheta(f, i);
        for (int i = 0; i < N; i++)
            for (int k = 0; k < f->nd[i]; k++) {
                int j = f->nb[i * MAXDEG + k];
                if (j < i) continue;
                double q = f->dt * f->gE * (f->E[i] - f->E[j]);
                dE[i] -= q; dE[j] += q;
            }
        for (int i = 0; i < N; i++) { f->th[i] = nth[i]; f->E[i] += dE[i]; }
        memcpy(f->pub, f->th, N * sizeof(double));
        for (int i = 0; i < N; i++) f->tau[i]++;
    }
    free(nth); free(dE);
}

/* Unified event loop over CELLS and CHANNELS. An event is eligible when
 * its own local time stays within LOOKAHEAD of the minimum over the cells
 * it touches. LOOKAHEAD is physical — how long a neighbour needs to reach
 * this cell at the rate limit — so there is no arbitrary K.
 * sel 0 canonical (smallest (t, kind, index)); 1 arrival order.
 * rot/dir = scan origin and direction: under sel=0 neither may change a
 * single bit. `skew` returns the largest instantaneous neighbour tick
 * difference seen — the quantity a mod-M counter would have to hold. */
static void run_local(Fab *f, double T, double look, long M,
                      int sel, int rot, int dir, long *skew)
{
    int N = f->N, NE = f->ne, NT = N + NE;
    long cap = 400L * NT * (long)(T / f->dt + 2), guard = 0;
    if (skew) *skew = 0;
    for (;;) {
        int pick = -1;          /* < N : cell;  >= N : channel */
        double bt = 0;
        for (int q = 0; q < NT; q++) {
            int u = dir ? (rot + NT - 1 - q) % NT : (rot + q) % NT;
            double tu, lim = 1e300;
            if (u < N) {
                if (f->t[u] >= T) continue;
                tu = f->t[u];
                for (int k = 0; k < f->nd[u]; k++) {
                    double tj = f->t[f->nb[u * MAXDEG + k]];
                    if (tj < lim) lim = tj;
                }
                if (tu + f->h[u] > lim + look) continue;
            } else {
                int e = u - N;
                if (f->et[e] >= T) continue;
                tu = f->et[e];
                double a = f->t[f->ei[e]], b = f->t[f->ej[e]];
                lim = a < b ? a : b;
                if (tu + f->eh[e] > lim + look) continue;
            }
            if (M > 0 && u < N) {
                long d = 0;
                for (int k = 0; k < f->nd[u]; k++) {
                    long dd = f->tau[u] - f->tau[f->nb[u * MAXDEG + k]];
                    if (dd < 0) dd = -dd;
                    if (dd > d) d = dd;
                }
                if (d >= M / 2) {
                    fprintf(stderr, "FATAL: neighbour tick skew %ld >= M/2 "
                            "(%ld) — the wrapped counter can no longer order "
                            "ticks\n", d, M / 2);
                    exit(3);
                }
            }
            if (sel == 1) { pick = u; break; }
            if (pick < 0 || tu < bt || (tu == bt && u < pick))
                { pick = u; bt = tu; }
        }
        if (pick < 0) break;
        if (pick < N) {
            advance(f, pick);
            if (skew) {
                for (int k = 0; k < f->nd[pick]; k++) {
                    long dd = f->tau[pick] - f->tau[f->nb[pick * MAXDEG + k]];
                    if (dd < 0) dd = -dd;
                    if (dd > *skew) *skew = dd;
                }
            }
        } else {
            advance_edge(f, pick - N);
        }
        if (++guard > cap) { f->cap_hits++; break; }
    }
}

static double mean_eligible(Fab *f, double T, double look)
{
    int N = f->N, NE = f->ne, NT = N + NE;
    double acc = 0;
    long n = 0, cap = 400L * NT * (long)(T / f->dt + 2), guard = 0;
    for (;;) {
        int cnt = 0, pick = -1;
        double bt = 0;
        for (int u = 0; u < NT; u++) {
            double tu, lim = 1e300;
            if (u < N) {
                if (f->t[u] >= T) continue;
                tu = f->t[u];
                for (int k = 0; k < f->nd[u]; k++) {
                    double tj = f->t[f->nb[u * MAXDEG + k]];
                    if (tj < lim) lim = tj;
                }
                if (tu + f->h[u] > lim + look) continue;
            } else {
                int e = u - N;
                if (f->et[e] >= T) continue;
                tu = f->et[e];
                double a = f->t[f->ei[e]], b = f->t[f->ej[e]];
                lim = a < b ? a : b;
                if (tu + f->eh[e] > lim + look) continue;
            }
            cnt++;
            if (pick < 0 || tu < bt || (tu == bt && u < pick))
                { pick = u; bt = tu; }
        }
        if (pick < 0) break;
        acc += cnt; n++;
        if (pick < N) advance(f, pick); else advance_edge(f, pick - N);
        if (++guard > cap) break;
    }
    return n ? acc / n : 0;
}

/* ---- BATCH-PARALLEL execution ----
 * The serial loop advances one event at a time, which caps throughput at
 * one event per step however many cells are eligible. R4 showed nearly
 * every event is eligible at once, so the width is there — it just needs
 * a conflict-free selection.
 *
 * Two events CONFLICT if they touch a common cell. Events that do not
 * share a cell commute EXACTLY (not just to rounding: they write
 * disjoint memory), so a conflict-free batch can be advanced in any
 * order — including simultaneously — and give bit-identical results.
 *
 * Selection rule, deterministic and computable locally: an eligible
 * event is in the batch iff it has the minimum (t, index) among all
 * eligible events touching any of its cells. That is a local maximal
 * independent set, a pure function of state, so it does not depend on
 * thread count or scheduling. Returns rounds executed; *width = mean
 * batch size. */
static long run_batch(Fab *f, double T, double look, double *width,
                      int nthreads)
{
    int N = f->N, NE = f->ne, NT = N + NE;
    char *elig = malloc(NT), *sel = malloc(NT);
    double *evt = malloc(NT * sizeof(double));
    long rounds = 0; double acc = 0;
    for (;;) {
        int any = 0;
        #pragma omp parallel for num_threads(nthreads) schedule(static) reduction(|:any)
        for (int u = 0; u < NT; u++) {
            elig[u] = 0; evt[u] = 0;
            double tu, lim = 1e300;
            if (u < N) {
                if (f->t[u] >= T) continue;
                tu = f->t[u];
                for (int k = 0; k < f->nd[u]; k++) {
                    double tj = f->t[f->nb[u * MAXDEG + k]];
                    if (tj < lim) lim = tj;
                }
                if (tu + f->h[u] > lim + look) continue;
            } else {
                int e = u - N;
                if (f->et[e] >= T) continue;
                tu = f->et[e];
                double a = f->t[f->ei[e]], b = f->t[f->ej[e]];
                lim = a < b ? a : b;
                if (tu + f->eh[e] > lim + look) continue;
            }
            elig[u] = 1; evt[u] = tu; any = 1;
        }
        if (!any) break;
        /* local minimum over the conflict neighbourhood */
        #pragma omp parallel for num_threads(nthreads) schedule(static)
        for (int u = 0; u < NT; u++) {
            sel[u] = 0;
            if (!elig[u]) continue;
            int win = 1;
            int cells[2], nc;
            if (u < N) { cells[0] = u; nc = 1; }
            else { cells[0] = f->ei[u - N]; cells[1] = f->ej[u - N]; nc = 2; }
            for (int c = 0; c < nc && win; c++) {
                int i = cells[c];
                if (elig[i] && (evt[i] < evt[u] || (evt[i] == evt[u] && i < u)))
                    win = 0;
                /* A cell event READS its neighbours' published phase, so
                 * two ADJACENT cell events race even though they write
                 * disjoint memory. Excluding only incident edges left
                 * that race in and it showed up as thread-count
                 * dependence at 8 and 16 threads. Neighbouring cell
                 * events are conflicts too. */
                for (int k = 0; k < f->nd[i] && win; k++) {
                    int v = f->nb[i * MAXDEG + k];
                    if (v != u && elig[v] &&
                        (evt[v] < evt[u] || (evt[v] == evt[u] && v < u))) win = 0;
                }
                for (int k = 0; k < f->cen[i] && win; k++) {
                    int v = N + f->ce[i * MAXDEG + k];
                    if (v != u && elig[v] &&
                        (evt[v] < evt[u] || (evt[v] == evt[u] && v < u))) win = 0;
                }
            }
            sel[u] = win;
        }
        int cnt = 0;
        for (int u = 0; u < NT; u++) if (sel[u]) cnt++;
        if (!cnt) break;
        /* conflict-free by construction: disjoint memory, so parallel */
        #pragma omp parallel for num_threads(nthreads) schedule(static)
        for (int u = 0; u < NT; u++) {
            if (!sel[u]) continue;
            if (u < N) advance(f, u); else advance_edge(f, u - N);
        }
        acc += cnt; rounds++;
    }
    free(elig); free(sel); free(evt);
    if (width) *width = rounds ? acc / rounds : 0;
    return rounds;
}

static double maxdiff_th(Fab *a, Fab *b)
{ double m = 0; for (int i = 0; i < a->N; i++)
  { double d = fabs(a->th[i] - b->th[i]); if (d > m) m = d; } return m; }
static double maxdiff_E(Fab *a, Fab *b)
{ double m = 0; for (int i = 0; i < a->N; i++)
  { double d = fabs(a->E[i] - b->E[i]); if (d > m) m = d; } return m; }

static Fab *mk(Fab *ic, double dt)
{
    Fab *f = fab_new(ic->N, ic->deg, ic->Kc, ic->gE, dt);
    fab_copy_ic(f, ic);
    for (int i = 0; i < f->N; i++) f->h[i] = ic->h[i] * (dt / ic->dt);
    /* the CHANNEL steps must rescale too. Missing this froze the async
     * transport at the dt0 discretisation, so its error against a fine
     * reference sat at exactly 1.200e-4 for every dt — a plateau that
     * looks precisely like 'different fixed point' and is not. */
    for (int e = 0; e < f->ne; e++) f->eh[e] = ic->eh[e] * (dt / ic->dt);
    return f;
}

int main(void)
{
    const int N = 96, DEG = 7;
    const double Kc = 0.5, gE = 0.30, dt0 = 0.02, T = 40.0;
    const double LOOK = 0.05;

    printf("# localclock — READINESS SUITE\n");
    printf("# N=%d degree=%d Kc=%.2f gE=%.2f T=%.1f lookahead=%.3f\n",
           N, DEG, Kc, gE, T, LOOK);
    printf("# skew bounded in LOCAL TIME; order canonical (t_i,i);\n");
    printf("# energy moves as paired antisymmetric transfers.\n\n");

    Fab *ic = fab_new(N, DEG, Kc, gE, dt0);
    fab_build(ic, 20260802);
    { int s = 0; for (int i = 0; i < N; i++) s += ic->nd[i];
      printf("graph: mean degree %.2f\n\n", (double)s / N); }

    printf("== R1. CONSERVATION — does the ledger close under async? ==\n");
    printf("                   before          after          drift\n");
    {
        Fab *S = mk(ic, dt0); double e0 = total_E(S);
        run_sync(S, T);
        printf("SYNC           %14.10f  %14.10f  %+.3e\n",
               e0, total_E(S), total_E(S) - e0);
        fab_free(S);
        double looks[] = {0.02, 0.05, 0.2, 1.0};
        for (unsigned q = 0; q < 4; q++) {
            Fab *L = mk(ic, dt0); double b0 = total_E(L);
            run_local(L, T, looks[q], 0, 0, 0, 0, NULL);
            printf("LOCAL look=%-5.2f %14.10f  %14.10f  %+.3e%s\n",
                   looks[q], b0, total_E(L), total_E(L) - b0,
                   L->cap_hits ? "  (CAP HIT)" : "");
            fab_free(L);
        }
    }
    {   /* is the diffusion actually finished in BOTH schemes? If E is
         * equilibrated the two must agree to the FP floor; a residual
         * means one of them did not finish, which is a scheduler bug
         * and not a discretisation difference. */
        Fab *S = mk(ic, dt0); run_sync(S, T);
        Fab *L = mk(ic, dt0); run_local(L, T, LOOK, 0, 0, 0, 0, NULL);
        double smn = S->E[0], smx = S->E[0], lmn = L->E[0], lmx = L->E[0];
        long emn = 0, emx = 0;
        for (int i = 0; i < N; i++) {
            if (S->E[i] < smn) smn = S->E[i];
            if (S->E[i] > smx) smx = S->E[i];
            if (L->E[i] < lmn) lmn = L->E[i];
            if (L->E[i] > lmx) lmx = L->E[i];
        }
        for (int e = 0; e < L->ne; e++) {
            long k = (long)(L->et[e] / L->eh[e] + 0.5);
            if (e == 0 || k < emn) emn = k;
            if (e == 0 || k > emx) emx = k;
        }
        printf("\nequilibration check (E should be uniform by T=%.0f):\n", T);
        printf("  SYNC  E spread = %.3e\n", smx - smn);
        printf("  LOCAL E spread = %.3e\n", lmx - lmn);
        printf("  LOCAL channel firings: min %ld max %ld\n", emn, emx);
        fab_free(S); fab_free(L);
    }

    printf("\nreading: paired transfers conserve by construction, so this\n");
    printf("must sit at the FP floor for EVERY lookahead. Anything else\n");
    printf("means the transfer structure is wrong, not the scheduler.\n\n");

    printf("== R2. CONVERGENCE — is it the SAME physics? ==\n");
    printf("Comparing async(dt) to sync(dt) conflates two independent O(dt)\n");
    printf("errors and cannot answer this. Both are measured instead against\n");
    printf("a common fine-step reference: sync at dt0/64.\n");
    printf("     dt    sync err     async err   sync ratio  async ratio\n");
    {
        Fab *R = mk(ic, dt0 / 64.0); run_sync(R, T);
        double ps = 0, pa = 0;
        double dts[] = {0.02, 0.01, 0.005, 0.0025, 0.00125};
        for (unsigned q = 0; q < 5; q++) {
            Fab *S = mk(ic, dts[q]); run_sync(S, T);
            Fab *L = mk(ic, dts[q]); run_local(L, T, LOOK, 0, 0, 0, 0, NULL);
            double es = maxdiff_E(R, S), ea = maxdiff_E(R, L);
            if (q) printf("%8.5f  %10.3e  %10.3e  %9.2fx  %10.2fx\n",
                          dts[q], es, ea, ps / es, pa / ea);
            else   printf("%8.5f  %10.3e  %10.3e       —            —\n",
                          dts[q], es, ea);
            ps = es; pa = ea;
            fab_free(S); fab_free(L);
        }
        fab_free(R);
    }
    printf("\nreading: BOTH columns must fall, at comparable order. That is\n");
    printf("what 'same physics, different integrator' looks like. If the\n");
    printf("async column plateaus while sync falls, the local scheme has a\n");
    printf("different fixed point and is not adoptable.\n\n");

    printf("== R3. DETERMINISM — bit-identity under reordering ==\n");
    printf("  variant                        max|dtheta|   verdict\n");
    {
        Fab *A = mk(ic, dt0); run_local(A, T, LOOK, 0, 0, 0, 0, NULL);
        struct { const char *name; int sel, rot, dir; } v[] = {
            {"canonical, reversed scan",      0,  0, 1},
            {"canonical, rotated origin 37",  0, 37, 0},
            {"canonical, rotated + reversed", 0, 53, 1},
            {"ARRIVAL order, reversed scan",  1,  0, 1},
        };
        for (unsigned q = 0; q < 4; q++) {
            Fab *B = mk(ic, dt0);
            run_local(B, T, LOOK, 0, v[q].sel, v[q].rot, v[q].dir, NULL);
            double d = maxdiff_th(A, B);
            printf("  %-29s  %11.3e   %s\n", v[q].name, d,
                   d == 0.0 ? "EXACT" : "differs");
            fab_free(B);
        }
        fab_free(A);
    }
    printf("\nreading: every canonical row must read EXACT — that property is\n");
    printf("what lets the battery compare byte-identical reruns.\n\n");

    printf("== R4. LOOKAHEAD — physical bound, and what it buys ==\n");
    printf(" lookahead  mean_eligible  neighbour skew  total dilation  bits\n");
    {
        double looks[] = {0.02, 0.05, 0.1, 0.2, 0.5, 1.0, 2.0};
        for (unsigned q = 0; q < 7; q++) {
            Fab *L = mk(ic, dt0);
            double me = mean_eligible(L, T, looks[q]);
            fab_free(L);
            Fab *L2 = mk(ic, dt0); long sk = 0;
            run_local(L2, T, looks[q], 0, 0, 0, 0, &sk);
            long mn = L2->tau[0], mx = L2->tau[0];
            for (int i = 0; i < N; i++) {
                if (L2->tau[i] < mn) mn = L2->tau[i];
                if (L2->tau[i] > mx) mx = L2->tau[i];
            }
            int bits = 1; while ((1L << bits) < 2 * (sk + 1)) bits++;
            printf("%10.3f  %13.2f  %14ld  %13ld  %8d\n",
                   looks[q], me, sk, mx - mn, bits);
            fab_free(L2);
        }
    }
    printf("\nreading: mean_eligible is the width available to run in parallel\n");
    printf("with no global barrier. 'counter bits' is what the K<M/2 rule\n");
    printf("demands at that lookahead — the answer to 'is a byte enough'.\n");

    printf("\n== R5. PARALLEL WIDTH — can this be batched? ==\n");
    printf("Conflict-free batch: events sharing no cell write disjoint\n");
    printf("memory, so they commute EXACTLY and may run simultaneously.\n");
    printf("Selection is a local minimum over the conflict neighbourhood —\n");
    printf("a pure function of state, so batches do not depend on threads.\n");
    printf(" lookahead   rounds   mean batch   speedup vs serial   matches serial\n");
    {
        Fab *ser = mk(ic, dt0);
        run_local(ser, T, LOOK, 0, 0, 0, 0, NULL);
        double looks[] = {0.05, 0.2, 1.0};
        for (unsigned q = 0; q < 3; q++) {
            Fab *B = mk(ic, dt0); double wd;
            long r = run_batch(B, T, looks[q], &wd, 1);
            long nev = 0;
            for (int i = 0; i < N; i++) nev += B->tau[i];
            for (int e = 0; e < B->ne; e++) nev += (long)(B->et[e] / B->eh[e] + 0.5);
            double d = (looks[q] == LOOK) ? maxdiff_th(ser, B) : -1;
            printf("%10.3f  %7ld  %11.2f  %17.2fx   %s\n",
                   looks[q], r, wd, (double)nev / r,
                   d < 0 ? "n/a" : (d == 0.0 ? "EXACT" : "differs"));
            fab_free(B);
        }
        fab_free(ser);
    }
    printf("\nreading: mean batch is how many events a GPU could retire per\n");
    printf("round. It is the ceiling on any parallel speedup here, and it\n");
    printf("must come with EXACT agreement against the serial run — a batch\n");
    printf("that is merely close is a different scheme.\n");

    printf("\n== R6. DOES THE WIDTH SCALE WITH N? (decides GPU viability) ==\n");
    printf("A GPU needs thousands of independent work items. A maximal\n");
    printf("independent set in a fixed-degree conflict graph is a constant\n");
    printf("FRACTION of the events, so width should grow linearly with N.\n");
    printf("If it saturates, no GPU can help and the answer is no.\n");
    printf("      N   events   rounds   mean batch   batch/event frac\n");
    {
        int Ns[] = {48, 96, 192, 384, 768, 1536};
        for (unsigned q = 0; q < 6; q++) {
            Fab *G = fab_new(Ns[q], DEG, Kc, gE, dt0);
            fab_build(G, 20260802);
            Fab *B = mk(G, dt0); double wd;
            long r = run_batch(B, 4.0, LOOK, &wd, 1);
            int NT = B->N + B->ne;
            printf("%7d  %7d  %7ld  %11.2f  %16.4f\n",
                   Ns[q], NT, r, wd, wd / NT);
            fab_free(B); fab_free(G);
        }
    }
    printf("\nreading: a constant batch/event fraction with mean batch rising\n");
    printf("linearly means the width is there at scale and a GPU port is\n");
    printf("worth building. A falling fraction means the conflict graph is\n");
    printf("too dense and the scheme is latency-bound whatever the hardware.\n");

    printf("\n== R7. THREAD-COUNT INVARIANCE (the ratchet requirement) ==\n");
    printf(" threads   bit-identical to 1 thread\n");
    {
        Fab *ref = mk(ic, dt0); double wd;
        run_batch(ref, T, LOOK, &wd, 1);
        int ths[] = {2, 4, 8, 16};
        for (unsigned q = 0; q < 4; q++) {
            Fab *B = mk(ic, dt0);
            run_batch(B, T, LOOK, &wd, ths[q]);
            printf("%8d   %s\n", ths[q],
                   maxdiff_th(ref, B) == 0.0 && maxdiff_E(ref, B) == 0.0
                   ? "YES" : "NO");
            fab_free(B);
        }
        fab_free(ref);
    }
    printf("\nreading: the selection rule is a pure function of state, so this\n");
    printf("must read YES at every thread count. It is the property that lets\n");
    printf("the battery compare byte-identical reruns on any machine.\n");

    fab_free(ic);
    return 0;
}
