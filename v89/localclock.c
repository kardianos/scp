/* localclock.c — v89 — does the global clock have to exist?
 *
 * Standalone. Does not touch cellfab.c. Build:
 *   gcc -O2 -Wall -Wextra -o localclock localclock.c -lm
 *
 * THE QUESTION
 * ------------
 * Idea (2b) of the 2026-08-02 design thread: amortise the global clock by
 * giving every cell its OWN tick counter — "a byte with a mod" — and
 * letting each advance at its own rate, relativistically. The v89 kernel
 * currently advances every cell together on one dt, which is a permanent
 * index set in time: something that persists and is merely re-valued,
 * which is the construction PRINCIPLE constraint 2 forbids in space. If
 * the argument is good against a spatial background it is good against a
 * temporal one.
 *
 * The idea is well-motivated beyond tidiness: v89 ALREADY has a local
 * rate. The pitch law x = (Em + fl*load)/cap flattens a loaded cell's
 * pitch, and pitch IS rate. So the physics is already local-rate; only
 * the integrator is global. Making the tick local formalises what the
 * pitch law already does, and makes time dilation structural rather than
 * an effect computed on top of a universal clock.
 *
 * WHAT THIS FILE MEASURES
 * -----------------------
 * A ring of coupled phase cells — the minimal stand-in for the dense
 * sector — integrated three ways on identical initial conditions:
 *
 *   SYNC   every cell advances on one global dt (what cellfab does now)
 *   LOCAL  every cell owns a counter tau_i; it may advance only while it
 *          is no more than K ticks ahead of its neighbours; it reads
 *          whatever its neighbours have PUBLISHED, which may be stale
 *   WRAP   the same, but with the counter stored in a narrow field of
 *          M values, so tau wraps
 *
 * Three questions, each with a number attached:
 *   Q1  does going local change the physics? (compare locked state)
 *   Q2  how much skew K can be tolerated before it does?
 *   Q3  how many bits does the counter need before wrap destroys
 *       causality? (the "is a byte enough" question, answered)
 *
 * and one that decides whether this is adoptable at all:
 *   Q4  is LOCAL still DETERMINISTIC — same answer regardless of the
 *       order cells happen to be visited in? The ratchet depends on
 *       byte-identical reruns; an update rule that loses that is not
 *       adoptable no matter how principled it is.
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define TWOPI 6.283185307179586

typedef struct {
    int N;
    double *th;      /* phase */
    double *pub;     /* last PUBLISHED phase (what neighbours can see) */
    double *w;       /* natural pitch = tick rate */
    long   *tau;     /* local tick counter */
    double Kc, dt;
} Ring;

static Ring *ring_new(int N, double Kc, double dt, unsigned long long seed)
{
    Ring *r = calloc(1, sizeof(Ring));
    r->N = N; r->Kc = Kc; r->dt = dt;
    r->th = malloc(N * sizeof(double));
    r->pub = malloc(N * sizeof(double));
    r->w = malloc(N * sizeof(double));
    r->tau = calloc(N, sizeof(long));
    unsigned long long s = seed ? seed : 1;
    for (int i = 0; i < N; i++) {
        s ^= s << 13; s ^= s >> 7; s ^= s << 17;
        double a = (double)((s >> 11) & 0xFFFFFFFFFFFFFULL) / 9007199254740992.0;
        s ^= s << 13; s ^= s >> 7; s ^= s << 17;
        double b = (double)((s >> 11) & 0xFFFFFFFFFFFFFULL) / 9007199254740992.0;
        r->th[i] = TWOPI * a;
        r->pub[i] = r->th[i];
        r->w[i] = 1.0 + 0.30 * (b - 0.5);   /* spread of local rates */
    }
    return r;
}
static void ring_free(Ring *r)
{ free(r->th); free(r->pub); free(r->w); free(r->tau); free(r); }

static void ring_copy_ic(Ring *d, const Ring *s)
{
    memcpy(d->th, s->th, s->N * sizeof(double));
    memcpy(d->pub, s->th, s->N * sizeof(double));
    memcpy(d->w, s->w, s->N * sizeof(double));
    memset(d->tau, 0, s->N * sizeof(long));
}

/* one cell's increment, from PUBLISHED neighbour phases */
static double dtheta(Ring *r, int i)
{
    int N = r->N, a = (i + N - 1) % N, b = (i + 1) % N;
    return r->w[i] + r->Kc * (sin(r->pub[a] - r->th[i])
                              + sin(r->pub[b] - r->th[i]));
}

/* order parameter of the phase differences: 1 = fully locked */
static double locked_R(Ring *r)
{
    double c = 0, s = 0;
    for (int i = 0; i < r->N; i++) {
        double d = r->th[(i + 1) % r->N] - r->th[i];
        c += cos(d); s += sin(d);
    }
    return sqrt(c * c + s * s) / r->N;
}

/* ---- SYNC: one global dt, read-old/write-new (Jacobi) ---- */
static void run_sync(Ring *r, double T)
{
    int N = r->N;
    double *nx = malloc(N * sizeof(double));
    int steps = (int)(T / r->dt + 0.5);
    for (int t = 0; t < steps; t++) {
        for (int i = 0; i < N; i++) nx[i] = r->th[i] + r->dt * dtheta(r, i);
        memcpy(r->th, nx, N * sizeof(double));
        memcpy(r->pub, nx, N * sizeof(double));
        for (int i = 0; i < N; i++) r->tau[i]++;
    }
    free(nx);
}

/* ---- LOCAL: event-driven, per-cell proper time ----
 * Every cell owns a local time t_i and a local step h_i. Cells with
 * different pitch take different steps, so after the same elapsed
 * coordinate time they have done DIFFERENT numbers of ticks — the
 * counters genuinely diverge, which is the whole point and is what the
 * wrap question is actually about.
 *
 * sel = 0 CANONICAL: always advance the cell with the smallest
 *        (t_i, i). The order is a function of STATE, so it does not
 *        matter which core gets there first.
 * sel = 1 ARRIVAL: advance whichever eligible cell the scan reaches
 *        first. This is what you get if you simply let threads run.
 * scan   direction of the scan; under sel=0 it must not matter, under
 *        sel=1 it will.                                                */
static void run_local(Ring *r, double T, int K, long M, int sel, int scan)
{
    int N = r->N;
    double *h = malloc(N * sizeof(double));
    double *t = calloc(N, sizeof(double));
    for (int i = 0; i < N; i++) h[i] = r->dt * (0.6 + 0.8 * (r->w[i] - 0.85));
    long guard = 0, cap = 4000L * N * (long)(T / r->dt + 2);
    for (;;) {
        int pick = -1;
        double bt = 0;
        for (int q = 0; q < N; q++) {
            int i = scan ? (N - 1 - q) : q;
            if (t[i] >= T) continue;
            int a = (i + N - 1) % N, b = (i + 1) % N;
            long d1, d2;
            if (M > 0) {
                d1 = ((r->tau[i] - r->tau[a]) % M + M) % M; if (d1 > M / 2) d1 -= M;
                d2 = ((r->tau[i] - r->tau[b]) % M + M) % M; if (d2 > M / 2) d2 -= M;
            } else {
                d1 = r->tau[i] - r->tau[a];
                d2 = r->tau[i] - r->tau[b];
            }
            if (d1 > K || d2 > K) continue;      /* too far ahead: wait */
            if (sel == 1) { pick = i; break; }   /* arrival order */
            /* the order must be a TOTAL function of state: ties broken
             * by index, never by who was scanned first. Without the tie
             * clause this is 'almost canonical' and leaks scan order back
             * in at ~1e-6 — small, but not byte-identical, and the
             * ratchet compares bytes. */
            if (pick < 0 || t[i] < bt || (t[i] == bt && i < pick))
                { pick = i; bt = t[i]; }
        }
        if (pick < 0) break;
        r->th[pick] += h[pick] * dtheta(r, pick);
        r->pub[pick] = r->th[pick];
        t[pick] += h[pick];
        r->tau[pick]++;
        if (++guard > cap) break;
    }
    free(h); free(t);
}

static double maxdiff(Ring *a, Ring *b)
{
    double m = 0;
    for (int i = 0; i < a->N; i++) {
        double d = fabs(a->th[i] - b->th[i]);
        if (d > m) m = d;
    }
    return m;
}

int main(void)
{
    const int N = 24;
    const double Kc = 0.6, dt = 0.02, T = 80.0;

    printf("# localclock — must the clock be global?\n");
    printf("# ring N=%d, Kc=%.2f, dt=%.3f, elapsed T=%.1f\n", N, Kc, dt, T);
    printf("# LOCAL cells own a local time and a local STEP, so their tick\n");
    printf("# counters genuinely diverge (that is the dilation).\n\n");

    Ring *ic = ring_new(N, Kc, dt, 20260802);
    Ring *S = ring_new(N, Kc, dt, 1); ring_copy_ic(S, ic);
    run_sync(S, T);
    printf("== reference ==\nSYNC  R=%.6f\n\n", locked_R(S));

    printf("== Q1/Q2. does the LOCK survive going local? ==\n");
    printf("   K    R(local)   dR vs SYNC   tick spread (max-min tau)\n");
    {
        int Ks[] = {1, 2, 4, 8, 16, 32, 64, 1000000};
        for (unsigned q = 0; q < sizeof(Ks) / sizeof(Ks[0]); q++) {
            Ring *L = ring_new(N, Kc, dt, 1); ring_copy_ic(L, ic);
            run_local(L, T, Ks[q], 0, 0, 0);
            long mn = L->tau[0], mx = L->tau[0];
            for (int i = 0; i < N; i++) {
                if (L->tau[i] < mn) mn = L->tau[i];
                if (L->tau[i] > mx) mx = L->tau[i];
            }
            printf("%4d    %8.6f   %+10.2e   %ld\n",
                   Ks[q], locked_R(L), locked_R(L) - locked_R(S), mx - mn);
            ring_free(L);
        }
    }
    printf("\nreading: the lock is the physics. If R survives, a global\n");
    printf("clock was never carrying it.\n\n");

    printf("== Q3. IS A BYTE ENOUGH? counter width M vs skew bound K ==\n");
    printf("wrapped comparison can only order ticks while K < M/2.\n");
    printf("    M     K   K<M/2    R(local)   tick spread\n");
    {
        long Ms[] = {8, 8, 16, 16, 64, 256, 256};
        int  Ks[] = {2, 6,  4, 12, 40,  64, 200};
        for (unsigned q = 0; q < sizeof(Ms) / sizeof(Ms[0]); q++) {
            Ring *L = ring_new(N, Kc, dt, 1); ring_copy_ic(L, ic);
            run_local(L, T, Ks[q], Ms[q], 0, 0);
            long mn = L->tau[0], mx = L->tau[0];
            for (int i = 0; i < N; i++) {
                if (L->tau[i] < mn) mn = L->tau[i];
                if (L->tau[i] > mx) mx = L->tau[i];
            }
            printf("%5ld  %4d   %-6s  %8.6f   %ld\n", Ms[q], Ks[q],
                   (Ks[q] < Ms[q] / 2) ? "yes" : "NO", locked_R(L), mx - mn);
            ring_free(L);
        }
    }
    printf("\n");

    printf("== Q4. DETERMINISM — the adoptability test ==\n");
    printf("Same run, opposite scan direction. CANONICAL advances the cell\n");
    printf("with the smallest (t_i, i); ARRIVAL advances whichever eligible\n");
    printf("cell the scan meets first (what threads give you for free).\n");
    printf("   K    canonical |fwd-rev|   arrival |fwd-rev|\n");
    {
        int Ks[] = {1, 2, 4, 8, 16, 64};
        for (unsigned q = 0; q < sizeof(Ks) / sizeof(Ks[0]); q++) {
            Ring *A = ring_new(N, Kc, dt, 1); ring_copy_ic(A, ic);
            Ring *B = ring_new(N, Kc, dt, 1); ring_copy_ic(B, ic);
            run_local(A, T, Ks[q], 0, 0, 0);
            run_local(B, T, Ks[q], 0, 0, 1);
            double dc = maxdiff(A, B);
            Ring *C = ring_new(N, Kc, dt, 1); ring_copy_ic(C, ic);
            Ring *D = ring_new(N, Kc, dt, 1); ring_copy_ic(D, ic);
            run_local(C, T, Ks[q], 0, 1, 0);
            run_local(D, T, Ks[q], 0, 1, 1);
            double da = maxdiff(C, D);
            printf("%4d   %8.2e %-8s  %8.2e %s\n", Ks[q],
                   dc, dc == 0.0 ? "(exact)" : "(DIFFERS)",
                   da, da == 0.0 ? "(exact)" : "(DIFFERS)");
            ring_free(A); ring_free(B); ring_free(C); ring_free(D);
        }
    }
    printf("\nreading: the ratchet compares byte-identical reruns. If\n");
    printf("CANONICAL is exact and ARRIVAL is not, then local clocks are\n");
    printf("adoptable but ONLY with an order defined by state, never by\n");
    printf("which core arrived first.\n");

    ring_free(S); ring_free(ic);
    return 0;
}
