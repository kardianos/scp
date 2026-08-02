/* localclock.cu — v89 — GPU batch execution of the local-clock scheme.
 *
 * Build: nvcc -O3 -arch=sm_70 -o localclock_cu localclock.cu -lm
 * Run:   ./localclock_cu
 *
 * WHY THIS CAN BE PARALLEL AT ALL
 * -------------------------------
 * The serial scheme advances one event at a time (smallest (t,index)),
 * which caps throughput at one event per step. But events that share no
 * CELL touch disjoint memory and commute exactly, so a conflict-free set
 * may be advanced simultaneously and gives bit-identical results.
 *
 * The conflict neighbourhood of an event is NOT just "events on the same
 * cell". A cell event READS its neighbours' published phase, so two
 * ADJACENT cell events race even though they write disjoint memory.
 * On CPU that omission showed up only at 8+ threads. The rule used here:
 *
 *   cell event u   conflicts with: incident channel events,
 *                                  and cell events on adjacent cells
 *   channel event e conflicts with: events on either endpoint,
 *                                  and channel events sharing an endpoint
 *
 * An eligible event enters the batch iff it is the minimum (t, index)
 * over its conflict neighbourhood — a local maximal independent set,
 * a pure function of state, hence independent of thread count, block
 * size and launch configuration.
 *
 * CPU measurement (localclock.c, R6) says the batch width grows with N
 * (5.2 at N=48 -> 83.5 at N=1536, ~N^0.8), so the width a GPU needs is
 * there at production size even though it is only ~8 at test size.
 * This file measures whether that translates into real acceleration.
 */

#include <cuda_runtime.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAXDEG 16
#define CK(x) do { cudaError_t e_ = (x); if (e_ != cudaSuccess) { \
    fprintf(stderr, "CUDA %s @%d: %s\n", #x, __LINE__, \
    cudaGetErrorString(e_)); exit(2); } } while (0)

/* ------------------------------------------------------------- host graph */

typedef struct {
    int N, ne;
    int *nb, *nd, *ei, *ej, *ce, *cen;
    double *th, *pub, *w, *E, *t, *h, *et, *eh;
    double Kc, gE, dt;
} Fab;

static unsigned long long RS;
static double rnd(void)
{ RS ^= RS << 13; RS ^= RS >> 7; RS ^= RS << 17;
  return (double)((RS >> 11) & 0xFFFFFFFFFFFFFULL) / 9007199254740992.0; }

static Fab *fab_new(int N)
{
    Fab *f = (Fab *)calloc(1, sizeof(Fab));
    f->N = N;
    f->nb  = (int *)calloc((size_t)N * MAXDEG, sizeof(int));
    f->nd  = (int *)calloc(N, sizeof(int));
    f->ei  = (int *)calloc((size_t)N * MAXDEG, sizeof(int));
    f->ej  = (int *)calloc((size_t)N * MAXDEG, sizeof(int));
    f->ce  = (int *)calloc((size_t)N * MAXDEG, sizeof(int));
    f->cen = (int *)calloc(N, sizeof(int));
    f->th  = (double *)calloc(N, sizeof(double));
    f->pub = (double *)calloc(N, sizeof(double));
    f->w   = (double *)calloc(N, sizeof(double));
    f->E   = (double *)calloc(N, sizeof(double));
    f->t   = (double *)calloc(N, sizeof(double));
    f->h   = (double *)calloc(N, sizeof(double));
    f->et  = (double *)calloc((size_t)N * MAXDEG, sizeof(double));
    f->eh  = (double *)calloc((size_t)N * MAXDEG, sizeof(double));
    return f;
}
static void fab_free(Fab *f)
{ free(f->nb); free(f->nd); free(f->ei); free(f->ej); free(f->ce);
  free(f->cen); free(f->th); free(f->pub); free(f->w); free(f->E);
  free(f->t); free(f->h); free(f->et); free(f->eh); free(f); }

static void fab_build(Fab *f, int deg, double Kc, double gE, double dt,
                      unsigned long long seed)
{
    int N = f->N;
    f->Kc = Kc; f->gE = gE; f->dt = dt;
    RS = seed;
    for (int i = 0; i < N; i++) {
        f->th[i] = 6.283185307179586 * rnd();
        f->pub[i] = f->th[i];
        f->w[i] = 1.0 + 0.30 * (rnd() - 0.5);
        f->E[i] = 1.0 + rnd();
        f->h[i] = dt * (0.6 + 0.8 * rnd());
    }
    f->ne = 0;
    memset(f->nd, 0, N * sizeof(int));
    for (int i = 0; i < N; i++) {
        int j = (i + 1) % N;
        f->nb[i * MAXDEG + f->nd[i]++] = j;
        f->nb[j * MAXDEG + f->nd[j]++] = i;
    }
    for (int a = 0; a < N * (deg - 2) / 2 * 20; a++) {
        int cnt = 0;
        for (int i = 0; i < N; i++) cnt += f->nd[i];
        if (cnt >= N * deg) break;
        int i = (int)(rnd() * N), j = (int)(rnd() * N);
        if (i == j || f->nd[i] >= deg || f->nd[j] >= deg) continue;
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

/* ------------------------------------------------------- CPU reference */

static void cpu_batch(Fab *f, double T, double look, long *rounds, double *width)
{
    int N = f->N, NE = f->ne, NT = N + NE;
    char *elig = (char *)malloc(NT), *sel = (char *)malloc(NT);
    double *evt = (double *)malloc(NT * sizeof(double));
    long r = 0; double acc = 0;
    for (;;) {
        int any = 0;
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
        for (int u = 0; u < NT; u++) {
            sel[u] = 0;
            if (!elig[u]) continue;
            int cells[2], nc;
            if (u < N) { cells[0] = u; nc = 1; }
            else { cells[0] = f->ei[u - N]; cells[1] = f->ej[u - N]; nc = 2; }
            int win = 1;
            for (int c = 0; c < nc && win; c++) {
                int i = cells[c];
                if (elig[i] && (evt[i] < evt[u] || (evt[i] == evt[u] && i < u)))
                    win = 0;
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
            sel[u] = (char)win;
        }
        int cnt = 0;
        for (int u = 0; u < NT; u++) if (sel[u]) cnt++;
        if (!cnt) break;
        for (int u = 0; u < NT; u++) {
            if (!sel[u]) continue;
            if (u < N) {
                double s = 0;
                for (int k = 0; k < f->nd[u]; k++)
                    s += sin(f->pub[f->nb[u * MAXDEG + k]] - f->th[u]);
                f->th[u] += f->h[u] * (f->w[u] + f->Kc * s);
                f->pub[u] = f->th[u];
                f->t[u] += f->h[u];
            } else {
                int e = u - N, i = f->ei[e], j = f->ej[e];
                double q = f->eh[e] * f->gE * (f->E[i] - f->E[j]);
                f->E[i] -= q; f->E[j] += q;
                f->et[e] += f->eh[e];
            }
        }
        acc += cnt; r++;
    }
    free(elig); free(sel); free(evt);
    *rounds = r; *width = r ? acc / r : 0;
}

/* ------------------------------------------------------------- device */

__global__ void k_elig(int N, int NE, double T, double look,
                       const int *nb, const int *nd, const int *ei, const int *ej,
                       const double *t, const double *h, const double *et,
                       const double *eh, char *elig, double *evt, int *any)
{
    int u = blockIdx.x * blockDim.x + threadIdx.x;
    int NT = N + NE;
    if (u >= NT) return;
    elig[u] = 0; evt[u] = 0.0;
    double tu, lim = 1e300;
    if (u < N) {
        if (t[u] >= T) return;
        tu = t[u];
        for (int k = 0; k < nd[u]; k++) {
            double tj = t[nb[u * MAXDEG + k]];
            if (tj < lim) lim = tj;
        }
        if (tu + h[u] > lim + look) return;
    } else {
        int e = u - N;
        if (et[e] >= T) return;
        tu = et[e];
        double a = t[ei[e]], b = t[ej[e]];
        lim = a < b ? a : b;
        if (tu + eh[e] > lim + look) return;
    }
    elig[u] = 1; evt[u] = tu;
    atomicOr(any, 1);
}

__global__ void k_select(int N, int NE,
                         const int *nb, const int *nd, const int *ei,
                         const int *ej, const int *ce, const int *cen,
                         const char *elig, const double *evt, char *sel, int *cnt)
{
    int u = blockIdx.x * blockDim.x + threadIdx.x;
    int NT = N + NE;
    if (u >= NT) return;
    sel[u] = 0;
    if (!elig[u]) return;
    int cells[2], nc;
    if (u < N) { cells[0] = u; nc = 1; }
    else { cells[0] = ei[u - N]; cells[1] = ej[u - N]; nc = 2; }
    int win = 1;
    for (int c = 0; c < nc && win; c++) {
        int i = cells[c];
        if (elig[i] && (evt[i] < evt[u] || (evt[i] == evt[u] && i < u))) win = 0;
        for (int k = 0; k < nd[i] && win; k++) {
            int v = nb[i * MAXDEG + k];
            if (v != u && elig[v] &&
                (evt[v] < evt[u] || (evt[v] == evt[u] && v < u))) win = 0;
        }
        for (int k = 0; k < cen[i] && win; k++) {
            int v = N + ce[i * MAXDEG + k];
            if (v != u && elig[v] &&
                (evt[v] < evt[u] || (evt[v] == evt[u] && v < u))) win = 0;
        }
    }
    sel[u] = (char)win;
    if (win) atomicAdd(cnt, 1);
}

__global__ void k_advance(int N, int NE, double Kc, double gE,
                          const int *nb, const int *nd, const int *ei,
                          const int *ej, const char *sel,
                          double *th, double *pub, const double *w,
                          double *E, double *t, const double *h,
                          double *et, const double *eh)
{
    int u = blockIdx.x * blockDim.x + threadIdx.x;
    int NT = N + NE;
    if (u >= NT || !sel[u]) return;
    if (u < N) {
        double s = 0;
        for (int k = 0; k < nd[u]; k++) s += sin(pub[nb[u * MAXDEG + k]] - th[u]);
        th[u] += h[u] * (w[u] + Kc * s);
        t[u] += h[u];
    } else {
        int e = u - N, i = ei[e], j = ej[e];
        double q = eh[e] * gE * (E[i] - E[j]);
        E[i] -= q; E[j] += q;   /* conflict-free: i,j touched by no other */
        et[e] += eh[e];
    }
}

/* publish is a separate pass so no cell reads a phase mid-update */
__global__ void k_publish(int N, const char *sel, const double *th, double *pub)
{
    int u = blockIdx.x * blockDim.x + threadIdx.x;
    if (u < N && sel[u]) pub[u] = th[u];
}

int main(void)
{
    const int DEG = 7;
    const double Kc = 0.5, gE = 0.30, dt = 0.02, T = 4.0, LOOK = 0.05;

    cudaDeviceProp p;
    CK(cudaGetDeviceProperties(&p, 0));
    printf("# localclock.cu — GPU batch execution\n");
    printf("# device: %s, SMs=%d, cc=%d.%d\n", p.name, p.multiProcessorCount,
           p.major, p.minor);
    printf("# degree=%d T=%.1f lookahead=%.3f\n\n", DEG, T, LOOK);
    printf("      N   events   rounds   mean_batch   cpu_s    gpu_s  speedup  agree\n");

    int Ns[] = {1536, 6144, 24576, 98304, 393216};
    for (unsigned q = 0; q < sizeof(Ns) / sizeof(Ns[0]); q++) {
        int N = Ns[q];
        Fab *C = fab_new(N); fab_build(C, DEG, Kc, gE, dt, 20260802);
        Fab *G = fab_new(N); fab_build(G, DEG, Kc, gE, dt, 20260802);
        int NE = C->ne, NT = N + NE;

        /* ---- CPU ---- */
        long rounds = 0; double width = 0;
        cudaEvent_t a, b; CK(cudaEventCreate(&a)); CK(cudaEventCreate(&b));
        float cpu_ms = 0, gpu_ms = 0;
        CK(cudaEventRecord(a));
        cpu_batch(C, T, LOOK, &rounds, &width);
        CK(cudaEventRecord(b)); CK(cudaEventSynchronize(b));
        CK(cudaEventElapsedTime(&cpu_ms, a, b));

        /* ---- GPU ---- */
        int *d_nb, *d_nd, *d_ei, *d_ej, *d_ce, *d_cen, *d_any, *d_cnt;
        double *d_th, *d_pub, *d_w, *d_E, *d_t, *d_h, *d_et, *d_eh, *d_evt;
        char *d_elig, *d_sel;
        size_t sN = (size_t)N * sizeof(double), sNT = (size_t)NT * sizeof(double);
        size_t sA = (size_t)N * MAXDEG * sizeof(int);
        CK(cudaMalloc(&d_nb, sA)); CK(cudaMalloc(&d_ce, sA));
        CK(cudaMalloc(&d_ei, (size_t)NE * sizeof(int)));
        CK(cudaMalloc(&d_ej, (size_t)NE * sizeof(int)));
        CK(cudaMalloc(&d_nd, N * sizeof(int)));
        CK(cudaMalloc(&d_cen, N * sizeof(int)));
        CK(cudaMalloc(&d_th, sN)); CK(cudaMalloc(&d_pub, sN));
        CK(cudaMalloc(&d_w, sN)); CK(cudaMalloc(&d_E, sN));
        CK(cudaMalloc(&d_t, sN)); CK(cudaMalloc(&d_h, sN));
        CK(cudaMalloc(&d_et, (size_t)NE * sizeof(double)));
        CK(cudaMalloc(&d_eh, (size_t)NE * sizeof(double)));
        CK(cudaMalloc(&d_evt, sNT));
        CK(cudaMalloc(&d_elig, NT)); CK(cudaMalloc(&d_sel, NT));
        CK(cudaMalloc(&d_any, sizeof(int))); CK(cudaMalloc(&d_cnt, sizeof(int)));
        CK(cudaMemcpy(d_nb, G->nb, sA, cudaMemcpyHostToDevice));
        CK(cudaMemcpy(d_ce, G->ce, sA, cudaMemcpyHostToDevice));
        CK(cudaMemcpy(d_ei, G->ei, NE * sizeof(int), cudaMemcpyHostToDevice));
        CK(cudaMemcpy(d_ej, G->ej, NE * sizeof(int), cudaMemcpyHostToDevice));
        CK(cudaMemcpy(d_nd, G->nd, N * sizeof(int), cudaMemcpyHostToDevice));
        CK(cudaMemcpy(d_cen, G->cen, N * sizeof(int), cudaMemcpyHostToDevice));
        CK(cudaMemcpy(d_th, G->th, sN, cudaMemcpyHostToDevice));
        CK(cudaMemcpy(d_pub, G->pub, sN, cudaMemcpyHostToDevice));
        CK(cudaMemcpy(d_w, G->w, sN, cudaMemcpyHostToDevice));
        CK(cudaMemcpy(d_E, G->E, sN, cudaMemcpyHostToDevice));
        CK(cudaMemcpy(d_t, G->t, sN, cudaMemcpyHostToDevice));
        CK(cudaMemcpy(d_h, G->h, sN, cudaMemcpyHostToDevice));
        CK(cudaMemcpy(d_et, G->et, NE * sizeof(double), cudaMemcpyHostToDevice));
        CK(cudaMemcpy(d_eh, G->eh, NE * sizeof(double), cudaMemcpyHostToDevice));

        int TPB = 256, BL = (NT + TPB - 1) / TPB, BLN = (N + TPB - 1) / TPB;
        long grounds = 0; double gacc = 0;
        CK(cudaEventRecord(a));
        for (;;) {
            int any = 0, cnt = 0;
            CK(cudaMemcpy(d_any, &any, sizeof(int), cudaMemcpyHostToDevice));
            k_elig<<<BL, TPB>>>(N, NE, T, LOOK, d_nb, d_nd, d_ei, d_ej,
                                d_t, d_h, d_et, d_eh, d_elig, d_evt, d_any);
            CK(cudaMemcpy(&any, d_any, sizeof(int), cudaMemcpyDeviceToHost));
            if (!any) break;
            CK(cudaMemcpy(d_cnt, &cnt, sizeof(int), cudaMemcpyHostToDevice));
            k_select<<<BL, TPB>>>(N, NE, d_nb, d_nd, d_ei, d_ej, d_ce, d_cen,
                                  d_elig, d_evt, d_sel, d_cnt);
            CK(cudaMemcpy(&cnt, d_cnt, sizeof(int), cudaMemcpyDeviceToHost));
            if (!cnt) break;
            k_advance<<<BL, TPB>>>(N, NE, Kc, gE, d_nb, d_nd, d_ei, d_ej,
                                   d_sel, d_th, d_pub, d_w, d_E, d_t, d_h,
                                   d_et, d_eh);
            k_publish<<<BLN, TPB>>>(N, d_sel, d_th, d_pub);
            gacc += cnt; grounds++;
        }
        CK(cudaEventRecord(b)); CK(cudaEventSynchronize(b));
        CK(cudaEventElapsedTime(&gpu_ms, a, b));

        double *gth = (double *)malloc(sN), *gE_ = (double *)malloc(sN);
        CK(cudaMemcpy(gth, d_th, sN, cudaMemcpyDeviceToHost));
        CK(cudaMemcpy(gE_, d_E, sN, cudaMemcpyDeviceToHost));
        double md = 0;
        for (int i = 0; i < N; i++) {
            double d1 = fabs(gth[i] - C->th[i]);
            double d2 = fabs(gE_[i] - C->E[i]);
            if (d1 > md) md = d1;
            if (d2 > md) md = d2;
        }
        printf("%7d  %7d  %7ld  %10.1f  %7.2f  %7.2f  %6.1fx  %s\n",
               N, NT, rounds, width, cpu_ms / 1000.0, gpu_ms / 1000.0,
               cpu_ms / gpu_ms, md == 0.0 ? "EXACT" : "differs");
        fflush(stdout);

        free(gth); free(gE_);
        cudaFree(d_nb); cudaFree(d_ce); cudaFree(d_ei); cudaFree(d_ej);
        cudaFree(d_nd); cudaFree(d_cen); cudaFree(d_th); cudaFree(d_pub);
        cudaFree(d_w); cudaFree(d_E); cudaFree(d_t); cudaFree(d_h);
        cudaFree(d_et); cudaFree(d_eh); cudaFree(d_evt);
        cudaFree(d_elig); cudaFree(d_sel); cudaFree(d_any); cudaFree(d_cnt);
        cudaEventDestroy(a); cudaEventDestroy(b);
        fab_free(C); fab_free(G);
    }
    printf("\nreading: 'agree' must read EXACT — the batch is a conflict-free\n");
    printf("set, so GPU and CPU execute the same schedule and must produce\n");
    printf("identical bits. Anything else means the conflict rule is wrong.\n");
    printf("Each round costs 2 device syncs (any/cnt readback), so speedup is\n");
    printf("bounded by batch width vs round latency: the win arrives only\n");
    printf("once the width is large.\n");
    return 0;
}
