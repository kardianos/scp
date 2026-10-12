/* c3pair.c — v103 lane c, step c3: two defects in the medium of bare quanta
 * that holds three dimensions (../PREREG.md).
 *
 * The medium is g1's, unchanged (../../../gate/g1/cas/geo3d.c; its code is
 * carried here verbatim where marked): N quanta and their contacts, six per
 * quantum, all paired, two quanta touching at most once.  No coordinates, no
 * cells, no slots.  "Where" is counted through the medium: the distance
 * between two quanta is the number of contacts on the shortest path.  The
 * medium re-forms by local flips (along i-j-k-l, ij and kl become ik and jl),
 * sampled at temperature t.  Candidates (all 8 combinations), 2t per unit
 * deviation squared:
 *   a  flat squares:     every contact lies in exactly four squares
 *   b  octahedral link:  at every quantum, each contact shares a square with
 *                        exactly four of the other five
 *   c  no triangles:     2t per triangle
 *
 * c3 makes two defects, gap with gap (GG), excess with excess (XX) and gap
 * with excess (GX), at a separation R counted through the medium, and reads
 * them against the same defects made alone in the same medium:
 *   - the interaction energy  E_int = E_pair - E_1 - E_2 + E_0;
 *   - the shared "where":     d_pair - d_1 - d_2 + d_0, over every pair of
 *     quanta present in all four media;
 *   - what the pair becomes:  two readable defects, one, or a healed medium.
 * Each defect is made with its own random stream, and every candidate list is
 * enumerated in the quanta's own identity order, so a defect made beside
 * another and the same defect made alone are the same closure whenever they
 * are out of each other's reach.  After making, the medium is quenched
 * (the strictly downhill flip of lowest energy change, ties by the quanta's
 * identities, until none is left), then sampled at t.
 *
 * Quantum identities (oid) are bookkeeping of quanta, never places.  The
 * dynamics reads them only to order candidate lists and to break ties.
 *
 * build: gcc -O3 -march=native -Wall -Wextra -o c3pair c3pair.c -lm
 * usage: ./c3pair mode=run L=12 a=1 b=0 c=1 k=0   (one registered process)
 *        ./c3pair mode=t0                         (tool gates, PREREG §6)
 */
#include <limits.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define DEG 6
static int OPT_A, OPT_B, OPT_C;
static const double TREG = 1.0; /* the registered temperature, t = 1 */

/* the registered schedule (PREREG §3) */
static const int M_SWEEPS = 500, MEAS_SWEEPS = 400, SAMPLE_EVERY = 100;
static const int RLIST[6] = {1, 2, 3, 4, 6, 8}; /* separations, those with 2R <= L */
#define QR 3 /* a strictly downhill flip has its middle quantum within 3 contacts of a mismatch (PREREG §3) */

/* stream roles (PREREG §3) */
enum { RO_M = 1, RO_P1 = 2, RO_P2 = 3, RO_G1 = 4, RO_X1 = 5, RO_G2 = 6, RO_X2 = 7, RO_W = 8 };
static const char *KIND[7] = {"G1", "X1", "G2", "X2", "GG", "XX", "GX"};

/* ---------------- rng: xoshiro256** seeded by splitmix64 (g1, verbatim) ---------------- */
typedef struct {
    uint64_t s[4];
} Rng;
static uint64_t splitmix64(uint64_t *x)
{
    uint64_t z = (*x += 0x9e3779b97f4a7c15ULL);
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
    return z ^ (z >> 31);
}
static void rng_seed(Rng *r, uint64_t seed)
{
    uint64_t x = seed;
    for (int i = 0; i < 4; i++) r->s[i] = splitmix64(&x);
}
static inline uint64_t rotl(uint64_t x, int k) { return (x << k) | (x >> (64 - k)); }
static inline uint64_t rnext(Rng *r)
{
    uint64_t *s = r->s, res = rotl(s[1] * 5, 7) * 9, t = s[1] << 17;
    s[2] ^= s[0];
    s[3] ^= s[1];
    s[1] ^= s[2];
    s[0] ^= s[3];
    s[2] ^= t;
    s[3] = rotl(s[3], 45);
    return res;
}
static inline double runif(Rng *r) { return (double)(rnext(r) >> 11) * 0x1.0p-53; }
static inline int rint_(Rng *r, int n) { return (int)(((rnext(r) >> 32) * (uint64_t)n) >> 32); }
static uint64_t stream_seed(uint64_t base, int k, int R, int role)
{
    return base * 1000003ULL + 100000ULL * (uint64_t)k + 100ULL * (uint64_t)R + (uint64_t)role;
}

/* ---------------- the medium: quanta and contacts (g1, plus identities) ---------------- */
typedef struct {
    int n, cap, W;
    int *deg, *nb; /* nb[i*DEG + k] */
    int *oid;      /* each quantum's own identity, carried through relabelling */
    uint64_t *adj; /* cap x W bits */
} G;

static G *g_new(int cap)
{
    G *g = calloc(1, sizeof *g);
    g->cap = cap;
    g->W = (cap + 63) / 64;
    g->deg = calloc((size_t)cap, sizeof(int));
    g->nb = calloc((size_t)cap * DEG, sizeof(int));
    g->oid = calloc((size_t)cap, sizeof(int));
    g->adj = calloc((size_t)cap * g->W, sizeof(uint64_t));
    if (!g->deg || !g->nb || !g->oid || !g->adj) {
        fprintf(stderr, "out of memory\n");
        exit(2);
    }
    return g;
}
static void g_free(G *g)
{
    free(g->deg);
    free(g->nb);
    free(g->oid);
    free(g->adj);
    free(g);
}
static void g_reset(G *g, int n)
{
    g->n = n;
    memset(g->deg, 0, sizeof(int) * (size_t)g->cap);
    memset(g->adj, 0, sizeof(uint64_t) * (size_t)g->cap * g->W);
}
static void g_copy(G *d, const G *s)
{
    if (d->cap != s->cap) {
        fprintf(stderr, "g_copy: cap mismatch\n");
        exit(2);
    }
    d->n = s->n;
    memcpy(d->deg, s->deg, sizeof(int) * (size_t)s->cap);
    memcpy(d->nb, s->nb, sizeof(int) * (size_t)s->cap * DEG);
    memcpy(d->oid, s->oid, sizeof(int) * (size_t)s->cap);
    memcpy(d->adj, s->adj, sizeof(uint64_t) * (size_t)s->cap * s->W);
}
static inline int ADJ(const G *g, int i, int j) { return (int)((g->adj[(size_t)i * g->W + (j >> 6)] >> (j & 63)) & 1u); }
static inline void setb(G *g, int i, int j) { g->adj[(size_t)i * g->W + (j >> 6)] |= 1ULL << (j & 63); }
static inline void clrb(G *g, int i, int j) { g->adj[(size_t)i * g->W + (j >> 6)] &= ~(1ULL << (j & 63)); }
static void add_edge(G *g, int i, int j)
{
    if (g->deg[i] >= DEG || g->deg[j] >= DEG || i == j || ADJ(g, i, j)) {
        fprintf(stderr, "add_edge: invalid %d-%d\n", i, j);
        exit(2);
    }
    g->nb[i * DEG + g->deg[i]++] = j;
    g->nb[j * DEG + g->deg[j]++] = i;
    setb(g, i, j);
    setb(g, j, i);
}
static void rm_nb(G *g, int i, int j)
{
    int *p = g->nb + i * DEG;
    for (int k = 0; k < g->deg[i]; k++)
        if (p[k] == j) {
            p[k] = p[--g->deg[i]];
            return;
        }
    fprintf(stderr, "rm_nb: %d not a contact of %d\n", j, i);
    exit(2);
}
static void del_edge(G *g, int i, int j)
{
    rm_nb(g, i, j);
    rm_nb(g, j, i);
    clrb(g, i, j);
    clrb(g, j, i);
}
/* drop quantum q (no contacts left): the last quantum's bookkeeping index becomes q; its identity goes with it */
static void g_drop(G *g, int q)
{
    int last = g->n - 1;
    if (g->deg[q]) {
        fprintf(stderr, "g_drop: quantum still has contacts\n");
        exit(2);
    }
    if (q != last) {
        for (int k = 0; k < g->deg[last]; k++) {
            int u = g->nb[last * DEG + k];
            int *p = g->nb + u * DEG;
            for (int m = 0; m < g->deg[u]; m++)
                if (p[m] == last) p[m] = q;
            clrb(g, u, last);
            setb(g, u, q);
            g->nb[q * DEG + k] = u;
            setb(g, q, u);
        }
        g->deg[q] = g->deg[last];
        g->deg[last] = 0;
        g->oid[q] = g->oid[last];
        memset(g->adj + (size_t)last * g->W, 0, sizeof(uint64_t) * g->W);
    }
    g->n--;
}
static int idx_of(const G *g, int o)
{
    for (int i = 0; i < g->n; i++)
        if (g->oid[i] == o) return i;
    return -1;
}

/* the flat start (bookkeeping coordinates, dropped at once); rev = the same contacts entered in reverse order (T0.2) */
static void build_cubic_order(G *g, int L, int rev)
{
    int N = L * L * L;
    g_reset(g, N);
    for (int i = 0; i < N; i++) g->oid[i] = i;
    for (int a = 0; a < N; a++) {
        int i = rev ? N - 1 - a : a, x = i % L, y = (i / L) % L, z = i / (L * L);
        int nbr[3] = {(x + 1) % L + L * (y + L * z), x + L * ((y + 1) % L + L * z), x + L * (y + L * ((z + 1) % L))};
        for (int m = 0; m < 3; m++) add_edge(g, i, nbr[rev ? 2 - m : m]);
    }
}

/* ---------------- local structure (g1, verbatim) ---------------- */
static int sq(const G *g, int u, int v) /* squares through contact u-v */
{
    int c = 0;
    const int *Nu = g->nb + u * DEG, *Nv = g->nb + v * DEG;
    for (int a = 0; a < g->deg[u]; a++) {
        int w = Nu[a];
        if (w == v) continue;
        for (int b = 0; b < g->deg[v]; b++) {
            int x = Nv[b];
            if (x == u || x == w) continue;
            if (ADJ(g, w, x)) c++;
        }
    }
    return c;
}
static int lk(const G *g, int i, int j) /* contacts k of i (k != j) that share a square with i and j */
{
    int c = 0;
    const int *Ni = g->nb + i * DEG, *Nj = g->nb + j * DEG;
    for (int a = 0; a < g->deg[i]; a++) {
        int k = Ni[a];
        if (k == j) continue;
        for (int b = 0; b < g->deg[j]; b++) {
            int m = Nj[b];
            if (m == i || m == k) continue;
            if (ADJ(g, m, k)) {
                c++;
                break;
            }
        }
    }
    return c;
}
static int cn(const G *g, int u, int v) /* common contacts of u and v (triangles through u-v) */
{
    int c = 0;
    const int *Nu = g->nb + u * DEG;
    for (int a = 0; a < g->deg[u]; a++)
        if (Nu[a] != v && ADJ(g, Nu[a], v)) c++;
    return c;
}

/* ---------------- energy in units of 2t (g1, verbatim) ---------------- */
static char *MARK;
static int *VL, SCAP;
static void scratch(int cap)
{
    if (cap <= SCAP) return;
    free(MARK);
    free(VL);
    MARK = calloc((size_t)cap, 1);
    VL = malloc(sizeof(int) * (size_t)cap);
    SCAP = cap;
}
static inline int anyopt(void) { return OPT_A || OPT_B || OPT_C; }
static long eloc_impl(const G *g, const int *V, int nv, int all)
{
    long e = 0;
    for (int a = 0; a < nv; a++) {
        int u = all ? a : V[a];
        const int *Nu = g->nb + u * DEG;
        for (int k = 0; k < g->deg[u]; k++) {
            int v = Nu[k];
            if (v <= u || !(all || MARK[v])) continue;
            if (OPT_A) {
                int s = sq(g, u, v) - 4;
                e += (long)s * s;
            }
            if (OPT_C)
                for (int m = 0; m < g->deg[u]; m++) {
                    int w = Nu[m];
                    if (w > v && (all || MARK[w]) && ADJ(g, v, w)) e++;
                }
        }
        if (OPT_B)
            for (int k = 0; k < g->deg[u]; k++) {
                int q = lk(g, u, Nu[k]) - 4;
                e += (long)q * q;
            }
    }
    return e;
}
static long eloc(const G *g, const int *V, int nv) { return eloc_impl(g, V, nv, 0); }
static long efull(const G *g) { return eloc_impl(g, NULL, g->n, 1); }
static int build_V1(const G *g, const int *V0, int n0)
{
    int nv = 0;
    for (int a = 0; a < n0; a++)
        if (!MARK[V0[a]]) {
            MARK[V0[a]] = 1;
            VL[nv++] = V0[a];
        }
    for (int a = 0; a < n0; a++)
        for (int k = 0; k < g->deg[V0[a]]; k++) {
            int y = g->nb[V0[a] * DEG + k];
            if (!MARK[y]) {
                MARK[y] = 1;
                VL[nv++] = y;
            }
        }
    return nv;
}
static void unmark(int nv)
{
    for (int a = 0; a < nv; a++) MARK[VL[a]] = 0;
}

/* ---------------- the re-forming: local flips (g1, verbatim) ---------------- */
typedef struct {
    int i, j, k, l;
} Flip;
static long NPROP, NACC;
static int pick_other(const G *g, Rng *r, int j, int excl)
{
    int a = rint_(r, g->deg[j] - 1), c = 0;
    const int *Nj = g->nb + j * DEG;
    for (int m = 0; m < g->deg[j]; m++) {
        if (Nj[m] == excl) continue;
        if (c++ == a) return Nj[m];
    }
    return -1;
}
static int propose(const G *g, Rng *r, Flip *f)
{
    int j = rint_(r, g->n), k = g->nb[j * DEG + rint_(r, g->deg[j])];
    int i = pick_other(g, r, j, k), l = pick_other(g, r, k, j);
    if (i < 0 || l < 0 || i == l || ADJ(g, i, k) || ADJ(g, j, l)) return 0;
    f->i = i;
    f->j = j;
    f->k = k;
    f->l = l;
    return 1;
}
static void flip_apply(G *g, const Flip *f)
{
    del_edge(g, f->i, f->j);
    del_edge(g, f->k, f->l);
    add_edge(g, f->i, f->k);
    add_edge(g, f->j, f->l);
}
static void flip_revert(G *g, const Flip *f)
{
    del_edge(g, f->i, f->k);
    del_edge(g, f->j, f->l);
    add_edge(g, f->i, f->j);
    add_edge(g, f->k, f->l);
}
static long flip_de(G *g, const Flip *f)
{
    int V0[4] = {f->i, f->j, f->k, f->l};
    int nv = build_V1(g, V0, 4);
    long e0 = eloc(g, VL, nv);
    flip_apply(g, f);
    long e1 = eloc(g, VL, nv);
    flip_revert(g, f);
    unmark(nv);
    return e1 - e0;
}
static int mc_step(G *g, Rng *r, double T, int free_)
{
    Flip f;
    NPROP++;
    if (!propose(g, r, &f)) return 0;
    if (free_ || !anyopt()) {
        flip_apply(g, &f);
        NACC++;
        return 1;
    }
    int V0[4] = {f.i, f.j, f.k, f.l};
    int nv = build_V1(g, V0, 4);
    long e0 = eloc(g, VL, nv);
    flip_apply(g, &f);
    long e1 = eloc(g, VL, nv);
    double dE = 2.0 * (double)(e1 - e0); /* in t */
    int acc = dE <= 0 || runif(r) < exp(-dE / T);
    if (acc)
        NACC++;
    else
        flip_revert(g, &f);
    unmark(nv);
    return acc;
}
static void sweeps(G *g, Rng *r, int ns, double T, int free_)
{
    for (int s = 0; s < ns; s++) {
        long m = 3L * g->n; /* one sweep = one proposal per contact */
        for (long p = 0; p < m; p++) mc_step(g, r, T, free_);
    }
}

/* ---------------- graph instruments ---------------- */
static int bfs(const G *g, const int *src, int ns, int *dist, int *que)
{
    for (int i = 0; i < g->n; i++) dist[i] = -1;
    int h = 0, t = 0;
    for (int a = 0; a < ns; a++)
        if (dist[src[a]] < 0) {
            dist[src[a]] = 0;
            que[t++] = src[a];
        }
    while (h < t) {
        int u = que[h++];
        for (int k = 0; k < g->deg[u]; k++) {
            int v = g->nb[u * DEG + k];
            if (dist[v] < 0) {
                dist[v] = dist[u] + 1;
                que[t++] = v;
            }
        }
    }
    return t;
}
static int pieces(const G *g)
{
    int *dist = malloc(sizeof(int) * (size_t)g->n), *que = malloc(sizeof(int) * (size_t)g->n), *seen = calloc((size_t)g->n, sizeof(int));
    int np = 0;
    for (int s = 0; s < g->n; s++) {
        if (seen[s]) continue;
        np++;
        int nr = bfs(g, &s, 1, dist, que);
        for (int a = 0; a < nr; a++) seen[que[a]] = 1;
    }
    free(dist);
    free(que);
    free(seen);
    return np;
}
/* g1's invariants, plus: every identity distinct */
static int check(const G *g, int nexp)
{
    int ok = g->n == nexp;
    for (int i = 0; i < g->n && ok; i++) {
        if (g->deg[i] != DEG) ok = 0;
        int pc = 0;
        for (int w = 0; w < g->W; w++) pc += __builtin_popcountll(g->adj[(size_t)i * g->W + w]);
        if (pc != g->deg[i]) ok = 0;
        for (int k = 0; k < g->deg[i]; k++) {
            int j = g->nb[i * DEG + k];
            if (j == i || j < 0 || j >= g->n || !ADJ(g, i, j) || !ADJ(g, j, i)) ok = 0;
            for (int m = k + 1; m < g->deg[i]; m++)
                if (g->nb[i * DEG + m] == j) ok = 0;
        }
    }
    if (ok) {
        char *seen = calloc((size_t)g->cap, 1);
        for (int i = 0; i < g->n && ok; i++) {
            int o = g->oid[i];
            if (o < 0 || o >= g->cap || seen[o]) ok = 0;
            else seen[o] = 1;
        }
        free(seen);
    }
    if (ok && pieces(g) != 1) ok = 0;
    return ok;
}
static int uf_find(int *p, int x)
{
    while (p[x] != x) {
        p[x] = p[p[x]];
        x = p[x];
    }
    return x;
}
static void uf_union(int *p, int a, int b)
{
    a = uf_find(p, a);
    b = uf_find(p, b);
    if (a != b) p[a] = b;
}

/* ---------------- readings (PREREG §4) ---------------- */
/* a contact is mismatched if s != 4, or it is in a triangle, or q != 4 at either end (g1).
 * isD[i] = 1 if quantum i touches one. Returns the number of mismatched contacts; *E = the energy in t. */
static int mismatch(const G *g, char *isD, double *E)
{
    long A2 = 0, B2 = 0, T3 = 0;
    int nmis = 0;
    memset(isD, 0, (size_t)g->n);
    for (int u = 0; u < g->n; u++)
        for (int a = 0; a < g->deg[u]; a++) {
            int v = g->nb[u * DEG + a];
            if (v <= u) continue;
            int s = sq(g, u, v), t = cn(g, u, v), l1 = lk(g, u, v), l2 = lk(g, v, u);
            A2 += (long)(s - 4) * (s - 4);
            B2 += (long)(l1 - 4) * (l1 - 4) + (long)(l2 - 4) * (l2 - 4);
            T3 += t;
            if (s != 4 || t > 0 || l1 != 4 || l2 != 4) {
                nmis++;
                isD[u] = isD[v] = 1;
            }
        }
    T3 /= 3;
    if (E) *E = 2.0 * (double)((OPT_A ? A2 : 0) + (OPT_B ? B2 : 0) + (OPT_C ? T3 : 0));
    return nmis;
}
typedef struct {
    double E, meand;
    int n, nmis, nD, ncl, dspread, diam, pieces;
} Reading;
static int NO; /* the original quanta, L^3: the distance matrices cover these */
static int NQ; /* L^3 as well: the excesses take identities NQ (first) and NQ + 1 (second) */
/* every reading of one state. If M, the distances between original quanta are written to M[o1*NO + o2]
 * (255 = not present). If cl, cl[i] = the cluster of a mismatched quantum, else -1. */
static void readm(const G *g, Reading *rd, uint8_t *M, int *cl)
{
    int n = g->n;
    char *isD = malloc((size_t)n);
    int *dist = malloc(sizeof(int) * (size_t)n), *que = malloc(sizeof(int) * (size_t)n), *uf = malloc(sizeof(int) * (size_t)n);
    rd->n = n;
    rd->nmis = mismatch(g, isD, &rd->E);
    for (int i = 0; i < n; i++) uf[i] = i;
    if (M) memset(M, 255, (size_t)NO * NO);
    double tot = 0;
    int diam = 0, dspread = -1, reached0 = 0;
    for (int s = 0; s < n; s++) {
        int nr = bfs(g, &s, 1, dist, que);
        if (s == 0) reached0 = nr;
        double sm = 0;
        uint8_t *row = (M && g->oid[s] < NO) ? M + (size_t)g->oid[s] * NO : NULL;
        for (int a = 0; a < nr; a++) {
            int y = que[a], d = dist[y];
            sm += d;
            if (d > diam) diam = d;
            if (row && g->oid[y] < NO) row[g->oid[y]] = (uint8_t)(d < 254 ? d : 254);
            if (isD[s] && isD[y]) {
                if (d > dspread) dspread = d;
                if (d <= 2) uf_union(uf, s, y);
            }
        }
        tot += nr > 1 ? sm / (nr - 1) : 0;
    }
    int nD = 0, ncl = 0;
    for (int i = 0; i < n; i++)
        if (isD[i]) {
            nD++;
            if (uf_find(uf, i) == i) ncl++;
        }
    if (cl)
        for (int i = 0; i < n; i++) cl[i] = isD[i] ? uf_find(uf, i) : -1;
    rd->meand = tot / n;
    rd->nD = nD;
    rd->ncl = ncl;
    rd->dspread = dspread;
    rd->diam = diam;
    rd->pieces = reached0 == n ? 1 : pieces(g);
    free(isD);
    free(dist);
    free(que);
    free(uf);
}
/* the medium distance between the centres of exactly two mismatch clusters (-1 otherwise). A cluster's centre
 * is the set of quanta whose largest distance to the cluster's members is smallest (read from the state). */
static int centre_dist(const G *g, const int *cl)
{
    int n = g->n, r1 = -1, r2 = -1;
    for (int i = 0; i < n; i++) {
        int c = cl[i];
        if (c < 0) continue;
        if (r1 < 0) r1 = c;
        else if (c != r1) {
            if (r2 < 0) r2 = c;
            else if (c != r2) return -1;
        }
    }
    if (r2 < 0) return -1;
    int *dist = malloc(sizeof(int) * (size_t)n), *que = malloc(sizeof(int) * (size_t)n);
    int *e1 = calloc((size_t)n, sizeof(int)), *e2 = calloc((size_t)n, sizeof(int)), *C = malloc(sizeof(int) * (size_t)n);
    for (int y = 0; y < n; y++) {
        if (cl[y] < 0) continue;
        int *e = cl[y] == r1 ? e1 : e2;
        bfs(g, &y, 1, dist, que);
        for (int x = 0; x < n; x++)
            if (dist[x] > e[x]) e[x] = dist[x];
    }
    int m1 = INT_MAX, m2 = INT_MAX, nc = 0, dcc = INT_MAX;
    for (int x = 0; x < n; x++) {
        if (e1[x] < m1) m1 = e1[x];
        if (e2[x] < m2) m2 = e2[x];
    }
    for (int x = 0; x < n; x++)
        if (e1[x] == m1) C[nc++] = x;
    bfs(g, C, nc, dist, que);
    for (int x = 0; x < n; x++)
        if (e2[x] == m2 && dist[x] >= 0 && dist[x] < dcc) dcc = dist[x];
    free(dist);
    free(que);
    free(e1);
    free(e2);
    free(C);
    return dcc == INT_MAX ? -1 : dcc;
}
/* the contacts as sorted identity pairs */
static int cmp_u64(const void *a, const void *b)
{
    uint64_t x = *(const uint64_t *)a, y = *(const uint64_t *)b;
    return (x > y) - (x < y);
}
static int clist(const G *g, uint64_t *out)
{
    int c = 0;
    for (int u = 0; u < g->n; u++)
        for (int k = 0; k < g->deg[u]; k++) {
            int v = g->nb[u * DEG + k];
            if (v <= u) continue;
            uint64_t a = (uint64_t)g->oid[u], b = (uint64_t)g->oid[v];
            out[c++] = a < b ? (a << 32 | b) : (b << 32 | a);
        }
    qsort(out, (size_t)c, sizeof(uint64_t), cmp_u64);
    return c;
}
/* the entries of b that are not in a (both sorted); written to out if given */
static int cdiff(const uint64_t *a, int na, const uint64_t *b, int nb, uint64_t *out)
{
    int i = 0, j = 0, c = 0;
    while (j < nb) {
        if (i < na && a[i] < b[j])
            i++;
        else if (i < na && a[i] == b[j]) {
            i++;
            j++;
        } else {
            if (out) out[c] = b[j];
            c++;
            j++;
        }
    }
    return c;
}

/* ---------------- making a gap and an excess (PREREG §3) ---------------- */
static int PM[15][6];
static void init_pm(void) /* the 15 perfect matchings of six (g1, verbatim) */
{
    int c = 0;
    for (int x = 1; x < 6; x++) {
        int rest[4], nr = 0;
        for (int y = 1; y < 6; y++)
            if (y != x) rest[nr++] = y;
        for (int z = 1; z < 4; z++) {
            int o[2], no = 0;
            for (int w = 1; w < 4; w++)
                if (w != z) o[no++] = rest[w];
            int *P = PM[c++];
            P[0] = 0;
            P[1] = x;
            P[2] = rest[0];
            P[3] = rest[z];
            P[4] = o[0];
            P[5] = o[1];
        }
    }
}
static int gibbs_pick(Rng *r, const long *dE, const int *valid, int n, double T, double *wsel) /* g1, verbatim */
{
    long emin = 0;
    int first = 1;
    for (int p = 0; p < n; p++)
        if (!valid || valid[p])
            if (first || dE[p] < emin) {
                emin = dE[p];
                first = 0;
            }
    double Z = 0;
    for (int p = 0; p < n; p++)
        if (!valid || valid[p]) Z += exp(-2.0 * (double)(dE[p] - emin) / T);
    double u = runif(r) * Z, c = 0;
    int pick = -1;
    for (int p = 0; p < n; p++)
        if (!valid || valid[p]) {
            pick = p;
            c += exp(-2.0 * (double)(dE[p] - emin) / T);
            if (u < c) break;
        }
    *wsel = exp(-2.0 * (double)(dE[pick] - emin) / T) / Z;
    return pick;
}
static const G *SORT_G;
static int cmp_oid(const void *a, const void *b)
{
    int x = SORT_G->oid[*(const int *)a], y = SORT_G->oid[*(const int *)b];
    return (x > y) - (x < y);
}
/* q's six contacts in identity order (J); for each of the 15 re-pairings, whether it is valid (no doubled contact)
 * and its energy change (2t units). If dfull, the full recomputation too (T0). Leaves the medium unchanged. */
static int gap_table(G *g, int q, int *J, int *valid, long *dl, long *dfull)
{
    memcpy(J, g->nb + q * DEG, sizeof(int) * DEG);
    SORT_G = g;
    qsort(J, DEG, sizeof(int), cmp_oid);
    int V0[DEG + 1];
    V0[0] = q;
    for (int m = 0; m < DEG; m++) V0[m + 1] = J[m];
    int nv = build_V1(g, V0, DEG + 1);
    long e0 = eloc(g, VL, nv), E0 = dfull ? efull(g) : 0;
    int nvalid = 0;
    for (int p = 0; p < 15; p++) {
        const int *P = PM[p];
        dl[p] = 0;
        valid[p] = !ADJ(g, J[P[0]], J[P[1]]) && !ADJ(g, J[P[2]], J[P[3]]) && !ADJ(g, J[P[4]], J[P[5]]);
        if (!valid[p]) continue;
        nvalid++;
        for (int m = 0; m < DEG; m++) del_edge(g, q, J[m]);
        for (int m = 0; m < 6; m += 2) add_edge(g, J[P[m]], J[P[m + 1]]);
        dl[p] = eloc(g, VL, nv) - e0;
        if (dfull) dfull[p] = efull(g) - E0;
        for (int m = 0; m < 6; m += 2) del_edge(g, J[P[m]], J[P[m + 1]]);
        for (int m = 0; m < DEG; m++) add_edge(g, q, J[m]);
    }
    unmark(nv);
    return nvalid;
}
/* remove quantum q: the medium closes over it at once, its six contacts re-paired among its neighbours by Gibbs
 * weight at t. Returns the energy change in t (NAN if no pairing is valid). pair_o: the identities of the three
 * new contacts' ends. */
static double make_gap_at(G *g, int q, Rng *r, int *pair_o)
{
    int J[DEG], valid[15];
    long dl[15];
    if (q < 0 || g->deg[q] != DEG || !gap_table(g, q, J, valid, dl, NULL)) return NAN;
    double w;
    int p = gibbs_pick(r, dl, valid, 15, TREG, &w);
    for (int m = 0; m < DEG; m++) del_edge(g, q, J[m]);
    for (int m = 0; m < 6; m += 2) add_edge(g, J[PM[p][m]], J[PM[p][m + 1]]);
    if (pair_o)
        for (int m = 0; m < 6; m++) pair_o[m] = g->oid[J[PM[p][m]]];
    g_drop(g, q);
    return 2.0 * (double)dl[p];
}
typedef struct {
    int u, v, k1, k2;
} Edge;
static int cmp_edge(const void *a, const void *b)
{
    const Edge *x = a, *y = b;
    if (x->k1 != y->k1) return (x->k1 > y->k1) - (x->k1 < y->k1);
    return (x->k2 > y->k2) - (x->k2 < y->k2);
}
typedef struct {
    int e[3];
} Tri;
/* the contacts with both ends within distance 2 of p, in identity order, and every triple of pairwise disjoint
 * ones with the energy change (2t units) of a new quantum taking all three. If df_o, the full recomputation (T0). */
static int excess_table(G *g, int p, Edge **E_o, Tri **tr_o, long **dl_o, long **df_o)
{
    int n = g->n, *dist = malloc(sizeof(int) * (size_t)n), *que = malloc(sizeof(int) * (size_t)n);
    bfs(g, &p, 1, dist, que);
    int ne = 0;
    Edge *E = malloc(sizeof(Edge) * 3 * (size_t)n);
    for (int u = 0; u < n; u++)
        if (dist[u] >= 0 && dist[u] <= 2)
            for (int k = 0; k < g->deg[u]; k++) {
                int v = g->nb[u * DEG + k];
                if (v > u && dist[v] >= 0 && dist[v] <= 2) {
                    int a = g->oid[u], b = g->oid[v];
                    E[ne].u = u;
                    E[ne].v = v;
                    E[ne].k1 = a < b ? a : b;
                    E[ne].k2 = a < b ? b : a;
                    ne++;
                }
            }
    qsort(E, (size_t)ne, sizeof(Edge), cmp_edge);
    int cap = 1024, nt = 0;
    Tri *tr = malloc(sizeof(Tri) * cap);
    for (int a = 0; a < ne; a++)
        for (int b = a + 1; b < ne; b++) {
            if (E[a].u == E[b].u || E[a].u == E[b].v || E[a].v == E[b].u || E[a].v == E[b].v) continue;
            for (int c = b + 1; c < ne; c++) {
                int x = E[c].u, y = E[c].v;
                if (x == E[a].u || x == E[a].v || y == E[a].u || y == E[a].v || x == E[b].u || x == E[b].v || y == E[b].u || y == E[b].v) continue;
                if (nt == cap) {
                    cap *= 2;
                    tr = realloc(tr, sizeof(Tri) * cap);
                }
                tr[nt].e[0] = a;
                tr[nt].e[1] = b;
                tr[nt].e[2] = c;
                nt++;
            }
        }
    long *dl = malloc(sizeof(long) * (nt ? nt : 1)), *df = df_o ? malloc(sizeof(long) * (nt ? nt : 1)) : NULL;
    int q = g->n;
    g->n++; /* the new quantum, no contacts yet */
    long E0 = df ? efull(g) : 0;
    for (int t = 0; t < nt; t++) {
        int V0[7];
        V0[0] = q;
        for (int m = 0; m < 3; m++) {
            V0[1 + 2 * m] = E[tr[t].e[m]].u;
            V0[2 + 2 * m] = E[tr[t].e[m]].v;
        }
        int nv = build_V1(g, V0, 7);
        long e0 = eloc(g, VL, nv);
        for (int m = 0; m < 3; m++) {
            int u = E[tr[t].e[m]].u, v = E[tr[t].e[m]].v;
            del_edge(g, u, v);
            add_edge(g, u, q);
            add_edge(g, q, v);
        }
        dl[t] = eloc(g, VL, nv) - e0;
        if (df) df[t] = efull(g) - E0;
        for (int m = 0; m < 3; m++) {
            int u = E[tr[t].e[m]].u, v = E[tr[t].e[m]].v;
            del_edge(g, u, q);
            del_edge(g, q, v);
            add_edge(g, u, v);
        }
        unmark(nv);
    }
    g->n--;
    free(dist);
    free(que);
    *E_o = E;
    *tr_o = tr;
    *dl_o = dl;
    if (df_o) *df_o = df;
    return nt;
}
/* a new quantum (identity newoid) splits the three contacts u[m]-v[m] into u-q-v */
static void insert_excess(G *g, const int *u, const int *v, int newoid)
{
    int q = g->n;
    g->n++;
    g->oid[q] = newoid;
    for (int m = 0; m < 3; m++) {
        del_edge(g, u[m], v[m]);
        add_edge(g, u[m], q);
        add_edge(g, q, v[m]);
    }
}
/* a new quantum takes three contacts near p by Gibbs weight at t. Returns the energy change in t (NAN if none). */
static double make_excess_at(G *g, int p, int newoid, Rng *r)
{
    if (p < 0) return NAN;
    Edge *E;
    Tri *tr;
    long *dl;
    int nt = excess_table(g, p, &E, &tr, &dl, NULL);
    double dE = NAN;
    if (nt) {
        double w;
        int t = gibbs_pick(r, dl, NULL, nt, TREG, &w), u[3], v[3];
        for (int m = 0; m < 3; m++) {
            u[m] = E[tr[t].e[m]].u;
            v[m] = E[tr[t].e[m]].v;
        }
        insert_excess(g, u, v, newoid);
        dE = 2.0 * (double)dl[t];
    }
    free(E);
    free(tr);
    free(dl);
    return dE;
}

/* ---------------- the quench (PREREG §3) ---------------- */
static void flip_key(const G *g, const Flip *f, int *key) /* a flip read from either end is the same flip */
{
    int a[4] = {g->oid[f->i], g->oid[f->j], g->oid[f->k], g->oid[f->l]}, b[4] = {a[3], a[2], a[1], a[0]}, rev = 0;
    for (int m = 0; m < 4; m++)
        if (a[m] != b[m]) {
            rev = b[m] < a[m];
            break;
        }
    memcpy(key, rev ? b : a, sizeof a);
}
static int key_less(const int *x, const int *y)
{
    for (int m = 0; m < 4; m++)
        if (x[m] != y[m]) return x[m] < y[m];
    return 0;
}
/* apply the strictly downhill flip with the lowest (energy change, identity key) until none is left.
 * full = 0: only flips whose middle quantum j is within QR contacts of a mismatched quantum (exact: a flip whose
 * neighbourhood holds no mismatch starts from zero local energy and cannot go down); full = 1: every flip (T0.3).
 * Returns the number of flips; *dsum = the energy change in 2t units. */
static int quench(G *g, int full, long *dsum)
{
    *dsum = 0;
    if (!anyopt()) return 0;
    int n = g->n, steps = 0;
    char *isD = malloc((size_t)n);
    int *dist = malloc(sizeof(int) * (size_t)n), *que = malloc(sizeof(int) * (size_t)n), *src = malloc(sizeof(int) * (size_t)n);
    for (;;) {
        int ns = 0;
        mismatch(g, isD, NULL);
        for (int i = 0; i < n; i++)
            if (isD[i]) src[ns++] = i;
        if (!ns) break;
        if (!full) bfs(g, src, ns, dist, que);
        Flip best = {0, 0, 0, 0};
        long bde = 0;
        int bkey[4] = {0, 0, 0, 0}, have = 0;
        for (int j = 0; j < n; j++) {
            if (!full && (dist[j] < 0 || dist[j] > QR)) continue;
            int Nj[DEG], dj = g->deg[j];
            memcpy(Nj, g->nb + j * DEG, sizeof Nj);
            for (int a = 0; a < dj; a++) {
                int k = Nj[a], Nk[DEG], dk = g->deg[k];
                memcpy(Nk, g->nb + k * DEG, sizeof Nk);
                for (int b = 0; b < dj; b++) {
                    int i = Nj[b];
                    if (i == k) continue;
                    for (int c = 0; c < dk; c++) {
                        int l = Nk[c];
                        if (l == j || i == l || ADJ(g, i, k) || ADJ(g, j, l)) continue;
                        Flip f = {i, j, k, l};
                        long de = flip_de(g, &f);
                        if (de >= 0) continue;
                        int key[4];
                        flip_key(g, &f, key);
                        if (!have || de < bde || (de == bde && key_less(key, bkey))) {
                            best = f;
                            bde = de;
                            memcpy(bkey, key, sizeof bkey);
                            have = 1;
                        }
                    }
                }
            }
        }
        if (!have) break;
        flip_apply(g, &best);
        *dsum += bde;
        steps++;
    }
    free(isD);
    free(dist);
    free(que);
    free(src);
    return steps;
}

/* ---------------- one medium: make, quench, read, sample (PREREG §3-§4) ---------------- */
static uint64_t BASE;
typedef struct {
    int ok, steps, changed, nmis_end, ncl_end, dcc;
    double Emade, Eafter, dEq, Ew, meandw, acc, secs;
    Reading q;
} Med;
/* one defect: a gap (isgap) or an excess at the quantum with identity site; second = the defect at p2 */
static double make_one(G *h, int isgap, int site, int second, int k, int R)
{
    int idx = idx_of(h, site);
    Rng r;
    if (isgap) {
        rng_seed(&r, stream_seed(BASE, k, second ? R : 0, second ? RO_G2 : RO_G1));
        return make_gap_at(h, idx, &r, NULL);
    }
    rng_seed(&r, stream_seed(BASE, k, second ? R : 0, second ? RO_X2 : RO_X1));
    return make_excess_at(h, idx, NQ + second, &r);
}
static const int ISG[7][2] = {{1, -1}, {0, -1}, {-1, 1}, {-1, 0}, {1, 1}, {0, 0}, {1, 0}}; /* per kind: defect at p1, at p2 (1 gap, 0 excess, -1 none) */
/* make the defects of `kind` in a copy of s0 (no quench, no window): T0 and do_medium */
static int make_kind(G *h, const G *s0, int kind, int k, int R, int p1, int p2, double *Em, int *nexp)
{
    g_copy(h, s0);
    *Em = 0;
    *nexp = s0->n;
    for (int d = 0; d < 2; d++) {
        int t = ISG[kind][d];
        if (t < 0) continue;
        double e = make_one(h, t, d ? p2 : p1, d, k, R);
        if (isnan(e)) return 0;
        *Em += e;
        *nexp += t ? -1 : 1;
    }
    return 1;
}
static void do_medium(G *h, const G *s0, int kind, int k, int R, int p1, int p2, uint8_t *M, int *cl, Med *m)
{
    clock_t c0 = clock();
    int nexp;
    memset(m, 0, sizeof *m);
    m->dcc = -1;
    if (!make_kind(h, s0, kind, k, R, p1, p2, &m->Emade, &nexp)) {
        printf("FAIL kind=%s R=%d k=%d what=make\n", KIND[kind], R, k);
        fflush(stdout);
        return;
    }
    m->Eafter = 2.0 * (double)efull(h);
    long dq = 0;
    m->steps = quench(h, 0, &dq);
    m->dEq = 2.0 * (double)dq;
    readm(h, &m->q, M, cl);
    if (cl && m->q.ncl == 2) m->dcc = centre_dist(h, cl);
    printf("MADE kind=%s R=%d k=%d n=%d E_made=%.6g E_after=%.6g quench_steps=%d dE_quench=%.6g\n", KIND[kind], R, k, h->n, m->Emade, m->Eafter, m->steps, m->dEq);
    printf("QUENCH kind=%s R=%d k=%d E=%.6f nmis=%d nD=%d ncl=%d dspread=%d meand=%.9f diam=%d pieces=%d dcc=%d\n", KIND[kind], R, k, m->q.E, m->q.nmis, m->q.nD, m->q.ncl,
           m->q.dspread, m->q.meand, m->q.diam, m->q.pieces, m->dcc);
    fflush(stdout);
    uint64_t *L0 = malloc(sizeof(uint64_t) * 3 * (size_t)h->cap), *L1 = malloc(sizeof(uint64_t) * 3 * (size_t)h->cap);
    int n0 = clist(h, L0), ns = MEAS_SWEEPS / SAMPLE_EVERY;
    Rng r;
    rng_seed(&r, stream_seed(BASE, k, kind >= 2 ? R : 0, RO_W + kind));
    double sE = 0, sd = 0;
    long np = 0, na = 0;
    Reading rd = m->q;
    for (int s = 1; s <= ns; s++) {
        NPROP = NACC = 0;
        sweeps(h, &r, SAMPLE_EVERY, TREG, 0);
        np += NPROP;
        na += NACC;
        readm(h, &rd, NULL, NULL);
        printf("SAMPLE kind=%s R=%d k=%d s=%d E=%.6f nmis=%d nD=%d ncl=%d dspread=%d meand=%.9f diam=%d pieces=%d\n", KIND[kind], R, k, s, rd.E, rd.nmis, rd.nD, rd.ncl, rd.dspread,
               rd.meand, rd.diam, rd.pieces);
        sE += rd.E;
        sd += rd.meand;
    }
    int n1 = clist(h, L1);
    m->changed = cdiff(L0, n0, L1, n1, NULL);
    m->Ew = sE / ns;
    m->meandw = sd / ns;
    m->acc = np ? (double)na / np : 0.0;
    m->nmis_end = rd.nmis;
    m->ncl_end = rd.ncl;
    m->secs = (double)(clock() - c0) / CLOCKS_PER_SEC;
    m->ok = 1;
    printf("WINDOW kind=%s R=%d k=%d Ew=%.6f meandw=%.9f acc=%.8g changed=%d nmis_end=%d ncl_end=%d secs=%.1f\n", KIND[kind], R, k, m->Ew, m->meandw, m->acc, m->changed, m->nmis_end,
           m->ncl_end, m->secs);
    printf("CHECK kind=%s R=%d k=%d ok=%d\n", KIND[kind], R, k, check(h, nexp));
    fflush(stdout);
    free(L0);
    free(L1);
}
/* the pair against its two singles and the medium (PREREG §4) */
static void pair_stats(int kind, int R, int k, int p1, int p2, const uint8_t *D0, const uint8_t *D1, const uint8_t *D2, const uint8_t *DP, const Med *m1, const Med *m2,
                       const Med *mp, double E0)
{
    double s1 = 0, s2 = 0, sp = 0, b0 = 0, b1 = 0, b2 = 0, bp = 0;
    long np = 0, nb = 0;
    char *sh1 = malloc((size_t)NO), *sh2 = malloc((size_t)NO);
    for (int u = 0; u < NO; u++) {
        sh1[u] = D0[(size_t)p1 * NO + u] <= 2;
        sh2[u] = D0[(size_t)p2 * NO + u] <= 2;
    }
    for (int u = 0; u < NO; u++) {
        const uint8_t *r0 = D0 + (size_t)u * NO, *r1 = D1 + (size_t)u * NO, *r2 = D2 + (size_t)u * NO, *rp = DP + (size_t)u * NO;
        if (r0[u] == 255 || r1[u] == 255 || r2[u] == 255 || rp[u] == 255) continue; /* u is not present in all four */
        for (int v = 0; v < NO; v++) {
            if (v == u || r0[v] == 255 || r1[v] == 255 || r2[v] == 255 || rp[v] == 255) continue;
            int d0 = r0[v], e1 = r1[v] - d0, e2 = r2[v] - d0, ep = rp[v] - d0;
            if (v > u) {
                s1 += e1;
                s2 += e2;
                sp += ep;
                np++;
            }
            if (sh1[u] && sh2[v]) {
                b0 += d0;
                b1 += e1;
                b2 += e2;
                bp += ep;
                nb++;
            }
        }
    }
    free(sh1);
    free(sh2);
    double S1 = np ? s1 / np : NAN, S2 = np ? s2 / np : NAN, SP = np ? sp / np : NAN;
    double B0 = nb ? b0 / nb : NAN, B1 = nb ? b1 / nb : NAN, B2 = nb ? b2 / nb : NAN, BP = nb ? bp / nb : NAN;
    double rho = (S1 + S2) != 0 ? SP / (S1 + S2) : NAN;
    printf("PAIR kind=%s R=%d k=%d p1=%d p2=%d Eint_made=%.6g Eint_q=%.6g Eint_w=%.6f healed_q=%d healed_end=%d ncl_q=%d ncl_end=%d dcc=%d S_P=%.9f S_1=%.9f S_2=%.9f I=%.9f "
           "rhoW=%.6f B_0=%.6f B_P=%.9f B_1=%.9f B_2=%.9f B_int=%.9f npairs=%ld nshell=%ld\n",
           KIND[kind], R, k, p1, p2, mp->Emade - m1->Emade - m2->Emade, mp->q.E - m1->q.E - m2->q.E + E0, mp->Ew - m1->Ew - m2->Ew + E0, mp->q.nmis == 0, mp->nmis_end == 0, mp->q.ncl,
           mp->ncl_end, mp->dcc, SP, S1, S2, SP - S1 - S2, rho, B0, BP, B1, B2, BP - B1 - B2, np, nb);
    fflush(stdout);
}

/* ---------------- the registered process (PREREG §3) ---------------- */
static const char *setname(void)
{
    static char s[8];
    int c = 0;
    if (OPT_A) s[c++] = 'a';
    if (OPT_B) s[c++] = 'b';
    if (OPT_C) s[c++] = 'c';
    s[c] = 0;
    return c ? s : "none";
}
static void run(int L, int k)
{
    int N = L * L * L, cap = N + 2;
    NO = NQ = N;
    scratch(cap);
    G *s0 = g_new(cap), *h = g_new(cap);
    size_t MS = (size_t)N * N;
    uint8_t *D0 = malloc(MS), *DG1 = malloc(MS), *DX1 = malloc(MS), *DG2 = malloc(MS), *DX2 = malloc(MS), *DP = malloc(MS);
    int *cl = malloc(sizeof(int) * (size_t)cap);
    uint64_t *La = malloc(sizeof(uint64_t) * 3 * (size_t)cap), *Lb = malloc(sizeof(uint64_t) * 3 * (size_t)cap);
    if (!D0 || !DG1 || !DX1 || !DG2 || !DX2 || !DP || !cl || !La || !Lb) {
        fprintf(stderr, "out of memory\n");
        exit(2);
    }
    printf("RUN set=%s L=%d N=%d k=%d seed=%llu T=%g msweeps=%d meas=%d every=%d qr=%d\n", setname(), L, N, k, (unsigned long long)BASE, TREG, M_SWEEPS, MEAS_SWEEPS, SAMPLE_EVERY, QR);
    fflush(stdout);
    /* the medium M: the flat start, M_SWEEPS at t (one stream for every placement) */
    build_cubic_order(s0, L, 0);
    Rng r;
    rng_seed(&r, stream_seed(BASE, 0, 0, RO_M));
    NPROP = NACC = 0;
    sweeps(s0, &r, M_SWEEPS, TREG, 0);
    double accM = NPROP ? (double)NACC / NPROP : 0.0;
    build_cubic_order(h, L, 0);
    int na = clist(h, La), nb = clist(s0, Lb), chg = cdiff(La, na, Lb, nb, NULL);
    Reading r0;
    readm(s0, &r0, D0, NULL);
    printf("MED set=%s L=%d k=%d E=%.6f nmis=%d nD=%d ncl=%d meand=%.9f diam=%d pieces=%d changed_vs_flat=%d acc=%.8g ok=%d\n", setname(), L, k, r0.E, r0.nmis, r0.nD, r0.ncl, r0.meand,
           r0.diam, r0.pieces, chg, accM, check(s0, N));
    /* the first site: uniform over the quanta */
    rng_seed(&r, stream_seed(BASE, k, 0, RO_P1));
    int p1 = rint_(&r, N);
    printf("SITE k=%d p1=%d\n", k, p1);
    fflush(stdout);
    Med mG1, mX1, mG2, mX2, mP;
    do_medium(h, s0, 0, k, 0, p1, -1, DG1, NULL, &mG1);
    do_medium(h, s0, 1, k, 0, p1, -1, DX1, NULL, &mX1);
    for (int a = 0; a < 6; a++) {
        int R = RLIST[a];
        if (2 * R > L) continue;
        /* the second site: uniform among the quanta at distance R from p1, counted through the medium */
        int cnt = 0, p2 = -1;
        for (int v = 0; v < N; v++)
            if (D0[(size_t)p1 * N + v] == R) cnt++;
        if (!cnt) {
            printf("SKIP k=%d R=%d what=no_quantum_at_R\n", k, R);
            continue;
        }
        rng_seed(&r, stream_seed(BASE, k, R, RO_P2));
        int pick = rint_(&r, cnt);
        for (int v = 0, c = 0; v < N; v++)
            if (D0[(size_t)p1 * N + v] == R && c++ == pick) {
                p2 = v;
                break;
            }
        printf("SITE k=%d R=%d p2=%d shell=%d\n", k, R, p2, cnt);
        fflush(stdout);
        do_medium(h, s0, 2, k, R, p1, p2, DG2, NULL, &mG2);
        do_medium(h, s0, 3, k, R, p1, p2, DX2, NULL, &mX2);
        for (int pk = 0; pk < 3; pk++) {
            int kind = 4 + pk;
            const Med *m1 = ISG[kind][0] ? &mG1 : &mX1, *m2 = ISG[kind][1] ? &mG2 : &mX2;
            const uint8_t *M1 = ISG[kind][0] ? DG1 : DX1, *M2 = ISG[kind][1] ? DG2 : DX2;
            do_medium(h, s0, kind, k, R, p1, p2, DP, cl, &mP);
            if (mP.ok && m1->ok && m2->ok)
                pair_stats(kind, R, k, p1, p2, D0, M1, M2, DP, m1, m2, &mP, r0.E);
            else
                printf("FAIL kind=%s R=%d k=%d what=pair_incomplete\n", KIND[kind], R, k);
        }
    }
    printf("DONE\n");
    fflush(stdout);
    g_free(s0);
    g_free(h);
    free(D0);
    free(DG1);
    free(DX1);
    free(DG2);
    free(DX2);
    free(DP);
    free(cl);
    free(La);
    free(Lb);
}

/* ---------------- tool gates (PREREG §6) ---------------- */
static void setopts(int m)
{
    OPT_A = (m >> 2) & 1;
    OPT_B = (m >> 1) & 1;
    OPT_C = m & 1;
}
static int first_at(const G *g, int src, int R) /* the lowest identity at distance R from src */
{
    int *dist = malloc(sizeof(int) * (size_t)g->n), *que = malloc(sizeof(int) * (size_t)g->n), best = -1;
    bfs(g, &src, 1, dist, que);
    for (int i = 0; i < g->n; i++)
        if (dist[i] == R && (best < 0 || g->oid[i] < g->oid[best])) best = i;
    free(dist);
    free(que);
    return best < 0 ? -1 : g->oid[best];
}
static int same_contacts(const G *x, const G *y)
{
    uint64_t *a = malloc(sizeof(uint64_t) * 3 * (size_t)x->cap), *b = malloc(sizeof(uint64_t) * 3 * (size_t)y->cap);
    int na = clist(x, a), nb = clist(y, b), same = x->n == y->n && na == nb && cdiff(a, na, b, nb, NULL) == 0;
    free(a);
    free(b);
    return same;
}
static int t0(void)
{
    int pass = 1;
    init_pm();
    scratch(4200);
    /* T0.1 the carried medium: invariants (with identities) and local = full, L = 6, every set, flat and melted */
    {
        long bad = 0, tested = 0;
        int inv = 0;
        NO = NQ = 216;
        for (int m = 0; m < 8; m++) {
            setopts(m);
            for (int start = 0; start < 2; start++) {
                G *g = g_new(218);
                build_cubic_order(g, 6, 0);
                Rng r;
                rng_seed(&r, 101 + 2 * m + start);
                if (start) sweeps(g, &r, 200, 1, 1);
                for (int p = 0; p < 2000; p++) {
                    Flip f;
                    if (!propose(g, &r, &f)) continue;
                    long E0 = efull(g), dl = flip_de(g, &f);
                    flip_apply(g, &f);
                    long E1 = efull(g);
                    tested++;
                    if (E1 - E0 != dl) bad++;
                    if (runif(&r) < 0.5) flip_revert(g, &f);
                }
                for (int p = 0; p < 5; p++) {
                    for (int s = 0; s < 10000; s++) mc_step(g, &r, TREG, 0);
                    if (!check(g, 216)) inv++;
                }
                if (!start) { /* identities survive a gap and an excess */
                    G *h = g_new(218);
                    g_copy(h, g);
                    double e1 = make_gap_at(h, idx_of(h, 7), &r, NULL);
                    double e2 = make_excess_at(h, idx_of(h, 100), 216, &r);
                    if (isnan(e1) || isnan(e2) || !check(h, 216) || idx_of(h, 7) >= 0 || idx_of(h, 216) < 0) inv++;
                    g_free(h);
                }
                g_free(g);
            }
        }
        printf("T0.1 invariants (six distinct contacts, no self contact, adjacency consistent, one piece, identities distinct) violations=%d %s\n", inv, inv ? "FAIL" : "PASS");
        printf("T0.1 local dE = full dE: tested=%ld mismatches=%ld %s\n", tested, bad, bad ? "FAIL" : "PASS");
        pass &= inv == 0 && bad == 0;
    }
    /* T0.2 candidate tables in identity order: counts, local = full, and the same table and closure whatever the
     * order in which the medium's contacts were entered */
    {
        int ok = 1;
        NO = NQ = 216;
        for (int m = 0; m < 8; m++) {
            setopts(m);
            G *g = g_new(218), *g2 = g_new(218);
            build_cubic_order(g, 6, 0);
            build_cubic_order(g2, 6, 1);
            int J[6], J2[6], v1[15], v2[15], bad = 0, diff = 0, xbad = 0;
            long d1[15], d2[15], dfu[15];
            int nv = gap_table(g, 0, J, v1, d1, dfu);
            gap_table(g2, 0, J2, v2, d2, NULL);
            for (int p = 0; p < 15; p++) {
                if (v1[p] && d1[p] != dfu[p]) bad++;
                if (v1[p] != v2[p] || (v1[p] && d1[p] != d2[p])) diff++;
            }
            for (int q = 0; q < 6; q++)
                if (g->oid[J[q]] != g2->oid[J2[q]]) diff++;
            Edge *E1, *E2;
            Tri *t1, *t2;
            long *x1, *x2, *xf;
            int nt = excess_table(g, 0, &E1, &t1, &x1, &xf), nt2 = excess_table(g2, 0, &E2, &t2, &x2, NULL);
            for (int t = 0; t < nt; t++)
                if (x1[t] != xf[t]) xbad++;
            if (nt != nt2)
                diff++;
            else
                for (int t = 0; t < nt; t++) {
                    if (x1[t] != x2[t]) diff++;
                    for (int e = 0; e < 3; e++)
                        if (E1[t1[t].e[e]].k1 != E2[t2[t].e[e]].k1 || E1[t1[t].e[e]].k2 != E2[t2[t].e[e]].k2) diff++;
                }
            free(E1);
            free(E2);
            free(t1);
            free(t2);
            free(x1);
            free(x2);
            free(xf);
            G *h1 = g_new(218), *h2 = g_new(218);
            g_copy(h1, g);
            g_copy(h2, g2);
            Rng ra, rb;
            rng_seed(&ra, 700 + m);
            rng_seed(&rb, 700 + m);
            double ea = make_gap_at(h1, idx_of(h1, 0), &ra, NULL), eb = make_gap_at(h2, idx_of(h2, 0), &rb, NULL);
            double fa = make_excess_at(h1, idx_of(h1, 50), 216, &ra), fb = make_excess_at(h2, idx_of(h2, 50), 216, &rb);
            int same = same_contacts(h1, h2) && ea == eb && fa == fb && check(h1, 216) && check(h2, 216);
            int oka = nv == 15 && nt == 3712 && bad == 0 && xbad == 0 && diff == 0 && same;
            ok &= oka;
            printf("T0.2 set=%s gap: valid=%d (expect 15) local!=full=%d; excess: triples=%d (expect 3712) local!=full=%d; tables differ by entry order=%d; same closure=%d %s\n",
                   setname(), nv, bad, nt, xbad, diff, same, oka ? "PASS" : "FAIL");
            g_free(g);
            g_free(g2);
            g_free(h1);
            g_free(h2);
        }
        pass &= ok;
    }
    /* T0.3 the restricted quench equals the full quench, and ends at a local minimum: L = 6 and 8, every set with
     * an energy, five configurations */
    {
        int bad = 0, notmin = 0, cases = 0;
        for (int m = 1; m < 8; m++) {
            setopts(m);
            for (int L = 6; L <= 8; L += 2) {
                int N = L * L * L;
                NO = NQ = N;
                scratch(N + 2);
                G *g = g_new(N + 2), *h1 = g_new(N + 2), *h2 = g_new(N + 2);
                build_cubic_order(g, L, 0);
                static const int ck[5] = {0, 1, 4, 5, 6}, cR[5] = {0, 0, 1, 1, 2};
                for (int cs = 0; cs < 5; cs++) {
                    int R = cR[cs], p2 = R ? first_at(g, 0, R) : -1;
                    double Em;
                    int nexp;
                    BASE = 300 + 10 * m + L;
                    if (!make_kind(h1, g, ck[cs], 0, R, 0, p2, &Em, &nexp)) {
                        bad++;
                        continue;
                    }
                    g_copy(h2, h1);
                    long d1, d2, d3;
                    int s1 = quench(h1, 0, &d1), s2 = quench(h2, 1, &d2);
                    if (s1 != s2 || d1 != d2 || !same_contacts(h1, h2)) bad++;
                    if (quench(h2, 1, &d3) != 0) notmin++;
                    cases++;
                }
                g_free(g);
                g_free(h1);
                g_free(h2);
            }
        }
        int ok = bad == 0 && notmin == 0;
        printf("T0.3 restricted quench = full quench (L = 6, 8; 7 sets; gap, excess, GG R=1, XX R=1, GX R=2): cases=%d differ=%d not_minimum=%d %s\n", cases, bad, notmin, ok ? "PASS" : "FAIL");
        pass &= ok;
    }
    /* T0.4 out of reach, a pair is its two singles exactly: L = 16, R = 16, every set, GG, XX and GX */
    {
        int ok = 1, L = 16, N = 4096;
        NO = NQ = N;
        scratch(N + 2);
        G *g = g_new(N + 2), *hp = g_new(N + 2), *h1 = g_new(N + 2), *h2 = g_new(N + 2);
        uint64_t *Lg = malloc(sizeof(uint64_t) * 3 * (N + 2)), *Lp = malloc(sizeof(uint64_t) * 3 * (N + 2)), *L1 = malloc(sizeof(uint64_t) * 3 * (N + 2)),
                 *L2 = malloc(sizeof(uint64_t) * 3 * (N + 2));
        uint64_t *ap = malloc(sizeof(uint64_t) * 64), *a1 = malloc(sizeof(uint64_t) * 64), *a2 = malloc(sizeof(uint64_t) * 64);
        uint64_t *rp = malloc(sizeof(uint64_t) * 64), *r1 = malloc(sizeof(uint64_t) * 64), *r2 = malloc(sizeof(uint64_t) * 64);
        for (int m = 0; m < 8; m++) {
            setopts(m);
            build_cubic_order(g, L, 0);
            int p2 = first_at(g, 0, 16), ng = clist(g, Lg);
            BASE = 400 + m;
            for (int pk = 0; pk < 3; pk++) {
                int kind = 4 + pk, k1 = ISG[kind][0] ? 0 : 1, k2 = ISG[kind][1] ? 2 : 3, ne;
                double Ep, E1, E2;
                long q;
                int okp = make_kind(hp, g, kind, 0, 16, 0, p2, &Ep, &ne), ok1 = make_kind(h1, g, k1, 0, 16, 0, p2, &E1, &ne), ok2 = make_kind(h2, g, k2, 0, 16, 0, p2, &E2, &ne);
                int madeok = okp && ok1 && ok2 && Ep == E1 + E2;
                quench(hp, 0, &q);
                quench(h1, 0, &q);
                quench(h2, 0, &q);
                long ep = efull(hp), e1 = efull(h1), e2 = efull(h2);
                int np = clist(hp, Lp), n1 = clist(h1, L1), n2 = clist(h2, L2);
                /* contacts added (not in the flat medium) and removed, for the pair and each single */
                int nap = cdiff(Lg, ng, Lp, np, NULL), na1 = cdiff(Lg, ng, L1, n1, NULL), na2 = cdiff(Lg, ng, L2, n2, NULL);
                int nrp = cdiff(Lp, np, Lg, ng, NULL), nr1 = cdiff(L1, n1, Lg, ng, NULL), nr2 = cdiff(L2, n2, Lg, ng, NULL);
                int dec = 0;
                if (nap <= 64 && na1 <= 64 && na2 <= 64 && nrp <= 64 && nr1 <= 64 && nr2 <= 64) {
                    cdiff(Lg, ng, Lp, np, ap);
                    cdiff(Lg, ng, L1, n1, a1);
                    cdiff(Lg, ng, L2, n2, a2);
                    cdiff(Lp, np, Lg, ng, rp);
                    cdiff(L1, n1, Lg, ng, r1);
                    cdiff(L2, n2, Lg, ng, r2);
                    dec = nap == na1 + na2 && nrp == nr1 + nr2 && cdiff(ap, nap, a1, na1, NULL) == 0 && cdiff(ap, nap, a2, na2, NULL) == 0 && cdiff(rp, nrp, r1, nr1, NULL) == 0 &&
                          cdiff(rp, nrp, r2, nr2, NULL) == 0;
                }
                int oka = madeok && ep == e1 + e2 && dec;
                ok &= oka;
                printf("T0.4 set=%s %s at R=16: made dE pair=%.6g singles=%.6g+%.6g; quenched E pair=%ld singles=%ld+%ld (2t units); contacts added %d=%d+%d removed %d=%d+%d; "
                       "decomposes=%d %s\n",
                       setname(), KIND[kind], Ep, E1, E2, ep, e1, e2, nap, na1, na2, nrp, nr1, nr2, dec, oka ? "PASS" : "FAIL");
            }
        }
        free(Lg);
        free(Lp);
        free(L1);
        free(L2);
        free(ap);
        free(a1);
        free(a2);
        free(rp);
        free(r1);
        free(r2);
        g_free(g);
        g_free(hp);
        g_free(h1);
        g_free(h2);
        pass &= ok;
    }
    /* T0.5 healing is reachable and readable: a gap, then an excess on exactly the three contacts its closure made,
     * restores the flat medium (L = 8, every set) */
    {
        int ok = 1, L = 8, N = 512;
        NO = NQ = N;
        scratch(N + 2);
        G *g = g_new(N + 2), *h = g_new(N + 2);
        uint8_t *M0 = malloc((size_t)N * N), *MP = malloc((size_t)N * N);
        double mex = 0.75 * L * N / (N - 1.0);
        for (int m = 0; m < 8; m++) {
            setopts(m);
            build_cubic_order(g, L, 0);
            g_copy(h, g);
            Rng r;
            rng_seed(&r, 500 + m);
            int pr[6], u[3], v[3];
            double e = make_gap_at(h, 0, &r, pr);
            for (int q = 0; q < 3; q++) {
                u[q] = idx_of(h, pr[2 * q]);
                v[q] = idx_of(h, pr[2 * q + 1]);
            }
            insert_excess(h, u, v, N);
            Reading r0, rd;
            readm(g, &r0, M0, NULL);
            readm(h, &rd, MP, NULL);
            long dsum = 0, ndif = 0;
            for (size_t x = 0; x < (size_t)N * N; x++)
                if (M0[x] != 255 && MP[x] != 255) {
                    dsum += (long)MP[x] - (long)M0[x];
                    if (MP[x] != M0[x]) ndif++;
                }
            int oka = !isnan(e) && rd.nmis == 0 && rd.E == 0 && fabs(rd.meand - mex) < 1e-9 && check(h, N) && ndif == 0 && dsum == 0;
            ok &= oka;
            printf("T0.5 set=%s gap dE=%.6g, then the excess on its closure: nmis=%d E=%.6g meand=%.9f (flat %.9f) distances changed=%ld %s\n", setname(), e, rd.nmis, rd.E, rd.meand, mex,
                   ndif, oka ? "PASS" : "FAIL");
        }
        free(M0);
        free(MP);
        g_free(g);
        g_free(h);
        pass &= ok;
    }
    /* T0.6 the instruments on known states (L = 8, every set): the flat medium; one gap (the matrix mean equals the
     * mean distance); two gaps at the far side of the medium read as two clusters at their separation */
    {
        int ok = 1, L = 8, N = 512;
        NO = NQ = N;
        scratch(N + 2);
        G *g = g_new(N + 2), *h = g_new(N + 2);
        uint8_t *MG = malloc((size_t)N * N);
        int *cl = malloc(sizeof(int) * (N + 2));
        double mex = 0.75 * L * N / (N - 1.0);
        for (int m = 0; m < 8; m++) {
            setopts(m);
            build_cubic_order(g, L, 0);
            Reading r0, rg, rp;
            readm(g, &r0, NULL, NULL);
            int okf = r0.nmis == 0 && r0.ncl == 0 && fabs(r0.meand - mex) < 1e-9 && r0.diam == 3 * L / 2;
            BASE = 600 + m;
            double Em;
            int nexp;
            make_kind(h, g, 0, 0, 0, 0, -1, &Em, &nexp);
            readm(h, &rg, MG, NULL);
            double s = 0;
            long c = 0;
            for (int u = 0; u < N; u++)
                for (int v = u + 1; v < N; v++)
                    if (MG[(size_t)u * N + v] != 255) {
                        s += MG[(size_t)u * N + v];
                        c++;
                    }
            int okg = c > 0 && fabs(s / c - rg.meand) < 1e-9 && rg.ncl == 1;
            int p2 = first_at(g, 0, 12);
            make_kind(h, g, 4, 0, 12, 0, p2, &Em, &nexp);
            readm(h, &rp, NULL, cl);
            int dcc = rp.ncl == 2 ? centre_dist(h, cl) : -1;
            int okp = rp.ncl == 2 && dcc >= 8 && dcc <= 12;
            int oka = okf && okg && okp;
            ok &= oka;
            printf("T0.6 set=%s flat: nmis=%d meand=%.9f (exact %.9f) diam=%d; one gap: matrix mean %.9f = meand %.9f, clusters=%d; two gaps at 12: clusters=%d centre distance=%d %s\n",
                   setname(), r0.nmis, r0.meand, mex, r0.diam, c ? s / c : NAN, rg.meand, rg.ncl, rp.ncl, dcc, oka ? "PASS" : "FAIL");
        }
        free(MG);
        free(cl);
        g_free(g);
        g_free(h);
        pass &= ok;
    }
    printf("T0 %s\n", pass ? "PASS" : "FAIL");
    return pass ? 0 : 1;
}

int main(int argc, char **argv)
{
    const char *mode = "run";
    int L = 8, a = 0, b = 0, c = 0, k = 0;
    long long seed = -1;
    for (int i = 1; i < argc; i++) {
        if (!strncmp(argv[i], "mode=", 5)) mode = argv[i] + 5;
        else if (!strncmp(argv[i], "L=", 2)) L = atoi(argv[i] + 2);
        else if (!strncmp(argv[i], "a=", 2)) a = atoi(argv[i] + 2);
        else if (!strncmp(argv[i], "b=", 2)) b = atoi(argv[i] + 2);
        else if (!strncmp(argv[i], "c=", 2)) c = atoi(argv[i] + 2);
        else if (!strncmp(argv[i], "k=", 2)) k = atoi(argv[i] + 2);
        else if (!strncmp(argv[i], "seed=", 5)) seed = atoll(argv[i] + 5);
        else {
            fprintf(stderr, "unknown argument %s\n", argv[i]);
            return 2;
        }
    }
    if (!strcmp(mode, "t0")) return t0();
    if (strcmp(mode, "run") || L < 6 || L > 16 || k < 0 || k > 9) {
        fprintf(stderr, "usage: ./c3pair mode=run L=12 a=1 b=0 c=1 k=0 | mode=t0\n");
        return 2;
    }
    OPT_A = a != 0;
    OPT_B = b != 0;
    OPT_C = c != 0;
    init_pm();
    /* the registered base seed (PREREG §3) */
    BASE = seed >= 0 ? (uint64_t)seed : (uint64_t)(2000003LL * L + 7LL * (4 * OPT_A + 2 * OPT_B + OPT_C) + 1);
    run(L, k);
    return 0;
}
