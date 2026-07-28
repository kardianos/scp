/*
 * construct_gauss.c — combinatorial Gauss law toy for v88/CONSTRUCTION.md §4
 *
 * No embedding, no ambient coordinates. Build abstract space graphs with a
 * designated interior holding dense energy M (parcel count removed from the
 * space graph and boundary rewired). Measure degree defect
 *
 *   K(I) = d_* * |boundary-adjacent space| contribution form:
 *   K     = sum_v (d_* - deg(v)) over all space vertices
 *
 * and compare to enclosed M (number of removed parcels * unit mass).
 *
 * Claim under test [P]: K ≈ κ M / ε_0  (linear relation, κ from foam).
 *
 * Build: gcc -O2 -o construct_gauss construct_gauss.c -lm
 * Run:   ./construct_gauss
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define VMAX 256
#define EMAX 2048

typedef struct {
    int n;
    int deg[VMAX];
    int adj[VMAX][16];
    int nadj[VMAX];
} Graph;

static void g_clear(Graph *g)
{
    memset(g, 0, sizeof(*g));
}

static void g_add_edge(Graph *g, int u, int v)
{
    if (u == v || u < 0 || v < 0 || u >= g->n || v >= g->n)
        return;
    for (int i = 0; i < g->nadj[u]; i++)
        if (g->adj[u][i] == v)
            return;
    if (g->nadj[u] >= 16 || g->nadj[v] >= 16)
        return;
    g->adj[u][g->nadj[u]++] = v;
    g->adj[v][g->nadj[v]++] = u;
    g->deg[u]++;
    g->deg[v]++;
}

/* Regular-ish foam patch: triangular lattice disk (combinatorial), degree
 * target d_*=6 in the bulk. No coordinates used for physics — only to
 * generate a plausible abstract graph once. */
static void build_tri_disk(Graph *g, int radius)
{
    /* Map axial hex coords (i,j) with |i|,|j|,|i+j| <= radius to vertices */
    int map[64][64];
    for (int i = 0; i < 64; i++)
        for (int j = 0; j < 64; j++)
            map[i][j] = -1;
    g_clear(g);
    int off = 32;
    for (int i = -radius; i <= radius; i++) {
        for (int j = -radius; j <= radius; j++) {
            int k = -i - j;
            if (abs(i) > radius || abs(j) > radius || abs(k) > radius)
                continue;
            map[i + off][j + off] = g->n++;
        }
    }
    int nbr[6][2] = {{1, 0}, {0, 1}, {-1, 1}, {-1, 0}, {0, -1}, {1, -1}};
    for (int i = -radius; i <= radius; i++) {
        for (int j = -radius; j <= radius; j++) {
            int u = map[i + off][j + off];
            if (u < 0)
                continue;
            for (int t = 0; t < 6; t++) {
                int i2 = i + nbr[t][0], j2 = j + nbr[t][1];
                if (i2 + off < 0 || i2 + off >= 64 || j2 + off < 0
                    || j2 + off >= 64)
                    continue;
                int v = map[i2 + off][j2 + off];
                if (v < 0)
                    continue;
                if (u < v)
                    g_add_edge(g, u, v);
            }
        }
    }
}

/* Remove a connected block of n_rem vertices starting from seed; rewire
 * boundary by connecting former neighbors of removed set to each other
 * in a cycle (abstract reseal). Returns M = n_rem (ε_0 = 1). */
static int convert_interior(Graph *g, int seed, int n_rem, int *removed)
{
    if (n_rem <= 0 || seed < 0 || seed >= g->n)
        return 0;
    int mark[VMAX];
    memset(mark, 0, sizeof(mark));
    int q[VMAX], qh = 0, qt = 0;
    q[qt++] = seed;
    mark[seed] = 1;
    int got = 0;
    while (qh < qt && got < n_rem) {
        int u = q[qh++];
        removed[got++] = u;
        for (int i = 0; i < g->nadj[u]; i++) {
            int v = g->adj[u][i];
            if (!mark[v]) {
                mark[v] = 1;
                q[qt++] = v;
            }
        }
    }
    n_rem = got;

    /* boundary: space vertices adjacent to removed */
    int is_rem[VMAX];
    memset(is_rem, 0, sizeof(is_rem));
    for (int i = 0; i < n_rem; i++)
        is_rem[removed[i]] = 1;

    int bnd[VMAX], nb = 0;
    int on_bnd[VMAX];
    memset(on_bnd, 0, sizeof(on_bnd));
    for (int i = 0; i < n_rem; i++) {
        int u = removed[i];
        for (int j = 0; j < g->nadj[u]; j++) {
            int v = g->adj[u][j];
            if (!is_rem[v] && !on_bnd[v]) {
                on_bnd[v] = 1;
                bnd[nb++] = v;
            }
        }
    }

    /* Build new graph on non-removed vertices */
    Graph h;
    g_clear(&h);
    int new_id[VMAX];
    for (int i = 0; i < VMAX; i++)
        new_id[i] = -1;
    for (int u = 0; u < g->n; u++) {
        if (is_rem[u])
            continue;
        new_id[u] = h.n++;
    }
    /* keep edges not touching removed */
    for (int u = 0; u < g->n; u++) {
        if (is_rem[u])
            continue;
        for (int j = 0; j < g->nadj[u]; j++) {
            int v = g->adj[u][j];
            if (is_rem[v] || v < u)
                continue;
            g_add_edge(&h, new_id[u], new_id[v]);
        }
    }
    /* reseal: connect boundary vertices in a cycle (abstract, no geometry) */
    for (int i = 0; i < nb; i++) {
        int a = new_id[bnd[i]];
        int b = new_id[bnd[(i + 1) % nb]];
        if (a >= 0 && b >= 0)
            g_add_edge(&h, a, b);
    }

    *g = h;
    return n_rem; /* M / ε_0 */
}

static double curvature(const Graph *g, int d_star)
{
    /* total combinatorial curvature: sum (d_* - deg)
     * On a closed surface this is topological; on a finite patch with
     * boundary it includes boundary deficit. We compare ΔK relative to
     * empty conversion (same outer size). */
    double K = 0;
    for (int v = 0; v < g->n; v++)
        K += (double)(d_star - g->deg[v]);
    return K;
}

int main(void)
{
    const int d_star = 6;
    const int radius = 5;
    printf("# construct_gauss — combinatorial K vs M (no embedding in physics)\n");
    printf("# d_*=%d tri-disk radius=%d\n", d_star, radius);
    printf("# M\tK\tK_minus_K0\tn_space\tmean_deg\n");

    Graph g0;
    build_tri_disk(&g0, radius);
    double K0 = curvature(&g0, d_star);
    printf("# baseline n=%d K0=%.6f mean_deg=%.4f\n", g0.n, K0,
           g0.n ? (2.0 * /*edges*/ 0) : 0.0);
    /* mean deg */
    long sdeg = 0;
    for (int i = 0; i < g0.n; i++)
        sdeg += g0.deg[i];
    printf("# baseline mean_deg=%.4f\n", g0.n ? (double)sdeg / g0.n : 0.0);

    /* seed = center-ish vertex 0 is rim; pick vertex with max deg as bulk */
    int seed = 0;
    for (int i = 1; i < g0.n; i++)
        if (g0.deg[i] > g0.deg[seed])
            seed = i;

    for (int M = 0; M <= 12; M++) {
        Graph g = g0;
        int removed[VMAX];
        int m = 0;
        if (M > 0)
            m = convert_interior(&g, seed, M, removed);
        double K = curvature(&g, d_star);
        long sd = 0;
        for (int i = 0; i < g.n; i++)
            sd += g.deg[i];
        double md = g.n ? (double)sd / g.n : 0.0;
        printf("%d\t%.6f\t%.6f\t%d\t%.4f\n", m, K, K - K0, g.n, md);
    }

    /* Linear fit K-K0 = kappa * M on M>=1 */
    double sum_x = 0, sum_y = 0, sum_xx = 0, sum_xy = 0;
    int npt = 0;
    /* recompute for fit */
    for (int M = 1; M <= 12; M++) {
        Graph g = g0;
        int removed[VMAX];
        int m = convert_interior(&g, seed, M, removed);
        double y = curvature(&g, d_star) - K0;
        double x = (double)m;
        sum_x += x;
        sum_y += y;
        sum_xx += x * x;
        sum_xy += x * y;
        npt++;
    }
    double denom = npt * sum_xx - sum_x * sum_x;
    double kappa = denom != 0 ? (npt * sum_xy - sum_x * sum_y) / denom : 0;
    double b = (sum_y - kappa * sum_x) / npt;
    /* residual */
    double ss_res = 0, ss_tot = 0, ymean = sum_y / npt;
    for (int M = 1; M <= 12; M++) {
        Graph g = g0;
        int removed[VMAX];
        int m = convert_interior(&g, seed, M, removed);
        double y = curvature(&g, d_star) - K0;
        double yhat = kappa * m + b;
        ss_res += (y - yhat) * (y - yhat);
        ss_tot += (y - ymean) * (y - ymean);
    }
    double r2 = ss_tot > 0 ? 1.0 - ss_res / ss_tot : 0;
    printf("# fit: (K-K0) = kappa*M + b  with kappa=%.6f b=%.6f r2=%.6f\n",
           kappa, b, r2);
    printf("# interpretation: linear K∝M is the combinatorial Gauss claim [P];\n");
    printf("# kappa is foam-dependent, not universal. No ambient metric used.\n");
    return 0;
}
