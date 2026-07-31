/* maxfab — EM5 PROTOTYPE (standalone; not the production kernel).
 *
 * Discrete source-free Maxwell on the foam's triangle complex, per the
 * EMF.md 2026-07-31 design: field 1-form E on links, 2-form B on
 * triangles, leapfrog pair
 *     B_t -= dt * C * (circulation of E around t)
 *     E_l += dt * C * (sum of B over triangles containing l, signed)
 * — the discrete exterior-calculus curl pair (unweighted hodge, v1).
 *
 * QUESTION (pre-registered): does this operator give a LINEAR cone —
 * v_g(k) ~ flat — where the production scalar-hop field sector gives
 * the measured band v_g = 0.610*sin(1.410 k) (EM1)? If yes, EM5's law
 * design (link-borne chiral pair) has a measured discrete operator to
 * build on. Standalone by doctrine: design first, production later.
 *
 * usage: maxfab foam.tsv kx [T dt amp sigma]
 * out:   per-diag rows "t  centroid_x  E_tot"; final "# RESULT vg=..."
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int NC = 0, NL = 0, NT = 0;
static double *cx, *cy, *cz, *cr;
static int *li, *lj;
static double *lmx, *lmy, *lmz, *lux, *luy, *luz;
static int *t1_, *t2_, *t3_;         /* triangle edge ids */
static signed char *s1_, *s2_, *s3_; /* orientation signs */
static double *E, *B;

int main(int argc, char **argv)
{
    if (argc < 3) { fprintf(stderr, "usage: %s foam.tsv kx [T dt amp sigma]\n", argv[0]); return 2; }
    double kx = atof(argv[2]);
    double T = argc > 3 ? atof(argv[3]) : 30.0;
    double dt = argc > 4 ? atof(argv[4]) : 0.005;
    double amp = argc > 5 ? atof(argv[5]) : 1.0;
    double sig = argc > 6 ? atof(argv[6]) : 2.5;
    double C = 1.0;

    /* ---- foam ---- */
    FILE *f = fopen(argv[1], "r");
    if (!f) { fprintf(stderr, "no foam\n"); return 2; }
    char ln[1024];
    if (!fgets(ln, sizeof(ln), f)) return 2;   /* header */
    int cap = 16384;
    cx = malloc(cap * sizeof(double)); cy = malloc(cap * sizeof(double));
    cz = malloc(cap * sizeof(double)); cr = malloc(cap * sizeof(double));
    while (fgets(ln, sizeof(ln), f)) {
        int id; double x, y, z, r;
        if (sscanf(ln, "%d %lf %lf %lf %lf", &id, &x, &y, &z, &r) != 5) continue;
        cx[NC] = x; cy[NC] = y; cz[NC] = z; cr[NC] = r; NC++;
        if (NC >= cap) { cap *= 2;
            cx = realloc(cx, cap * sizeof(double)); cy = realloc(cy, cap * sizeof(double));
            cz = realloc(cz, cap * sizeof(double)); cr = realloc(cr, cap * sizeof(double)); }
    }
    fclose(f);

    /* ---- links (the production rule) + per-cell adjacency ---- */
    int lcap = NC * 12;
    li = malloc(lcap * sizeof(int)); lj = malloc(lcap * sizeof(int));
    for (int i = 0; i < NC; i++)
        for (int j = i + 1; j < NC; j++) {
            double dx = cx[j] - cx[i], dy = cy[j] - cy[i], dz = cz[j] - cz[i];
            if (fabs(dx) > 2.1 || fabs(dy) > 2.1 || fabs(dz) > 2.1) continue;
            double d2 = dx * dx + dy * dy + dz * dz;
            double cut = 1.15 * (cr[i] + cr[j]);
            if (d2 >= cut * cut) continue;
            if (NL >= lcap) { lcap *= 2;
                li = realloc(li, lcap * sizeof(int)); lj = realloc(lj, lcap * sizeof(int)); }
            li[NL] = i; lj[NL] = j; NL++;
        }
    lmx = malloc(NL * sizeof(double)); lmy = malloc(NL * sizeof(double));
    lmz = malloc(NL * sizeof(double));
    lux = malloc(NL * sizeof(double)); luy = malloc(NL * sizeof(double));
    luz = malloc(NL * sizeof(double));
    /* link id lookup: per-cell sorted neighbor lists */
    int *deg = calloc(NC, sizeof(int));
    for (int l = 0; l < NL; l++) { deg[li[l]]++; deg[lj[l]]++; }
    int *off = malloc((NC + 1) * sizeof(int));
    off[0] = 0;
    for (int i = 0; i < NC; i++) off[i + 1] = off[i] + deg[i];
    int *adjc = malloc(2 * (size_t)NL * sizeof(int));
    int *adjl = malloc(2 * (size_t)NL * sizeof(int));
    int *fill = calloc(NC, sizeof(int));
    for (int l = 0; l < NL; l++) {
        int a = li[l], b = lj[l];
        adjc[off[a] + fill[a]] = b; adjl[off[a] + fill[a]] = l; fill[a]++;
        adjc[off[b] + fill[b]] = a; adjl[off[b] + fill[b]] = l; fill[b]++;
        double d = 0;
        double dx = cx[b] - cx[a], dy = cy[b] - cy[a], dz = cz[b] - cz[a];
        d = sqrt(dx * dx + dy * dy + dz * dz);
        lmx[l] = 0.5 * (cx[a] + cx[b]); lmy[l] = 0.5 * (cy[a] + cy[b]);
        lmz[l] = 0.5 * (cz[a] + cz[b]);
        lux[l] = dx / d; luy[l] = dy / d; luz[l] = dz / d;
    }

    /* ---- triangles: i<j<k pairwise linked; boundary (ij)+(jk)-(ik) ---- */
    int tcap = NL * 4;
    t1_ = malloc(tcap * sizeof(int)); t2_ = malloc(tcap * sizeof(int));
    t3_ = malloc(tcap * sizeof(int));
    s1_ = malloc(tcap); s2_ = malloc(tcap); s3_ = malloc(tcap);
    for (int l = 0; l < NL; l++) {
        int i = li[l], j = lj[l];          /* i<j by construction */
        for (int q = off[j]; q < off[j + 1]; q++) {
            int k = adjc[q];
            if (k <= j) continue;          /* enforce i<j<k via k>j */
            /* need link (i,k) */
            int lik = -1;
            for (int p = off[i]; p < off[i + 1]; p++)
                if (adjc[p] == k) { lik = adjl[p]; break; }
            if (lik < 0) continue;
            if (NT >= tcap) { tcap *= 2;
                t1_ = realloc(t1_, tcap * sizeof(int)); t2_ = realloc(t2_, tcap * sizeof(int));
                t3_ = realloc(t3_, tcap * sizeof(int));
                s1_ = realloc(s1_, tcap); s2_ = realloc(s2_, tcap); s3_ = realloc(s3_, tcap); }
            t1_[NT] = l;        s1_[NT] = +1;   /* (i,j) */
            t2_[NT] = adjl[q];  s2_[NT] = +1;   /* (j,k) */
            t3_[NT] = lik;      s3_[NT] = -1;   /* -(i,k) */
            NT++;
        }
    }
    fprintf(stderr, "# maxfab: NC=%d NL=%d NT=%d kx=%.2f dt=%g\n", NC, NL, NT, kx, dt);

    /* ---- packet: E along y, gaussian at x0, carrier kx along x ---- */
    E = calloc(NL, sizeof(double));
    B = calloc(NT, sizeof(double));
    double x0 = 6.0, yc = 12.0, zc = 12.0;
    for (int l = 0; l < NL; l++) {
        double g = exp(-((lmx[l] - x0) * (lmx[l] - x0)
                       + (lmy[l] - yc) * (lmy[l] - yc)
                       + (lmz[l] - zc) * (lmz[l] - zc)) / (2 * sig * sig));
        E[l] = amp * g * cos(kx * lmx[l]) * luy[l];
    }

    /* ---- probes: sample links near the packet center; the frozen
     * (curl-kernel) component sits at omega=0, the propagating part at
     * omega(k) — spectral peak separates them cleanly ---- */
    int npr = 0, prl[24];
    for (int l = 0; l < NL && npr < 24; l++) {
        double dx = lmx[l] - x0, dy = lmy[l] - yc, dz = lmz[l] - zc;
        if (dx * dx + dy * dy + dz * dz < 4.0 && fabs(luy[l]) > 0.7)
            prl[npr++] = l;
    }
    long NS = (long)(T / dt + 0.5);
    double *rec = malloc((NS + 1) * sizeof(double));
    long dg = (long)(1.0 / dt + 0.5);
    printf("# t\tcx\tEtot\n");
    for (long s = 0; s <= NS; s++) {
        if (s % dg == 0) {
            double et = 0, cxw = 0;
            for (int l = 0; l < NL; l++) { double e = E[l] * E[l]; et += e; cxw += e * lmx[l]; }
            double bt = 0;
            for (int t = 0; t < NT; t++) bt += B[t] * B[t];
            printf("%.2f\t%.4f\t%.6g\n", s * dt, et > 0 ? cxw / et : 0, 0.5 * (et + bt));
            fflush(stdout);
        }
        {
            double pv = 0;
            for (int p = 0; p < npr; p++) pv += E[prl[p]];
            rec[s] = pv;
        }
        for (int t = 0; t < NT; t++)
            B[t] -= dt * C * (s1_[t] * E[t1_[t]] + s2_[t] * E[t2_[t]] + s3_[t] * E[t3_[t]]);
        for (int t = 0; t < NT; t++) {
            double b = dt * C * B[t];
            E[t1_[t]] += b * s1_[t];
            E[t2_[t]] += b * s2_[t];
            E[t3_[t]] += b * s3_[t];
        }
    }
    /* spectral peak: DFT of the probe record over a frequency grid,
     * skipping omega < 0.05 (the frozen kernel line at DC) */
    {
        double best_w = 0, best_p = -1;
        for (double w = 0.05; w <= 3.5; w += 0.01) {
            double cr_ = 0, ci_ = 0;
            for (long s = 0; s <= NS; s++) {
                double t = s * dt;
                cr_ += rec[s] * cos(w * t);
                ci_ += rec[s] * sin(w * t);
            }
            double p = cr_ * cr_ + ci_ * ci_;
            if (p > best_p) { best_p = p; best_w = w; }
        }
        printf("# RESULT maxfab kx=%.2f omega_peak=%.3f (probes=%d)\n",
               kx, best_w, npr);
    }
    return 0;
}
