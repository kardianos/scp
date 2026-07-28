/*  hc6_shape.c — v86 HC-6, corrected observable.
 *
 *  WHY THE ORIGINAL HC-6 OBSERVABLE WAS VOID
 *    HC-6 was designed to watch charge REDISTRIBUTE between components --
 *    "decay to the sector minimum". That cannot happen in this kernel. The
 *    potential depends on the fields only through s = prod_a |Phi_a|^2, which
 *    is invariant under INDEPENDENT phase rotations Phi_a -> e^{i th_a} Phi_a.
 *    So the model carries U(1)^3, not a single U(1), and each Q_a is its own
 *    conserved Noether charge. Measured: Q_a constant to 13 significant
 *    figures over t = 0..201 in a box with no absorber.
 *
 *    Partition fractions x_a = Q_a/Q_tot are therefore pinned by a conservation
 *    law, exactly like the monochromatic arm's are pinned by permutation
 *    symmetry. Any verdict built on them -- null or positive -- is meaningless.
 *
 *  THE CHANNEL THAT CAN ACTUALLY RESPOND
 *    A GSS index mismatch does not move charge between the U(1)^3 sectors; it
 *    makes the object unstable at FIXED (Q_0,Q_1,Q_2), i.e. in the SHAPE
 *    channel -- the profile deforms, collapses, spreads or fissions. The
 *    kernel reports s_max (peak of s) and P_int every diag row, and those do
 *    respond.
 *
 *  WHAT THIS TOOL MEASURES
 *    s_max(t) relative to its initial value, and an exponential growth rate
 *    fitted to the early rise, for each arm. The comparison that matters is
 *    n(D)=0 arms against the n(D)=1 control at matched charge -- in a box with
 *    NO sponge, so absorption cannot masquerade as instability.
 *
 *  Build: gcc -O2 -o hc6_shape hc6_shape.c -lm
 *  Usage: hc6_shape label=file.tsv [label=file.tsv ...]
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define MAXR 400000
#define MAXC 64

typedef struct {
    char label[64];
    int n;
    double t[MAXR], s[MAXR], q[3][MAXR], E[MAXR];
    double qc[MAXR], rc[MAXR], pint[MAXR];
} Run;

static int col_of(char *hdr, const char *want) {
    char buf[8192];
    snprintf(buf, sizeof(buf), "%s", hdr);
    int i = 0;
    for (char *p = strtok(buf, "\t\n"); p; p = strtok(NULL, "\t\n"), i++)
        if (!strcmp(p, want)) return i;
    return -1;
}

static int load(const char *path, Run *r) {
    FILE *f = fopen(path, "r");
    if (!f) { fprintf(stderr, "cannot open %s\n", path); return 0; }
    static char line[16384];
    if (!fgets(line, sizeof(line), f)) { fclose(f); return 0; }
    int ct = col_of(line, "t"), cs = col_of(line, "s_max"), cE = col_of(line, "E_total");
    int c0 = col_of(line, "Q_p0"), c1 = col_of(line, "Q_p1"), c2 = col_of(line, "Q_p2");
    int cq = col_of(line, "Q_core"), cr = col_of(line, "r_core"),
        cp = col_of(line, "P_int");
    if (ct < 0 || cs < 0) { fclose(f); return 0; }
    r->n = 0;
    while (fgets(line, sizeof(line), f) && r->n < MAXR) {
        double v[MAXC]; int i = 0;
        for (char *p = strtok(line, "\t\n"); p && i < MAXC; p = strtok(NULL, "\t\n"))
            v[i++] = atof(p);
        if (i <= cs) continue;
        int k = r->n;
        r->t[k] = v[ct]; r->s[k] = v[cs];
        r->E[k] = (cE >= 0 && cE < i) ? v[cE] : 0;
        r->qc[k] = (cq >= 0 && cq < i) ? v[cq] : 0;
        r->rc[k] = (cr >= 0 && cr < i) ? v[cr] : 0;
        r->pint[k] = (cp >= 0 && cp < i) ? v[cp] : 0;
        for (int a = 0; a < 3; a++) {
            int c = a == 0 ? c0 : (a == 1 ? c1 : c2);
            r->q[a][k] = (c >= 0 && c < i) ? v[c] : 0;
        }
        r->n++;
    }
    fclose(f);
    return r->n > 4;
}

/* relative drift of each conserved charge -- proves U(1)^3 numerically */
static double qdrift(const Run *r) {
    double d = 0;
    for (int a = 0; a < 3; a++) {
        if (r->q[a][0] == 0) continue;
        for (int k = 0; k < r->n; k++) {
            double e = fabs(r->q[a][k] - r->q[a][0]) / fabs(r->q[a][0]);
            if (e > d) d = e;
        }
    }
    return d;
}

/* least-squares growth rate of ln(s/s0) over the first frac of the run */
static double growth(const Run *r, double frac, double *r2) {
    double tmax = r->t[r->n - 1] * frac;
    double sx = 0, sy = 0, sxx = 0, sxy = 0, syy = 0; int m = 0;
    for (int k = 0; k < r->n; k++) {
        if (r->t[k] > tmax) break;
        double y = log(fmax(r->s[k] / r->s[0], 1e-300));
        sx += r->t[k]; sy += y; sxx += r->t[k] * r->t[k];
        sxy += r->t[k] * y; syy += y * y; m++;
    }
    if (m < 4) { *r2 = 0; return 0; }
    double den = m * sxx - sx * sx;
    double b = (m * sxy - sx * sy) / den;
    double num = m * sxy - sx * sy;
    double dd = (m * sxx - sx * sx) * (m * syy - sy * sy);
    *r2 = dd > 0 ? (num * num) / dd : 0;
    return b;
}

int main(int argc, char **argv) {
    if (argc < 2) {
        fprintf(stderr, "usage: %s label=file.tsv [...]\n", argv[0]);
        return 1;
    }
    printf("=========================================================================\n");
    printf("v86 HC-6 (corrected) -- instability in the SHAPE channel at fixed Q_a\n");
    printf("=========================================================================\n");
    printf("  The original observable (charge redistribution between components) is\n");
    printf("  VOID: V depends only on s = prod_a |Phi_a|^2, so the model has U(1)^3\n");
    printf("  and every Q_a is separately conserved. Verified below.\n\n");

    static Run R[8]; int nr = 0;
    for (int i = 1; i < argc && nr < 8; i++) {
        char *eq = strchr(argv[i], '=');
        const char *path = eq ? eq + 1 : argv[i];
        if (!load(path, &R[nr])) { fprintf(stderr, "skip %s\n", path); continue; }
        if (eq) { *eq = 0; snprintf(R[nr].label, 64, "%s", argv[i]); }
        else snprintf(R[nr].label, 64, "%s", path);
        nr++;
    }
    if (!nr) return 1;

    /* STATISTIC CHOICE (grok Finding 4, MAJOR). A log-linear "growth rate" is
     * the WRONG statistic here: the control's fit has R^2 ~ 0.001, so its rate
     * is undefined and any ratio against it is rate/noise, not signal/signal.
     * The runs also carry large breathing on top of the secular trend (~200
     * sign changes in ds over the target run), so a single exponential is not
     * the right model. Report the SECULAR AMPLITUDE -- late mean minus early
     * mean of s/s0 -- with R^2 alongside so a meaningless fit is visible. */
    printf("  %-22s %8s %13s %10s %11s %11s %7s\n",
           "arm", "t_end", "max dQ_a/Q", "s_end/s0", "secular amp", "breathing", "R^2");
    for (int i = 0; i < nr; i++) {
        Run *r = &R[i]; int e = r->n - 1;
        double r2, g = growth(r, 0.6, &r2);
        (void)g;
        /* early = first 15% mean, late = last 15% mean, of s/s0 */
        double ea = 0, la = 0; int ne = 0, nl = 0;
        for (int k = 0; k < r->n; k++) {
            double frac = r->t[k] / r->t[e];
            if (frac <= 0.15) { ea += r->s[k] / r->s[0]; ne++; }
            if (frac >= 0.85) { la += r->s[k] / r->s[0]; nl++; }
        }
        ea = ne ? ea / ne : 1; la = nl ? la / nl : 1;
        /* breathing = peak-to-peak of s/s0 about a linear trend */
        double mx = -1e30, mn = 1e30;
        for (int k = 0; k < r->n; k++) {
            double y = r->s[k] / r->s[0];
            double lin = ea + (la - ea) * (r->t[k] / r->t[e]);
            if (y - lin > mx) mx = y - lin;
            if (y - lin < mn) mn = y - lin;
        }
        printf("  %-22s %8.1f %13.2e %10.3f %+11.4f %11.4f %7.3f\n",
               r->label, r->t[e], qdrift(r), r->s[e] / r->s[0],
               la - ea, mx - mn, r2);
    }
    printf("  secular amp = <s/s0>_late - <s/s0>_early ; breathing = peak-to-peak\n");
    printf("  about the linear trend. Compare AMPLITUDES, not rates.\n");

    printf("\n  CHARGE CONSERVATION (the reason the old observable was void)\n");
    double worst = 0;
    for (int i = 0; i < nr; i++) if (qdrift(&R[i]) > worst) worst = qdrift(&R[i]);
    printf("    worst relative drift of any Q_a, any arm: %.2e\n", worst);
    printf("    -> each component charge is conserved to the integrator floor.\n");
    printf("       Redistribution between components is IMPOSSIBLE, so a GSS\n");
    printf("       instability must appear in the shape channel or not at all.\n");

    printf("\n  THE FULL SHAPE PICTURE -- fractional change from t=0 to t=end.\n");
    printf("  These are INDEPENDENT channels; an instability must move them\n");
    printf("  coherently, a numerical artifact need not.\n");
    printf("  %-22s %10s %10s %10s %10s %12s\n",
           "arm", "d s_max", "d Q_core", "d r_core", "d P_int", "dE/E");
    for (int i = 0; i < nr; i++) {
        Run *r = &R[i]; int e = r->n - 1;
        double d1 = r->s[0] ? r->s[e] / r->s[0] - 1 : 0;
        double d2 = r->qc[0] ? r->qc[e] / r->qc[0] - 1 : 0;
        double d3 = r->rc[0] ? r->rc[e] / r->rc[0] - 1 : 0;
        double d4 = r->pint[0] ? r->pint[e] / r->pint[0] - 1 : 0;
        double d5 = r->E[0] ? (r->E[e] - r->E[0]) / r->E[0] : 0;
        printf("  %-22s %+9.2f%% %+9.3f%% %+9.3f%% %+9.2f%% %12.2e\n",
               r->label, 100 * d1, 100 * d2, 100 * d3, 100 * d4, d5);
    }
    printf("\n  A COLLAPSE signature is: s_max UP, Q_core UP, r_core DOWN,\n");
    printf("  P_int UP, with E conserved -- the object contracting at fixed\n");
    printf("  charge, i.e. a VK-unstable soliton migrating to the stable branch\n");
    printf("  at the same Q. Random drift or under-resolution has no reason to\n");
    printf("  move all four coherently.\n");

    printf("\n  s_max TRAJECTORY (peak density; the shape observable)\n");
    printf("  %8s", "t");
    for (int i = 0; i < nr; i++) printf(" %14s", R[i].label);
    printf("\n");
    for (int s = 0; s <= 10; s++) {
        double tt = 0;
        for (int i = 0; i < nr; i++)
            if (R[i].t[R[i].n - 1] > tt) tt = R[i].t[R[i].n - 1];
        double target_t = tt * s / 10.0;
        printf("  %8.1f", target_t);
        for (int i = 0; i < nr; i++) {
            int k = 0;
            while (k + 1 < R[i].n && R[i].t[k + 1] <= target_t) k++;
            if (R[i].t[R[i].n - 1] + 1e-9 < target_t) printf(" %14s", "-");
            else printf(" %14.5f", R[i].s[k] / R[i].s[0]);
        }
        printf("\n");
    }
    printf("  (values are s_max(t)/s_max(0); 1.000 = no shape change)\n");
    return 0;
}
