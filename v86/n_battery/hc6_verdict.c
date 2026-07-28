/*  hc6_verdict.c — v86 HC-6: read the converse-decay runs and apply the
 *  pre-registered verdict, with the confound control the grok-4.5 review
 *  (Finding B1, MAJOR) demanded.
 *
 *  THE CONFOUND, AND WHY THE THIRD ARM EXISTS
 *    The target (flavoured, n(D)=0) and control (flavoured, n(D)=1) differ in
 *    THREE ways that all bias the target to look worse: its mean effective
 *    frequency sits at the VK turn (~1.482) while the control is below it; its
 *    tail reaches r=15 with only ~4 of margin to the sponge against the
 *    control's ~9.5; and it is the softer, more radiative object. So "target
 *    decays, control does not" would NOT establish a flavoured GSS effect.
 *
 *    The discriminator is a MONOCHROMATIC seed past the same VK turn, matched
 *    in charge and extent to the target (omega=1.4955: Q=115.43 vs 116.44,
 *    E=175.64 vs 177.40, extent 14.64 vs 15.0) and carrying NO flavour
 *    structure at all. Then:
 *
 *      target ~ mono, both differ from control  -> the effect is VK + box,
 *                                                  NOT flavoured GSS
 *      target differs from BOTH mono and control -> flavour is doing something
 *      nothing moves anywhere                    -> report a bound only
 *
 *  THE OBSERVABLE THAT SEPARATES THEM
 *    "Decay to the sector minimum" means charge REDISTRIBUTES between
 *    components: the partition fractions x_a = Q_a/Q_tot must move. Uniform
 *    loss of Q_tot at FROZEN fractions is radiation into the sponge, which is
 *    a different phenomenon and must not be reported as decay. This tool
 *    therefore reports dQ_tot/Q and max|dx_a| separately and never conflates
 *    them.
 *
 *  Build: gcc -O2 -o hc6_verdict hc6_verdict.c -lm
 *  Usage: hc6_verdict target.tsv control.tsv [mono.tsv]
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define MAXC 64
#define MAXR 200000

typedef struct {
    char name[64];
    int  n;
    double t[MAXR], Qt[MAXR], E[MAXR], smax[MAXR];
    double q[3][MAXR];
} Run;

static int col_of(char *hdr, const char *want) {
    char buf[4096];
    strncpy(buf, hdr, sizeof(buf) - 1);
    buf[sizeof(buf) - 1] = 0;
    int i = 0;
    for (char *p = strtok(buf, "\t\n"); p; p = strtok(NULL, "\t\n"), i++)
        if (!strcmp(p, want)) return i;
    return -1;
}

static int load(const char *path, Run *r) {
    FILE *f = fopen(path, "r");
    if (!f) { fprintf(stderr, "cannot open %s\n", path); return 0; }
    static char line[8192];
    if (!fgets(line, sizeof(line), f)) { fclose(f); return 0; }
    int ct = col_of(line, "t"),  cE = col_of(line, "E_total");
    int cs = col_of(line, "s_max");
    int c0 = col_of(line, "Q_p0"), c1 = col_of(line, "Q_p1"), c2 = col_of(line, "Q_p2");
    if (ct < 0 || c0 < 0 || c1 < 0 || c2 < 0) {
        fprintf(stderr, "%s: missing required columns\n", path);
        fclose(f); return 0;
    }
    const char *b = strrchr(path, '/');
    snprintf(r->name, sizeof(r->name), "%s", b ? b + 1 : path);
    r->n = 0;
    while (fgets(line, sizeof(line), f) && r->n < MAXR) {
        double v[MAXC]; int i = 0;
        for (char *p = strtok(line, "\t\n"); p && i < MAXC; p = strtok(NULL, "\t\n"))
            v[i++] = atof(p);
        if (i <= c2) continue;
        int k = r->n;
        r->t[k] = v[ct];
        r->E[k] = (cE >= 0 && cE < i) ? v[cE] : 0.0;
        r->smax[k] = (cs >= 0 && cs < i) ? v[cs] : 0.0;
        r->q[0][k] = v[c0]; r->q[1][k] = v[c1]; r->q[2][k] = v[c2];
        r->Qt[k] = v[c0] + v[c1] + v[c2];
        r->n++;
    }
    fclose(f);
    return r->n > 1;
}

/* max |x_a(t) - x_a(0)| over components: the redistribution signal */
static double frac_drift(const Run *r, int k) {
    double d = 0;
    for (int a = 0; a < 3; a++) {
        double x0 = r->q[a][0] / r->Qt[0], xk = r->q[a][k] / r->Qt[k];
        double e = fabs(xk - x0);
        if (e > d) d = e;
    }
    return d;
}

static void report(const Run *r) {
    int last = r->n - 1;
    printf("\n--- %s ---\n", r->name);
    printf("  t range 0 .. %.2f   (%d diag rows)\n", r->t[last], r->n);
    printf("  %8s %12s %10s %12s %12s %12s %11s %10s\n",
           "t", "Q_tot", "dQ/Q %", "x_0", "x_1", "x_2", "max|dx|", "s_max");
    int steps = 8;
    for (int s = 0; s <= steps; s++) {
        int k = (int)((double)s / steps * last);
        double x[3];
        for (int a = 0; a < 3; a++) x[a] = r->q[a][k] / r->Qt[k];
        printf("  %8.2f %12.4f %10.4f %12.6f %12.6f %12.6f %11.2e %10.3e\n",
               r->t[k], r->Qt[k], 100.0 * (r->Qt[k] - r->Qt[0]) / r->Qt[0],
               x[0], x[1], x[2], frac_drift(r, k), r->smax[k]);
    }
}

int main(int argc, char **argv) {
    if (argc < 3) {
        fprintf(stderr, "usage: %s target.tsv control.tsv [mono.tsv]\n", argv[0]);
        return 1;
    }
    static Run R[3]; int nr = 0;
    const char *lbl[3] = { "TARGET  (flavoured, n(D)=0)",
                           "CONTROL (flavoured, n(D)=1)",
                           "MONO    (monochromatic past the VK turn, n(D)=0)" };
    for (int i = 1; i < argc && i <= 3; i++)
        if (load(argv[i], &R[nr])) nr++;
    if (nr < 2) { fprintf(stderr, "need at least target and control\n"); return 1; }

    printf("================================================================\n");
    printf("v86 HC-6 -- converse decay, with the VK/box confound controlled\n");
    printf("================================================================\n");
    for (int i = 0; i < nr; i++) {
        printf("\n[%s]", lbl[i]);
        report(&R[i]);
    }

    printf("\n================================================================\n");
    printf("VERDICT\n");
    printf("================================================================\n");
    /* compare on the COMMON time window -- the runs are not equally advanced */
    double tmin = R[0].t[R[0].n - 1];
    for (int i = 1; i < nr; i++)
        if (R[i].t[R[i].n - 1] < tmin) tmin = R[i].t[R[i].n - 1];
    printf("  common window: t = 0 .. %.2f\n\n", tmin);
    printf("  %-52s %12s %12s\n", "run", "dQ_tot/Q %", "max|dx_a|");
    double dq[3], dx[3];
    for (int i = 0; i < nr; i++) {
        int k = 0;
        while (k + 1 < R[i].n && R[i].t[k + 1] <= tmin) k++;
        dq[i] = 100.0 * (R[i].Qt[k] - R[i].Qt[0]) / R[i].Qt[0];
        dx[i] = frac_drift(&R[i], k);
        printf("  %-52s %12.4f %12.3e\n", lbl[i], dq[i], dx[i]);
    }

    /* CRITICAL: the mono arm's fractions are pinned at 1/3 by EXACT permutation
     * symmetry of its initial data -- three identical components evolve
     * identically forever. Its max|dx| is therefore ~1e-17 BY CONSTRUCTION and
     * carries ZERO information about redistribution. Using it as the
     * discriminator would manufacture a conclusion. The mono arm's job is the
     * CHARGE-LOSS comparison only. */
    const double DX_SIG = 1e-3;     /* pre-registered: redistribution threshold */
    printf("\n  pre-registered redistribution threshold: max|dx_a| > %.0e\n", DX_SIG);
    int t_redist = dx[0] > DX_SIG, c_redist = dx[1] > DX_SIG;
    int m_redist = (nr > 2) ? dx[2] > DX_SIG : 0;

    if (!t_redist && !c_redist) {
        printf("\n  NEITHER target nor control redistributes charge between\n");
        printf("  components on this window. The target's charge loss is UNIFORM\n");
        printf("  (fractions frozen), which is radiation into the sponge, not\n");
        printf("  decay to a sector minimum.\n");
        printf("\n  => REPORT A BOUND, NOT A CONFIRMATION. n(D)=0 is not sufficient\n");
        printf("     for observable flavour redistribution on t <= %.0f at this\n", tmin);
        printf("     box and resolution. The GSS converse is NOT confirmed here.\n");
    } else if (t_redist && !c_redist) {
        printf("\n  Target's partition fractions move (%.2e) and the control's do\n", dx[0]);
        printf("  not (%.2e). But that comparison is NOT yet a result:\n", dx[1]);
        /* normalise the drift by how much charge was actually lost -- the
         * components sit at different omega, so they radiate at different
         * rates, and ANY object losing charge will show some fraction drift. */
        double n0 = dx[0] / fmax(fabs(dq[0]), 1e-30);
        double n1 = dx[1] / fmax(fabs(dq[1]), 1e-30);
        printf("\n  DRIFT PER UNIT CHARGE LOST (the fair comparison, since the\n");
        printf("  components sit at different omega and radiate at different rates):\n");
        printf("    target : %.3e per %% lost\n", n0);
        printf("    control: %.3e per %% lost\n", n1);
        if (n0 <= n1)
            printf("    -> the target drifts LESS per unit charge lost than the\n"
                   "       control. Its fraction drift is fully accounted for by\n"
                   "       differential radiation, NOT by redistribution to a\n"
                   "       sector minimum. NO GSS converse signal.\n");
        else
            printf("    -> the target drifts MORE per unit charge lost (%.1fx);\n"
                   "       that excess is the only candidate GSS signal.\n", n0 / n1);
    } else {
        printf("\n  Control also moves: the seeding or the box is the instability,\n");
        printf("  not the partition. Test void, redesign required.\n");
    }

    if (nr > 2) {
        printf("\n  CONFOUND CHECK (the reason the mono arm exists):\n");
        printf("    target dQ/Q = %+.4f %%   mono dQ/Q = %+.4f %%   control = %+.4f %%\n",
               dq[0], dq[2], dq[1]);
        printf("    (mono max|dx| = %.1e is ZERO BY SYMMETRY, not a measurement --\n", dx[2]);
        printf("     three identical components stay identical. Ignore it.)\n");
        if (fabs(dq[2]) >= fabs(dq[0]))
            printf("\n    -> the MONOCHROMATIC object past the same turn, matched in\n"
                   "       charge and extent, loses MORE (%.2f%% vs %.2f%%). The\n"
                   "       target's loss is therefore entirely within what 'past the\n"
                   "       VK turn, in this box' produces -- if anything the\n"
                   "       flavoured object is the MORE stable one. No flavoured\n"
                   "       destabilisation is present to explain.\n",
                   fabs(dq[2]), fabs(dq[0]));
        else
            printf("\n    -> the target loses MORE than the matched monochromatic\n"
                   "       object (%.2f%% vs %.2f%%); that excess is a candidate\n"
                   "       flavoured effect.\n", fabs(dq[0]), fabs(dq[2]));
    }
    return 0;
}
