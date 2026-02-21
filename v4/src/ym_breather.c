/*
 * ym_breather.c — SU(2) Yang-Mills quasi-breather search
 *
 * Tests whether the gauge self-interaction [A,A] can create
 * long-lived oscillating localized configurations WITHOUT
 * any mass term. This is the purest realization of
 * "the interaction constrains it to itself."
 *
 * Physics: SU(2) gauge field in hedgehog ansatz
 *   A_i^a = epsilon_{ija} x_a/r^2 (1 - w(r,t))
 *
 * Reduced PDE: w_tt = w_rr - w(w^2 - 1) / r^2
 * Energy: E = integral [w_t^2 + w_r^2 + (w^2-1)^2/(2r^2)] dr
 * (Angular integration absorbed r^2 — no friction term.)
 *
 * Effective potential V(w) = (w^2-1)^2/(4r^2):
 *   - Degenerate minima at w = +1 and w = -1 (both vacuum)
 *   - Unstable maximum at w = 0
 *
 * BCs: w(0,t) = 1 (gauge regularity), w(R_max,t) = 1 (vacuum)
 * IC: w(r,0) = 1 - A*exp(-(r-R0)^2/sigma^2), w_t(r,0) = 0
 *
 * For A = 2, w dips from +1 to -1: a "bubble" of the other vacuum.
 * The bubble oscillates and may be self-confined by the gauge
 * nonlinearity for many oscillation periods.
 *
 * Usage:
 *   ym_breather -A 1.5 -outdir data
 *   ym_breather -scan -outdir data
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define NR_MAX 20001

static double w[NR_MAX], wdot[NR_MAX];

/*
 * Compute energy and core energy (r < r_core).
 */
static void compute_energy(int N, double h, double r_core,
                           double *E_total, double *E_core)
{
    double et = 0, ec = 0;
    for (int i = 1; i < N; i++) {
        double r = i * h;
        double wr = (w[i + 1] - w[i - 1]) / (2.0 * h);
        double wm = w[i] * w[i] - 1.0;
        double e = wdot[i] * wdot[i] + wr * wr + wm * wm / (2.0 * r * r);
        et += e * h;
        if (r < r_core) ec += e * h;
    }
    *E_total = et;
    *E_core = ec;
}

/*
 * Evolve: leapfrog scheme.
 *   wdot^{n+1/2} = wdot^{n-1/2} + dt * F(w^n)
 *   w^{n+1} = w^n + dt * wdot^{n+1/2}
 * where F_i = (w_{i+1} - 2w_i + w_{i-1})/h^2 - w_i*(w_i^2-1)/r_i^2
 */
static void step(int N, double h, double dt)
{
    /* Update wdot (half-step ahead) */
    for (int i = 1; i < N; i++) {
        double r = i * h;
        double wpp = (w[i + 1] - 2.0 * w[i] + w[i - 1]) / (h * h);
        double force = wpp - w[i] * (w[i] * w[i] - 1.0) / (r * r);
        wdot[i] += dt * force;
    }
    /* Update w */
    for (int i = 1; i < N; i++) {
        w[i] += dt * wdot[i];
    }
    /* BCs: w[0] = 1, w[N] = 1 */
}

static void run_single(double A, double R0, double sigma,
                        double rmax, int N, double tmax,
                        const char *outdir)
{
    double h = rmax / N;
    double dt = 0.4 * h; /* CFL: dt < h for wave eq */

    /* Initialize */
    for (int i = 0; i <= N; i++) {
        double r = i * h;
        double dr = r - R0;
        w[i] = 1.0 - A * exp(-dr * dr / (sigma * sigma));
        wdot[i] = 0.0;
    }
    w[0] = 1.0;
    w[N] = 1.0;

    /* Half-step for leapfrog initialization */
    for (int i = 1; i < N; i++) {
        double r = i * h;
        double wpp = (w[i + 1] - 2.0 * w[i] + w[i - 1]) / (h * h);
        double force = wpp - w[i] * (w[i] * w[i] - 1.0) / (r * r);
        wdot[i] = 0.5 * dt * force;
    }

    double r_core = 2.0 * R0; /* energy inside this radius */
    double E0, Ec0;
    compute_energy(N, h, r_core, &E0, &Ec0);

    printf("# YM quasi-breather: A=%.2f R0=%.1f sigma=%.1f\n", A, R0, sigma);
    printf("# Rmax=%.1f N=%d h=%.4f dt=%.5f\n", rmax, N, h, dt);
    printf("# E_init=%.6f E_core_init=%.6f (r<%.1f)\n", E0, Ec0, r_core);
    printf("# %10s %14s %14s %10s %12s\n",
           "t", "E_total", "E_core", "E_core/E0", "w(R0)");

    /* Time series file */
    char fname[512];
    snprintf(fname, sizeof(fname), "%s/ym_ts_A%.2f.dat", outdir, A);
    FILE *fp = fopen(fname, "w");
    if (fp) fprintf(fp, "# t E_total E_core E_frac w_R0\n");

    /* Snapshot file */
    char sname[512];
    snprintf(sname, sizeof(sname), "%s/ym_snap_A%.2f.dat", outdir, A);
    FILE *sfp = fopen(sname, "w");

    int nsteps = (int)(tmax / dt);
    int out_freq = nsteps / 2000; /* ~2000 output points */
    if (out_freq < 1) out_freq = 1;
    int snap_freq = nsteps / 20; /* ~20 snapshots */
    if (snap_freq < 1) snap_freq = 1;

    int iR0 = (int)(R0 / h); /* grid index near R0 */

    for (int n = 0; n < nsteps; n++) {
        step(N, h, dt);
        double t = (n + 1) * dt;

        if ((n + 1) % out_freq == 0) {
            double Et, Ec;
            compute_energy(N, h, r_core, &Et, &Ec);
            double frac = (E0 > 0) ? Ec / E0 : 0;
            double wR0 = (iR0 <= N) ? w[iR0] : 1.0;
            printf("  %10.4f %14.6e %14.6e %10.4f %12.6f\n",
                   t, Et, Ec, frac, wR0);
            if (fp) fprintf(fp, "%.6f %.8e %.8e %.6f %.8f\n",
                            t, Et, Ec, frac, wR0);
        }

        /* Snapshots of w(r) */
        if ((n + 1) % snap_freq == 0 && sfp) {
            double t_snap = (n + 1) * dt;
            fprintf(sfp, "# t = %.4f\n", t_snap);
            int sstep = (N > 500) ? N / 500 : 1;
            for (int i = 0; i <= N; i += sstep)
                fprintf(sfp, "%.4f %.8f\n", i * h, w[i]);
            fprintf(sfp, "\n");
        }
    }

    if (fp) { fclose(fp); printf("# Time series: %s\n", fname); }
    if (sfp) { fclose(sfp); printf("# Snapshots: %s\n", sname); }

    /* Final energy */
    double Ef, Ecf;
    compute_energy(N, h, r_core, &Ef, &Ecf);
    printf("# Final: E=%.6e E_core=%.6e E_core/E0=%.4f\n",
           Ef, Ecf, (E0 > 0) ? Ecf / E0 : 0);

    /* Assess result */
    double retention = (E0 > 0) ? Ecf / E0 : 0;
    if (retention > 0.5)
        printf("# RESULT: POSITIVE — %.0f%% energy retained in core after t=%.0f\n",
               retention * 100, tmax);
    else if (retention > 0.1)
        printf("# RESULT: PARTIAL — %.0f%% energy retained (quasi-breather)\n",
               retention * 100);
    else
        printf("# RESULT: NULL — energy dispersed (%.0f%% retained)\n",
               retention * 100);
}

static void run_scan(double R0, double sigma, double rmax, int N,
                     double tmax, const char *outdir)
{
    printf("# YM quasi-breather scan: R0=%.1f sigma=%.1f\n", R0, sigma);
    printf("# %6s %12s %12s %10s %s\n",
           "A", "E_init", "E_core_fin", "retention", "result");

    double A_vals[] = {0.5, 0.8, 1.0, 1.2, 1.5, 1.8, 2.0, 2.2, 2.5, 3.0};
    int nA = sizeof(A_vals) / sizeof(A_vals[0]);

    for (int j = 0; j < nA; j++) {
        double A = A_vals[j];
        double h = rmax / N;
        double dt = 0.4 * h;

        for (int i = 0; i <= N; i++) {
            double r = i * h;
            double dr = r - R0;
            w[i] = 1.0 - A * exp(-dr * dr / (sigma * sigma));
            wdot[i] = 0.0;
        }
        w[0] = 1.0;
        w[N] = 1.0;

        /* Half-step init */
        for (int i = 1; i < N; i++) {
            double r = i * h;
            double wpp = (w[i + 1] - 2.0 * w[i] + w[i - 1]) / (h * h);
            double force = wpp - w[i] * (w[i] * w[i] - 1.0) / (r * r);
            wdot[i] = 0.5 * dt * force;
        }

        double r_core = 2.0 * R0;
        double E0, Ec0;
        compute_energy(N, h, r_core, &E0, &Ec0);

        int nsteps = (int)(tmax / dt);
        for (int n = 0; n < nsteps; n++)
            step(N, h, dt);

        double Ef, Ecf;
        compute_energy(N, h, r_core, &Ef, &Ecf);
        double ret = (E0 > 0) ? Ecf / E0 : 0;

        const char *res = (ret > 0.5) ? "POSITIVE" :
                          (ret > 0.1) ? "PARTIAL" : "NULL";
        printf("  %6.2f %12.4f %12.4f %10.2f%% %s\n",
               A, E0, Ecf, ret * 100, res);
        fflush(stdout);
    }
}

int main(int argc, char **argv)
{
    double A = 2.0, R0 = 3.0, sigma = 1.0;
    double rmax = 40.0;
    int N = 8000;
    double tmax = 100.0;
    char outdir[256] = "data";
    int do_scan = 0;

    for (int i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "-A") && i + 1 < argc)
            A = atof(argv[++i]);
        else if (!strcmp(argv[i], "-R0") && i + 1 < argc)
            R0 = atof(argv[++i]);
        else if (!strcmp(argv[i], "-sigma") && i + 1 < argc)
            sigma = atof(argv[++i]);
        else if (!strcmp(argv[i], "-rmax") && i + 1 < argc)
            rmax = atof(argv[++i]);
        else if (!strcmp(argv[i], "-N") && i + 1 < argc) {
            N = atoi(argv[++i]);
            if (N >= NR_MAX) N = NR_MAX - 1;
        }
        else if (!strcmp(argv[i], "-tmax") && i + 1 < argc)
            tmax = atof(argv[++i]);
        else if (!strcmp(argv[i], "-scan"))
            do_scan = 1;
        else if (!strcmp(argv[i], "-outdir") && i + 1 < argc)
            strncpy(outdir, argv[++i], sizeof(outdir) - 1);
    }

    if (do_scan) {
        run_scan(R0, sigma, rmax, N, tmax, outdir);
    } else {
        run_single(A, R0, sigma, rmax, N, tmax, outdir);
    }

    return 0;
}
