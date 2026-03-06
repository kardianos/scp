/*
 * static1d.c — Gradient flow search for static 1D soliton
 *
 * Energy functional (symmetric ansatz phi_1 = phi_2 = phi_3 = phi):
 *   E = integral [ (1/2)(phi')^2 + (beta/2)(phi'')^2 + (m^2/2)phi^2
 *                  + (mu/2) phi^6 / (1 + kappa phi^6) ] dx
 *
 * The mass term m^2 phi^2/2 ensures phi -> 0 at infinity (exponential decay).
 * The phi^6 attractive potential (mu<0) creates a potential well.
 * The biharmonic (beta) prevents collapse to a delta function.
 *
 * EOM: phi'' - beta phi'''' - m^2 phi - dV/dphi = 0
 *   dV/dphi = 3 mu phi^5 / (1 + kappa phi^6)^2
 *
 * Gradient flow: d phi/dtau = phi'' - beta phi'''' - m^2 phi - dV/dphi
 *
 * Derrick virial (1D): E_grad + 3 E_biharm - E_mass - E_pot = 0
 *   (E_mass and E_pot scale as lambda^{-1} under x -> lambda x)
 *
 * Compile: gcc -O3 -Wall -o static1d v21/src/static1d.c -lm
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static double mu     = -50.0;
static double kappa  = 1.0;
static double beta   = 1.0;
static double mass   = 1.0;
static double A_init = 2.0;
static double sigma  = 1.5;
static int    Nx     = 800;
static double xmax   = 20.0;
static int    niter  = 100000000;
static char   outdir[512] = "v21/data";

static void parse_args(int argc, char **argv)
{
    for (int i = 1; i < argc - 1; i += 2) {
        if      (!strcmp(argv[i], "-mu"))     mu     = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-kappa"))  kappa  = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-beta"))   beta   = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-mass"))   mass   = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-A"))      A_init = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-sigma"))  sigma  = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-Nx"))     Nx     = atoi(argv[i+1]);
        else if (!strcmp(argv[i], "-xmax"))   xmax   = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-niter"))  niter  = atoi(argv[i+1]);
        else if (!strcmp(argv[i], "-o"))      strncpy(outdir, argv[i+1], sizeof(outdir)-1);
        else { fprintf(stderr, "Unknown flag: %s\n", argv[i]); exit(1); }
    }
}

int main(int argc, char **argv)
{
    parse_args(argc, argv);

    double dx = 2.0 * xmax / (Nx - 1);
    double dx2 = dx * dx;
    double dx4 = dx2 * dx2;
    double m2 = mass * mass;

    /* CFL: dt < min(dx^2/(4+2m^2 dx^2), dx^4/(16 beta)) */
    double dt_lapl = dx2 / (4.0 + 2.0 * m2 * dx2);
    double dt_biharm = (beta > 0) ? dx4 / (16.0 * beta) : 1e10;
    double dtau = 0.5 * ((dt_lapl < dt_biharm) ? dt_lapl : dt_biharm);

    printf("static1d: mu=%.3f kappa=%.4f beta=%.4f mass=%.4f A=%.3f sigma=%.3f\n",
           mu, kappa, beta, mass, A_init, sigma);
    printf("  Nx=%d xmax=%.1f dx=%.5f dtau=%.6e niter=%d\n",
           Nx, xmax, dx, dtau, niter);

    double *phi  = calloc(Nx, sizeof(double));
    double *force = calloc(Nx, sizeof(double));

    /* initialize: Gaussian */
    for (int i = 0; i < Nx; i++) {
        double x = -xmax + i * dx;
        phi[i] = A_init * exp(-x * x / (2.0 * sigma * sigma));
    }

    /* timeseries */
    char tspath[600];
    snprintf(tspath, sizeof(tspath), "%s/static1d_timeseries.tsv", outdir);
    FILE *fts = fopen(tspath, "w");
    if (!fts) { fprintf(stderr, "Cannot open %s\n", tspath); return 1; }
    fprintf(fts, "iter\tpeak\tE_grad\tE_biharm\tE_mass\tE_pot\tE_total\tvirial\tmax_force\n");

    int rec_every = niter / 10000;
    if (rec_every < 1) rec_every = 1;
    int print_every = niter / 40;
    if (print_every < 1) print_every = 1;

    double prev_E = 1e30;
    int converged = 0;

    for (int n = 0; n <= niter; n++) {
        int do_rec   = (n % rec_every == 0);
        int do_print = (n % print_every == 0) || (n == niter);

        if (do_rec || do_print) {
            double E_grad = 0, E_bih = 0, E_mass_tot = 0, E_pot = 0;
            double peak = 0, max_f = 0;

            for (int i = 2; i < Nx - 2; i++) {
                double dphi_dx = (phi[i+1] - phi[i-1]) / (2.0 * dx);
                E_grad += 0.5 * dphi_dx * dphi_dx * dx;

                double d2phi = (phi[i+1] - 2.0*phi[i] + phi[i-1]) / dx2;
                E_bih += 0.5 * beta * d2phi * d2phi * dx;

                E_mass_tot += 0.5 * m2 * phi[i] * phi[i] * dx;

                double p = phi[i];
                double p6 = p*p*p*p*p*p;
                double V = 0.5 * mu * p6 / (1.0 + kappa * p6);
                E_pot += V * dx;

                if (fabs(phi[i]) > peak) peak = fabs(phi[i]);
            }

            /* compute force */
            for (int i = 2; i < Nx - 2; i++) {
                double lapl = (phi[i+1] - 2.0*phi[i] + phi[i-1]) / dx2;
                double biharm = (phi[i-2] - 4.0*phi[i-1] + 6.0*phi[i]
                                - 4.0*phi[i+1] + phi[i+2]) / dx4;
                double p = phi[i];
                double p5 = p*p*p*p*p;
                double p6 = p5*p;
                double denom2 = (1.0 + kappa * p6) * (1.0 + kappa * p6);
                double dVdp = 3.0 * mu * p5 / denom2;
                force[i] = lapl - beta * biharm - m2 * p - dVdp;
                if (fabs(force[i]) > max_f) max_f = fabs(force[i]);
            }

            double E_total = E_grad + E_bih + E_mass_tot + E_pot;
            /* Derrick: E_grad + 3 E_bih - E_mass - E_pot = 0 */
            double virial = E_grad + 3.0 * E_bih - E_mass_tot - E_pot;

            if (do_rec) {
                fprintf(fts, "%d\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\n",
                        n, peak, E_grad, E_bih, E_mass_tot, E_pot, E_total, virial, max_f);
            }
            if (do_print) {
                printf("  iter=%9d  pk=%.4f  E=%+.4f  Eg=%.3f  Eb=%.3f  "
                       "Em=%.3f  Ep=%+.3f  vir=%+.3f  |F|=%.2e\n",
                       n, peak, E_total, E_grad, E_bih, E_mass_tot, E_pot,
                       virial, max_f);
            }

            if (n > 0 && max_f < 1e-8) {
                printf("\n*** CONVERGED at iter %d (max_force = %.2e) ***\n", n, max_f);
                converged = 1;
            }
            if (n > 0 && fabs(E_total - prev_E) < 1e-14 * (fabs(E_total)+1) && max_f < 1e-6) {
                printf("\n*** CONVERGED at iter %d (energy stalled, |F| = %.2e) ***\n", n, max_f);
                converged = 1;
            }
            if (peak < 1e-10) {
                printf("\n*** DISPERSED to vacuum at iter %d ***\n", n);
                converged = 1;
            }
            if (peak > 1e6) {
                printf("\n*** BLOWUP at iter %d (peak=%.2e) ***\n", n, peak);
                converged = 1;
            }
            prev_E = E_total;
        }

        if (converged || n == niter) break;

        /* gradient flow step */
        for (int i = 2; i < Nx - 2; i++) {
            double lapl = (phi[i+1] - 2.0*phi[i] + phi[i-1]) / dx2;
            double biharm = (phi[i-2] - 4.0*phi[i-1] + 6.0*phi[i]
                            - 4.0*phi[i+1] + phi[i+2]) / dx4;
            double p = phi[i];
            double p5 = p*p*p*p*p;
            double p6 = p5*p;
            double denom2 = (1.0 + kappa * p6) * (1.0 + kappa * p6);
            double dVdp = 3.0 * mu * p5 / denom2;
            force[i] = lapl - beta * biharm - m2 * p - dVdp;
        }
        for (int i = 2; i < Nx - 2; i++)
            phi[i] += dtau * force[i];

        phi[0] = phi[1] = 0.0;
        phi[Nx-1] = phi[Nx-2] = 0.0;
    }

    fclose(fts);

    /* write profile */
    char profpath[600];
    snprintf(profpath, sizeof(profpath), "%s/static1d_profile.tsv", outdir);
    FILE *fprof = fopen(profpath, "w");
    if (fprof) {
        fprintf(fprof, "x\tphi\n");
        for (int i = 0; i < Nx; i++) {
            double x = -xmax + i * dx;
            fprintf(fprof, "%.6f\t%.10e\n", x, phi[i]);
        }
        fclose(fprof);
        printf("Profile written to %s\n", profpath);
    }
    printf("Timeseries written to %s\n", tspath);

    free(phi); free(force);
    return 0;
}
