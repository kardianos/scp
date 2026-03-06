/*
 * biharm1d.c — 1D three-scalar field with biharmonic stabilizer + saturating potential
 *
 * Lagrangian density:
 *   L = sum_a [ (1/2)(dt phi_a)^2 - (1/2)(dx phi_a)^2 - (beta/2)(d^2x phi_a)^2 ] - V
 *
 * Potential:
 *   V = (mu/2) P^2 / (1 + kappa P^2),  P = phi_1 phi_2 phi_3
 *
 * EOM:
 *   d^2t phi_a = d^2x phi_a - beta * d^4x phi_a - dV/dphi_a
 *
 * Derrick scaling (1D, x -> lambda*x):
 *   E_biharm ~ lambda^3,  E_grad ~ lambda^1,  E_pot ~ lambda^{-1}
 *   -> stable minimum at E_grad + 3*E_biharm + E_pot = 0
 *
 * CFL: dt < C * dx^2 / sqrt(beta)  (biharmonic dominates at high k)
 *
 * Compile: gcc -O3 -Wall -o biharm1d v21/src/biharm1d.c -lm
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* parameters */
static double mu     = -10.0;
static double kappa  = 1.0;
static double beta   = 0.01;
static double A_init = 1.0;
static double sigma  = 1.5;
static double sep    = 0.0;
static int    Nx     = 2000;
static double xmax   = 20.0;
static double tfinal = 200.0;
static char   outdir[512] = "v21/data";

static void parse_args(int argc, char **argv)
{
    for (int i = 1; i < argc - 1; i += 2) {
        if      (!strcmp(argv[i], "-mu"))     mu     = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-kappa"))  kappa  = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-beta"))   beta   = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-A"))      A_init = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-sigma"))  sigma  = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-sep"))    sep    = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-Nx"))     Nx     = atoi(argv[i+1]);
        else if (!strcmp(argv[i], "-xmax"))   xmax   = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-tfinal")) tfinal = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-o"))      strncpy(outdir, argv[i+1], sizeof(outdir)-1);
        else { fprintf(stderr, "Unknown flag: %s\n", argv[i]); exit(1); }
    }
}

/* V = (mu/2) P^2 / (1 + kappa P^2) */
static double potential(double p1, double p2, double p3)
{
    double P  = p1 * p2 * p3;
    double P2 = P * P;
    return 0.5 * mu * P2 / (1.0 + kappa * P2);
}

/* -dV/dphi_a */
static double force_pot(double p1, double p2, double p3, int a)
{
    double P  = p1 * p2 * p3;
    double P2 = P * P;
    double denom2 = (1.0 + kappa * P2) * (1.0 + kappa * P2);

    double dP;
    switch (a) {
        case 0: dP = p2 * p3; break;
        case 1: dP = p1 * p3; break;
        case 2: dP = p1 * p2; break;
        default: dP = 0.0;
    }

    double dV = mu * P * dP / denom2;
    return -dV;
}

int main(int argc, char **argv)
{
    parse_args(argc, argv);

    double dx = 2.0 * xmax / (Nx - 1);

    /* CFL: dt < min( 0.5*dx (wave), C*dx^2/sqrt(beta) (biharm) ) */
    double dt_wave = 0.5 * dx;
    double dt_biharm = 0.25 * dx * dx / sqrt(beta);  /* conservative factor */
    double dt = (dt_wave < dt_biharm) ? dt_wave : dt_biharm;
    int    Nt = (int)(tfinal / dt) + 1;

    printf("biharm1d: mu=%.3f kappa=%.4f beta=%.6f A=%.3f sigma=%.3f sep=%.3f\n",
           mu, kappa, beta, A_init, sigma, sep);
    printf("  Nx=%d xmax=%.1f dx=%.5f dt=%.6f (wave=%.5f biharm=%.5f) Nt=%d\n",
           Nx, xmax, dx, dt, dt_wave, dt_biharm, Nt);

    /* allocate fields */
    double *phi[3], *vel[3], *acc[3];
    for (int a = 0; a < 3; a++) {
        phi[a] = calloc(Nx, sizeof(double));
        vel[a] = calloc(Nx, sizeof(double));
        acc[a] = calloc(Nx, sizeof(double));
    }

    /* absorbing boundary damping */
    double *damp = malloc(Nx * sizeof(double));
    double absorb_width = 0.15 * xmax;
    double x_inner = xmax - absorb_width;
    for (int i = 0; i < Nx; i++) {
        double x = -xmax + i * dx;
        double ax = fabs(x);
        if (ax > x_inner) {
            double frac = (ax - x_inner) / absorb_width;
            damp[i] = 1.0 - 0.95 * frac;
        } else {
            damp[i] = 1.0;
        }
    }

    /* initialize: three Gaussians */
    double centers[3] = { -sep, 0.0, sep };
    for (int a = 0; a < 3; a++)
        for (int i = 0; i < Nx; i++) {
            double x = -xmax + i * dx;
            double r = x - centers[a];
            phi[a][i] = A_init * exp(-r * r / (2.0 * sigma * sigma));
        }

    /* compute_accel: d^2x phi - beta * d^4x phi - dV/dphi */
    /* d^4x phi uses 5-point stencil: (phi[i-2] - 4phi[i-1] + 6phi[i] - 4phi[i+1] + phi[i+2]) / dx^4 */
    double dx2 = dx * dx;
    double dx4 = dx2 * dx2;

    #define COMPUTE_ACCEL() do { \
        for (int a = 0; a < 3; a++) { \
            acc[a][0] = acc[a][1] = 0.0; \
            acc[a][Nx-1] = acc[a][Nx-2] = 0.0; \
            for (int i = 2; i < Nx - 2; i++) { \
                double lapl = (phi[a][i+1] - 2.0*phi[a][i] + phi[a][i-1]) / dx2; \
                double biharm = (phi[a][i-2] - 4.0*phi[a][i-1] + 6.0*phi[a][i] \
                                - 4.0*phi[a][i+1] + phi[a][i+2]) / dx4; \
                double fp = force_pot(phi[0][i], phi[1][i], phi[2][i], a); \
                acc[a][i] = lapl - beta * biharm + fp; \
            } \
        } \
    } while(0)

    COMPUTE_ACCEL();

    /* timeseries output */
    char tspath[600];
    snprintf(tspath, sizeof(tspath), "%s/biharm1d_timeseries.tsv", outdir);
    FILE *fts = fopen(tspath, "w");
    if (!fts) { fprintf(stderr, "Cannot open %s\n", tspath); return 1; }
    fprintf(fts, "time\tphi1_peak\tphi2_peak\tphi3_peak\t"
                 "E_kin\tE_grad\tE_biharm\tE_pot\tE_total\tf_core\n");

    int rec_every = Nt / 4000;
    if (rec_every < 1) rec_every = 1;
    int print_every = Nt / 20;
    if (print_every < 1) print_every = 1;

    double core_radius = 3.0 * sigma;

    /* snapshot output at specific times */
    char snappath[600];
    snprintf(snappath, sizeof(snappath), "%s/biharm1d_snapshot.tsv", outdir);

    for (int n = 0; n <= Nt; n++) {
        double t = n * dt;

        int do_rec   = (n % rec_every == 0);
        int do_print = (n % print_every == 0);

        if (do_rec || do_print) {
            double E_kin = 0, E_grad = 0, E_biharm_tot = 0, E_pot_tot = 0;
            double E_core = 0, E_all = 0;
            double phi_peak[3] = {0};

            for (int i = 2; i < Nx - 2; i++) {
                double x = -xmax + i * dx;

                for (int a = 0; a < 3; a++) {
                    E_kin  += 0.5 * vel[a][i] * vel[a][i] * dx;
                    double dphi = (phi[a][i+1] - phi[a][i-1]) / (2.0 * dx);
                    E_grad += 0.5 * dphi * dphi * dx;
                    /* biharmonic energy: (beta/2)(d^2x phi)^2 */
                    double d2phi = (phi[a][i+1] - 2.0*phi[a][i] + phi[a][i-1]) / dx2;
                    E_biharm_tot += 0.5 * beta * d2phi * d2phi * dx;
                    double ap = fabs(phi[a][i]);
                    if (ap > phi_peak[a]) phi_peak[a] = ap;
                }

                double Vloc = potential(phi[0][i], phi[1][i], phi[2][i]);
                E_pot_tot += Vloc * dx;

                double e_loc = Vloc;
                for (int a = 0; a < 3; a++) {
                    e_loc += 0.5 * vel[a][i] * vel[a][i];
                    double dphi = (phi[a][i+1] - phi[a][i-1]) / (2.0 * dx);
                    e_loc += 0.5 * dphi * dphi;
                    double d2phi = (phi[a][i+1] - 2.0*phi[a][i] + phi[a][i-1]) / dx2;
                    e_loc += 0.5 * beta * d2phi * d2phi;
                }
                E_all += e_loc * dx;
                if (fabs(x) < core_radius)
                    E_core += e_loc * dx;
            }

            double E_total = E_kin + E_grad + E_biharm_tot + E_pot_tot;
            double f_core = (E_all > 1e-30) ? E_core / E_all : 0.0;

            if (do_rec) {
                fprintf(fts, "%.6f\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6f\n",
                        t, phi_peak[0], phi_peak[1], phi_peak[2],
                        E_kin, E_grad, E_biharm_tot, E_pot_tot, E_total, f_core);
            }
            if (do_print) {
                printf("  t=%7.2f  E_tot=%10.4f  E_grad=%8.4f  E_bih=%8.4f  E_pot=%10.4f  "
                       "peaks=(%.4f,%.4f,%.4f)  f_core=%.4f\n",
                       t, E_total, E_grad, E_biharm_tot, E_pot_tot,
                       phi_peak[0], phi_peak[1], phi_peak[2], f_core);
            }
        }

        if (n == Nt) break;

        /* velocity Verlet */
        for (int a = 0; a < 3; a++)
            for (int i = 2; i < Nx - 2; i++)
                vel[a][i] += 0.5 * dt * acc[a][i];

        for (int a = 0; a < 3; a++)
            for (int i = 2; i < Nx - 2; i++)
                phi[a][i] += dt * vel[a][i];

        COMPUTE_ACCEL();

        for (int a = 0; a < 3; a++)
            for (int i = 2; i < Nx - 2; i++)
                vel[a][i] += 0.5 * dt * acc[a][i];

        /* absorbing boundary */
        for (int a = 0; a < 3; a++)
            for (int i = 2; i < Nx - 2; i++)
                vel[a][i] *= damp[i];
    }

    fclose(fts);

    /* write final spatial snapshot */
    FILE *fsnap = fopen(snappath, "w");
    if (fsnap) {
        fprintf(fsnap, "x\tphi1\tphi2\tphi3\n");
        for (int i = 0; i < Nx; i++) {
            double x = -xmax + i * dx;
            fprintf(fsnap, "%.6f\t%.6e\t%.6e\t%.6e\n", x, phi[0][i], phi[1][i], phi[2][i]);
        }
        fclose(fsnap);
        printf("Snapshot written to %s\n", snappath);
    }

    printf("Timeseries written to %s\n", tspath);

    /* Derrick virial check */
    double E_grad_f = 0, E_bih_f = 0, E_pot_f = 0;
    for (int i = 2; i < Nx - 2; i++) {
        for (int a = 0; a < 3; a++) {
            double dphi = (phi[a][i+1] - phi[a][i-1]) / (2.0 * dx);
            E_grad_f += 0.5 * dphi * dphi * dx;
            double d2phi = (phi[a][i+1] - 2.0*phi[a][i] + phi[a][i-1]) / dx2;
            E_bih_f += 0.5 * beta * d2phi * d2phi * dx;
        }
        E_pot_f += potential(phi[0][i], phi[1][i], phi[2][i]) * dx;
    }
    printf("\n=== DERRICK VIRIAL (1D: E_grad + 3*E_biharm + E_pot = 0 for equilibrium) ===\n");
    printf("  E_grad  = %+.6f\n", E_grad_f);
    printf("  E_biharm= %+.6f  (3x = %+.6f)\n", E_bih_f, 3*E_bih_f);
    printf("  E_pot   = %+.6f\n", E_pot_f);
    printf("  Virial residual: %+.6f (should be ~0 for equilibrium)\n",
           E_grad_f + 3*E_bih_f + E_pot_f);
    printf("===\n");

    for (int a = 0; a < 3; a++) {
        free(phi[a]); free(vel[a]); free(acc[a]);
    }
    free(damp);
    return 0;
}
