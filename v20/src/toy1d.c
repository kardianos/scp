/*
 * toy1d.c — 1D three-scalar field simulator with saturating triple-product potential
 *
 * Lagrangian: L = sum_a [ (1/2)(dt phi_a)^2 - (1/2)(dx phi_a)^2 ] - V
 * Potential:  V = (mu/2) P^2 / (1 + kappa P^2),  P = phi_1 phi_2 phi_3
 *
 * Time integration: velocity Verlet (symplectic)
 * Boundary: Dirichlet (phi=0 at edges) + absorbing damping layer in outer 15%
 *
 * Compile: gcc -O3 -Wall -o toy1d v20/src/toy1d.c -lm
 * Usage:   ./toy1d [-mu -10] [-kappa 0.1] [-A 1.0] [-sigma 1.5] [-sep 0]
 *                  [-Nx 4000] [-xmax 40] [-tfinal 200] [-o v20/data]
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* ─── parameters ─── */
static double mu     = -10.0;
static double kappa  = 0.1;
static double A_init = 1.0;
static double sigma  = 1.5;
static double sep    = 0.0;
static int    Nx     = 4000;
static double xmax   = 40.0;
static double tfinal = 200.0;
static char   outdir[512] = "v20/data";

/* ─── helpers ─── */

static void parse_args(int argc, char **argv)
{
    for (int i = 1; i < argc - 1; i += 2) {
        if      (!strcmp(argv[i], "-mu"))     mu     = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-kappa"))  kappa  = atof(argv[i+1]);
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

/* Compute potential V = (mu/2) P^2 / (1 + kappa P^2) */
static double potential(double p1, double p2, double p3)
{
    double P  = p1 * p2 * p3;
    double P2 = P * P;
    return 0.5 * mu * P2 / (1.0 + kappa * P2);
}

/* Compute -dV/dphi_a.  Returns force (negative gradient). */
static double force_pot(double p1, double p2, double p3, int a)
{
    double P  = p1 * p2 * p3;
    double P2 = P * P;
    double denom = 1.0 + kappa * P2;
    double denom2 = denom * denom;

    /* dP/dphi_a */
    double dP;
    switch (a) {
        case 0: dP = p2 * p3; break;
        case 1: dP = p1 * p3; break;
        case 2: dP = p1 * p2; break;
        default: dP = 0.0;
    }

    /* dV/dphi_a = mu * P * dP / (1 + kappa P^2)^2 */
    double dV = mu * P * dP / denom2;
    return -dV;  /* force = -dV/dphi_a */
}

/* ─── main ─── */

int main(int argc, char **argv)
{
    parse_args(argc, argv);

    double dx = 2.0 * xmax / (Nx - 1);
    double dt = 0.5 * dx;
    int    Nt = (int)(tfinal / dt) + 1;

    printf("toy1d: mu=%.3f kappa=%.4f A=%.3f sigma=%.3f sep=%.3f\n",
           mu, kappa, A_init, sigma, sep);
    printf("  Nx=%d xmax=%.1f dx=%.5f dt=%.5f Nt=%d tfinal=%.1f\n",
           Nx, xmax, dx, dt, Nt, tfinal);

    /* allocate fields: phi[3][Nx], vel[3][Nx], acc[3][Nx] */
    double *phi[3], *vel[3], *acc[3];
    for (int a = 0; a < 3; a++) {
        phi[a] = calloc(Nx, sizeof(double));
        vel[a] = calloc(Nx, sizeof(double));
        acc[a] = calloc(Nx, sizeof(double));
        if (!phi[a] || !vel[a] || !acc[a]) {
            fprintf(stderr, "Allocation failed\n");
            return 1;
        }
    }

    /* precompute damping coefficients for absorbing boundary */
    double *damp = malloc(Nx * sizeof(double));
    double absorb_width = 0.15 * xmax;
    double x_inner = xmax - absorb_width;
    for (int i = 0; i < Nx; i++) {
        double x = -xmax + i * dx;
        double ax = fabs(x);
        if (ax > x_inner) {
            double frac = (ax - x_inner) / absorb_width;  /* 0 at inner edge, 1 at boundary */
            damp[i] = 1.0 - 0.95 * frac;  /* ramps from 1.0 to 0.05 */
        } else {
            damp[i] = 1.0;
        }
    }

    /* initialize: three Gaussians centered at -sep, 0, +sep */
    double centers[3] = { -sep, 0.0, sep };
    for (int a = 0; a < 3; a++) {
        for (int i = 0; i < Nx; i++) {
            double x = -xmax + i * dx;
            double r = x - centers[a];
            phi[a][i] = A_init * exp(-r * r / (2.0 * sigma * sigma));
        }
    }

    /* compute initial acceleration */
    for (int a = 0; a < 3; a++) {
        acc[a][0] = 0.0;
        acc[a][Nx-1] = 0.0;
        for (int i = 1; i < Nx - 1; i++) {
            double lapl = (phi[a][i+1] - 2.0*phi[a][i] + phi[a][i-1]) / (dx*dx);
            double fp = force_pot(phi[0][i], phi[1][i], phi[2][i], a);
            acc[a][i] = lapl + fp;
        }
    }

    /* open timeseries file */
    char tspath[600];
    snprintf(tspath, sizeof(tspath), "%s/toy1d_timeseries.tsv", outdir);
    FILE *fts = fopen(tspath, "w");
    if (!fts) { fprintf(stderr, "Cannot open %s\n", tspath); return 1; }
    fprintf(fts, "time\tphi1_peak\tphi2_peak\tphi3_peak\t"
                 "E_kin\tE_grad\tE_pot\tE_total\tE1\tE2\tE3\tf_core\n");

    int rec_every = Nt / 2000;
    if (rec_every < 1) rec_every = 1;
    int print_every = Nt / 20;
    if (print_every < 1) print_every = 1;

    /* core region: |x| < 3*sigma */
    double core_radius = 3.0 * sigma;

    /* ─── time integration (velocity Verlet) ─── */
    for (int n = 0; n <= Nt; n++) {
        double t = n * dt;

        /* diagnostics */
        int do_rec   = (n % rec_every == 0);
        int do_print = (n % print_every == 0);

        if (do_rec || do_print) {
            double E_kin_a[3] = {0}, E_grad_a[3] = {0};
            double E_pot_total = 0.0;
            double phi_peak[3] = {0};
            double E_core = 0.0, E_all = 0.0;

            for (int i = 1; i < Nx - 1; i++) {
                double x = -xmax + i * dx;

                /* per-field kinetic and gradient */
                for (int a = 0; a < 3; a++) {
                    E_kin_a[a]  += 0.5 * vel[a][i] * vel[a][i] * dx;
                    double dphi = (phi[a][i+1] - phi[a][i-1]) / (2.0 * dx);
                    E_grad_a[a] += 0.5 * dphi * dphi * dx;
                    double ap = fabs(phi[a][i]);
                    if (ap > phi_peak[a]) phi_peak[a] = ap;
                }

                /* potential */
                double Vloc = potential(phi[0][i], phi[1][i], phi[2][i]);
                E_pot_total += Vloc * dx;

                /* core energy density */
                double e_loc = Vloc;
                for (int a = 0; a < 3; a++) {
                    e_loc += 0.5 * vel[a][i] * vel[a][i];
                    double dphi = (phi[a][i+1] - phi[a][i-1]) / (2.0 * dx);
                    e_loc += 0.5 * dphi * dphi;
                }
                E_all += e_loc * dx;
                if (fabs(x) < core_radius)
                    E_core += e_loc * dx;
            }

            double E_kin  = E_kin_a[0] + E_kin_a[1] + E_kin_a[2];
            double E_grad = E_grad_a[0] + E_grad_a[1] + E_grad_a[2];
            double E_total = E_kin + E_grad + E_pot_total;
            double f_core = (E_all > 1e-30) ? E_core / E_all : 0.0;
            double Ea[3];
            for (int a = 0; a < 3; a++)
                Ea[a] = E_kin_a[a] + E_grad_a[a];

            if (do_rec) {
                fprintf(fts, "%.6f\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t"
                             "%.6e\t%.6e\t%.6e\t%.6f\n",
                        t, phi_peak[0], phi_peak[1], phi_peak[2],
                        E_kin, E_grad, E_pot_total, E_total,
                        Ea[0], Ea[1], Ea[2], f_core);
            }
            if (do_print) {
                printf("  t=%7.2f  E_tot=%10.4f  E_pot=%10.4f  "
                       "peaks=(%.4f,%.4f,%.4f)  f_core=%.4f\n",
                       t, E_total, E_pot_total,
                       phi_peak[0], phi_peak[1], phi_peak[2], f_core);
            }
        }

        if (n == Nt) break;  /* last step: diagnostics only */

        /* --- Velocity Verlet --- */

        /* 1. half-kick */
        for (int a = 0; a < 3; a++)
            for (int i = 1; i < Nx - 1; i++)
                vel[a][i] += 0.5 * dt * acc[a][i];

        /* 2. drift */
        for (int a = 0; a < 3; a++)
            for (int i = 1; i < Nx - 1; i++)
                phi[a][i] += dt * vel[a][i];

        /* 3. recompute acceleration */
        for (int a = 0; a < 3; a++) {
            for (int i = 1; i < Nx - 1; i++) {
                double lapl = (phi[a][i+1] - 2.0*phi[a][i] + phi[a][i-1]) / (dx*dx);
                double fp = force_pot(phi[0][i], phi[1][i], phi[2][i], a);
                acc[a][i] = lapl + fp;
            }
        }

        /* 4. half-kick */
        for (int a = 0; a < 3; a++)
            for (int i = 1; i < Nx - 1; i++)
                vel[a][i] += 0.5 * dt * acc[a][i];

        /* 5. absorbing boundary damping */
        for (int a = 0; a < 3; a++)
            for (int i = 1; i < Nx - 1; i++)
                vel[a][i] *= damp[i];
    }

    fclose(fts);
    printf("\nTimeseries written to %s\n", tspath);

    /* ─── final summary ─── */
    /* recompute final-state quantities */
    double E_pot_final = 0.0;
    double E_core_final = 0.0, E_all_final = 0.0;
    double P_center = 0.0;
    int i_center = Nx / 2;

    P_center = phi[0][i_center] * phi[1][i_center] * phi[2][i_center];

    for (int i = 1; i < Nx - 1; i++) {
        double x = -xmax + i * dx;
        double Vloc = potential(phi[0][i], phi[1][i], phi[2][i]);
        E_pot_final += Vloc * dx;

        double e_loc = Vloc;
        for (int a = 0; a < 3; a++) {
            e_loc += 0.5 * vel[a][i] * vel[a][i];
            double dphi = (phi[a][i+1] - phi[a][i-1]) / (2.0 * dx);
            e_loc += 0.5 * dphi * dphi;
        }
        E_all_final += e_loc * dx;
        if (fabs(x) < core_radius)
            E_core_final += e_loc * dx;
    }

    double f_core_final = (E_all_final > 1e-30) ? E_core_final / E_all_final : 0.0;

    /* initial kinetic + gradient energy (no potential contribution at t=0 with sep>>sigma) */
    double E_grad_init = 0.0;
    for (int a = 0; a < 3; a++) {
        for (int i = 1; i < Nx - 1; i++) {
            double x = -xmax + i * dx;
            double r = x - centers[a];
            double g = A_init * exp(-r * r / (2.0 * sigma * sigma));
            double dg = -r / (sigma * sigma) * g;
            E_grad_init += 0.5 * dg * dg * dx;
        }
    }

    printf("\n=== FINAL SUMMARY (t=%.1f) ===\n", tfinal);
    printf("  E_pot at t_final:   %+.6f\n", E_pot_final);
    printf("  f_core at t_final:  %.6f\n", f_core_final);
    printf("  |P| at center:      %.6e\n", fabs(P_center));
    printf("  Binding (E_pot<0):  %s\n", (E_pot_final < 0) ? "YES" : "NO");
    printf("  Persistent (f>0.5): %s\n", (f_core_final > 0.5) ? "YES" : "NO");
    if (E_pot_final < 0.0) {
        printf("  Binding energy:     %.6f (%.2f%% of E_grad_init=%.4f)\n",
               -E_pot_final, 100.0 * (-E_pot_final) / E_grad_init, E_grad_init);
    }
    printf("===\n");

    /* cleanup */
    for (int a = 0; a < 3; a++) {
        free(phi[a]); free(vel[a]); free(acc[a]);
    }
    free(damp);
    return 0;
}
