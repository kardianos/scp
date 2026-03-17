/*
 * phonon1d.c — V23-D Phase 1: Inter-oscillon potential measurement
 *
 * Phase 1a: Equilibrate a single oscillon (mu=-20, kappa=20, m=1.0) from
 *   Gaussian init for t=5000. Save the equilibrated profile at breathing max.
 *
 * Phase 1b: For each separation D in {8,10,12,14,16,18,20,25,30}, init TWO
 *   oscillons from saved equilibrated profile (in-phase), evolve for t=2000,
 *   measure force F(D) from center-of-energy acceleration at early times.
 *
 * Output: D, F(D), V(D) table; Yukawa fit F = F0*exp(-D/lambda).
 *
 * Compile: gcc -O3 -Wall -o phonon1d src/phonon1d.c -lm
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* --- Parameters --- */
static double mu    = -20.0;
static double kappa = 20.0;
static double mass  = 1.0;
static double A_init = 1.0;
static double sig   = 3.0;
static char   outdir[512] = "v23/phonon/data";

/* Phase 1a: single oscillon equilibration */
static int    Nx_eq   = 4000;
static double xmax_eq = 100.0;
static double t_equil = 5000.0;

/* Phase 1b: two-oscillon measurement */
static int    Nx_2   = 8000;
static double xmax_2 = 200.0;
static double t_meas = 2000.0;
static double t_fit  = 200.0;  /* fit window for force measurement */

static double sep_list[] = {8, 10, 12, 14, 16, 18, 20, 25, 30};
static int    n_sep = 9;

/* Saved equilibrated profile */
static double *eq_phi[3], *eq_vel[3];
static int    eq_Nx = 0;
static double eq_dx = 0;
static double eq_xmax = 0;
static double eq_mass_osc = 0;  /* oscillon mass (energy) */

/* -dV/dphi_a where V = (mu/2)P^2/(1+kappa P^2), P = phi1 phi2 phi3 */
static inline double force_pot(double p1, double p2, double p3, int a)
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
    return -mu * P * dP / denom2;
}

/* Compute energy density at grid point i (needs phi, vel arrays and dx) */
static double energy_density(double *phi0, double *phi1, double *phi2,
                             double *v0, double *v1, double *v2,
                             int i, int Nx, double dx)
{
    double m2 = mass * mass;
    double e = 0;
    double *phi_a[3] = {phi0, phi1, phi2};
    double *vel_a[3] = {v0, v1, v2};

    for (int a = 0; a < 3; a++) {
        e += 0.5 * vel_a[a][i] * vel_a[a][i];
        if (i > 0 && i < Nx - 1) {
            double dp = (phi_a[a][i+1] - phi_a[a][i-1]) / (2.0 * dx);
            e += 0.5 * dp * dp;
        }
        e += 0.5 * m2 * phi_a[a][i] * phi_a[a][i];
    }

    double P = phi0[i] * phi1[i] * phi2[i];
    double P2 = P * P;
    e += 0.5 * mu * P2 / (1.0 + kappa * P2);

    return e;
}

/* ===================================================================
 *  Phase 1a: Equilibrate single oscillon
 * =================================================================== */
static void run_equilibrate(void)
{
    printf("===== Phase 1a: Single oscillon equilibration =====\n");
    printf("  mu=%.0f kappa=%.0f mass=%.2f A=%.2f sigma=%.2f\n",
           mu, kappa, mass, A_init, sig);

    int Nx = Nx_eq;
    double xmax = xmax_eq;
    double dx  = 2.0 * xmax / (Nx - 1);
    double dx2 = dx * dx;
    double m2  = mass * mass;

    double kmax = M_PI / dx;
    double dt   = 0.8 * 2.0 / sqrt(kmax * kmax + m2);
    int    Nt   = (int)(t_equil / dt) + 1;

    printf("  Nx=%d xmax=%.1f dx=%.5f dt=%.6f Nt=%d tfinal=%.0f\n",
           Nx, xmax, dx, dt, Nt, t_equil);

    double *phi[3], *vel[3], *acc[3];
    for (int a = 0; a < 3; a++) {
        phi[a] = calloc(Nx, sizeof(double));
        vel[a] = calloc(Nx, sizeof(double));
        acc[a] = calloc(Nx, sizeof(double));
    }

    /* Absorbing boundary: outer 25% */
    double *damp = malloc(Nx * sizeof(double));
    double x_abs = xmax * 0.75;
    for (int i = 0; i < Nx; i++) {
        double x = -xmax + i * dx;
        double ax = fabs(x);
        if (ax > x_abs) {
            double f = (ax - x_abs) / (xmax - x_abs);
            damp[i] = 1.0 - 0.98 * f * f;
        } else {
            damp[i] = 1.0;
        }
    }

    /* Initialize: symmetric Gaussians */
    int ic = Nx / 2;
    for (int a = 0; a < 3; a++)
        for (int i = 0; i < Nx; i++) {
            double x = -xmax + i * dx;
            phi[a][i] = A_init * exp(-x * x / (2.0 * sig * sig));
        }

    #define COMPUTE_ACC_EQ() do { \
        for (int a = 0; a < 3; a++) { \
            acc[a][0] = acc[a][1] = acc[a][Nx-2] = acc[a][Nx-1] = 0; \
            for (int ii = 1; ii < Nx - 1; ii++) { \
                double lapl = (phi[a][ii+1] - 2.0*phi[a][ii] + phi[a][ii-1]) / dx2; \
                double fp = force_pot(phi[0][ii], phi[1][ii], phi[2][ii], a); \
                acc[a][ii] = lapl - m2*phi[a][ii] + fp; \
            } \
        } \
    } while(0)

    COMPUTE_ACC_EQ();

    /* Track phi(0) to find breathing maxima */
    double prev_phi0 = phi[0][ic];
    double prev_prev_phi0 = phi[0][ic];
    int print_every = Nt / 40;
    if (print_every < 1) print_every = 1;

    /* We want to save the profile at the LAST breathing maximum before t_equil.
     * A breathing maximum is when phi(0) transitions from increasing to decreasing
     * AND phi(0) > 0 (positive phase of breathing). */
    double *save_phi[3], *save_vel[3];
    for (int a = 0; a < 3; a++) {
        save_phi[a] = calloc(Nx, sizeof(double));
        save_vel[a] = calloc(Nx, sizeof(double));
    }
    double save_t = -1;
    int n_maxima = 0;

    /* Only start looking for maxima after t > 0.5*t_equil (well equilibrated) */
    double t_start_save = 0.5 * t_equil;

    for (int n = 0; n <= Nt; n++) {
        double t = n * dt;

        /* Check for breathing maximum: prev > prev_prev AND prev > current, prev > 0 */
        if (n > 2 && t > t_start_save) {
            double cur = phi[0][ic];
            if (prev_phi0 > prev_prev_phi0 && prev_phi0 > cur && prev_phi0 > 0.1) {
                /* Save profile from previous step (we're one step past the max,
                 * but close enough — the breathing oscillation is smooth) */
                for (int a = 0; a < 3; a++) {
                    memcpy(save_phi[a], phi[a], Nx * sizeof(double));
                    memcpy(save_vel[a], vel[a], Nx * sizeof(double));
                }
                save_t = t - dt;
                n_maxima++;
            }
        }

        if (n % print_every == 0) {
            /* Compute energy */
            double Etot = 0;
            for (int i = 1; i < Nx - 1; i++) {
                Etot += energy_density(phi[0], phi[1], phi[2],
                                       vel[0], vel[1], vel[2], i, Nx, dx) * dx;
            }
            printf("  t=%7.1f  phi0=%+.4f  E=%.4f  maxima_found=%d\n",
                   t, phi[0][ic], Etot, n_maxima);
        }

        prev_prev_phi0 = prev_phi0;
        prev_phi0 = phi[0][ic];

        if (n == Nt) break;

        /* Velocity Verlet */
        for (int a = 0; a < 3; a++)
            for (int i = 1; i < Nx - 1; i++)
                vel[a][i] += 0.5 * dt * acc[a][i];
        for (int a = 0; a < 3; a++)
            for (int i = 1; i < Nx - 1; i++)
                phi[a][i] += dt * vel[a][i];
        COMPUTE_ACC_EQ();
        for (int a = 0; a < 3; a++)
            for (int i = 1; i < Nx - 1; i++)
                vel[a][i] += 0.5 * dt * acc[a][i];

        /* Absorbing boundary */
        for (int a = 0; a < 3; a++)
            for (int i = 0; i < Nx; i++) {
                vel[a][i] *= damp[i];
                phi[a][i] *= damp[i];
            }
    }

    if (n_maxima == 0) {
        printf("ERROR: No breathing maxima found! Using final profile.\n");
        for (int a = 0; a < 3; a++) {
            memcpy(save_phi[a], phi[a], Nx * sizeof(double));
            memcpy(save_vel[a], vel[a], Nx * sizeof(double));
        }
        save_t = t_equil;
    }

    printf("  Saved profile at t=%.1f (%d maxima found in second half)\n",
           save_t, n_maxima);

    /* Compute oscillon mass from saved profile */
    eq_mass_osc = 0;
    for (int i = 1; i < Nx - 1; i++) {
        eq_mass_osc += energy_density(save_phi[0], save_phi[1], save_phi[2],
                                       save_vel[0], save_vel[1], save_vel[2],
                                       i, Nx, dx) * dx;
    }
    printf("  Oscillon mass (energy) = %.4f\n", eq_mass_osc);

    /* Save profile to file */
    {
        char path[600];
        snprintf(path, sizeof(path), "%s/equilibrated_profile.tsv", outdir);
        FILE *fp = fopen(path, "w");
        if (fp) {
            fprintf(fp, "x\tphi1\tphi2\tphi3\tvel1\tvel2\tvel3\n");
            for (int i = 0; i < Nx; i++) {
                double x = -xmax + i * dx;
                fprintf(fp, "%.6f\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\n",
                        x, save_phi[0][i], save_phi[1][i], save_phi[2][i],
                        save_vel[0][i], save_vel[1][i], save_vel[2][i]);
            }
            fclose(fp);
            printf("  Profile saved to %s\n", path);
        }
    }

    /* Store in global arrays for Phase 1b */
    eq_Nx = Nx;
    eq_dx = dx;
    eq_xmax = xmax;
    for (int a = 0; a < 3; a++) {
        eq_phi[a] = save_phi[a];
        eq_vel[a] = save_vel[a];
    }

    /* Cleanup (don't free save_phi/save_vel — stored globally) */
    for (int a = 0; a < 3; a++) { free(phi[a]); free(vel[a]); free(acc[a]); }
    free(damp);

    printf("  Phase 1a complete.\n\n");
}

/* ===================================================================
 *  Phase 1b: Two-oscillon force measurement
 * =================================================================== */

/* Find center-of-energy for x > x_split (right oscillon) */
static double centroid_energy_right(double *phi0, double *phi1, double *phi2,
                                     double *v0, double *v1, double *v2,
                                     int Nx, double dx, double xmax, double x_split)
{
    double num = 0, den = 0;
    for (int i = 1; i < Nx - 1; i++) {
        double x = -xmax + i * dx;
        if (x > x_split) {
            double e = energy_density(phi0, phi1, phi2, v0, v1, v2, i, Nx, dx);
            if (e > 0) {
                num += x * e;
                den += e;
            }
        }
    }
    return (den > 1e-30) ? num / den : 0;
}

/* Find center-of-energy for x < x_split (left oscillon) */
static double centroid_energy_left(double *phi0, double *phi1, double *phi2,
                                    double *v0, double *v1, double *v2,
                                    int Nx, double dx, double xmax, double x_split)
{
    double num = 0, den = 0;
    for (int i = 1; i < Nx - 1; i++) {
        double x = -xmax + i * dx;
        if (x < x_split) {
            double e = energy_density(phi0, phi1, phi2, v0, v1, v2, i, Nx, dx);
            if (e > 0) {
                num += x * e;
                den += e;
            }
        }
    }
    return (den > 1e-30) ? num / den : 0;
}

/* Initialize two-oscillon field from equilibrated single-oscillon profile.
 * Places oscillons at +D/2 and -D/2. Both in-phase (same sign). */
static void init_two_oscillons(double *phi[3], double *vel[3],
                                int Nx, double dx, double xmax, double D)
{
    /* The equilibrated profile is centered at x=0 in the eq grid.
     * We interpolate it onto the two-oscillon grid shifted to +/-D/2. */
    for (int a = 0; a < 3; a++) {
        for (int i = 0; i < Nx; i++) {
            double x = -xmax + i * dx;

            /* Contribution from right oscillon at +D/2 */
            double x_r = x - D / 2.0;
            /* Contribution from left oscillon at -D/2 */
            double x_l = x + D / 2.0;

            double phi_r = 0, vel_r = 0;
            double phi_l = 0, vel_l = 0;

            /* Interpolate from equilibrated profile */
            /* eq profile: x_eq = -eq_xmax + j * eq_dx */
            double j_r = (x_r + eq_xmax) / eq_dx;
            if (j_r >= 0 && j_r < eq_Nx - 1) {
                int j0 = (int)j_r;
                double frac = j_r - j0;
                phi_r = (1 - frac) * eq_phi[a][j0] + frac * eq_phi[a][j0 + 1];
                vel_r = (1 - frac) * eq_vel[a][j0] + frac * eq_vel[a][j0 + 1];
            }

            double j_l = (x_l + eq_xmax) / eq_dx;
            if (j_l >= 0 && j_l < eq_Nx - 1) {
                int j0 = (int)j_l;
                double frac = j_l - j0;
                phi_l = (1 - frac) * eq_phi[a][j0] + frac * eq_phi[a][j0 + 1];
                vel_l = (1 - frac) * eq_vel[a][j0] + frac * eq_vel[a][j0 + 1];
            }

            phi[a][i] = phi_r + phi_l;
            vel[a][i] = vel_r + vel_l;
        }
    }
}

typedef struct {
    double D;
    double F;        /* force (positive = repulsive) */
    double a_coeff;  /* quadratic coefficient from fit */
    double D_final;  /* separation at end */
    int    merged;   /* 1 if oscillons merged */
} force_result_t;

static force_result_t measure_force(double D)
{
    force_result_t res = {0};
    res.D = D;

    int Nx = Nx_2;
    double xmax = xmax_2;
    double dx  = 2.0 * xmax / (Nx - 1);
    double dx2 = dx * dx;
    double m2  = mass * mass;

    double kmax_l = M_PI / dx;
    double dt   = 0.8 * 2.0 / sqrt(kmax_l * kmax_l + m2);
    int    Nt   = (int)(t_meas / dt) + 1;

    printf("  D=%.0f: Nx=%d xmax=%.0f dx=%.5f dt=%.6f Nt=%d\n",
           D, Nx, xmax, dx, dt, Nt);

    double *phi[3], *vel[3], *acc[3];
    for (int a = 0; a < 3; a++) {
        phi[a] = calloc(Nx, sizeof(double));
        vel[a] = calloc(Nx, sizeof(double));
        acc[a] = calloc(Nx, sizeof(double));
    }

    /* Absorbing boundary: outer 25% */
    double *damp = malloc(Nx * sizeof(double));
    double x_abs = xmax * 0.75;
    for (int i = 0; i < Nx; i++) {
        double x = -xmax + i * dx;
        double ax = fabs(x);
        if (ax > x_abs) {
            double f = (ax - x_abs) / (xmax - x_abs);
            damp[i] = 1.0 - 0.98 * f * f;
        } else {
            damp[i] = 1.0;
        }
    }

    /* Initialize from equilibrated profile */
    init_two_oscillons(phi, vel, Nx, dx, xmax, D);

    /* Compute initial acceleration */
    for (int a = 0; a < 3; a++) {
        acc[a][0] = acc[a][1] = acc[a][Nx-2] = acc[a][Nx-1] = 0;
        for (int i = 1; i < Nx - 1; i++) {
            double lapl = (phi[a][i+1] - 2.0*phi[a][i] + phi[a][i-1]) / dx2;
            double fp = force_pot(phi[0][i], phi[1][i], phi[2][i], a);
            acc[a][i] = lapl - m2*phi[a][i] + fp;
        }
    }

    /* Track separation vs time */
    int max_track = 20000;
    double *sep_hist = malloc(max_track * sizeof(double));
    double *t_track  = malloc(max_track * sizeof(double));
    int n_track = 0;
    int track_every = Nt / max_track;
    if (track_every < 1) track_every = 1;

    /* Output time series */
    char tspath[600];
    snprintf(tspath, sizeof(tspath), "%s/force_D%.0f_ts.tsv", outdir, D);
    FILE *fts = fopen(tspath, "w");
    if (fts) fprintf(fts, "time\tseparation\tx_right\tx_left\tE_total\n");

    int print_every = Nt / 20;
    if (print_every < 1) print_every = 1;

    for (int n = 0; n <= Nt; n++) {
        double t = n * dt;

        if (n % track_every == 0 && n_track < max_track) {
            double xr = centroid_energy_right(phi[0], phi[1], phi[2],
                                               vel[0], vel[1], vel[2],
                                               Nx, dx, xmax, 0.0);
            double xl = centroid_energy_left(phi[0], phi[1], phi[2],
                                              vel[0], vel[1], vel[2],
                                              Nx, dx, xmax, 0.0);
            double sep = xr - xl;

            sep_hist[n_track] = sep;
            t_track[n_track]  = t;
            n_track++;

            if (fts) {
                /* Quick energy */
                double Etot = 0;
                for (int i = 1; i < Nx - 1; i++)
                    Etot += energy_density(phi[0], phi[1], phi[2],
                                           vel[0], vel[1], vel[2], i, Nx, dx) * dx;
                fprintf(fts, "%.6f\t%.6f\t%.6f\t%.6f\t%.6e\n",
                        t, sep, xr, xl, Etot);
            }

            if (n % print_every == 0)
                printf("    t=%7.1f  sep=%.4f  xR=%.3f  xL=%.3f\n",
                       t, sep, xr, xl);
        }

        if (n == Nt) break;

        /* Velocity Verlet */
        for (int a = 0; a < 3; a++)
            for (int i = 1; i < Nx - 1; i++)
                vel[a][i] += 0.5 * dt * acc[a][i];
        for (int a = 0; a < 3; a++)
            for (int i = 1; i < Nx - 1; i++)
                phi[a][i] += dt * vel[a][i];

        for (int a = 0; a < 3; a++) {
            acc[a][0] = acc[a][1] = acc[a][Nx-2] = acc[a][Nx-1] = 0;
            for (int i = 1; i < Nx - 1; i++) {
                double lapl = (phi[a][i+1] - 2.0*phi[a][i] + phi[a][i-1]) / dx2;
                double fp = force_pot(phi[0][i], phi[1][i], phi[2][i], a);
                acc[a][i] = lapl - m2*phi[a][i] + fp;
            }
        }

        for (int a = 0; a < 3; a++)
            for (int i = 1; i < Nx - 1; i++)
                vel[a][i] += 0.5 * dt * acc[a][i];

        for (int a = 0; a < 3; a++)
            for (int i = 0; i < Nx; i++) {
                vel[a][i] *= damp[i];
                phi[a][i] *= damp[i];
            }
    }

    if (fts) fclose(fts);

    /* Measure force from quadratic fit to sep(t) over t in [0, t_fit].
     * sep(t) = s0 + v0*t + (1/2)*a*t^2
     * F = M_osc * a  (each oscillon feels F, separation acc = 2*a_each = a_sep)
     * So F_per_oscillon = (M/2) * a_sep */

    /* Find how many track points are in [0, t_fit] */
    int n_fit = 0;
    for (int j = 0; j < n_track; j++) {
        if (t_track[j] <= t_fit) n_fit = j + 1;
    }
    if (n_fit < 10) n_fit = n_track < 10 ? n_track : 10;

    if (n_fit >= 5) {
        /* Quadratic fit: sep(t) = a0 + a1*t + a2*t^2 */
        double S[5] = {0};
        double Sy[3] = {0};
        for (int j = 0; j < n_fit; j++) {
            double tj = t_track[j];
            double sj = sep_hist[j];
            double tk = 1.0;
            for (int k = 0; k < 5; k++) { S[k] += tk; tk *= tj; }
            Sy[0] += sj;
            Sy[1] += sj * tj;
            Sy[2] += sj * tj * tj;
        }

        /* 3x3 normal equations */
        double M[3][4] = {
            {S[0], S[1], S[2], Sy[0]},
            {S[1], S[2], S[3], Sy[1]},
            {S[2], S[3], S[4], Sy[2]}
        };

        /* Gaussian elimination */
        for (int col = 0; col < 3; col++) {
            int piv = col;
            for (int r = col+1; r < 3; r++)
                if (fabs(M[r][col]) > fabs(M[piv][col])) piv = r;
            if (piv != col)
                for (int c = 0; c < 4; c++) {
                    double tmp = M[col][c]; M[col][c] = M[piv][c]; M[piv][c] = tmp;
                }
            if (fabs(M[col][col]) < 1e-30) continue;
            for (int r = col+1; r < 3; r++) {
                double fac = M[r][col] / M[col][col];
                for (int c = col; c < 4; c++)
                    M[r][c] -= fac * M[col][c];
            }
        }
        double coeff[3] = {0};
        for (int r = 2; r >= 0; r--) {
            if (fabs(M[r][r]) < 1e-30) continue;
            coeff[r] = M[r][3];
            for (int c = r+1; c < 3; c++)
                coeff[r] -= M[r][c] * coeff[c];
            coeff[r] /= M[r][r];
        }

        /* a_sep = 2*coeff[2], F_per_oscillon = (M_osc/2)*a_sep = M_osc*coeff[2] */
        res.a_coeff = 2.0 * coeff[2];  /* separation acceleration */
        res.F = eq_mass_osc * coeff[2]; /* force on each oscillon */

        printf("    Fit: s0=%.4f v0=%.4e a=%.4e -> F=%.4e\n",
               coeff[0], coeff[1], 2.0*coeff[2], res.F);
    }

    res.D_final = (n_track > 0) ? sep_hist[n_track - 1] : D;
    res.merged = (res.D_final < 2.0);

    for (int a = 0; a < 3; a++) { free(phi[a]); free(vel[a]); free(acc[a]); }
    free(damp); free(sep_hist); free(t_track);

    return res;
}

/* =================================================================== */
int main(int argc, char **argv)
{
    /* Simple arg parsing */
    for (int i = 1; i < argc - 1; i += 2) {
        if      (!strcmp(argv[i], "-mu"))     mu     = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-kappa"))  kappa  = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-mass"))   mass   = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-A"))      A_init = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-sigma"))  sig    = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-o"))      strncpy(outdir, argv[i+1], sizeof(outdir)-1);
        else if (!strcmp(argv[i], "-t_equil")) t_equil = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-t_meas"))  t_meas  = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-t_fit"))    t_fit   = atof(argv[i+1]);
        else { fprintf(stderr, "Unknown flag: %s\n", argv[i]); exit(1); }
    }

    printf("=== V23-D Phase 1: Inter-Oscillon Potential ===\n\n");

    /* ===== Phase 1a: Equilibrate ===== */
    run_equilibrate();

    /* ===== Phase 1b: Force measurement ===== */
    printf("===== Phase 1b: Two-oscillon force measurement =====\n");
    printf("  Oscillon mass M = %.4f\n", eq_mass_osc);
    printf("  Fit window: t in [0, %.0f]\n", t_fit);
    printf("  Separations:");
    for (int s = 0; s < n_sep; s++) printf(" %.0f", sep_list[s]);
    printf("\n\n");

    force_result_t results[20];

    for (int s = 0; s < n_sep; s++) {
        printf("\n--- D = %.0f ---\n", sep_list[s]);
        results[s] = measure_force(sep_list[s]);
    }

    /* ===== Summary table ===== */
    printf("\n\n===== FORCE TABLE =====\n");
    printf("%-6s  %-12s  %-12s  %-8s  %-8s\n",
           "D", "F(D)", "a_sep", "D_final", "merged");
    printf("------  ------------  ------------  --------  --------\n");

    for (int s = 0; s < n_sep; s++) {
        force_result_t *r = &results[s];
        printf("%-6.0f  %+12.4e  %+12.4e  %-8.3f  %-8s\n",
               r->D, r->F, r->a_coeff, r->D_final,
               r->merged ? "YES" : "NO");
    }

    /* Save force table to file */
    {
        char path[600];
        snprintf(path, sizeof(path), "%s/force_table.tsv", outdir);
        FILE *fp = fopen(path, "w");
        if (fp) {
            fprintf(fp, "D\tF\ta_sep\tD_final\tmerged\n");
            for (int s = 0; s < n_sep; s++) {
                force_result_t *r = &results[s];
                fprintf(fp, "%.1f\t%.6e\t%.6e\t%.4f\t%d\n",
                        r->D, r->F, r->a_coeff, r->D_final, r->merged);
            }
            fclose(fp);
            printf("\nForce table saved to %s\n", path);
        }
    }

    /* ===== Yukawa fit: F = F0 * exp(-D/lambda) ===== */
    /* Use only non-merged results with |F| > 0 */
    printf("\n===== YUKAWA FIT =====\n");
    {
        double sx = 0, sy = 0, sxx = 0, sxy = 0;
        int nfit = 0;
        int sign_positive = 0, sign_negative = 0;

        for (int s = 0; s < n_sep; s++) {
            force_result_t *r = &results[s];
            if (r->merged) continue;
            if (fabs(r->F) < 1e-20) continue;

            double lnF = log(fabs(r->F));
            sx  += r->D;
            sy  += lnF;
            sxx += r->D * r->D;
            sxy += r->D * lnF;
            nfit++;

            if (r->F > 0) sign_positive++;
            else sign_negative++;
        }

        if (nfit >= 2) {
            double det = nfit * sxx - sx * sx;
            if (fabs(det) > 1e-30) {
                double slope = (nfit * sxy - sx * sy) / det;
                double intercept = (sy - slope * sx) / nfit;
                double lambda = -1.0 / slope;
                double F0 = exp(intercept);

                printf("  F(D) ~ F0 * exp(-D/lambda)\n");
                printf("  F0     = %.4e\n", F0);
                printf("  lambda = %.4f\n", lambda);
                printf("  Fit points: %d\n", nfit);

                /* Print fit comparison */
                printf("\n  %-6s  %-12s  %-12s  %-10s\n",
                       "D", "F_meas", "F_fit", "ratio");
                printf("  ------  ------------  ------------  ----------\n");
                for (int s = 0; s < n_sep; s++) {
                    force_result_t *r = &results[s];
                    if (r->merged) continue;
                    double Ffit = F0 * exp(-r->D / lambda);
                    if (r->F < 0) Ffit = -Ffit;  /* sign */
                    printf("  %-6.0f  %+12.4e  %+12.4e  %-10.4f\n",
                           r->D, r->F, Ffit,
                           (fabs(r->F) > 1e-20) ? fabs(Ffit / r->F) : 0.0);
                }
            }
        }

        /* Check for sign change (equilibrium spacing) */
        printf("\n===== EQUILIBRIUM SPACING =====\n");
        if (sign_positive > 0 && sign_negative > 0) {
            printf("  SIGN CHANGE DETECTED: F changes from attractive to repulsive!\n");
            /* Find where */
            for (int s = 0; s < n_sep - 1; s++) {
                if (results[s].merged || results[s+1].merged) continue;
                if (results[s].F * results[s+1].F < 0) {
                    /* Linear interpolation for zero crossing */
                    double D1 = results[s].D, F1 = results[s].F;
                    double D2 = results[s+1].D, F2 = results[s+1].F;
                    double d_eq = D1 - F1 * (D2 - D1) / (F2 - F1);
                    printf("  Equilibrium at D_eq ~ %.2f (between D=%.0f and D=%.0f)\n",
                           d_eq, D1, D2);
                    printf("  F(%.0f) = %+.4e, F(%.0f) = %+.4e\n",
                           D1, F1, D2, F2);
                }
            }
        } else if (sign_positive > 0) {
            printf("  All forces REPULSIVE (positive) -- no equilibrium.\n");
        } else if (sign_negative > 0) {
            printf("  All forces ATTRACTIVE (negative) -- oscillons always attract.\n");
            printf("  No equilibrium spacing -> no lattice -> no phonons.\n");
        } else {
            printf("  No valid force measurements.\n");
        }
    }

    /* ===== Integrate F(D) to get potential V(D) ===== */
    printf("\n===== POTENTIAL V(D) (integrated from infinity) =====\n");
    {
        char path[600];
        snprintf(path, sizeof(path), "%s/potential.tsv", outdir);
        FILE *fp = fopen(path, "w");
        if (fp) fprintf(fp, "D\tF\tV\n");

        printf("  %-6s  %-12s  %-12s\n", "D", "F(D)", "V(D)");
        printf("  ------  ------------  ------------\n");

        /* V(D) = -integral from D to infinity of F(D') dD'
         * Since F ~ F0*exp(-D/lambda), V(D) = -lambda*F(D)
         * But let's do numerical trapezoid from largest D down. */

        /* Sort results by D (they already are) and integrate from large D.
         * V(D_max) = 0 (reference). V(D) = V(D+dD) + F_avg * dD */
        double V_vals[20];
        V_vals[n_sep - 1] = 0;  /* reference */

        for (int s = n_sep - 2; s >= 0; s--) {
            if (results[s].merged || results[s+1].merged) {
                V_vals[s] = V_vals[s+1];
                continue;
            }
            double dD = results[s+1].D - results[s].D;
            double F_avg = 0.5 * (results[s].F + results[s+1].F);
            /* V(D) = V(D+dD) - integral_D^{D+dD} F dD' = V(D+dD) - F_avg*dD
             * (force is -dV/dD, so V(D) = V(D+dD) + integral of (-F) from D to D+dD)
             * Wait: F = -dV/dD -> dV = -F*dD
             * V(D) = V(D_max) - integral_D_max^D (-F) dD' = integral_D^D_max (-F) dD'
             * Going from large D down: V(D) = V(D+dD) + (-F_avg)*dD */
            V_vals[s] = V_vals[s+1] + (-F_avg) * dD;
        }

        for (int s = 0; s < n_sep; s++) {
            printf("  %-6.0f  %+12.4e  %+12.4e\n",
                   results[s].D, results[s].F, V_vals[s]);
            if (fp && !results[s].merged)
                fprintf(fp, "%.1f\t%.6e\t%.6e\n", results[s].D, results[s].F, V_vals[s]);
        }
        if (fp) fclose(fp);
        printf("  Potential saved to %s\n", path);
    }

    /* Cleanup global arrays */
    for (int a = 0; a < 3; a++) {
        if (eq_phi[a]) free(eq_phi[a]);
        if (eq_vel[a]) free(eq_vel[a]);
    }

    printf("\n=== V23-D Phase 1 Complete ===\n");
    return 0;
}
