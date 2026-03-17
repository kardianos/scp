/*
 * critical2.c — V23-F: Critical Gravity Phase 2 Redo with Equilibrated Profiles
 *
 * Two modes:
 *   -mode equil    : Equilibrate a single oscillon, save full profile at peak phase
 *   -mode interact : Load equilibrated profile, run two-oscillon interaction
 *   -mode full     : Run full protocol (equil + interact for all mu, D)
 *
 * Lagrangian:
 *   L = sum_a [ (1/2)(dt phi_a)^2 - (1/2)(dx phi_a)^2 - (m^2/2)phi_a^2 ]
 *     - (mu/2) P^2 / (1 + kappa P^2)
 *   P = phi_1 phi_2 phi_3
 *
 * Compile: gcc -O3 -Wall -o critical2 src/critical2.c -lm
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* --- Parameters --- */
static double kappa    = 20.0;
static double mass     = 1.0;
static double A_init   = 1.0;
static double sig      = 3.0;
static char   outdir[512] = "data";
static char   mode[64] = "full";
static double mu_single = -10.0;
static double D_single  = 20.0;

/* Equilibration grid */
static int    Nx_eq    = 4000;
static double xmax_eq  = 80.0;
static double t_equil  = 10000.0;

/* Interaction grid */
static int    Nx_int   = 8000;
static double xmax_int = 200.0;
static double t_run    = 5000.0;

/* Protocol scan */
static double mu_list[]  = {-20.0, -14.0, -10.0, -8.0, -6.0};
static int    n_mu       = 5;
static double D_list[]   = {15.0, 20.0, 25.0, 30.0, 40.0};
static int    n_D        = 5;

static void parse_args(int argc, char **argv)
{
    for (int i = 1; i < argc - 1; i += 2) {
        if      (!strcmp(argv[i], "-mode"))   strncpy(mode, argv[i+1], sizeof(mode)-1);
        else if (!strcmp(argv[i], "-mu"))     mu_single = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-D"))      D_single  = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-kappa"))  kappa     = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-mass"))   mass      = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-A"))      A_init    = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-sigma"))  sig       = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-o"))      strncpy(outdir, argv[i+1], sizeof(outdir)-1);
        else if (!strcmp(argv[i], "-Nx_eq"))  Nx_eq     = atoi(argv[i+1]);
        else if (!strcmp(argv[i], "-Nx_int")) Nx_int    = atoi(argv[i+1]);
        else if (!strcmp(argv[i], "-xmax_eq"))  xmax_eq  = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-xmax_int")) xmax_int = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-t_equil"))  t_equil  = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-t_run"))    t_run    = atof(argv[i+1]);
        else { fprintf(stderr, "Unknown flag: %s\n", argv[i]); exit(1); }
    }
}

/* -dV/dphi_a where V = (mu/2)P^2/(1+kappa P^2), P = phi1 phi2 phi3 */
static inline double force_pot(double p1, double p2, double p3, int a, double mu)
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

/* ===================================================================
 *  Saved profile data
 * =================================================================== */
typedef struct {
    int    N;          /* number of grid points */
    double xmax;       /* half-domain */
    double dx;
    double omega;      /* breathing frequency */
    double A_peak;     /* peak amplitude at save time */
    double E_total;    /* total energy */
    double *phi[3];    /* field values */
    double *vel[3];    /* velocity values */
} profile_t;

static void profile_free(profile_t *p)
{
    for (int a = 0; a < 3; a++) {
        if (p->phi[a]) free(p->phi[a]);
        if (p->vel[a]) free(p->vel[a]);
    }
}

static int profile_save(const char *path, const profile_t *p)
{
    FILE *f = fopen(path, "w");
    if (!f) { fprintf(stderr, "Cannot open %s for writing\n", path); return -1; }
    fprintf(f, "# N=%d xmax=%.6f dx=%.6e omega=%.6f A_peak=%.6f E_total=%.6f\n",
            p->N, p->xmax, p->dx, p->omega, p->A_peak, p->E_total);
    fprintf(f, "# x  phi1  phi2  phi3  vel1  vel2  vel3\n");
    for (int i = 0; i < p->N; i++) {
        double x = -p->xmax + i * p->dx;
        fprintf(f, "%.8e  %.12e  %.12e  %.12e  %.12e  %.12e  %.12e\n",
                x, p->phi[0][i], p->phi[1][i], p->phi[2][i],
                p->vel[0][i], p->vel[1][i], p->vel[2][i]);
    }
    fclose(f);
    return 0;
}

static int profile_load(const char *path, profile_t *p)
{
    FILE *f = fopen(path, "r");
    if (!f) { fprintf(stderr, "Cannot open %s for reading\n", path); return -1; }

    /* Read header */
    char line[1024];
    if (!fgets(line, sizeof(line), f)) { fclose(f); return -1; }
    sscanf(line, "# N=%d xmax=%lf dx=%*s omega=%lf A_peak=%lf E_total=%lf",
           &p->N, &p->xmax, &p->omega, &p->A_peak, &p->E_total);
    p->dx = 2.0 * p->xmax / (p->N - 1);

    /* Skip second header line */
    if (!fgets(line, sizeof(line), f)) { fclose(f); return -1; }

    /* Allocate */
    for (int a = 0; a < 3; a++) {
        p->phi[a] = calloc(p->N, sizeof(double));
        p->vel[a] = calloc(p->N, sizeof(double));
    }

    /* Read data */
    for (int i = 0; i < p->N; i++) {
        double x;
        if (fscanf(f, "%lf %lf %lf %lf %lf %lf %lf",
                   &x, &p->phi[0][i], &p->phi[1][i], &p->phi[2][i],
                   &p->vel[0][i], &p->vel[1][i], &p->vel[2][i]) != 7) {
            fprintf(stderr, "Read error at line %d in %s\n", i+3, path);
            fclose(f);
            return -1;
        }
    }

    fclose(f);
    return 0;
}

/* ===================================================================
 *  MODE 1: Equilibrate a single oscillon
 * =================================================================== */
static profile_t run_equilibrate(double mu, int save_files)
{
    profile_t prof = {0};
    prof.N = Nx_eq;
    prof.xmax = xmax_eq;
    prof.dx = 2.0 * xmax_eq / (Nx_eq - 1);

    double dx  = prof.dx;
    double dx2 = dx * dx;
    double m2  = mass * mass;
    int    Nx  = Nx_eq;
    int    ic  = Nx / 2;

    /* CFL */
    double kmax = M_PI / dx;
    double dt   = 0.8 * 2.0 / sqrt(kmax * kmax + m2);
    int    Nt   = (int)(t_equil / dt) + 1;

    printf("  Equilibrate: mu=%.0f Nx=%d xmax=%.0f dx=%.5f dt=%.6f Nt=%d\n",
           mu, Nx, xmax_eq, dx, dt, Nt);

    /* Allocate fields */
    double *phi[3], *vel[3], *acc[3];
    for (int a = 0; a < 3; a++) {
        phi[a] = calloc(Nx, sizeof(double));
        vel[a] = calloc(Nx, sizeof(double));
        acc[a] = calloc(Nx, sizeof(double));
    }

    /* Absorbing boundary: outer 25% */
    double *damp = malloc(Nx * sizeof(double));
    double x_abs = xmax_eq * 0.75;
    for (int i = 0; i < Nx; i++) {
        double x = -xmax_eq + i * dx;
        double ax = fabs(x);
        if (ax > x_abs) {
            double f = (ax - x_abs) / (xmax_eq - x_abs);
            damp[i] = 1.0 - 0.98 * f * f;
        } else {
            damp[i] = 1.0;
        }
    }

    /* Initialize: symmetric Gaussians */
    for (int a = 0; a < 3; a++)
        for (int i = 0; i < Nx; i++) {
            double x = -xmax_eq + i * dx;
            phi[a][i] = A_init * exp(-x * x / (2.0 * sig * sig));
        }

    /* Compute acceleration */
    #define COMPUTE_ACC_EQ() do { \
        for (int a = 0; a < 3; a++) { \
            acc[a][0] = acc[a][1] = acc[a][Nx-2] = acc[a][Nx-1] = 0; \
            for (int ii = 1; ii < Nx - 1; ii++) { \
                double lapl = (phi[a][ii+1] - 2.0*phi[a][ii] + phi[a][ii-1]) / dx2; \
                double fp = force_pot(phi[0][ii], phi[1][ii], phi[2][ii], a, mu); \
                acc[a][ii] = lapl - m2*phi[a][ii] + fp; \
            } \
        } \
    } while(0)

    COMPUTE_ACC_EQ();

    /* DFT storage for frequency measurement */
    int max_dft = 100000;
    double *phi0_hist = malloc(max_dft * sizeof(double));
    double *t_hist    = malloc(max_dft * sizeof(double));
    int n_dft = 0;
    int dft_every = Nt / max_dft;
    if (dft_every < 1) dft_every = 1;

    int print_every = Nt / 20;
    if (print_every < 1) print_every = 1;

    /* Save state at positive peaks in the last 5% of evolution.
     * We keep a rolling snapshot: overwrite every time we find a new peak.
     * This avoids expensive re-evolution. */
    double *snap_phi[3], *snap_vel[3];
    for (int a = 0; a < 3; a++) {
        snap_phi[a] = calloc(Nx, sizeof(double));
        snap_vel[a] = calloc(Nx, sizeof(double));
    }
    int have_snap = 0;
    double prev_phi0 = 0, prev_prev_phi0 = 0;
    double snap_phi0_val = 0;

    double last_E = 0;

    for (int n = 0; n <= Nt; n++) {
        double t = n * dt;

        /* Record for DFT */
        if (n % dft_every == 0 && n_dft < max_dft) {
            phi0_hist[n_dft] = phi[0][ic];
            t_hist[n_dft] = t;
            n_dft++;
        }

        /* Detect positive peaks in last 5% and save full state */
        if (t > t_equil * 0.95) {
            double cur = phi[0][ic];
            if (prev_phi0 > prev_prev_phi0 && prev_phi0 > cur && prev_phi0 > 0.01) {
                /* prev was a local maximum — save state from PREVIOUS step.
                 * Since we already advanced, save current state (close enough
                 * — one dt from peak). For cleaner results, we save at n-1
                 * but we don't have it. Instead, just note the peak value
                 * and save current state which is within 1 dt of peak. */
                for (int a = 0; a < 3; a++) {
                    memcpy(snap_phi[a], phi[a], Nx * sizeof(double));
                    memcpy(snap_vel[a], vel[a], Nx * sizeof(double));
                }
                have_snap = 1;
                snap_phi0_val = prev_phi0;
            }
            prev_prev_phi0 = prev_phi0;
            prev_phi0 = cur;
        }

        /* Energy + printing */
        if (n % print_every == 0) {
            double Ek = 0, Eg = 0, Em = 0, Ep = 0;
            double peak = 0;
            for (int i = 1; i < Nx - 1; i++) {
                for (int a = 0; a < 3; a++) {
                    Ek += 0.5 * vel[a][i] * vel[a][i] * dx;
                    double dp = (phi[a][i+1] - phi[a][i-1]) / (2.0*dx);
                    Eg += 0.5 * dp * dp * dx;
                    Em += 0.5 * m2 * phi[a][i] * phi[a][i] * dx;
                    if (fabs(phi[a][i]) > peak) peak = fabs(phi[a][i]);
                }
                double P  = phi[0][i] * phi[1][i] * phi[2][i];
                double P2 = P * P;
                Ep += 0.5 * mu * P2 / (1.0 + kappa * P2) * dx;
            }
            last_E = Ek + Eg + Em + Ep;
            printf("    t=%7.0f  phi0=%+.4f  pk=%.4f  E=%+.4f\n",
                   t, phi[0][ic], peak, last_E);
        }

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

    /* Use the saved snapshot (state at last positive peak) */
    if (have_snap) {
        printf("  Using saved peak snapshot: phi1(0)=%.6f\n", snap_phi0_val);
        for (int a = 0; a < 3; a++) {
            memcpy(phi[a], snap_phi[a], Nx * sizeof(double));
            memcpy(vel[a], snap_vel[a], Nx * sizeof(double));
        }
    } else {
        printf("  WARNING: no positive peak found in last 5%%, using final state\n");
    }

    /* Compute energy at this state */
    double E_save = 0;
    double A_save = 0;
    for (int i = 1; i < Nx - 1; i++) {
        for (int a = 0; a < 3; a++) {
            E_save += 0.5 * vel[a][i] * vel[a][i] * dx;
            double dp = (phi[a][i+1] - phi[a][i-1]) / (2.0*dx);
            E_save += 0.5 * dp * dp * dx;
            E_save += 0.5 * m2 * phi[a][i] * phi[a][i] * dx;
            if (fabs(phi[a][i]) > A_save) A_save = fabs(phi[a][i]);
        }
        double P  = phi[0][i] * phi[1][i] * phi[2][i];
        double P2 = P * P;
        E_save += 0.5 * mu * P2 / (1.0 + kappa * P2) * dx;
    }

    /* Measure frequency from DFT (second half) */
    double peak_omega = 0;
    double peak_pow_val = 0;
    int dft_start = n_dft / 2;
    if (n_dft - dft_start > 100) {
        double T = t_hist[n_dft-1] - t_hist[dft_start];
        int nf = 2000;
        for (int k = 1; k < nf; k++) {
            double omega = 2.0 * mass * k / nf;
            double re = 0, im = 0;
            for (int j = dft_start; j < n_dft; j++) {
                double dtj = (j > dft_start) ?
                    (t_hist[j]-t_hist[j-1]) : (t_hist[dft_start+1]-t_hist[dft_start]);
                re += phi0_hist[j] * cos(omega * t_hist[j]) * dtj;
                im += phi0_hist[j] * sin(omega * t_hist[j]) * dtj;
            }
            double pw = (re*re + im*im) / (T*T);
            if (pw > peak_pow_val) { peak_pow_val = pw; peak_omega = omega; }
        }
    }

    /* Fill profile struct */
    for (int a = 0; a < 3; a++) {
        prof.phi[a] = malloc(Nx * sizeof(double));
        prof.vel[a] = malloc(Nx * sizeof(double));
        memcpy(prof.phi[a], phi[a], Nx * sizeof(double));
        memcpy(prof.vel[a], vel[a], Nx * sizeof(double));
    }
    prof.omega  = peak_omega;
    prof.A_peak = A_save;
    prof.E_total = E_save;

    printf("  omega=%.4f  A_peak=%.4f  E_total=%.4f\n", peak_omega, A_save, E_save);

    /* Save profile */
    if (save_files) {
        char path[600];
        snprintf(path, sizeof(path), "%s/profile_mu%d.dat", outdir, (int)fabs(mu));
        profile_save(path, &prof);
        printf("  Saved profile to %s\n", path);

        /* Also save summary */
        snprintf(path, sizeof(path), "%s/equil_mu%d_summary.txt", outdir, (int)fabs(mu));
        FILE *fs = fopen(path, "w");
        if (fs) {
            fprintf(fs, "mu=%.2f\n", mu);
            fprintf(fs, "omega=%.6f\n", peak_omega);
            fprintf(fs, "A_peak=%.6f\n", A_save);
            fprintf(fs, "E_total=%.6f\n", E_save);
            double k2 = mass*mass - peak_omega*peak_omega;
            double xi_pred = (k2 > 0) ? 1.0/sqrt(k2) : 9999.0;
            fprintf(fs, "xi_pred=%.6f\n", xi_pred);
            fprintf(fs, "gap_margin=%.6f\n", (mass - peak_omega)/mass);
            fclose(fs);
        }
    }

    /* Cleanup */
    for (int a = 0; a < 3; a++) { free(phi[a]); free(vel[a]); free(acc[a]); }
    for (int a = 0; a < 3; a++) { free(snap_phi[a]); free(snap_vel[a]); }
    free(damp); free(phi0_hist); free(t_hist);

    return prof;
}

/* ===================================================================
 *  Interpolate profile onto interaction grid
 * =================================================================== */
static double interp_profile(const double *data, int N_src, double dx_src, double xmax_src, double x_query)
{
    /* x_query is relative to profile center (0) */
    double xi = (x_query + xmax_src) / dx_src;
    int i0 = (int)floor(xi);
    if (i0 < 0 || i0 >= N_src - 1) return 0.0;
    double frac = xi - i0;
    return data[i0] * (1.0 - frac) + data[i0 + 1] * frac;
}

/* ===================================================================
 *  Energy centroid using energy density (more robust than phi^2)
 * =================================================================== */
static void compute_energy_centroids(double *phi[], double *vel[], int Nx, double dx,
                                      double xmax_l, double m2, double mu,
                                      double *x_left, double *x_right,
                                      double *E_left, double *E_right)
{
    double num_L = 0, den_L = 0, num_R = 0, den_R = 0;
    int ic = Nx / 2;

    for (int i = 1; i < Nx - 1; i++) {
        double x = -xmax_l + i * dx;
        double e = 0;
        for (int a = 0; a < 3; a++) {
            e += 0.5 * vel[a][i] * vel[a][i];
            double dp = (phi[a][i+1] - phi[a][i-1]) / (2.0*dx);
            e += 0.5 * dp * dp;
            e += 0.5 * m2 * phi[a][i] * phi[a][i];
        }
        double P  = phi[0][i] * phi[1][i] * phi[2][i];
        double P2 = P * P;
        e += 0.5 * mu * P2 / (1.0 + kappa * P2);

        double w = fabs(e) * dx;  /* use |e| to avoid sign issues from negative V */
        if (i < ic) {
            num_L += x * w;
            den_L += w;
        } else {
            num_R += x * w;
            den_R += w;
        }
    }

    *x_left  = (den_L > 1e-30) ? num_L / den_L : 0;
    *x_right = (den_R > 1e-30) ? num_R / den_R : 0;
    *E_left  = den_L;
    *E_right = den_R;
}

/* ===================================================================
 *  MODE 2: Two-oscillon interaction
 * =================================================================== */
typedef struct {
    double F_fit;      /* fitted acceleration coefficient (force/mass = 2*c in s=a+bt+ct^2) */
    double D_final;
    int    valid;
} interact_result_t;

static interact_result_t run_interact(double mu, double D, const profile_t *prof)
{
    interact_result_t res = {0};

    int    Nx  = Nx_int;
    double xm  = xmax_int;
    double dx  = 2.0 * xm / (Nx - 1);
    double dx2 = dx * dx;
    double m2  = mass * mass;

    double kmax = M_PI / dx;
    double dt   = 0.8 * 2.0 / sqrt(kmax * kmax + m2);
    int    Nt   = (int)(t_run / dt) + 1;

    printf("    Interact: mu=%.0f D=%.0f Nx=%d xmax=%.0f dx=%.5f dt=%.6f Nt=%d\n",
           mu, D, Nx, xm, dx, dt, Nt);

    double *phi[3], *vel[3], *acc[3];
    for (int a = 0; a < 3; a++) {
        phi[a] = calloc(Nx, sizeof(double));
        vel[a] = calloc(Nx, sizeof(double));
        acc[a] = calloc(Nx, sizeof(double));
    }

    /* Absorbing boundary */
    double *damp = malloc(Nx * sizeof(double));
    double x_abs = xm * 0.75;
    for (int i = 0; i < Nx; i++) {
        double x = -xm + i * dx;
        double ax = fabs(x);
        if (ax > x_abs) {
            double f = (ax - x_abs) / (xm - x_abs);
            damp[i] = 1.0 - 0.98 * f * f;
        } else {
            damp[i] = 1.0;
        }
    }

    /* Initialize from saved profile: two copies at +/- D/2 */
    for (int a = 0; a < 3; a++) {
        for (int i = 0; i < Nx; i++) {
            double x = -xm + i * dx;
            double f1 = interp_profile(prof->phi[a], prof->N, prof->dx, prof->xmax, x - D/2);
            double f2 = interp_profile(prof->phi[a], prof->N, prof->dx, prof->xmax, x + D/2);
            phi[a][i] = f1 + f2;

            double v1 = interp_profile(prof->vel[a], prof->N, prof->dx, prof->xmax, x - D/2);
            double v2 = interp_profile(prof->vel[a], prof->N, prof->dx, prof->xmax, x + D/2);
            vel[a][i] = v1 + v2;
        }
    }

    /* Compute initial acceleration */
    #define COMPUTE_ACC_INT(NX) do { \
        for (int a = 0; a < 3; a++) { \
            acc[a][0] = acc[a][1] = acc[a][(NX)-2] = acc[a][(NX)-1] = 0; \
            for (int ii = 1; ii < (NX) - 1; ii++) { \
                double lapl = (phi[a][ii+1] - 2.0*phi[a][ii] + phi[a][ii-1]) / dx2; \
                double fp = force_pot(phi[0][ii], phi[1][ii], phi[2][ii], a, mu); \
                acc[a][ii] = lapl - m2*phi[a][ii] + fp; \
            } \
        } \
    } while(0)

    COMPUTE_ACC_INT(Nx);

    /* Track separation vs time */
    int max_track = 50000;
    double *sep_hist = malloc(max_track * sizeof(double));
    double *t_track  = malloc(max_track * sizeof(double));
    int n_track = 0;
    int track_every = Nt / max_track;
    if (track_every < 1) track_every = 1;

    /* Output TSV */
    char path[600];
    snprintf(path, sizeof(path), "%s/interact_mu%d_D%d_ts.tsv",
             outdir, (int)fabs(mu), (int)D);
    FILE *fts = fopen(path, "w");
    if (fts) fprintf(fts, "time\tseparation\tx_right\tx_left\tE_total\n");

    int print_every = Nt / 20;
    if (print_every < 1) print_every = 1;

    for (int n = 0; n <= Nt; n++) {
        double t = n * dt;

        if (n % track_every == 0) {
            double xL, xR, EL, ER;
            compute_energy_centroids(phi, vel, Nx, dx, xm, m2, mu, &xL, &xR, &EL, &ER);
            double sep = xR - xL;

            /* Total energy for conservation check */
            double Et = EL + ER;

            if (n_track < max_track) {
                sep_hist[n_track] = sep;
                t_track[n_track]  = t;
                n_track++;
            }

            if (fts)
                fprintf(fts, "%.6f\t%.6f\t%.6f\t%.6f\t%.6f\n", t, sep, xR, xL, Et);

            if (n % print_every == 0)
                printf("      t=%7.0f  sep=%.4f  xR=%.4f  xL=%.4f  E=%.4f\n",
                       t, sep, xR, xL, Et);
        }

        if (n == Nt) break;

        /* Velocity Verlet */
        for (int a = 0; a < 3; a++)
            for (int i = 1; i < Nx - 1; i++)
                vel[a][i] += 0.5 * dt * acc[a][i];
        for (int a = 0; a < 3; a++)
            for (int i = 1; i < Nx - 1; i++)
                phi[a][i] += dt * vel[a][i];
        COMPUTE_ACC_INT(Nx);
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

    /* --- Fit sep(t) using cycle-averaged separation --- */
    /* First: compute cycle-averaged separation by averaging over breathing period.
     * Breathing period T ~ 2*pi/omega ~ 7-8. Use a boxcar average over T_avg=10. */
    if (n_track > 100) {
        double T_avg = 10.0;  /* averaging window (> breathing period) */

        /* Compute averaged separation */
        int n_avg = 0;
        int max_avg = n_track;
        double *t_avg_arr = malloc(max_avg * sizeof(double));
        double *s_avg_arr = malloc(max_avg * sizeof(double));

        for (int j = 0; j < n_track; j++) {
            double tc = t_track[j];
            if (tc > t_run - T_avg) break;  /* don't average past end */

            /* Find average of sep_hist over [tc - T_avg/2, tc + T_avg/2] */
            double sum_s = 0;
            int cnt = 0;
            for (int k = j; k < n_track; k++) {
                if (t_track[k] > tc + T_avg/2) break;
                if (t_track[k] >= tc - T_avg/2) {
                    sum_s += sep_hist[k];
                    cnt++;
                }
            }
            /* Also look backwards */
            for (int k = j - 1; k >= 0; k--) {
                if (t_track[k] < tc - T_avg/2) break;
                sum_s += sep_hist[k];
                cnt++;
            }

            if (cnt > 0 && tc >= T_avg/2) {
                t_avg_arr[n_avg] = tc;
                s_avg_arr[n_avg] = sum_s / cnt;
                n_avg++;
            }
        }

        /* Fit averaged separation over t in [20, 500] to s = a + b*t + c*t^2 */
        double S[5] = {0};
        double Sy[3] = {0};
        int nfit = 0;

        for (int j = 0; j < n_avg; j++) {
            if (t_avg_arr[j] < 20.0) continue;   /* skip initial transient */
            if (t_avg_arr[j] > 500.0) break;
            double tj = t_avg_arr[j];
            double sj = s_avg_arr[j];
            double tk = 1.0;
            for (int k = 0; k < 5; k++) { S[k] += tk; tk *= tj; }
            Sy[0] += sj;
            Sy[1] += sj * tj;
            Sy[2] += sj * tj * tj;
            nfit++;
        }

        if (nfit >= 10) {
            /* 3x3 normal equations */
            double M[3][4] = {
                {S[0], S[1], S[2], Sy[0]},
                {S[1], S[2], S[3], Sy[1]},
                {S[2], S[3], S[4], Sy[2]}
            };

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

            /* sep(t) = coeff[0] + coeff[1]*t + coeff[2]*t^2
             * d^2(sep)/dt^2 = 2*coeff[2]
             * Force on each oscillon = (M/2) * d^2(sep)/dt^2 = M * coeff[2] */
            res.F_fit = coeff[2];
            res.valid = 1;

            printf("    Fit(avg): a=%.4f b=%.4e c=%.4e  => accel=%.4e\n",
                   coeff[0], coeff[1], coeff[2], coeff[2]);
        }

        /* Also compute simple average velocity over [200, 1000] as fallback */
        {
            int j_start = -1, j_end = -1;
            for (int j = 0; j < n_avg; j++) {
                if (j_start < 0 && t_avg_arr[j] >= 200.0) j_start = j;
                if (t_avg_arr[j] <= 1000.0) j_end = j;
            }
            if (j_start >= 0 && j_end > j_start) {
                double v_avg = (s_avg_arr[j_end] - s_avg_arr[j_start]) /
                               (t_avg_arr[j_end] - t_avg_arr[j_start]);
                printf("    v_avg[200,1000] = %.4e\n", v_avg);
            }
        }

        free(t_avg_arr);
        free(s_avg_arr);
    }

    res.D_final = (n_track > 0) ? sep_hist[n_track - 1] : D;

    for (int a = 0; a < 3; a++) { free(phi[a]); free(vel[a]); free(acc[a]); }
    free(damp); free(sep_hist); free(t_track);

    return res;
}

/* ===================================================================
 *  Yukawa fit: ln|F| = ln|F0| - D/lambda  =>  linear fit in D
 * =================================================================== */
static void fit_yukawa(double *D_vals, double *F_vals, int nf,
                       double *F0_out, double *lambda_out)
{
    double sx = 0, sy = 0, sxx = 0, sxy = 0;
    int nfit = 0;
    for (int j = 0; j < nf; j++) {
        if (fabs(F_vals[j]) < 1e-30) continue;
        double lnF = log(fabs(F_vals[j]));
        sx  += D_vals[j];
        sy  += lnF;
        sxx += D_vals[j] * D_vals[j];
        sxy += D_vals[j] * lnF;
        nfit++;
    }
    if (nfit < 2) { *F0_out = 0; *lambda_out = 0; return; }
    double det = nfit * sxx - sx * sx;
    if (fabs(det) < 1e-30) { *F0_out = 0; *lambda_out = 0; return; }
    double slope     = (nfit * sxy - sx * sy) / det;
    double intercept = (sy * sxx - sx * sxy) / det;
    *lambda_out = -1.0 / slope;
    *F0_out     = exp(intercept);
}

/* =================================================================== */
int main(int argc, char **argv)
{
    parse_args(argc, argv);

    printf("=== V23-F: Critical Gravity Phase 2 Redo ===\n\n");
    printf("Mode: %s  kappa=%.1f  mass=%.3f\n\n", mode, kappa, mass);

    if (!strcmp(mode, "equil")) {
        /* Single equilibration */
        printf("--- Equilibrating mu=%.0f ---\n", mu_single);
        profile_t prof = run_equilibrate(mu_single, 1);
        profile_free(&prof);
    }
    else if (!strcmp(mode, "interact")) {
        /* Single interaction: load profile and run */
        char path[600];
        snprintf(path, sizeof(path), "%s/profile_mu%d.dat", outdir, (int)fabs(mu_single));
        printf("--- Loading profile from %s ---\n", path);
        profile_t prof = {0};
        if (profile_load(path, &prof) != 0) {
            fprintf(stderr, "Failed to load profile!\n");
            return 1;
        }
        printf("  omega=%.4f  A_peak=%.4f  E_total=%.4f\n",
               prof.omega, prof.A_peak, prof.E_total);

        printf("\n--- Running interaction mu=%.0f D=%.0f ---\n", mu_single, D_single);
        interact_result_t ir = run_interact(mu_single, D_single, &prof);
        printf("  F_fit=%.4e  D_final=%.4f  valid=%d\n", ir.F_fit, ir.D_final, ir.valid);

        profile_free(&prof);
    }
    else if (!strcmp(mode, "full")) {
        /* Full protocol */

        /* Step 1: Equilibrate all mu values */
        printf("========== STEP 1: EQUILIBRATION ==========\n\n");
        profile_t profiles[10];
        double omegas[10], E_totals[10], A_peaks[10];
        int alive[10];
        memset(alive, 0, sizeof(alive));

        for (int mi = 0; mi < n_mu; mi++) {
            double mu_i = mu_list[mi];
            printf("\n--- Equilibrating mu=%.0f ---\n", mu_i);
            profiles[mi] = run_equilibrate(mu_i, 1);
            omegas[mi]   = profiles[mi].omega;
            E_totals[mi] = profiles[mi].E_total;
            A_peaks[mi]  = profiles[mi].A_peak;
            alive[mi]    = (A_peaks[mi] > 0.05 && omegas[mi] > 0.01 && omegas[mi] < mass);
            printf("  => omega=%.4f  E=%.4f  A=%.4f  alive=%s\n",
                   omegas[mi], E_totals[mi], A_peaks[mi], alive[mi] ? "YES" : "NO");
        }

        /* Summary table */
        printf("\n\n========== EQUILIBRATION SUMMARY ==========\n\n");
        printf("%-6s  %-8s  %-10s  %-10s  %-10s  %-5s\n",
               "mu", "omega", "E_total", "A_peak", "xi_pred", "alive");
        printf("------  --------  ----------  ----------  ----------  -----\n");
        for (int mi = 0; mi < n_mu; mi++) {
            double k2 = mass*mass - omegas[mi]*omegas[mi];
            double xi_pred = (k2 > 0) ? 1.0/sqrt(k2) : 9999.0;
            printf("%-6.0f  %-8.4f  %-10.4f  %-10.4f  %-10.4f  %-5s\n",
                   mu_list[mi], omegas[mi], E_totals[mi], A_peaks[mi],
                   xi_pred, alive[mi] ? "YES" : "NO");
        }

        /* Step 2 & 3: Interactions */
        printf("\n\n========== STEP 2: TWO-OSCILLON INTERACTIONS ==========\n\n");

        /* Results storage */
        double F_table[10][10];   /* F_table[mi][di] */
        double D_table[10][10];
        int    valid_table[10][10];
        memset(valid_table, 0, sizeof(valid_table));

        for (int mi = 0; mi < n_mu; mi++) {
            if (!alive[mi]) {
                printf("\n--- Skipping mu=%.0f (not alive) ---\n", mu_list[mi]);
                continue;
            }

            printf("\n--- mu=%.0f: interactions ---\n", mu_list[mi]);
            for (int di = 0; di < n_D; di++) {
                double D = D_list[di];
                printf("\n  D=%.0f:\n", D);
                interact_result_t ir = run_interact(mu_list[mi], D, &profiles[mi]);
                F_table[mi][di] = ir.F_fit;
                D_table[mi][di] = D;
                valid_table[mi][di] = ir.valid;
                printf("    F=%.4e  D_final=%.4f\n", ir.F_fit, ir.D_final);
            }
        }

        /* Step 3: Yukawa fits */
        printf("\n\n========== STEP 3: YUKAWA FITS ==========\n\n");
        printf("%-6s  %-8s  %-10s  %-10s  %-10s  %-10s  %-10s\n",
               "mu", "omega", "xi_pred", "lambda", "F0", "ratio", "sign");
        printf("------  --------  ----------  ----------  ----------  ----------  ----------\n");

        /* Master results file */
        char rpath[600];
        snprintf(rpath, sizeof(rpath), "%s/yukawa_fits.tsv", outdir);
        FILE *fr = fopen(rpath, "w");
        if (fr) fprintf(fr, "mu\tomega\txi_pred\tlambda\tF0\tE_total\n");

        for (int mi = 0; mi < n_mu; mi++) {
            if (!alive[mi]) continue;

            /* Collect valid F(D) data */
            double Dv[10], Fv[10];
            int nf = 0;
            int sign_pos = 0, sign_neg = 0;
            for (int di = 0; di < n_D; di++) {
                if (valid_table[mi][di]) {
                    Dv[nf] = D_table[mi][di];
                    Fv[nf] = F_table[mi][di];
                    if (Fv[nf] > 0) sign_pos++;
                    else sign_neg++;
                    nf++;
                }
            }

            double F0, lambda;
            fit_yukawa(Dv, Fv, nf, &F0, &lambda);

            double k2 = mass*mass - omegas[mi]*omegas[mi];
            double xi_pred = (k2 > 0) ? 1.0/sqrt(k2) : 9999.0;

            const char *sign_str = (sign_pos > sign_neg) ? "repuls" :
                                   (sign_neg > sign_pos) ? "attract" : "mixed";

            printf("%-6.0f  %-8.4f  %-10.4f  %-10.4f  %-10.4e  %-10.4f  %-10s\n",
                   mu_list[mi], omegas[mi], xi_pred, lambda, F0,
                   (xi_pred > 0.01 && xi_pred < 9000) ? lambda / xi_pred : 0.0,
                   sign_str);

            if (fr) fprintf(fr, "%.2f\t%.6f\t%.6f\t%.6f\t%.6e\t%.6f\n",
                            mu_list[mi], omegas[mi], xi_pred, lambda, F0, E_totals[mi]);

            /* Print individual F(D) values */
            printf("  F(D): ");
            for (int j = 0; j < nf; j++)
                printf("D=%.0f:%.2e  ", Dv[j], Fv[j]);
            printf("\n");
        }

        if (fr) fclose(fr);

        /* Cleanup */
        for (int mi = 0; mi < n_mu; mi++)
            profile_free(&profiles[mi]);
    }
    else {
        fprintf(stderr, "Unknown mode: %s (use equil, interact, or full)\n", mode);
        return 1;
    }

    printf("\n=== Done ===\n");
    return 0;
}
