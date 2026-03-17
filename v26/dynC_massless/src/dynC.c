/*
 * dynC.c — V26-DynC: Massless Propagating Braid
 *
 * Tests whether a propagating helical wave (m=0, omega=k) can be confined
 * by the triple product potential, preventing the collapse seen in static m=0.
 *
 * Two run modes:
 *   mode=0: Propagating helical wave (vel != 0)
 *   mode=1: Static control (vel = 0, same spatial init)
 *
 * Lagrangian (triple product only, m=0):
 *   L = Sum_a [1/2(dt phi_a)^2 - 1/2|grad phi_a|^2] - V(P)
 *   V(P) = (mu/2) P^2 / (1 + kappa P^2),  P = phi_1 phi_2 phi_3
 *
 * Compile: gcc -O3 -fopenmp -Wall -o dynC src/dynC.c -lm
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>

/* --- Parameters --- */
static double mu_pot    = -20.0;
static double kappa     = 20.0;
static double mass      = 0.0;     /* MASSLESS */
static double A0        = 0.8;
static double R_tube    = 3.0;

static int    N         = 128;
static double L         = 20.0;
static double tfinal    = 500.0;
static double cfl_frac  = 0.20;
static int    mode      = 0;       /* 0=propagating, 1=static control */
static char   outdir[512] = "data";

static void parse_args(int argc, char **argv)
{
    for (int i = 1; i < argc - 1; i += 2) {
        if      (!strcmp(argv[i], "-mu"))     mu_pot   = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-kappa"))  kappa    = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-mass"))   mass     = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-A"))      A0       = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-N"))      N        = atoi(argv[i+1]);
        else if (!strcmp(argv[i], "-L"))      L        = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-tfinal")) tfinal   = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-cfl"))    cfl_frac = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-Rtube"))  R_tube   = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-mode"))   mode     = atoi(argv[i+1]);
        else if (!strcmp(argv[i], "-o"))      strncpy(outdir, argv[i+1], sizeof(outdir)-1);
        else { fprintf(stderr, "Unknown flag: %s\n", argv[i]); exit(1); }
    }
}

/* --- Index helpers --- */
#define IDX(i,j,k) ((long)(i)*N*N + (long)(j)*N + (long)(k))

/* --- Globals --- */
static double *phi[3], *vel[3], *acc[3];
static double *damp;
static double dx, dx2, m2, dt;
static long Ngrid;

/* --- Mode names --- */
static const char *mode_names[] = {
    "propagating",
    "static_control"
};

/* --- Periodic wrap in z --- */
static inline int wrap_z(int k)
{
    if (k < 0)  return k + N;
    if (k >= N) return k - N;
    return k;
}

/* --- Initialization: Propagating helical wave --- */
static void init_helical(int propagating)
{
    double Lz = 2.0 * L;
    double k = 2.0 * M_PI / Lz;
    double omega = k;  /* m=0 -> omega = k (massless dispersion) */

    printf("  k = %.6f, omega = %.6f (Lz = %.1f)\n", k, omega, Lz);

    #pragma omp parallel for schedule(static)
    for (long idx = 0; idx < Ngrid; idx++) {
        int ii = idx / (N * N);
        int jj = (idx / N) % N;
        int kk = idx % N;
        double x = -L + ii * dx;
        double y = -L + jj * dx;
        double z = -L + kk * dx;

        double r_perp = sqrt(x * x + y * y);
        double envelope = A0 * exp(-r_perp * r_perp / (2.0 * R_tube * R_tube));

        for (int a = 0; a < 3; a++) {
            double phase = k * z + 2.0 * M_PI * a / 3.0;
            phi[a][idx] = envelope * cos(phase);
            if (propagating) {
                /* d/dt [A*cos(kz + 2pi*a/3 - omega*t)] at t=0
                 * = omega * A * sin(kz + 2pi*a/3) */
                vel[a][idx] = omega * envelope * sin(phase);
            } else {
                vel[a][idx] = 0.0;
            }
        }
    }
}

/* --- Compute acceleration (triple product only, periodic z) --- */
static void compute_acc(void)
{
    for (int a = 0; a < 3; a++) {
        #pragma omp parallel for schedule(static)
        for (int i = 2; i < N-2; i++) {
            for (int j = 2; j < N-2; j++) {
                for (int k = 0; k < N; k++) {
                    long idx = IDX(i, j, k);

                    /* 7-point Laplacian with periodic z */
                    int km1 = wrap_z(k-1);
                    int kp1 = wrap_z(k+1);

                    double lapl = (phi[a][IDX(i+1,j,k)] + phi[a][IDX(i-1,j,k)]
                                 + phi[a][IDX(i,j+1,k)] + phi[a][IDX(i,j-1,k)]
                                 + phi[a][IDX(i,j,kp1)] + phi[a][IDX(i,j,km1)]
                                 - 6.0 * phi[a][idx]) / dx2;

                    /* Triple-product potential force */
                    double p0 = phi[0][idx], p1 = phi[1][idx], p2 = phi[2][idx];
                    double P = p0 * p1 * p2;
                    double P2 = P * P;
                    double denom2 = (1.0 + kappa * P2);
                    denom2 *= denom2;

                    double dP;
                    switch (a) {
                        case 0: dP = p1 * p2; break;
                        case 1: dP = p0 * p2; break;
                        default: dP = p0 * p1; break;
                    }
                    double dVdphi = mu_pot * P * dP / denom2;

                    acc[a][idx] = lapl - m2 * phi[a][idx] - dVdphi;
                }
            }
        }

        /* Boundary: zero acceleration for x,y boundaries */
        #pragma omp parallel for schedule(static)
        for (long idx = 0; idx < Ngrid; idx++) {
            int ii = idx / (N * N);
            int jj = (idx / N) % N;
            if (ii < 2 || ii >= N-2 || jj < 2 || jj >= N-2) {
                acc[a][idx] = 0.0;
            }
        }
    }
}

/* --- Velocity Verlet step --- */
static void verlet_step(void)
{
    /* 1. Half-kick */
    for (int a = 0; a < 3; a++) {
        #pragma omp parallel for schedule(static)
        for (long idx = 0; idx < Ngrid; idx++)
            vel[a][idx] += 0.5 * dt * acc[a][idx];
    }

    /* 2. Drift */
    for (int a = 0; a < 3; a++) {
        #pragma omp parallel for schedule(static)
        for (long idx = 0; idx < Ngrid; idx++)
            phi[a][idx] += dt * vel[a][idx];
    }

    /* 3. Recompute acceleration */
    compute_acc();

    /* 4. Half-kick */
    for (int a = 0; a < 3; a++) {
        #pragma omp parallel for schedule(static)
        for (long idx = 0; idx < Ngrid; idx++)
            vel[a][idx] += 0.5 * dt * acc[a][idx];
    }

    /* 5. Absorbing boundary damping (transverse only) */
    for (int a = 0; a < 3; a++) {
        #pragma omp parallel for schedule(static)
        for (long idx = 0; idx < Ngrid; idx++) {
            vel[a][idx] *= damp[idx];
            phi[a][idx] *= damp[idx];
        }
    }
}

/* --- Diagnostics --- */
typedef struct {
    double Ek, Eg, Em, Ep, Et;
    double peak[3];
    double fc;            /* core energy fraction */
    double peak_P;        /* max |triple product| */
    double rho_center;    /* energy density at center */
    double Pz;            /* z-momentum (for propagation check) */
} Diag;

static Diag compute_diag(double core_radius)
{
    Diag d;
    memset(&d, 0, sizeof(d));

    /* Center energy density */
    {
        long ic = IDX(N/2, N/2, N/2);
        int ci = N/2, cj = N/2, ck = N/2;
        double rho = 0;
        for (int a = 0; a < 3; a++) {
            rho += 0.5 * vel[a][ic] * vel[a][ic];
            double gx = (phi[a][IDX(ci+1,cj,ck)] - phi[a][IDX(ci-1,cj,ck)]) / (2*dx);
            double gy = (phi[a][IDX(ci,cj+1,ck)] - phi[a][IDX(ci,cj-1,ck)]) / (2*dx);
            double gz = (phi[a][IDX(ci,cj,ck+1)] - phi[a][IDX(ci,cj,ck-1)]) / (2*dx);
            rho += 0.5 * (gx*gx + gy*gy + gz*gz);
            rho += 0.5 * m2 * phi[a][ic] * phi[a][ic];
        }
        double p0c = phi[0][ic], p1c = phi[1][ic], p2c = phi[2][ic];
        double Pc = p0c * p1c * p2c;
        double Pc2 = Pc * Pc;
        rho += 0.5 * mu_pot * Pc2 / (1.0 + kappa * Pc2);
        d.rho_center = rho;
    }

    double Ek=0, Eg=0, Em=0, Ep=0;
    double Ecore=0, Eall=0;
    double peak[3] = {0,0,0};
    double peak_P = 0;
    double Pz = 0;

    #pragma omp parallel
    {
        double lEk=0, lEg=0, lEm=0, lEp=0, lEc=0, lEa=0;
        double lpk[3] = {0,0,0};
        double lpkP = 0;
        double lPz = 0;

        #pragma omp for schedule(static) nowait
        for (int i = 2; i < N-2; i++) {
            double x = -L + i * dx;
            for (int j = 2; j < N-2; j++) {
                double y = -L + j * dx;
                for (int k = 0; k < N; k++) {
                    long idx = IDX(i,j,k);
                    double dV = dx * dx * dx;
                    double e_loc = 0;

                    int km1 = wrap_z(k-1);
                    int kp1 = wrap_z(k+1);

                    for (int a = 0; a < 3; a++) {
                        double v2 = vel[a][idx] * vel[a][idx];
                        lEk += 0.5 * v2 * dV;
                        e_loc += 0.5 * v2;

                        double gx = (phi[a][IDX(i+1,j,k)] - phi[a][IDX(i-1,j,k)]) / (2*dx);
                        double gy = (phi[a][IDX(i,j+1,k)] - phi[a][IDX(i,j-1,k)]) / (2*dx);
                        double gz = (phi[a][IDX(i,j,kp1)] - phi[a][IDX(i,j,km1)]) / (2*dx);
                        double grad2 = gx*gx + gy*gy + gz*gz;
                        lEg += 0.5 * grad2 * dV;
                        e_loc += 0.5 * grad2;

                        double mass_e = 0.5 * m2 * phi[a][idx] * phi[a][idx];
                        lEm += mass_e * dV;
                        e_loc += mass_e;

                        double ap = fabs(phi[a][idx]);
                        if (ap > lpk[a]) lpk[a] = ap;

                        /* z-momentum: -dot(phi_a) * d_z(phi_a) */
                        lPz -= vel[a][idx] * gz * dV;
                    }

                    double p0 = phi[0][idx], p1 = phi[1][idx], p2 = phi[2][idx];
                    double P = p0 * p1 * p2;
                    double P2 = P * P;
                    double Vloc = 0.5 * mu_pot * P2 / (1.0 + kappa * P2);
                    lEp += Vloc * dV;
                    e_loc += Vloc;

                    double absP = fabs(P);
                    if (absP > lpkP) lpkP = absP;

                    lEa += e_loc * dV;
                    double r_perp = sqrt(x*x + y*y);
                    if (r_perp < core_radius) lEc += e_loc * dV;
                }
            }
        }

        #pragma omp critical
        {
            Ek += lEk; Eg += lEg; Em += lEm; Ep += lEp;
            Ecore += lEc; Eall += lEa;
            Pz += lPz;
            for (int a = 0; a < 3; a++)
                if (lpk[a] > peak[a]) peak[a] = lpk[a];
            if (lpkP > peak_P) peak_P = lpkP;
        }
    }

    d.Ek = Ek; d.Eg = Eg; d.Em = Em; d.Ep = Ep;
    d.Et = Ek + Eg + Em + Ep;
    for (int a = 0; a < 3; a++) d.peak[a] = peak[a];
    d.fc = (Eall > 1e-20) ? Ecore / Eall : 0.0;
    d.peak_P = peak_P;
    d.Pz = Pz;

    return d;
}

/* --- Fast z-averaged radial profile --- */
static void write_radial_profile(const char *fname)
{
    int nr = 60;
    double dr = L / nr;
    double *amp_sum = calloc(nr, sizeof(double));
    int *amp_cnt = calloc(nr, sizeof(int));

    for (int i = 2; i < N-2; i++) {
        double x = -L + i * dx;
        for (int j = 2; j < N-2; j++) {
            double y = -L + j * dx;
            double r = sqrt(x*x + y*y);
            int ir = (int)(r / dr);
            if (ir >= nr) continue;
            /* Average amplitude over z */
            for (int k = 0; k < N; k++) {
                long idx = IDX(i,j,k);
                double a2 = 0;
                for (int a = 0; a < 3; a++)
                    a2 += phi[a][idx] * phi[a][idx];
                amp_sum[ir] += sqrt(a2);
                amp_cnt[ir]++;
            }
        }
    }

    FILE *fp = fopen(fname, "w");
    if (fp) {
        fprintf(fp, "r_perp\tamp_avg\n");
        for (int ir = 0; ir < nr; ir++) {
            double r = (ir + 0.5) * dr;
            double avg = (amp_cnt[ir] > 0) ? amp_sum[ir] / amp_cnt[ir] : 0.0;
            fprintf(fp, "%.4f\t%.6e\n", r, avg);
        }
        fclose(fp);
    }
    free(amp_sum);
    free(amp_cnt);
}

/* --- Run one configuration --- */
static void run_config(int run_mode)
{
    const char *name = mode_names[run_mode];
    printf("\n================================================================\n");
    printf("  V26-DynC: %s (mode=%d)\n", name, run_mode);
    printf("================================================================\n\n");

    m2 = mass * mass;

    /* Zero fields */
    for (int a = 0; a < 3; a++) {
        memset(phi[a], 0, Ngrid * sizeof(double));
        memset(vel[a], 0, Ngrid * sizeof(double));
        memset(acc[a], 0, Ngrid * sizeof(double));
    }

    /* Initialize */
    int propagating = (run_mode == 0);
    init_helical(propagating);

    /* Damping: transverse absorbing layer only (periodic z) */
    double R_abs_inner = L * 0.70;
    double R_abs_outer = L * 0.95;
    #pragma omp parallel for schedule(static)
    for (long idx = 0; idx < Ngrid; idx++) {
        int ii = idx / (N * N);
        int jj = (idx / N) % N;
        double x = -L + ii * dx;
        double y = -L + jj * dx;
        double rperp = sqrt(x*x + y*y);
        if (rperp > R_abs_inner) {
            double f = (rperp - R_abs_inner) / (R_abs_outer - R_abs_inner);
            if (f > 1.0) f = 1.0;
            damp[idx] = 1.0 - 0.98 * f * f;
        } else {
            damp[idx] = 1.0;
        }
    }

    double core_radius = 8.0;
    int Nt = (int)(tfinal / dt) + 1;

    /* DFT storage */
    int max_dft = 50000;
    double *rho_hist  = malloc(max_dft * sizeof(double));
    double *phi0_hist = malloc(max_dft * sizeof(double));
    double *t_hist    = malloc(max_dft * sizeof(double));
    int n_dft = 0;

    int diag_every  = Nt / 5000;  if (diag_every < 1) diag_every = 1;
    int print_every = Nt / 30;    if (print_every < 1) print_every = 1;
    int dft_every   = Nt / max_dft; if (dft_every < 1) dft_every = 1;

    /* Output file */
    char path[600];
    snprintf(path, sizeof(path), "%s/dynC_%s.tsv", outdir, name);
    FILE *f1 = fopen(path, "w");
    if (!f1) { fprintf(stderr, "Cannot open %s\n", path); return; }
    fprintf(f1, "time\tE_total\tE_kin\tE_grad\tE_mass\tE_pot\t"
                "fc\tpeak0\tpeak1\tpeak2\tpeak_P\trho_center\tPz\n");

    compute_acc();

    double wall_start = omp_get_wtime();

    /* Initial diagnostics */
    Diag d0 = compute_diag(core_radius);
    double E0 = d0.Et;
    printf("  t=%7.1f  E=%.4f  fc=%.4f  pk=(%.4f,%.4f,%.4f)  |P|=%.6f  Pz=%.4f\n",
           0.0, d0.Et, d0.fc, d0.peak[0], d0.peak[1], d0.peak[2], d0.peak_P, d0.Pz);

    /* Write initial radial profile */
    snprintf(path, sizeof(path), "%s/dynC_%s_profile_t0.tsv", outdir, name);
    write_radial_profile(path);

    for (int n = 0; n <= Nt; n++) {
        double t = n * dt;

        /* DFT history */
        if (n % dft_every == 0 && n_dft < max_dft) {
            long ic = IDX(N/2, N/2, N/2);
            double rho = 0;
            int ci = N/2, cj = N/2, ck = N/2;
            for (int a = 0; a < 3; a++) {
                rho += 0.5 * vel[a][ic] * vel[a][ic];
                double gx = (phi[a][IDX(ci+1,cj,ck)] - phi[a][IDX(ci-1,cj,ck)]) / (2*dx);
                double gy = (phi[a][IDX(ci,cj+1,ck)] - phi[a][IDX(ci,cj-1,ck)]) / (2*dx);
                double gz = (phi[a][IDX(ci,cj,ck+1)] - phi[a][IDX(ci,cj,ck-1)]) / (2*dx);
                rho += 0.5 * (gx*gx + gy*gy + gz*gz);
            }
            rho_hist[n_dft] = rho;
            phi0_hist[n_dft] = phi[0][ic];
            t_hist[n_dft] = t;
            n_dft++;
        }

        int do_diag  = (n % diag_every == 0);
        int do_print = (n % print_every == 0);

        if (do_diag || do_print) {
            Diag dg = compute_diag(core_radius);

            if (do_diag)
                fprintf(f1, "%.4f\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t"
                            "%.6f\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\n",
                        t, dg.Et, dg.Ek, dg.Eg, dg.Em, dg.Ep,
                        dg.fc, dg.peak[0], dg.peak[1], dg.peak[2],
                        dg.peak_P, dg.rho_center, dg.Pz);

            if (do_print) {
                double elapsed = omp_get_wtime() - wall_start;
                double frac = (double)n / Nt;
                double eta_t = (frac > 0.001) ? elapsed * (1.0-frac)/frac : 0;
                printf("  t=%7.1f  E=%.4f (%.1f%%)  fc=%.4f  |P|=%.6f  Pz=%.4f  [%.0fs, ETA %.0fs]\n",
                       t, dg.Et, 100.0*dg.Et/E0, dg.fc, dg.peak_P, dg.Pz,
                       elapsed, eta_t);
                fflush(stdout);
            }
        }

        if (n == Nt) break;
        verlet_step();
    }

    fclose(f1);

    /* Final diagnostics */
    Diag dfinal = compute_diag(core_radius);
    double elapsed1 = omp_get_wtime() - wall_start;

    printf("\nEvolution complete (%.1f sec)\n", elapsed1);
    printf("  Final: E=%.4f (%.1f%% of initial)  fc=%.4f  |P|=%.6f  Pz=%.4f\n",
           dfinal.Et, 100.0*dfinal.Et/E0, dfinal.fc, dfinal.peak_P, dfinal.Pz);

    int survived = (dfinal.fc > 0.01 && dfinal.peak_P > 1e-6);
    printf("  Survived? %s (fc=%.4f, |P|max=%.6f)\n",
           survived ? "YES" : "NO", dfinal.fc, dfinal.peak_P);

    /* Write final radial profile */
    snprintf(path, sizeof(path), "%s/dynC_%s_profile_tfinal.tsv", outdir, name);
    write_radial_profile(path);

    /* DFT analysis */
    printf("\n--- DFT Analysis ---\n");
    if (n_dft >= 100) {
        int start = n_dft / 2;
        double mean_rho = 0;
        for (int j = start; j < n_dft; j++) mean_rho += rho_hist[j];
        mean_rho /= (n_dft - start);
        double var_rho = 0;
        for (int j = start; j < n_dft; j++) {
            double dr = rho_hist[j] - mean_rho;
            var_rho += dr * dr;
        }
        var_rho /= (n_dft - start);
        printf("  rho(center): mean=%.4e, variance=%.4e, rel_var=%.4e\n",
               mean_rho, var_rho,
               (fabs(mean_rho) > 1e-30) ? var_rho / (mean_rho * mean_rho) : 0.0);

        /* Find peak omega */
        int nf = 500;
        double omega_max = 5.0;
        double peak_pow = 0, peak_om = 0;
        for (int kk = 1; kk < nf; kk++) {
            double omega = omega_max * kk / nf;
            double re = 0, im = 0;
            for (int j = start; j < n_dft; j++) {
                double dtj = (j > start) ?
                    (t_hist[j]-t_hist[j-1]) : (t_hist[start+1]-t_hist[start]);
                double val = rho_hist[j] - mean_rho;
                re += val * cos(omega * t_hist[j]) * dtj;
                im += val * sin(omega * t_hist[j]) * dtj;
            }
            double pw = re*re + im*im;
            if (pw > peak_pow) { peak_pow = pw; peak_om = omega; }
        }
        printf("  Peak omega(rho_center) = %.4f (power=%.4e)\n", peak_om, peak_pow);

        /* Write DFT */
        snprintf(path, sizeof(path), "%s/dynC_%s_dft.tsv", outdir, name);
        FILE *fdft = fopen(path, "w");
        if (fdft) {
            fprintf(fdft, "omega\tpower_rho\tpower_phi\n");
            double mean_phi0 = 0;
            for (int j = start; j < n_dft; j++) mean_phi0 += phi0_hist[j];
            mean_phi0 /= (n_dft - start);
            for (int kk = 1; kk < nf; kk++) {
                double omega = omega_max * kk / nf;
                double re_r=0, im_r=0, re_p=0, im_p=0;
                for (int j = start; j < n_dft; j++) {
                    double dtj = (j > start) ?
                        (t_hist[j]-t_hist[j-1]) : (t_hist[start+1]-t_hist[start]);
                    double c = cos(omega * t_hist[j]);
                    double s = sin(omega * t_hist[j]);
                    re_r += (rho_hist[j] - mean_rho) * c * dtj;
                    im_r += (rho_hist[j] - mean_rho) * s * dtj;
                    re_p += (phi0_hist[j] - mean_phi0) * c * dtj;
                    im_p += (phi0_hist[j] - mean_phi0) * s * dtj;
                }
                fprintf(fdft, "%.6f\t%.6e\t%.6e\n", omega,
                        re_r*re_r+im_r*im_r, re_p*re_p+im_p*im_p);
            }
            fclose(fdft);
        }

        /* Write time history */
        snprintf(path, sizeof(path), "%s/dynC_%s_history.tsv", outdir, name);
        FILE *fh = fopen(path, "w");
        if (fh) {
            fprintf(fh, "time\trho_center\tphi0_center\n");
            for (int j = 0; j < n_dft; j++)
                fprintf(fh, "%.4f\t%.6e\t%.6e\n", t_hist[j], rho_hist[j], phi0_hist[j]);
            fclose(fh);
        }
    } else {
        printf("  Not enough DFT samples (%d)\n", n_dft);
    }

    free(rho_hist);
    free(phi0_hist);
    free(t_hist);
}

/* --- Main --- */
int main(int argc, char **argv)
{
    parse_args(argc, argv);

    dx  = 2.0 * L / (N - 1);
    dx2 = dx * dx;
    m2  = mass * mass;
    dt  = cfl_frac * dx;
    Ngrid = (long)N * N * N;

    printf("=== V26-DynC: Massless Propagating Braid ===\n");
    printf("Parameters:\n");
    printf("  mu=%.1f  kappa=%.1f  mass=%.3f  A0=%.3f  R_tube=%.1f\n",
           mu_pot, kappa, mass, A0, R_tube);
    printf("  N=%d  L=%.1f  dx=%.4f  dt=%.5f\n", N, L, dx, dt);
    printf("  Ngrid=%ld (%.1f M)  Memory: %.1f MB\n",
           Ngrid, Ngrid/1e6, Ngrid*8.0*10/1e6);
    printf("  Threads: %d\n", omp_get_max_threads());
    printf("  mode=%d (%s)  tfinal=%.0f\n", mode,
           (mode >= 0 && mode <= 1) ? mode_names[mode] : "unknown", tfinal);
    fflush(stdout);

    /* Allocate */
    for (int a = 0; a < 3; a++) {
        phi[a] = calloc(Ngrid, sizeof(double));
        vel[a] = calloc(Ngrid, sizeof(double));
        acc[a] = calloc(Ngrid, sizeof(double));
        if (!phi[a] || !vel[a] || !acc[a]) {
            fprintf(stderr, "Allocation failed\n"); return 1;
        }
    }
    damp = malloc(Ngrid * sizeof(double));
    if (!damp) { fprintf(stderr, "Allocation failed for damp\n"); return 1; }

    if (mode == 0 || mode == 1) {
        run_config(mode);
    } else if (mode == 99) {
        run_config(0);
        run_config(1);
    } else {
        fprintf(stderr, "Unknown mode %d. Use 0 (propagating), 1 (static), or 99 (both).\n", mode);
        return 1;
    }

    /* Cleanup */
    for (int a = 0; a < 3; a++) {
        free(phi[a]); free(vel[a]); free(acc[a]);
    }
    free(damp);

    printf("\n=== V26-DynC Complete ===\n");
    return 0;
}
