/*
 * m3.c — V27-M3: Electromagnetic Structure from Torsion
 *
 * Uses the m=0, mu=-50, kappa=50 propagating braid (V27-M4 optimal config).
 * ALL 3D, N=128, L=20, periodic z, absorbing x,y.
 *
 * Lagrangian:
 *   L = Sum_a [1/2(dt phi_a)^2 - 1/2(di phi_a)^2] - (mu/2)P^2/(1+kappa*P^2)
 *   m=0 (massless — confinement from triple-product alone)
 *
 * Tests:
 *   M3a: Map the torsion field (vorticity Omega_k = eps_{ijk} omega_{ij})
 *   M3b: Torsion flux Phi_T = integral Omega_z dx dy at various z-slices
 *   M3c: Two braids — torsion-mediated force (same vs opposite twist)
 *
 * Compile: gcc -O3 -fopenmp -Wall -o m3 src/m3.c -lm
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>

/* ─── Parameters ─── */
static double mu_pot    = -50.0;
static double kappa     = 50.0;
static double mass      = 0.0;
static double A0        = 0.8;
static double R_tube    = 3.0;

static int    N         = 128;
static double L         = 20.0;
static double tfinal    = 200.0;
static double cfl_frac  = 0.20;
static char   outdir[512] = "data";

/* ─── Index helpers ─── */
#define IDX(i,j,k) ((long)(i)*N*N + (long)(j)*N + (long)(k))

/* ─── Globals ─── */
static double *phi[3], *vel[3], *acc[3];
static double *damp;
static double dx, dx2, m2, dt;
static long Ngrid;

/* ─── Periodic z wrap ─── */
static inline int wrap_z(int k)
{
    if (k < 0)   return k + N;
    if (k >= N)  return k - N;
    return k;
}

/* ─── Initialization: Propagating Helical Braid ─── */
static void init_braid(double x_off, double y_off, double twist_sign)
{
    double Lz = 2.0 * L;
    double k_wave = twist_sign * 2.0 * M_PI / Lz;
    double omega = sqrt(k_wave * k_wave + m2);

    #pragma omp parallel for schedule(static)
    for (long idx = 0; idx < Ngrid; idx++) {
        int ii = idx / (N * N);
        int jj = (idx / N) % N;
        int kk = idx % N;
        double x = -L + ii * dx - x_off;
        double y = -L + jj * dx - y_off;
        double z = -L + kk * dx;

        double r_perp = sqrt(x * x + y * y);
        double envelope = A0 * exp(-r_perp * r_perp / (2.0 * R_tube * R_tube));

        for (int a = 0; a < 3; a++) {
            double phase = k_wave * z + 2.0 * M_PI * a / 3.0;
            phi[a][idx] += envelope * cos(phase);
            vel[a][idx] += omega * envelope * sin(phase);
        }
    }
}

/* ─── Compute acceleration ─── */
static void compute_acc(void)
{
    for (int a = 0; a < 3; a++) {
        #pragma omp parallel for schedule(static)
        for (int i = 2; i < N-2; i++) {
            for (int j = 2; j < N-2; j++) {
                for (int k = 0; k < N; k++) {
                    long idx = IDX(i, j, k);

                    int km1 = wrap_z(k - 1);
                    int kp1 = wrap_z(k + 1);

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

        /* Boundary: zero acceleration for x,y edges */
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

/* ─── Velocity Verlet step ─── */
static void verlet_step(void)
{
    for (int a = 0; a < 3; a++) {
        #pragma omp parallel for schedule(static)
        for (long idx = 0; idx < Ngrid; idx++)
            vel[a][idx] += 0.5 * dt * acc[a][idx];
    }

    for (int a = 0; a < 3; a++) {
        #pragma omp parallel for schedule(static)
        for (long idx = 0; idx < Ngrid; idx++)
            phi[a][idx] += dt * vel[a][idx];
    }

    compute_acc();

    for (int a = 0; a < 3; a++) {
        #pragma omp parallel for schedule(static)
        for (long idx = 0; idx < Ngrid; idx++)
            vel[a][idx] += 0.5 * dt * acc[a][idx];
    }

    /* Absorbing boundary damping (x,y only) */
    for (int a = 0; a < 3; a++) {
        #pragma omp parallel for schedule(static)
        for (long idx = 0; idx < Ngrid; idx++) {
            vel[a][idx] *= damp[idx];
            phi[a][idx] *= damp[idx];
        }
    }
}

/* ─── Setup damping layer ─── */
static void setup_damping(void)
{
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
}

/* ─── Diagnostics ─── */
typedef struct {
    double Ek, Eg, Em, Ep, Et;
    double peak_P;
    double Pz;
} Diag;

static Diag compute_diag(void)
{
    Diag d;
    memset(&d, 0, sizeof(d));

    double Ek = 0, Eg = 0, Em = 0, Ep = 0;
    double peak_P = 0;
    double Pz_total = 0;

    #pragma omp parallel
    {
        double lEk = 0, lEg = 0, lEm = 0, lEp = 0;
        double lpkP = 0, lPz = 0;

        #pragma omp for schedule(static) nowait
        for (int i = 2; i < N-2; i++) {
            for (int j = 2; j < N-2; j++) {
                for (int k = 0; k < N; k++) {
                    long idx = IDX(i, j, k);
                    double dV = dx * dx * dx;

                    int km1 = wrap_z(k - 1);
                    int kp1 = wrap_z(k + 1);

                    for (int a = 0; a < 3; a++) {
                        lEk += 0.5 * vel[a][idx] * vel[a][idx] * dV;

                        double gx = (phi[a][IDX(i+1,j,k)] - phi[a][IDX(i-1,j,k)]) / (2*dx);
                        double gy = (phi[a][IDX(i,j+1,k)] - phi[a][IDX(i,j-1,k)]) / (2*dx);
                        double gz = (phi[a][IDX(i,j,kp1)] - phi[a][IDX(i,j,km1)]) / (2*dx);
                        lEg += 0.5 * (gx*gx + gy*gy + gz*gz) * dV;
                        lEm += 0.5 * m2 * phi[a][idx] * phi[a][idx] * dV;
                        lPz += -vel[a][idx] * gz * dV;
                    }

                    double p0 = phi[0][idx], p1 = phi[1][idx], p2 = phi[2][idx];
                    double P_val = p0 * p1 * p2;
                    double P2 = P_val * P_val;
                    lEp += 0.5 * mu_pot * P2 / (1.0 + kappa * P2) * dV;

                    double absP = fabs(P_val);
                    if (absP > lpkP) lpkP = absP;
                }
            }
        }

        #pragma omp critical
        {
            Ek += lEk; Eg += lEg; Em += lEm; Ep += lEp;
            Pz_total += lPz;
            if (lpkP > peak_P) peak_P = lpkP;
        }
    }

    d.Ek = Ek; d.Eg = Eg; d.Em = Em; d.Ep = Ep;
    d.Et = Ek + Eg + Em + Ep;
    d.peak_P = peak_P;
    d.Pz = Pz_total;
    return d;
}

/* ─── Evolve braid to settle ─── */
static void evolve(int Nt, int print_interval)
{
    compute_acc();
    double wall_start = omp_get_wtime();

    for (int n = 1; n <= Nt; n++) {
        verlet_step();

        if (print_interval > 0 && (n % print_interval == 0 || n == Nt)) {
            double t = n * dt;
            Diag d = compute_diag();
            double elapsed = omp_get_wtime() - wall_start;
            double frac = (double)n / Nt;
            double eta = (frac > 0.001) ? elapsed * (1.0 - frac) / frac : 0;
            printf("  t=%7.1f  E=%.2f  |P|=%.4f  Pz=%.4f  [%.0fs, ETA %.0fs]\n",
                   t, d.Et, d.peak_P, d.Pz, elapsed, eta);
            fflush(stdout);
        }
    }
}

/* ─── Compute torsion (vorticity) at a single point ─── */
static void compute_vorticity(int i, int j, int k,
                              double *Ox, double *Oy, double *Oz)
{
    /*
     * omega_{ij} = (1/2)(d_i phi_j - d_j phi_i)
     *   Here phi_a is the a-th field, and i,j are spatial indices.
     *   The "displacement vector" at each point is (phi_0, phi_1, phi_2).
     *   So d_i phi_j means the spatial-i derivative of the j-th field component.
     *
     * omega_{xy} = (1/2)(dx phi_y - dy phi_x) = (1/2)(dx phi_1 - dy phi_0)
     * omega_{yz} = (1/2)(dy phi_z - dz phi_y) = (1/2)(dy phi_2 - dz phi_1)
     * omega_{zx} = (1/2)(dz phi_x - dx phi_z) = (1/2)(dz phi_0 - dx phi_2)
     *
     * Vorticity (dual):
     * Omega_x = 2*omega_{yz} = dy phi_2 - dz phi_1
     * Omega_y = 2*omega_{zx} = dz phi_0 - dx phi_2
     * Omega_z = 2*omega_{xy} = dx phi_1 - dy phi_0
     */

    int km1 = wrap_z(k - 1);
    int kp1 = wrap_z(k + 1);

    /* Central differences for spatial derivatives */
    /* d_x phi_a = (phi_a[i+1,j,k] - phi_a[i-1,j,k]) / (2*dx) */
    double dx_phi[3], dy_phi[3], dz_phi[3];
    for (int a = 0; a < 3; a++) {
        dx_phi[a] = (phi[a][IDX(i+1,j,k)] - phi[a][IDX(i-1,j,k)]) / (2.0 * dx);
        dy_phi[a] = (phi[a][IDX(i,j+1,k)] - phi[a][IDX(i,j-1,k)]) / (2.0 * dx);
        dz_phi[a] = (phi[a][IDX(i,j,kp1)] - phi[a][IDX(i,j,km1)]) / (2.0 * dx);
    }

    /* Omega_x = dy phi_2 - dz phi_1 */
    *Ox = dy_phi[2] - dz_phi[1];
    /* Omega_y = dz phi_0 - dx phi_2 */
    *Oy = dz_phi[0] - dx_phi[2];
    /* Omega_z = dx phi_1 - dy phi_0 */
    *Oz = dx_phi[1] - dy_phi[0];
}

/* ================================================================== */
/*  M3a: Map the torsion field                                         */
/* ================================================================== */
static void run_M3a(void)
{
    printf("\n================================================================\n");
    printf("  M3a: Map the Torsion Field\n");
    printf("================================================================\n\n");

    /* Output vorticity on z=0 plane (k=N/2) */
    {
        char path[600];
        snprintf(path, sizeof(path), "%s/m3a_vorticity_z0.tsv", outdir);
        FILE *f = fopen(path, "w");
        if (!f) { fprintf(stderr, "Cannot open %s\n", path); return; }
        fprintf(f, "# Vorticity on z=0 plane (k=%d)\n", N/2);
        fprintf(f, "x\ty\tOmega_x\tOmega_y\tOmega_z\t|Omega|\n");

        int kk = N / 2;
        double Oz_max = 0, Omega_max = 0;
        double rmax_sig = 0;

        for (int i = 3; i < N-3; i++) {
            double x = -L + i * dx;
            for (int j = 3; j < N-3; j++) {
                double y = -L + j * dx;
                double Ox, Oy, Oz;
                compute_vorticity(i, j, kk, &Ox, &Oy, &Oz);
                double mag = sqrt(Ox*Ox + Oy*Oy + Oz*Oz);
                fprintf(f, "%.4f\t%.4f\t%.6e\t%.6e\t%.6e\t%.6e\n",
                        x, y, Ox, Oy, Oz, mag);

                if (fabs(Oz) > Oz_max) Oz_max = fabs(Oz);
                if (mag > Omega_max) Omega_max = mag;
                double rperp = sqrt(x*x + y*y);
                if (mag > 0.01 * Omega_max && rperp > rmax_sig)
                    rmax_sig = rperp;
            }
        }
        fclose(f);
        printf("  z=0 plane: |Omega|_max = %.6e, |Omega_z|_max = %.6e\n",
               Omega_max, Oz_max);
        printf("  Output: %s\n", path);
    }

    /* Output vorticity on x=0 plane (i=N/2) */
    {
        char path[600];
        snprintf(path, sizeof(path), "%s/m3a_vorticity_x0.tsv", outdir);
        FILE *f = fopen(path, "w");
        if (!f) { fprintf(stderr, "Cannot open %s\n", path); return; }
        fprintf(f, "# Vorticity on x=0 plane (i=%d)\n", N/2);
        fprintf(f, "y\tz\tOmega_x\tOmega_y\tOmega_z\t|Omega|\n");

        int ii = N / 2;
        double Omega_max = 0;

        for (int j = 3; j < N-3; j++) {
            double y = -L + j * dx;
            for (int k = 0; k < N; k++) {
                double z = -L + k * dx;
                double Ox, Oy, Oz;
                compute_vorticity(ii, j, k, &Ox, &Oy, &Oz);
                double mag = sqrt(Ox*Ox + Oy*Oy + Oz*Oz);
                fprintf(f, "%.4f\t%.4f\t%.6e\t%.6e\t%.6e\t%.6e\n",
                        y, z, Ox, Oy, Oz, mag);
                if (mag > Omega_max) Omega_max = mag;
            }
        }
        fclose(f);
        printf("  x=0 plane: |Omega|_max = %.6e\n", Omega_max);
        printf("  Output: %s\n", path);
    }

    /* Radial profile of |Omega| on z=0 plane */
    {
        char path[600];
        snprintf(path, sizeof(path), "%s/m3a_vorticity_radial.tsv", outdir);
        FILE *f = fopen(path, "w");
        if (!f) { fprintf(stderr, "Cannot open %s\n", path); return; }
        fprintf(f, "# Radial profile of vorticity magnitude (z=0 plane, azimuthal average)\n");
        fprintf(f, "r\t|Omega|_avg\tOmega_z_avg\tn_samples\n");

        int kk = N / 2;
        int n_bins = 60;
        double dr = L / n_bins;
        double *sum_mag = calloc(n_bins, sizeof(double));
        double *sum_Oz  = calloc(n_bins, sizeof(double));
        int    *count   = calloc(n_bins, sizeof(int));

        for (int i = 3; i < N-3; i++) {
            double x = -L + i * dx;
            for (int j = 3; j < N-3; j++) {
                double y = -L + j * dx;
                double rperp = sqrt(x*x + y*y);
                int bin = (int)(rperp / dr);
                if (bin >= 0 && bin < n_bins) {
                    double Ox, Oy, Oz;
                    compute_vorticity(i, j, kk, &Ox, &Oy, &Oz);
                    sum_mag[bin] += sqrt(Ox*Ox + Oy*Oy + Oz*Oz);
                    sum_Oz[bin]  += Oz;
                    count[bin]++;
                }
            }
        }

        for (int b = 0; b < n_bins; b++) {
            if (count[b] > 0) {
                double r = (b + 0.5) * dr;
                fprintf(f, "%.4f\t%.6e\t%.6e\t%d\n",
                        r, sum_mag[b] / count[b], sum_Oz[b] / count[b], count[b]);
            }
        }
        fclose(f);
        free(sum_mag); free(sum_Oz); free(count);
        printf("  Output: %s\n", path);
    }

    /* Localization check: what fraction of vorticity is within r < R_tube? */
    {
        int kk = N / 2;
        double total_sq = 0, core_sq = 0;
        for (int i = 3; i < N-3; i++) {
            double x = -L + i * dx;
            for (int j = 3; j < N-3; j++) {
                double y = -L + j * dx;
                double Ox, Oy, Oz;
                compute_vorticity(i, j, kk, &Ox, &Oy, &Oz);
                double sq = Ox*Ox + Oy*Oy + Oz*Oz;
                total_sq += sq;
                double rperp = sqrt(x*x + y*y);
                if (rperp < 2.0 * R_tube)
                    core_sq += sq;
            }
        }
        double frac = (total_sq > 0) ? core_sq / total_sq : 0;
        printf("  Vorticity localization: %.1f%% within r < 2*R_tube = %.1f\n",
               100.0 * frac, 2.0 * R_tube);
    }
}

/* ================================================================== */
/*  M3b: Torsion flux through cross-sections                           */
/* ================================================================== */
static void run_M3b(void)
{
    printf("\n================================================================\n");
    printf("  M3b: Torsion Flux Through Cross-Sections\n");
    printf("================================================================\n\n");

    char path[600];
    snprintf(path, sizeof(path), "%s/m3b_flux.tsv", outdir);
    FILE *f = fopen(path, "w");
    if (!f) { fprintf(stderr, "Cannot open %s\n", path); return; }
    fprintf(f, "# Torsion flux Phi_T = integral Omega_z dx dy\n");
    fprintf(f, "z_frac\tz_val\tk_idx\tPhi_T\tPhi_T_core\n");

    /* Compute flux at z = 0, L/4, L/2, 3L/4, and a fine z-scan */
    int n_slices = N;  /* scan all z slices for a smooth profile */
    double *flux_total = calloc(n_slices, sizeof(double));
    double *flux_core  = calloc(n_slices, sizeof(double));
    double core_r = 2.0 * R_tube;

    printf("  Computing Omega_z flux at all %d z-slices...\n", n_slices);

    #pragma omp parallel for schedule(dynamic, 4)
    for (int kk = 0; kk < n_slices; kk++) {
        double fl_tot = 0, fl_core = 0;
        double dA = dx * dx;

        for (int i = 3; i < N-3; i++) {
            double x = -L + i * dx;
            for (int j = 3; j < N-3; j++) {
                double y = -L + j * dx;
                double Ox, Oy, Oz;
                compute_vorticity(i, j, kk, &Ox, &Oy, &Oz);
                fl_tot += Oz * dA;
                double rperp = sqrt(x*x + y*y);
                if (rperp < core_r)
                    fl_core += Oz * dA;
            }
        }
        flux_total[kk] = fl_tot;
        flux_core[kk]  = fl_core;
    }

    /* Output all slices */
    for (int kk = 0; kk < n_slices; kk++) {
        double z = -L + kk * dx;
        fprintf(f, "%.4f\t%.4f\t%d\t%.6e\t%.6e\n",
                (double)kk / n_slices, z, kk, flux_total[kk], flux_core[kk]);
    }
    fclose(f);

    /* Print key slices */
    int k_slices[] = {0, N/4, N/2, 3*N/4};
    const char *labels[] = {"z=-L", "z=-L/2", "z=0", "z=L/2"};

    double sum = 0, sum2 = 0;
    for (int kk = 0; kk < n_slices; kk++) {
        sum  += flux_total[kk];
        sum2 += flux_total[kk] * flux_total[kk];
    }
    double mean_flux = sum / n_slices;
    double std_flux  = sqrt(sum2 / n_slices - mean_flux * mean_flux);

    printf("\n  %-10s  %-6s  %-14s  %-14s\n", "Slice", "z", "Phi_T(total)", "Phi_T(core)");
    printf("  %-10s  %-6s  %-14s  %-14s\n", "-----", "---", "------------", "-----------");
    for (int s = 0; s < 4; s++) {
        int kk = k_slices[s];
        double z = -L + kk * dx;
        printf("  %-10s  %6.2f  %14.6e  %14.6e\n",
               labels[s], z, flux_total[kk], flux_core[kk]);
    }
    printf("\n  Mean Phi_T = %.6e\n", mean_flux);
    printf("  Std  Phi_T = %.6e  (%.2f%% variation)\n",
           std_flux, (fabs(mean_flux) > 1e-30) ? 100.0*std_flux/fabs(mean_flux) : 0.0);
    printf("  Phi_T / (2*pi) = %.6f\n", mean_flux / (2.0 * M_PI));
    printf("  Output: %s\n", path);

    free(flux_total);
    free(flux_core);
}

/* ================================================================== */
/*  M3c: Two braids — torsion-mediated force                           */
/* ================================================================== */
static void run_M3c(void)
{
    printf("\n================================================================\n");
    printf("  M3c: Two Braids — Torsion-Mediated Force\n");
    printf("================================================================\n\n");

    /*
     * Test A: same twist (both +k_z) — should they repel?
     * Test B: opposite twist (+k_z and -k_z) — should they attract?
     *
     * Two braids separated by D=30 along x-axis: one at x=-15, one at x=+15.
     * We need L large enough: L=40 so domain is [-40,40].
     */

    double L_save = L;
    double D = 30.0;
    L = 40.0;
    dx  = 2.0 * L / (N - 1);
    dx2 = dx * dx;

    /* Re-allocate for larger grid */
    Ngrid = (long)N * N * N;

    printf("  Two-braid test: D=%.1f, L=%.1f, dx=%.4f\n", D, L, dx);
    printf("  Braid 1 at x=%.1f, Braid 2 at x=+%.1f\n", -D/2, D/2);
    dt = cfl_frac * dx;

    const char *config_labels[] = {"same_twist", "opposite_twist"};
    double twist_signs[][2] = {{+1.0, +1.0}, {+1.0, -1.0}};

    for (int cfg = 0; cfg < 2; cfg++) {
        printf("\n  --- Config: %s ---\n", config_labels[cfg]);

        /* Zero fields */
        for (int a = 0; a < 3; a++) {
            memset(phi[a], 0, Ngrid * sizeof(double));
            memset(vel[a], 0, Ngrid * sizeof(double));
            memset(acc[a], 0, Ngrid * sizeof(double));
        }

        /* Place two braids. Note: init_braid uses += so we can superpose. */
        /* Braid 1 at x = -D/2, y = 0 */
        init_braid(-D/2.0, 0.0, twist_signs[cfg][0]);
        /* Braid 2 at x = +D/2, y = 0 */
        init_braid(+D/2.0, 0.0, twist_signs[cfg][1]);

        setup_damping();

        /* Measure centroid of each braid over time */
        double t_evolve = 200.0;
        int Nt = (int)(t_evolve / dt) + 1;
        int meas_every = Nt / 100;
        if (meas_every < 1) meas_every = 1;

        char tpath[600];
        snprintf(tpath, sizeof(tpath), "%s/m3c_%s.tsv", outdir, config_labels[cfg]);
        FILE *ft = fopen(tpath, "w");
        if (!ft) { fprintf(stderr, "Cannot open %s\n", tpath); continue; }
        fprintf(ft, "time\tcx1\tcx2\tsep\tE_total\tpeak_P\n");

        compute_acc();
        double wall_start = omp_get_wtime();
        int print_every = Nt / 10;
        if (print_every < 1) print_every = 1;

        double cx1_init = 0, cx2_init = 0;

        for (int n = 0; n <= Nt; n++) {
            double t = n * dt;

            if (n % meas_every == 0 || n == Nt) {
                /* Compute energy-weighted centroid for left (x<0) and right (x>0) halves */
                double sum_xE_left = 0, sum_E_left = 0;
                double sum_xE_right = 0, sum_E_right = 0;
                double total_E = 0;
                double peak_P = 0;

                #pragma omp parallel
                {
                    double lxEl = 0, lEl = 0, lxEr = 0, lEr = 0, ltE = 0, lpP = 0;

                    #pragma omp for schedule(static) nowait
                    for (int i = 2; i < N-2; i++) {
                        double x = -L + i * dx;
                        for (int j = 2; j < N-2; j++) {
                            for (int k = 0; k < N; k++) {
                                long idx = IDX(i, j, k);
                                int km1 = wrap_z(k - 1);
                                int kp1 = wrap_z(k + 1);

                                double e_loc = 0;
                                for (int a = 0; a < 3; a++) {
                                    e_loc += 0.5 * vel[a][idx] * vel[a][idx];
                                    double gx = (phi[a][IDX(i+1,j,k)] - phi[a][IDX(i-1,j,k)]) / (2*dx);
                                    double gy = (phi[a][IDX(i,j+1,k)] - phi[a][IDX(i,j-1,k)]) / (2*dx);
                                    double gz = (phi[a][IDX(i,j,kp1)] - phi[a][IDX(i,j,km1)]) / (2*dx);
                                    e_loc += 0.5 * (gx*gx + gy*gy + gz*gz);
                                }
                                double p0 = phi[0][idx], p1 = phi[1][idx], p2 = phi[2][idx];
                                double Pv = p0 * p1 * p2;
                                (void)Pv;
                                /* Use |e_loc| for centroid weighting (potential can be negative) */
                                double w = fabs(e_loc);

                                if (x < 0) {
                                    lxEl += x * w;
                                    lEl  += w;
                                } else {
                                    lxEr += x * w;
                                    lEr  += w;
                                }
                                ltE += e_loc;
                                double aP = fabs(Pv);
                                if (aP > lpP) lpP = aP;
                            }
                        }
                    }

                    #pragma omp critical
                    {
                        sum_xE_left += lxEl;   sum_E_left += lEl;
                        sum_xE_right += lxEr;  sum_E_right += lEr;
                        total_E += ltE;
                        if (lpP > peak_P) peak_P = lpP;
                    }
                }

                double cx1 = (sum_E_left > 1e-20)  ? sum_xE_left / sum_E_left   : -D/2;
                double cx2 = (sum_E_right > 1e-20) ? sum_xE_right / sum_E_right : +D/2;
                double sep = cx2 - cx1;

                if (n == 0) {
                    cx1_init = cx1;
                    cx2_init = cx2;
                }

                fprintf(ft, "%.4f\t%.6f\t%.6f\t%.6f\t%.6e\t%.6e\n",
                        t, cx1, cx2, sep, total_E * dx * dx * dx, peak_P);

                if (n % print_every == 0 || n == Nt) {
                    double elapsed = omp_get_wtime() - wall_start;
                    printf("  t=%6.1f  cx1=%7.3f  cx2=%7.3f  sep=%7.3f  |P|=%.4f  [%.0fs]\n",
                           t, cx1, cx2, sep, peak_P, elapsed);
                    fflush(stdout);
                }
            }

            if (n == Nt) break;
            verlet_step();
        }

        fclose(ft);

        /* Measure final centroid shift */
        /* Re-compute final state */
        double sum_xE_left = 0, sum_E_left = 0;
        double sum_xE_right = 0, sum_E_right = 0;

        #pragma omp parallel
        {
            double lxEl = 0, lEl = 0, lxEr = 0, lEr = 0;

            #pragma omp for schedule(static) nowait
            for (int i = 2; i < N-2; i++) {
                double x = -L + i * dx;
                for (int j = 2; j < N-2; j++) {
                    for (int k = 0; k < N; k++) {
                        long idx = IDX(i, j, k);
                        int km1 = wrap_z(k - 1);
                        int kp1 = wrap_z(k + 1);

                        double e_loc = 0;
                        for (int a = 0; a < 3; a++) {
                            e_loc += 0.5 * vel[a][idx] * vel[a][idx];
                            double gx = (phi[a][IDX(i+1,j,k)] - phi[a][IDX(i-1,j,k)]) / (2*dx);
                            double gy = (phi[a][IDX(i,j+1,k)] - phi[a][IDX(i,j-1,k)]) / (2*dx);
                            double gz = (phi[a][IDX(i,j,kp1)] - phi[a][IDX(i,j,km1)]) / (2*dx);
                            e_loc += 0.5 * (gx*gx + gy*gy + gz*gz);
                        }
                        double w = fabs(e_loc);
                        if (x < 0) { lxEl += x * w; lEl += w; }
                        else        { lxEr += x * w; lEr += w; }
                    }
                }
            }

            #pragma omp critical
            {
                sum_xE_left += lxEl;   sum_E_left += lEl;
                sum_xE_right += lxEr;  sum_E_right += lEr;
            }
        }

        double cx1_final = (sum_E_left > 1e-20)  ? sum_xE_left / sum_E_left   : -D/2;
        double cx2_final = (sum_E_right > 1e-20) ? sum_xE_right / sum_E_right : +D/2;
        double sep_init  = cx2_init - cx1_init;
        double sep_final = cx2_final - cx1_final;
        double delta_sep = sep_final - sep_init;

        printf("\n  %s: sep_init=%.3f  sep_final=%.3f  delta=%.3f  (%s)\n",
               config_labels[cfg], sep_init, sep_final, delta_sep,
               delta_sep > 0.1 ? "REPULSION" : (delta_sep < -0.1 ? "ATTRACTION" : "NEUTRAL"));
        printf("  Output: %s\n", tpath);
    }

    /* Restore */
    L = L_save;
    dx  = 2.0 * L / (N - 1);
    dx2 = dx * dx;
    dt  = cfl_frac * dx;
}

/* ─── Main ─── */
int main(int argc, char **argv)
{
    int test = -1;  /* -1=all, 0=M3a, 1=M3b, 2=M3c */
    for (int i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "-test") && i+1 < argc)    { test = atoi(argv[++i]); }
        else if (!strcmp(argv[i], "-N") && i+1 < argc)  { N = atoi(argv[++i]); }
        else if (!strcmp(argv[i], "-L") && i+1 < argc)  { L = atof(argv[++i]); }
        else if (!strcmp(argv[i], "-mu") && i+1 < argc)  { mu_pot = atof(argv[++i]); }
        else if (!strcmp(argv[i], "-kappa") && i+1 < argc) { kappa = atof(argv[++i]); }
        else if (!strcmp(argv[i], "-tfinal") && i+1 < argc) { tfinal = atof(argv[++i]); }
        else if (!strcmp(argv[i], "-cfl") && i+1 < argc) { cfl_frac = atof(argv[++i]); }
        else if (!strcmp(argv[i], "-o") && i+1 < argc)  { strncpy(outdir, argv[++i], sizeof(outdir)-1); }
    }

    /* Use L=40 for M3c allocation (largest domain needed) */
    double L_alloc = 40.0;
    dx  = 2.0 * L_alloc / (N - 1);
    dx2 = dx * dx;
    m2  = mass * mass;
    dt  = cfl_frac * dx;
    Ngrid = (long)N * N * N;

    printf("=== V27-M3: Electromagnetic Structure from Torsion ===\n");
    printf("  m=%.1f  mu=%.1f  kappa=%.1f  A0=%.1f  R_tube=%.1f\n",
           mass, mu_pot, kappa, A0, R_tube);
    printf("  N=%d  L=%.1f (single) / 40.0 (two-braid)\n", N, L);
    printf("  Ngrid=%.1fM  Memory=%.1f MB\n", Ngrid/1e6, Ngrid*8.0*10/1e6);
    printf("  Threads: %d\n", omp_get_max_threads());
    printf("  test=%d (-1=all, 0=M3a, 1=M3b, 2=M3c)\n\n", test);
    fflush(stdout);

    /* Allocate for largest grid */
    for (int a = 0; a < 3; a++) {
        phi[a] = calloc(Ngrid, sizeof(double));
        vel[a] = calloc(Ngrid, sizeof(double));
        acc[a] = calloc(Ngrid, sizeof(double));
        if (!phi[a] || !vel[a] || !acc[a]) {
            fprintf(stderr, "Allocation failed\n");
            return 1;
        }
    }
    damp = malloc(Ngrid * sizeof(double));
    if (!damp) { fprintf(stderr, "Damp allocation failed\n"); return 1; }

    /* ─── M3a + M3b: Single braid, evolve then measure ─── */
    if (test == -1 || test == 0 || test == 1) {
        /* Reset to single-braid domain */
        dx  = 2.0 * L / (N - 1);
        dx2 = dx * dx;
        dt  = cfl_frac * dx;

        printf("================================================================\n");
        printf("  Evolving single braid to t=%.0f (settling)...\n", tfinal);
        printf("  dx=%.4f  dt=%.5f\n", dx, dt);
        printf("================================================================\n\n");

        /* Zero and init */
        for (int a = 0; a < 3; a++) {
            memset(phi[a], 0, Ngrid * sizeof(double));
            memset(vel[a], 0, Ngrid * sizeof(double));
            memset(acc[a], 0, Ngrid * sizeof(double));
        }
        init_braid(0.0, 0.0, +1.0);
        setup_damping();

        int Nt_settle = (int)(tfinal / dt) + 1;
        int pe = Nt_settle / 10;
        if (pe < 1) pe = 1;
        evolve(Nt_settle, pe);

        Diag d = compute_diag();
        printf("\n  After settling: E=%.2f  |P|=%.4f  Pz=%.4f\n\n",
               d.Et, d.peak_P, d.Pz);

        if (test == -1 || test == 0)
            run_M3a();
        if (test == -1 || test == 1)
            run_M3b();
    }

    /* ─── M3c: Two braids ─── */
    if (test == -1 || test == 2) {
        run_M3c();
    }

    /* Cleanup */
    for (int a = 0; a < 3; a++) {
        free(phi[a]); free(vel[a]); free(acc[a]);
    }
    free(damp);

    printf("\n=== V27-M3 Complete ===\n");
    return 0;
}
