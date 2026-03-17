/*
 * m4.c — V27-M4: Mass from Dynamics
 *
 * Based on V26-DynA propagating braid (dynA.c).
 *
 * Tests:
 *   M4a: Mass reduction scan m={1.0,0.8,0.6,0.4,0.2,0.1,0.0}
 *   M4b: Strong binding at m=0, mu={-20,-50,-100,-200}, kappa=|mu|
 *   M4c: Pairwise coupling at low m with lambda_pw=0.5
 *   M4d: Dispersion relation inside vs outside braid at m=1
 *
 * Compile: gcc -O3 -fopenmp -Wall -o m4 src/m4.c -lm
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>

/* ─── Parameters ─── */
static double mu_pot    = -20.0;
static double kappa     = 20.0;
static double mass      = 1.0;
static double A0        = 0.8;
static double R_tube    = 3.0;
static double lambda_pw = 0.0;   /* pairwise coupling */

static int    N         = 128;
static double L         = 20.0;
static double tfinal    = 500.0;
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
static void init_braid(void)
{
    double Lz = 2.0 * L;
    double k_wave = 2.0 * M_PI / Lz;
    double omega = sqrt(k_wave * k_wave + mass * mass);

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
            double phase = k_wave * z + 2.0 * M_PI * a / 3.0;
            phi[a][idx] = envelope * cos(phase);
            vel[a][idx] = omega * envelope * sin(phase);
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
                    double dVdphi_triple = mu_pot * P * dP / denom2;

                    /* Pairwise coupling force: -d/dphi_a [lambda_pw*(phi1*phi2+phi2*phi3+phi3*phi1)] */
                    double dVdphi_pw = 0.0;
                    if (lambda_pw != 0.0) {
                        switch (a) {
                            case 0: dVdphi_pw = lambda_pw * (p1 + p2); break;
                            case 1: dVdphi_pw = lambda_pw * (p0 + p2); break;
                            default: dVdphi_pw = lambda_pw * (p0 + p1); break;
                        }
                    }

                    acc[a][idx] = lapl - m2 * phi[a][idx] - dVdphi_triple - dVdphi_pw;
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

/* ─── Diagnostics structure ─── */
typedef struct {
    double Ek, Eg, Em, Ep, Epw, Et;
    double fc;
    double peak_P;
    double Pz;
} Diag;

/* ─── Compute diagnostics ─── */
static Diag compute_diag(double core_radius)
{
    Diag d;
    memset(&d, 0, sizeof(d));

    double Ek = 0, Eg = 0, Em = 0, Ep = 0, Epw = 0;
    double Ecore = 0, Eall = 0;
    double peak_P = 0;
    double Pz_total = 0;

    #pragma omp parallel
    {
        double lEk = 0, lEg = 0, lEm = 0, lEp = 0, lEpw = 0, lEc = 0, lEa = 0;
        double lpkP = 0;
        double lPz = 0;

        #pragma omp for schedule(static) nowait
        for (int i = 2; i < N-2; i++) {
            double x = -L + i * dx;
            for (int j = 2; j < N-2; j++) {
                double y = -L + j * dx;
                for (int k = 0; k < N; k++) {
                    double z = -L + k * dx;
                    long idx = IDX(i, j, k);
                    double dV = dx * dx * dx;
                    double e_loc = 0;

                    int km1 = wrap_z(k - 1);
                    int kp1 = wrap_z(k + 1);

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

                        lPz += -vel[a][idx] * gz * dV;
                    }

                    double p0 = phi[0][idx], p1 = phi[1][idx], p2 = phi[2][idx];
                    double P_val = p0 * p1 * p2;
                    double P2 = P_val * P_val;
                    double Vloc = 0.5 * mu_pot * P2 / (1.0 + kappa * P2);
                    lEp += Vloc * dV;
                    e_loc += Vloc;

                    /* Pairwise energy */
                    if (lambda_pw != 0.0) {
                        double pw_e = lambda_pw * (p0*p1 + p1*p2 + p2*p0);
                        lEpw += pw_e * dV;
                        e_loc += pw_e;
                    }

                    double absP = fabs(P_val);
                    if (absP > lpkP) lpkP = absP;

                    lEa += e_loc * dV;
                    double r = sqrt(x*x + y*y + z*z);
                    if (r < core_radius) lEc += e_loc * dV;
                }
            }
        }

        #pragma omp critical
        {
            Ek += lEk; Eg += lEg; Em += lEm; Ep += lEp; Epw += lEpw;
            Ecore += lEc; Eall += lEa;
            Pz_total += lPz;
            if (lpkP > peak_P) peak_P = lpkP;
        }
    }

    d.Ek = Ek; d.Eg = Eg; d.Em = Em; d.Ep = Ep; d.Epw = Epw;
    d.Et = Ek + Eg + Em + Ep + Epw;
    d.fc = (fabs(Eall) > 1e-20) ? Ecore / Eall : 0.0;
    d.peak_P = peak_P;
    d.Pz = Pz_total;

    return d;
}

/* ─── Run one evolution ─── */
typedef struct {
    double fc_final, peakP_final, E_final, Pz_final;
} RunResult;

static RunResult run_evolution(const char *label, double run_tfinal)
{
    RunResult res = {0};

    m2 = mass * mass;

    /* Zero fields */
    for (int a = 0; a < 3; a++) {
        memset(phi[a], 0, Ngrid * sizeof(double));
        memset(vel[a], 0, Ngrid * sizeof(double));
        memset(acc[a], 0, Ngrid * sizeof(double));
    }

    init_braid();
    setup_damping();

    double core_radius = 8.0;
    int Nt = (int)(run_tfinal / dt) + 1;
    int print_every = Nt / 20;
    if (print_every < 1) print_every = 1;
    int diag_every = Nt / 2000;
    if (diag_every < 1) diag_every = 1;

    /* Output file */
    char path[600];
    snprintf(path, sizeof(path), "%s/m4_%s.tsv", outdir, label);
    FILE *fout = fopen(path, "w");
    if (!fout) { fprintf(stderr, "Cannot open %s\n", path); return res; }
    fprintf(fout, "time\tE_total\tE_kin\tE_grad\tE_mass\tE_pot\tE_pw\t"
                  "fc\tpeak_P\tPz\n");

    compute_acc();
    double wall_start = omp_get_wtime();

    Diag d0 = compute_diag(core_radius);
    printf("  t=%7.1f  E=%.2f  fc=%.4f  |P|=%.6f  Pz=%.4f\n",
           0.0, d0.Et, d0.fc, d0.peak_P, d0.Pz);

    for (int n = 0; n <= Nt; n++) {
        double t = n * dt;

        int do_diag  = (n % diag_every == 0);
        int do_print = (n % print_every == 0);

        if (do_diag || do_print) {
            Diag d = compute_diag(core_radius);

            if (do_diag)
                fprintf(fout, "%.4f\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t"
                              "%.6f\t%.6e\t%.6e\n",
                        t, d.Et, d.Ek, d.Eg, d.Em, d.Ep, d.Epw,
                        d.fc, d.peak_P, d.Pz);

            if (do_print) {
                double elapsed = omp_get_wtime() - wall_start;
                double frac = (double)n / Nt;
                double eta_t = (frac > 0.001) ? elapsed * (1.0-frac)/frac : 0;
                printf("  t=%7.1f  E=%.2f  fc=%.4f  |P|=%.6f  Pz=%.4f  [%.0fs, ETA %.0fs]\n",
                       t, d.Et, d.fc, d.peak_P, d.Pz, elapsed, eta_t);
                fflush(stdout);
            }
        }

        if (n == Nt) break;
        verlet_step();
    }

    fclose(fout);

    Diag dfinal = compute_diag(core_radius);
    res.fc_final = dfinal.fc;
    res.peakP_final = dfinal.peak_P;
    res.E_final = dfinal.Et;
    res.Pz_final = dfinal.Pz;

    printf("  Final: E=%.2f  fc=%.4f  |P|=%.6f  Pz=%.4f  (%.1fs)\n",
           dfinal.Et, dfinal.fc, dfinal.peak_P, dfinal.Pz,
           omp_get_wtime() - wall_start);

    return res;
}

/* ─── M4d: Add perturbation and track dispersion ─── */
static void run_dispersion(void)
{
    printf("\n================================================================\n");
    printf("  M4d: Dispersion Relation Inside vs Outside Braid\n");
    printf("================================================================\n\n");

    mass = 1.0;
    m2 = mass * mass;
    mu_pot = -20.0;
    kappa = 20.0;
    lambda_pw = 0.0;

    /* Zero fields */
    for (int a = 0; a < 3; a++) {
        memset(phi[a], 0, Ngrid * sizeof(double));
        memset(vel[a], 0, Ngrid * sizeof(double));
        memset(acc[a], 0, Ngrid * sizeof(double));
    }

    init_braid();
    setup_damping();

    /* Store the base braid state */
    double *phi_base[3];
    for (int a = 0; a < 3; a++) {
        phi_base[a] = malloc(Ngrid * sizeof(double));
        memcpy(phi_base[a], phi[a], Ngrid * sizeof(double));
    }

    /*
     * We'll track perturbations at two probe points:
     *   Inside:  (0, 0, 0) — center of the braid tube
     *   Outside: (12, 0, 0) — well outside the tube (R_tube=3)
     *
     * Add small Gaussian perturbation at each location, run separately.
     * Track phi(probe, t) to extract oscillation frequency.
     * The perturbation is localized in z, so it excites a range of k modes.
     * The dominant k is set by the perturbation width.
     */

    /* Perturbation parameters */
    double pert_amp = 0.01;  /* small compared to A0=0.8 */
    (void)0; /* pert localization width = 1.5 (applied via grid-scale Gaussian) */

    /* Probe locations (grid indices) */
    int probe_in_i = N/2, probe_in_j = N/2;     /* (0,0,z) */
    int probe_out_i = (int)((12.0 + L)/dx + 0.5), probe_out_j = N/2; /* (12,0,z) */

    /* We'll do k-scan: add perturbation with definite k, measure omega */
    int n_k = 8;
    double k_vals[8] = {0.2, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0};

    char path[600];
    snprintf(path, sizeof(path), "%s/m4_dispersion.tsv", outdir);
    FILE *fdisp = fopen(path, "w");
    if (!fdisp) { fprintf(stderr, "Cannot open dispersion output\n"); return; }
    fprintf(fdisp, "# Dispersion relation: omega vs k inside and outside braid\n");
    fprintf(fdisp, "k\tomega_in\tomega_out\tomega_theory_out\n");

    double t_disp = 80.0;  /* short evolution for frequency extraction */
    int Nt_disp = (int)(t_disp / dt) + 1;

    /* DFT storage */
    int max_pts = 20000;
    double *phi_hist_in  = malloc(max_pts * sizeof(double));
    double *phi_hist_out = malloc(max_pts * sizeof(double));
    double *t_hist       = malloc(max_pts * sizeof(double));

    for (int ik = 0; ik < n_k; ik++) {
        double kz = k_vals[ik];
        double omega_theory = sqrt(kz*kz + m2);

        printf("  k=%.2f (omega_theory=%.4f)...", kz, omega_theory);
        fflush(stdout);

        /* Reset to base braid */
        for (int a = 0; a < 3; a++) {
            memcpy(phi[a], phi_base[a], Ngrid * sizeof(double));
            memset(vel[a], 0, Ngrid * sizeof(double));
            memset(acc[a], 0, Ngrid * sizeof(double));
        }

        /* Re-initialize velocities for propagating braid */
        {
            double Lz = 2.0 * L;
            double k_wave = 2.0 * M_PI / Lz;
            double omega = sqrt(k_wave * k_wave + m2);
            #pragma omp parallel for schedule(static)
            for (long idx = 0; idx < Ngrid; idx++) {
                int ii = idx / (N * N);
                int jj = (idx / N) % N;
                int kk = idx % N;
                double x = -L + ii * dx;
                double y = -L + jj * dx;
                double z = -L + kk * dx;
                double r_perp = sqrt(x*x + y*y);
                double envelope = A0 * exp(-r_perp*r_perp / (2.0*R_tube*R_tube));
                for (int a = 0; a < 3; a++) {
                    double phase = k_wave * z + 2.0 * M_PI * a / 3.0;
                    vel[a][idx] = omega * envelope * sin(phase);
                }
            }
        }

        /* Add perturbation on phi[0] only, with definite kz */
        /* Inside: along z-axis at (0,0,z) */
        /* Outside: along z-axis at (12,0,z) */
        #pragma omp parallel for schedule(static)
        for (int k = 0; k < N; k++) {
            double z = -L + k * dx;
            double pert = pert_amp * sin(kz * z);

            /* Inside perturbation */
            long idx_in = IDX(probe_in_i, probe_in_j, k);
            phi[0][idx_in] += pert;

            /* Outside perturbation (spread over a few points for smoothness) */
            for (int di = -2; di <= 2; di++) {
                for (int dj = -2; dj <= 2; dj++) {
                    int ii = probe_out_i + di;
                    int jj = probe_out_j + dj;
                    if (ii >= 2 && ii < N-2 && jj >= 2 && jj < N-2) {
                        double w = exp(-(di*di + dj*dj) / 2.0);
                        long idx_out = IDX(ii, jj, k);
                        phi[0][idx_out] += pert * w;
                    }
                }
            }
        }

        compute_acc();

        /* Evolve and record probe signals */
        int n_pts = 0;
        int rec_every = Nt_disp / max_pts;
        if (rec_every < 1) rec_every = 1;

        for (int n = 0; n <= Nt_disp; n++) {
            if (n % rec_every == 0 && n_pts < max_pts) {
                /* Record phi[0] at probe points, z=N/2 (middle of domain) */
                long idx_in  = IDX(probe_in_i, probe_in_j, N/2);
                long idx_out = IDX(probe_out_i, probe_out_j, N/2);
                /* Subtract base to get perturbation */
                phi_hist_in[n_pts]  = phi[0][idx_in]  - phi_base[0][idx_in];
                phi_hist_out[n_pts] = phi[0][idx_out] - phi_base[0][idx_out];
                t_hist[n_pts] = n * dt;
                n_pts++;
            }
            if (n == Nt_disp) break;
            verlet_step();
        }

        /* DFT to find peak omega for inside and outside */
        int start = n_pts / 3;  /* skip transient */
        double omega_in = 0, omega_out = 0;

        /* Scan omega from 0 to 8 */
        int nf = 400;
        double omega_max = 8.0;
        double best_pow_in = 0, best_pow_out = 0;

        double mean_in = 0, mean_out = 0;
        for (int j = start; j < n_pts; j++) {
            mean_in  += phi_hist_in[j];
            mean_out += phi_hist_out[j];
        }
        mean_in  /= (n_pts - start);
        mean_out /= (n_pts - start);

        for (int kk = 1; kk < nf; kk++) {
            double om = omega_max * kk / nf;
            double re_in = 0, im_in = 0, re_out = 0, im_out = 0;
            for (int j = start; j < n_pts; j++) {
                double dtj = (j > start) ?
                    (t_hist[j]-t_hist[j-1]) : (t_hist[start+1]-t_hist[start]);
                double c = cos(om * t_hist[j]);
                double s = sin(om * t_hist[j]);
                re_in  += (phi_hist_in[j]  - mean_in)  * c * dtj;
                im_in  += (phi_hist_in[j]  - mean_in)  * s * dtj;
                re_out += (phi_hist_out[j] - mean_out) * c * dtj;
                im_out += (phi_hist_out[j] - mean_out) * s * dtj;
            }
            double pw_in  = re_in*re_in + im_in*im_in;
            double pw_out = re_out*re_out + im_out*im_out;
            if (pw_in  > best_pow_in)  { best_pow_in  = pw_in;  omega_in  = om; }
            if (pw_out > best_pow_out) { best_pow_out = pw_out; omega_out = om; }
        }

        printf(" omega_in=%.3f  omega_out=%.3f  (theory=%.3f)\n",
               omega_in, omega_out, omega_theory);

        fprintf(fdisp, "%.4f\t%.6f\t%.6f\t%.6f\n", kz, omega_in, omega_out, omega_theory);
        fflush(fdisp);
    }

    fclose(fdisp);

    /* Compute effective mass from dispersion fit */
    /* omega^2 = k^2 + m_eff^2, so m_eff^2 = omega^2 - k^2 */
    printf("\n  Effective mass extraction (omega^2 - k^2):\n");
    printf("  %-6s  %-10s  %-10s  %-10s  %-10s\n",
           "k", "m2_eff_in", "m2_eff_out", "m_eff_in", "m_eff_out");

    /* Re-read the file we just wrote */
    snprintf(path, sizeof(path), "%s/m4_dispersion.tsv", outdir);
    FILE *fread = fopen(path, "r");
    if (fread) {
        char line[256];
        /* Write a summary file */
        char path2[600];
        snprintf(path2, sizeof(path2), "%s/m4_meff.tsv", outdir);
        FILE *fmeff = fopen(path2, "w");
        if (fmeff)
            fprintf(fmeff, "k\tm2_eff_in\tm2_eff_out\tm_eff_in\tm_eff_out\n");

        while (fgets(line, sizeof(line), fread)) {
            if (line[0] == '#') continue;
            double kk, oi, oo, ot;
            if (sscanf(line, "%lf\t%lf\t%lf\t%lf", &kk, &oi, &oo, &ot) == 4) {
                double m2_in  = oi*oi - kk*kk;
                double m2_out = oo*oo - kk*kk;
                double meff_in  = (m2_in  > 0) ? sqrt(m2_in)  : -sqrt(-m2_in);
                double meff_out = (m2_out > 0) ? sqrt(m2_out) : -sqrt(-m2_out);
                printf("  %-6.2f  %-10.4f  %-10.4f  %-10.4f  %-10.4f\n",
                       kk, m2_in, m2_out, meff_in, meff_out);
                if (fmeff)
                    fprintf(fmeff, "%.4f\t%.6f\t%.6f\t%.6f\t%.6f\n",
                            kk, m2_in, m2_out, meff_in, meff_out);
            }
        }
        fclose(fread);
        if (fmeff) fclose(fmeff);
    }

    free(phi_hist_in);
    free(phi_hist_out);
    free(t_hist);
    for (int a = 0; a < 3; a++) free(phi_base[a]);
}

/* ─── Main ─── */
int main(int argc, char **argv)
{
    /* Parse args */
    int test = -1;  /* -1 = all, 0=M4a, 1=M4b, 2=M4c, 3=M4d */
    for (int i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "-test") && i+1 < argc) { test = atoi(argv[++i]); }
        else if (!strcmp(argv[i], "-N")  && i+1 < argc) { N = atoi(argv[++i]); }
        else if (!strcmp(argv[i], "-L")  && i+1 < argc) { L = atof(argv[++i]); }
        else if (!strcmp(argv[i], "-o")  && i+1 < argc) { strncpy(outdir, argv[++i], sizeof(outdir)-1); }
        else if (!strcmp(argv[i], "-tfinal") && i+1 < argc) { tfinal = atof(argv[++i]); }
        else if (!strcmp(argv[i], "-cfl") && i+1 < argc) { cfl_frac = atof(argv[++i]); }
    }

    dx  = 2.0 * L / (N - 1);
    dx2 = dx * dx;
    dt  = cfl_frac * dx;
    Ngrid = (long)N * N * N;

    printf("=== V27-M4: Mass from Dynamics ===\n");
    printf("  N=%d  L=%.1f  dx=%.4f  dt=%.5f  Ngrid=%.1fM\n",
           N, L, dx, dt, Ngrid/1e6);
    printf("  Threads: %d\n", omp_get_max_threads());
    printf("  tfinal=%.0f  test=%d (-1=all)\n\n", tfinal, test);
    fflush(stdout);

    /* Allocate */
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

    /* ============================================================ */
    /* M4a: Mass reduction scan                                      */
    /* ============================================================ */
    if (test == -1 || test == 0) {
        printf("================================================================\n");
        printf("  M4a: Mass Reduction Scan\n");
        printf("================================================================\n\n");

        double m_vals[] = {1.0, 0.8, 0.6, 0.4, 0.2, 0.1, 0.0};
        int n_m = 7;

        char path[600];
        snprintf(path, sizeof(path), "%s/m4a_summary.tsv", outdir);
        FILE *fsum = fopen(path, "w");
        if (fsum)
            fprintf(fsum, "mass\tfc\tpeak_P\tE_total\tPz\tsurvived\n");

        mu_pot = -20.0;
        kappa = 20.0;
        lambda_pw = 0.0;

        double m_critical = -1.0;

        for (int im = 0; im < n_m; im++) {
            mass = m_vals[im];
            m2 = mass * mass;

            char label[64];
            snprintf(label, sizeof(label), "m4a_m%.2f", mass);
            printf("\n--- m=%.2f ---\n", mass);

            RunResult r = run_evolution(label, tfinal);

            int survived = (r.peakP_final > 0.01);
            if (survived && (m_critical < 0 || mass < m_critical))
                m_critical = mass;

            printf("  Survived (|P|>0.01)? %s\n\n", survived ? "YES" : "NO");

            if (fsum)
                fprintf(fsum, "%.2f\t%.6f\t%.6e\t%.6e\t%.6e\t%d\n",
                        mass, r.fc_final, r.peakP_final, r.E_final, r.Pz_final, survived);
        }

        if (fsum) fclose(fsum);

        printf("\n  >>> M4a RESULT: m_critical = %.2f (minimum m for |P|>0.01 at t=%.0f)\n\n",
               m_critical, tfinal);
    }

    /* ============================================================ */
    /* M4b: Strong binding compensation at m=0                       */
    /* ============================================================ */
    if (test == -1 || test == 1) {
        printf("================================================================\n");
        printf("  M4b: Strong Binding Compensation (m=0)\n");
        printf("================================================================\n\n");

        double mu_vals[] = {-20.0, -50.0, -100.0, -200.0};
        int n_mu = 4;

        char path[600];
        snprintf(path, sizeof(path), "%s/m4b_summary.tsv", outdir);
        FILE *fsum = fopen(path, "w");
        if (fsum)
            fprintf(fsum, "mu\tkappa\tfc\tpeak_P\tE_total\tPz\tsurvived\n");

        mass = 0.0;
        m2 = 0.0;
        lambda_pw = 0.0;

        for (int im = 0; im < n_mu; im++) {
            mu_pot = mu_vals[im];
            kappa = fabs(mu_pot);  /* V_max = mu^2/(4*kappa) = |mu|/4 when kappa=|mu| */

            char label[64];
            snprintf(label, sizeof(label), "m4b_mu%.0f", fabs(mu_pot));
            printf("\n--- m=0, mu=%.0f, kappa=%.0f ---\n", mu_pot, kappa);

            RunResult r = run_evolution(label, tfinal);

            int survived = (r.peakP_final > 0.01);
            printf("  Survived (|P|>0.01)? %s\n\n", survived ? "YES" : "NO");

            if (fsum)
                fprintf(fsum, "%.1f\t%.1f\t%.6f\t%.6e\t%.6e\t%.6e\t%d\n",
                        mu_pot, kappa, r.fc_final, r.peakP_final, r.E_final, r.Pz_final, survived);
        }

        if (fsum) fclose(fsum);
    }

    /* ============================================================ */
    /* M4c: Pairwise + propagation at low mass                       */
    /* ============================================================ */
    if (test == -1 || test == 2) {
        printf("================================================================\n");
        printf("  M4c: Pairwise Coupling at Low Mass\n");
        printf("================================================================\n\n");

        /* Use a few low mass values around the critical region */
        double m_vals[] = {0.4, 0.2, 0.1, 0.0};
        int n_m = 4;

        char path[600];
        snprintf(path, sizeof(path), "%s/m4c_summary.tsv", outdir);
        FILE *fsum = fopen(path, "w");
        if (fsum)
            fprintf(fsum, "mass\tlambda_pw\tfc\tpeak_P\tE_total\tPz\tsurvived\n");

        mu_pot = -20.0;
        kappa = 20.0;
        lambda_pw = 0.5;

        for (int im = 0; im < n_m; im++) {
            mass = m_vals[im];
            m2 = mass * mass;

            char label[64];
            snprintf(label, sizeof(label), "m4c_m%.2f_pw0.5", mass);
            printf("\n--- m=%.2f, lambda_pw=%.1f ---\n", mass, lambda_pw);

            RunResult r = run_evolution(label, tfinal);

            int survived = (r.peakP_final > 0.01);
            printf("  Survived (|P|>0.01)? %s\n\n", survived ? "YES" : "NO");

            if (fsum)
                fprintf(fsum, "%.2f\t%.1f\t%.6f\t%.6e\t%.6e\t%.6e\t%d\n",
                        mass, lambda_pw, r.fc_final, r.peakP_final, r.E_final, r.Pz_final, survived);
        }

        if (fsum) fclose(fsum);

        lambda_pw = 0.0;  /* reset */
    }

    /* ============================================================ */
    /* M4d: Dispersion relation inside vs outside                    */
    /* ============================================================ */
    if (test == -1 || test == 3) {
        run_dispersion();
    }

    /* Cleanup */
    for (int a = 0; a < 3; a++) {
        free(phi[a]); free(vel[a]); free(acc[a]);
    }
    free(damp);

    printf("\n=== V27-M4 Complete ===\n");
    return 0;
}
