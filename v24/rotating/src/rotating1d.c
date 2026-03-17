/*
 * rotating1d.c — Rotating triad: frequency-split three-body oscillon
 *
 * V24-DG: Split frequencies ω₁=ω₀+δ, ω₂=ω₀, ω₃=ω₀-δ
 * The triad rotates in internal field space at rate δ.
 *
 * Protocol:
 *   Phase 1: Equilibrate standard oscillon for t_equil (all fields equal)
 *   Phase 2: Apply phase perturbation v₁ += δ·φ₁, v₃ -= δ·φ₃
 *   Phase 3: Evolve for t_run, measure energy loss per field
 *
 * Key question: does δ>0 produce lower dE/dt than δ=0?
 *
 * Compile: gcc -O3 -Wall -o rotating1d src/rotating1d.c -lm
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Parameters (proposal defaults) */
static double mu     = -20.0;
static double kappa  = 20.0;
static double mass   = 1.0;
static double A_init = 0.8;
static double sigma  = 3.0;
static int    Nx     = 4000;
static double xmax   = 100.0;
static double t_equil = 5000.0;
static double t_run   = 10000.0;
static double delta   = 0.0;     /* frequency splitting */
static char   outdir[512] = "v24/rotating/data";

static void parse_args(int argc, char **argv)
{
    for (int i = 1; i < argc - 1; i += 2) {
        if      (!strcmp(argv[i], "-mu"))     mu     = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-kappa"))  kappa  = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-mass"))   mass   = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-A"))      A_init = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-sigma"))  sigma  = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-Nx"))     Nx     = atoi(argv[i+1]);
        else if (!strcmp(argv[i], "-xmax"))   xmax   = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-t_equil")) t_equil = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-t_run"))  t_run  = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-delta"))  delta  = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-o"))      strncpy(outdir, argv[i+1], sizeof(outdir)-1);
        else { fprintf(stderr, "Unknown flag: %s\n", argv[i]); exit(1); }
    }
}

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

/* Compute energy in a region [0, rcut] (symmetric, so half-domain from center) */
static void compute_energies(double *phi[3], double *vel[3], int Nx, double dx,
                             double m2, double xmax_val,
                             double *Ek_out, double *Eg_out, double *Em_out,
                             double *Ep_out, double *Et_out, double *Ecore_out,
                             double *Etot_out, double peak_out[3])
{
    double Ek = 0, Eg = 0, Em = 0, Ep = 0;
    double Ecore = 0, Eall = 0;
    double core_r = 3.0 * sigma;
    peak_out[0] = peak_out[1] = peak_out[2] = 0;

    for (int i = 1; i < Nx - 1; i++) {
        double x = -xmax_val + i * dx;
        double e_loc = 0;
        for (int a = 0; a < 3; a++) {
            double v2 = vel[a][i] * vel[a][i];
            Ek += 0.5 * v2 * dx;
            double dp = (phi[a][i+1] - phi[a][i-1]) / (2.0*dx);
            Eg += 0.5 * dp * dp * dx;
            Em += 0.5 * m2 * phi[a][i] * phi[a][i] * dx;
            e_loc += 0.5*v2 + 0.5*dp*dp + 0.5*m2*phi[a][i]*phi[a][i];
            if (fabs(phi[a][i]) > peak_out[a]) peak_out[a] = fabs(phi[a][i]);
        }
        double P = phi[0][i] * phi[1][i] * phi[2][i];
        double P2 = P * P;
        double V = 0.5 * mu * P2 / (1.0 + kappa * P2);
        Ep += V * dx;
        e_loc += V;
        Eall += e_loc * dx;
        if (fabs(x) < core_r) Ecore += e_loc * dx;
    }

    *Ek_out = Ek; *Eg_out = Eg; *Em_out = Em; *Ep_out = Ep;
    *Et_out = Ek + Eg + Em + Ep;
    *Ecore_out = Ecore;
    *Etot_out = Eall;
}

int main(int argc, char **argv)
{
    parse_args(argc, argv);

    double dx = 2.0 * xmax / (Nx - 1);
    double dx2 = dx * dx;
    double m2 = mass * mass;

    /* CFL */
    double kmax = M_PI / dx;
    double dt = 0.8 * 2.0 / sqrt(kmax*kmax + m2);
    int Nt_equil = (int)(t_equil / dt) + 1;
    int Nt_run   = (int)(t_run / dt) + 1;
    (void)(Nt_equil + Nt_run); /* total steps, used only for reference */

    printf("rotating1d: delta=%.4f\n", delta);
    printf("  mu=%.1f kappa=%.1f mass=%.4f A=%.3f sigma=%.3f\n",
           mu, kappa, mass, A_init, sigma);
    printf("  Nx=%d xmax=%.1f dx=%.5f dt=%.6f\n", Nx, xmax, dx, dt);
    printf("  t_equil=%.0f (Nt=%d), t_run=%.0f (Nt=%d)\n",
           t_equil, Nt_equil, t_run, Nt_run);

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

    /* Initialize: Gaussians (all three identical) */
    for (int a = 0; a < 3; a++)
        for (int i = 0; i < Nx; i++) {
            double x = -xmax + i * dx;
            phi[a][i] = A_init * exp(-x * x / (2.0 * sigma * sigma));
        }

    /* Acceleration macro */
    #define COMPUTE_ACC() do { \
        for (int a = 0; a < 3; a++) { \
            acc[a][0] = acc[a][1] = acc[a][Nx-2] = acc[a][Nx-1] = 0; \
            for (int i = 1; i < Nx - 1; i++) { \
                double lapl = (phi[a][i+1] - 2.0*phi[a][i] + phi[a][i-1]) / dx2; \
                double fp = force_pot(phi[0][i], phi[1][i], phi[2][i], a); \
                acc[a][i] = lapl - m2*phi[a][i] + fp; \
            } \
        } \
    } while(0)

    COMPUTE_ACC();

    int ic = Nx / 2;  /* center index */

    /* ============================================================
     * PHASE 1: Equilibrate standard oscillon
     * ============================================================ */
    printf("\n--- Phase 1: Equilibration (t=0..%.0f) ---\n", t_equil);

    int print_every_eq = Nt_equil / 10;
    if (print_every_eq < 1) print_every_eq = 1;

    for (int n = 0; n < Nt_equil; n++) {
        if (n % print_every_eq == 0) {
            double Ek, Eg, Em, Ep, Et, Ec, Ea;
            double pk[3];
            compute_energies(phi, vel, Nx, dx, m2, xmax,
                             &Ek, &Eg, &Em, &Ep, &Et, &Ec, &Ea, pk);
            double fc = (Ea > 1e-20) ? Ec / Ea : 0.0;
            printf("  t=%7.1f  phi0=(%+.4f,%+.4f,%+.4f)  E=%.4f  fc=%.4f\n",
                   n*dt, phi[0][ic], phi[1][ic], phi[2][ic], Et, fc);
        }

        /* Velocity Verlet */
        for (int a = 0; a < 3; a++)
            for (int i = 1; i < Nx - 1; i++)
                vel[a][i] += 0.5 * dt * acc[a][i];
        for (int a = 0; a < 3; a++)
            for (int i = 1; i < Nx - 1; i++)
                phi[a][i] += dt * vel[a][i];
        COMPUTE_ACC();
        for (int a = 0; a < 3; a++)
            for (int i = 1; i < Nx - 1; i++)
                vel[a][i] += 0.5 * dt * acc[a][i];

        /* absorbing boundary */
        for (int a = 0; a < 3; a++)
            for (int i = 0; i < Nx; i++) {
                vel[a][i] *= damp[i];
                phi[a][i] *= damp[i];
            }
    }

    /* Record pre-perturbation state */
    double Ek0, Eg0, Em0, Ep0, Et0, Ec0, Ea0;
    double pk0[3];
    compute_energies(phi, vel, Nx, dx, m2, xmax,
                     &Ek0, &Eg0, &Em0, &Ep0, &Et0, &Ec0, &Ea0, pk0);
    double fc0 = (Ea0 > 1e-20) ? Ec0 / Ea0 : 0.0;

    printf("\n--- Post-equilibration state ---\n");
    printf("  E=%.6f  Ep=%.6f  fc=%.4f\n", Et0, Ep0, fc0);
    printf("  peaks: (%.4f, %.4f, %.4f)\n", pk0[0], pk0[1], pk0[2]);

    /* ============================================================
     * PHASE 2: Apply frequency perturbation
     * ============================================================ */
    if (delta > 0.0) {
        printf("\n--- Applying phase perturbation: delta=%.4f ---\n", delta);
        /* v₁ += δ·φ₁ (advance phase), v₃ -= δ·φ₃ (retard phase) */
        for (int i = 0; i < Nx; i++) {
            vel[0][i] += delta * phi[0][i];
            vel[2][i] -= delta * phi[2][i];
        }
        /* Recompute energy after perturbation */
        double Ekp, Egp, Emp, Epp, Etp, Ecp, Eap;
        double pkp[3];
        compute_energies(phi, vel, Nx, dx, m2, xmax,
                         &Ekp, &Egp, &Emp, &Epp, &Etp, &Ecp, &Eap, pkp);
        printf("  E after perturbation: %.6f (delta_E = %+.6f)\n",
               Etp, Etp - Et0);
    }

    COMPUTE_ACC();

    /* ============================================================
     * PHASE 3: Evolution with measurement
     * ============================================================ */
    printf("\n--- Phase 3: Evolution (t=%.0f..%.0f) ---\n", t_equil, t_equil + t_run);

    /* Time series file */
    char tspath[600];
    snprintf(tspath, sizeof(tspath), "%s/rotating_d%.4f_ts.tsv", outdir, delta);
    FILE *fts = fopen(tspath, "w");
    if (!fts) { fprintf(stderr, "Cannot open %s\n", tspath); return 1; }
    fprintf(fts, "time\tphi1_0\tphi2_0\tphi3_0\tpeak1\tpeak2\tpeak3\t"
                 "E_kin\tE_grad\tE_mass\tE_pot\tE_total\tf_core\n");

    /* DFT storage: per-field histories at center */
    int max_dft = 50000;
    double *phi1_hist = malloc(max_dft * sizeof(double));
    double *phi2_hist = malloc(max_dft * sizeof(double));
    double *phi3_hist = malloc(max_dft * sizeof(double));
    double *t_hist    = malloc(max_dft * sizeof(double));
    int n_dft = 0;

    int rec_every = Nt_run / 20000;
    if (rec_every < 1) rec_every = 1;
    int print_every = Nt_run / 20;
    if (print_every < 1) print_every = 1;
    int dft_every = Nt_run / max_dft;
    if (dft_every < 1) dft_every = 1;

    /* Energy tracking for dE/dt measurement */
    #define N_ETRACK 200
    double E_track[N_ETRACK], t_track[N_ETRACK];
    int n_etrack = 0;
    int etrack_every = Nt_run / N_ETRACK;
    if (etrack_every < 1) etrack_every = 1;

    for (int n = 0; n <= Nt_run; n++) {
        double t = t_equil + n * dt;

        /* DFT sampling */
        if (n % dft_every == 0 && n_dft < max_dft) {
            phi1_hist[n_dft] = phi[0][ic];
            phi2_hist[n_dft] = phi[1][ic];
            phi3_hist[n_dft] = phi[2][ic];
            t_hist[n_dft] = t;
            n_dft++;
        }

        int do_rec = (n % rec_every == 0);
        int do_print = (n % print_every == 0);
        int do_etrack = (n % etrack_every == 0 && n_etrack < N_ETRACK);

        if (do_rec || do_print || do_etrack) {
            double Ek, Eg, Em, Ep, Et, Ec, Ea;
            double pk[3];
            compute_energies(phi, vel, Nx, dx, m2, xmax,
                             &Ek, &Eg, &Em, &Ep, &Et, &Ec, &Ea, pk);
            double fc = (Ea > 1e-20) ? Ec / Ea : 0.0;

            if (do_etrack) {
                E_track[n_etrack] = Et;
                t_track[n_etrack] = t;
                n_etrack++;
            }

            if (do_rec)
                fprintf(fts, "%.6f\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t"
                             "%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6f\n",
                        t, phi[0][ic], phi[1][ic], phi[2][ic],
                        pk[0], pk[1], pk[2], Ek, Eg, Em, Ep, Et, fc);
            if (do_print)
                printf("  t=%7.1f  p0=(%+.4f,%+.4f,%+.4f)  pk=(%.3f,%.3f,%.3f)  "
                       "E=%.4f  fc=%.4f\n",
                       t, phi[0][ic], phi[1][ic], phi[2][ic],
                       pk[0], pk[1], pk[2], Et, fc);
        }

        if (n == Nt_run) break;

        /* Velocity Verlet */
        for (int a = 0; a < 3; a++)
            for (int i = 1; i < Nx - 1; i++)
                vel[a][i] += 0.5 * dt * acc[a][i];
        for (int a = 0; a < 3; a++)
            for (int i = 1; i < Nx - 1; i++)
                phi[a][i] += dt * vel[a][i];
        COMPUTE_ACC();
        for (int a = 0; a < 3; a++)
            for (int i = 1; i < Nx - 1; i++)
                vel[a][i] += 0.5 * dt * acc[a][i];

        /* absorbing boundary */
        for (int a = 0; a < 3; a++)
            for (int i = 0; i < Nx; i++) {
                vel[a][i] *= damp[i];
                phi[a][i] *= damp[i];
            }
    }

    fclose(fts);
    printf("  Time series: %s\n", tspath);

    /* ============================================================
     * PHASE 4: DFT analysis — per-field spectra
     * ============================================================ */
    printf("\n--- Spectral Analysis ---\n");

    /* Use second half for steady-state */
    int dft_start = n_dft / 2;
    int dft_len = n_dft - dft_start;

    if (dft_len > 100) {
        char dftpath[600];
        snprintf(dftpath, sizeof(dftpath), "%s/rotating_d%.4f_spectrum.tsv", outdir, delta);
        FILE *fdft = fopen(dftpath, "w");
        fprintf(fdft, "omega\tpower1\tpower2\tpower3\n");

        double T = t_hist[n_dft-1] - t_hist[dft_start];
        int nf = 1000;
        double peak_pow[3] = {0}, peak_om[3] = {0};

        for (int k = 0; k < nf; k++) {
            double omega = 3.0 * mass * k / nf;
            double re[3] = {0}, im[3] = {0};
            for (int j = dft_start; j < n_dft; j++) {
                double dtj = (j > dft_start) ?
                    (t_hist[j]-t_hist[j-1]) : (t_hist[dft_start+1]-t_hist[dft_start]);
                double c = cos(omega * t_hist[j]) * dtj;
                double s = sin(omega * t_hist[j]) * dtj;
                re[0] += phi1_hist[j] * c; im[0] += phi1_hist[j] * s;
                re[1] += phi2_hist[j] * c; im[1] += phi2_hist[j] * s;
                re[2] += phi3_hist[j] * c; im[2] += phi3_hist[j] * s;
            }
            double pw[3];
            for (int a = 0; a < 3; a++) {
                pw[a] = (re[a]*re[a] + im[a]*im[a]) / (T*T);
                if (pw[a] > peak_pow[a]) { peak_pow[a] = pw[a]; peak_om[a] = omega; }
            }
            fprintf(fdft, "%.6f\t%.6e\t%.6e\t%.6e\n", omega, pw[0], pw[1], pw[2]);
        }
        fclose(fdft);

        printf("  Per-field peak frequencies:\n");
        printf("    phi1: omega=%.4f (power=%.4e)\n", peak_om[0], peak_pow[0]);
        printf("    phi2: omega=%.4f (power=%.4e)\n", peak_om[1], peak_pow[1]);
        printf("    phi3: omega=%.4f (power=%.4e)\n", peak_om[2], peak_pow[2]);
        printf("  Frequency splitting: omega1-omega3 = %.4f\n",
               peak_om[0] - peak_om[2]);
        printf("  Spectrum: %s\n", dftpath);
    }

    /* ============================================================
     * PHASE 5: Energy loss rate analysis
     * ============================================================ */
    printf("\n--- Energy Loss Rate ---\n");

    /* Fit dE/dt using linear regression on second half of energy track */
    int eh = n_etrack / 2;
    if (eh > 5) {
        /* Early period: first quarter */
        int q1s = 0, q1e = n_etrack / 4;
        double sum_t = 0, sum_E = 0, sum_tE = 0, sum_t2 = 0;
        for (int j = q1s; j < q1e; j++) {
            sum_t += t_track[j]; sum_E += E_track[j];
            sum_tE += t_track[j]*E_track[j]; sum_t2 += t_track[j]*t_track[j];
        }
        int nq = q1e - q1s;
        double dEdt_early = (nq*sum_tE - sum_t*sum_E) / (nq*sum_t2 - sum_t*sum_t);

        /* Late period: last quarter */
        int q3s = 3 * n_etrack / 4, q3e = n_etrack;
        sum_t = sum_E = sum_tE = sum_t2 = 0;
        for (int j = q3s; j < q3e; j++) {
            sum_t += t_track[j]; sum_E += E_track[j];
            sum_tE += t_track[j]*E_track[j]; sum_t2 += t_track[j]*t_track[j];
        }
        nq = q3e - q3s;
        double dEdt_late = (nq*sum_tE - sum_t*sum_E) / (nq*sum_t2 - sum_t*sum_t);

        /* Overall */
        sum_t = sum_E = sum_tE = sum_t2 = 0;
        for (int j = eh; j < n_etrack; j++) {
            sum_t += t_track[j]; sum_E += E_track[j];
            sum_tE += t_track[j]*E_track[j]; sum_t2 += t_track[j]*t_track[j];
        }
        nq = n_etrack - eh;
        double dEdt_overall = (nq*sum_tE - sum_t*sum_E) / (nq*sum_t2 - sum_t*sum_t);

        printf("  E(start_run) = %.6f\n", E_track[0]);
        printf("  E(end_run)   = %.6f\n", E_track[n_etrack-1]);
        printf("  dE/dt (early quarter):  %.6e\n", dEdt_early);
        printf("  dE/dt (late quarter):   %.6e\n", dEdt_late);
        printf("  dE/dt (second half):    %.6e\n", dEdt_overall);

        /* Write to summary line (append mode) */
        char sumpath[600];
        snprintf(sumpath, sizeof(sumpath), "%s/rotating_summary.tsv", outdir);
        int first_write = 0;
        {
            FILE *ft = fopen(sumpath, "r");
            if (!ft) first_write = 1;
            else fclose(ft);
        }
        FILE *fsum = fopen(sumpath, "a");
        if (first_write)
            fprintf(fsum, "delta\tE_start\tE_end\tdEdt_early\tdEdt_late\t"
                          "dEdt_overall\tomega1\tomega2\tomega3\tdom13\n");
        double om1 = 0, om2 = 0, om3 = 0;
        if (dft_len > 100) {
            /* Re-read from above (captured in peak_om) — need to recompute */
            /* Actually just use the printed values */
        }
        /* Recompute peak frequencies inline */
        if (dft_len > 100) {
            double T = t_hist[n_dft-1] - t_hist[dft_start];
            int nf = 1000;
            double ppow[3] = {0}, pom[3] = {0};
            for (int k = 0; k < nf; k++) {
                double omega = 3.0 * mass * k / nf;
                double re[3] = {0}, im[3] = {0};
                for (int j = dft_start; j < n_dft; j++) {
                    double dtj = (j > dft_start) ?
                        (t_hist[j]-t_hist[j-1]) : (t_hist[dft_start+1]-t_hist[dft_start]);
                    double c = cos(omega * t_hist[j]) * dtj;
                    double s = sin(omega * t_hist[j]) * dtj;
                    re[0] += phi1_hist[j]*c; im[0] += phi1_hist[j]*s;
                    re[1] += phi2_hist[j]*c; im[1] += phi2_hist[j]*s;
                    re[2] += phi3_hist[j]*c; im[2] += phi3_hist[j]*s;
                }
                for (int a = 0; a < 3; a++) {
                    double pw = (re[a]*re[a]+im[a]*im[a])/(T*T);
                    if (pw > ppow[a]) { ppow[a] = pw; pom[a] = omega; }
                }
            }
            om1 = pom[0]; om2 = pom[1]; om3 = pom[2];
        }
        fprintf(fsum, "%.4f\t%.6f\t%.6f\t%.6e\t%.6e\t%.6e\t%.4f\t%.4f\t%.4f\t%.4f\n",
                delta, E_track[0], E_track[n_etrack-1],
                dEdt_early, dEdt_late, dEdt_overall,
                om1, om2, om3, om1 - om3);
        fclose(fsum);
    }

    printf("\nDone. delta=%.4f\n", delta);

    for (int a = 0; a < 3; a++) { free(phi[a]); free(vel[a]); free(acc[a]); }
    free(damp); free(phi1_hist); free(phi2_hist); free(phi3_hist); free(t_hist);
    return 0;
}
