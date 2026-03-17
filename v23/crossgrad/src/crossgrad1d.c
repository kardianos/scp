/*
 * crossgrad1d.c — 1D cross-gradient oscillon dynamics (V23-B Phase 1)
 *
 * Three massive scalars with triple-product coupling plus cross-gradient term.
 * In 1D (x-direction), field 1 gets Laplacian coefficient (1+eta), fields 2,3 get 1.
 *
 * Modified equations of motion:
 *   d^2 phi_1/dt^2 = (1+eta) d^2 phi_1/dx^2 - m^2 phi_1 - dV/dphi_1
 *   d^2 phi_2/dt^2 =         d^2 phi_2/dx^2 - m^2 phi_2 - dV/dphi_2
 *   d^2 phi_3/dt^2 =         d^2 phi_3/dx^2 - m^2 phi_3 - dV/dphi_3
 *
 * where V = (mu/2) P^2 / (1 + kappa P^2), P = phi_1 phi_2 phi_3.
 *
 * Scans eta = 0.0, 0.2, 0.5, 1.0, 1.5, 2.0 automatically.
 * For each eta: equilibrate for t_equil, then add antisymmetric perturbation
 * and evolve for t_pert more to measure stability of shear mode.
 *
 * Based on v21/src/triad1d.c.
 *
 * Compile: gcc -O3 -Wall -o crossgrad1d src/crossgrad1d.c -lm
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Parameters */
static double mu     = -20.0;
static double kappa  = 20.0;
static double mass   = 1.0;
static double A_init = 0.8;
static double sigma  = 3.0;
static int    Nx     = 2000;
static double xmax   = 60.0;
static double tfinal_equil = 5000.0;
static double tfinal_pert  = 1000.0;
static double eta_single = -1.0;  /* if >=0, run single eta instead of scan */
static double pert_eps = 1e-3;    /* antisymmetric perturbation amplitude */
static double pert_sigma = 3.0;   /* perturbation Gaussian width */
static char   outdir[512] = "v23/crossgrad/data";

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
        else if (!strcmp(argv[i], "-tequil")) tfinal_equil = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-tpert"))  tfinal_pert  = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-eta"))    eta_single   = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-eps"))    pert_eps     = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-o"))      strncpy(outdir, argv[i+1], sizeof(outdir)-1);
        else { fprintf(stderr, "Unknown flag: %s\n", argv[i]); exit(1); }
    }
}

/* -dV/dphi_a where V = (mu/2)P^2/(1+kappa P^2), P = phi1 phi2 phi3 */
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
    return -mu * P * dP / denom2;
}

/* DFT: find peak frequency in second half of time series */
static double find_peak_omega(double *phi_hist, double *t_hist, int n_pts,
                              double *out_power)
{
    int start = n_pts / 2;
    if (n_pts - start < 100) { *out_power = 0; return 0; }
    double T = t_hist[n_pts-1] - t_hist[start];
    int nf = 500;
    double peak_pow = 0, peak_om = 0;
    for (int k = 1; k < nf; k++) {
        double omega = 3.0 * mass * k / nf;
        double re = 0, im = 0;
        for (int j = start; j < n_pts; j++) {
            double dtj = (j > start) ?
                (t_hist[j] - t_hist[j-1]) : (t_hist[start+1] - t_hist[start]);
            re += phi_hist[j] * cos(omega * t_hist[j]) * dtj;
            im += phi_hist[j] * sin(omega * t_hist[j]) * dtj;
        }
        double pw = (re*re + im*im) / (T*T);
        if (pw > peak_pow) { peak_pow = pw; peak_om = omega; }
    }
    *out_power = peak_pow;
    return peak_om;
}

/* Result structure for one eta run */
typedef struct {
    double eta;
    int    survived;
    double omega;
    double peak_amp;
    double gap_margin;
    double energy;
    double pert_growth_rate;
    double pert_final_amp;
    int    pert_grew;
    double pert_extent;
} Result;

static Result run_one_eta(double eta)
{
    Result res;
    memset(&res, 0, sizeof(res));
    res.eta = eta;

    double dx = 2.0 * xmax / (Nx - 1);
    double dx2 = dx * dx;
    double m2 = mass * mass;

    /* CFL: max wave speed is sqrt(1+eta) for field 1 */
    double c_max = sqrt(1.0 + fmax(eta, 0.0));
    double kmax = M_PI / dx;
    double dt = 0.8 * 2.0 / sqrt(c_max * c_max * kmax * kmax + m2);
    int Nt_equil = (int)(tfinal_equil / dt) + 1;
    int Nt_pert  = (int)(tfinal_pert / dt) + 1;

    printf("\n=== eta = %.2f ===\n", eta);
    printf("  dx=%.5f dt=%.6f c_max=%.3f Nt_equil=%d Nt_pert=%d\n",
           dx, dt, c_max, Nt_equil, Nt_pert);

    /* Allocate fields */
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

    /* Initialize: symmetric Gaussian oscillon (all three fields equal) */
    for (int a = 0; a < 3; a++)
        for (int i = 0; i < Nx; i++) {
            double x = -xmax + i * dx;
            phi[a][i] = A_init * exp(-x * x / (2.0 * sigma * sigma));
        }

    /* Laplacian coefficients: field 0 gets (1+eta), fields 1,2 get 1 */
    double lap_coeff[3] = {1.0 + eta, 1.0, 1.0};

    /* Compute acceleration macro */
    #define COMPUTE_ACC_ETA() do { \
        for (int a = 0; a < 3; a++) { \
            acc[a][0] = acc[a][1] = acc[a][Nx-2] = acc[a][Nx-1] = 0; \
            for (int i = 1; i < Nx - 1; i++) { \
                double lapl = (phi[a][i+1] - 2.0*phi[a][i] + phi[a][i-1]) / dx2; \
                double fp = force_pot(phi[0][i], phi[1][i], phi[2][i], a); \
                acc[a][i] = lap_coeff[a] * lapl - m2*phi[a][i] + fp; \
            } \
        } \
    } while(0)

    COMPUTE_ACC_ETA();

    /* Output file for time series */
    char tspath[600];
    snprintf(tspath, sizeof(tspath), "%s/crossgrad_eta%.1f_ts.tsv", outdir, eta);
    FILE *fts = fopen(tspath, "w");
    if (!fts) { fprintf(stderr, "Cannot open %s\n", tspath); exit(1); }
    fprintf(fts, "time\tphi1_0\tphi2_0\tphi3_0\tpeak1\tpeak2\tpeak3\t"
                 "E_kin\tE_grad\tE_mass\tE_pot\tE_total\tf_core\t"
                 "delta_phi\tphase\n");

    /* DFT storage */
    int max_dft = 50000;
    double *phi0_hist = malloc(max_dft * sizeof(double));
    double *t_hist = malloc(max_dft * sizeof(double));
    int n_dft = 0;

    /* Antisymmetric mode history (for perturbation phase) */
    int max_pert_hist = 20000;
    double *asym_hist = malloc(max_pert_hist * sizeof(double));
    double *asym_t = malloc(max_pert_hist * sizeof(double));
    int n_asym = 0;

    int rec_every_equil = Nt_equil / 20000;
    if (rec_every_equil < 1) rec_every_equil = 1;
    int print_every_equil = Nt_equil / 20;
    if (print_every_equil < 1) print_every_equil = 1;
    int dft_every = Nt_equil / max_dft;
    if (dft_every < 1) dft_every = 1;

    double core_r = 3.0 * sigma;
    int ic = Nx / 2;

    /* Track late-time peak amplitude for survival check (last 20% of equil) */
    double late_peak_max = 0;
    double late_energy = 0;
    double late_fc_min = 1.0;
    int late_start = (int)(0.8 * Nt_equil);

    /* ============ PHASE 1: Equilibration ============ */
    printf("  Phase 1: Equilibration (t=0..%.0f)\n", tfinal_equil);

    for (int n = 0; n <= Nt_equil; n++) {
        double t = n * dt;

        /* Record DFT history */
        if (n % dft_every == 0 && n_dft < max_dft) {
            phi0_hist[n_dft] = phi[0][ic];
            t_hist[n_dft] = t;
            n_dft++;
        }

        int do_rec = (n % rec_every_equil == 0);
        int do_print = (n % print_every_equil == 0);

        if (do_rec || do_print) {
            double Ek = 0, Eg = 0, Em = 0, Ep = 0;
            double Ecore = 0, Eall = 0;
            double peak[3] = {0};

            for (int i = 1; i < Nx - 1; i++) {
                double x = -xmax + i * dx;
                for (int a = 0; a < 3; a++) {
                    Ek += 0.5 * vel[a][i] * vel[a][i] * dx;
                    double dp = (phi[a][i+1] - phi[a][i-1]) / (2.0*dx);
                    /* Gradient energy includes cross-gradient for field 0 */
                    double g_coeff = (a == 0) ? (1.0 + eta) : 1.0;
                    Eg += 0.5 * g_coeff * dp * dp * dx;
                    Em += 0.5 * m2 * phi[a][i] * phi[a][i] * dx;
                    if (fabs(phi[a][i]) > peak[a]) peak[a] = fabs(phi[a][i]);
                }

                double P = phi[0][i] * phi[1][i] * phi[2][i];
                double P2 = P * P;
                double V = 0.5 * mu * P2 / (1.0 + kappa * P2);
                Ep += V * dx;

                double e = V;
                for (int a = 0; a < 3; a++) {
                    e += 0.5*vel[a][i]*vel[a][i];
                    double dp = (phi[a][i+1] - phi[a][i-1]) / (2.0*dx);
                    double g_coeff = (a == 0) ? (1.0 + eta) : 1.0;
                    e += 0.5*g_coeff*dp*dp + 0.5*m2*phi[a][i]*phi[a][i];
                }
                Eall += e * dx;
                if (fabs(x) < core_r) Ecore += e * dx;
            }

            double Et = Ek + Eg + Em + Ep;
            double fc = (Eall > 1e-20) ? Ecore / Eall : 0.0;
            double delta_phi = phi[0][ic] - phi[1][ic];

            /* Track late-time peak for survival */
            if (n >= late_start) {
                if (peak[0] > late_peak_max) late_peak_max = peak[0];
                late_energy = Et;
                if (fc < late_fc_min) late_fc_min = fc;
            }

            if (do_rec)
                fprintf(fts, "%.6f\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t"
                             "%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6f\t%.6e\t%s\n",
                        t, phi[0][ic], phi[1][ic], phi[2][ic],
                        peak[0], peak[1], peak[2], Ek, Eg, Em, Ep, Et, fc,
                        delta_phi, "equil");
            if (do_print)
                printf("  t=%7.1f  p=(%+.4f,%+.4f,%+.4f)  pk=(%.3f,%.3f,%.3f)  "
                       "E=%+.4f  fc=%.3f\n",
                       t, phi[0][ic], phi[1][ic], phi[2][ic],
                       peak[0], peak[1], peak[2], Et, fc);
        }

        if (n == Nt_equil) break;

        /* Velocity Verlet */
        for (int a = 0; a < 3; a++)
            for (int i = 1; i < Nx - 1; i++)
                vel[a][i] += 0.5 * dt * acc[a][i];
        for (int a = 0; a < 3; a++)
            for (int i = 1; i < Nx - 1; i++)
                phi[a][i] += dt * vel[a][i];
        COMPUTE_ACC_ETA();
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

    /* Survival from late-time tracking */
    res.peak_amp = late_peak_max;
    res.energy = late_energy;
    res.survived = (late_fc_min > 0.3 && late_peak_max > 0.05) ? 1 : 0;

    /* DFT analysis */
    double dft_power;
    res.omega = find_peak_omega(phi0_hist, t_hist, n_dft, &dft_power);
    res.gap_margin = (mass - res.omega) / mass;

    printf("  Equilibrium: omega=%.4f gap_margin=%.4f peak=%.4f fc_min=%.3f survived=%d E=%.4f\n",
           res.omega, res.gap_margin, res.peak_amp, late_fc_min, res.survived, res.energy);

    /* ============ PHASE 2: Perturbation ============ */
    if (res.survived) {
        printf("  Phase 2: Antisymmetric perturbation (eps=%.1e, t=%.0f..%.0f)\n",
               pert_eps, tfinal_equil, tfinal_equil + tfinal_pert);

        /* Add antisymmetric perturbation: delta_phi1 = +eps*G, delta_phi2 = -eps*G */
        for (int i = 0; i < Nx; i++) {
            double x = -xmax + i * dx;
            double g = exp(-x * x / (2.0 * pert_sigma * pert_sigma));
            phi[0][i] += pert_eps * g;
            phi[1][i] -= pert_eps * g;
            /* phi[2] unchanged */
        }

        COMPUTE_ACC_ETA();

        int rec_every_pert = Nt_pert / 10000;
        if (rec_every_pert < 1) rec_every_pert = 1;
        int print_every_pert = Nt_pert / 10;
        if (print_every_pert < 1) print_every_pert = 1;
        int asym_every = Nt_pert / max_pert_hist;
        if (asym_every < 1) asym_every = 1;

        double asym_initial = 0;

        for (int n = 0; n <= Nt_pert; n++) {
            double t = tfinal_equil + n * dt;

            /* Record antisymmetric amplitude at center */
            if (n % asym_every == 0 && n_asym < max_pert_hist) {
                double delta = phi[0][ic] - phi[1][ic];
                asym_hist[n_asym] = delta;
                asym_t[n_asym] = n * dt;  /* time from perturbation start */
                if (n_asym == 0) asym_initial = fabs(delta);
                n_asym++;
            }

            int do_rec = (n % rec_every_pert == 0);
            int do_print = (n % print_every_pert == 0);

            if (do_rec || do_print) {
                double Ek = 0, Eg = 0, Em = 0, Ep = 0;
                double Ecore = 0, Eall = 0;
                double peak[3] = {0};

                for (int i = 1; i < Nx - 1; i++) {
                    double x = -xmax + i * dx;
                    for (int a = 0; a < 3; a++) {
                        Ek += 0.5 * vel[a][i] * vel[a][i] * dx;
                        double dp = (phi[a][i+1] - phi[a][i-1]) / (2.0*dx);
                        double g_coeff = (a == 0) ? (1.0 + eta) : 1.0;
                        Eg += 0.5 * g_coeff * dp * dp * dx;
                        Em += 0.5 * m2 * phi[a][i] * phi[a][i] * dx;
                        if (fabs(phi[a][i]) > peak[a]) peak[a] = fabs(phi[a][i]);
                    }

                    double P = phi[0][i] * phi[1][i] * phi[2][i];
                    double P2 = P * P;
                    double V = 0.5 * mu * P2 / (1.0 + kappa * P2);
                    Ep += V * dx;

                    double e = V;
                    for (int a = 0; a < 3; a++) {
                        e += 0.5*vel[a][i]*vel[a][i];
                        double dp = (phi[a][i+1] - phi[a][i-1]) / (2.0*dx);
                        double g_coeff = (a == 0) ? (1.0 + eta) : 1.0;
                        e += 0.5*g_coeff*dp*dp + 0.5*m2*phi[a][i]*phi[a][i];
                    }
                    Eall += e * dx;
                    if (fabs(x) < core_r) Ecore += e * dx;
                }

                double Et = Ek + Eg + Em + Ep;
                double fc = (Eall > 1e-20) ? Ecore / Eall : 0.0;
                double delta_phi = phi[0][ic] - phi[1][ic];

                if (do_rec)
                    fprintf(fts, "%.6f\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t"
                                 "%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6f\t%.6e\t%s\n",
                            t, phi[0][ic], phi[1][ic], phi[2][ic],
                            peak[0], peak[1], peak[2], Ek, Eg, Em, Ep, Et, fc,
                            delta_phi, "pert");
                if (do_print)
                    printf("  t=%7.1f  delta=%.3e  pk=(%.3f,%.3f,%.3f)  E=%+.4f  fc=%.3f\n",
                           t, delta_phi, peak[0], peak[1], peak[2], Et, fc);
            }

            if (n == Nt_pert) break;

            /* Velocity Verlet */
            for (int a = 0; a < 3; a++)
                for (int i = 1; i < Nx - 1; i++)
                    vel[a][i] += 0.5 * dt * acc[a][i];
            for (int a = 0; a < 3; a++)
                for (int i = 1; i < Nx - 1; i++)
                    phi[a][i] += dt * vel[a][i];
            COMPUTE_ACC_ETA();
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

        /* Analyze antisymmetric perturbation evolution */
        if (n_asym > 10) {
            /* Find envelope: max |delta| in sliding windows */
            double final_amp = 0;
            for (int j = n_asym - n_asym/10; j < n_asym; j++)
                if (fabs(asym_hist[j]) > final_amp) final_amp = fabs(asym_hist[j]);
            res.pert_final_amp = final_amp;
            res.pert_grew = (final_amp > 1.5 * asym_initial) ? 1 : 0;

            /* Growth/decay rate: fit log(envelope) vs t */
            /* Use peak envelope in windows */
            int n_win = 20;
            int win_sz = n_asym / n_win;
            if (win_sz < 2) win_sz = 2;
            double sum_t = 0, sum_lnA = 0, sum_t2 = 0, sum_tlnA = 0;
            int n_fit = 0;
            for (int w = 0; w < n_win; w++) {
                int j0 = w * win_sz;
                int j1 = j0 + win_sz;
                if (j1 > n_asym) j1 = n_asym;
                double env = 0;
                double t_mid = 0;
                for (int j = j0; j < j1; j++) {
                    if (fabs(asym_hist[j]) > env) env = fabs(asym_hist[j]);
                    t_mid += asym_t[j];
                }
                t_mid /= (j1 - j0);
                if (env > 1e-20) {
                    double lnA = log(env);
                    sum_t += t_mid;
                    sum_lnA += lnA;
                    sum_t2 += t_mid * t_mid;
                    sum_tlnA += t_mid * lnA;
                    n_fit++;
                }
            }
            if (n_fit > 2) {
                double denom = n_fit * sum_t2 - sum_t * sum_t;
                if (fabs(denom) > 1e-30)
                    res.pert_growth_rate = (n_fit * sum_tlnA - sum_t * sum_lnA) / denom;
            }

            /* Spatial extent: measure RMS width of (phi1-phi2) at final time */
            double sum_d2x2 = 0, sum_d2 = 0;
            for (int i = 1; i < Nx - 1; i++) {
                double x = -xmax + i * dx;
                double d = phi[0][i] - phi[1][i];
                sum_d2 += d * d * dx;
                sum_d2x2 += d * d * x * x * dx;
            }
            res.pert_extent = (sum_d2 > 1e-30) ? sqrt(sum_d2x2 / sum_d2) : 0;

            printf("  Perturbation: initial=%.3e final=%.3e rate=%.4e grew=%d extent=%.2f\n",
                   asym_initial, final_amp, res.pert_growth_rate, res.pert_grew,
                   res.pert_extent);
        }
    } else {
        printf("  Oscillon did not survive — skipping perturbation phase.\n");
        res.pert_growth_rate = 0;
        res.pert_final_amp = 0;
        res.pert_grew = 0;
        res.pert_extent = 0;
    }

    fclose(fts);
    printf("  Output: %s\n", tspath);

    /* Cleanup */
    for (int a = 0; a < 3; a++) { free(phi[a]); free(vel[a]); free(acc[a]); }
    free(damp); free(phi0_hist); free(t_hist);
    free(asym_hist); free(asym_t);

    #undef COMPUTE_ACC_ETA
    return res;
}

int main(int argc, char **argv)
{
    parse_args(argc, argv);

    printf("crossgrad1d — V23-B Phase 1: Cross-Gradient Oscillon Dynamics\n");
    printf("  mu=%.3f kappa=%.4f mass=%.4f A=%.3f sigma=%.3f\n",
           mu, kappa, mass, A_init, sigma);
    printf("  Nx=%d xmax=%.1f tequil=%.0f tpert=%.0f eps=%.1e\n",
           Nx, xmax, tfinal_equil, tfinal_pert, pert_eps);

    /* Eta values to scan */
    double eta_vals[] = {0.0, 0.2, 0.5, 1.0, 1.5, 2.0};
    int n_eta = sizeof(eta_vals) / sizeof(eta_vals[0]);

    /* If single eta requested, override */
    if (eta_single >= 0.0) {
        eta_vals[0] = eta_single;
        n_eta = 1;
    }

    Result *results = malloc(n_eta * sizeof(Result));

    for (int ie = 0; ie < n_eta; ie++) {
        results[ie] = run_one_eta(eta_vals[ie]);
    }

    /* Print summary table */
    printf("\n");
    printf("============================================================\n");
    printf("  SUMMARY TABLE: Cross-Gradient Oscillon Scan\n");
    printf("============================================================\n");
    printf("  eta    surv  omega   gap%%   peak    E_total  rate      grew  extent\n");
    printf("  -----  ----  ------  -----  ------  -------  --------  ----  ------\n");
    for (int ie = 0; ie < n_eta; ie++) {
        Result *r = &results[ie];
        printf("  %5.2f  %4s  %6.4f  %5.2f  %6.4f  %7.3f  %+8.5f  %4s  %6.2f\n",
               r->eta,
               r->survived ? "YES" : "NO",
               r->omega,
               100.0 * r->gap_margin,
               r->peak_amp,
               r->energy,
               r->pert_growth_rate,
               r->pert_grew ? "YES" : "NO",
               r->pert_extent);
    }
    printf("============================================================\n\n");

    /* Write RESULTS.md */
    char respath[600];
    snprintf(respath, sizeof(respath), "v23/crossgrad/RESULTS.md");
    FILE *fres = fopen(respath, "w");
    if (fres) {
        fprintf(fres, "# V23-B Cross-Gradient Oscillon Results\n\n");
        fprintf(fres, "## Parameters\n\n");
        fprintf(fres, "- mu = %.1f, kappa = %.1f, m = %.1f\n", mu, kappa, mass);
        fprintf(fres, "- A_init = %.1f, sigma = %.1f\n", A_init, sigma);
        fprintf(fres, "- Nx = %d, xmax = %.1f\n", Nx, xmax);
        fprintf(fres, "- t_equil = %.0f, t_pert = %.0f\n", tfinal_equil, tfinal_pert);
        fprintf(fres, "- Perturbation: eps = %.1e, sigma_pert = %.1f\n\n",
                pert_eps, pert_sigma);

        fprintf(fres, "## Summary Table\n\n");
        fprintf(fres, "| eta  | survived | omega  | gap %%  | peak  | E_total | growth rate | grew | extent |\n");
        fprintf(fres, "|------|----------|--------|--------|-------|---------|-------------|------|--------|\n");
        for (int ie = 0; ie < n_eta; ie++) {
            Result *r = &results[ie];
            fprintf(fres, "| %.2f | %-8s | %.4f | %5.2f | %.4f | %7.3f | %+.5f | %-4s | %.2f |\n",
                    r->eta,
                    r->survived ? "YES" : "NO",
                    r->omega,
                    100.0 * r->gap_margin,
                    r->peak_amp,
                    r->energy,
                    r->pert_growth_rate,
                    r->pert_grew ? "YES" : "NO",
                    r->pert_extent);
        }

        fprintf(fres, "\n## Analysis\n\n");
        fprintf(fres, "### Oscillon Survival\n\n");
        int max_survived_idx = -1;
        for (int ie = n_eta - 1; ie >= 0; ie--) {
            if (results[ie].survived) { max_survived_idx = ie; break; }
        }
        if (max_survived_idx >= 0) {
            fprintf(fres, "Oscillon survives up to eta = %.2f.\n\n",
                    results[max_survived_idx].eta);
        } else {
            fprintf(fres, "Oscillon does not survive at any tested eta value.\n\n");
        }

        fprintf(fres, "### Breathing Frequency\n\n");
        fprintf(fres, "The symmetric oscillon breathing frequency omega should increase with eta\n");
        fprintf(fres, "since the effective wave speed for the symmetric mode is sqrt(1+eta/3).\n\n");

        fprintf(fres, "### Antisymmetric Perturbation\n\n");
        fprintf(fres, "The antisymmetric mode (delta_phi1 = +eps, delta_phi2 = -eps) breaks the\n");
        fprintf(fres, "phi_1 <-> phi_2 symmetry. With eta > 0, field 1 has faster wave propagation\n");
        fprintf(fres, "than field 2, so the shear perturbation experiences a different effective\n");
        fprintf(fres, "potential than the compression mode.\n\n");
        fprintf(fres, "Growth rate > 0 means the antisymmetric mode is unstable (shear instability).\n");
        fprintf(fres, "Growth rate < 0 means the perturbation decays (stable oscillon).\n\n");

        fprintf(fres, "## Data Files\n\n");
        for (int ie = 0; ie < n_eta; ie++) {
            fprintf(fres, "- `data/crossgrad_eta%.1f_ts.tsv` — time series for eta=%.2f\n",
                    results[ie].eta, results[ie].eta);
        }
        fprintf(fres, "\n");

        fclose(fres);
        printf("Results written to %s\n", respath);
    }

    free(results);
    return 0;
}
