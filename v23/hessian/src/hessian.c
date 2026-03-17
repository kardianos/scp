/*
 * hessian.c — Hessian splitting analysis for the symmetric triad oscillon
 *
 * Computes the eigenvalues of the 3x3 Hessian matrix ∂²V/∂φ_a∂φ_b
 * evaluated on the symmetric background φ₁=φ₂=φ₃=f(r), where f(r)
 * is the radial profile of the v21 1D oscillon.
 *
 * V(P) = (μ/2)P²/(1+κP²),  P = φ₁φ₂φ₃
 *
 * Hessian eigenvalues on symmetric background:
 *   λ_sym(r)  = μf⁴(5 - 7κf⁶)/(1+κf⁶)³     (compression mode, (1,1,1)/√3)
 *   λ_anti(r) = -μf⁴/(1+κf⁶)²                (shear modes, 2-fold degenerate)
 *
 * Effective masses:
 *   m²_sym(r)  = m² + λ_sym(r)
 *   m²_anti(r) = m² + λ_anti(r)
 *
 * The oscillon oscillates as f(r)cos(ωt), so time-averaged quantities
 * involve ⟨f⁴cos⁴(ωt)⟩ = (3/8)f⁴_peak.
 *
 * Method: evolve 1D symmetric triad to equilibrium, extract envelope,
 * compute Hessian eigenvalues.
 *
 * Compile: gcc -O3 -Wall -o hessian src/hessian.c -lm
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* --- Parameters (v21 production) --- */
static double mu    = -20.0;
static double kappa = 20.0;
static double mass  = 1.0;
static double A_init = 0.8;
static double sigma_init = 3.0;

/* Grid */
static int    Nx    = 4000;
static double xmax  = 80.0;
static double tfinal = 5000.0;

static char outdir[512] = "v23/hessian/data";

/* --- Potential force: -dV/dφ_a for symmetric mode φ₁=φ₂=φ₃=φ --- */
/* V = (μ/2)P²/(1+κP²), P = φ³ */
/* dV/dφ = dV/dP · dP/dφ = μP/(1+κP²)² · 3φ² = 3μφ⁵/(1+κφ⁶)² · φ² */
/* Actually: P = φ³, dP/dφ = 3φ² (for symmetric mode where all three are the same) */
/* But we evolve the FULL 3-field system to be safe. */

/* -dV/dphi_a where V = (mu/2)P^2/(1+kappa P^2), P = phi1 phi2 phi3 */
static inline double force_pot_sym(double phi)
{
    /* For symmetric mode: phi1=phi2=phi3=phi, P=phi^3 */
    double P  = phi * phi * phi;
    double P2 = P * P;  /* phi^6 */
    double denom2 = (1.0 + kappa * P2);
    denom2 = denom2 * denom2;
    double dP = phi * phi;  /* dP/dphi_a = phi_b * phi_c = phi^2 */
    return -mu * P * dP / denom2;
}

int main(int argc, char **argv)
{
    /* Parse command-line args */
    for (int i = 1; i < argc - 1; i += 2) {
        if      (!strcmp(argv[i], "-mu"))     mu     = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-kappa"))  kappa  = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-mass"))   mass   = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-A"))      A_init = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-sigma"))  sigma_init = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-Nx"))     Nx     = atoi(argv[i+1]);
        else if (!strcmp(argv[i], "-xmax"))   xmax   = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-tfinal")) tfinal = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-o"))      strncpy(outdir, argv[i+1], sizeof(outdir)-1);
    }

    double dx  = 2.0 * xmax / (Nx - 1);
    double dx2 = dx * dx;
    double m2  = mass * mass;

    /* CFL condition */
    double kmax = M_PI / dx;
    double dt = 0.8 * 2.0 / sqrt(kmax * kmax + m2);
    int Nt = (int)(tfinal / dt) + 1;

    printf("Hessian splitting analysis\n");
    printf("  mu=%.1f kappa=%.1f mass=%.3f A=%.3f sigma=%.3f\n",
           mu, kappa, mass, A_init, sigma_init);
    printf("  Nx=%d xmax=%.1f dx=%.5f dt=%.6f Nt=%d tfinal=%.0f\n",
           Nx, xmax, dx, dt, Nt, tfinal);

    /* Allocate: evolve ONE field (symmetric mode) */
    double *phi = calloc(Nx, sizeof(double));
    double *vel = calloc(Nx, sizeof(double));
    double *acc = calloc(Nx, sizeof(double));

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

    /* Initialize: Gaussian */
    for (int i = 0; i < Nx; i++) {
        double x = -xmax + i * dx;
        phi[i] = A_init * exp(-x * x / (2.0 * sigma_init * sigma_init));
    }

    /* Compute acceleration macro */
    #define COMPUTE_ACC() do { \
        acc[0] = acc[Nx-1] = 0; \
        for (int i = 1; i < Nx - 1; i++) { \
            double lapl = (phi[i+1] - 2.0*phi[i] + phi[i-1]) / dx2; \
            double fp = force_pot_sym(phi[i]); \
            acc[i] = lapl - m2*phi[i] + fp; \
        } \
    } while(0)

    COMPUTE_ACC();

    int ic = Nx / 2;  /* center index */

    /* --- Envelope tracking --- */
    /* Track the peak amplitude envelope by recording local maxima over time */
    double *envelope = calloc(Nx, sizeof(double));  /* peak |phi| seen */
    double *phi_snapshot = calloc(Nx, sizeof(double));  /* snapshot at peak */

    /* We'll record the envelope during the last portion of the simulation */
    double t_env_start = tfinal * 0.8;  /* start tracking envelope at 80% */
    int envelope_active = 0;

    /* DFT storage for frequency measurement */
    int max_dft = 50000;
    double *phi0_hist = malloc(max_dft * sizeof(double));
    double *t_hist = malloc(max_dft * sizeof(double));
    int n_dft = 0;
    int dft_every = Nt / max_dft;
    if (dft_every < 1) dft_every = 1;

    int print_every = Nt / 40;
    if (print_every < 1) print_every = 1;

    /* --- Time evolution (Velocity Verlet) --- */
    printf("\nEvolving 1D symmetric oscillon...\n");

    for (int n = 0; n <= Nt; n++) {
        double t = n * dt;

        /* DFT recording */
        if (n % dft_every == 0 && n_dft < max_dft) {
            phi0_hist[n_dft] = phi[ic];
            t_hist[n_dft] = t;
            n_dft++;
        }

        /* Envelope tracking */
        if (t >= t_env_start) {
            if (!envelope_active) {
                envelope_active = 1;
                printf("  Envelope tracking started at t=%.1f\n", t);
            }
            for (int i = 0; i < Nx; i++) {
                double ap = fabs(phi[i]);
                if (ap > envelope[i]) {
                    envelope[i] = ap;
                }
            }
        }

        /* Print progress */
        if (n % print_every == 0) {
            double peak = 0;
            for (int i = 0; i < Nx; i++)
                if (fabs(phi[i]) > peak) peak = fabs(phi[i]);
            printf("  t=%7.1f  phi(0)=%+.4f  peak=%.4f\n", t, phi[ic], peak);
        }

        if (n == Nt) break;

        /* Velocity Verlet */
        for (int i = 1; i < Nx - 1; i++)
            vel[i] += 0.5 * dt * acc[i];
        for (int i = 1; i < Nx - 1; i++)
            phi[i] += dt * vel[i];
        COMPUTE_ACC();
        for (int i = 1; i < Nx - 1; i++)
            vel[i] += 0.5 * dt * acc[i];

        /* Absorbing boundary */
        for (int i = 0; i < Nx; i++) {
            vel[i] *= damp[i];
            phi[i] *= damp[i];
        }
    }

    /* --- Measure frequency from DFT --- */
    int dft_start = n_dft / 2;
    double omega_peak = 0;
    if (n_dft - dft_start > 100) {
        double T = t_hist[n_dft-1] - t_hist[dft_start];
        int nf = 1000;
        double peak_pow = 0;
        for (int k = 1; k < nf; k++) {
            double omega = 3.0 * mass * k / nf;
            double re = 0, im = 0;
            for (int j = dft_start; j < n_dft; j++) {
                double dtj = (j > dft_start) ?
                    (t_hist[j]-t_hist[j-1]) : (t_hist[dft_start+1]-t_hist[dft_start]);
                re += phi0_hist[j] * cos(omega * t_hist[j]) * dtj;
                im += phi0_hist[j] * sin(omega * t_hist[j]) * dtj;
            }
            double pw = (re*re + im*im) / (T*T);
            if (pw > peak_pow) { peak_pow = pw; omega_peak = omega; }
        }
        printf("\nOscillon frequency: omega = %.4f (mass gap m = %.4f)\n",
               omega_peak, mass);
        printf("Sub-gap? %s (omega/m = %.4f)\n",
               (omega_peak < mass) ? "YES" : "NO", omega_peak / mass);
    }

    /* --- Use ONLY positive-x half (r >= 0) since symmetric --- */
    int Nr = Nx / 2 + 1;  /* number of radial points: x=0 to x=xmax */
    double *r_arr  = malloc(Nr * sizeof(double));
    double *f_arr  = malloc(Nr * sizeof(double));  /* envelope = peak amplitude */

    for (int j = 0; j < Nr; j++) {
        int i = ic + j;  /* index from center outward */
        r_arr[j] = j * dx;
        f_arr[j] = (i < Nx) ? envelope[i] : 0.0;
    }

    printf("\nEnvelope: f(0)=%.6f, f(5)=%.6f, f(10)=%.6f\n",
           f_arr[0],
           f_arr[(int)(5.0/dx)],
           f_arr[(int)(10.0/dx)]);

    /* --- Compute Hessian eigenvalues --- */
    printf("\nComputing Hessian eigenvalues...\n");

    /* Output file: instantaneous Hessian (using peak envelope) */
    char path1[600];
    snprintf(path1, sizeof(path1), "%s/hessian_instant.tsv", outdir);
    FILE *f1 = fopen(path1, "w");
    if (!f1) { fprintf(stderr, "Cannot open %s\n", path1); return 1; }
    fprintf(f1, "r\tf\tP\tlambda_sym\tlambda_anti\tm2_sym\tm2_anti\tdelta_m2\n");

    /* Output file: time-averaged Hessian */
    char path2[600];
    snprintf(path2, sizeof(path2), "%s/hessian_tavg.tsv", outdir);
    FILE *f2 = fopen(path2, "w");
    if (!f2) { fprintf(stderr, "Cannot open %s\n", path2); return 1; }
    fprintf(f2, "r\tf\tlambda_sym_avg\tlambda_anti_avg\tm2_sym_avg\tm2_anti_avg\tdelta_m2_avg\n");

    double m2_anti_min = 1e30;
    double r_anti_min = 0;
    double m2_anti_avg_min = 1e30;
    double r_anti_avg_min = 0;
    double f_at_min = 0;

    double m2_sym_min = 1e30;
    double r_sym_min = 0;
    double m2_sym_avg_min = 1e30;
    double r_sym_avg_min = 0;

    for (int j = 0; j < Nr; j++) {
        double r = r_arr[j];
        double f = f_arr[j];
        double f2v = f * f;
        double f4 = f2v * f2v;
        double f6 = f4 * f2v;
        double P = f * f * f;  /* f³ */

        double kf6 = kappa * f6;
        double denom1 = 1.0 + kf6;
        double denom2 = denom1 * denom1;
        double denom3 = denom2 * denom1;

        /* Instantaneous Hessian eigenvalues (at peak amplitude) */
        double lam_sym  = mu * f4 * (5.0 - 7.0 * kf6) / denom3;
        double lam_anti = -mu * f4 / denom2;

        double m2_sym  = m2 + lam_sym;
        double m2_anti = m2 + lam_anti;
        double delta   = m2_sym - m2_anti;

        fprintf(f1, "%.6f\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\n",
                r, f, P, lam_sym, lam_anti, m2_sym, m2_anti, delta);

        if (m2_anti < m2_anti_min) {
            m2_anti_min = m2_anti;
            r_anti_min = r;
            f_at_min = f;
        }
        if (m2_sym < m2_sym_min) {
            m2_sym_min = m2_sym;
            r_sym_min = r;
        }

        /* Time-averaged: ⟨f⁴cos⁴(ωt)⟩ = (3/8)f⁴ */
        /* The Hessian eigenvalues are proportional to f⁴ times functions of f⁶.
         * For the time average, the oscillon profile is f(r)cos(ωt), so
         * f → f·cos(ωt) in the expressions. We need ⟨f⁴cos⁴⟩ and ⟨f⁶cos⁶⟩.
         *
         * More precisely, with φ = A(r)cos(ωt):
         *   φ⁴ → A⁴⟨cos⁴⟩ = (3/8)A⁴
         *   φ⁶ → A⁶⟨cos⁶⟩ = (5/16)A⁶
         *
         * But the Hessian has f⁴/(1+κf⁶)^n which is NOT simply f⁴ times const.
         * We need ⟨H(A cos(ωt))⟩ over a period.
         *
         * For proper time averaging, do a numerical average over one period.
         */
        int Ntheta = 1000;
        double dtheta = 2.0 * M_PI / Ntheta;
        double avg_lam_sym = 0, avg_lam_anti = 0;

        for (int k = 0; k < Ntheta; k++) {
            double theta = k * dtheta;
            double ct = cos(theta);
            double phi_t = f * ct;  /* f is the peak envelope */
            double pt2 = phi_t * phi_t;
            double pt4 = pt2 * pt2;
            double pt6 = pt4 * pt2;
            double kpt6 = kappa * pt6;
            double d1 = 1.0 + kpt6;
            double d2 = d1 * d1;
            double d3 = d2 * d1;

            avg_lam_sym  += mu * pt4 * (5.0 - 7.0 * kpt6) / d3;
            avg_lam_anti += -mu * pt4 / d2;
        }
        avg_lam_sym  /= Ntheta;
        avg_lam_anti /= Ntheta;

        double m2_sym_avg  = m2 + avg_lam_sym;
        double m2_anti_avg = m2 + avg_lam_anti;
        double delta_avg   = m2_sym_avg - m2_anti_avg;

        fprintf(f2, "%.6f\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\n",
                r, f, avg_lam_sym, avg_lam_anti, m2_sym_avg, m2_anti_avg, delta_avg);

        if (m2_anti_avg < m2_anti_avg_min) {
            m2_anti_avg_min = m2_anti_avg;
            r_anti_avg_min = r;
        }
        if (m2_sym_avg < m2_sym_avg_min) {
            m2_sym_avg_min = m2_sym_avg;
            r_sym_avg_min = r;
        }
    }

    fclose(f1);
    fclose(f2);

    /* --- Summary output --- */
    printf("\n=== HESSIAN SPLITTING RESULTS ===\n\n");

    printf("Parameters: mu=%.1f, kappa=%.1f, m=%.3f\n", mu, kappa, mass);
    printf("Oscillon frequency: omega=%.4f (omega/m=%.4f)\n", omega_peak, omega_peak/mass);
    printf("Peak envelope amplitude: f(0) = %.6f\n\n", f_arr[0]);

    /* Instantaneous values at center */
    {
        double f = f_arr[0];
        double f4 = f*f*f*f;
        double f6 = f4*f*f;
        double kf6 = kappa * f6;
        double d1 = 1.0 + kf6;
        double d2 = d1*d1;
        double d3 = d2*d1;

        double ls = mu * f4 * (5.0 - 7.0*kf6) / d3;
        double la = -mu * f4 / d2;

        printf("--- Instantaneous (at peak amplitude) ---\n");
        printf("  f(0)           = %.6f\n", f);
        printf("  f⁶(0)          = %.6f\n", f6);
        printf("  κf⁶(0)         = %.6f\n", kf6);
        printf("  λ_sym(0)       = %+.6f\n", ls);
        printf("  λ_anti(0)      = %+.6f\n", la);
        printf("  m²_sym(0)      = %+.6f\n", m2 + ls);
        printf("  m²_anti(0)     = %+.6f\n", m2 + la);
        printf("  Δm²(0)         = %+.6f\n", (m2+ls) - (m2+la));
        printf("  m²_sym min     = %+.6f at r=%.3f\n", m2_sym_min, r_sym_min);
        printf("  m²_anti min    = %+.6f at r=%.3f (f=%.4f)\n", m2_anti_min, r_anti_min, f_at_min);
        printf("  Sym tachyonic? %s\n", (m2_sym_min < 0) ? "YES" : "NO");
        printf("  Anti tachyonic? %s\n", (m2_anti_min < 0) ? "YES" : "NO");
        printf("\n");
    }

    /* Time-averaged values at center (recompute for display) */
    {
        double f = f_arr[0];
        int Ntheta = 10000;
        double dtheta = 2.0 * M_PI / Ntheta;
        double avg_ls = 0, avg_la = 0;
        for (int k = 0; k < Ntheta; k++) {
            double theta = k * dtheta;
            double ct = cos(theta);
            double pt = f * ct;
            double pt2 = pt*pt;
            double pt4 = pt2*pt2;
            double pt6 = pt4*pt2;
            double kpt6 = kappa*pt6;
            double d1 = 1.0 + kpt6;
            double d2 = d1*d1;
            double d3 = d2*d1;
            avg_ls += mu * pt4 * (5.0 - 7.0*kpt6) / d3;
            avg_la += -mu * pt4 / d2;
        }
        avg_ls /= Ntheta;
        avg_la /= Ntheta;

        printf("--- Time-averaged (⟨H(f cos ωt)⟩) ---\n");
        printf("  ⟨λ_sym⟩(0)     = %+.6f\n", avg_ls);
        printf("  ⟨λ_anti⟩(0)    = %+.6f\n", avg_la);
        printf("  ⟨m²_sym⟩(0)    = %+.6f\n", m2 + avg_ls);
        printf("  ⟨m²_anti⟩(0)   = %+.6f\n", m2 + avg_la);
        printf("  ⟨Δm²⟩(0)       = %+.6f\n", (m2+avg_ls) - (m2+avg_la));
        printf("  ⟨m²_sym⟩ min   = %+.6f at r=%.3f\n", m2_sym_avg_min, r_sym_avg_min);
        printf("  ⟨m²_anti⟩ min  = %+.6f at r=%.3f\n", m2_anti_avg_min, r_anti_avg_min);
        printf("  Sym tachyonic (avg)? %s\n", (m2_sym_avg_min < 0) ? "YES" : "NO");
        printf("  Anti tachyonic (avg)? %s\n", (m2_anti_avg_min < 0) ? "YES" : "NO");
        printf("\n");
    }

    /* Splitting analysis: scan for sign changes and key features */
    printf("--- Splitting Profile ---\n");
    double r_half = 0;
    {
        double delta_0 = 0;
        /* Compute Δm² at r=0 */
        double f = f_arr[0];
        double f4 = f*f*f*f; double f6 = f4*f*f;
        double kf6 = kappa*f6;
        double d1 = 1.0+kf6; double d2=d1*d1; double d3=d2*d1;
        double ls = mu*f4*(5.0-7.0*kf6)/d3;
        double la = -mu*f4/d2;
        delta_0 = ls - la;

        /* Find half-width of splitting */
        for (int j = 1; j < Nr; j++) {
            double fj = f_arr[j];
            double fj4 = fj*fj*fj*fj; double fj6 = fj4*fj*fj;
            double kfj6 = kappa*fj6;
            double dj1 = 1.0+kfj6; double dj2=dj1*dj1; double dj3=dj2*dj1;
            double lsj = mu*fj4*(5.0-7.0*kfj6)/dj3;
            double laj = -mu*fj4/dj2;
            double dj = lsj - laj;
            if (fabs(dj) < 0.5 * fabs(delta_0) && r_half == 0) {
                r_half = r_arr[j];
                break;
            }
        }
        printf("  Δm²(0)     = %+.6f\n", delta_0);
        printf("  Δm² half-width r₁/₂ = %.3f\n", r_half);
    }

    /* Check where 7κf⁶ = 5 (sign change of λ_sym) */
    {
        double f_crit = pow(5.0 / (7.0 * kappa), 1.0/6.0);
        printf("  λ_sym sign change at f = %.6f (7κf⁶=5)\n", f_crit);
        /* Find corresponding r */
        double r_crit = -1;
        for (int j = 0; j < Nr - 1; j++) {
            if (f_arr[j] >= f_crit && f_arr[j+1] < f_crit) {
                double frac = (f_arr[j] - f_crit) / (f_arr[j] - f_arr[j+1]);
                r_crit = r_arr[j] + frac * dx;
                break;
            }
        }
        if (r_crit >= 0)
            printf("  λ_sym sign change at r = %.3f\n", r_crit);
        else
            printf("  λ_sym does not change sign (f(0) < f_crit or f(0) > f_crit everywhere)\n");
        printf("  f(0) %s f_crit → λ_sym(0) is %s\n",
               (f_arr[0] > f_crit) ? ">" : "<=",
               (f_arr[0] > f_crit) ?
                   (mu < 0 ? "POSITIVE (saturation dominates)" : "NEGATIVE") :
                   (mu < 0 ? "NEGATIVE" : "POSITIVE"));
    }

    printf("\nOutput files:\n");
    printf("  %s\n", path1);
    printf("  %s\n", path2);

    printf("\nDone. Write RESULTS.md manually from the above output.\n");

    /* Cleanup */
    free(phi); free(vel); free(acc); free(damp);
    free(envelope); free(phi_snapshot);
    free(phi0_hist); free(t_hist);
    free(r_arr); free(f_arr);

    return 0;
}
