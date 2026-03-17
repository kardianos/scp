/*
 * maxwell_a.c — Three-field oscillon with U(1) gauge coupling
 *
 * Fields: φ₁, φ₂ (charged under U(1)), φ₃ (neutral), A (gauge field)
 * Temporal gauge: A₀ = 0
 *
 * Lagrangian:
 *   L = ½|D_t Ψ|² - ½|D_x Ψ|² - ½m²|Ψ|²
 *     + ½(∂_t φ₃)² - ½(∂_x φ₃)² - ½m²φ₃²
 *     - V_coupling - ¼F_{μν}F^{μν}
 *
 * where Ψ = φ₁ + iφ₂, D_μ = ∂_μ - ieA_μ
 *
 * Coupling: V = (μ/2) P² / (1 + κ P²),  P = φ₁ φ₂ φ₃
 *   (original triple product — NOT gauge-invariant, but oscillon is proven)
 *
 * NOTE: The triple product P = φ₁φ₂φ₃ breaks U(1) gauge invariance.
 * This means Q is NOT exactly conserved even classically.
 * We study this as a diagnostic: does the gauge field develop structure?
 * Does Q remain approximately conserved?
 *
 * EOMs (component form, A₀=0):
 *   ∂²φ₁/∂t² = ∂²φ₁/∂x² + e(∂_xA)φ₂ + 2eA(∂_xφ₂) + e²A²φ₁
 *               - m²φ₁ + fp₁
 *   ∂²φ₂/∂t² = ∂²φ₂/∂x² - e(∂_xA)φ₁ - 2eA(∂_xφ₁) + e²A²φ₂
 *               - m²φ₂ + fp₂
 *   ∂²φ₃/∂t² = ∂²φ₃/∂x² - m²φ₃ + fp₃
 *   ∂²A/∂t²  = ∂²A/∂x² - j
 *
 * where fp_a = -∂V/∂φ_a (triple-product force)
 *       j = e(φ₁∂_xφ₂ - φ₂∂_xφ₁) - e²A(φ₁²+φ₂²)
 *
 * Charge: Q = ∫(φ₁ v₂ - φ₂ v₁) dx
 *
 * Compile: gcc -O3 -Wall -o maxwell_a src/maxwell_a.c -lm
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Parameters — v21 optimized: μ=-20, κ=20, m=1.0, A=0.8, σ=3.0 */
static double mu      = -20.0;
static double kappa   = 20.0;
static double mass    = 1.0;
static double A_init  = 0.8;
static double sig     = 3.0;
static double echarge = 0.0;   /* gauge coupling */
static double twist   = 0.3;   /* v₂ kick amplitude (fraction of v₁) */
static int    Nx      = 4000;
static double xmax    = 100.0;
static double tfinal  = 10000.0;
static char   outdir[512] = "v24/maxwell_a/data";

static void parse_args(int argc, char **argv)
{
    for (int i = 1; i < argc - 1; i += 2) {
        if      (!strcmp(argv[i], "-mu"))     mu      = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-kappa"))  kappa   = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-mass"))   mass    = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-A"))      A_init  = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-sigma"))  sig     = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-e"))      echarge = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-twist"))  twist   = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-Nx"))     Nx      = atoi(argv[i+1]);
        else if (!strcmp(argv[i], "-xmax"))   xmax    = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-tfinal")) tfinal  = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-o"))      strncpy(outdir, argv[i+1], sizeof(outdir)-1);
        else { fprintf(stderr, "Unknown: %s\n", argv[i]); exit(1); }
    }
}

/* Triple product force: -dV/dphi_a where V = (mu/2)P^2/(1+kappa P^2) */
static inline double force_pot(double p1, double p2, double p3, int a)
{
    double P = p1 * p2 * p3;
    double P2 = P * P;
    double denom2 = (1.0 + kappa * P2); denom2 *= denom2;
    double dP;
    switch (a) {
        case 0: dP = p2 * p3; break;
        case 1: dP = p1 * p3; break;
        case 2: dP = p1 * p2; break;
        default: dP = 0.0;
    }
    return -mu * P * dP / denom2;
}

int main(int argc, char **argv)
{
    parse_args(argc, argv);

    double dx = 2.0 * xmax / (Nx - 1);
    double dx2 = dx * dx;
    double m2 = mass * mass;
    double e = echarge;
    double e2 = e * e;

    /* CFL: use 0.4 safety factor */
    double kmax = M_PI / dx;
    double dt = 0.4 * 2.0 / sqrt(kmax*kmax + m2);
    int Nt = (int)(tfinal / dt) + 1;

    printf("maxwell_a: triad oscillon + U(1) gauge (triple product coupling)\n");
    printf("  mu=%.1f kappa=%.1f mass=%.3f A=%.3f sigma=%.3f e=%.4f twist=%.3f\n",
           mu, kappa, mass, A_init, sig, e, twist);
    printf("  Nx=%d xmax=%.1f dx=%.5f dt=%.6f Nt=%d tfinal=%.0f\n",
           Nx, xmax, dx, dt, Nt, tfinal);

    /* Allocate fields */
    double *phi1 = calloc(Nx, sizeof(double));
    double *phi2 = calloc(Nx, sizeof(double));
    double *phi3 = calloc(Nx, sizeof(double));
    double *Af   = calloc(Nx, sizeof(double));

    double *v1 = calloc(Nx, sizeof(double));
    double *v2 = calloc(Nx, sizeof(double));
    double *v3 = calloc(Nx, sizeof(double));
    double *vA = calloc(Nx, sizeof(double));

    double *a1 = calloc(Nx, sizeof(double));
    double *a2 = calloc(Nx, sizeof(double));
    double *a3 = calloc(Nx, sizeof(double));
    double *aA = calloc(Nx, sizeof(double));

    /* Absorbing boundary */
    double *damp = malloc(Nx * sizeof(double));
    double x_abs = xmax * 0.75;
    for (int i = 0; i < Nx; i++) {
        double x = -xmax + i * dx;
        double ax = fabs(x);
        if (ax > x_abs)
            damp[i] = 1.0 - 0.98 * ((ax-x_abs)/(xmax-x_abs)) * ((ax-x_abs)/(xmax-x_abs));
        else
            damp[i] = 1.0;
    }

    /* Initialize: standard triad oscillon (all three fields equal Gaussians).
     * Add small velocity kick to φ₂ to seed U(1) charge.
     * φ₁ = φ₂ = φ₃ = A·g(x), v₁ = v₃ = 0
     * v₂ = twist · ω₀ · A · g(x) where ω₀ ≈ 0.76 is the oscillon frequency
     *
     * This makes φ₂ slightly out of phase with φ₁, creating a rotating
     * Ψ = φ₁ + iφ₂ that carries U(1) charge.
     * Q(0) = ∫ φ₁ v₂ dx = twist·ω₀·A² ∫g² dx > 0
     *
     * Gauss constraint: at t=0, ρ = e(φ₁v₂ - φ₂v₁) = e·φ₁·twist·ω₀·φ₁
     * We set E(x) = ∫₋∞ˣ ρ dx' to satisfy ∂_x E = ρ.
     */
    double omega0 = 0.76 * mass;
    for (int i = 0; i < Nx; i++) {
        double x = -xmax + i * dx;
        double g = A_init * exp(-x*x / (2.0*sig*sig));
        phi1[i] = g;
        phi2[i] = g;  /* all three equal for triad binding */
        phi3[i] = g;
        Af[i] = 0.0;
        v1[i] = 0.0;
        v2[i] = twist * omega0 * g;  /* velocity kick on φ₂ only */
        v3[i] = 0.0;
        vA[i] = 0.0;
    }

    /* Fix Gauss constraint: vA = -∫₋∞ˣ ρ dx' where ρ = e(φ₁v₂ - φ₂v₁) */
    if (fabs(e) > 1e-15) {
        double Esum = 0.0;
        for (int i = 0; i < Nx; i++) {
            double rho_i = e * (phi1[i]*v2[i] - phi2[i]*v1[i]);
            Esum += rho_i * dx;
            vA[i] = -Esum;
        }
        printf("  Gauss fix: net charge = %.6e, E_boundary = %.6e\n",
               Esum/e, Esum);
    }

    /* Initial Q */
    {
        double Q0 = 0;
        for (int i = 0; i < Nx; i++)
            Q0 += (phi1[i]*v2[i] - phi2[i]*v1[i]) * dx;
        printf("  Initial Q = %.6f\n", Q0);
    }

    /* Compute accelerations */
    #define COMPUTE_ACC() do { \
        for (int b = 0; b < 2; b++) { \
            a1[b]=a2[b]=a3[b]=aA[b]=0; \
            a1[Nx-1-b]=a2[Nx-1-b]=a3[Nx-1-b]=aA[Nx-1-b]=0; \
        } \
        for (int i = 2; i < Nx - 2; i++) { \
            double L1 = (phi1[i+1] - 2.0*phi1[i] + phi1[i-1]) / dx2; \
            double L2 = (phi2[i+1] - 2.0*phi2[i] + phi2[i-1]) / dx2; \
            double L3 = (phi3[i+1] - 2.0*phi3[i] + phi3[i-1]) / dx2; \
            double LA = (Af[i+1] - 2.0*Af[i] + Af[i-1]) / dx2; \
            double dxA = (Af[i+1] - Af[i-1]) / (2.0*dx); \
            double dxp1 = (phi1[i+1] - phi1[i-1]) / (2.0*dx); \
            double dxp2 = (phi2[i+1] - phi2[i-1]) / (2.0*dx); \
            double Ai = Af[i]; \
            double psi2 = phi1[i]*phi1[i] + phi2[i]*phi2[i]; \
            double fp1 = force_pot(phi1[i], phi2[i], phi3[i], 0); \
            double fp2 = force_pot(phi1[i], phi2[i], phi3[i], 1); \
            double fp3 = force_pot(phi1[i], phi2[i], phi3[i], 2); \
            /* φ₁: Laplacian + gauge terms + mass + potential */ \
            a1[i] = L1 + e*dxA*phi2[i] + 2.0*e*Ai*dxp2 + e2*Ai*Ai*phi1[i] \
                    - m2*phi1[i] + fp1; \
            /* φ₂: Laplacian + gauge terms (opposite sign) + mass + potential */ \
            a2[i] = L2 - e*dxA*phi1[i] - 2.0*e*Ai*dxp1 + e2*Ai*Ai*phi2[i] \
                    - m2*phi2[i] + fp2; \
            /* φ₃: standard KG + potential */ \
            a3[i] = L3 - m2*phi3[i] + fp3; \
            /* A: Maxwell with current */ \
            double j_curr = e*(phi1[i]*dxp2 - phi2[i]*dxp1) - e2*Ai*psi2; \
            aA[i] = LA - j_curr; \
        } \
    } while(0)

    COMPUTE_ACC();

    /* Output files */
    char tspath[600];
    snprintf(tspath, sizeof(tspath), "%s/maxwell_e%.4f_ts.tsv", outdir, e);
    FILE *fts = fopen(tspath, "w");
    if (!fts) { fprintf(stderr, "Cannot open %s\n", tspath); return 1; }
    fprintf(fts, "time\tphi1_0\tphi2_0\tphi3_0\tA_0\tpeak1\tpeak2\tpeak3\tpeak_A\t"
                 "E_kin\tE_grad\tE_mass\tE_pot\tE_EM\tE_total\tQ_charge\t"
                 "f_core\tgauss_err\n");

    /* DFT storage */
    int max_dft = 50000;
    double *phi0_hist = malloc(max_dft * sizeof(double));
    double *t_hist = malloc(max_dft * sizeof(double));
    int n_dft = 0;

    int rec_every = Nt / 20000;
    if (rec_every < 1) rec_every = 1;
    int print_every = Nt / 40;
    if (print_every < 1) print_every = 1;
    int dft_every = Nt / max_dft;
    if (dft_every < 1) dft_every = 1;

    double core_r = 3.0 * sig;
    int ic = Nx / 2;

    for (int n = 0; n <= Nt; n++) {
        double t = n * dt;

        if (n % dft_every == 0 && n_dft < max_dft) {
            phi0_hist[n_dft] = phi1[ic];
            t_hist[n_dft] = t;
            n_dft++;
        }

        int do_rec = (n % rec_every == 0);
        int do_print = (n % print_every == 0);

        if (do_rec || do_print) {
            double Ek=0, Eg=0, Em=0, Ep=0, Eem=0;
            double Ecore=0, Eall=0;
            double peak[3]={0}, peak_A=0;
            double Q = 0, gauss_max = 0;

            for (int i = 1; i < Nx - 1; i++) {
                double x = -xmax + i * dx;

                /* Kinetic */
                double ek_i = 0.5*(v1[i]*v1[i] + v2[i]*v2[i] + v3[i]*v3[i]);
                Ek += ek_i * dx;

                /* Gradient: use covariant for charged fields */
                double dxp1 = (phi1[i+1]-phi1[i-1])/(2.0*dx);
                double dxp2 = (phi2[i+1]-phi2[i-1])/(2.0*dx);
                double dxp3 = (phi3[i+1]-phi3[i-1])/(2.0*dx);
                double Dx_re = dxp1 - e*Af[i]*phi2[i];
                double Dx_im = dxp2 + e*Af[i]*phi1[i];
                double eg_i = 0.5*(Dx_re*Dx_re + Dx_im*Dx_im) + 0.5*dxp3*dxp3;
                Eg += eg_i * dx;

                /* Mass */
                double em_i = 0.5*m2*(phi1[i]*phi1[i]+phi2[i]*phi2[i]+phi3[i]*phi3[i]);
                Em += em_i * dx;

                /* Triple product potential */
                double P = phi1[i]*phi2[i]*phi3[i];
                double P2 = P*P;
                double V = 0.5*mu*P2/(1.0+kappa*P2);
                Ep += V * dx;

                /* EM energy: ½E² = ½vA² */
                double eem_i = 0.5*vA[i]*vA[i];
                Eem += eem_i * dx;

                /* Peaks */
                for (int a = 0; a < 3; a++) {
                    double *phi_a = (a==0)?phi1:(a==1)?phi2:phi3;
                    if (fabs(phi_a[i]) > peak[a]) peak[a] = fabs(phi_a[i]);
                }
                if (fabs(Af[i]) > peak_A) peak_A = fabs(Af[i]);

                /* Charge */
                Q += (phi1[i]*v2[i] - phi2[i]*v1[i]) * dx;

                /* Gauss check */
                if (i > 1 && i < Nx-2) {
                    double dxE = -(vA[i+1]-vA[i-1])/(2.0*dx);
                    double rho = e*(phi1[i]*v2[i] - phi2[i]*v1[i]);
                    double ge = fabs(dxE - rho);
                    if (ge > gauss_max) gauss_max = ge;
                }

                /* Total for core fraction */
                double etot_i = ek_i + eg_i + em_i + V + eem_i;
                Eall += etot_i * dx;
                if (fabs(x) < core_r) Ecore += etot_i * dx;
            }

            double Et = Ek + Eg + Em + Ep + Eem;
            double fc = (Eall > 1e-20) ? Ecore/Eall : 0.0;

            if (do_rec)
                fprintf(fts, "%.6f\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t"
                             "%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6f\t%.3e\n",
                        t, phi1[ic], phi2[ic], phi3[ic], Af[ic],
                        peak[0], peak[1], peak[2], peak_A,
                        Ek, Eg, Em, Ep, Eem, Et, Q, fc, gauss_max);

            if (do_print)
                printf("  t=%7.1f  pk=(%.3f,%.3f,%.3f) A=%.5f  "
                       "E=%+.3f Ep=%+.3f Eem=%.4f  Q=%+.5f  fc=%.3f  g=%.1e\n",
                       t, peak[0], peak[1], peak[2], peak_A,
                       Et, Ep, Eem, Q, fc, gauss_max);
        }

        if (n == Nt) break;

        /* Velocity Verlet */
        for (int i = 2; i < Nx-2; i++) {
            v1[i] += 0.5*dt*a1[i]; v2[i] += 0.5*dt*a2[i];
            v3[i] += 0.5*dt*a3[i]; vA[i] += 0.5*dt*aA[i];
        }
        for (int i = 2; i < Nx-2; i++) {
            phi1[i] += dt*v1[i]; phi2[i] += dt*v2[i];
            phi3[i] += dt*v3[i]; Af[i] += dt*vA[i];
        }
        COMPUTE_ACC();
        for (int i = 2; i < Nx-2; i++) {
            v1[i] += 0.5*dt*a1[i]; v2[i] += 0.5*dt*a2[i];
            v3[i] += 0.5*dt*a3[i]; vA[i] += 0.5*dt*aA[i];
        }

        /* Absorbing boundary */
        for (int i = 0; i < Nx; i++) {
            phi1[i] *= damp[i]; v1[i] *= damp[i];
            phi2[i] *= damp[i]; v2[i] *= damp[i];
            phi3[i] *= damp[i]; v3[i] *= damp[i];
            Af[i] *= damp[i]; vA[i] *= damp[i];
        }
    }

    fclose(fts);

    /* DFT */
    int dft_start = n_dft / 2;
    if (n_dft - dft_start > 100) {
        char dftpath[600];
        snprintf(dftpath, sizeof(dftpath), "%s/maxwell_e%.4f_spectrum.tsv", outdir, e);
        FILE *fdft = fopen(dftpath, "w");
        fprintf(fdft, "omega\tpower\n");
        double T = t_hist[n_dft-1] - t_hist[dft_start];
        int nf = 500;
        double peak_pow=0, peak_om=0;
        for (int k = 0; k < nf; k++) {
            double omega = 3.0*mass*k/nf;
            double re=0, im=0;
            for (int j = dft_start; j < n_dft; j++) {
                double dtj = (j>dft_start) ?
                    (t_hist[j]-t_hist[j-1]) : (t_hist[dft_start+1]-t_hist[dft_start]);
                re += phi0_hist[j]*cos(omega*t_hist[j])*dtj;
                im += phi0_hist[j]*sin(omega*t_hist[j])*dtj;
            }
            double pw = (re*re+im*im)/(T*T);
            fprintf(fdft, "%.6f\t%.6e\n", omega, pw);
            if (pw > peak_pow) { peak_pow = pw; peak_om = omega; }
        }
        fclose(fdft);
        printf("\nSpectrum: peak omega = %.4f (mass = %.4f)\n", peak_om, mass);
        printf("Oscillon? %s\n", (peak_om>0.01 && peak_om<mass) ? "YES" : "NO");
    }

    /* Profile snapshot */
    {
        char pp[600];
        snprintf(pp, sizeof(pp), "%s/maxwell_e%.4f_profile.tsv", outdir, e);
        FILE *fp = fopen(pp, "w");
        fprintf(fp, "x\tphi1\tphi2\tphi3\tA\n");
        for (int i = 0; i < Nx; i += 4) {
            double x = -xmax + i*dx;
            fprintf(fp, "%.4f\t%.6e\t%.6e\t%.6e\t%.6e\n",
                    x, phi1[i], phi2[i], phi3[i], Af[i]);
        }
        fclose(fp);
    }

    printf("Output: %s\n", tspath);

    free(phi1); free(phi2); free(phi3); free(Af);
    free(v1); free(v2); free(v3); free(vA);
    free(a1); free(a2); free(a3); free(aA);
    free(damp); free(phi0_hist); free(t_hist);
    return 0;
}
