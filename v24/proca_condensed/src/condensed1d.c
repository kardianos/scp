/*
 * condensed1d.c — Condensed phase Goldstone search
 *
 * Three massive scalars with pairwise + triple-product + quartic self-coupling:
 *   L = sum_a [ (1/2)(dt phi_a)^2 - (1/2)(dx phi_a)^2 - (m^2/2)phi_a^2 ]
 *     - lambda(phi_1 phi_2 + phi_2 phi_3 + phi_3 phi_1)
 *     - (mu/2) P^2 / (1 + kappa P^2)    [P = phi_1 phi_2 phi_3]
 *     - (g4/4) (sum_a phi_a^2)^2         [O(3)-symmetric quartic stabilization]
 *
 * When lambda > m^2: antisymmetric modes are tachyonic in phi=0 vacuum.
 * The quartic g4 stabilizes the condensate at finite VEV.
 * Condensate breaks SO(2) in the antisymmetric subspace -> Goldstone mode.
 *
 * EOM:  d^2 phi_a/dt^2 = d^2 phi_a/dx^2 - m^2 phi_a - lambda(phi_b + phi_c)
 *                        - g4 phi_a (sum phi^2) - mu P dP/dphi_a / (1+kappa P^2)^2
 *
 * Phase 1: Find condensed vacuum (gradient descent for uniform fields)
 * Phase 2: Fluctuation spectrum around condensate (Goldstone vs massive)
 * Phase 3: Oscillon on top of condensate
 *
 * Compile: gcc -O3 -Wall -o condensed1d src/condensed1d.c -lm
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* ---- Parameters ---- */
static double mu      = -20.0;
static double kappa   = 20.0;
static double mass    = 1.0;
static double lambda  = 1.1;
static double g4      = 1.0;   /* quartic self-coupling for condensate stabilization */
static double A_init  = 0.8;
static double sigma_w = 3.0;
static int    Nx      = 8000;
static double xmax    = 200.0;
static double tfinal  = 10000.0;
static int    phase   = 0;     /* 0=all, 1=vacuum, 2=spectrum, 3=oscillon */
static char   outdir[512] = "v24/proca_condensed/data";

static void parse_args(int argc, char **argv)
{
    for (int i = 1; i < argc - 1; i += 2) {
        if      (!strcmp(argv[i], "-mu"))     mu      = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-kappa"))  kappa   = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-mass"))   mass    = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-lambda")) lambda  = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-g4"))     g4      = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-A"))      A_init  = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-sigma"))  sigma_w = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-Nx"))     Nx      = atoi(argv[i+1]);
        else if (!strcmp(argv[i], "-xmax"))   xmax    = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-tfinal")) tfinal  = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-phase"))  phase   = atoi(argv[i+1]);
        else if (!strcmp(argv[i], "-o"))      strncpy(outdir, argv[i+1], sizeof(outdir)-1);
        else { fprintf(stderr, "Unknown flag: %s\n", argv[i]); exit(1); }
    }
}

/* -dV_triple/dphi_a where V = (mu/2)P^2/(1+kappa P^2), P = phi1 phi2 phi3 */
static inline double force_triple(double p1, double p2, double p3, int a)
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

/* Total potential energy density at a point */
static inline double pot_density(double p1, double p2, double p3,
                                  double m2, double lam)
{
    double S2 = p1*p1 + p2*p2 + p3*p3;
    double V = 0.5 * m2 * S2;
    V += lam * (p1*p2 + p2*p3 + p3*p1);
    V += (g4/4.0) * S2 * S2;
    double P = p1 * p2 * p3;
    V += 0.5 * mu * P * P / (1.0 + kappa * P * P);
    return V;
}

/* ======================================================================
 * Find the condensed vacuum numerically for uniform fields.
 * Minimize V(phi1,phi2,phi3) with gradient descent.
 * Antisymmetric seed: phi1=v, phi2=-v/2, phi3=-v/2
 * ====================================================================== */
static void find_vacuum(double lam, double vev_out[3], double *V_min)
{
    double m2 = mass * mass;
    double phi[3];

    /* Seed: small antisymmetric perturbation */
    double v0 = 0.3;
    phi[0] = v0;
    phi[1] = -v0 / 2.0;
    phi[2] = -v0 / 2.0;

    /* Gradient descent with adaptive step */
    double eta = 0.001;
    for (int iter = 0; iter < 500000; iter++) {
        double grad[3];
        for (int a = 0; a < 3; a++) {
            int b = (a+1)%3, c = (a+2)%3;
            double S2 = phi[0]*phi[0] + phi[1]*phi[1] + phi[2]*phi[2];
            grad[a] = m2 * phi[a] + lam * (phi[b] + phi[c])
                     + g4 * phi[a] * S2;
            grad[a] -= force_triple(phi[0], phi[1], phi[2], a);
        }
        double gnorm = 0;
        for (int a = 0; a < 3; a++) gnorm += grad[a] * grad[a];
        if (gnorm < 1e-28) break;

        /* Adaptive step: smaller step for larger gradients */
        double step = eta / (1.0 + sqrt(gnorm));
        for (int a = 0; a < 3; a++)
            phi[a] -= step * grad[a];
    }

    for (int a = 0; a < 3; a++) vev_out[a] = phi[a];
    *V_min = pot_density(phi[0], phi[1], phi[2], m2, lam);
}

/* ======================================================================
 * Compute Hessian eigenvalues at the vacuum (mass spectrum of fluctuations)
 * 3x3 matrix: H_{ab} = d^2V / dphi_a dphi_b
 * ====================================================================== */
static void compute_hessian(double lam, double phi[3], double eigenvals[3])
{
    double m2 = mass * mass;
    double H[3][3];

    /* Numerical Hessian */
    double eps = 1e-5;
    for (int a = 0; a < 3; a++) {
        for (int b = a; b < 3; b++) {
            double pp[3], pm[3], mp[3], mm[3];
            for (int c = 0; c < 3; c++) pp[c] = pm[c] = mp[c] = mm[c] = phi[c];
            pp[a] += eps; pp[b] += eps;
            pm[a] += eps; pm[b] -= eps;
            mp[a] -= eps; mp[b] += eps;
            mm[a] -= eps; mm[b] -= eps;
            H[a][b] = (pot_density(pp[0],pp[1],pp[2], m2, lam)
                       - pot_density(pm[0],pm[1],pm[2], m2, lam)
                       - pot_density(mp[0],mp[1],mp[2], m2, lam)
                       + pot_density(mm[0],mm[1],mm[2], m2, lam)) / (4.0*eps*eps);
            H[b][a] = H[a][b];
        }
    }

    /* Print Hessian for debug */
    printf("  Hessian matrix:\n");
    for (int a = 0; a < 3; a++)
        printf("    [%.6f  %.6f  %.6f]\n", H[a][0], H[a][1], H[a][2]);

    /* Eigenvalues of 3x3 symmetric matrix — Cardano's formula */
    double p = -(H[0][0] + H[1][1] + H[2][2]);
    double q = H[0][0]*H[1][1] + H[1][1]*H[2][2] + H[2][2]*H[0][0]
             - H[0][1]*H[0][1] - H[1][2]*H[1][2] - H[2][0]*H[2][0];
    double r = -H[0][0]*H[1][1]*H[2][2] - 2*H[0][1]*H[1][2]*H[2][0]
             + H[0][0]*H[1][2]*H[1][2] + H[1][1]*H[2][0]*H[2][0] + H[2][2]*H[0][1]*H[0][1];

    /* Solve x^3 + px^2 + qx + r = 0 */
    double a2 = q - p*p/3.0;
    double b2 = r - p*q/3.0 + 2.0*p*p*p/27.0;
    double disc = b2*b2/4.0 + a2*a2*a2/27.0;

    if (disc < 0) {
        double R2 = sqrt(-a2*a2*a2/27.0);
        double theta = acos(fmax(-1.0, fmin(1.0, -b2/(2.0*R2))));
        double R = cbrt(R2);
        eigenvals[0] = 2.0*R*cos(theta/3.0) - p/3.0;
        eigenvals[1] = 2.0*R*cos((theta+2*M_PI)/3.0) - p/3.0;
        eigenvals[2] = 2.0*R*cos((theta+4*M_PI)/3.0) - p/3.0;
    } else {
        double sd = sqrt(disc);
        eigenvals[0] = cbrt(-b2/2.0+sd) + cbrt(-b2/2.0-sd) - p/3.0;
        eigenvals[1] = eigenvals[2] = eigenvals[0];
    }

    /* Sort */
    for (int i = 0; i < 2; i++)
        for (int j = i+1; j < 3; j++)
            if (eigenvals[j] < eigenvals[i]) {
                double tmp = eigenvals[i];
                eigenvals[i] = eigenvals[j];
                eigenvals[j] = tmp;
            }
}

/* ======================================================================
 * PHASE 1: Find condensed vacuum for each lambda
 * ====================================================================== */
static void phase1_vacuum(void)
{
    printf("\n=== PHASE 1: Condensed Vacuum Structure ===\n");
    printf("  mu=%.1f kappa=%.1f mass=%.4f g4=%.4f\n", mu, kappa, mass, g4);

    double lam_vals[] = {1.01, 1.05, 1.1, 1.2, 1.5};
    int n_lam = sizeof(lam_vals) / sizeof(lam_vals[0]);

    char path[600];
    snprintf(path, sizeof(path), "%s/phase1_vacuum.tsv", outdir);
    FILE *fp = fopen(path, "w");
    if (!fp) { fprintf(stderr, "Cannot open %s\n", path); return; }
    fprintf(fp, "lambda\tm2_anti\tvev1\tvev2\tvev3\tvev_norm\tV_vac\tV_origin\t"
                "eig0\teig1\teig2\tgoldstone\n");

    for (int k = 0; k < n_lam; k++) {
        double lam = lam_vals[k];
        double m2 = mass * mass;
        double m2_anti = m2 - lam;

        printf("\n--- lambda=%.4f  m2_anti=%.4f (tachyonic=%s) ---\n",
               lam, m2_anti, m2_anti < 0 ? "YES" : "NO");

        double vev[3], V_min;
        find_vacuum(lam, vev, &V_min);

        double V_origin = pot_density(0, 0, 0, m2, lam);

        double vev_norm = sqrt(vev[0]*vev[0] + vev[1]*vev[1] + vev[2]*vev[2]);
        printf("  VEV = (%.6f, %.6f, %.6f)  |VEV|=%.6f\n",
               vev[0], vev[1], vev[2], vev_norm);
        printf("  V(vacuum) = %.6e,  V(origin) = %.6e\n", V_min, V_origin);
        printf("  Sum = %.6e (should be ~0 for antisymmetric)\n",
               vev[0]+vev[1]+vev[2]);

        /* Hessian eigenvalues */
        double eig[3];
        compute_hessian(lam, vev, eig);
        printf("  Hessian eigenvalues: %.6f, %.6f, %.6f\n", eig[0], eig[1], eig[2]);

        int goldstone = (fabs(eig[0]) < 0.01);
        printf("  Goldstone (eigenvalue ~0)? %s (min_eig=%.6e)\n",
               goldstone ? "YES" : "NO", eig[0]);

        /* Analytical prediction for VEV:
         * V_eff = (3/4)(m2-lam)v^2 + (9g4/16)v^4 + O(v^6)
         * Minimum: v^2 = 2(lam-m2)/(3g4)
         * For phi=(v,-v/2,-v/2): |phi|^2 = 3v^2/2 */
        double v2_pred = 2.0 * (lam - m2) / (3.0 * g4);
        double norm_pred = sqrt(1.5 * v2_pred);
        printf("  Analytical |VEV| ~ %.6f (quadratic approx)\n", norm_pred);

        fprintf(fp, "%.4f\t%.4f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6e\t%.6e\t"
                    "%.6f\t%.6f\t%.6f\t%d\n",
                lam, m2_anti, vev[0], vev[1], vev[2], vev_norm, V_min, V_origin,
                eig[0], eig[1], eig[2], goldstone);
        fflush(fp);
    }

    fclose(fp);
    printf("\nOutput: %s\n", path);
}

/* ======================================================================
 * PHASE 2: Fluctuation spectrum around condensate
 * Multi-k excitation: excite many wavelengths, measure dispersion omega(k).
 * ====================================================================== */
static void phase2_spectrum(void)
{
    double m2 = mass * mass;
    double m2_anti = m2 - lambda;

    printf("\n=== PHASE 2: Fluctuation Spectrum at lambda=%.4f ===\n", lambda);
    printf("  m2_anti=%.4f  g4=%.4f\n", m2_anti, g4);

    /* Find vacuum */
    double vev[3], V_min;
    find_vacuum(lambda, vev, &V_min);
    double vev_norm = sqrt(vev[0]*vev[0] + vev[1]*vev[1] + vev[2]*vev[2]);
    printf("  VEV = (%.6f, %.6f, %.6f)  |VEV|=%.6f\n",
           vev[0], vev[1], vev[2], vev_norm);

    /* Hessian for reference */
    double eig[3];
    compute_hessian(lambda, vev, eig);
    printf("  Hessian eigenvalues (=mass^2): %.6f, %.6f, %.6f\n", eig[0], eig[1], eig[2]);

    /* Goldstone direction: perpendicular to VEV in the constraint surface
     * For VEV = (v, -v/2, -v/2), Goldstone ~ (0, 1, -1)/sqrt(2) */
    double gold_dir[3] = {0.0, 1.0/sqrt(2.0), -1.0/sqrt(2.0)};
    /* Radial direction: along vev */
    double rad_dir[3] = {vev[0]/vev_norm, vev[1]/vev_norm, vev[2]/vev_norm};

    /* Grid */
    double dx  = 2.0 * xmax / (Nx - 1);
    double dx2 = dx * dx;
    double kmax_grid = M_PI / dx;
    double m2_max = fabs(eig[2]) + 2.0;
    double dt = 0.4 * 2.0 / sqrt(kmax_grid * kmax_grid + m2_max);

    /* Spectrum measurement time */
    double t_spec = fmin(tfinal, 3000.0);
    int Nt = (int)(t_spec / dt) + 1;

    printf("  Nx=%d xmax=%.1f dx=%.5f dt=%.6f Nt=%d t_spec=%.0f\n",
           Nx, xmax, dx, dt, Nt, t_spec);

    double *phi[3], *vel[3], *acc[3];
    for (int a = 0; a < 3; a++) {
        phi[a] = calloc(Nx, sizeof(double));
        vel[a] = calloc(Nx, sizeof(double));
        acc[a] = calloc(Nx, sizeof(double));
    }

    /* Absorbing boundary */
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

    /* Initialize: condensate + localized perturbation in Goldstone direction.
     * Use a narrow Gaussian to excite many k modes simultaneously. */
    double eps_pert = 0.02;
    double sig_pert = 2.0;
    for (int a = 0; a < 3; a++)
        for (int i = 0; i < Nx; i++) {
            double x = -xmax + i * dx;
            double pert_g = eps_pert * gold_dir[a] * exp(-x*x/(2.0*sig_pert*sig_pert));
            double pert_r = eps_pert * 0.5 * rad_dir[a] * exp(-x*x/(2.0*sig_pert*sig_pert));
            phi[a][i] = vev[a] + pert_g + pert_r;
        }

    #define COMPUTE_ACC2() do { \
        for (int a = 0; a < 3; a++) { \
            acc[a][0] = acc[a][1] = acc[a][Nx-2] = acc[a][Nx-1] = 0; \
            int b = (a+1)%3, c_idx = (a+2)%3; \
            for (int i = 1; i < Nx - 1; i++) { \
                double lapl = (phi[a][i+1] - 2.0*phi[a][i] + phi[a][i-1]) / dx2; \
                double fp = force_triple(phi[0][i], phi[1][i], phi[2][i], a); \
                acc[a][i] = lapl - m2*phi[a][i] - lambda*(phi[b][i]+phi[c_idx][i]) \
                           - g4*phi[a][i]*(phi[0][i]*phi[0][i]+phi[1][i]*phi[1][i]+phi[2][i]*phi[2][i]) + fp; \
            } \
        } \
    } while(0)

    COMPUTE_ACC2();

    /* Probe at multiple distances to measure dispersion */
    int n_probe = 8;
    int probe_idx[8];
    double probe_x[8] = {0.0, 5.0, 10.0, 20.0, 30.0, 50.0, 80.0, 120.0};
    for (int p = 0; p < n_probe; p++)
        probe_idx[p] = (int)((probe_x[p] + xmax) / dx);

    int max_rec = 40000;
    double *gold_hist = calloc(max_rec * n_probe, sizeof(double));
    double *rad_hist  = calloc(max_rec * n_probe, sizeof(double));
    double *t_hist    = calloc(max_rec, sizeof(double));
    int n_rec = 0;
    int rec_every = Nt / max_rec;
    if (rec_every < 1) rec_every = 1;

    int print_every = Nt / 20;
    if (print_every < 1) print_every = 1;

    for (int n = 0; n <= Nt; n++) {
        double t = n * dt;

        if (n % rec_every == 0 && n_rec < max_rec) {
            t_hist[n_rec] = t;
            for (int p = 0; p < n_probe; p++) {
                int idx = probe_idx[p];
                if (idx < 1 || idx >= Nx-1) continue;
                double dphi[3];
                for (int a = 0; a < 3; a++) dphi[a] = phi[a][idx] - vev[a];
                double g = 0, r = 0;
                for (int a = 0; a < 3; a++) {
                    g += dphi[a] * gold_dir[a];
                    r += dphi[a] * rad_dir[a];
                }
                gold_hist[n_rec * n_probe + p] = g;
                rad_hist[n_rec * n_probe + p] = r;
            }
            n_rec++;
        }

        if (n % print_every == 0) {
            double dphi0[3];
            for (int a = 0; a < 3; a++) dphi0[a] = phi[a][Nx/2] - vev[a];
            double g0 = 0, r0 = 0;
            for (int a = 0; a < 3; a++) {
                g0 += dphi0[a] * gold_dir[a];
                r0 += dphi0[a] * rad_dir[a];
            }
            printf("  t=%7.1f  gold(0)=%+.4e  rad(0)=%+.4e\n", t, g0, r0);
        }

        if (n == Nt) break;

        /* Velocity Verlet */
        for (int a = 0; a < 3; a++)
            for (int i = 1; i < Nx - 1; i++)
                vel[a][i] += 0.5 * dt * acc[a][i];
        for (int a = 0; a < 3; a++)
            for (int i = 1; i < Nx - 1; i++)
                phi[a][i] += dt * vel[a][i];
        COMPUTE_ACC2();
        for (int a = 0; a < 3; a++)
            for (int i = 1; i < Nx - 1; i++)
                vel[a][i] += 0.5 * dt * acc[a][i];

        /* Absorbing boundary: damp only the fluctuation */
        for (int a = 0; a < 3; a++)
            for (int i = 0; i < Nx; i++) {
                vel[a][i] *= damp[i];
                phi[a][i] = vev[a] + (phi[a][i] - vev[a]) * damp[i];
            }
    }
    #undef COMPUTE_ACC2

    /* DFT analysis */
    printf("\n--- Spectrum analysis ---\n");

    char specpath[600];
    snprintf(specpath, sizeof(specpath), "%s/phase2_spectrum_lam%.2f.tsv",
             outdir, lambda);
    FILE *fspec = fopen(specpath, "w");
    fprintf(fspec, "omega");
    for (int p = 0; p < n_probe; p++)
        fprintf(fspec, "\tgold_x%.0f\trad_x%.0f", probe_x[p], probe_x[p]);
    fprintf(fspec, "\n");

    int dft_start = n_rec / 4;
    double T = t_hist[n_rec-1] - t_hist[dft_start];
    int nf = 1000;
    double dom = 5.0 / nf;

    double gold_peak_om[8] = {0}, gold_peak_pw[8] = {0};
    double rad_peak_om[8] = {0}, rad_peak_pw[8] = {0};

    for (int f = 0; f < nf; f++) {
        double omega = f * dom;
        fprintf(fspec, "%.6f", omega);
        for (int p = 0; p < n_probe; p++) {
            double re_g = 0, im_g = 0, re_r = 0, im_r = 0;
            for (int j = dft_start; j < n_rec; j++) {
                double dtj = (j > dft_start) ?
                    (t_hist[j] - t_hist[j-1]) : (t_hist[dft_start+1] - t_hist[dft_start]);
                double gv = gold_hist[j * n_probe + p];
                double rv = rad_hist[j * n_probe + p];
                double c = cos(omega * t_hist[j]);
                double s = sin(omega * t_hist[j]);
                re_g += gv * c * dtj;
                im_g += gv * s * dtj;
                re_r += rv * c * dtj;
                im_r += rv * s * dtj;
            }
            double pw_g = (re_g*re_g + im_g*im_g) / (T*T);
            double pw_r = (re_r*re_r + im_r*im_r) / (T*T);
            fprintf(fspec, "\t%.6e\t%.6e", pw_g, pw_r);
            if (pw_g > gold_peak_pw[p] && omega > 0.02) {
                gold_peak_pw[p] = pw_g;
                gold_peak_om[p] = omega;
            }
            if (pw_r > rad_peak_pw[p] && omega > 0.02) {
                rad_peak_pw[p] = pw_r;
                rad_peak_om[p] = omega;
            }
        }
        fprintf(fspec, "\n");
    }
    fclose(fspec);

    printf("  Probe results:\n");
    printf("  %-10s  %-14s  %-14s\n", "x", "gold_omega", "rad_omega");
    for (int p = 0; p < n_probe; p++)
        printf("  %-10.1f  %-14.6f  %-14.6f\n",
               probe_x[p], gold_peak_om[p], rad_peak_om[p]);

    /* Arrival time analysis: measure when wavefront arrives at each probe */
    printf("\n  Wavefront arrival (Goldstone):\n");
    printf("  %-10s  %-14s  %-14s\n", "x", "t_arrive", "v_group");
    double t_prev = 0, x_prev = 0;
    for (int p = 0; p < n_probe; p++) {
        double thresh = 1e-5;
        double t_arr = -1;
        for (int j = 0; j < n_rec; j++) {
            if (fabs(gold_hist[j * n_probe + p]) > thresh) {
                t_arr = t_hist[j];
                break;
            }
        }
        double v_grp = 0;
        if (t_arr > 0 && p > 0 && t_arr > t_prev) {
            v_grp = (probe_x[p] - x_prev) / (t_arr - t_prev);
        }
        printf("  %-10.1f  %-14.4f  %-14.4f\n", probe_x[p], t_arr, v_grp);
        if (t_arr > 0) { t_prev = t_arr; x_prev = probe_x[p]; }
    }

    printf("\n  Hessian mass^2 eigenvalues: %.6f, %.6f, %.6f\n", eig[0], eig[1], eig[2]);
    printf("  If Goldstone: omega(k->0) -> 0\n");
    printf("  Predicted Goldstone mass = %.6f (should be ~0)\n", sqrt(fabs(eig[0])));
    printf("  Predicted radial mass = %.6f\n", sqrt(fabs(eig[1])));
    printf("  Predicted symmetric mass = %.6f\n", sqrt(fabs(eig[2])));

    printf("\nOutput: %s\n", specpath);

    free(gold_hist); free(rad_hist); free(t_hist);
    for (int a = 0; a < 3; a++) { free(phi[a]); free(vel[a]); free(acc[a]); }
    free(damp);
}

/* ======================================================================
 * PHASE 3: Oscillon on top of condensate
 * ====================================================================== */
static void phase3_oscillon(void)
{
    double m2 = mass * mass;

    printf("\n=== PHASE 3: Oscillon in Condensed Vacuum at lambda=%.4f ===\n", lambda);

    /* Find vacuum */
    double vev[3], V_min;
    find_vacuum(lambda, vev, &V_min);
    double vev_norm = sqrt(vev[0]*vev[0] + vev[1]*vev[1] + vev[2]*vev[2]);
    printf("  VEV = (%.6f, %.6f, %.6f)  |VEV|=%.6f\n",
           vev[0], vev[1], vev[2], vev_norm);

    /* Mass spectrum around condensate */
    double eig[3];
    compute_hessian(lambda, vev, eig);
    printf("  Hessian eigenvalues: %.6f, %.6f, %.6f\n", eig[0], eig[1], eig[2]);

    /* Grid */
    double dx  = 2.0 * xmax / (Nx - 1);
    double dx2 = dx * dx;
    double kmax_grid = M_PI / dx;
    double m2_max = fabs(eig[2]) + 2.0;
    double dt = 0.4 * 2.0 / sqrt(kmax_grid * kmax_grid + m2_max);
    int Nt = (int)(tfinal / dt) + 1;

    printf("  Nx=%d xmax=%.1f dx=%.5f dt=%.6f Nt=%d tfinal=%.0f\n",
           Nx, xmax, dx, dt, Nt, tfinal);

    double *phi[3], *vel[3], *acc[3];
    for (int a = 0; a < 3; a++) {
        phi[a] = calloc(Nx, sizeof(double));
        vel[a] = calloc(Nx, sizeof(double));
        acc[a] = calloc(Nx, sizeof(double));
    }

    /* Absorbing boundary */
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

    /* Initialize: condensate + symmetric Gaussian oscillon */
    for (int a = 0; a < 3; a++)
        for (int i = 0; i < Nx; i++) {
            double x = -xmax + i * dx;
            double osc = A_init * exp(-x*x / (2.0 * sigma_w * sigma_w));
            phi[a][i] = vev[a] + osc;
        }

    #define COMPUTE_ACC3() do { \
        for (int a = 0; a < 3; a++) { \
            acc[a][0] = acc[a][1] = acc[a][Nx-2] = acc[a][Nx-1] = 0; \
            int b = (a+1)%3, c_idx = (a+2)%3; \
            for (int i = 1; i < Nx - 1; i++) { \
                double lapl = (phi[a][i+1] - 2.0*phi[a][i] + phi[a][i-1]) / dx2; \
                double fp = force_triple(phi[0][i], phi[1][i], phi[2][i], a); \
                acc[a][i] = lapl - m2*phi[a][i] - lambda*(phi[b][i]+phi[c_idx][i]) \
                           - g4*phi[a][i]*(phi[0][i]*phi[0][i]+phi[1][i]*phi[1][i]+phi[2][i]*phi[2][i]) + fp; \
            } \
        } \
    } while(0)

    COMPUTE_ACC3();

    /* Time series output */
    char tspath[600];
    snprintf(tspath, sizeof(tspath), "%s/phase3_oscillon_lam%.2f.tsv", outdir, lambda);
    FILE *fts = fopen(tspath, "w");
    fprintf(fts, "time\tphi1_0\tphi2_0\tphi3_0\tsym_0\tgold_0\t"
                 "peak_sym\tE_core\tE_total\tfc\n");

    /* DFT storage */
    int max_dft = 50000;
    double *sym_hist = calloc(max_dft, sizeof(double));
    double *gold_hist_c = calloc(max_dft, sizeof(double));
    double *t_hist = calloc(max_dft, sizeof(double));
    int n_dft = 0;

    int rec_every = Nt / max_dft;
    if (rec_every < 1) rec_every = 1;
    int print_every = Nt / 40;
    if (print_every < 1) print_every = 1;
    double core_r = 4.0 * sigma_w;
    int ic = Nx / 2;

    double gold_dir[3] = {0.0, 1.0/sqrt(2.0), -1.0/sqrt(2.0)};
    double sym_dir[3] = {1.0/sqrt(3.0), 1.0/sqrt(3.0), 1.0/sqrt(3.0)};

    for (int n = 0; n <= Nt; n++) {
        double t = n * dt;

        if (n % rec_every == 0 && n_dft < max_dft) {
            double dphi[3];
            for (int a = 0; a < 3; a++) dphi[a] = phi[a][ic] - vev[a];
            double sv = 0, gv = 0;
            for (int a = 0; a < 3; a++) {
                sv += dphi[a] * sym_dir[a];
                gv += dphi[a] * gold_dir[a];
            }
            sym_hist[n_dft] = sv;
            gold_hist_c[n_dft] = gv;
            t_hist[n_dft] = t;
            n_dft++;
        }

        int do_rec = (n % rec_every == 0);
        int do_print = (n % print_every == 0);

        if (do_rec || do_print) {
            double Ecore = 0, Eall = 0;
            double peak_sym = 0;

            for (int i = 1; i < Nx - 1; i++) {
                double x = -xmax + i * dx;
                double e = 0;
                for (int a = 0; a < 3; a++) {
                    e += 0.5 * vel[a][i] * vel[a][i];
                    double dp = (phi[a][i+1] - phi[a][i-1]) / (2.0*dx);
                    e += 0.5 * dp * dp + 0.5 * m2 * phi[a][i] * phi[a][i];
                }
                e += lambda * (phi[0][i]*phi[1][i] + phi[1][i]*phi[2][i]
                              + phi[2][i]*phi[0][i]);
                {
                    double S2i = phi[0][i]*phi[0][i]+phi[1][i]*phi[1][i]+phi[2][i]*phi[2][i];
                    e += (g4/4.0) * S2i * S2i;
                }
                double P = phi[0][i] * phi[1][i] * phi[2][i];
                e += 0.5 * mu * P * P / (1.0 + kappa * P * P);
                e -= V_min;

                Eall += e * dx;
                if (fabs(x) < core_r) Ecore += e * dx;

                double dphi[3];
                for (int a = 0; a < 3; a++) dphi[a] = phi[a][i] - vev[a];
                double sv = 0;
                for (int a = 0; a < 3; a++) sv += dphi[a] * sym_dir[a];
                if (fabs(sv) > peak_sym) peak_sym = fabs(sv);
            }
            double fc = (fabs(Eall) > 1e-20) ? Ecore / Eall : 0.0;

            double dphi0[3];
            for (int a = 0; a < 3; a++) dphi0[a] = phi[a][ic] - vev[a];
            double s0 = 0, g0 = 0;
            for (int a = 0; a < 3; a++) {
                s0 += dphi0[a] * sym_dir[a];
                g0 += dphi0[a] * gold_dir[a];
            }

            if (do_rec)
                fprintf(fts, "%.6f\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t"
                             "%.6e\t%.6e\t%.6e\t%.6f\n",
                        t, phi[0][ic], phi[1][ic], phi[2][ic],
                        s0, g0, peak_sym, Ecore, Eall, fc);

            if (do_print)
                printf("  t=%7.1f  sym=%+.4e  gold=%+.4e  pk_s=%.4e  fc=%.4f  E=%.4f\n",
                       t, s0, g0, peak_sym, fc, Eall);
        }

        if (n == Nt) break;

        /* Velocity Verlet */
        for (int a = 0; a < 3; a++)
            for (int i = 1; i < Nx - 1; i++)
                vel[a][i] += 0.5 * dt * acc[a][i];
        for (int a = 0; a < 3; a++)
            for (int i = 1; i < Nx - 1; i++)
                phi[a][i] += dt * vel[a][i];
        COMPUTE_ACC3();
        for (int a = 0; a < 3; a++)
            for (int i = 1; i < Nx - 1; i++)
                vel[a][i] += 0.5 * dt * acc[a][i];

        /* Absorbing boundary: damp fluctuations around condensate */
        for (int a = 0; a < 3; a++)
            for (int i = 0; i < Nx; i++) {
                vel[a][i] *= damp[i];
                phi[a][i] = vev[a] + (phi[a][i] - vev[a]) * damp[i];
            }
    }
    #undef COMPUTE_ACC3

    fclose(fts);

    /* DFT of symmetric and Goldstone modes at center */
    printf("\n--- Oscillon spectrum analysis ---\n");
    char dftpath[600];
    snprintf(dftpath, sizeof(dftpath), "%s/phase3_spectrum_lam%.2f.tsv", outdir, lambda);
    FILE *fdft = fopen(dftpath, "w");
    fprintf(fdft, "omega\tsym_power\tgold_power\n");

    int dft_start = n_dft / 2;
    double T_dft = t_hist[n_dft-1] - t_hist[dft_start];
    int nf = 500;
    double sym_peak_om = 0, sym_peak_pw = 0;
    double gold_peak_om_c = 0, gold_peak_pw_c = 0;

    for (int f = 0; f < nf; f++) {
        double omega = 8.0 * f / nf;
        double re_s = 0, im_s = 0, re_g = 0, im_g = 0;
        for (int j = dft_start; j < n_dft; j++) {
            double dtj = (j > dft_start) ?
                (t_hist[j] - t_hist[j-1]) : (t_hist[dft_start+1] - t_hist[dft_start]);
            double c = cos(omega * t_hist[j]);
            double s = sin(omega * t_hist[j]);
            re_s += sym_hist[j] * c * dtj;
            im_s += sym_hist[j] * s * dtj;
            re_g += gold_hist_c[j] * c * dtj;
            im_g += gold_hist_c[j] * s * dtj;
        }
        double pw_s = (re_s*re_s + im_s*im_s) / (T_dft*T_dft);
        double pw_g = (re_g*re_g + im_g*im_g) / (T_dft*T_dft);
        fprintf(fdft, "%.6f\t%.6e\t%.6e\n", omega, pw_s, pw_g);
        if (pw_s > sym_peak_pw && omega > 0.02) {
            sym_peak_pw = pw_s; sym_peak_om = omega;
        }
        if (pw_g > gold_peak_pw_c && omega > 0.02) {
            gold_peak_pw_c = pw_g; gold_peak_om_c = omega;
        }
    }
    fclose(fdft);

    double m_eff_sym = sqrt(fabs(eig[2]));
    printf("  Symmetric mode peak: omega=%.4f (mass=%.4f)\n", sym_peak_om, m_eff_sym);
    printf("  Goldstone mode peak: omega=%.4f\n", gold_peak_om_c);
    printf("  Oscillon (sym_omega < m_sym_eff)? %s\n",
           (sym_peak_om > 0.01 && sym_peak_om < m_eff_sym) ? "YES" : "NO");

    printf("\nOutput: %s\n", tspath);
    printf("Output: %s\n", dftpath);

    free(sym_hist); free(gold_hist_c); free(t_hist);
    for (int a = 0; a < 3; a++) { free(phi[a]); free(vel[a]); free(acc[a]); }
    free(damp);
}

/* ====================================================================== */
int main(int argc, char **argv)
{
    parse_args(argc, argv);

    printf("condensed1d: lambda=%.4f mu=%.1f kappa=%.1f mass=%.4f g4=%.4f\n",
           lambda, mu, kappa, mass, g4);
    printf("  Nx=%d xmax=%.1f tfinal=%.0f\n", Nx, xmax, tfinal);

    if (phase == 0 || phase == 1) phase1_vacuum();
    if (phase == 0 || phase == 2) phase2_spectrum();
    if (phase == 0 || phase == 3) phase3_oscillon();

    return 0;
}
