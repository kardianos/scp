/*  blv_cosserat.c — BLV effective metric for 6-field Cosserat braid
 *
 *  Level 1: scalar perturbation (phi only, no theta)
 *  Level 2: coupled phi-theta perturbation dispersion
 *
 *  Reads the measured radial energy density profile from the V34 phonon test,
 *  reconstructs the background field configuration, and computes:
 *    - Effective mass m_eff^2(r) from V(P) coupling
 *    - Group velocity v_g(r) for test waves
 *    - Gravitational potential analog Phi(r)/c^2
 *    - Deflection angle theta(b) via Born integral
 *
 *  Build: gcc -O3 -o blv_cosserat src/blv_cosserat.c -lm
 *  Run:   ./blv_cosserat [profile.tsv]
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define PI 3.14159265358979323846

/* Physics parameters */
static const double MU     = -41.345;
static const double KAPPA  = 50.0;
static const double MASS2  = 2.25;    /* phi mass^2 */
static const double MTHETA2 = 0.0;   /* theta mass^2 (massless) */
static const double ETA    = 0.5;    /* phi-theta coupling */
static const double A_BG   = 0.1;   /* background amplitude */

/* Braid parameters */
static const double DELTA[3] = {0.0, 3.0005, 4.4325};  /* phase offsets */
static const double A_CORE = 0.80;  /* field amplitude at braid core */

/* ================================================================
   V(P) and derivatives
   ================================================================ */

static double V_func(double P) {
    double P2 = P * P;
    return (MU / 2.0) * P2 / (1.0 + KAPPA * P2);
}

static double Vp_func(double P) {
    /* dV/dP */
    double P2 = P * P;
    double den = 1.0 + KAPPA * P2;
    return MU * P / (den * den);
}

static double Vpp_func(double P) {
    /* d^2V/dP^2 */
    double P2 = P * P;
    double den = 1.0 + KAPPA * P2;
    return MU * (1.0 - 3.0 * KAPPA * P2) / (den * den * den);
}

/* ================================================================
   Time-averaged background quantities at radius r

   The braid fields oscillate as:
     phi_a(theta) = A(r) * cos(theta + delta_a)
   where theta = kz - omega*t sweeps 0..2pi.

   We numerically average all V(P)-derived quantities over one period.
   ================================================================ */

typedef struct {
    double A;         /* field amplitude at this r */
    double rho_phi2;  /* <sum phi_a^2> = 3A^2/2 */
    double P_rms;     /* sqrt(<P^2>) */
    double V_avg;     /* <V(P)> */
    double W_diag;    /* <V''(P)*(dP/dphi_a)^2 + V'(P)*d^2P/dphi_a^2>, averaged over a */
    double W_offdiag; /* off-diagonal W_ab average (for cross-coupling) */
} BgAvg;

static BgAvg compute_bg_avg(double A) {
    BgAvg bg = {0};
    bg.A = A;
    bg.rho_phi2 = 1.5 * A * A;

    int N = 1000;  /* integration points */
    double dth = 2.0 * PI / N;

    double sum_P2 = 0, sum_V = 0;
    double sum_W00 = 0, sum_W11 = 0, sum_W22 = 0;
    double sum_W01 = 0;

    for (int i = 0; i < N; i++) {
        double th = (i + 0.5) * dth;

        double p0 = A * cos(th + DELTA[0]);
        double p1 = A * cos(th + DELTA[1]);
        double p2 = A * cos(th + DELTA[2]);

        double P = p0 * p1 * p2;
        double P2 = P * P;

        sum_P2 += P2;
        sum_V += V_func(P);

        double Vp = Vp_func(P);
        double Vpp = Vpp_func(P);

        /* dP/dphi_a = product of other two */
        double dPd0 = p1 * p2;
        double dPd1 = p0 * p2;
        double dPd2 = p0 * p1;

        /* d^2P/(dphi_a dphi_b) for a != b: the third field */
        double d2P_01 = p2;
        double d2P_02 = p1;
        double d2P_12 = p0;
        /* d^2P/(dphi_a)^2 = 0 for triple product */

        /* W_aa = V''(P) * (dP/dphi_a)^2 + V'(P) * d^2P/dphi_a^2
                = V''(P) * (dP/dphi_a)^2 + 0 */
        sum_W00 += Vpp * dPd0 * dPd0;
        sum_W11 += Vpp * dPd1 * dPd1;
        sum_W22 += Vpp * dPd2 * dPd2;

        /* W_01 = V''(P) * dPd0 * dPd1 + V'(P) * d2P_01 */
        sum_W01 += Vpp * dPd0 * dPd1 + Vp * d2P_01;
    }

    bg.P_rms = sqrt(sum_P2 / N);
    bg.V_avg = sum_V / N;
    /* Average diagonal element over the three components */
    bg.W_diag = (sum_W00 + sum_W11 + sum_W22) / (3.0 * N);
    bg.W_offdiag = sum_W01 / N;

    return bg;
}

/* ================================================================
   Read depletion profile
   ================================================================ */

#define MAX_PROFILE 200

typedef struct {
    double r[MAX_PROFILE];
    double rho[MAX_PROFILE];
    double delta_rho[MAX_PROFILE];
    int n;
    double rho_bg;
} Profile;

static Profile read_profile(const char *path) {
    Profile p = {0};
    p.rho_bg = 0.01794307;  /* default from phonon test header */

    FILE *fp = fopen(path, "r");
    if (!fp) {
        fprintf(stderr, "Cannot open %s, using analytic model\n", path);
        /* Generate analytic power-law profile */
        p.n = 80;
        for (int i = 0; i < p.n; i++) {
            p.r[i] = 0.5 + i * 0.5;
            if (p.r[i] < 4.0) {
                p.rho[i] = 0.08;  /* flat core */
            } else {
                p.rho[i] = p.rho_bg + 0.35 / pow(p.r[i], 1.2);
            }
            p.delta_rho[i] = p.rho[i] - p.rho_bg;
        }
        return p;
    }

    char line[512];
    while (fgets(line, sizeof(line), fp)) {
        if (line[0] == '#' || line[0] == 'r') continue;
        if (p.n >= MAX_PROFILE) break;
        double r, rho, dr;
        int cnt;
        if (sscanf(line, "%lf %lf %lf %d", &r, &rho, &dr, &cnt) >= 3) {
            p.r[p.n] = r;
            p.rho[p.n] = rho;
            p.delta_rho[p.n] = dr;
            p.n++;
        }
    }
    fclose(fp);
    printf("Loaded profile: %d points, rho_bg=%.6e\n", p.n, p.rho_bg);
    return p;
}

/* Linear interpolation on profile */
static double interp_rho(const Profile *p, double r) {
    if (r <= p->r[0]) return p->rho[0];
    if (r >= p->r[p->n-1]) return p->rho_bg;
    for (int i = 0; i < p->n - 1; i++) {
        if (r >= p->r[i] && r < p->r[i+1]) {
            double t = (r - p->r[i]) / (p->r[i+1] - p->r[i]);
            return p->rho[i] + t * (p->rho[i+1] - p->rho[i]);
        }
    }
    return p->rho_bg;
}

/* ================================================================
   Convert energy density rho to field amplitude A

   rho = (3/4)*A^2*(omega^2 + k^2 + m^2) + <V(P)>
       ≈ (3/2)*A^2*m^2  (dominant term for k << m)

   More precisely, we calibrate: at the core, A_core=0.8, rho_core~0.08
   So A(r) = A_core * sqrt(rho(r) / rho_core)
   ================================================================ */

static double rho_to_amplitude(double rho, double rho_core) {
    if (rho <= 0) return 0;
    return A_CORE * sqrt(rho / rho_core);
}

/* ================================================================
   Level 1: Scalar BLV metric (phi only)

   Effective mass: m_eff^2(r) = m^2 + W_diag(r)
   Group velocity: v_g = k/omega = 1/sqrt(1 + m_eff^2/k^2)
   Potential: Phi(r)/c^2 = -(m_eff^2(r) - m_eff^2(inf)) / (2*k_ref^2)
   ================================================================ */

typedef struct {
    double r;
    double A;
    double P_rms;
    double W_diag;
    double m_eff2;          /* m^2 + W_diag */
    double v_g;             /* group velocity at k_ref */
    double phi_over_c2;     /* gravitational potential */
    /* Level 2 additions */
    double m_eff2_coupled;  /* effective mass with theta coupling */
    double v_g_coupled;
    double phi_coupled;
} MetricPoint;

/* ================================================================
   Level 2: Coupled phi-theta dispersion

   For z-propagation, the coupled dispersion for transverse modes is:
     (alpha - m^2 - W)(alpha - m_theta^2) = eta^2 * k^2
   where alpha = omega^2 - k^2.

   The phi-dominant solution (alpha > 0):
     alpha = [(m^2+W+mt^2) + sqrt((m^2+W-mt^2)^2 + 4*eta^2*k^2)] / 2
   ================================================================ */

static double coupled_alpha(double W, double k2) {
    double mW = MASS2 + W;
    double diff = mW - MTHETA2;
    double sum = mW + MTHETA2;
    double disc = diff * diff + 4.0 * ETA * ETA * k2;
    /* phi-dominant branch (larger alpha) */
    return 0.5 * (sum + sqrt(disc));
}

/* ================================================================
   Deflection angle via Born integral

   theta(b) = -b * integral_b^inf (dn/dr) / sqrt(r^2 - b^2) dr
   where n(r) = 1/v_g(r) is the refractive index.

   Use cosh substitution: r = b*cosh(t), dr = b*sinh(t)*dt
   sqrt(r^2-b^2) = b*sinh(t), so integrand simplifies.
   ================================================================ */

static double deflection_angle(MetricPoint *pts, int n, double b, int use_coupled) {
    double theta = 0;
    int nt = 2000;
    double t_max = 5.0;  /* cosh(5) ~ 74, so integrates to r ~ 74*b */
    double dt = t_max / nt;

    for (int i = 0; i < nt; i++) {
        double t = (i + 0.5) * dt;
        double r = b * cosh(t);

        /* Find n(r) and dn/dr by finite difference */
        double dr_fd = 0.05;
        double r_p = r + dr_fd, r_m = r - dr_fd;
        if (r_m < 0.3) r_m = 0.3;

        /* Interpolate metric at r, r+dr, r-dr */
        double n_r = 0, n_p = 0, n_m = 0;
        /* Simple linear search (profile is small) */
        for (int j = 0; j < n - 1; j++) {
            if (pts[j].r <= r && r < pts[j+1].r) {
                double f = (r - pts[j].r) / (pts[j+1].r - pts[j].r);
                n_r = use_coupled ?
                    1.0 / (pts[j].v_g_coupled + f*(pts[j+1].v_g_coupled - pts[j].v_g_coupled)) :
                    1.0 / (pts[j].v_g + f*(pts[j+1].v_g - pts[j].v_g));
                break;
            }
        }
        for (int j = 0; j < n - 1; j++) {
            if (pts[j].r <= r_p && r_p < pts[j+1].r) {
                double f = (r_p - pts[j].r) / (pts[j+1].r - pts[j].r);
                n_p = use_coupled ?
                    1.0 / (pts[j].v_g_coupled + f*(pts[j+1].v_g_coupled - pts[j].v_g_coupled)) :
                    1.0 / (pts[j].v_g + f*(pts[j+1].v_g - pts[j].v_g));
                break;
            }
        }
        for (int j = 0; j < n - 1; j++) {
            if (pts[j].r <= r_m && r_m < pts[j+1].r) {
                double f = (r_m - pts[j].r) / (pts[j+1].r - pts[j].r);
                n_m = use_coupled ?
                    1.0 / (pts[j].v_g_coupled + f*(pts[j+1].v_g_coupled - pts[j].v_g_coupled)) :
                    1.0 / (pts[j].v_g + f*(pts[j+1].v_g - pts[j].v_g));
                break;
            }
        }

        if (n_r == 0 || n_p == 0 || n_m == 0) continue;

        double dn_dr = (n_p - n_m) / (r_p - r_m);
        /* In cosh coords: integrand = dn/dr * b*sinh(t)*dt / (b*sinh(t)) = dn/dr * dt */
        theta += dn_dr * dt;
    }
    return -2.0 * theta;  /* factor 2 from symmetry */
}

/* ================================================================
   Main
   ================================================================ */

int main(int argc, char **argv) {
    const char *profile_path = "../v34/phonon_test/data/depletion_t0100.tsv";
    if (argc > 1) profile_path = argv[1];

    Profile prof = read_profile(profile_path);

    /* Reference wavenumber for dispersion evaluation */
    double k_ref = sqrt(MASS2);  /* k = m: Compton wavelength */
    double k2_ref = k_ref * k_ref;

    /* Core rho (for amplitude scaling) */
    double rho_core = prof.rho[0];
    for (int i = 0; i < prof.n && i < 10; i++)
        if (prof.rho[i] > rho_core) rho_core = prof.rho[i];

    printf("=== BLV Effective Metric for 6-field Cosserat Braid ===\n\n");
    printf("Parameters: m^2=%.4f, m_theta^2=%.4f, eta=%.3f\n", MASS2, MTHETA2, ETA);
    printf("            mu=%.3f, kappa=%.1f\n", MU, KAPPA);
    printf("Phase offsets: delta = {%.4f, %.4f, %.4f}\n", DELTA[0], DELTA[1], DELTA[2]);
    printf("Reference k = sqrt(m^2) = %.4f  (Compton wavelength)\n", k_ref);
    printf("Core rho = %.6f, A_core = %.3f\n", rho_core, A_CORE);
    printf("Background rho = %.6f, A_bg = %.3f\n\n", prof.rho_bg, A_BG);

    /* Compute background averages at infinity (vacuum) */
    BgAvg bg_inf = compute_bg_avg(A_BG);
    double m_eff2_inf = MASS2 + bg_inf.W_diag;
    double alpha_inf = coupled_alpha(bg_inf.W_diag, k2_ref);

    printf("Background (r->inf):\n");
    printf("  A=%.4f, <P^2>^(1/2)=%.6e, <V>=%.6e\n", bg_inf.A, bg_inf.P_rms, bg_inf.V_avg);
    printf("  W_diag=%.6e, m_eff^2=%.6f\n", bg_inf.W_diag, m_eff2_inf);
    printf("  v_g(inf) = %.6f (at k=m)\n", k_ref / sqrt(k2_ref + m_eff2_inf));
    printf("  alpha_coupled(inf) = %.6f\n\n", alpha_inf);

    /* Compute metric at each profile radius */
    int n_pts = prof.n;
    MetricPoint *pts = calloc(n_pts, sizeof(MetricPoint));

    printf("%-8s %-8s %-10s %-12s %-12s %-12s | %-12s %-12s %-12s\n",
           "r", "A(r)", "P_rms", "W_diag", "Phi_L1/c^2", "v_g_L1",
           "W_coupled", "Phi_L2/c^2", "v_g_L2");
    printf("%-8s %-8s %-10s %-12s %-12s %-12s | %-12s %-12s %-12s\n",
           "----", "----", "------", "------", "----------", "------",
           "---------", "----------", "------");

    for (int i = 0; i < n_pts; i++) {
        double r = prof.r[i];
        double rho = prof.rho[i];
        double A = rho_to_amplitude(rho, rho_core);

        BgAvg bg = compute_bg_avg(A);

        /* Level 1: scalar m_eff */
        double m_eff2 = MASS2 + bg.W_diag;
        double omega2_L1 = k2_ref + m_eff2;
        double v_g_L1 = (omega2_L1 > 0) ? k_ref / sqrt(omega2_L1) : 0;
        double phi_L1 = -(m_eff2 - m_eff2_inf) / (2.0 * k2_ref);

        /* Level 2: coupled dispersion */
        double alpha = coupled_alpha(bg.W_diag, k2_ref);
        double omega2_L2 = k2_ref + alpha;
        double v_g_L2 = (omega2_L2 > 0) ? k_ref / sqrt(omega2_L2) : 0;
        double phi_L2 = -(alpha - alpha_inf) / (2.0 * k2_ref);

        pts[i].r = r;
        pts[i].A = A;
        pts[i].P_rms = bg.P_rms;
        pts[i].W_diag = bg.W_diag;
        pts[i].m_eff2 = m_eff2;
        pts[i].v_g = v_g_L1;
        pts[i].phi_over_c2 = phi_L1;
        pts[i].m_eff2_coupled = alpha;
        pts[i].v_g_coupled = v_g_L2;
        pts[i].phi_coupled = phi_L2;

        if (i % 2 == 0 || r <= 5.0) {
            printf("%-8.2f %-8.4f %-10.4e %-12.4e %-12.6f %-12.6f | %-12.6f %-12.6f %-12.6f\n",
                   r, A, bg.P_rms, bg.W_diag, phi_L1, v_g_L1,
                   alpha, phi_L2, v_g_L2);
        }
    }

    /* ============================================================
       Summary of potential well
       ============================================================ */

    double phi_min_L1 = 0, phi_min_L2 = 0;
    double r_min_L1 = 0, r_min_L2 = 0;
    for (int i = 0; i < n_pts; i++) {
        if (pts[i].phi_over_c2 < phi_min_L1) {
            phi_min_L1 = pts[i].phi_over_c2;
            r_min_L1 = pts[i].r;
        }
        if (pts[i].phi_coupled < phi_min_L2) {
            phi_min_L2 = pts[i].phi_coupled;
            r_min_L2 = pts[i].r;
        }
    }

    printf("\n=== POTENTIAL WELL SUMMARY ===\n\n");
    printf("Level 1 (scalar, no theta):\n");
    printf("  Phi_min/c^2 = %.6f  at r = %.2f\n", phi_min_L1, r_min_L1);
    printf("  Potential depth = %.2f MeV  (1 code energy = 9.098 MeV)\n",
           fabs(phi_min_L1) * 103.13 * 9.098);  /* E_sol * MeV/code */
    printf("  v_g_min/c = %.6f  (%.2f%% slower than vacuum)\n",
           pts[0].v_g, 100.0 * (1.0 - pts[0].v_g / pts[n_pts-1].v_g));
    printf("  delta_c/c at core = %.6f\n",
           (pts[0].v_g - pts[n_pts-1].v_g) / pts[n_pts-1].v_g);

    printf("\nLevel 2 (coupled phi-theta, eta=%.2f):\n", ETA);
    printf("  Phi_min/c^2 = %.6f  at r = %.2f\n", phi_min_L2, r_min_L2);
    printf("  Potential depth = %.2f MeV\n",
           fabs(phi_min_L2) * 103.13 * 9.098);
    printf("  v_g_min/c = %.6f  (%.2f%% slower than vacuum)\n",
           pts[0].v_g_coupled,
           100.0 * (1.0 - pts[0].v_g_coupled / pts[n_pts-1].v_g_coupled));
    printf("  delta_c/c at core = %.6f\n",
           (pts[0].v_g_coupled - pts[n_pts-1].v_g_coupled) / pts[n_pts-1].v_g_coupled);

    double ratio = (phi_min_L1 != 0) ? phi_min_L2 / phi_min_L1 : 0;
    printf("\n  Theta coupling effect: Phi_L2/Phi_L1 = %.4f (%.1f%% %s)\n",
           ratio, 100.0 * fabs(1.0 - ratio),
           (fabs(phi_min_L2) > fabs(phi_min_L1)) ? "deeper" : "shallower");

    /* ============================================================
       Deflection angles
       ============================================================ */

    printf("\n=== DEFLECTION ANGLES (Born approximation) ===\n\n");
    printf("%-8s %-14s %-14s %-14s\n",
           "b", "theta_L1 (rad)", "theta_L2 (rad)", "L2/L1 ratio");
    printf("%-8s %-14s %-14s %-14s\n",
           "----", "-----------", "-----------", "-----------");

    double b_values[] = {1.0, 2.0, 3.0, 5.0, 8.0, 10.0, 15.0, 20.0, 30.0};
    int n_b = sizeof(b_values) / sizeof(b_values[0]);

    for (int ib = 0; ib < n_b; ib++) {
        double b = b_values[ib];
        double th_L1 = deflection_angle(pts, n_pts, b, 0);
        double th_L2 = deflection_angle(pts, n_pts, b, 1);
        double ratio_b = (th_L1 != 0) ? th_L2 / th_L1 : 0;
        printf("%-8.1f %-14.6e %-14.6e %-14.4f\n", b, th_L1, th_L2, ratio_b);
    }

    /* ============================================================
       Physical comparison
       ============================================================ */

    printf("\n=== PHYSICAL COMPARISON ===\n\n");
    /* Convert to physical units */
    double fm_per_code = 0.5624;
    double MeV_per_code = 9.098;
    double E_sol = 103.13;  /* braid energy in code units */
    double M_sol_MeV = E_sol * MeV_per_code;

    printf("Braid mass: %.1f MeV (= %.3f GeV, proton-scale)\n", M_sol_MeV, M_sol_MeV/1000);
    printf("Core radius: ~3 code = %.2f fm\n", 3.0 * fm_per_code);
    printf("Potential half-width: ~5 code = %.2f fm\n\n", 5.0 * fm_per_code);

    printf("Level 1 potential depth: %.2f MeV\n", fabs(phi_min_L1) * M_sol_MeV);
    printf("Level 2 potential depth: %.2f MeV\n", fabs(phi_min_L2) * M_sol_MeV);

    /* GR comparison: Schwarzschild Phi/c^2 = -GM/(rc^2) */
    /* At r=3 (core): Phi_GR = -G M_p / (r_p c^2) ≈ -M_p/(M_Planck^2 * r_p) */
    /* G = 6.674e-11, M_p = 1.67e-27 kg, c = 3e8, r_p = 0.84 fm = 0.84e-15 m */
    double G_N = 6.674e-11;
    double M_proton = 1.6726e-27;
    double c_light = 2.998e8;
    double r_proton = 0.84e-15;
    double Phi_GR = -G_N * M_proton / (r_proton * c_light * c_light);

    printf("\nGR Schwarzschild at proton surface: Phi/c^2 = %.2e\n", Phi_GR);
    printf("Our Level 1:                        Phi/c^2 = %.2e\n", phi_min_L1);
    printf("Our Level 2:                        Phi/c^2 = %.2e\n", phi_min_L2);
    printf("Ratio (ours/GR): %.2e (nuclear vs gravitational)\n",
           phi_min_L1 / Phi_GR);

    /* N braids needed */
    printf("\n=== N-BRAID SCALING ===\n\n");
    printf("Potential at distance r from N-braid cluster (superposition):\n");
    printf("  Phi_N(r) ~ N * Phi_1(r)\n\n");
    printf("%-8s %-10s %-10s %-10s\n", "r", "Phi_1/c^2", "N for 0.01", "N for GR");
    for (int i = 0; i < n_pts; i += 4) {
        double phi1 = fabs(pts[i].phi_over_c2);
        if (phi1 < 1e-15) continue;
        double N_001 = 0.01 / phi1;
        double N_GR = fabs(Phi_GR) / phi1;
        printf("%-8.1f %-10.2e %-10.0f %-10.2e\n",
               pts[i].r, pts[i].phi_over_c2, N_001, N_GR);
    }

    free(pts);
    return 0;
}
