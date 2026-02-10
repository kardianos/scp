/*
 * angular_modes.c — K=1 angular mode analysis of B=1 hedgehog Skyrmion
 *
 * Computes the K=1 (grand spin 1) perturbation spectrum to check for
 * bound states below the pion mass threshold.
 *
 * Three K=1 channels are analyzed:
 *
 * Channel A: ISOROTATIONAL (L=0, I=1)
 *   δq = g(r) × iσ_a q₀  (isospin rotation by angle g(r))
 *   - Kinetic: T_iso = (1/2c²) ∫ Ω(r) ġ² 4πr² dr
 *   - Potential: δ²E = ∫ Ω(r) g'² 4πr² dr  (NO pion mass — isospin invariant)
 *   - Ω(r) = ρ₀² sin²f [1 + c₄(f'² + sin²f/r²)]
 *   - Zero mode at ω=0 (g=const), continuum at ω≥0
 *   - Moment of inertia: Λ = (8π/3c²) ∫ r² Ω(r) dr
 *
 * Channel B: TRANSLATIONAL (L=1, I=0)
 *   δq = g(r) × ∂q₀/∂z  (displacement along z)
 *   - Has centrifugal barrier 2/r² from L=1
 *   - Zero mode at ω=0 (g=const → uniform translation)
 *
 * Channel C: MIXED (L=1, I=1, K=1)
 *   Both orbital and isospin angular momentum.
 *   - Centrifugal barrier 2/r² from L=1
 *   - Isospin coupling modifies the potential
 *   - Pion mass contributes: W_pi(r) = m_π² cos(f) (can be attractive)
 *   - Question: does the attractive well overcome the centrifugal barrier?
 *
 * For each channel, the Sturm-Liouville equation is:
 *   -(P g')' + W g = (ω²/c²) m g
 *
 * Solved by shooting with momentum formulation (p = Pg').
 */

#define _USE_MATH_DEFINES
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ========== Profile I/O (from normal_modes.c) ========== */

typedef struct {
    int n;
    double *r, *f, *fp;
    double dr;
} Profile;

static int read_profile(const char *fname, Profile *prof)
{
    FILE *fp = fopen(fname, "r");
    if (!fp) { fprintf(stderr, "Cannot open %s\n", fname); return -1; }

    int cap = 4096;
    prof->r  = malloc(cap * sizeof(double));
    prof->f  = malloc(cap * sizeof(double));
    prof->fp = malloc(cap * sizeof(double));
    int n = 0;
    char line[512];

    while (fgets(line, sizeof(line), fp)) {
        if (line[0] == '#') continue;
        double rv, fval, fpval = 0;
        double d1, d2, d3;
        int nc = sscanf(line, "%lf %lf %lf %lf %lf %lf",
                        &rv, &fval, &fpval, &d1, &d2, &d3);
        if (nc < 2) continue;
        if (n >= cap) {
            cap *= 2;
            prof->r  = realloc(prof->r,  cap * sizeof(double));
            prof->f  = realloc(prof->f,  cap * sizeof(double));
            prof->fp = realloc(prof->fp, cap * sizeof(double));
        }
        prof->r[n] = rv; prof->f[n] = fval; prof->fp[n] = fpval;
        n++;
    }
    fclose(fp);
    prof->n = n;
    prof->dr = (n > 1) ? prof->r[1] - prof->r[0] : 0.001;

    int has_fp = 0;
    for (int i = 0; i < n; i++)
        if (fabs(prof->fp[i]) > 1e-30) { has_fp = 1; break; }
    if (!has_fp) {
        for (int i = 0; i < n; i++) {
            if (i == 0) prof->fp[i] = (prof->f[1]-prof->f[0])/(prof->r[1]-prof->r[0]);
            else if (i == n-1) prof->fp[i] = (prof->f[n-1]-prof->f[n-2])/(prof->r[n-1]-prof->r[n-2]);
            else prof->fp[i] = (prof->f[i+1]-prof->f[i-1])/(prof->r[i+1]-prof->r[i-1]);
        }
    }
    return 0;
}

static void interp_profile(const Profile *prof, double r,
                            double *f_out, double *fp_out)
{
    if (r <= prof->r[0]) {
        *f_out = prof->f[0]; *fp_out = prof->fp[0]; return;
    }
    if (r >= prof->r[prof->n-1]) {
        *f_out = prof->f[prof->n-1]; *fp_out = prof->fp[prof->n-1]; return;
    }
    int lo = 0, hi = prof->n - 1;
    while (hi - lo > 1) { int mid = (lo+hi)/2; if (prof->r[mid] <= r) lo = mid; else hi = mid; }
    double t = (r - prof->r[lo]) / (prof->r[hi] - prof->r[lo]);
    *f_out  = prof->f[lo]  + t * (prof->f[hi]  - prof->f[lo]);
    *fp_out = prof->fp[lo] + t * (prof->fp[hi] - prof->fp[lo]);
}

static double compute_fpp(double r, double f, double fp, double c4, double m_pi2)
{
    double sf = sin(f), sf2 = sf*sf, s2f = sin(2*f);
    double denom = r*r + 2*c4*sf2;
    if (denom < 1e-30) return 0.0;
    double num = s2f*(1.0 + c4*sf2/(r*r + 1e-30))
               + m_pi2*r*r*sf
               - 2*fp*r - c4*fp*fp*s2f;
    return num / denom;
}

/* ========== Potential construction ========== */

typedef struct {
    int nr;
    double dr, r_max;
    double *r, *P, *W, *m;
} AMPotential;

/*
 * Channel types:
 *   0 = Isorotational (L=0, I=1)
 *   1 = Translational (L=1, I=0)
 *   2 = Mixed pion-like (L=1, I=1, K=1)
 */
static AMPotential *build_potential(const Profile *prof, double c4,
                                     double m_pi2, double rho0,
                                     int nr, double r_max, int channel)
{
    AMPotential *pot = calloc(1, sizeof(AMPotential));
    pot->nr = nr;
    pot->dr = r_max / nr;
    pot->r_max = r_max;
    pot->r = calloc(nr+1, sizeof(double));
    pot->P = calloc(nr+1, sizeof(double));
    pot->W = calloc(nr+1, sizeof(double));
    pot->m = calloc(nr+1, sizeof(double));

    for (int ir = 0; ir <= nr; ir++) {
        double r = (ir == 0) ? 1e-6 : ir * pot->dr;
        pot->r[ir] = r;

        double f, fp;
        interp_profile(prof, r, &f, &fp);
        double fpp = compute_fpp(r, f, fp, c4, m_pi2);

        double sf = sin(f), sf2 = sf*sf;
        double cf = cos(f);
        double c2f = cos(2*f), s2f = sin(2*f);
        double r2 = r*r;

        /* Isospin inertia density: Ω(r) = ρ₀² sin²f [1 + c₄(f'² + sin²f/r²)] */
        double Omega = rho0*rho0 * sf2 * (1.0 + c4*(fp*fp + sf2/(r2+1e-30)));

        switch (channel) {
        case 0: /* Isorotational (L=0, I=1) */
            /* -(r²Ω g')' = ω² r²Ω g  →  P = r²Ω, W = 0, m = r²Ω */
            pot->P[ir] = r2 * Omega;
            pot->W[ir] = 0.0;  /* Pion mass is isospin-invariant! */
            pot->m[ir] = r2 * Omega;
            break;

        case 1: { /* Translational (L=1, I=0) */
            /* Uses breathing mode coefficients + centrifugal barrier 2/r² */
            double P_br = 2*r2 + 4*c4*sf2;
            double m_br = r2 + 2*c4*sf2;

            /* Q and B' from second variation of E₂ + E₄ */
            double Q = 4*c2f*(1 + c4*fp*fp) + 4*c4*sf2*(1 + 2*c2f)/(r2 + 1e-30);
            double Bp = 4*c4*(fpp*s2f + 2*fp*fp*c2f);
            double W_br = Q - Bp + m_pi2 * r2 * cf;

            /* Add centrifugal barrier: L(L+1)/r² × (coefficient) */
            /* For L=1: 2/r² acts on the "kinetic" part. In the SL equation,
             * the centrifugal term adds 2 × P_br / r² to W */
            double W_cent = 2.0 * P_br / (r2 + 1e-30);

            pot->P[ir] = P_br;
            pot->W[ir] = W_br + W_cent;
            pot->m[ir] = m_br;
            break;
        }

        case 2: { /* Mixed (L=1, I=1, K=1) */
            /* The perturbation has both L=1 orbital and I=1 isospin.
             *
             * Stiffness: contributions from E₂ (gradient) and E₄ (Skyrme).
             * E₂ contribution: P_E2 = r² + centrifugal from L=1
             * E₄ contribution: involves angular integrals of cross-products
             *
             * Effective potential includes:
             * - Centrifugal barrier: 2/r² × P(r)  [repulsive]
             * - Skyrme coupling: angular-averaged E₄ [repulsive]
             * - Pion mass: m_π² r² sin²f cosf  [attractive near core]
             *
             * We use an approximation: the mixed channel SL equation uses
             * the isospin-weighted stiffness with centrifugal barrier.
             */

            /* P from isospin + orbital gradient:
             * P_mix = r²[1 + c₄(f'² + sin²f/r²)] for the isospin part
             * + radial stiffness from the orbital L=1 component */
            double P_mix = r2 * (1.0 + c4*(fp*fp + sf2/(r2+1e-30)));

            /* Centrifugal barrier from L=1 */
            double W_cent = 2.0 * P_mix / (r2 + 1e-30);

            /* Isospin coupling from E₄ cross-term: angular integral gives
             * W_E4 ~ c₄ × sin²f × (f'² + sin²f/r²) / r²
             * This is repulsive (stabilizing). */
            double W_E4 = c4 * sf2 * (fp*fp + sf2/(r2+1e-30)) / (r2 + 1e-30);

            /* Pion mass contribution for the mixed mode:
             * δ²V_pi = m_π² ρ₀² cosf × g² (same form as breathing mode).
             * At large r (f→0, cosf→1): W_pi → m_π²r² = mass gap.
             * Near core (f→π, cosf→-1): W_pi < 0 (attractive well). */
            double W_pi = m_pi2 * r2 * cf;

            pot->P[ir] = P_mix;
            pot->W[ir] = W_cent + W_E4 + W_pi;
            pot->m[ir] = P_mix;
            break;
        }
        }
    }

    return pot;
}

static void free_potential(AMPotential *p) {
    free(p->r); free(p->P); free(p->W); free(p->m); free(p);
}

/* ========== Shooting solver (from normal_modes.c) ========== */

typedef struct {
    double lambda;
    const AMPotential *pot;
    int L;  /* orbital angular momentum (affects BC) */
} ShootParams;

static void interp_pot(const AMPotential *pot, double r,
                        double *P_out, double *W_out, double *m_out)
{
    double idx_f = r / pot->dr;
    int i0 = (int)idx_f;
    if (i0 < 0) i0 = 0;
    if (i0 >= pot->nr) i0 = pot->nr - 1;
    int i1 = (i0 < pot->nr) ? i0+1 : i0;
    double t = idx_f - i0;
    if (t < 0) t = 0; if (t > 1) t = 1;
    *P_out = (1-t)*pot->P[i0] + t*pot->P[i1];
    *W_out = (1-t)*pot->W[i0] + t*pot->W[i1];
    *m_out = (1-t)*pot->m[i0] + t*pot->m[i1];
}

static void gp_rhs(double r, double g, double p,
                    const ShootParams *sp, double *dg, double *dp)
{
    double P, W, m;
    interp_pot(sp->pot, r, &P, &W, &m);
    if (P < 1e-15) P = 1e-15;
    *dg = p / P;
    *dp = (W - sp->lambda * m) * g;
}

/* Compute indicial exponent at r=0 for given L.
 * For L=0: g ~ r^ν where ν = (-1+√(1+2W₀))/2
 * For L=1: g ~ r near origin (ν=1) */
static double compute_nu(const AMPotential *pot, int L)
{
    if (L >= 1) return (double)L;
    double W0 = pot->W[0];
    double disc = 1.0 + 2.0 * W0;
    if (disc < 0) return 0;
    return (-1.0 + sqrt(disc)) / 2.0;
}

static double shoot(const ShootParams *sp, double *g_out, int nr_out)
{
    const AMPotential *pot = sp->pot;
    double dr = pot->dr;
    int nr = pot->nr;

    double nu = compute_nu(pot, sp->L);
    double delta = pot->r[0];

    /* BC: g ~ r^ν, g' = ν r^{ν-1}, p = P g' */
    double g = pow(delta, nu);
    double P0 = pot->P[0];
    double p = P0 * nu * pow(delta, nu - 1.0);

    if (g_out && nr_out > 0) g_out[0] = g;

    for (int ir = 0; ir < nr; ir++) {
        double r = pot->r[ir];
        if (r < 1e-10) r = 1e-10;

        double k1g, k1p, k2g, k2p, k3g, k3p, k4g, k4p;
        double rmid = r + 0.5*dr, rend = r + dr;

        gp_rhs(r,    g,              p,              sp, &k1g, &k1p);
        gp_rhs(rmid, g+0.5*dr*k1g,  p+0.5*dr*k1p,  sp, &k2g, &k2p);
        gp_rhs(rmid, g+0.5*dr*k2g,  p+0.5*dr*k2p,  sp, &k3g, &k3p);
        gp_rhs(rend, g+dr*k3g,      p+dr*k3p,       sp, &k4g, &k4p);

        g += (dr/6.0)*(k1g + 2*k2g + 2*k3g + k4g);
        p += (dr/6.0)*(k1p + 2*k2p + 2*k3p + k4p);

        if (g_out && ir+1 < nr_out) g_out[ir+1] = g;

        if (fabs(g) > 1e30 || fabs(p) > 1e30) {
            double s = 1e-20;
            g *= s; p *= s;
            if (g_out)
                for (int j = 0; j <= ir+1 && j < nr_out; j++) g_out[j] *= s;
        }
    }
    return g;
}

static int find_eigenvalues(const AMPotential *pot, int L,
                            double lam_lo, double lam_hi,
                            double tol, double *eigenvalues, int max_ev)
{
    ShootParams sp = { .lambda = 0, .pot = pot, .L = L };
    int n_scan = 5000;
    double d_lam = (lam_hi - lam_lo) / n_scan;

    double brackets_lo[100], brackets_hi[100];
    int n_brackets = 0;

    sp.lambda = lam_lo;
    double g_prev = shoot(&sp, NULL, 0);

    for (int is = 1; is <= n_scan; is++) {
        double lam = lam_lo + is * d_lam;
        sp.lambda = lam;
        double g_cur = shoot(&sp, NULL, 0);

        if (g_prev * g_cur < 0 && n_brackets < 100) {
            brackets_lo[n_brackets] = lam - d_lam;
            brackets_hi[n_brackets] = lam;
            n_brackets++;
        }
        g_prev = g_cur;
    }

    int n_ev = 0;
    for (int ib = 0; ib < n_brackets && n_ev < max_ev; ib++) {
        double lo = brackets_lo[ib], hi = brackets_hi[ib];
        sp.lambda = lo;
        double g_lo = shoot(&sp, NULL, 0);

        for (int bi = 0; bi < 200; bi++) {
            double mid = 0.5*(lo + hi);
            sp.lambda = mid;
            double g_mid = shoot(&sp, NULL, 0);
            if (g_mid * g_lo > 0) { lo = mid; g_lo = g_mid; }
            else hi = mid;
            if (hi - lo < tol) break;
        }
        eigenvalues[n_ev++] = 0.5*(lo + hi);
    }
    return n_ev;
}

static int count_nodes(const double *g, int n) {
    int nodes = 0;
    for (int i = 1; i < n; i++)
        if (g[i-1]*g[i] < 0) nodes++;
    return nodes;
}

/* ========== Main ========== */

int main(int argc, char *argv[])
{
    const char *profile_file = NULL;
    double e_skyrme = 1.0;
    double rho0 = 1.0;
    double m_pi = 0.0;
    double r_max = 15.0;
    int nr = 4000;
    int print_potential = 0;

    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-profile") == 0 && i+1 < argc) profile_file = argv[++i];
        else if (strcmp(argv[i], "-e") == 0 && i+1 < argc) e_skyrme = atof(argv[++i]);
        else if (strcmp(argv[i], "-rho0") == 0 && i+1 < argc) rho0 = atof(argv[++i]);
        else if (strcmp(argv[i], "-mpi") == 0 && i+1 < argc) m_pi = atof(argv[++i]);
        else if (strcmp(argv[i], "-rmax") == 0 && i+1 < argc) r_max = atof(argv[++i]);
        else if (strcmp(argv[i], "-nr") == 0 && i+1 < argc) nr = atoi(argv[++i]);
        else if (strcmp(argv[i], "-pot") == 0) print_potential = 1;
        else {
            fprintf(stderr, "Usage: %s -profile <file> [-e val] [-rho0 val] "
                    "[-mpi pion_mass] [-rmax val] [-nr val] [-pot]\n", argv[0]);
            return 1;
        }
    }

    if (!profile_file) {
        fprintf(stderr, "Must specify -profile <file>\n");
        return 1;
    }

    double c4 = 2.0 * rho0 * rho0 / (e_skyrme * e_skyrme);
    double m_pi2 = m_pi * m_pi;
    double code_to_MeV = 938.272 * e_skyrme / (103.13 * rho0*rho0*rho0);

    printf("===== K=1 Angular Mode Analysis =====\n");
    printf("Profile: %s\n", profile_file);
    printf("e=%.4f, rho0=%.4f, c4=%.6f\n", e_skyrme, rho0, c4);
    printf("Pion mass: m_pi=%.4f (m_pi²=%.6f)\n", m_pi, m_pi2);
    printf("Energy conversion: 1 code unit = %.3f MeV\n", code_to_MeV);

    Profile prof;
    if (read_profile(profile_file, &prof) != 0) return 1;
    printf("Profile: %d points, r_max=%.3f\n\n", prof.n, prof.r[prof.n-1]);

    const char *channel_names[] = {
        "Isorotational (L=0, I=1)",
        "Translational (L=1, I=0)",
        "Mixed pion-like (L=1, I=1, K=1)"
    };
    int channel_L[] = {0, 1, 1};  /* orbital angular momentum */

    double lam_lo = -50.0, lam_hi = 200.0;
    double tol = 1e-10;

    for (int ch = 0; ch < 3; ch++) {
        printf("===== Channel %c: %s =====\n", 'A'+ch, channel_names[ch]);

        AMPotential *pot = build_potential(&prof, c4, m_pi2, rho0,
                                           nr, r_max, ch);

        /* Print W/m ratio range */
        double Wm_min = 1e30, Wm_max = -1e30;
        for (int i = 1; i <= pot->nr; i++) {
            if (pot->m[i] < 1e-10) continue;
            double ratio = pot->W[i] / pot->m[i];
            if (ratio < Wm_min) Wm_min = ratio;
            if (ratio > Wm_max) Wm_max = ratio;
        }
        printf("W/m range: [%.4f, %.4f]\n", Wm_min, Wm_max);

        if (print_potential) {
            printf("# r  P(r)  W(r)  m(r)  W/m\n");
            int step = nr / 200;
            if (step < 1) step = 1;
            for (int i = 0; i <= pot->nr; i += step) {
                double ratio = (pot->m[i] > 1e-20) ? pot->W[i]/pot->m[i] : 0;
                printf("%.6f  %.6e  %.6e  %.6e  %.6f\n",
                       pot->r[i], pot->P[i], pot->W[i], pot->m[i], ratio);
            }
        }

        /* Continuum threshold */
        double threshold = m_pi2;
        if (ch == 0) threshold = 0;  /* isorotational: massless continuum */

        printf("Continuum threshold: lambda = %.6f", threshold);
        if (threshold > 0) printf(" (%.1f MeV)", sqrt(threshold) * code_to_MeV);
        printf("\n");

        /* Search for eigenvalues */
        double eigenvalues[50];
        int L = channel_L[ch];
        double nu = compute_nu(pot, L);
        printf("Indicial exponent: nu = %.4f (g ~ r^nu near r=0)\n", nu);

        int n_ev = find_eigenvalues(pot, L, lam_lo, lam_hi, tol, eigenvalues, 50);

        /* Count bound states (below threshold) */
        int n_bound = 0;
        for (int ie = 0; ie < n_ev; ie++)
            if (eigenvalues[ie] < threshold - 1e-6) n_bound++;

        printf("Found %d eigenvalues, %d below threshold\n", n_ev, n_bound);

        /* Print eigenvalues */
        double *g_wf = malloc((nr+1) * sizeof(double));
        ShootParams sp = { .lambda = 0, .pot = pot, .L = L };

        for (int ie = 0; ie < n_ev && ie < 15; ie++) {
            sp.lambda = eigenvalues[ie];
            shoot(&sp, g_wf, nr+1);

            double gmax = 0;
            for (int i = 0; i <= nr; i++)
                if (fabs(g_wf[i]) > gmax) gmax = fabs(g_wf[i]);
            if (gmax > 0) for (int i = 0; i <= nr; i++) g_wf[i] /= gmax;

            int nodes = count_nodes(g_wf, nr+1);
            double omega2 = eigenvalues[ie];
            double omega = (omega2 > 0) ? sqrt(omega2) : -sqrt(-omega2);
            double omega_MeV = omega * code_to_MeV;

            printf("  n=%2d: lam=%.6f  omega=%.4f  E=%.1f MeV  nodes=%d",
                   ie, omega2, omega, omega_MeV, nodes);

            if (omega2 < threshold - 1e-6)
                printf("  ** BOUND **");
            else if (ie == 0 && fabs(omega2) < 0.01)
                printf("  ** ZERO MODE **");
            printf("\n");
        }

        free(g_wf);

        /* Summary for this channel */
        if (n_bound == 0) {
            printf("→ No bound states in channel %c.\n", 'A'+ch);
            if (ch == 0)
                printf("  (Isorotational zero mode: pion from collective quantization)\n");
        } else {
            printf("→ %d BOUND STATE(S) found in channel %c!\n", n_bound, 'A'+ch);
        }
        printf("\n");

        free_potential(pot);
    }

    /* ===== Summary ===== */
    printf("===== Summary =====\n");
    printf("K=1 angular mode spectrum for B=1 hedgehog Skyrmion:\n");
    if (m_pi > 0) {
        printf("  Pion mass: m_pi = %.4f code = %.1f MeV\n",
               m_pi, m_pi * code_to_MeV);
        printf("  Continuum threshold: omega = m_pi = %.4f\n", m_pi);
    }
    printf("  Channel A (isorotation): zero mode → pion from quantization\n");
    printf("  Channel B (translation): zero mode → momentum\n");
    printf("  Channel C (mixed L=1,I=1): check for pion-like bound states\n");
    printf("\n");
    printf("Physical pion mass from collective quantization:\n");

    /* Compute moment of inertia */
    double Lambda_sum = 0;
    for (int i = 1; i < prof.n - 1; i++) {
        double r = prof.r[i], f = prof.f[i], fp = prof.fp[i];
        double sf = sin(f), sf2 = sf*sf, r2 = r*r;
        double integrand = r2 * sf2 * (1.0 + c4*(fp*fp + sf2/(r2+1e-30)));
        Lambda_sum += integrand * prof.dr;
    }
    double Lambda = (8.0*M_PI*rho0*rho0/3.0) * Lambda_sum;
    double delta_N = 1.5 / Lambda;
    double delta_N_MeV = delta_N * code_to_MeV;

    printf("  Moment of inertia: Lambda = %.4f\n", Lambda);
    printf("  Delta-N splitting: 3/(2*Lambda) = %.6f code = %.1f MeV\n",
           delta_N, delta_N_MeV);
    printf("  (Experiment: Delta-N = 293.7 MeV, m_pi = 139.6 MeV)\n");

    free(prof.r); free(prof.f); free(prof.fp);
    return 0;
}
