/*
 * normal_modes.c — Normal mode analysis of B=1 hedgehog Skyrmion
 *
 * Computes:
 *   1. Breathing mode (K=0) eigenfrequencies
 *   2. Moment of inertia Λ → Delta-N splitting
 *
 * BREATHING MODE:
 *
 * For f(r,t) = f₀(r) + ε g(r) cos(ωt), the linearized EOM gives:
 *
 *   -(P(r) g')' + W(r) g = (ω²/c²) m(r) g
 *
 * where the coefficients come from the second variation of E = E₂ + E₄ + E_V:
 *
 *   P(r) = 2r² + 4c₄ sin²f₀            (stiffness)
 *   m(r) = r² + 2c₄ sin²f₀ = P/2       (inertia)
 *
 *   W(r) = Q(r) - B'(r)                 (effective potential)
 *   Q(r) = 4cos(2f₀)(1 + c₄f₀'²) + 4c₄ sin²f₀(1+2cos(2f₀))/r²
 *   B(r) = 4c₄ f₀' sin(2f₀)
 *   B'(r) = 4c₄[f₀'' sin(2f₀) + 2f₀'² cos(2f₀)]
 *
 * With pion mass V = m_π² ∫(1-cosf) d³x:
 *   W → W + m_π² r² cos(f₀)
 *
 * Continuum threshold: ω² = 0 (massless) or ω² = m_π² c² (massive)
 *
 * MOMENT OF INERTIA:
 *
 *   Λ = (8πρ₀²)/(3c²) ∫ r² sin²f₀ [1 + c₄(f₀'² + sin²f₀/r²)] dr
 *   Delta-N splitting: ΔE = 3ℏ²c²/(2Λ) = 3/(2Λ) in natural units
 */

#define _USE_MATH_DEFINES
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ========== Profile I/O ========== */

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

    /* Compute f' numerically if not provided */
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

/* Interpolate f, f' at given r (linear) */
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

/* Compute f'' from the ODE (with optional pion mass):
 * f'' = [sin(2f)(1+c₄sin²f/r²) + m_π²r²sinf - 2f'r - c₄f'²sin(2f)] / (r²+2c₄sin²f) */
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
} NMPotential;

static NMPotential *build_potential(const Profile *prof, double c4,
                                     double m_pi2, int nr, double r_max)
{
    NMPotential *pot = calloc(1, sizeof(NMPotential));
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
        double c2f = cos(2*f), s2f = sin(2*f);
        double r2 = r*r;

        /* Stiffness and inertia */
        pot->P[ir] = 2*r2 + 4*c4*sf2;
        pot->m[ir] = r2 + 2*c4*sf2;   /* = P/2 */

        /* Q from second variation of E₂ + E₄ */
        double Q = 4*c2f*(1 + c4*fp*fp) + 4*c4*sf2*(1 + 2*c2f)/(r2 + 1e-30);

        /* B' = 4c₄[f₀''sin(2f₀) + 2f₀'²cos(2f₀)] */
        double Bp = 4*c4*(fpp*s2f + 2*fp*fp*c2f);

        /* W = Q - B' + pion mass contribution */
        pot->W[ir] = Q - Bp + m_pi2 * r2 * cos(f);
    }

    return pot;
}

static void free_potential(NMPotential *p) {
    free(p->r); free(p->P); free(p->W); free(p->m); free(p);
}

/* ========== Shooting solver ========== */

/*
 * The ODE in standard form (g'' = ...):
 *   -(Pg')' + Wg = λ m g
 *   -Pg'' - P'g' + Wg = λ m g
 *   g'' = -(P'/P) g' + (W - λm)/P · g
 *
 * We use a momentum formulation (p = P g') to avoid needing P':
 *   g' = p / P
 *   p' = (W - λm) g
 *
 * This is equivalent to -(Pg')' + Wg = λmg.
 */

typedef struct {
    double lambda;      /* ω²/c² trial eigenvalue */
    const NMPotential *pot;
} ShootParams;

/* Interpolate potential at given r */
static void interp_pot(const NMPotential *pot, double r,
                        double *P_out, double *W_out, double *m_out)
{
    double idx = r / pot->dr;
    int i0 = (int)idx;
    if (i0 < 0) i0 = 0;
    if (i0 >= pot->nr) i0 = pot->nr - 1;
    int i1 = (i0 < pot->nr) ? i0+1 : i0;
    double t = idx - i0;
    if (t < 0) t = 0;
    if (t > 1) t = 1;
    *P_out = (1-t)*pot->P[i0] + t*pot->P[i1];
    *W_out = (1-t)*pot->W[i0] + t*pot->W[i1];
    *m_out = (1-t)*pot->m[i0] + t*pot->m[i1];
}

/* RHS for (g, p=Pg') system */
static void gp_rhs(double r, double g, double p,
                    const ShootParams *sp, double *dg, double *dp)
{
    double P, W, m;
    interp_pot(sp->pot, r, &P, &W, &m);
    if (P < 1e-15) P = 1e-15;
    *dg = p / P;
    *dp = (W - sp->lambda * m) * g;
}

/* Compute the indicial exponent ν at r=0.
 * Near r=0: P ≈ 2r², m ≈ r², W → W₀.
 * Regular solution: g ~ r^ν where ν(ν+1) = W₀/2.
 * ν = (-1 + √(1 + 2W₀)) / 2. */
static double compute_nu(const NMPotential *pot)
{
    double W0 = pot->W[0];
    double disc = 1.0 + 2.0 * W0;
    if (disc < 0) return 0;   /* shouldn't happen for physical potential */
    return (-1.0 + sqrt(disc)) / 2.0;
}

/* Shoot from r=delta to r_max; return g(r_max), optionally fill wavefunction */
static double shoot(const ShootParams *sp, double *g_out, int nr_out)
{
    const NMPotential *pot = sp->pot;
    double dr = pot->dr;
    int nr = pot->nr;

    /* Regular BC at r=δ: g ~ r^ν, g' = ν r^{ν-1}, p = P g' = 2r²·ν r^{ν-1} = 2ν r^{ν+1} */
    double nu = compute_nu(pot);
    double delta = pot->r[0];
    double g = pow(delta, nu);
    double p = 2.0 * nu * pow(delta, nu + 1.0);  /* p = P(δ)·g'(δ) ≈ 2δ²·νδ^{ν-1} */

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

        /* Overflow protection */
        if (fabs(g) > 1e30 || fabs(p) > 1e30) {
            double s = 1e-20;
            g *= s; p *= s;
            if (g_out)
                for (int j = 0; j <= ir+1 && j < nr_out; j++) g_out[j] *= s;
        }
    }
    return g;
}

/* Find eigenvalues by scanning + bisection */
static int find_eigenvalues(const NMPotential *pot,
                            double lam_lo, double lam_hi,
                            double tol, double *eigenvalues, int max_ev)
{
    ShootParams sp = { .lambda = 0, .pot = pot };
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

/* ========== Moment of inertia ========== */

/*
 * Λ = (8πρ₀²)/(3c²) ∫₀^∞ r² sin²f [1 + c₄(f'² + sin²f/r²)] dr
 *
 * Delta-N splitting: ΔE = 3/(2Λ) (in code units where ℏ=c=1)
 */
static double compute_moment_of_inertia(const Profile *prof, double rho0, double c4)
{
    double sum = 0;
    for (int i = 1; i < prof->n - 1; i++) {
        double r = prof->r[i];
        double f = prof->f[i], fp = prof->fp[i];
        double sf = sin(f), sf2 = sf*sf;
        double r2 = r*r;

        double integrand = r2 * sf2 * (1.0 + c4*(fp*fp + sf2/(r2+1e-30)));
        double dr = (i == 0 || i == prof->n-1) ? prof->dr/2 : prof->dr;
        sum += integrand * dr;
    }
    return (8.0*M_PI*rho0*rho0/3.0) * sum;
}

/* ========== Variational consistency check ========== */

/* Compute ω² = E[g] / T[g] by direct quadrature */
static double variational_omega2(const NMPotential *pot, const double *g, int n)
{
    double E_num = 0, T_num = 0;
    double dr = pot->dr;

    for (int i = 0; i <= n && i <= pot->nr; i++) {
        double w = dr * ((i == 0 || i == pot->nr) ? 0.5 : 1.0);

        /* g' by central differences */
        double gp;
        if (i == 0)       gp = (g[1] - g[0]) / dr;
        else if (i >= n)  gp = (g[n] - g[n-1]) / dr;
        else              gp = (g[i+1] - g[i-1]) / (2*dr);

        E_num += w * (pot->P[i]*gp*gp + pot->W[i]*g[i]*g[i]);
        T_num += w * pot->m[i] * g[i]*g[i];
    }

    if (T_num < 1e-30) return 0;
    return E_num / T_num;
}

/* ========== Count nodes ========== */

static int count_nodes(const double *g, int n) {
    int nodes = 0;
    for (int i = 1; i < n; i++)
        if (g[i-1]*g[i] < 0) nodes++;
    return nodes;
}

/* ========== Main ========== */

int main(int argc, char *argv[])
{
    /* Default parameters */
    const char *profile_file = NULL;
    double e_skyrme = 4.0;
    double rho0 = 1.0;
    double m_pi = 0.0;      /* pion mass in code units */
    double r_max = 15.0;
    int nr = 4000;
    int print_potential = 0;
    int print_wavefunction = 0;

    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-profile") == 0 && i+1 < argc) profile_file = argv[++i];
        else if (strcmp(argv[i], "-e") == 0 && i+1 < argc) e_skyrme = atof(argv[++i]);
        else if (strcmp(argv[i], "-rho0") == 0 && i+1 < argc) rho0 = atof(argv[++i]);
        else if (strcmp(argv[i], "-mpi") == 0 && i+1 < argc) m_pi = atof(argv[++i]);
        else if (strcmp(argv[i], "-rmax") == 0 && i+1 < argc) r_max = atof(argv[++i]);
        else if (strcmp(argv[i], "-nr") == 0 && i+1 < argc) nr = atoi(argv[++i]);
        else if (strcmp(argv[i], "-pot") == 0) print_potential = 1;
        else if (strcmp(argv[i], "-wf") == 0) print_wavefunction = 1;
        else {
            fprintf(stderr, "Usage: %s -profile <file> [-e val] [-rho0 val] "
                    "[-mpi pion_mass] [-rmax val] [-nr val] [-pot] [-wf]\n", argv[0]);
            return 1;
        }
    }

    if (!profile_file) {
        fprintf(stderr, "Must specify -profile <file>\n");
        return 1;
    }

    double c4 = 2.0 * rho0 * rho0 / (e_skyrme * e_skyrme);
    double m_pi2 = m_pi * m_pi;

    printf("===== Hedgehog Normal Mode Analysis =====\n");
    printf("Profile: %s\n", profile_file);
    printf("e=%.4f, rho0=%.4f, c4=%.6f\n", e_skyrme, rho0, c4);
    printf("Pion mass: m_pi=%.4f (m_pi²=%.6f)\n", m_pi, m_pi2);
    printf("Grid: nr=%d, r_max=%.1f, dr=%.6f\n", nr, r_max, r_max/nr);

    /* Read profile */
    Profile prof;
    if (read_profile(profile_file, &prof) != 0) return 1;
    printf("Profile: %d points, r_max=%.3f\n", prof.n, prof.r[prof.n-1]);

    /* ===== Moment of Inertia ===== */
    double Lambda = compute_moment_of_inertia(&prof, rho0, c4);
    double delta_N = 1.5 / Lambda;   /* ΔE = 3/(2Λ) in natural units */
    printf("\n===== Moment of Inertia =====\n");
    printf("Lambda = %.6f\n", Lambda);
    printf("Delta-N splitting = 3/(2*Lambda) = %.6f code units\n", delta_N);

    /* Physical units (conversion from paramfit_analysis.md at e=1, rho0=1) */
    double Lambda_E = 9.098;  /* MeV per code energy unit */
    /* Scale by actual e, rho0: E_code scales as rho0^3/e, so Lambda_E scales too */
    /* Actually the conversion is universal: 1 code unit = 9.098 MeV at e=1, rho0=1 */
    /* For general e, rho0: energy scale is rho0^3/e relative to e=1,rho0=1 */
    double energy_scale = rho0*rho0*rho0 / e_skyrme;
    double delta_N_MeV = delta_N * Lambda_E * energy_scale / (1.0/e_skyrme); /* need to think about this */
    /* Actually, all code units are at the specified e, rho0 already.
     * The profile was solved at those parameters. The moment of inertia
     * integral uses the profile. So ΔE is in "code units" for those params.
     * The physical conversion is: E_phys = E_code × (K × rho0^3/e) × (M_p / E_sol_code)
     * But E_sol_code = K = 103.13 at e=1,rho0=1. More simply:
     * 1 code energy = M_p / E_sol = 938.3 / 103.13 = 9.098 MeV when e=1, rho0=1.
     * For general e, rho0: E_sol = 103.13 × rho0^3/e, and M_p is fixed,
     * so 1 code energy = M_p × e / (103.13 × rho0^3) MeV.
     */
    double code_to_MeV = 938.272 * e_skyrme / (103.13 * rho0*rho0*rho0);
    delta_N_MeV = delta_N * code_to_MeV;
    printf("Energy conversion: 1 code unit = %.3f MeV (at e=%.2f, rho0=%.2f)\n",
           code_to_MeV, e_skyrme, rho0);
    printf("Delta-N splitting = %.1f MeV (experiment: 293.7 MeV)\n", delta_N_MeV);

    /* ===== Build Potential ===== */
    NMPotential *pot = build_potential(&prof, c4, m_pi2, nr, r_max);

    if (print_potential) {
        printf("\n# Potential profile:\n");
        printf("# r  P(r)  m(r)  W(r)  W/m\n");
        int step = nr / 500;
        if (step < 1) step = 1;
        for (int i = 0; i <= pot->nr; i += step) {
            double ratio = (pot->m[i] > 1e-20) ? pot->W[i]/pot->m[i] : 0;
            printf("%.6f  %.6f  %.6f  %.6f  %.6f\n",
                   pot->r[i], pot->P[i], pot->m[i], pot->W[i], ratio);
        }
    }

    /* ===== Eigenvalue Search ===== */
    printf("\n===== Breathing Mode Eigenvalues =====\n");

    /* Determine search range.
     * Continuum threshold: ω² = m_π² c² (massive) or 0 (massless).
     * For box quantization (massless), eigenvalues are at ω² > 0.
     * For massive pions: bound states at 0 < ω² < m_π².
     *
     * The lowest breathing mode is typically at large ω²
     * (the breathing frequency is high). Search a wide range.
     *
     * W/m at large r → 0 (sigma model) or m_π² (massive).
     * W/m near origin can be large positive or negative.
     */

    /* Estimate W/m range for scan limits */
    double Wm_min = 1e30, Wm_max = -1e30;
    for (int i = 1; i <= pot->nr; i++) {
        if (pot->m[i] < 1e-10) continue;
        double ratio = pot->W[i] / pot->m[i];
        if (ratio < Wm_min) Wm_min = ratio;
        if (ratio > Wm_max) Wm_max = ratio;
    }
    printf("W/m range: [%.4f, %.4f]\n", Wm_min, Wm_max);

    /* For box quantization (g(R_max)=0), scan λ = ω²/c² from near-zero upward.
     * The "box modes" have λ_n ≈ (nπ/R_max)² × 2 at large n (free particle).
     * The breathing mode is the one with lowest λ that has special structure
     * (concentrated near the soliton).
     *
     * Scan range: [0, λ_max] where λ_max covers several box modes */
    double lam_lo = -50.0;    /* allow negative for massive case */
    double lam_hi = 200.0;    /* should be enough to see first few modes */
    double tol = 1e-10;

    double eigenvalues[50];
    int n_ev = find_eigenvalues(pot, lam_lo, lam_hi, tol, eigenvalues, 50);

    printf("\nFound %d box eigenvalues (g(R_max)=0) in [%.1f, %.1f]:\n",
           n_ev, lam_lo, lam_hi);

    if (m_pi > 0)
        printf("Continuum threshold: lambda = m_pi^2 = %.6f\n", m_pi2);
    else
        printf("Continuum threshold: lambda = 0 (massless pions)\n");

    /* Compute and print wavefunctions for each eigenvalue */
    double *g_wf = malloc((nr+1) * sizeof(double));
    ShootParams sp = { .lambda = 0, .pot = pot };

    double nu = compute_nu(pot);
    printf("Indicial exponent nu = %.4f (g ~ r^nu near r=0)\n", nu);

    for (int ie = 0; ie < n_ev && ie < 20; ie++) {
        sp.lambda = eigenvalues[ie];
        shoot(&sp, g_wf, nr+1);

        /* Normalize max|g|=1 before analysis */
        double gmax = 0;
        for (int i = 0; i <= nr; i++)
            if (fabs(g_wf[i]) > gmax) gmax = fabs(g_wf[i]);
        if (gmax > 0) for (int i = 0; i <= nr; i++) g_wf[i] /= gmax;

        int nodes = count_nodes(g_wf, nr+1);
        double omega2 = eigenvalues[ie];  /* = ω²/c² with c=1 */
        double omega = (omega2 > 0) ? sqrt(omega2) : -sqrt(-omega2);
        double omega_MeV = omega * code_to_MeV;
        double var_check = variational_omega2(pot, g_wf, nr);

        printf("  n=%2d: lam=%.6f  omega=%.4f  E=%.1f MeV  nodes=%d  var=%.4f",
               ie, eigenvalues[ie], omega, omega_MeV, nodes, var_check);

        if (m_pi > 0 && omega2 < m_pi2)
            printf("  ** BOUND (below m_pi=%.4f) **", m_pi);
        printf("\n");
    }

    /* Print wavefunction of lowest mode */
    if (print_wavefunction && n_ev > 0) {
        printf("\n# Wavefunction of lowest mode (lambda=%.6f):\n", eigenvalues[0]);
        printf("# r  g(r)  P(r)  W(r)  m(r)\n");
        sp.lambda = eigenvalues[0];
        shoot(&sp, g_wf, nr+1);
        /* Normalize: max|g| = 1 */
        double gmax = 0;
        for (int i = 0; i <= nr; i++)
            if (fabs(g_wf[i]) > gmax) gmax = fabs(g_wf[i]);
        if (gmax > 0) for (int i = 0; i <= nr; i++) g_wf[i] /= gmax;

        int step = nr / 500;
        if (step < 1) step = 1;
        for (int i = 0; i <= nr; i += step)
            printf("%.6f  %.6e  %.6f  %.6f  %.6f\n",
                   pot->r[i], g_wf[i], pot->P[i], pot->W[i], pot->m[i]);
    }

    /* ===== Breathing mode frequency in physical units ===== */
    if (n_ev > 0) {
        printf("\n===== Physical Interpretation =====\n");
        double omega0 = sqrt(fabs(eigenvalues[0]));
        double freq0 = omega0 / (2*M_PI);
        double period0 = 1.0 / freq0;
        printf("Lowest breathing mode:\n");
        printf("  omega = %.4f code  (%.1f MeV)\n", omega0, omega0*code_to_MeV);
        printf("  Period = %.4f code time units\n", period0);
        printf("  (For comparison: light-crossing time 2R ≈ 3.0 code)\n");

        /* Check consistency with observed 19 MeV spectral peak */
        printf("\n  Observed spectral peak from B+Bbar scattering: ~19 MeV\n");
        printf("  Breathing mode prediction: %.1f MeV\n", omega0*code_to_MeV);
    }

    free(g_wf);
    free_potential(pot);
    free(prof.r); free(prof.f); free(prof.fp);

    return 0;
}
