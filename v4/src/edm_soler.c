/*
 * edm_soler.c — Einstein-Dirac-Maxwell(-Soler) solver
 *
 * Self-consistent field (SCF) iteration for coupled gravity + EM:
 *   1. Fix EM potential φ(r), solve Dirac (+Einstein if G>0) via bisection
 *   2. From spinor, compute new φ via Green's function / forward integration
 *   3. Under-relax φ, repeat until convergence
 *
 * Dirac+Einstein (4 coupled ODEs, evolved with fixed φ):
 *   √A·α' = (κ/r)α - ((ω-eφ)T + m_eff)β
 *   √A·β' = ((ω-eφ)T - m_eff)α - (κ/r)β
 *   rA' = 1-A - 8πGN(ω-eφ)T²(α²+β²) - GΨ²/r² - 4πGr²λV²
 *   2rA(T'/T) = A-1 - 8πGNT√A(αβ'-βα') + GΨ²/r² + 4πGr²λV²
 *
 * Maxwell (solved separately via Green's function, φ(∞)=0 gauge):
 *   V(r) = eN[M_enc(r)/r + I_ext(r)]    (Coulomb self-potential, V > 0)
 *   φ(r) = -V(r)                          (FSY convention, φ ≤ 0)
 *   ω_eff = ω - eφ = ω + eV ≥ ω          (repulsive Coulomb)
 *
 * Convention:
 *   β = -r·g (minus sign vs soler.c)
 *   λ_edm = λ_soler/N (default 0.5 for N=2 matches soler.c λ=1)
 *   N=2 fermions → total charge = 2eQ_norm per orbital
 *   Gaussian CGS units for EM (4π in Gauss's law)
 *
 * Refs: Finster, Smoller, Yau, Phys. Lett. A 259, 431-436 (1999)
 *       gr-qc/9802012  (EDM particlelike solutions)
 *       Zhang & He, arXiv:2508.11714 (2025)  (ED-Soler)
 *
 * Usage:
 *   edm_soler -omega 0.9 -e 0.1                (flat + charge)
 *   edm_soler -omega 0.9 -e 0.1 -G 0.001       (gravity + charge)
 *   edm_soler -omega 0.9 -e 0 -G 0.001         (ED-Soler, no EM)
 *   edm_soler -scan -e 0.1                       (scan omega)
 *   edm_soler -scan_e -omega 0.9                 (scan charge)
 */

#define _DEFAULT_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define NMAX 200001
#define NEQN 4    /* Dirac+Einstein: α, β, A, T */

static double al_arr[NMAX];    /* α(r) */
static double be_arr[NMAX];    /* β(r) */
static double A_arr[NMAX];     /* A(r) metric */
static double T_arr[NMAX];     /* T(r) metric */
static double phi_arr[NMAX];   /* φ(r) EM potential (fixed during Dirac solve) */
static double phi_new[NMAX];   /* φ(r) from Maxwell solve */
static double psi_arr[NMAX];   /* Ψ(r) EM flux (fixed, for GR EM stress-energy) */

/* Parameters */
static double par_m      = 1.0;
static double par_omega  = 0.9;
static double par_lambda = 0.5;   /* Soler coupling (= λ_soler/N) */
static double par_G      = 0.0;   /* Newton's constant (0 = flat) */
static double par_e      = 0.0;   /* EM charge coupling */
static int    par_N      = 2;     /* fermions in filled shell */

/*
 * RHS of the 4 coupled ODEs for Dirac+Einstein.
 * y[0]=α, y[1]=β, y[2]=A, y[3]=T
 *
 * EM potential φ and flux Ψ are read from phi_arr[gi] and psi_arr[gi].
 * The caller is responsible for setting phi_arr[gi] to the appropriate
 * value (interpolated for RK4 midpoints).
 */
static void rhs(double r, const double *y, int gi, double *dy)
{
    double alpha = y[0], beta = y[1], A = y[2], T = y[3];
    double sqrtA = sqrt(A);
    double N = par_N;
    double kappa = N * 0.5;  /* = 1 for N=2 (j=1/2 ground state) */

    /* Scalar bilinear and effective mass */
    double al2 = alpha * alpha, be2 = beta * beta;
    double V = N * T * (al2 - be2) / (r * r);
    double m_eff = par_m - par_lambda * V;

    /* Gauge-shifted frequency: ω_eff = ω - eφ (= ω + eV since φ = -V) */
    double omega_eff = par_omega - par_e * phi_arr[gi];

    /* Dirac equations */
    double dalpha = ((kappa / r) * alpha
                     - (omega_eff * T + m_eff) * beta) / sqrtA;
    double dbeta  = ((omega_eff * T - m_eff) * alpha
                     - (kappa / r) * beta) / sqrtA;
    dy[0] = dalpha;
    dy[1] = dbeta;

    /* EM flux contribution: GΨ²/r² (from previous SCF iteration) */
    double Psi = psi_arr[gi];
    double em_term = par_G * Psi * Psi / (r * r);

    /* Einstein equations */
    double G8piN = 8.0 * M_PI * par_G * N;
    double Vsq = V * V;

    /* rA' = 1-A - 8πGN(ω-eφ)T²(α²+β²) - GΨ²/r² - 4πGr²λV² */
    dy[2] = (1.0 - A
             - G8piN * omega_eff * T * T * (al2 + be2)
             - em_term
             - 4.0 * M_PI * par_G * r * r * par_lambda * Vsq) / r;

    /* 2rA(T'/T) = A-1 - 8πGNT√A(αβ'-βα') + GΨ²/r² + 4πGr²λV² */
    double cross = alpha * dbeta - beta * dalpha;
    dy[3] = T * (A - 1.0
                 - G8piN * T * sqrtA * cross
                 + em_term
                 + 4.0 * M_PI * par_G * r * r * par_lambda * Vsq)
            / (2.0 * r * A);
}

/*
 * RK4 integration from r_start to rmax.
 * Evolves 4 variables (α, β, A, T) with φ read from phi_arr[].
 * For RK4 midpoints, φ is interpolated between grid points.
 */
static int integrate(double alpha1, int Ngrid, double rmax)
{
    double r_start = 1e-5;
    double h = (rmax - r_start) / Ngrid;
    double T0 = 1.0;
    int N = par_N;

    /* Taylor expansion coefficients for Dirac */
    double a1sq = alpha1 * alpha1;
    double V0 = N * T0 * a1sq;
    double m_eff0 = par_m - par_lambda * V0;

    /* ω_eff at origin (uses phi_arr[0]) */
    double omega_eff0 = par_omega - par_e * phi_arr[0];
    double beta2 = (omega_eff0 * T0 - m_eff0) * alpha1 / (N + 1.0);

    /* Metric corrections (gravity only) */
    double A2 = 0, T2 = 0;
    if (par_G > 0) {
        A2 = -(8.0 * M_PI * par_G * N * omega_eff0 * T0 * T0 * a1sq
               + 4.0 * M_PI * par_G * par_lambda * V0 * V0) / 3.0;
        double term_cross = 8.0 * M_PI * par_G * N * T0 * alpha1 * beta2;
        double term_V = 4.0 * M_PI * par_G * par_lambda * V0 * V0;
        T2 = (A2 - term_cross + term_V) / 4.0;
    }

    double rs = r_start;
    al_arr[0] = alpha1 * rs;
    be_arr[0] = beta2 * rs * rs;
    A_arr[0]  = 1.0 + A2 * rs * rs;
    T_arr[0]  = T0 + T2 * rs * rs;

    for (int i = 0; i < Ngrid; i++) {
        double r = r_start + i * h;
        double y[NEQN] = {al_arr[i], be_arr[i], A_arr[i], T_arr[i]};
        double k1[NEQN], k2[NEQN], k3[NEQN], k4[NEQN], yt[NEQN];

        /* Save phi/psi for midpoint interpolation */
        double phi_save = phi_arr[i];
        double psi_save = psi_arr[i];
        double phi_mid = 0.5 * (phi_arr[i] + phi_arr[i + 1]);
        double psi_mid = 0.5 * (psi_arr[i] + psi_arr[i + 1]);

        /* k1: at (r, i) */
        rhs(r, y, i, k1);

        /* k2, k3: at (r+h/2) with interpolated phi */
        phi_arr[i] = phi_mid;
        psi_arr[i] = psi_mid;

        for (int j = 0; j < NEQN; j++) yt[j] = y[j] + 0.5 * h * k1[j];
        rhs(r + 0.5 * h, yt, i, k2);

        for (int j = 0; j < NEQN; j++) yt[j] = y[j] + 0.5 * h * k2[j];
        rhs(r + 0.5 * h, yt, i, k3);

        /* Restore phi/psi */
        phi_arr[i] = phi_save;
        psi_arr[i] = psi_save;

        /* k4: at (r+h, i+1) */
        for (int j = 0; j < NEQN; j++) yt[j] = y[j] + h * k3[j];
        rhs(r + h, yt, i + 1, k4);

        for (int j = 0; j < NEQN; j++)
            yt[j] = y[j] + (h / 6.0) * (k1[j] + 2*k2[j] + 2*k3[j] + k4[j]);

        al_arr[i+1] = yt[0];
        be_arr[i+1] = yt[1];
        A_arr[i+1]  = yt[2];
        T_arr[i+1]  = yt[3];

        /* Bail on divergence */
        if (fabs(al_arr[i+1]) > 1e10 || fabs(be_arr[i+1]) > 1e10 ||
            A_arr[i+1] < 1e-6 || T_arr[i+1] < 1e-6 || T_arr[i+1] > 1e6)
            return i + 1;
    }
    return Ngrid;
}

/*
 * Classify solution: -1 if f=α/r crosses zero in physical region
 * (α₁ too large), +1 if f stays positive (α₁ too small).
 */
static int classify(int nv, double rmax)
{
    double r_start = 1e-5;
    double h = (rmax - r_start) / nv;

    /* Find minimum of f²+g² = (α/r)²+(β/r)² */
    double amp_min = 1e30;
    int imin = nv;
    for (int i = 1; i <= nv; i++) {
        double r = r_start + i * h;
        double f = al_arr[i] / r, g = be_arr[i] / r;
        double amp = f * f + g * g;
        if (amp < amp_min) {
            amp_min = amp;
            imin = i;
        }
    }

    /* Check for α crossing zero in physical region */
    for (int i = 1; i <= imin; i++) {
        if (al_arr[i] < 0) return -1;
    }
    return +1;
}

/*
 * Find α₁ by bisection (with current fixed φ).
 */
static double find_alpha1(int Ngrid, double rmax)
{
    double last_pos = -1, first_neg = -1;
    int found_pos = 0;

    for (double a_scan = 0.001; a_scan <= 30.0; a_scan *= 1.01) {
        int nv = integrate(a_scan, Ngrid, rmax);
        int s = classify(nv, rmax);
        if (s == +1) {
            last_pos = a_scan;
            found_pos = 1;
        } else if (s == -1 && found_pos) {
            first_neg = a_scan;
            break;
        }
    }

    if (last_pos < 0 || first_neg < 0) return -1;

    double a_lo = last_pos, a_hi = first_neg;

    for (int iter = 0; iter < 80; iter++) {
        double a_mid = 0.5 * (a_lo + a_hi);
        int nv = integrate(a_mid, Ngrid, rmax);
        int s = classify(nv, rmax);
        if (s == +1)
            a_lo = a_mid;
        else
            a_hi = a_mid;
    }

    double a1_best = 0.5 * (a_lo + a_hi);
    integrate(a1_best, Ngrid, rmax);

    /* Clip growing tail */
    double r_start = 1e-5;
    double h = (rmax - r_start) / Ngrid;
    double amp_min = 1e30;
    int imin = Ngrid;
    for (int i = 1; i <= Ngrid; i++) {
        double r = r_start + i * h;
        double f = al_arr[i] / r, g = be_arr[i] / r;
        double amp = f * f + g * g;
        if (amp < amp_min) {
            amp_min = amp;
            imin = i;
        }
    }
    for (int i = imin + 1; i <= Ngrid; i++) {
        al_arr[i] = 0;
        be_arr[i] = 0;
    }

    return a1_best;
}

/*
 * Solve Maxwell equation via Green's function.
 *
 * Given spinor (α, β) and metric (A, T), compute the Coulomb self-potential.
 *
 * In flat space (A=T=1):
 *   V(r) = eN[M_enc(r)/r + I_ext(r)]
 *   M_enc(r) = 4π∫₀ʳ (α²+β²) dr'
 *   I_ext(r) = 4π∫ᵣ^∞ (α²+β²)/r' dr'
 *
 * In GR: use curved-space integrals
 *   Ψ(r) = ∫₀ʳ 4πeN(α²+β²)T/√A dr'          (forward)
 *   φ(r) = -∫ᵣ^∞ Ψ(r')/(r'²√A·T) dr'          (backward, φ(∞)=0)
 *
 * Also stores Ψ in psi_arr for EM stress-energy in Einstein equations.
 *
 * Output: phi_new[] (in φ(∞)=0 gauge, φ ≤ 0 for e > 0)
 */
static void solve_maxwell(int Ngrid, double rmax)
{
    double r_start = 1e-5;
    double h = (rmax - r_start) / Ngrid;
    int N = par_N;
    static double q_arr[NMAX];  /* Ψ(r) = enclosed flux */

    /* Forward integration: Ψ(r) = ∫₀ʳ 4πeN(α²+β²)T/√A dr' */
    q_arr[0] = 0;
    for (int i = 1; i <= Ngrid; i++) {
        double r_prev = r_start + (i - 1) * h;
        double r_curr = r_start + i * h;

        double al2p = al_arr[i-1] * al_arr[i-1] + be_arr[i-1] * be_arr[i-1];
        double al2c = al_arr[i] * al_arr[i] + be_arr[i] * be_arr[i];

        double Ap = A_arr[i-1], Tp = T_arr[i-1];
        double Ac = A_arr[i], Tc = T_arr[i];

        double fp = 4.0 * M_PI * par_e * N * al2p * Tp / sqrt(Ap);
        double fc = 4.0 * M_PI * par_e * N * al2c * Tc / sqrt(Ac);

        q_arr[i] = q_arr[i-1] + 0.5 * (fp + fc) * h;
        (void)r_prev; (void)r_curr;
    }

    /* Store Ψ for EM stress-energy */
    for (int i = 0; i <= Ngrid; i++)
        psi_arr[i] = q_arr[i];

    /* Backward integration: φ(r) = -∫ᵣ^∞ Ψ/(r²√A·T) dr' */
    phi_new[Ngrid] = 0;
    for (int i = Ngrid - 1; i >= 0; i--) {
        double r_curr = r_start + i * h;
        double r_next = r_start + (i + 1) * h;

        double Ac = A_arr[i], Tc = T_arr[i];
        double An = A_arr[i+1], Tn = T_arr[i+1];

        double fc = q_arr[i] / (r_curr * r_curr * sqrt(Ac) * Tc);
        double fn = q_arr[i+1] / (r_next * r_next * sqrt(An) * Tn);

        phi_new[i] = phi_new[i+1] - 0.5 * (fc + fn) * h;
    }
}

/*
 * Compute physical observables.
 *
 * Mass decomposition (from A equation integrated):
 *   M_ADM = M_Dirac + E_EM + E_self
 * where:
 *   M_Dirac = 4πN ∫ (ω-eφ)T²(α²+β²) dr
 *   E_EM    = (1/2) ∫ Ψ²/r² dr
 *   E_self  = 2πλ ∫ r²V² dr
 */
static void compute_props(int Ngrid, double rmax,
                          double *Q_out, double *R_out,
                          double *M_ADM_out, double *Phi_c_out,
                          double *E_self_out, double *E_EM_out,
                          double *M_Dirac_out)
{
    double r_start = 1e-5;
    double h = (rmax - r_start) / Ngrid;
    double Q = 0, R2w = 0, Eself = 0, E_em = 0, M_dir = 0;

    for (int i = 0; i <= Ngrid; i++) {
        double r = r_start + i * h;
        double al2 = al_arr[i] * al_arr[i];
        double be2 = be_arr[i] * be_arr[i];

        if (al2 + be2 < 1e-30) continue;

        double A = A_arr[i], T = T_arr[i];
        if (A < 1e-10 || !isfinite(A) || !isfinite(T)) continue;
        double sqrtA = sqrt(A);
        double rho = (al2 + be2) * T / sqrtA;

        /* Soler self-energy: V = N·T·(α²-β²)/r² */
        double V_sc = par_N * T * (al2 - be2) / (r * r);
        double eself = r * r * V_sc * V_sc;

        /* EM field energy: (1/2)Ψ²/r² */
        double Psi = psi_arr[i];
        double eem = Psi * Psi / (r * r);

        /* Dirac mass contribution: 4πN(ω-eφ)T²(α²+β²) */
        double omega_eff = par_omega - par_e * phi_arr[i];
        double mdir = par_N * omega_eff * T * T * (al2 + be2);

        double w;
        if (i == 0 || i == Ngrid)
            w = 1.0 / 3.0;
        else if (i % 2 == 1)
            w = 4.0 / 3.0;
        else
            w = 2.0 / 3.0;

        Q     += w * rho * h;
        R2w   += w * rho * r * r * h;
        Eself += w * eself * h;
        E_em  += w * eem * h;
        M_dir += w * mdir * h;
    }

    *Q_out      = 4.0 * M_PI * Q;
    *R_out      = (Q > 0) ? sqrt(R2w / Q) : 0;
    *E_self_out = 2.0 * M_PI * par_lambda * Eself;
    *E_EM_out   = 0.5 * E_em;
    *M_Dirac_out = 4.0 * M_PI * M_dir;

    /* ADM mass: A → 1 - 2GM/r */
    if (par_G > 0 && Ngrid > 0) {
        double A_end = A_arr[Ngrid];
        *M_ADM_out = (1.0 - A_end) * rmax / (2.0 * par_G);
    } else {
        *M_ADM_out = 0;
    }

    /* Gravitational redshift */
    double T_inf = T_arr[Ngrid];
    *Phi_c_out = 1.0 - (T_inf * T_inf) / (T_arr[0] * T_arr[0]);
}

/*
 * SCF iteration for charged soliton.
 *
 * 1. Init φ = 0
 * 2. Solve Dirac(+Einstein) with fixed φ → spinor (α, β)
 * 3. Solve Maxwell from spinor → φ_new
 * 4. Under-relax: φ = (1-η)φ + η·φ_new
 * 5. Check convergence, repeat
 *
 * Returns converged α₁, or -1 on failure.
 */
static double find_alpha1_scf(int Ngrid, double rmax)
{
    /* If no charge, just do single bisection (no SCF needed) */
    if (par_e == 0)
        return find_alpha1(Ngrid, rmax);

    /* Initialize φ = 0, Ψ = 0 */
    for (int i = 0; i <= Ngrid; i++) {
        phi_arr[i] = 0;
        psi_arr[i] = 0;
    }

    double relax = 0.4;
    int max_iter = 200;
    double alpha1_prev = 0;

    fprintf(stderr, "# SCF iteration: e=%.6f, G=%.6e\n", par_e, par_G);
    fprintf(stderr, "# %4s %14s %14s %14s\n",
            "iter", "alpha1", "phi(0)", "dphi_max");

    for (int iter = 0; iter < max_iter; iter++) {
        /* Solve Dirac(+Einstein) with current φ */
        double alpha1 = find_alpha1(Ngrid, rmax);
        if (alpha1 < 0) {
            fprintf(stderr, "# SCF iter %d: no bracket found\n", iter);
            /* Try reducing relaxation */
            if (relax > 0.1) {
                relax *= 0.5;
                fprintf(stderr, "# Reducing relaxation to %.3f\n", relax);
                continue;
            }
            return -1;
        }

        /* Solve Maxwell: compute φ from spinor */
        solve_maxwell(Ngrid, rmax);

        /* Under-relax: φ = (1-η)φ_old + η·φ_new */
        double dphi_max = 0;
        for (int i = 0; i <= Ngrid; i++) {
            double dphi = phi_new[i] - phi_arr[i];
            if (fabs(dphi) > dphi_max) dphi_max = fabs(dphi);
            phi_arr[i] = (1.0 - relax) * phi_arr[i] + relax * phi_new[i];
        }

        fprintf(stderr, "  %4d %14.8f %14.8e %14.8e\n",
                iter, alpha1, phi_arr[0], dphi_max);

        /* Check convergence: relative α₁ change and absolute φ change */
        if (iter > 3 && fabs(alpha1 - alpha1_prev) < 1e-7 * fabs(alpha1)
            && dphi_max < 1e-6) {
            fprintf(stderr, "# Converged at iter %d\n", iter);
            /* Final solve with converged φ */
            alpha1 = find_alpha1(Ngrid, rmax);
            solve_maxwell(Ngrid, rmax);
            /* Update phi_arr to exact converged value */
            for (int i = 0; i <= Ngrid; i++)
                phi_arr[i] = phi_new[i];
            return alpha1;
        }
        alpha1_prev = alpha1;
    }

    fprintf(stderr, "# WARNING: SCF did not converge in %d iterations\n", max_iter);
    /* Return last result anyway */
    return find_alpha1(Ngrid, rmax);
}

static void run_single(double omega, double rmax, int Ngrid, const char *outdir)
{
    par_omega = omega;
    double kappa = sqrt(par_m * par_m - par_omega * par_omega);

    if (10.0 / kappa > rmax) rmax = 10.0 / kappa;
    if (rmax > 500) rmax = 500;

    printf("# Einstein-Dirac-Maxwell-Soler: omega=%.6f, m=%.6f, "
           "lambda=%.6f, e=%.6f, G=%.6e, N=%d\n",
           par_omega, par_m, par_lambda, par_e, par_G, par_N);
    printf("# rmax=%.1f, Ngrid=%d, h=%.6f, kappa=%.6f\n",
           rmax, Ngrid, (rmax - 1e-5) / Ngrid, kappa);

    double alpha1 = find_alpha1_scf(Ngrid, rmax);
    if (alpha1 < 0) {
        printf("# FAILED: No bracket found for omega=%.6f\n", par_omega);
        return;
    }

    printf("#\n");

    double Q, R, M_ADM, Phi_c, E_self, E_EM, M_Dirac;
    compute_props(Ngrid, rmax, &Q, &R, &M_ADM, &Phi_c,
                  &E_self, &E_EM, &M_Dirac);

    double T_inf = T_arr[Ngrid];
    double omega_phys = par_omega * T_inf;
    double M_total = M_Dirac + E_self + E_EM;
    double Q_charge = psi_arr[Ngrid];  /* total electric charge = Ψ(∞) → eNQ_norm */
    double phi_0 = phi_arr[0];

    printf("# RESULT: EDMS soliton\n");
    printf("#   alpha1     = %.8f  (spinor amplitude at origin)\n", alpha1);
    printf("#   omega      = %.8f  (coordinate freq, T₀=1)\n", par_omega);
    printf("#   omega_phys = %.8f  (physical freq, T_∞=1)\n", omega_phys);
    printf("#   E_bind     = %.8f  (ω_phys - m)\n", omega_phys - par_m);
    printf("#   Q_norm     = %.6f  (per-orbital normalization)\n", Q);
    printf("#   R_rms      = %.6f\n", R);

    printf("#\n# Mass decomposition:\n");
    printf("#   M_Dirac    = %.6f  (4πN∫(ω-eφ)T²(α²+β²)dr)\n", M_Dirac);
    printf("#   E_self     = %.6f  (Soler: 2πλ∫r²V²dr)\n", E_self);
    printf("#   E_EM       = %.6f  (EM: (1/2)∫Ψ²/r²dr)\n", E_EM);
    printf("#   M_total    = %.6f  (M_Dirac + E_self + E_EM)\n", M_total);

    if (par_G > 0) {
        printf("#   M_ADM      = %.8f  ((1-A)r/(2G))\n", M_ADM);
        printf("#   dM/M       = %.6e  ((M_ADM-M_total)/M_total)\n",
               M_total > 0 ? (M_ADM - M_total) / M_total : 0);
        printf("#   Phi_c/c2   = %.8e  (1 - T_∞²/T₀²)\n", Phi_c);
        printf("#   GM/R       = %.6e  (compactness)\n", par_G * M_ADM / R);
    }

    if (par_e != 0) {
        printf("#\n# Electromagnetic:\n");
        printf("#   e          = %.8f  (charge coupling)\n", par_e);
        printf("#   Q_charge   = %.6f  (total charge = Ψ(∞) = eNQ)\n", Q_charge);
        printf("#   Q_expected = %.6f  (eNQ_norm = %.4f×%d×%.4f)\n",
               par_e * par_N * Q, par_e, par_N, Q);
        printf("#   phi(0)     = %.8e  (EM potential at center, φ(∞)=0 gauge)\n", phi_0);
        printf("#   V(0)       = %.8e  (Coulomb: -phi(0), should be > 0)\n", -phi_0);
        printf("#   Psi(rmax)  = %.8e  (EM flux, should → eNQ)\n", Q_charge);
    }

    printf("#\n# Metric at boundary:\n");
    printf("#   A(rmax)    = %.10f\n", A_arr[Ngrid]);
    printf("#   T(rmax)    = %.10f\n", T_inf);
    printf("#   T(0)       = %.10f\n", T_arr[0]);

    /* Save profile */
    if (outdir) {
        char fname[256];
        snprintf(fname, sizeof(fname), "%s/edm_om%.3f_e%.3f_G%.1e.dat",
                 outdir, par_omega, par_e, par_G);
        FILE *fp = fopen(fname, "w");
        if (fp) {
            double r_start = 1e-5;
            double h = (rmax - r_start) / Ngrid;
            fprintf(fp, "# EDMS: omega=%.6f m=%.6f lambda=%.6f e=%.6f "
                    "G=%.6e N=%d alpha1=%.8f\n",
                    par_omega, par_m, par_lambda, par_e, par_G, par_N, alpha1);
            fprintf(fp, "# r  alpha  beta  A  T  phi  Psi  rho\n");
            int step = (Ngrid > 5000) ? Ngrid / 5000 : 1;
            for (int i = 0; i <= Ngrid; i += step) {
                double r = r_start + i * h;
                double al = al_arr[i], be = be_arr[i];
                fprintf(fp, "%.6f %.10e %.10e %.10e %.10e %.10e %.10e %.10e\n",
                        r, al, be, A_arr[i], T_arr[i],
                        phi_arr[i], psi_arr[i],
                        al * al + be * be);
            }
            fclose(fp);
            printf("#   Profile: %s\n", fname);
        }
    }
}

static void run_scan(double rmax, int Ngrid, const char *outdir)
{
    double omega_vals[] = {0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.99};
    int nvals = sizeof(omega_vals) / sizeof(omega_vals[0]);

    printf("# EDMS scan: m=%.6f, lambda=%.6f, e=%.6f, G=%.6e, N=%d\n",
           par_m, par_lambda, par_e, par_G, par_N);
    printf("# %10s  %12s  %12s  %8s  %12s  %12s  %12s  %12s  %12s\n",
           "omega", "alpha1", "Q_norm", "R_rms",
           "M_Dirac", "E_self", "E_EM", "M_total", "M_ADM");

    for (int k = 0; k < nvals; k++) {
        par_omega = omega_vals[k];
        double rmax_k = rmax;
        double kappa = sqrt(par_m * par_m - par_omega * par_omega);
        if (kappa > 0 && 10.0 / kappa > rmax_k)
            rmax_k = 10.0 / kappa;
        if (rmax_k > 500) rmax_k = 500;

        double alpha1 = find_alpha1_scf(Ngrid, rmax_k);
        if (alpha1 < 0) {
            printf("  %10.4f  FAILED\n", par_omega);
            continue;
        }

        double Q, R, M_ADM, Phi_c, E_self, E_EM, M_Dirac;
        compute_props(Ngrid, rmax_k, &Q, &R, &M_ADM, &Phi_c,
                      &E_self, &E_EM, &M_Dirac);
        double M_total = M_Dirac + E_self + E_EM;

        printf("  %10.4f  %12.6f  %12.6f  %8.4f  %12.6f  %12.6f  "
               "%12.6e  %12.6f  %12.6e\n",
               par_omega, alpha1, Q, R,
               M_Dirac, E_self, E_EM, M_total, M_ADM);
    }
}

static void run_scan_e(double omega, double rmax, int Ngrid, const char *outdir)
{
    double e_vals[] = {0.0, 0.001, 0.005, 0.01, 0.02, 0.05, 0.08, 0.1, 0.15, 0.2, 0.3};
    int nvals = sizeof(e_vals) / sizeof(e_vals[0]);

    printf("# EDMS scan_e: omega=%.6f, m=%.6f, lambda=%.6f, G=%.6e, N=%d\n",
           omega, par_m, par_lambda, par_G, par_N);
    printf("# %10s  %12s  %12s  %8s  %12s  %12s  %12s  %12s\n",
           "e", "alpha1", "Q_norm", "R_rms",
           "M_total", "E_EM", "phi(0)", "Q_charge");

    for (int k = 0; k < nvals; k++) {
        par_e = e_vals[k];
        par_omega = omega;
        double rmax_k = rmax;
        double kappa = sqrt(par_m * par_m - par_omega * par_omega);
        if (kappa > 0 && 10.0 / kappa > rmax_k)
            rmax_k = 10.0 / kappa;
        if (rmax_k > 500) rmax_k = 500;

        double alpha1 = find_alpha1_scf(Ngrid, rmax_k);
        if (alpha1 < 0) {
            printf("  %10.4f  FAILED (charge too large or no bracket)\n", par_e);
            continue;
        }

        double Q, R, M_ADM, Phi_c, E_self, E_EM, M_Dirac;
        compute_props(Ngrid, rmax_k, &Q, &R, &M_ADM, &Phi_c,
                      &E_self, &E_EM, &M_Dirac);
        double M_total = M_Dirac + E_self + E_EM;
        double Q_charge = psi_arr[Ngrid];

        printf("  %10.4f  %12.6f  %12.6f  %8.4f  %12.6f  %12.6e  "
               "%12.6e  %12.6f\n",
               par_e, alpha1, Q, R,
               M_total, E_EM, phi_arr[0], Q_charge);
        fflush(stdout);
    }
}

static void usage(const char *prog)
{
    fprintf(stderr,
            "Usage: %s [-omega ω] [-m mass] [-lambda λ] [-e charge]\n"
            "          [-G G] [-Nferm N] [-rmax R] [-N n]\n"
            "          [-scan] [-scan_e] [-outdir dir]\n",
            prog);
}

int main(int argc, char *argv[])
{
    double omega = 0.9;
    double rmax = 50.0;
    int Ngrid = 100000;
    int scan = 0, scan_e = 0;
    const char *outdir = NULL;

    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-omega") == 0 && i + 1 < argc)
            omega = atof(argv[++i]);
        else if (strcmp(argv[i], "-m") == 0 && i + 1 < argc)
            par_m = atof(argv[++i]);
        else if (strcmp(argv[i], "-lambda") == 0 && i + 1 < argc)
            par_lambda = atof(argv[++i]);
        else if (strcmp(argv[i], "-e") == 0 && i + 1 < argc)
            par_e = atof(argv[++i]);
        else if (strcmp(argv[i], "-G") == 0 && i + 1 < argc)
            par_G = atof(argv[++i]);
        else if (strcmp(argv[i], "-Nferm") == 0 && i + 1 < argc)
            par_N = atoi(argv[++i]);
        else if (strcmp(argv[i], "-rmax") == 0 && i + 1 < argc)
            rmax = atof(argv[++i]);
        else if (strcmp(argv[i], "-N") == 0 && i + 1 < argc)
            Ngrid = atoi(argv[++i]);
        else if (strcmp(argv[i], "-scan") == 0)
            scan = 1;
        else if (strcmp(argv[i], "-scan_e") == 0)
            scan_e = 1;
        else if (strcmp(argv[i], "-outdir") == 0 && i + 1 < argc)
            outdir = argv[++i];
        else {
            usage(argv[0]);
            return 1;
        }
    }

    if (Ngrid > NMAX - 1) Ngrid = NMAX - 1;

    if (scan) {
        run_scan(rmax, Ngrid, outdir);
    } else if (scan_e) {
        run_scan_e(omega, rmax, Ngrid, outdir);
    } else {
        run_single(omega, rmax, Ngrid, outdir);
    }

    return 0;
}
