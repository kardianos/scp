/*
 * boson_star.c — Newtonian boson star (Schrödinger-Poisson system)
 *
 * Demonstrates self-confinement by gravity: the scalar field's own
 * gravitational potential traps it. No self-interaction needed.
 *
 * System (natural units ħ = μ = 1):
 *   u'' + 2u'/r + 2(E - Φ)u = 0      (Schrödinger, fixed Φ)
 *   Φ'' + 2Φ'/r = 4πG u²             (Poisson)
 *
 * Method: self-consistent field (SCF) iteration:
 *   1. Bootstrap Φ from Gaussian u profile
 *   2. Shoot Schrödinger for eigenvalue E (linear, fixed Φ)
 *   3. Clip growing tail of u beyond its minimum |u|
 *   4. Solve Poisson for new Φ from u²
 *   5. Under-relax: Φ ← α·Φ_new + (1-α)·Φ_old
 *   6. Repeat until E converges
 *
 * Scale invariance: E ∝ uc, R ∝ 1/√uc. Solutions exist for all uc > 0.
 *
 * Usage:
 *   boson_star -uc 1.0              (single central amplitude)
 *   boson_star -scan                (scan uc from 0.01 to 10)
 */

#define _DEFAULT_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define NR_MAX 100001

static double u_arr[NR_MAX], v_arr[NR_MAX];
static double phi_arr[NR_MAX];
static double phi_new[NR_MAX];

/*
 * Shoot Schrödinger u'' + 2u'/r + 2(E-Φ)u = 0 with FIXED Φ.
 * Returns last valid index (N if complete, < N if |u| > 1e10).
 */
static int shoot(double uc, double E, int N, double rmax)
{
    double h = rmax / N;
    double eps0 = E - phi_arr[0];
    double u2c = -eps0 * uc / 3.0;  /* u'' = -2(E-Φ₀)u/3 at r=0 */

    u_arr[0] = uc;        v_arr[0] = 0;
    u_arr[1] = uc + u2c * h * h;
    v_arr[1] = 2.0 * u2c * h;

    for (int i = 1; i < N; i++) {
        double r  = i * h;
        double re = r + h;
        double rm = r + 0.5 * h;
        double u = u_arr[i], v = v_arr[i];

        double Phi_i = phi_arr[i];
        double Phi_m = 0.5 * (phi_arr[i] + phi_arr[i + 1]);
        double Phi_e = phi_arr[i + 1];

        double k1u = v;
        double k1v = -2.0 * v / r - 2.0 * (E - Phi_i) * u;

        double um = u + 0.5 * h * k1u;
        double vm = v + 0.5 * h * k1v;
        double k2u = vm;
        double k2v = -2.0 * vm / rm - 2.0 * (E - Phi_m) * um;

        um = u + 0.5 * h * k2u;
        vm = v + 0.5 * h * k2v;
        double k3u = vm;
        double k3v = -2.0 * vm / rm - 2.0 * (E - Phi_m) * um;

        um = u + h * k3u;
        vm = v + h * k3v;
        double k4u = vm;
        double k4v = -2.0 * vm / re - 2.0 * (E - Phi_e) * um;

        u_arr[i + 1] = u + (h / 6.0) * (k1u + 2*k2u + 2*k3u + k4u);
        v_arr[i + 1] = v + (h / 6.0) * (k1v + 2*k2v + 2*k3v + k4v);

        if (fabs(u_arr[i + 1]) > 1e10) return i + 1;
    }
    return N;
}

/*
 * Classify shooting result for ground state (no nodes):
 *   -1  = overshoot (u crosses zero — E too negative)
 *   +1  = undershoot (u stays positive, diverges — E too close to 0)
 */
static int classify(int nv)
{
    for (int i = 1; i <= nv; i++)
        if (u_arr[i] < 0) return -1;
    return +1;
}

/*
 * Find ground-state eigenvalue E by scanning + bisection with fixed Φ.
 * Bound state has Φ(0) < E < 0.
 *
 * As E increases from Φ(0) toward 0:
 *   E < E₁: 0 nodes → undershoot (s = +1)
 *   E > E₁: ≥1 node → overshoot  (s = -1)
 * Ground state is at the FIRST undershoot→overshoot transition.
 *
 * Bisection: E_lo (undershoot) < E₁ < E_hi (overshoot).
 */
static double find_E(double uc, int N, double rmax, double E_hint)
{
    double Phi0 = phi_arr[0];
    if (Phi0 >= 0) return 0;

    double E_bot = Phi0 + 0.005 * fabs(Phi0);  /* just above well bottom */
    double E_top = -1e-10;                      /* near zero */
    double E_lo, E_hi;

    /* If we have a hint from previous iteration, try tight bracket first */
    if (E_hint < -1e-12 && E_hint > Phi0) {
        double dE = 0.15 * fabs(E_hint) + 0.01;
        double try_lo = E_hint - dE;
        double try_hi = E_hint + dE;
        if (try_lo < Phi0) try_lo = E_bot;
        if (try_hi > -1e-10) try_hi = E_top;

        int n1 = shoot(uc, try_lo, N, rmax);
        int s1 = classify(n1);
        int n2 = shoot(uc, try_hi, N, rmax);
        int s2 = classify(n2);

        if (s1 == +1 && s2 == -1) {
            E_lo = try_lo;
            E_hi = try_hi;
            goto bisect;
        }
    }

    /* Full scan: find first +1 → -1 transition (ground state) */
    {
        int nscan = 300;
        int prev_s = 0;
        double prev_E = E_bot;
        int found = 0;
        E_lo = E_bot;
        E_hi = E_top;

        for (int j = 0; j <= nscan; j++) {
            double E = E_bot + j * (E_top - E_bot) / nscan;
            int nv = shoot(uc, E, N, rmax);
            int sj = classify(nv);

            if (j > 0 && prev_s == +1 && sj == -1) {
                E_lo = prev_E;   /* undershoot */
                E_hi = E;        /* overshoot */
                found = 1;
                break;  /* first transition = ground state */
            }
            prev_E = E;
            prev_s = sj;
        }

        if (!found) return 0;
    }

bisect:
    for (int iter = 0; iter < 200; iter++) {
        double E_mid = 0.5 * (E_lo + E_hi);
        if (fabs(E_hi - E_lo) < 1e-13 * fabs(E_mid)) break;

        int nv = shoot(uc, E_mid, N, rmax);
        int s = classify(nv);

        if (s > 0)  E_lo = E_mid;   /* undershoot → raise lower bound */
        else        E_hi = E_mid;   /* overshoot  → lower upper bound */
    }

    double E_best = 0.5 * (E_lo + E_hi);
    int nv = shoot(uc, E_best, N, rmax);
    /* Find crossover (minimum |u| in valid range) and clip growing tail.
     * This removes both the exponentially growing mode AND stale data. */
    double umin = fabs(u_arr[0]);
    int imin = 0;
    for (int i = 1; i <= nv; i++) {
        if (fabs(u_arr[i]) < umin) {
            umin = fabs(u_arr[i]);
            imin = i;
        }
    }
    for (int i = imin + 1; i <= N; i++) {
        u_arr[i] = 0;
        v_arr[i] = 0;
    }
    return E_best;
}

/*
 * Clip growing tail: find global minimum of |u| and zero out beyond it.
 * For well-converged eigenvalue, this removes only the exponentially
 * growing numerical contamination.
 */
static void clip_tail(int N)
{
    double umin = fabs(u_arr[0]);
    int imin = 0;
    for (int i = 1; i <= N; i++) {
        if (fabs(u_arr[i]) < umin) {
            umin = fabs(u_arr[i]);
            imin = i;
        }
    }
    for (int i = imin + 1; i <= N; i++) {
        u_arr[i] = 0;
        v_arr[i] = 0;
    }
}

/*
 * Solve Poisson for Φ(r) from u(r) via Green's function:
 *   Φ(r) = -(1/r) M(r) - I(r)
 * where M(r) = 4πG ∫₀ʳ u²r'² dr', I(r) = 4πG ∫ᵣ^∞ u²r' dr'
 */
static void solve_poisson(double G, int N, double rmax)
{
    double h = rmax / N;
    static double M_enc[NR_MAX], I_ext[NR_MAX];

    M_enc[0] = 0;
    for (int i = 1; i <= N; i++) {
        double r = i * h, rp = (i - 1) * h;
        double u2p = u_arr[i - 1] * u_arr[i - 1];
        double u2c = u_arr[i] * u_arr[i];
        M_enc[i] = M_enc[i - 1]
                  + 0.5 * 4.0 * M_PI * G * (u2p * rp * rp + u2c * r * r) * h;
    }

    I_ext[N] = 0;
    for (int i = N - 1; i >= 0; i--) {
        double r = i * h, rn = (i + 1) * h;
        double u2c = u_arr[i] * u_arr[i];
        double u2n = u_arr[i + 1] * u_arr[i + 1];
        I_ext[i] = I_ext[i + 1]
                  + 0.5 * 4.0 * M_PI * G * (u2c * r + u2n * rn) * h;
    }

    phi_new[0] = -I_ext[0];
    for (int i = 1; i <= N; i++) {
        double r = i * h;
        phi_new[i] = -M_enc[i] / r - I_ext[i];
    }
}

static void run_single(double uc, double rmax, int N, double G,
                        double relax, const char *outdir)
{
    double h = rmax / N;

    printf("# Boson star: SCF Schrödinger-Poisson\n");
    printf("# uc=%.6f  G=%.6f  rmax=%.1f  N=%d  h=%.6f  alpha=%.2f\n",
           uc, G, rmax, N, h, relax);

    /* Bootstrap: Gaussian u → Poisson Φ.
     * sigma ~ 1/√uc gives Φ(0) ∝ uc, matching scaling. */
    double sigma = 1.0 / sqrt(uc);
    if (sigma < 0.3) sigma = 0.3;
    if (sigma > rmax / 5.0) sigma = rmax / 5.0;
    for (int i = 0; i <= N; i++) {
        double r = i * h;
        u_arr[i] = uc * exp(-r * r / (2.0 * sigma * sigma));
    }
    solve_poisson(G, N, rmax);
    for (int i = 0; i <= N; i++)
        phi_arr[i] = phi_new[i];

    printf("# Bootstrap: sigma=%.2f  Phi(0)=%.6f\n", sigma, phi_arr[0]);
    printf("# %4s %14s %14s %12s %12s\n",
           "iter", "E", "Phi(0)", "M_grav", "R_rms");

    double E = 0, E_prev = 1;
    int max_iter = 150;

    for (int iter = 0; iter < max_iter; iter++) {
        E = find_E(uc, N, rmax, E_prev != 1 ? E_prev : 0);
        if (E == 0) {
            printf("# Failed to find eigenvalue at iter %d\n", iter);
            printf("# Phi(0)=%.8f  Phi(rmax)=%.8f\n", phi_arr[0], phi_arr[N]);
            return;
        }

        /* Clip growing tail */
        clip_tail(N);

        /* Solve Poisson for new Φ */
        solve_poisson(G, N, rmax);

        /* Under-relax */
        for (int i = 0; i <= N; i++)
            phi_arr[i] = relax * phi_new[i] + (1.0 - relax) * phi_arr[i];

        /* Compute properties */
        double M_grav = 0, norm = 0, mom4 = 0;
        for (int i = 0; i <= N; i++) {
            double r = i * h;
            double u2 = u_arr[i] * u_arr[i];
            double w = (i == 0 || i == N) ? 0.5 : 1.0;
            M_grav += w * 4.0 * M_PI * G * u2 * r * r * h;
            norm += w * u2 * r * r * h;
            mom4 += w * u2 * r * r * r * r * h;
        }
        double Rrms = (norm > 0) ? sqrt(mom4 / norm) : 0;

        printf("  %4d %14.8f %14.8f %12.6f %12.4f\n",
               iter, E, phi_arr[0], M_grav, Rrms);
        fflush(stdout);

        if (iter > 3 && fabs(E - E_prev) < 1e-11 * fabs(E)) {
            printf("# Converged at iter %d\n", iter);
            break;
        }
        E_prev = E;
    }

    /* Final properties */
    double M_grav = 0, Npart = 0, norm = 0, mom4 = 0;
    for (int i = 0; i <= N; i++) {
        double r = i * h;
        double u2 = u_arr[i] * u_arr[i];
        double w = (i == 0 || i == N) ? 0.5 : 1.0;
        M_grav += w * 4.0 * M_PI * G * u2 * r * r * h;
        Npart += w * 4.0 * M_PI * u2 * r * r * h;
        norm += w * u2 * r * r * h;
        mom4 += w * u2 * r * r * r * r * h;
    }
    double Rrms = (norm > 0) ? sqrt(mom4 / norm) : 0;
    double kappa = (E < 0) ? sqrt(-2.0 * E) : 0;

    printf("\n# RESULT: Newtonian boson star\n");
    printf("#   u(0)      = %.8f\n", uc);
    printf("#   E         = %.10f  (binding energy eigenvalue)\n", E);
    printf("#   Phi(0)    = %.10f  (central gravitational potential)\n",
           phi_arr[0]);
    printf("#   M_grav    = %.8f  (gravitational mass)\n", M_grav);
    printf("#   N_part    = %.8f  (particle number)\n", Npart);
    printf("#   R_rms     = %.6f\n", Rrms);
    printf("#   kappa     = %.6f  (decay const sqrt(-2E))\n", kappa);
    printf("#   E/Phi(0)  = %.6f  (ratio, virial ~0.5)\n",
           fabs(phi_arr[0]) > 1e-30 ? E / phi_arr[0] : 0);

    if (E < 0 && M_grav > 0) {
        printf("#\n# POSITIVE: Gravity provides self-confinement.\n");
        printf("# The scalar field IS bound by its own gravitational potential.\n");
        printf("# Without gravity (G=0), no bound state exists.\n");
    } else {
        printf("#\n# NULL: No self-consistent bound state found.\n");
    }

    /* Output profile */
    char fname[512];
    snprintf(fname, sizeof(fname), "%s/boson_star_uc%.3f.dat", outdir, uc);
    FILE *fp = fopen(fname, "w");
    if (fp) {
        fprintf(fp, "# r  u  u'  Phi\n");
        int step = (N > 2000) ? N / 2000 : 1;
        for (int i = 0; i <= N; i += step)
            fprintf(fp, "%.6f %.12e %.12e %.12e\n",
                    i * h, u_arr[i], v_arr[i], phi_arr[i]);
        fclose(fp);
        printf("#   Profile: %s\n", fname);
    }
}

static void run_scan(double rmax, int N, double G, double relax)
{
    printf("# Boson star scan: SCF Schrödinger-Poisson, G=%.6f\n", G);
    printf("# %10s %14s %14s %12s %10s\n",
           "u_c", "E", "Phi(0)", "M_grav", "R_rms");

    double uc_vals[] = {0.01, 0.02, 0.05, 0.1, 0.2, 0.5,
                        1.0, 2.0, 5.0, 10.0};
    int nuc = sizeof(uc_vals) / sizeof(uc_vals[0]);

    for (int j = 0; j < nuc; j++) {
        double uc = uc_vals[j];

        /* Scale rmax with 1/sqrt(uc) for diffuse solutions */
        double rm = rmax / sqrt(uc);
        if (rm < rmax) rm = rmax;
        if (rm > 500) rm = 500;
        double h = rm / N;

        /* Bootstrap */
        double sigma = 1.0 / sqrt(uc);
        if (sigma < 0.3) sigma = 0.3;
        if (sigma > rm / 5.0) sigma = rm / 5.0;

        for (int i = 0; i <= N; i++) {
            double r = i * h;
            u_arr[i] = uc * exp(-r * r / (2.0 * sigma * sigma));
        }
        solve_poisson(G, N, rm);
        for (int i = 0; i <= N; i++) phi_arr[i] = phi_new[i];

        /* SCF iteration */
        double E = 0, E_prev = 1;
        int converged = 0;
        for (int iter = 0; iter < 200; iter++) {
            E = find_E(uc, N, rm, E_prev != 1 ? E_prev : 0);
            if (E == 0) break;

            clip_tail(N);
            solve_poisson(G, N, rm);
            for (int i = 0; i <= N; i++)
                phi_arr[i] = relax * phi_new[i] + (1.0 - relax) * phi_arr[i];

            if (iter > 3 && fabs(E - E_prev) < 1e-10 * fabs(E)) {
                converged = 1;
                break;
            }
            E_prev = E;
        }

        if (E == 0 || !converged) {
            printf("  %10.4f   FAILED\n", uc);
            continue;
        }

        double M_grav = 0, norm = 0, mom4 = 0;
        for (int i = 0; i <= N; i++) {
            double r = i * h;
            double u2 = u_arr[i] * u_arr[i];
            double w = (i == 0 || i == N) ? 0.5 : 1.0;
            M_grav += w * 4.0 * M_PI * G * u2 * r * r * h;
            norm += w * u2 * r * r * h;
            mom4 += w * u2 * r * r * r * r * h;
        }
        double Rrms = (norm > 0) ? sqrt(mom4 / norm) : 0;

        printf("  %10.4f %14.8f %14.8f %12.6f %10.4f  BOUND\n",
               uc, E, phi_arr[0], M_grav, Rrms);
        fflush(stdout);
    }
}

int main(int argc, char **argv)
{
    double uc = 1.0;
    double rmax = 50.0;
    int N = 20000;
    double G = 1.0;
    double relax = 0.2;
    char outdir[256] = "data";
    int do_scan = 0;

    for (int i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "-uc") && i + 1 < argc)
            uc = atof(argv[++i]);
        else if (!strcmp(argv[i], "-rmax") && i + 1 < argc)
            rmax = atof(argv[++i]);
        else if (!strcmp(argv[i], "-N") && i + 1 < argc) {
            N = atoi(argv[++i]);
            if (N >= NR_MAX) N = NR_MAX - 1;
        }
        else if (!strcmp(argv[i], "-G") && i + 1 < argc)
            G = atof(argv[++i]);
        else if (!strcmp(argv[i], "-alpha") && i + 1 < argc)
            relax = atof(argv[++i]);
        else if (!strcmp(argv[i], "-scan"))
            do_scan = 1;
        else if (!strcmp(argv[i], "-outdir") && i + 1 < argc)
            strncpy(outdir, argv[++i], sizeof(outdir) - 1);
    }

    if (do_scan)
        run_scan(rmax, N, G, relax);
    else
        run_single(uc, rmax, N, G, relax, outdir);

    return 0;
}
