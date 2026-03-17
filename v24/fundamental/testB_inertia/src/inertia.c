/*
 * inertia.c — Test B: Oscillon Deformation Under Acceleration
 *
 * Based on v21/src/triad1d.c (1D three-body oscillon).
 *
 * Method: Equilibrate oscillon for t_equil, then apply a constant external
 * force F by injecting momentum dP = F*dt at each timestep. The momentum
 * injection is done by boosting field velocities:
 *   vel[a][i] += (F*dt / (2*E_grad)) * (phi[a][i+1] - phi[a][i-1]) / (2*dx)
 * where E_grad = sum_a ∫ (∂φ_a/∂x)² dx is the gradient energy (the natural
 * "inertia" for momentum coupling). This ensures dP/dt = F exactly.
 *
 * Alternatively, uses a linear gravitational potential V_ext = -(F/2)*x*Σφ²
 * implemented as acc[a][i] += F*x*φ_a (controlled by -method flag).
 *
 * Measures:
 *   - Center of mass trajectory x_c(t), velocity v_c(t)
 *   - Leading/trailing width asymmetry σ+ vs σ-
 *   - Energy density profile in co-moving frame
 *   - Deformation δρ = ρ_accel - ρ_rest
 *   - Linearity of deformation in F
 *
 * Compile: gcc -O3 -Wall -o inertia src/inertia.c -lm
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Parameters from PROPOSAL: mu=-20, kappa=20, m=1.0 (3D-stable oscillon) */
static double mu     = -20.0;
static double kappa  = 20.0;
static double mass   = 1.0;
static double A_init = 0.8;
static double sigma  = 3.0;
static int    Nx     = 4000;
static double xmax   = 100.0;
static double t_equil = 5000.0;
static double t_accel = 3000.0;
static double F_ext  = 0.01;
static int    method = 1;   /* 1=momentum injection, 2=linear potential */
static char   outdir[512] = "v24/fundamental/testB_inertia/data";

static void parse_args(int argc, char **argv)
{
    for (int i = 1; i < argc - 1; i += 2) {
        if      (!strcmp(argv[i], "-mu"))      mu      = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-kappa"))   kappa   = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-mass"))    mass    = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-A"))       A_init  = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-sigma"))   sigma   = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-Nx"))      Nx      = atoi(argv[i+1]);
        else if (!strcmp(argv[i], "-xmax"))    xmax    = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-t_equil")) t_equil = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-t_accel")) t_accel = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-F"))       F_ext   = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-method"))  method  = atoi(argv[i+1]);
        else if (!strcmp(argv[i], "-o"))       strncpy(outdir, argv[i+1], sizeof(outdir)-1);
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

/* Compute energy density at grid point i */
static double energy_density(double *phi[3], double *vel[3], double dx, double m2, int i)
{
    double e = 0.0;
    for (int a = 0; a < 3; a++) {
        e += 0.5 * vel[a][i] * vel[a][i];
        double dp = (phi[a][i+1] - phi[a][i-1]) / (2.0 * dx);
        e += 0.5 * dp * dp;
        e += 0.5 * m2 * phi[a][i] * phi[a][i];
    }
    double P = phi[0][i] * phi[1][i] * phi[2][i];
    double P2 = P * P;
    e += 0.5 * mu * P2 / (1.0 + kappa * P2);
    return e;
}

/* Find center of energy */
static double find_center(double *phi[3], double *vel[3], double dx, double m2,
                          int N, double xmin, double *E_total)
{
    double xc = 0.0, Etot = 0.0;
    for (int i = 1; i < N - 1; i++) {
        double x = xmin + i * dx;
        double e = energy_density(phi, vel, dx, m2, i);
        xc += x * e * dx;
        Etot += e * dx;
    }
    *E_total = Etot;
    return (Etot > 1e-20) ? xc / Etot : 0.0;
}

/* Measure half-widths on leading (+x) and trailing (-x) sides of center */
static void measure_widths(double *phi[3], double *vel[3], double dx, double m2,
                           int N, double xmin, double xc,
                           double *sigma_plus, double *sigma_minus)
{
    double s2_plus = 0, w_plus = 0;
    double s2_minus = 0, w_minus = 0;
    for (int i = 1; i < N - 1; i++) {
        double x = xmin + i * dx;
        double e = energy_density(phi, vel, dx, m2, i);
        double dxc = x - xc;
        if (dxc >= 0) {
            s2_plus += dxc * dxc * e * dx;
            w_plus += e * dx;
        } else {
            s2_minus += dxc * dxc * e * dx;
            w_minus += e * dx;
        }
    }
    *sigma_plus  = (w_plus  > 1e-20) ? sqrt(s2_plus  / w_plus)  : 0;
    *sigma_minus = (w_minus > 1e-20) ? sqrt(s2_minus / w_minus) : 0;
}

/* Compute field momentum P = -sum_a ∫ (∂_t φ_a)(∂_x φ_a) dx */
static double field_momentum(double *phi[3], double *vel[3], double dx, int N)
{
    double P = 0.0;
    for (int a = 0; a < 3; a++)
        for (int i = 1; i < N - 1; i++) {
            double grad = (phi[a][i+1] - phi[a][i-1]) / (2.0 * dx);
            P -= vel[a][i] * grad * dx;
        }
    return P;
}

int main(int argc, char **argv)
{
    parse_args(argc, argv);

    double dx = 2.0 * xmax / (Nx - 1);
    double dx2 = dx * dx;
    double m2 = mass * mass;
    double xmin = -xmax;
    double tfinal = t_equil + t_accel;

    /* CFL */
    double kmax = M_PI / dx;
    double dt = 0.8 * 2.0 / sqrt(kmax * kmax + m2);
    int Nt = (int)(tfinal / dt) + 1;
    int N_equil = (int)(t_equil / dt);

    const char *method_desc = (method == 1) ? "momentum injection" : "linear potential";
    printf("inertia: Oscillon deformation under acceleration\n");
    printf("  mu=%.1f kappa=%.1f mass=%.3f A=%.3f sigma=%.3f\n", mu, kappa, mass, A_init, sigma);
    printf("  Nx=%d xmax=%.1f dx=%.5f dt=%.6f\n", Nx, xmax, dx, dt);
    printf("  t_equil=%.0f t_accel=%.0f F=%.6f method=%d (%s)\n",
           t_equil, t_accel, F_ext, method, method_desc);
    printf("  Nt=%d  N_equil=%d\n", Nt, N_equil);

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
        double x = xmin + i * dx;
        double ax = fabs(x);
        if (ax > x_abs) {
            double f = (ax - x_abs) / (xmax - x_abs);
            damp[i] = 1.0 - 0.98 * f * f;
        } else {
            damp[i] = 1.0;
        }
    }

    /* Initialize: symmetric triad Gaussians centered at x=0 */
    for (int a = 0; a < 3; a++)
        for (int i = 0; i < Nx; i++) {
            double x = xmin + i * dx;
            phi[a][i] = A_init * exp(-x * x / (2.0 * sigma * sigma));
        }

    /* Compute acceleration (no external force in this macro — force handled separately) */
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

    /* Method 2: add linear potential force to acceleration */
    #define ADD_LINEAR_FORCE() do { \
        for (int a = 0; a < 3; a++) \
            for (int i = 1; i < Nx - 1; i++) { \
                double x_i = xmin + i * dx; \
                acc[a][i] += F_ext * x_i * phi[a][i]; \
            } \
    } while(0)

    COMPUTE_ACC();

    /* --- Storage for rest-frame profile snapshot --- */
    int prof_N = 800;       /* profile bins: [-prof_R, +prof_R] */
    double prof_R = 20.0;   /* profile half-width */
    double prof_dx = 2.0 * prof_R / prof_N;
    double *rho_rest = calloc(prof_N, sizeof(double));
    int rest_snap_count = 0;

    /* Time series output */
    char tspath[600];
    snprintf(tspath, sizeof(tspath), "%s/inertia_F%.6f_ts.tsv", outdir, F_ext);
    FILE *fts = fopen(tspath, "w");
    if (!fts) { fprintf(stderr, "Cannot open %s\n", tspath); return 1; }
    fprintf(fts, "time\txc\tvx\tP_field\tE_total\tsigma_plus\tsigma_minus\tasymmetry\tphase\n");

    int rec_every = Nt / 40000;
    if (rec_every < 1) rec_every = 1;
    int print_every = Nt / 40;
    if (print_every < 1) print_every = 1;

    /* Snapshot the rest profile: average over last 500 t.u. of equilibration */
    int snap_start = (int)((t_equil - 500.0) / dt);
    int snap_end   = N_equil;
    int snap_count_target = 200;
    int snap_every = (snap_end - snap_start) / snap_count_target;
    if (snap_every < 1) snap_every = 1;

    /* Accel profile snapshot: average over last 500 t.u. of acceleration */
    double *rho_accel = calloc(prof_N, sizeof(double));
    int accel_snap_count = 0;
    int accel_snap_start = (int)((tfinal - 500.0) / dt);
    int accel_snap_every = (Nt - accel_snap_start) / snap_count_target;
    if (accel_snap_every < 1) accel_snap_every = 1;

    double xc_prev = 0.0;
    double vx_est = 0.0;

    /* Measure rest-frame energy and "charge" Q = ∫Σφ² for method 2 */
    double E_rest = 0, Q_rest = 0;

    for (int n = 0; n <= Nt; n++) {
        double t = n * dt;
        int force_on = (n >= N_equil);
        const char *phase = force_on ? "accel" : "equil";

        /* --- Record time series --- */
        if (n % rec_every == 0) {
            double Etot;
            double xc = find_center(phi, vel, dx, m2, Nx, xmin, &Etot);
            double sp, sm;
            measure_widths(phi, vel, dx, m2, Nx, xmin, xc, &sp, &sm);
            double asym = sp - sm;
            double Pfield = field_momentum(phi, vel, dx, Nx);

            if (n > 0) vx_est = (xc - xc_prev) / (rec_every * dt);
            xc_prev = xc;

            fprintf(fts, "%.4f\t%.8f\t%.8e\t%.8e\t%.6f\t%.6f\t%.6f\t%.6e\t%s\n",
                    t, xc, vx_est, Pfield, Etot, sp, sm, asym, phase);

            if (n % print_every == 0)
                printf("  t=%7.1f  xc=%+8.4f  vx=%+.3e  P=%+.3e  E=%.4f  "
                       "s+/s-=%.3f/%.3f  asym=%+.2e  [%s]\n",
                       t, xc, vx_est, Pfield, Etot, sp, sm, asym, phase);
        }

        /* Record rest-frame energy just before turning on force */
        if (n == N_equil - 1) {
            double Etot;
            find_center(phi, vel, dx, m2, Nx, xmin, &Etot);
            E_rest = Etot;
            for (int i = 1; i < Nx - 1; i++)
                for (int a = 0; a < 3; a++)
                    Q_rest += phi[a][i] * phi[a][i] * dx;
        }

        /* --- Rest-frame profile snapshot (average over window) --- */
        if (n >= snap_start && n < snap_end && (n - snap_start) % snap_every == 0) {
            double Etot;
            double xc = find_center(phi, vel, dx, m2, Nx, xmin, &Etot);
            for (int b = 0; b < prof_N; b++) {
                double xi = -prof_R + (b + 0.5) * prof_dx;
                double xabs = xc + xi;
                double fi = (xabs - xmin) / dx;
                int il = (int)fi;
                if (il < 1 || il >= Nx - 2) continue;
                double w = fi - il;
                double e = (1.0 - w) * energy_density(phi, vel, dx, m2, il)
                         +        w  * energy_density(phi, vel, dx, m2, il + 1);
                rho_rest[b] += e;
            }
            rest_snap_count++;
        }

        /* --- Accelerated profile snapshot --- */
        if (force_on && n >= accel_snap_start) {
            if ((n - accel_snap_start) % accel_snap_every == 0) {
                double Etot;
                double xc = find_center(phi, vel, dx, m2, Nx, xmin, &Etot);
                for (int b = 0; b < prof_N; b++) {
                    double xi = -prof_R + (b + 0.5) * prof_dx;
                    double xabs = xc + xi;
                    double fi = (xabs - xmin) / dx;
                    int il = (int)fi;
                    if (il < 1 || il >= Nx - 2) continue;
                    double w = fi - il;
                    double e = (1.0 - w) * energy_density(phi, vel, dx, m2, il)
                             +        w  * energy_density(phi, vel, dx, m2, il + 1);
                    rho_accel[b] += e;
                }
                accel_snap_count++;
            }
        }

        if (n == Nt) break;

        /* ==================== TIME STEP ==================== */

        if (method == 1) {
            /* Method 1: Velocity Verlet + momentum injection
             *
             * Momentum injection: at each dt, add dP = F*dt to total field momentum.
             * Field momentum: P = -Σ_a ∫ (∂_t φ_a)(∂_x φ_a) dx
             * To add dP, boost velocities: δ(vel[a][i]) = -dP * (∂φ_a/∂x) / Σ_a∫(∂φ_a/∂x)² dx
             * This is the minimal (least-squares) perturbation that adds exactly dP.
             */
            /* Half-kick */
            for (int a = 0; a < 3; a++)
                for (int i = 1; i < Nx - 1; i++)
                    vel[a][i] += 0.5 * dt * acc[a][i];

            /* Momentum injection (after half-kick, before drift) */
            if (force_on) {
                /* Compute gradient norm */
                double grad_norm2 = 0;
                for (int a = 0; a < 3; a++)
                    for (int i = 1; i < Nx - 1; i++) {
                        double g = (phi[a][i+1] - phi[a][i-1]) / (2.0 * dx);
                        grad_norm2 += g * g * dx;
                    }
                if (grad_norm2 > 1e-30) {
                    double dP = F_ext * dt;
                    double coeff = -dP / grad_norm2;
                    for (int a = 0; a < 3; a++)
                        for (int i = 1; i < Nx - 1; i++) {
                            double g = (phi[a][i+1] - phi[a][i-1]) / (2.0 * dx);
                            vel[a][i] += coeff * g;
                        }
                }
            }

            /* Drift */
            for (int a = 0; a < 3; a++)
                for (int i = 1; i < Nx - 1; i++)
                    phi[a][i] += dt * vel[a][i];

            COMPUTE_ACC();

            /* Half-kick */
            for (int a = 0; a < 3; a++)
                for (int i = 1; i < Nx - 1; i++)
                    vel[a][i] += 0.5 * dt * acc[a][i];

        } else {
            /* Method 2: Velocity Verlet with linear potential V = -(F/2)x Σφ²
             * Force: acc[a] += F*x*phi[a]
             */
            for (int a = 0; a < 3; a++)
                for (int i = 1; i < Nx - 1; i++)
                    vel[a][i] += 0.5 * dt * acc[a][i];
            for (int a = 0; a < 3; a++)
                for (int i = 1; i < Nx - 1; i++)
                    phi[a][i] += dt * vel[a][i];
            COMPUTE_ACC();
            if (force_on) ADD_LINEAR_FORCE();
            for (int a = 0; a < 3; a++)
                for (int i = 1; i < Nx - 1; i++)
                    vel[a][i] += 0.5 * dt * acc[a][i];
        }

        /* Absorbing boundary */
        for (int a = 0; a < 3; a++)
            for (int i = 0; i < Nx; i++) {
                vel[a][i] *= damp[i];
                phi[a][i] *= damp[i];
            }
    }

    fclose(fts);
    printf("\nTime series: %s\n", tspath);

    /* --- Normalize and write profiles --- */
    if (rest_snap_count > 0)
        for (int b = 0; b < prof_N; b++) rho_rest[b] /= rest_snap_count;
    if (accel_snap_count > 0)
        for (int b = 0; b < prof_N; b++) rho_accel[b] /= accel_snap_count;

    char profpath[600];
    snprintf(profpath, sizeof(profpath), "%s/inertia_F%.6f_profile.tsv", outdir, F_ext);
    FILE *fp = fopen(profpath, "w");
    fprintf(fp, "xi\trho_rest\trho_accel\tdelta_rho\n");
    for (int b = 0; b < prof_N; b++) {
        double xi = -prof_R + (b + 0.5) * prof_dx;
        fprintf(fp, "%.6f\t%.8e\t%.8e\t%.8e\n",
                xi, rho_rest[b], rho_accel[b], rho_accel[b] - rho_rest[b]);
    }
    fclose(fp);
    printf("Profiles: %s  (rest_snaps=%d, accel_snaps=%d)\n",
           profpath, rest_snap_count, accel_snap_count);

    /* --- Compute summary statistics --- */
    double odd_L1 = 0, even_L1 = 0, rho_rest_L1 = 0;
    double odd_peak = 0, even_peak = 0;
    for (int b = 0; b < prof_N / 2; b++) {
        int b2 = prof_N - 1 - b;
        double dr_p = rho_accel[b2] - rho_rest[b2];  /* leading side (+xi) */
        double dr_m = rho_accel[b]  - rho_rest[b];    /* trailing side (-xi) */
        double odd  = 0.5 * (dr_p - dr_m);   /* antisymmetric deformation */
        double even = 0.5 * (dr_p + dr_m);   /* symmetric deformation */
        odd_L1  += fabs(odd)  * prof_dx;
        even_L1 += fabs(even) * prof_dx;
        if (fabs(odd) > fabs(odd_peak)) odd_peak = odd;
        if (fabs(even) > fabs(even_peak)) even_peak = even;
        rho_rest_L1 += (fabs(rho_rest[b]) + fabs(rho_rest[b2])) * 0.5 * prof_dx;
    }

    /* Final state */
    double Etot_final;
    double xc_final = find_center(phi, vel, dx, m2, Nx, xmin, &Etot_final);
    double Pf = field_momentum(phi, vel, dx, Nx);
    double v_final = (Etot_final > 1e-20) ? Pf / Etot_final : 0;

    printf("\n=== DEFORMATION SUMMARY (F=%.6f, method=%d) ===\n", F_ext, method);
    printf("  Rest energy E_rest = %.6f\n", E_rest);
    printf("  Rest 'charge' Q = ∫Σφ² = %.6f\n", Q_rest);
    printf("  Final xc = %.6f\n", xc_final);
    printf("  Final E  = %.6f\n", Etot_final);
    printf("  Final P  = %.6e\n", Pf);
    printf("  Final v  = P/E = %.6e\n", v_final);
    printf("  Expected a = F/E_rest = %.6e\n", F_ext / (E_rest + 1e-30));
    printf("  Expected v(t_accel) = a*t = %.6e\n", F_ext * t_accel / (E_rest + 1e-30));
    printf("  Expected x(t_accel) = ½at² = %.6f\n",
           0.5 * F_ext * t_accel * t_accel / (E_rest + 1e-30));
    printf("  ---\n");
    printf("  Odd  (asymmetric) |δρ|₁ = %.6e  peak = %.6e\n", odd_L1, odd_peak);
    printf("  Even (symmetric)  |δρ|₁ = %.6e  peak = %.6e\n", even_L1, even_peak);
    printf("  Rest profile      |ρ|₁  = %.6e\n", rho_rest_L1);
    printf("  Relative odd  deformation: %.6e\n", odd_L1 / (rho_rest_L1 + 1e-30));
    printf("  Relative even deformation: %.6e\n", even_L1 / (rho_rest_L1 + 1e-30));
    printf("  Odd/Even ratio: %.4f\n", odd_L1 / (even_L1 + 1e-30));

    for (int a = 0; a < 3; a++) { free(phi[a]); free(vel[a]); free(acc[a]); }
    free(damp); free(rho_rest); free(rho_accel);
    return 0;
}
