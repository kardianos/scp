/*
 * thermal1d.c — V24 Test 3: Thermal Bath Stabilization of Oscillon
 *
 * Tests whether stochastic thermal noise stabilizes or destabilizes the
 * three-field oscillon. Langevin dynamics with spatially-selective damping:
 *   - Noise eta_a applied EVERYWHERE
 *   - Damping gamma * dphi_a/dt applied ONLY for |x| > R_damp
 *
 * Key diagnostic: OSCILLON OVERLAP integral
 *   O(t) = integral[ phi_sym(x,t) * f_eq(x) dx ] / integral[ f_eq(x)^2 dx ]
 * where phi_sym = (phi1+phi2+phi3)/3 is the symmetric mode and f_eq(x) is
 * the equilibrated profile saved at t=t_equil. O=1 means perfect alignment
 * with equilibrium shape. O=0 means no coherent oscillon.
 *
 * The envelope of |O(t)| over breathing cycles tracks oscillon survival.
 * Noise is spatially uncorrelated, so <O_noise> = 0 on average.
 *
 * Compile: gcc -O3 -Wall -o thermal1d v24/src/thermal1d.c -lm
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Physics parameters */
static double mu      = -20.0;
static double kappa   = 20.0;
static double mass    = 1.0;
static double A_init  = 0.8;
static double sigma_w = 3.0;

/* Grid */
static int    Nx      = 4000;
static double xmax    = 80.0;

/* Time */
static double t_equil   = 5000.0;
static double t_thermal = 10000.0;

/* Thermal bath */
static double T_bath     = 0.0;
static double gamma_damp = 0.01;
static double R_damp     = 20.0;

/* Core measurement window */
static double R_core    = 15.0;

/* RNG seed */
static unsigned long long rng_state = 123456789ULL;

static char outdir[512] = "v24/data";

static void parse_args(int argc, char **argv)
{
    for (int i = 1; i < argc - 1; i += 2) {
        if      (!strcmp(argv[i], "-mu"))       mu       = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-kappa"))    kappa    = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-mass"))     mass     = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-A"))        A_init   = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-sigma"))    sigma_w  = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-Nx"))       Nx       = atoi(argv[i+1]);
        else if (!strcmp(argv[i], "-xmax"))     xmax     = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-t_equil"))  t_equil  = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-t_thermal"))t_thermal= atof(argv[i+1]);
        else if (!strcmp(argv[i], "-T"))        T_bath   = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-gamma"))    gamma_damp = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-R_damp"))   R_damp   = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-R_core"))   R_core   = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-seed"))     rng_state= (unsigned long long)atol(argv[i+1]);
        else if (!strcmp(argv[i], "-o"))        strncpy(outdir, argv[i+1], sizeof(outdir)-1);
        else { fprintf(stderr, "Unknown flag: %s\n", argv[i]); exit(1); }
    }
}

/* xorshift64 RNG */
static double rng_uniform(void)
{
    rng_state ^= rng_state << 13;
    rng_state ^= rng_state >> 7;
    rng_state ^= rng_state << 17;
    return (rng_state & 0x1FFFFFFFFFFFFFULL) / (double)0x20000000000000ULL;
}

/* Box-Muller */
static void rng_gauss(double *z1, double *z2)
{
    double u1 = rng_uniform();
    double u2 = rng_uniform();
    if (u1 < 1e-300) u1 = 1e-300;
    double r = sqrt(-2.0 * log(u1));
    *z1 = r * cos(2.0 * M_PI * u2);
    *z2 = r * sin(2.0 * M_PI * u2);
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

int main(int argc, char **argv)
{
    parse_args(argc, argv);

    double dx  = 2.0 * xmax / (Nx - 1);
    double dx2 = dx * dx;
    double m2  = mass * mass;

    /* CFL */
    double kmax = M_PI / dx;
    double dt = 0.8 * 2.0 / sqrt(kmax * kmax + m2);
    int Nt_equil   = (int)(t_equil / dt) + 1;
    int Nt_thermal = (int)(t_thermal / dt) + 1;
    int Nt_total   = Nt_equil + Nt_thermal;

    printf("thermal1d: T=%.6f gamma=%.3f R_damp=%.1f\n", T_bath, gamma_damp, R_damp);
    printf("  mu=%.1f kappa=%.1f mass=%.3f A=%.3f sigma=%.3f\n",
           mu, kappa, mass, A_init, sigma_w);
    printf("  Nx=%d xmax=%.1f dx=%.5f dt=%.6f\n", Nx, xmax, dx, dt);
    printf("  t_equil=%.0f (%d steps) t_thermal=%.0f (%d steps)\n",
           t_equil, Nt_equil, t_thermal, Nt_thermal);

    double noise_amp = (T_bath > 0) ? sqrt(2.0 * T_bath / (dx * dt)) : 0.0;
    printf("  noise_amp=%.6e\n", noise_amp);

    /* Allocate fields */
    double *phi[3], *vel[3], *acc[3];
    for (int a = 0; a < 3; a++) {
        phi[a] = calloc(Nx, sizeof(double));
        vel[a] = calloc(Nx, sizeof(double));
        acc[a] = calloc(Nx, sizeof(double));
    }

    /* Equilibrium profile template (saved at t_equil) */
    double *f_eq = calloc(Nx, sizeof(double));
    double f_eq_norm = 0.0; /* integral f_eq^2 dx */

    /* Absorbing boundary: outer 10% of domain */
    double *bdamp = malloc(Nx * sizeof(double));
    double x_abs = xmax * 0.90;
    for (int i = 0; i < Nx; i++) {
        double x = -xmax + i * dx;
        double ax = fabs(x);
        if (ax > x_abs) {
            double f = (ax - x_abs) / (xmax - x_abs);
            bdamp[i] = 1.0 - 0.98 * f * f;
        } else {
            bdamp[i] = 1.0;
        }
    }

    /* Damping mask: gamma only for |x| > R_damp, smooth transition */
    double *gdamp = calloc(Nx, sizeof(double));
    for (int i = 0; i < Nx; i++) {
        double x = -xmax + i * dx;
        double ax = fabs(x);
        if (ax > R_damp + 5.0) {
            gdamp[i] = gamma_damp;
        } else if (ax > R_damp) {
            double f = (ax - R_damp) / 5.0;
            gdamp[i] = gamma_damp * f * f * (3.0 - 2.0 * f);
        }
    }

    /* Initialize: Gaussians (symmetric triad) */
    for (int a = 0; a < 3; a++)
        for (int i = 0; i < Nx; i++) {
            double x = -xmax + i * dx;
            phi[a][i] = A_init * exp(-x * x / (2.0 * sigma_w * sigma_w));
        }

    #define COMPUTE_ACC() do { \
        for (int a = 0; a < 3; a++) { \
            acc[a][0] = acc[a][1] = acc[a][Nx-2] = acc[a][Nx-1] = 0; \
            for (int i = 1; i < Nx - 1; i++) { \
                double lapl = (phi[a][i+1] - 2.0*phi[a][i] + phi[a][i-1]) / dx2; \
                double fp = force_pot(phi[0][i], phi[1][i], phi[2][i], a); \
                acc[a][i] = lapl - m2 * phi[a][i] + fp; \
            } \
        } \
    } while(0)

    COMPUTE_ACC();

    /* Output file */
    char tspath[600];
    snprintf(tspath, sizeof(tspath), "%s/thermal_T%.6f_ts.tsv", outdir, T_bath);
    FILE *fts = fopen(tspath, "w");
    if (!fts) { fprintf(stderr, "Cannot open %s\n", tspath); return 1; }
    fprintf(fts, "time\tphase\toverlap\toverlap_env\tE_core\tpeak_amp\tphi1_0\n");

    int rec_every = Nt_total / 20000;
    if (rec_every < 1) rec_every = 1;
    int print_every = Nt_total / 40;
    if (print_every < 1) print_every = 1;

    /* Overlap envelope tracking: running max of |O| over ~1 breathing period */
    double overlap_env = 0.0;    /* current envelope value */
    double overlap_env_init = 0.0; /* envelope at start of thermal phase */
    int thermal_on = 0;
    double lifetime = t_thermal;
    int lifetime_found = 0;

    /* Ring buffer for running-max envelope (over ~1 period, T_breath ~ 7.2) */
    int env_buf_len = (int)(8.0 / (dt * rec_every)) + 1;
    if (env_buf_len < 10) env_buf_len = 10;
    if (env_buf_len > 2000) env_buf_len = 2000;
    double *env_buf = calloc(env_buf_len, sizeof(double));
    int env_idx = 0;
    int env_filled = 0;

    /* For computing running envelope after thermal start */
    int thermal_rec_count = 0;
    /* Track envelope at specific checkpoints */
    double env_at_1000 = -1, env_at_2000 = -1, env_at_5000 = -1;

    double *eta[3];
    for (int a = 0; a < 3; a++)
        eta[a] = calloc(Nx, sizeof(double));

    int ic = Nx / 2;

    for (int n = 0; n <= Nt_total; n++) {
        double t = n * dt;
        int phase_thermal = (n >= Nt_equil) ? 1 : 0;

        /* Save equilibrium profile at t_equil */
        if (n == Nt_equil) {
            /* Use the PEAK of the breathing cycle as reference.
             * Actually, save current phi_sym profile and normalize. */
            f_eq_norm = 0.0;
            for (int i = 0; i < Nx; i++) {
                f_eq[i] = (phi[0][i] + phi[1][i] + phi[2][i]) / 3.0;
                f_eq_norm += f_eq[i] * f_eq[i] * dx;
            }
            if (f_eq_norm < 1e-20) f_eq_norm = 1e-20;

            thermal_on = 1;
            printf("\n--- Thermal phase starts (t=%.1f) ---\n", t);
            printf("  f_eq_norm = %.6f  phi_center = %.6f\n\n",
                   f_eq_norm, phi[0][ic]);

            /* Calibrate envelope over next breathing cycles */
        }

        /* Record */
        int do_rec = (n % rec_every == 0);
        int do_print = (n % print_every == 0);

        if (do_rec || do_print) {
            /* Compute overlap: O = int(phi_sym * f_eq dx) / int(f_eq^2 dx) */
            double overlap = 0.0;
            double Ecore = 0.0;
            double pk = 0.0;

            for (int i = 1; i < Nx - 1; i++) {
                double x = -xmax + i * dx;
                double phi_sym = (phi[0][i] + phi[1][i] + phi[2][i]) / 3.0;
                overlap += phi_sym * f_eq[i] * dx;

                if (fabs(x) < R_core) {
                    double e = 0;
                    for (int a = 0; a < 3; a++) {
                        e += 0.5 * vel[a][i] * vel[a][i];
                        double dp = (phi[a][i+1] - phi[a][i-1]) / (2.0*dx);
                        e += 0.5*dp*dp + 0.5*m2*phi[a][i]*phi[a][i];
                        if (fabs(phi[a][i]) > pk) pk = fabs(phi[a][i]);
                    }
                    double P = phi[0][i]*phi[1][i]*phi[2][i];
                    double P2 = P*P;
                    e += 0.5*mu*P2/(1.0+kappa*P2);
                    Ecore += e * dx;
                }
            }
            overlap /= f_eq_norm;

            /* Update running-max envelope */
            if (thermal_on) {
                env_buf[env_idx] = fabs(overlap);
                env_idx = (env_idx + 1) % env_buf_len;
                if (!env_filled && env_idx == 0) env_filled = 1;

                int len = env_filled ? env_buf_len : env_idx;
                double mx = 0;
                for (int j = 0; j < len; j++)
                    if (env_buf[j] > mx) mx = env_buf[j];
                overlap_env = mx;

                thermal_rec_count++;

                /* After first full buffer, record initial envelope */
                if (thermal_rec_count == env_buf_len)
                    overlap_env_init = overlap_env;

                /* Track envelope at checkpoints */
                double t_therm = (n - Nt_equil) * dt;
                if (env_at_1000 < 0 && t_therm >= 1000) env_at_1000 = overlap_env;
                if (env_at_2000 < 0 && t_therm >= 2000) env_at_2000 = overlap_env;
                if (env_at_5000 < 0 && t_therm >= 5000) env_at_5000 = overlap_env;
            }

            if (do_rec)
                fprintf(fts, "%.4f\t%d\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\n",
                        t, phase_thermal, overlap, overlap_env, Ecore, pk, phi[0][ic]);

            if (do_print) {
                printf("  t=%8.1f [%s]  O=%+.4f  env=%.4f  Ec=%.3f  pk=%.4f  p0=%+.4f",
                       t, phase_thermal ? "THERM" : "equil",
                       overlap, overlap_env, Ecore, pk, phi[0][ic]);
                if (thermal_on && overlap_env_init > 1e-10)
                    printf("  env/env0=%.4f", overlap_env / overlap_env_init);
                printf("\n");
            }

            /* Lifetime: envelope drops below 50% of initial */
            if (thermal_on && !lifetime_found && overlap_env_init > 1e-10) {
                if (overlap_env < 0.5 * overlap_env_init) {
                    lifetime = (n - Nt_equil) * dt;
                    lifetime_found = 1;
                    printf("  *** LIFETIME at t_therm=%.1f (env/env0=%.4f) ***\n",
                           lifetime, overlap_env / overlap_env_init);
                }
            }
        }

        if (n == Nt_total) break;

        /* Generate noise */
        if (thermal_on && noise_amp > 0) {
            for (int a = 0; a < 3; a++) {
                for (int i = 0; i < Nx - 1; i += 2) {
                    double z1, z2;
                    rng_gauss(&z1, &z2);
                    eta[a][i]   = noise_amp * z1;
                    eta[a][i+1] = noise_amp * z2;
                }
                if (Nx % 2 == 1) {
                    double z1, z2;
                    rng_gauss(&z1, &z2);
                    eta[a][Nx-1] = noise_amp * z1;
                }
            }
        }

        /* Velocity Verlet */
        for (int a = 0; a < 3; a++)
            for (int i = 1; i < Nx - 1; i++)
                vel[a][i] += 0.5 * dt * acc[a][i];
        for (int a = 0; a < 3; a++)
            for (int i = 1; i < Nx - 1; i++)
                phi[a][i] += dt * vel[a][i];
        COMPUTE_ACC();
        for (int a = 0; a < 3; a++)
            for (int i = 1; i < Nx - 1; i++)
                vel[a][i] += 0.5 * dt * acc[a][i];

        /* Langevin: damping (far-field) + noise (everywhere) */
        if (thermal_on) {
            for (int a = 0; a < 3; a++) {
                for (int i = 1; i < Nx - 1; i++) {
                    vel[a][i] -= dt * gdamp[i] * vel[a][i];
                    if (noise_amp > 0)
                        vel[a][i] += dt * eta[a][i];
                }
            }
        }

        /* Absorbing boundary */
        for (int a = 0; a < 3; a++)
            for (int i = 0; i < Nx; i++) {
                vel[a][i] *= bdamp[i];
                phi[a][i] *= bdamp[i];
            }
    }

    fclose(fts);

    /* Final summary */
    double env_ratio = (overlap_env_init > 1e-10) ?
        overlap_env / overlap_env_init : 0.0;

    printf("\n=== SUMMARY ===\n");
    printf("T=%.6f  gamma=%.3f  R_damp=%.1f\n", T_bath, gamma_damp, R_damp);
    printf("Overlap envelope: init=%.6f final=%.6f ratio=%.4f\n",
           overlap_env_init, overlap_env, env_ratio);
    if (env_at_1000 >= 0 && overlap_env_init > 0)
        printf("  env(t=1k)/env0=%.4f\n", env_at_1000/overlap_env_init);
    if (env_at_2000 >= 0 && overlap_env_init > 0)
        printf("  env(t=2k)/env0=%.4f\n", env_at_2000/overlap_env_init);
    if (env_at_5000 >= 0 && overlap_env_init > 0)
        printf("  env(t=5k)/env0=%.4f\n", env_at_5000/overlap_env_init);
    printf("Lifetime: %.1f (max=%.1f)\n", lifetime, t_thermal);
    printf("Output: %s\n", tspath);

    /* Append to summary */
    char sumpath[600];
    snprintf(sumpath, sizeof(sumpath), "%s/thermal_summary.tsv", outdir);
    FILE *fsum = fopen(sumpath, "a");
    if (fsum) {
        fseek(fsum, 0, SEEK_END);
        if (ftell(fsum) == 0)
            fprintf(fsum, "T\tgamma\tR_damp\tenv_init\tenv_final\tenv_ratio\t"
                          "env_1k\tenv_2k\tenv_5k\tlifetime\tt_thermal\n");
        fprintf(fsum, "%.6f\t%.3f\t%.1f\t%.6f\t%.6f\t%.4f\t%.4f\t%.4f\t%.4f\t%.1f\t%.1f\n",
                T_bath, gamma_damp, R_damp, overlap_env_init, overlap_env, env_ratio,
                (env_at_1000 >= 0 && overlap_env_init > 0) ? env_at_1000/overlap_env_init : -1.0,
                (env_at_2000 >= 0 && overlap_env_init > 0) ? env_at_2000/overlap_env_init : -1.0,
                (env_at_5000 >= 0 && overlap_env_init > 0) ? env_at_5000/overlap_env_init : -1.0,
                lifetime, t_thermal);
        fclose(fsum);
    }

    for (int a = 0; a < 3; a++) {
        free(phi[a]); free(vel[a]); free(acc[a]); free(eta[a]);
    }
    free(bdamp); free(gdamp); free(f_eq); free(env_buf);
    return 0;
}
