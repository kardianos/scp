/*
 * radial1d.c — 1D radial Skyrme PDE solver for frequency calibration
 * Phase 0 of the v19 three-mode resonance experiment.
 *
 * Evolves a non-topological B=0 bump  f(r,0) = A*exp(-r^2/(2*sigma^2))
 * under the hedgehog-ansatz S^3 sigma-model wave equation:
 *
 *   f_tt = f'' + (2/r)f' - sin(2f)/(2r^2)
 *        + c4 * [...Skyrme terms...]
 *
 * The Skyrme (L_4) terms contain products of 1/r^2 factors that create
 * O(1/r^4) stiffness near the origin for non-topological configurations
 * with f(0) != n*pi, making explicit time integration unstable. Two modes:
 *
 *   Mode 1 (default): Pure sigma model (c4_eff=0). Captures the natural
 *   frequencies of the S^3 target space, which is what Phase 0 needs.
 *
 *   Mode 2 (-skyrme): Include Skyrme term only for r > r_skyrme (default 1.0)
 *   where the angular terms are mild. This gives an approximate Skyrme
 *   correction without the origin instability.
 *
 * Uses a cell-centered (staggered) grid: r_i = (i+0.5)*dr, avoiding r=0.
 * Boundary conditions: even symmetry at r=0 (ghost f[-1]=f[0]),
 *                      f -> 0 at r=rmax.
 * Time integration: Velocity Verlet (symplectic).
 * Absorbing boundary in outer 15% of grid.
 *
 * Compile: gcc -O3 -o radial1d v19/src/radial1d.c -lm
 * Usage:   ./radial1d [-A 0.8] [-sigma 1.5] [-c4 2.0] [-skyrme]
 *                     [-Nr 2000] [-rmax 30] [-tfinal 500] [-o v19/data]
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* ------------------------------------------------------------------ */
/*  Compute acceleration f_tt at each grid point.                     */
/*  Staggered grid: r_i = (i+0.5)*dr.                                */
/*  Returns max|acc| for diagnostics.                                 */
/* ------------------------------------------------------------------ */
static double compute_accel(const double *f, int Nr, double dr, double c4,
                            double r_skyrme, double *acc)
{
    double amax = 0.0;

    for (int i = 0; i < Nr; i++) {
        double r  = (i + 0.5) * dr;
        double r2 = r * r;

        double fi = f[i];

        /* Even-symmetry ghost at i=0: f[-1] = f[0] */
        double fm    = (i == 0)      ? f[0] : f[i-1];
        double fp_val = (i == Nr - 1) ? 0.0  : f[i+1];

        double fp  = (fp_val - fm) / (2.0 * dr);
        double fpp = (fp_val - 2.0 * fi + fm) / (dr * dr);

        double s2f = sin(2.0 * fi);

        /* L_2: f'' + 2f'/r - sin(2f)/(2r^2) */
        double a = fpp + 2.0 * fp / r - s2f / (2.0 * r2);

        /* L_4: Skyrme term, only for r > r_skyrme */
        if (c4 > 0.0 && r > r_skyrme) {
            double sf    = sin(fi);
            double sin2f = sf * sf;
            double A_term = fpp - s2f / (2.0 * r2);
            double L4 = c4 * (
                2.0 * fp * sin2f / r2 * A_term
              + fp * fp * s2f / r2
              - 2.0 * fp * sin2f / (r2 * r)
            );
            a += L4;
        }

        acc[i] = a;
        if (fabs(a) > amax) amax = fabs(a);
    }

    return amax;
}

/* ------------------------------------------------------------------ */
/*  Energy integrals.                                                 */
/* ------------------------------------------------------------------ */
static void compute_energy(const double *f, const double *v, int Nr,
                           double dr, double c4,
                           double *Ekin, double *E2, double *E4)
{
    double ek = 0.0, e2 = 0.0, e4 = 0.0;

    for (int i = 0; i < Nr; i++) {
        double r  = (i + 0.5) * dr;
        double r2 = r * r;
        double fi = f[i];

        double fm    = (i == 0)      ? f[0] : f[i-1];
        double fp_val = (i == Nr - 1) ? 0.0  : f[i+1];
        double fp = (fp_val - fm) / (2.0 * dr);

        double sf    = sin(fi);
        double sin2f = sf * sf;
        double sin4f = sin2f * sin2f;

        ek += v[i] * v[i] * r2;
        e2 += fp * fp * r2 + 2.0 * sin2f;
        e4 += c4 * (2.0 * fp * fp * sin2f + sin4f / r2);
    }

    *Ekin = 2.0 * M_PI * ek * dr;
    *E2   = 2.0 * M_PI * e2 * dr;
    *E4   = 2.0 * M_PI * e4 * dr;
}

/* ------------------------------------------------------------------ */
/*  Absorbing boundary: smoothstep damping in outer 15% of grid.      */
/* ------------------------------------------------------------------ */
static void apply_damping(double *v, int Nr, double dr, double rmax)
{
    double r_damp = 0.85 * rmax;
    double width  = rmax - r_damp;

    for (int i = 0; i < Nr; i++) {
        double r = (i + 0.5) * dr;
        if (r > r_damp) {
            double xi = (r - r_damp) / width;
            double alpha = 1.0 - 0.95 * xi * xi * (3.0 - 2.0 * xi);
            v[i] *= alpha;
        }
    }
}

/* ------------------------------------------------------------------ */
/*  Simple DFT power spectrum of a real signal.                       */
/* ------------------------------------------------------------------ */
static void compute_spectrum(const double *signal, int Nsig, double dt_sig,
                             double **omega_out, double **power_out,
                             int *Nfreq_out)
{
    int Nfreq = Nsig / 2;
    double *omega = (double *)malloc(Nfreq * sizeof(double));
    double *power = (double *)malloc(Nfreq * sizeof(double));
    double T = Nsig * dt_sig;

    for (int k = 1; k <= Nfreq; k++) {
        double w = 2.0 * M_PI * k / T;
        double re = 0.0, im = 0.0;
        for (int n = 0; n < Nsig; n++) {
            double phase = -w * n * dt_sig;
            re += signal[n] * cos(phase);
            im += signal[n] * sin(phase);
        }
        omega[k-1] = w;
        power[k-1] = (re * re + im * im) / ((double)Nsig * Nsig);
    }

    *omega_out = omega;
    *power_out = power;
    *Nfreq_out = Nfreq;
}

/* ------------------------------------------------------------------ */
/*  Main                                                              */
/* ------------------------------------------------------------------ */
int main(int argc, char **argv)
{
    /* Default parameters */
    double A         = 0.8;
    double sigma     = 1.5;
    double c4        = 2.0;
    int    Nr        = 2000;
    double rmax      = 30.0;
    double tfinal    = 500.0;
    int    use_skyrme = 0;       /* default: pure sigma model */
    double r_skyrme  = 1.0;     /* Skyrme cutoff radius */
    char   outdir[512] = "v19/data";

    /* Parse command line */
    for (int i = 1; i < argc; i++) {
        if      (!strcmp(argv[i], "-A")       && i+1 < argc) A         = atof(argv[++i]);
        else if (!strcmp(argv[i], "-sigma")   && i+1 < argc) sigma     = atof(argv[++i]);
        else if (!strcmp(argv[i], "-c4")      && i+1 < argc) c4        = atof(argv[++i]);
        else if (!strcmp(argv[i], "-Nr")      && i+1 < argc) Nr        = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-rmax")    && i+1 < argc) rmax      = atof(argv[++i]);
        else if (!strcmp(argv[i], "-tfinal")  && i+1 < argc) tfinal    = atof(argv[++i]);
        else if (!strcmp(argv[i], "-skyrme"))                 use_skyrme = 1;
        else if (!strcmp(argv[i], "-rskyrme") && i+1 < argc) { r_skyrme = atof(argv[++i]); use_skyrme = 1; }
        else if (!strcmp(argv[i], "-o")       && i+1 < argc) strncpy(outdir, argv[++i], 511);
        else {
            fprintf(stderr, "Usage: %s [-A 0.8] [-sigma 1.5] [-c4 2.0] [-skyrme] [-rskyrme 1.0]\n"
                    "       [-Nr 2000] [-rmax 30] [-tfinal 500] [-o v19/data]\n", argv[0]);
            return 1;
        }
    }

    double c4_eff = use_skyrme ? c4 : 0.0;
    double dr = rmax / Nr;            /* staggered: r_i = (i+0.5)*dr */
    double dt = 0.3 * dr;             /* CFL = 0.3, stable for sigma model */
    int    Nsteps = (int)(tfinal / dt + 0.5);
    int    rec_every = 10;
    int    Nrec = Nsteps / rec_every + 2;

    printf("=== radial1d: Phase 0 frequency calibration ===\n");
    printf("A=%.3f  sigma=%.3f  c4=%.3f  mode=%s\n",
           A, sigma, c4, use_skyrme ? "sigma+Skyrme" : "pure sigma model");
    if (use_skyrme)
        printf("Skyrme active for r > %.3f\n", r_skyrme);
    printf("Nr=%d  rmax=%.1f  dr=%.6f  dt=%.6f  (CFL=%.3f)\n", Nr, rmax, dr, dt, dt/dr);
    printf("tfinal=%.1f  Nsteps=%d  rec_every=%d\n", tfinal, Nsteps, rec_every);
    printf("Grid: staggered, r_0=%.6f  r_max=%.6f\n", 0.5*dr, (Nr-0.5)*dr);
    printf("Output: %s/\n\n", outdir);

    /* Allocate */
    double *f   = (double *)calloc(Nr, sizeof(double));
    double *v   = (double *)calloc(Nr, sizeof(double));
    double *acc = (double *)calloc(Nr, sizeof(double));

    double *ts_time  = (double *)malloc(Nrec * sizeof(double));
    double *ts_f0    = (double *)malloc(Nrec * sizeof(double));
    double *ts_ekin  = (double *)malloc(Nrec * sizeof(double));
    double *ts_e2    = (double *)malloc(Nrec * sizeof(double));
    double *ts_e4    = (double *)malloc(Nrec * sizeof(double));
    double *ts_etot  = (double *)malloc(Nrec * sizeof(double));

    if (!f || !v || !acc || !ts_time || !ts_f0 || !ts_ekin || !ts_e2 || !ts_e4 || !ts_etot) {
        fprintf(stderr, "Memory allocation failed\n");
        return 1;
    }

    /* Initialize: Gaussian bump */
    for (int i = 0; i < Nr; i++) {
        double r = (i + 0.5) * dr;
        f[i] = A * exp(-r * r / (2.0 * sigma * sigma));
    }

    /* Initial acceleration */
    compute_accel(f, Nr, dr, c4_eff, r_skyrme, acc);

    /* ---- Velocity Verlet time loop ---- */
    int irec = 0;
    int print_interval = Nsteps / 20 + 1;

    for (int step = 0; step <= Nsteps; step++) {
        double t = step * dt;

        /* Record */
        if (step % rec_every == 0 && irec < Nrec) {
            double Ekin, E2, E4;
            compute_energy(f, v, Nr, dr, c4, &Ekin, &E2, &E4);

            ts_time[irec] = t;
            ts_f0[irec]   = f[0];  /* f at r = dr/2, closest to origin */
            ts_ekin[irec]  = Ekin;
            ts_e2[irec]    = E2;
            ts_e4[irec]    = E4;
            ts_etot[irec]  = Ekin + E2 + E4;

            if (step % print_interval == 0) {
                printf("t=%8.2f  f(r0)=%+.6f  Ekin=%8.4f  E2=%8.4f  E4=%8.4f  Etot=%8.4f\n",
                       t, f[0], Ekin, E2, E4, Ekin + E2 + E4);
            }

            /* NaN check */
            if (isnan(f[0]) || isnan(Ekin)) {
                fprintf(stderr, "ERROR: NaN at step %d (t=%.4f). Aborting.\n", step, t);
                irec++;
                break;
            }

            irec++;
        }

        if (step == Nsteps) break;

        /* Velocity Verlet: half-kick */
        for (int i = 0; i < Nr; i++)
            v[i] += 0.5 * dt * acc[i];

        /* Drift */
        for (int i = 0; i < Nr; i++)
            f[i] += dt * v[i];

        /* New acceleration */
        compute_accel(f, Nr, dr, c4_eff, r_skyrme, acc);

        /* Second half-kick */
        for (int i = 0; i < Nr; i++)
            v[i] += 0.5 * dt * acc[i];

        /* Absorbing boundary */
        apply_damping(v, Nr, dr, rmax);
    }
    int Nrec_actual = irec;

    printf("\nSimulation complete. %d time records.\n", Nrec_actual);

    /* ---- Write timeseries ---- */
    {
        char fname[600];
        snprintf(fname, sizeof(fname), "%s/phase0_timeseries.tsv", outdir);
        FILE *fp = fopen(fname, "w");
        if (!fp) {
            fprintf(stderr, "Cannot open %s for writing\n", fname);
            return 1;
        }
        fprintf(fp, "time\tf_at_origin\tE_kin\tE_2\tE_4\tE_total\n");
        for (int i = 0; i < Nrec_actual; i++) {
            fprintf(fp, "%.6f\t%.10e\t%.10e\t%.10e\t%.10e\t%.10e\n",
                    ts_time[i], ts_f0[i], ts_ekin[i], ts_e2[i], ts_e4[i], ts_etot[i]);
        }
        fclose(fp);
        printf("Wrote %s\n", fname);
    }

    /* ---- Compute and write power spectrum ---- */
    if (Nrec_actual > 20) {
        int skip = Nrec_actual / 10;
        int Nsig = Nrec_actual - skip;
        double dt_sig = ts_time[1] - ts_time[0];

        /* Hann window, mean-subtracted */
        double *windowed = (double *)malloc(Nsig * sizeof(double));
        double mean = 0.0;
        for (int i = 0; i < Nsig; i++) mean += ts_f0[skip + i];
        mean /= Nsig;
        for (int i = 0; i < Nsig; i++) {
            double w = 0.5 * (1.0 - cos(2.0 * M_PI * i / (Nsig - 1)));
            windowed[i] = (ts_f0[skip + i] - mean) * w;
        }

        double *omega, *power;
        int Nfreq;
        compute_spectrum(windowed, Nsig, dt_sig, &omega, &power, &Nfreq);
        free(windowed);

        /* Output omega < 20 */
        int Nout = 0;
        for (int i = 0; i < Nfreq; i++) {
            if (omega[i] > 20.0) break;
            Nout++;
        }

        char fname[600];
        snprintf(fname, sizeof(fname), "%s/phase0_spectrum.tsv", outdir);
        FILE *fp = fopen(fname, "w");
        if (!fp) {
            fprintf(stderr, "Cannot open %s for writing\n", fname);
            return 1;
        }
        fprintf(fp, "omega\tpower\n");
        for (int i = 0; i < Nout; i++)
            fprintf(fp, "%.8f\t%.10e\n", omega[i], power[i]);
        fclose(fp);
        printf("Wrote %s (%d frequency bins up to omega=20)\n", fname, Nout);

        /* Find top 3 peaks */
        typedef struct { double omega; double power; } peak_t;
        peak_t peaks[3] = {{0,0},{0,0},{0,0}};
        for (int i = 1; i < Nout - 1; i++) {
            if (power[i] > power[i-1] && power[i] > power[i+1]) {
                /* Local maximum */
                for (int p = 0; p < 3; p++) {
                    if (power[i] > peaks[p].power) {
                        /* Shift down */
                        for (int q = 2; q > p; q--) peaks[q] = peaks[q-1];
                        peaks[p].omega = omega[i];
                        peaks[p].power = power[i];
                        break;
                    }
                }
            }
        }

        printf("\n=== RESULTS ===\n");
        printf("Top spectral peaks:\n");
        for (int p = 0; p < 3; p++) {
            if (peaks[p].power > 0)
                printf("  Peak %d: omega = %.4f  power = %.6e\n", p+1, peaks[p].omega, peaks[p].power);
        }

        double peak_omega = peaks[0].omega;

        /* Check oscillon: E_kin contrast in late half */
        double ekin_max = 0.0, ekin_min = 1e30;
        for (int i = Nrec_actual / 2; i < Nrec_actual; i++) {
            if (ts_ekin[i] > ekin_max) ekin_max = ts_ekin[i];
            if (ts_ekin[i] < ekin_min) ekin_min = ts_ekin[i];
        }
        double contrast = (ekin_max + ekin_min > 1e-20) ?
            (ekin_max - ekin_min) / (ekin_max + ekin_min) : 0.0;
        int oscillon = (contrast > 0.1);

        printf("\nE_kin contrast (late half): %.4f  -> %s\n",
               contrast, oscillon ? "OSCILLON DETECTED" : "decayed / no oscillon");

        printf("\nDominant frequency: omega* = %.4f\n", peak_omega);
        printf("Suggested triad frequencies:\n");
        printf("  omega_1 = omega*/3   = %.4f\n", peak_omega / 3.0);
        printf("  omega_2 = 2*omega*/3 = %.4f\n", 2.0 * peak_omega / 3.0);
        printf("  omega_3 = omega*     = %.4f\n", peak_omega);

        free(omega);
        free(power);
    } else {
        printf("\nToo few records for spectrum analysis.\n");
    }

    /* Cleanup */
    free(f); free(v); free(acc);
    free(ts_time); free(ts_f0); free(ts_ekin); free(ts_e2); free(ts_e4); free(ts_etot);

    return 0;
}
