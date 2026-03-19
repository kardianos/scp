/*  v33.c — Clean simulation: standard equation, single alloc, no modifications
 *
 *  The standard equation already produces long-range attraction (F ∝ 1/D^2.83).
 *  No gradient coupling, no c(ρ), no S/B split, no smoothing.
 *
 *  ∂²φ_a/∂t² = ∇²φ_a - m²φ_a - ∂V/∂φ_a
 *  V(P) = (μ/2)P²/(1+κP²), P = φ₀φ₁φ₂
 *
 *  ONE malloc for entire simulation. Periodic BC. Symplectic Verlet.
 *
 *  Build: gcc -O3 -march=native -fopenmp -o v33 src/v33.c -lm
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <omp.h>
#include <sys/stat.h>

#define NFIELDS 3
#define PI 3.14159265358979323846

/* Physics parameters */
static double MU    = -41.345;
static double KAPPA = 50.0;
static double MASS2 = 2.25;
static double A_BG  = 0.1;

/* ================================================================
   Grid: ONE allocation for everything
   ================================================================ */

typedef struct {
    double *mem;           /* single allocation */
    double *phi[NFIELDS];  /* pointers into mem */
    double *vel[NFIELDS];
    double *acc[NFIELDS];
    int N;
    long N3;
    double L, dx, dt;
} Grid;

static Grid *grid_alloc(int N, double L) {
    Grid *g = calloc(1, sizeof(Grid));
    g->N  = N;
    g->N3 = (long)N * N * N;
    g->L  = L;
    g->dx = 2.0 * L / (N - 1);
    g->dt = 0.12 * g->dx;

    /* ONE allocation */
    long total = 9 * g->N3;
    double bytes = total * sizeof(double);
    printf("Allocating %.2f GB (%ld doubles, N=%d)\n", bytes / 1e9, total, N);
    g->mem = malloc(total * sizeof(double));
    if (!g->mem) {
        fprintf(stderr, "FATAL: malloc failed for %.2f GB\n", bytes / 1e9);
        exit(1);
    }
    memset(g->mem, 0, total * sizeof(double));

    /* Pointers into the single block */
    for (int a = 0; a < NFIELDS; a++) {
        g->phi[a] = g->mem + (0 + a) * g->N3;
        g->vel[a] = g->mem + (3 + a) * g->N3;
        g->acc[a] = g->mem + (6 + a) * g->N3;
    }

    return g;
}

static void grid_free(Grid *g) {
    free(g->mem);  /* ONE free */
    free(g);
}

/* ================================================================
   Forces: standard equation, fully periodic
   ================================================================ */

static void compute_forces(Grid *g) {
    const int N = g->N, NN = N * N;
    const long N3 = g->N3;
    const double idx2 = 1.0 / (g->dx * g->dx);

    #pragma omp parallel for schedule(static)
    for (long idx = 0; idx < N3; idx++) {
        int i = (int)(idx / NN);
        int j = (int)((idx / N) % N);
        int k = (int)(idx % N);

        /* Fully periodic neighbors */
        int ip = (i+1)%N, im = (i-1+N)%N;
        int jp = (j+1)%N, jm = (j-1+N)%N;
        int kp = (k+1)%N, km = (k-1+N)%N;

        long n_ip = (long)ip*NN + j*N + k;
        long n_im = (long)im*NN + j*N + k;
        long n_jp = (long)i*NN + jp*N + k;
        long n_jm = (long)i*NN + jm*N + k;
        long n_kp = (long)i*NN + j*N + kp;
        long n_km = (long)i*NN + j*N + km;

        double p0 = g->phi[0][idx];
        double p1 = g->phi[1][idx];
        double p2 = g->phi[2][idx];
        double P  = p0 * p1 * p2;
        double den = 1.0 + KAPPA * P * P;
        double mPd2 = MU * P / (den * den);

        for (int a = 0; a < NFIELDS; a++) {
            double lap = (g->phi[a][n_ip] + g->phi[a][n_im]
                        + g->phi[a][n_jp] + g->phi[a][n_jm]
                        + g->phi[a][n_kp] + g->phi[a][n_km]
                        - 6.0 * g->phi[a][idx]) * idx2;

            double dPda = (a==0) ? p1*p2 : (a==1) ? p0*p2 : p0*p1;

            g->acc[a][idx] = lap - MASS2 * g->phi[a][idx] - mPd2 * dPda;
        }
    }
}

/* ================================================================
   Symplectic Velocity Verlet
   ================================================================ */

static void verlet_step(Grid *g) {
    const long N3 = g->N3;
    const double hdt = 0.5 * g->dt;
    const double dt  = g->dt;

    /* Half-kick */
    for (int a = 0; a < NFIELDS; a++) {
        double *v = g->vel[a], *ac = g->acc[a];
        for (long idx = 0; idx < N3; idx++)
            v[idx] += hdt * ac[idx];
    }
    /* Drift */
    for (int a = 0; a < NFIELDS; a++) {
        double *p = g->phi[a], *v = g->vel[a];
        for (long idx = 0; idx < N3; idx++)
            p[idx] += dt * v[idx];
    }
    /* Recompute forces */
    compute_forces(g);
    /* Half-kick */
    for (int a = 0; a < NFIELDS; a++) {
        double *v = g->vel[a], *ac = g->acc[a];
        for (long idx = 0; idx < N3; idx++)
            v[idx] += hdt * ac[idx];
    }
}

/* ================================================================
   Initialization: bimodal braid + background
   ================================================================ */

static void init_braid(Grid *g, double x_cen, double y_cen) {
    const int N = g->N, NN = N * N;
    const double dx = g->dx, L = g->L;
    const double A[3] = {0.8, 0.8, 0.8};
    const double delta[3] = {0, 3.0005, 4.4325};
    const double R_tube = 3.0, ellip = 0.3325;
    const double kw = PI / L;
    const double omega = sqrt(kw*kw + MASS2);
    const double sx = 1+ellip, sy = 1-ellip;
    const double inv2R2 = 1.0 / (2*R_tube*R_tube);
    const double k_bg = PI / L;
    const double omega_bg = sqrt(k_bg*k_bg + MASS2);

    for (int i = 0; i < N; i++) {
        double x = -L + i*dx;
        for (int j = 0; j < N; j++) {
            double y = -L + j*dx;
            for (int kk = 0; kk < N; kk++) {
                double z = -L + kk*dx;
                long idx = (long)i*NN + j*N + kk;

                double xc = x - x_cen, yc = y - y_cen;
                double r2e = xc*xc/(sx*sx) + yc*yc/(sy*sy);
                double env = exp(-r2e * inv2R2);

                for (int a = 0; a < NFIELDS; a++) {
                    double ph = kw*z + delta[a];
                    double ph_bg = k_bg*z + 2*PI*a/3.0;
                    g->phi[a][idx] += A[a]*env*cos(ph) + A_BG*cos(ph_bg);
                    g->vel[a][idx] += omega*A[a]*env*sin(ph)
                                    + omega_bg*A_BG*sin(ph_bg);
                }
            }
        }
    }
}

/* ================================================================
   Diagnostics
   ================================================================ */

/* Energy components */
static void compute_energy(Grid *g, double *E_kin, double *E_grad,
                           double *E_mass, double *E_pot, double *E_total) {
    const int N = g->N, NN = N*N;
    const long N3 = g->N3;
    const double dx = g->dx, dV = dx*dx*dx;
    double ek=0, eg=0, em=0, ep=0;

    #pragma omp parallel for reduction(+:ek,eg,em,ep) schedule(static)
    for (long idx = 0; idx < N3; idx++) {
        int i=(int)(idx/NN), j=(int)((idx/N)%N), k=(int)(idx%N);
        int ip=(i+1)%N,im=(i-1+N)%N,jp=(j+1)%N,jm=(j-1+N)%N;
        int kp=(k+1)%N,km=(k-1+N)%N;

        for (int a = 0; a < NFIELDS; a++) {
            ek += 0.5 * g->vel[a][idx] * g->vel[a][idx] * dV;
            double gx = (g->phi[a][(long)ip*NN+j*N+k]-g->phi[a][(long)im*NN+j*N+k])/(2*dx);
            double gy = (g->phi[a][(long)i*NN+jp*N+k]-g->phi[a][(long)i*NN+jm*N+k])/(2*dx);
            double gz = (g->phi[a][(long)i*NN+j*N+kp]-g->phi[a][(long)i*NN+j*N+km])/(2*dx);
            eg += 0.5*(gx*gx+gy*gy+gz*gz) * dV;
            em += 0.5*MASS2*g->phi[a][idx]*g->phi[a][idx] * dV;
        }
        double P = g->phi[0][idx]*g->phi[1][idx]*g->phi[2][idx];
        ep += (MU/2.0)*P*P/(1.0+KAPPA*P*P) * dV;
    }
    *E_kin=ek; *E_grad=eg; *E_mass=em; *E_pot=ep; *E_total=ek+eg+em+ep;
}

/* Radial energy density profile */
static void compute_radial_profile(Grid *g, double *rho_bins, int *counts,
                                    int nbins, double dr) {
    const int N = g->N, NN = N*N;
    const double dx = g->dx, L = g->L;
    for (int b = 0; b < nbins; b++) { rho_bins[b] = 0; counts[b] = 0; }

    for (int i = 0; i < N; i++) {
        double x = -L + i*dx;
        for (int j = 0; j < N; j++) {
            double y = -L + j*dx;
            double rp = sqrt(x*x + y*y);
            int b = (int)(rp / dr);
            if (b >= nbins) continue;
            for (int kk = 0; kk < N; kk++) {
                long idx = (long)i*NN + j*N + kk;
                int kp=(kk+1)%N, km=(kk-1+N)%N;
                double e = 0;
                for (int a = 0; a < NFIELDS; a++) {
                    e += 0.5*g->vel[a][idx]*g->vel[a][idx];
                    double gx=(g->phi[a][(long)((i+1)%N)*NN+j*N+kk]-g->phi[a][(long)((i-1+N)%N)*NN+j*N+kk])/(2*dx);
                    double gy=(g->phi[a][(long)i*NN+((j+1)%N)*N+kk]-g->phi[a][(long)i*NN+((j-1+N)%N)*N+kk])/(2*dx);
                    double gz=(g->phi[a][(long)i*NN+j*N+kp]-g->phi[a][(long)i*NN+j*N+km])/(2*dx);
                    e += 0.5*(gx*gx+gy*gy+gz*gz);
                    e += 0.5*MASS2*g->phi[a][idx]*g->phi[a][idx];
                }
                double P = g->phi[0][idx]*g->phi[1][idx]*g->phi[2][idx];
                e += (MU/2.0)*P*P/(1.0+KAPPA*P*P);
                rho_bins[b] += e;
                counts[b]++;
            }
        }
    }
    for (int b = 0; b < nbins; b++)
        if (counts[b] > 0) rho_bins[b] /= counts[b];
}

/* Energy flux through a shell at radius R */
static double compute_energy_flux(Grid *g, double R, double dR) {
    const int N = g->N, NN = N*N;
    const double dx = g->dx, L = g->L;
    double flux = 0;

    for (int i = 0; i < N; i++) {
        double x = -L + i*dx;
        for (int j = 0; j < N; j++) {
            double y = -L + j*dx;
            double rp = sqrt(x*x + y*y);
            if (rp < R - dR/2 || rp > R + dR/2) continue;
            double rx = x/(rp+1e-30), ry = y/(rp+1e-30);
            for (int kk = 0; kk < N; kk++) {
                long idx = (long)i*NN + j*N + kk;
                /* Radial energy flux: Σ_a v_a × (r̂·∇φ_a) */
                for (int a = 0; a < NFIELDS; a++) {
                    double dpx = (g->phi[a][(long)((i+1)%N)*NN+j*N+kk]-g->phi[a][(long)((i-1+N)%N)*NN+j*N+kk])/(2*dx);
                    double dpy = (g->phi[a][(long)i*NN+((j+1)%N)*N+kk]-g->phi[a][(long)i*NN+((j-1+N)%N)*N+kk])/(2*dx);
                    flux += g->vel[a][idx] * (rx*dpx + ry*dpy) * dx*dx*dx;
                }
            }
        }
    }
    return flux;
}

/* Braid separation (high-energy centroid tracking) */
static double measure_separation(Grid *g) {
    const long N3 = g->N3;
    const int NN = g->N * g->N, N = g->N;
    const double dx = g->dx, L = g->L;

    /* Compute average energy density for threshold */
    double avg_e = 0;
    for (long idx = 0; idx < N3; idx++) {
        double p2 = 0;
        for (int a = 0; a < NFIELDS; a++)
            p2 += g->phi[a][idx]*g->phi[a][idx];
        avg_e += p2;
    }
    avg_e /= N3;
    double thresh = 3.0 * avg_e;

    double wL=0, wR=0, xL=0, xR=0;
    for (long idx = 0; idx < N3; idx++) {
        double p2 = 0;
        for (int a = 0; a < NFIELDS; a++)
            p2 += g->phi[a][idx]*g->phi[a][idx];
        if (p2 < thresh) continue;
        int i = (int)(idx / NN);
        double x = -L + i*dx;
        if (x < 0) { wL += p2; xL += x*p2; }
        else       { wR += p2; xR += x*p2; }
    }
    if (wL > 0) xL /= wL;
    if (wR > 0) xR /= wR;
    return xR - xL;
}

/* Save field snapshot */
static void save_field(Grid *g, double t, const char *dir) {
    char fn[512];
    snprintf(fn, sizeof(fn), "%s/field_t%04d.bin", dir, (int)(t+0.5));
    FILE *fp = fopen(fn, "wb");
    if (!fp) return;
    int n = g->N; double l = g->L;
    fwrite(&n, sizeof(int), 1, fp);
    fwrite(&l, sizeof(double), 1, fp);
    fwrite(&t, sizeof(double), 1, fp);
    for (int a = 0; a < NFIELDS; a++)
        fwrite(g->phi[a], sizeof(double), g->N3, fp);
    fclose(fp);
    printf("  [SNAP] %s (%.1f GB)\n", fn, 3.0*g->N3*8.0/1e9);
}

/* ================================================================
   Main
   ================================================================ */

int main(int argc, char **argv) {
    int N = 128;
    double L = 30.0, T = 300.0;
    int n_braids = 1;
    double D = 20.0;
    double diag_dt = 5.0, snap_dt = 100.0, prof_dt = 50.0;
    char outdir[256] = "data/run";

    for (int i = 1; i < argc; i++) {
        if      (!strcmp(argv[i],"-N"))      N = atoi(argv[++i]);
        else if (!strcmp(argv[i],"-L"))      L = atof(argv[++i]);
        else if (!strcmp(argv[i],"-T"))      T = atof(argv[++i]);
        else if (!strcmp(argv[i],"-braids")) n_braids = atoi(argv[++i]);
        else if (!strcmp(argv[i],"-D"))      D = atof(argv[++i]);
        else if (!strcmp(argv[i],"-bg"))     A_BG = atof(argv[++i]);
        else if (!strcmp(argv[i],"-m"))      MASS2 = atof(argv[++i]) * atof(argv[i]);
        else if (!strcmp(argv[i],"-diag"))   diag_dt = atof(argv[++i]);
        else if (!strcmp(argv[i],"-snap"))   snap_dt = atof(argv[++i]);
        else if (!strcmp(argv[i],"-prof"))   prof_dt = atof(argv[++i]);
        else if (!strcmp(argv[i],"-o"))      strncpy(outdir, argv[++i], 255);
    }

    omp_set_num_threads(16);

    printf("=== V33: Clean Standard Equation ===\n");
    printf("∂²φ/∂t² = ∇²φ - m²φ - V'(φ)\n");
    printf("No modifications. Single alloc. Symplectic Verlet.\n\n");
    printf("N=%d L=%.0f T=%.0f braids=%d D=%.0f bg=%.2f m²=%.2f\n",
           N, L, T, n_braids, D, A_BG, MASS2);

    mkdir("data", 0755);
    mkdir(outdir, 0755);

    Grid *g = grid_alloc(N, L);
    printf("dx=%.4f dt=%.5f\n\n", g->dx, g->dt);

    /* Initialize */
    if (n_braids == 1) {
        init_braid(g, 0, 0);
    } else {
        for (int b = 0; b < n_braids; b++) {
            double angle = (n_braids == 2) ? PI * b : 2*PI*b/n_braids;
            double bx = (n_braids == 2) ? (b == 0 ? -D/2 : D/2) : D/2*cos(angle);
            double by = (n_braids == 2) ? 0 : D/2*sin(angle);
            printf("  Braid %d at (%.1f, %.1f)\n", b, bx, by);
            init_braid(g, bx, by);
        }
    }

    compute_forces(g);

    /* Output files */
    char tspath[512], profpath[512];
    snprintf(tspath, sizeof(tspath), "%s/timeseries.tsv", outdir);
    snprintf(profpath, sizeof(profpath), "%s/profiles.tsv", outdir);
    FILE *fp_ts = fopen(tspath, "w");
    FILE *fp_prof = fopen(profpath, "w");
    fprintf(fp_ts, "t\tE_kin\tE_grad\tE_mass\tE_pot\tE_total\tD\n");
    fprintf(fp_prof, "t\tr\trho\n");

    int n_steps  = (int)(T / g->dt);
    int diag_every = (int)(diag_dt / g->dt); if (diag_every < 1) diag_every = 1;
    int snap_every = (int)(snap_dt / g->dt); if (snap_every < 1) snap_every = 1;
    int prof_every = (int)(prof_dt / g->dt); if (prof_every < 1) prof_every = 1;

    printf("Steps=%d diag=%d snap=%d prof=%d\n\n", n_steps, diag_every, snap_every, prof_every);

    /* Radial profile bins */
    int nbins = (int)(L / 0.5) + 1;
    double dr = 0.5;
    double *rho_bins = calloc(nbins, sizeof(double));
    int *counts = calloc(nbins, sizeof(int));

    double wall0 = omp_get_wtime();
    double E0 = 0;  /* initial energy for conservation check */
    save_field(g, 0, outdir);

    for (int step = 0; step <= n_steps; step++) {
        if (step > 0) verlet_step(g);
        double t = step * g->dt;

        /* Diagnostics */
        if (step % diag_every == 0) {
            double ek, eg, em, ep, et;
            compute_energy(g, &ek, &eg, &em, &ep, &et);
            if (step == 0) E0 = et;
            double Dsep = (n_braids > 1) ? measure_separation(g) : 0;

            fprintf(fp_ts, "%.2f\t%.4e\t%.4e\t%.4e\t%.4e\t%.4e\t%.2f\n",
                    t, ek, eg, em, ep, et, Dsep);
            fflush(fp_ts);

            if (step % (diag_every * 10) == 0) {
                double wall = omp_get_wtime() - wall0;
                double pct = 100.0 * step / n_steps;
                double drift_pct = 100.0 * (et - E0) / (fabs(E0) + 1e-30);
                printf("t=%7.1f E=%.4e (drift %+.3f%%) D=%.1f [%.0f%% %.0fs]\n",
                       t, et, drift_pct, Dsep, pct, wall);
                fflush(stdout);
            }
        }

        /* Radial profile */
        if (step % prof_every == 0) {
            compute_radial_profile(g, rho_bins, counts, nbins, dr);
            for (int b = 0; b < nbins; b++) {
                if (counts[b] > 0)
                    fprintf(fp_prof, "%.2f\t%.2f\t%.6e\n", t, (b+0.5)*dr, rho_bins[b]);
            }
            fflush(fp_prof);
        }

        /* Snapshots */
        if (step > 0 && step % snap_every == 0)
            save_field(g, t, outdir);
    }

    save_field(g, T, outdir);
    fclose(fp_ts);
    fclose(fp_prof);

    double wall = omp_get_wtime() - wall0;
    printf("\n=== Complete: %.0fs (%.1f min) ===\n", wall, wall/60);
    printf("Avg: %.4f s/step\n", wall / n_steps);

    free(rho_bins);
    free(counts);
    grid_free(g);
    return 0;
}
