/*
 * m2.c — V27-M2: Gravity Between Two Massless Propagating Braids
 *
 * Uses the optimal m=0, mu=-50, kappa=50 configuration from M4.
 *
 * Lagrangian: L = sum_a [1/2(dt phi_a)^2 - 1/2(di phi_a)^2]
 *                 - (mu/2) P^2/(1+kappa*P^2)
 * where P = phi_1 * phi_2 * phi_3.  NO mass term.
 *
 * Tests:
 *   M2a: Strain field multipole decomposition at R={5,8,12,16}
 *   M2b: Two-braid interaction at D=30 along x
 *   M2c: Angular pattern of force (2nd braid at different angles)
 *
 * Compile: gcc -O3 -fopenmp -Wall -o m2 src/m2.c -lm
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>

/* ─── Parameters ─── */
static double mu_pot    = -50.0;
static double kappa     = 50.0;
static double A0        = 0.8;
static double R_tube    = 3.0;

static int    N         = 128;
static double L         = 20.0;   /* half-box: domain is [-L, +L]^3 */
static double tfinal    = 300.0;
static double cfl_frac  = 0.20;
static char   outdir[512] = "data";

/* ─── Index helpers ─── */
#define IDX(i,j,k) ((long)(i)*N*N + (long)(j)*N + (long)(k))

/* ─── Globals ─── */
static double *phi[3], *vel[3], *acc[3];
static double *damp;
static double dx, dx2, dt;
static long Ngrid;

/* ─── Periodic z wrap ─── */
static inline int wrap_z(int k)
{
    if (k < 0)   return k + N;
    if (k >= N)  return k - N;
    return k;
}

/* ─── Grid coords ─── */
static inline double xcoord(int i) { return -L + i * dx; }
static inline double ycoord(int j) { return -L + j * dx; }
static inline double zcoord(int k) { return -L + k * dx; }

/* ─── Initialize single braid centered at (cx, cy) ─── */
static void add_braid(double cx, double cy, double sign)
{
    /* sign = +1 or -1 for same/opposite helicity */
    double Lz = 2.0 * L;
    double k_wave = 2.0 * M_PI / Lz;
    double omega = k_wave;  /* m=0: omega = k */

    #pragma omp parallel for schedule(static)
    for (long idx = 0; idx < Ngrid; idx++) {
        int ii = idx / (N * N);
        int jj = (idx / N) % N;
        int kk = idx % N;
        double x = xcoord(ii) - cx;
        double y = ycoord(jj) - cy;
        double z = zcoord(kk);

        double r_perp = sqrt(x * x + y * y);
        double envelope = A0 * exp(-r_perp * r_perp / (2.0 * R_tube * R_tube));

        for (int a = 0; a < 3; a++) {
            double phase = sign * k_wave * z + 2.0 * M_PI * a / 3.0;
            phi[a][idx] += envelope * cos(phase);
            vel[a][idx] += omega * envelope * sin(phase);
        }
    }
}

/* ─── Compute acceleration (m=0, no mass term) ─── */
static void compute_acc(void)
{
    for (int a = 0; a < 3; a++) {
        #pragma omp parallel for schedule(static)
        for (int i = 2; i < N-2; i++) {
            for (int j = 2; j < N-2; j++) {
                for (int k = 0; k < N; k++) {
                    long idx = IDX(i, j, k);

                    int km1 = wrap_z(k - 1);
                    int kp1 = wrap_z(k + 1);

                    double lapl = (phi[a][IDX(i+1,j,k)] + phi[a][IDX(i-1,j,k)]
                                 + phi[a][IDX(i,j+1,k)] + phi[a][IDX(i,j-1,k)]
                                 + phi[a][IDX(i,j,kp1)] + phi[a][IDX(i,j,km1)]
                                 - 6.0 * phi[a][idx]) / dx2;

                    /* Triple-product potential force */
                    double p0 = phi[0][idx], p1 = phi[1][idx], p2 = phi[2][idx];
                    double P = p0 * p1 * p2;
                    double P2 = P * P;
                    double denom2 = (1.0 + kappa * P2);
                    denom2 *= denom2;

                    double dP;
                    switch (a) {
                        case 0: dP = p1 * p2; break;
                        case 1: dP = p0 * p2; break;
                        default: dP = p0 * p1; break;
                    }
                    double dVdphi = mu_pot * P * dP / denom2;

                    acc[a][idx] = lapl - dVdphi;
                }
            }
        }

        /* Boundary: zero acceleration at x,y edges */
        #pragma omp parallel for schedule(static)
        for (long idx = 0; idx < Ngrid; idx++) {
            int ii = idx / (N * N);
            int jj = (idx / N) % N;
            if (ii < 2 || ii >= N-2 || jj < 2 || jj >= N-2) {
                acc[a][idx] = 0.0;
            }
        }
    }
}

/* ─── Velocity Verlet step ─── */
static void verlet_step(void)
{
    for (int a = 0; a < 3; a++) {
        #pragma omp parallel for schedule(static)
        for (long idx = 0; idx < Ngrid; idx++)
            vel[a][idx] += 0.5 * dt * acc[a][idx];
    }

    for (int a = 0; a < 3; a++) {
        #pragma omp parallel for schedule(static)
        for (long idx = 0; idx < Ngrid; idx++)
            phi[a][idx] += dt * vel[a][idx];
    }

    compute_acc();

    for (int a = 0; a < 3; a++) {
        #pragma omp parallel for schedule(static)
        for (long idx = 0; idx < Ngrid; idx++)
            vel[a][idx] += 0.5 * dt * acc[a][idx];
    }

    /* Absorbing boundary damping (x,y only) */
    for (int a = 0; a < 3; a++) {
        #pragma omp parallel for schedule(static)
        for (long idx = 0; idx < Ngrid; idx++) {
            vel[a][idx] *= damp[idx];
            phi[a][idx] *= damp[idx];
        }
    }
}

/* ─── Setup damping layer ─── */
static void setup_damping(void)
{
    double R_abs_inner = L * 0.70;
    double R_abs_outer = L * 0.95;
    #pragma omp parallel for schedule(static)
    for (long idx = 0; idx < Ngrid; idx++) {
        int ii = idx / (N * N);
        int jj = (idx / N) % N;
        double x = xcoord(ii);
        double y = ycoord(jj);
        double rperp = sqrt(x*x + y*y);
        if (rperp > R_abs_inner) {
            double f = (rperp - R_abs_inner) / (R_abs_outer - R_abs_inner);
            if (f > 1.0) f = 1.0;
            damp[idx] = 1.0 - 0.98 * f * f;
        } else {
            damp[idx] = 1.0;
        }
    }
}

/* ─── Setup damping for two-braid (rectangular, damps x and y edges) ─── */
static void setup_damping_rect(void)
{
    double margin = L * 0.15;  /* damping layer thickness */
    #pragma omp parallel for schedule(static)
    for (long idx = 0; idx < Ngrid; idx++) {
        int ii = idx / (N * N);
        int jj = (idx / N) % N;
        double x = xcoord(ii);
        double y = ycoord(jj);

        double fx = 0.0, fy = 0.0;
        /* Left/right x edges */
        if (x < -L + margin) {
            fx = (-L + margin - x) / margin;
        } else if (x > L - margin) {
            fx = (x - (L - margin)) / margin;
        }
        /* Top/bottom y edges */
        if (y < -L + margin) {
            fy = (-L + margin - y) / margin;
        } else if (y > L - margin) {
            fy = (y - (L - margin)) / margin;
        }
        double f = fx > fy ? fx : fy;
        if (f > 1.0) f = 1.0;
        damp[idx] = 1.0 - 0.98 * f * f;
    }
}

/* ──────────────────────────────────────────────────────────────────────
 *  M2a: Strain field multipole decomposition
 *  Compute strain tensor eps_{ij} = 1/2(d_i phi_j + d_j phi_i)
 *  at shells of radius R from braid center, decompose into l=0,1,2
 * ─────────────────────────────────────────────────────────────────── */

/* Spherical harmonics for l=0,2 decomposition of symmetric traceless tensor */
/* We decompose the energy density on spherical shells */
static void run_m2a(void)
{
    printf("================================================================\n");
    printf("  M2a: Strain Field Multipole Decomposition (Single Braid)\n");
    printf("================================================================\n\n");

    /* Single braid at origin */
    for (int a = 0; a < 3; a++) {
        memset(phi[a], 0, Ngrid * sizeof(double));
        memset(vel[a], 0, Ngrid * sizeof(double));
        memset(acc[a], 0, Ngrid * sizeof(double));
    }
    add_braid(0.0, 0.0, +1.0);
    setup_damping();

    /* Evolve to let the braid settle */
    printf("  Evolving single braid for t=100 to settle...\n");
    compute_acc();
    int Nsettle = (int)(100.0 / dt);
    for (int n = 0; n < Nsettle; n++) {
        verlet_step();
        if (n % (Nsettle/5) == 0) {
            printf("    t=%.1f\n", n * dt);
            fflush(stdout);
        }
    }

    /* Now compute strain at several radii */
    double R_shells[] = {5.0, 8.0, 12.0, 16.0};
    int n_shells = 4;
    double shell_width = 1.5 * dx;  /* shell thickness for sampling */

    char path[600];
    snprintf(path, sizeof(path), "%s/m2a_strain.tsv", outdir);
    FILE *fout = fopen(path, "w");
    if (!fout) { fprintf(stderr, "Cannot open %s\n", path); return; }
    fprintf(fout, "R\tl0_power\tl1_power\tl2_power\tl2_over_l0\ttotal_energy\n");

    printf("  Computing strain multipoles...\n");
    for (int iR = 0; iR < n_shells; iR++) {
        double R = R_shells[iR];

        /*
         * Strategy: sample energy density on the shell.
         * Decompose into spherical harmonics using projection.
         *
         * l=0: <e> (isotropic average)
         * l=1: <e * cos(theta)> etc  (dipole)
         * l=2: <e * (3cos^2(theta)-1)/2> etc (quadrupole)
         *
         * For symmetric strain tensor, the l=2 content comes from
         * the traceless symmetric part of eps_{ij} at each point.
         */

        /* Accumulate multipole moments on the shell */
        double Y00_sum = 0, Y00_count = 0;
        double Y10_sum = 0, Y11c_sum = 0, Y11s_sum = 0;
        double Y20_sum = 0, Y21c_sum = 0, Y21s_sum = 0;
        double Y22c_sum = 0, Y22s_sum = 0;
        double E_shell = 0;

        for (int i = 2; i < N-2; i++) {
            double x = xcoord(i);
            for (int j = 2; j < N-2; j++) {
                double y = ycoord(j);
                for (int k = 0; k < N; k++) {
                    double z = zcoord(k);
                    double r = sqrt(x*x + y*y + z*z);
                    if (fabs(r - R) > shell_width) continue;

                    long idx = IDX(i, j, k);
                    int km1 = wrap_z(k-1), kp1 = wrap_z(k+1);

                    /* Compute strain energy density:
                     * e_strain = 1/2 * eps_{ij} * eps_{ij}
                     * where eps_{ij} = 1/2(d_i phi_j + d_j phi_i)
                     * Note: phi has 3 components, spatial has 3 dims
                     * eps is 3x3 symmetric tensor */
                    double dphi[3][3]; /* dphi[a][i] = d_i phi_a */
                    for (int a = 0; a < 3; a++) {
                        dphi[a][0] = (phi[a][IDX(i+1,j,k)] - phi[a][IDX(i-1,j,k)]) / (2*dx);
                        dphi[a][1] = (phi[a][IDX(i,j+1,k)] - phi[a][IDX(i,j-1,k)]) / (2*dx);
                        dphi[a][2] = (phi[a][IDX(i,j,kp1)]  - phi[a][IDX(i,j,km1)])  / (2*dx);
                    }

                    /* eps_{ij} = 1/2(d_i phi_j + d_j phi_i) */
                    double eps[3][3];
                    for (int ii2 = 0; ii2 < 3; ii2++)
                        for (int jj2 = 0; jj2 < 3; jj2++)
                            eps[ii2][jj2] = 0.5 * (dphi[jj2][ii2] + dphi[ii2][jj2]);

                    /* Strain energy density */
                    double e_strain = 0;
                    for (int ii2 = 0; ii2 < 3; ii2++)
                        for (int jj2 = 0; jj2 < 3; jj2++)
                            e_strain += eps[ii2][jj2] * eps[ii2][jj2];
                    e_strain *= 0.5;

                    /* Also compute total energy density (KE + grad + pot) */
                    double e_total = 0;
                    for (int a = 0; a < 3; a++) {
                        e_total += 0.5 * vel[a][idx] * vel[a][idx];
                        e_total += 0.5 * (dphi[a][0]*dphi[a][0] + dphi[a][1]*dphi[a][1] + dphi[a][2]*dphi[a][2]);
                    }
                    double p0 = phi[0][idx], p1 = phi[1][idx], p2 = phi[2][idx];
                    double P = p0*p1*p2;
                    double P2 = P*P;
                    e_total += 0.5 * mu_pot * P2 / (1.0 + kappa * P2);

                    /* Use strain energy density for multipole decomposition */
                    double e = e_strain;
                    if (r < 1e-10) continue;

                    double ct = z / r;  /* cos(theta) */
                    double st = sqrt(x*x + y*y) / r;
                    double cp = (st > 1e-10) ? x / (r * st) : 1.0;
                    double sp = (st > 1e-10) ? y / (r * st) : 0.0;

                    /* Real spherical harmonics (unnormalized for power) */
                    /* l=0 */
                    Y00_sum += e;

                    /* l=1 */
                    Y10_sum  += e * ct;
                    Y11c_sum += e * st * cp;
                    Y11s_sum += e * st * sp;

                    /* l=2 */
                    Y20_sum  += e * (3*ct*ct - 1) / 2.0;
                    Y21c_sum += e * st * ct * cp;
                    Y21s_sum += e * st * ct * sp;
                    Y22c_sum += e * st * st * (cp*cp - sp*sp);
                    Y22s_sum += e * st * st * 2*cp*sp;

                    E_shell += e_total;
                    Y00_count += 1;
                }
            }
        }

        if (Y00_count < 10) {
            printf("  R=%.1f: too few points (%.0f), skipping\n", R, Y00_count);
            continue;
        }

        /* Normalize */
        double norm = Y00_count;
        double l0_pow = (Y00_sum / norm) * (Y00_sum / norm);
        double l1_pow = (Y10_sum*Y10_sum + Y11c_sum*Y11c_sum + Y11s_sum*Y11s_sum) / (norm*norm);
        double l2_pow = (Y20_sum*Y20_sum + Y21c_sum*Y21c_sum + Y21s_sum*Y21s_sum
                       + Y22c_sum*Y22c_sum + Y22s_sum*Y22s_sum) / (norm*norm);

        double ratio = (l0_pow > 1e-30) ? l2_pow / l0_pow : 0;

        printf("  R=%5.1f: l=0: %.6e  l=1: %.6e  l=2: %.6e  l2/l0=%.4f  (N_pts=%.0f)\n",
               R, l0_pow, l1_pow, l2_pow, ratio, Y00_count);

        fprintf(fout, "%.2f\t%.8e\t%.8e\t%.8e\t%.6f\t%.8e\n",
                R, l0_pow, l1_pow, l2_pow, ratio, E_shell * dx*dx*dx);
    }

    fclose(fout);
    printf("  M2a complete, results in %s/m2a_strain.tsv\n\n", outdir);
}


/* ──────────────────────────────────────────────────────────────────────
 *  M2b: Two-braid interaction
 *  Two braids separated by D along x, both propagating along z.
 *  Track separation vs time, measure force from acceleration.
 * ─────────────────────────────────────────────────────────────────── */

/* Find centroid of energy density in a half of the domain */
static void find_centroid(double *cx, double *cy, double *cz,
                          double *E_half, double x_min, double x_max)
{
    double sx = 0, sy = 0, sz = 0, se = 0;

    #pragma omp parallel
    {
        double lsx = 0, lsy = 0, lsz = 0, lse = 0;
        #pragma omp for schedule(static) nowait
        for (int i = 2; i < N-2; i++) {
            double x = xcoord(i);
            if (x < x_min || x > x_max) continue;
            for (int j = 2; j < N-2; j++) {
                double y = ycoord(j);
                for (int k = 0; k < N; k++) {
                    double z = zcoord(k);
                    long idx = IDX(i, j, k);
                    int km1 = wrap_z(k-1), kp1 = wrap_z(k+1);

                    double e = 0;
                    for (int a = 0; a < 3; a++) {
                        e += 0.5 * vel[a][idx] * vel[a][idx];
                        double gx = (phi[a][IDX(i+1,j,k)] - phi[a][IDX(i-1,j,k)]) / (2*dx);
                        double gy = (phi[a][IDX(i,j+1,k)] - phi[a][IDX(i,j-1,k)]) / (2*dx);
                        double gz = (phi[a][IDX(i,j,kp1)]  - phi[a][IDX(i,j,km1)])  / (2*dx);
                        e += 0.5 * (gx*gx + gy*gy + gz*gz);
                    }
                    double p0 = phi[0][idx], p1 = phi[1][idx], p2 = phi[2][idx];
                    double P = p0*p1*p2;
                    double P2 = P*P;
                    e += 0.5 * mu_pot * P2 / (1.0 + kappa * P2);
                    if (e < 0) e = -e;  /* use |e| for centroid weight */

                    lsx += e * x;
                    lsy += e * y;
                    lsz += e * z;
                    lse += e;
                }
            }
        }
        #pragma omp critical
        {
            sx += lsx; sy += lsy; sz += lsz; se += lse;
        }
    }

    *cx = (se > 1e-30) ? sx / se : 0;
    *cy = (se > 1e-30) ? sy / se : 0;
    *cz = (se > 1e-30) ? sz / se : 0;
    *E_half = se * dx * dx * dx;
}

/* Compute peak |P| */
static double compute_peakP(void)
{
    double pk = 0;
    #pragma omp parallel
    {
        double lpk = 0;
        #pragma omp for schedule(static) nowait
        for (long idx = 0; idx < Ngrid; idx++) {
            double P = fabs(phi[0][idx] * phi[1][idx] * phi[2][idx]);
            if (P > lpk) lpk = P;
        }
        #pragma omp critical
        { if (lpk > pk) pk = lpk; }
    }
    return pk;
}

/* Compute total energy */
static double compute_total_energy(void)
{
    double E = 0;
    #pragma omp parallel
    {
        double lE = 0;
        #pragma omp for schedule(static) nowait
        for (int i = 2; i < N-2; i++) {
            for (int j = 2; j < N-2; j++) {
                for (int k = 0; k < N; k++) {
                    long idx = IDX(i, j, k);
                    int km1 = wrap_z(k-1), kp1 = wrap_z(k+1);
                    double e = 0;
                    for (int a = 0; a < 3; a++) {
                        e += 0.5 * vel[a][idx] * vel[a][idx];
                        double gx = (phi[a][IDX(i+1,j,k)] - phi[a][IDX(i-1,j,k)]) / (2*dx);
                        double gy = (phi[a][IDX(i,j+1,k)] - phi[a][IDX(i,j-1,k)]) / (2*dx);
                        double gz = (phi[a][IDX(i,j,kp1)]  - phi[a][IDX(i,j,km1)])  / (2*dx);
                        e += 0.5 * (gx*gx + gy*gy + gz*gz);
                    }
                    double p0 = phi[0][idx], p1 = phi[1][idx], p2 = phi[2][idx];
                    double P = p0*p1*p2;
                    double P2 = P*P;
                    e += 0.5 * mu_pot * P2 / (1.0 + kappa * P2);
                    lE += e;
                }
            }
        }
        #pragma omp critical
        { E += lE; }
    }
    return E * dx * dx * dx;
}

static void run_m2b(void)
{
    printf("================================================================\n");
    printf("  M2b: Two-Braid Interaction (D=30 along x)\n");
    printf("================================================================\n\n");

    /* Need larger domain for two braids at D=30 */
    double L_save = L;
    L = 40.0;
    dx = 2.0 * L / (N - 1);
    dx2 = dx * dx;
    dt = cfl_frac * dx;

    printf("  Two-braid domain: L=%.0f, dx=%.4f, dt=%.5f\n", L, dx, dt);

    /* Re-allocate damping for new domain */
    setup_damping_rect();

    double D = 30.0;  /* initial separation */
    double cx1 = -D/2.0, cx2 = +D/2.0;

    for (int a = 0; a < 3; a++) {
        memset(phi[a], 0, Ngrid * sizeof(double));
        memset(vel[a], 0, Ngrid * sizeof(double));
        memset(acc[a], 0, Ngrid * sizeof(double));
    }

    add_braid(cx1, 0.0, +1.0);
    add_braid(cx2, 0.0, +1.0);

    char path[600];
    snprintf(path, sizeof(path), "%s/m2b_twobraid.tsv", outdir);
    FILE *fout = fopen(path, "w");
    if (!fout) { fprintf(stderr, "Cannot open %s\n", path); return; }
    fprintf(fout, "time\tx1\ty1\tx2\ty2\tsep\tE_total\tpeak_P\n");

    compute_acc();

    int Nt = (int)(tfinal / dt) + 1;
    int diag_every = Nt / 1000;
    if (diag_every < 1) diag_every = 1;
    int print_every = Nt / 20;
    if (print_every < 1) print_every = 1;

    double wall_start = omp_get_wtime();

    printf("  Evolving two braids for t=%.0f...\n", tfinal);
    fflush(stdout);

    for (int n = 0; n <= Nt; n++) {
        double t = n * dt;
        int do_diag  = (n % diag_every == 0);
        int do_print = (n % print_every == 0);

        if (do_diag || do_print) {
            double x1, y1, z1, E1;
            double x2, y2, z2, E2;
            find_centroid(&x1, &y1, &z1, &E1, -L, 0.0);
            find_centroid(&x2, &y2, &z2, &E2, 0.0, L);
            double sep = sqrt((x2-x1)*(x2-x1) + (y2-y1)*(y2-y1));
            double peakP = compute_peakP();
            double Etot = E1 + E2;  /* approximate */

            if (do_diag)
                fprintf(fout, "%.4f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6e\t%.6e\n",
                        t, x1, y1, x2, y2, sep, Etot, peakP);

            if (do_print) {
                double elapsed = omp_get_wtime() - wall_start;
                double frac = (double)n / Nt;
                double eta = (frac > 0.001) ? elapsed * (1.0-frac)/frac : 0;
                printf("  t=%7.1f  sep=%.4f  x1=%.3f x2=%.3f  |P|=%.4f  [%.0fs, ETA %.0fs]\n",
                       t, sep, x1, x2, peakP, elapsed, eta);
                fflush(stdout);
            }
        }

        if (n == Nt) break;
        verlet_step();
    }

    fclose(fout);
    printf("  M2b complete, results in %s/m2b_twobraid.tsv\n\n", outdir);

    /* Restore domain */
    L = L_save;
    dx = 2.0 * L / (N - 1);
    dx2 = dx * dx;
    dt = cfl_frac * dx;
}


/* ──────────────────────────────────────────────────────────────────────
 *  M2c: Angular pattern of force
 *  Place second braid at different angles relative to first.
 *  Track separation change (attraction/repulsion) for each angle.
 * ─────────────────────────────────────────────────────────────────── */

static void run_m2c_angle(double angle_deg, FILE *fsum)
{
    double L_save = L;
    L = 40.0;
    dx = 2.0 * L / (N - 1);
    dx2 = dx * dx;
    dt = cfl_frac * dx;

    setup_damping_rect();

    double D = 30.0;
    double angle_rad = angle_deg * M_PI / 180.0;
    double cx2 = D * cos(angle_rad);
    double cy2 = D * sin(angle_rad);

    for (int a = 0; a < 3; a++) {
        memset(phi[a], 0, Ngrid * sizeof(double));
        memset(vel[a], 0, Ngrid * sizeof(double));
        memset(acc[a], 0, Ngrid * sizeof(double));
    }

    add_braid(0.0, 0.0, +1.0);
    add_braid(cx2, cy2, +1.0);

    char label[64];
    snprintf(label, sizeof(label), "m2c_angle%.0f", angle_deg);
    char path[600];
    snprintf(path, sizeof(path), "%s/%s.tsv", outdir, label);
    FILE *fout = fopen(path, "w");
    if (!fout) { fprintf(stderr, "Cannot open %s\n", path); return; }
    fprintf(fout, "time\tsep\tpeak_P\n");

    compute_acc();

    double t_run = 200.0;  /* shorter run for angular scan */
    int Nt = (int)(t_run / dt) + 1;
    int diag_every = Nt / 500;
    if (diag_every < 1) diag_every = 1;
    int print_every = Nt / 10;
    if (print_every < 1) print_every = 1;

    double wall_start = omp_get_wtime();

    /* Track separation at start and end */
    double sep_initial = -1, sep_final = -1;

    printf("  angle=%.0f deg (braid2 at %.1f, %.1f)...\n", angle_deg, cx2, cy2);
    fflush(stdout);

    for (int n = 0; n <= Nt; n++) {
        double t = n * dt;
        int do_diag = (n % diag_every == 0);
        int do_print = (n % print_every == 0);

        if (do_diag || do_print || n == 0 || n == Nt) {
            /* Find centroids: braid 1 near origin, braid 2 near (cx2, cy2) */
            /* Split domain along the perpendicular bisector of the line connecting them */
            double mx = cx2 / 2.0, my = cy2 / 2.0;

            /* Use simple half-plane split: dot product with direction vector */
            double dx_dir = cx2, dy_dir = cy2;
            double dnorm = sqrt(dx_dir*dx_dir + dy_dir*dy_dir);
            if (dnorm < 1e-10) { dx_dir = 1; dy_dir = 0; dnorm = 1; }
            dx_dir /= dnorm; dy_dir /= dnorm;

            double sx1 = 0, sy1 = 0, se1 = 0;
            double sx2 = 0, sy2 = 0, se2 = 0;

            #pragma omp parallel
            {
                double l1x = 0, l1y = 0, l1e = 0;
                double l2x = 0, l2y = 0, l2e = 0;
                #pragma omp for schedule(static) nowait
                for (int i = 2; i < N-2; i++) {
                    double x = xcoord(i);
                    for (int j = 2; j < N-2; j++) {
                        double y = ycoord(j);
                        /* Which side of midpoint? */
                        double proj = (x - mx) * dx_dir + (y - my) * dy_dir;
                        for (int k = 0; k < N; k++) {
                            long idx = IDX(i, j, k);
                            int km1 = wrap_z(k-1), kp1 = wrap_z(k+1);
                            double e = 0;
                            for (int a = 0; a < 3; a++) {
                                e += 0.5 * vel[a][idx] * vel[a][idx];
                                double gx2 = (phi[a][IDX(i+1,j,k)] - phi[a][IDX(i-1,j,k)]) / (2*dx);
                                double gy2 = (phi[a][IDX(i,j+1,k)] - phi[a][IDX(i,j-1,k)]) / (2*dx);
                                double gz2 = (phi[a][IDX(i,j,kp1)]  - phi[a][IDX(i,j,km1)])  / (2*dx);
                                e += 0.5 * (gx2*gx2 + gy2*gy2 + gz2*gz2);
                            }
                            double p0 = phi[0][idx], p1 = phi[1][idx], p2 = phi[2][idx];
                            double P = p0*p1*p2;
                            double P2 = P*P;
                            e += 0.5 * mu_pot * P2 / (1.0 + kappa * P2);
                            if (e < 0) e = -e;

                            if (proj < 0) {
                                l1x += e * x; l1y += e * y; l1e += e;
                            } else {
                                l2x += e * x; l2y += e * y; l2e += e;
                            }
                        }
                    }
                }
                #pragma omp critical
                {
                    sx1 += l1x; sy1 += l1y; se1 += l1e;
                    sx2 += l2x; sy2 += l2y; se2 += l2e;
                }
            }

            double xx1 = (se1 > 1e-30) ? sx1/se1 : 0;
            double yy1 = (se1 > 1e-30) ? sy1/se1 : 0;
            double xx2 = (se2 > 1e-30) ? sx2/se2 : cx2;
            double yy2 = (se2 > 1e-30) ? sy2/se2 : cy2;
            double sep = sqrt((xx2-xx1)*(xx2-xx1) + (yy2-yy1)*(yy2-yy1));
            double peakP = compute_peakP();

            if (n == 0) sep_initial = sep;
            sep_final = sep;

            if (do_diag)
                fprintf(fout, "%.4f\t%.6f\t%.6e\n", t, sep, peakP);

            if (do_print) {
                double elapsed = omp_get_wtime() - wall_start;
                double frac = (double)n / Nt;
                double eta = (frac > 0.001) ? elapsed * (1.0-frac)/frac : 0;
                printf("    t=%6.1f  sep=%.4f  |P|=%.4f  [%.0fs, ETA %.0fs]\n",
                       t, sep, peakP, elapsed, eta);
                fflush(stdout);
            }
        }

        if (n == Nt) break;
        verlet_step();
    }

    fclose(fout);

    double dsep = sep_final - sep_initial;
    double force_proxy = -dsep / t_run;  /* positive = attraction */

    printf("  angle=%.0f: sep_i=%.4f -> sep_f=%.4f, dsep=%.4f, force_proxy=%.6f\n",
           angle_deg, sep_initial, sep_final, dsep, force_proxy);

    if (fsum)
        fprintf(fsum, "%.1f\t%.6f\t%.6f\t%.6f\t%.6e\n",
                angle_deg, sep_initial, sep_final, dsep, force_proxy);

    /* Restore domain */
    L = L_save;
    dx = 2.0 * L / (N - 1);
    dx2 = dx * dx;
    dt = cfl_frac * dx;
}

static void run_m2c(void)
{
    printf("================================================================\n");
    printf("  M2c: Angular Pattern of Force\n");
    printf("================================================================\n\n");

    char path[600];
    snprintf(path, sizeof(path), "%s/m2c_angular.tsv", outdir);
    FILE *fsum = fopen(path, "w");
    if (fsum)
        fprintf(fsum, "angle_deg\tsep_initial\tsep_final\tdsep\tforce_proxy\n");

    /* Test angles: 0, 45, 90 degrees */
    double angles[] = {0.0, 45.0, 90.0};
    int n_angles = 3;

    for (int ia = 0; ia < n_angles; ia++) {
        run_m2c_angle(angles[ia], fsum);
    }

    if (fsum) fclose(fsum);
    printf("  M2c complete, results in %s/m2c_angular.tsv\n\n", outdir);
}


/* ─── Main ─── */
int main(int argc, char **argv)
{
    int test = -1;  /* -1=all, 0=M2a, 1=M2b, 2=M2c */

    for (int i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "-test") && i+1 < argc)   test = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-N")  && i+1 < argc) N = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-L")  && i+1 < argc) L = atof(argv[++i]);
        else if (!strcmp(argv[i], "-o")  && i+1 < argc) strncpy(outdir, argv[++i], sizeof(outdir)-1);
        else if (!strcmp(argv[i], "-tfinal") && i+1 < argc) tfinal = atof(argv[++i]);
        else if (!strcmp(argv[i], "-cfl") && i+1 < argc) cfl_frac = atof(argv[++i]);
        else if (!strcmp(argv[i], "-mu")  && i+1 < argc) mu_pot = atof(argv[++i]);
        else if (!strcmp(argv[i], "-kappa") && i+1 < argc) kappa = atof(argv[++i]);
        else if (!strcmp(argv[i], "-A0") && i+1 < argc) A0 = atof(argv[++i]);
        else if (!strcmp(argv[i], "-R")  && i+1 < argc) R_tube = atof(argv[++i]);
    }

    dx  = 2.0 * L / (N - 1);
    dx2 = dx * dx;
    dt  = cfl_frac * dx;
    Ngrid = (long)N * N * N;

    printf("=== V27-M2: Gravity Between Massless Propagating Braids ===\n");
    printf("  N=%d  L=%.1f  dx=%.4f  dt=%.5f  Ngrid=%.1fM\n",
           N, L, dx, dt, Ngrid/1e6);
    printf("  mu=%.1f  kappa=%.1f  A0=%.2f  R_tube=%.1f  m=0 (massless)\n",
           mu_pot, kappa, A0, R_tube);
    printf("  Threads: %d\n", omp_get_max_threads());
    printf("  tfinal=%.0f  test=%d (-1=all)\n\n", tfinal, test);
    fflush(stdout);

    /* Allocate */
    for (int a = 0; a < 3; a++) {
        phi[a] = calloc(Ngrid, sizeof(double));
        vel[a] = calloc(Ngrid, sizeof(double));
        acc[a] = calloc(Ngrid, sizeof(double));
        if (!phi[a] || !vel[a] || !acc[a]) {
            fprintf(stderr, "Allocation failed\n"); return 1;
        }
    }
    damp = malloc(Ngrid * sizeof(double));
    if (!damp) { fprintf(stderr, "Damp alloc failed\n"); return 1; }

    if (test == -1 || test == 0) run_m2a();
    if (test == -1 || test == 1) run_m2b();
    if (test == -1 || test == 2) run_m2c();

    /* Cleanup */
    for (int a = 0; a < 3; a++) {
        free(phi[a]); free(vel[a]); free(acc[a]);
    }
    free(damp);

    printf("\n=== V27-M2 Complete ===\n");
    return 0;
}
