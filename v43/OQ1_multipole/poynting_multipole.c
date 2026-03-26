/*  poynting_multipole.c — Spherical harmonic decomposition of time-averaged
 *  Poynting flux |S| = |E × B| on spherical shells.
 *
 *  E = -dtheta/dt  (finite difference between consecutive frames)
 *  B = curl(theta)  (centered spatial differences)
 *  S = E × B        (Cosserat Poynting vector)
 *
 *  Accumulates <|S|> over all frame pairs, then decomposes into Y_l^m
 *  on spherical shells at specified radii.
 *
 *  Build: gcc -O3 -march=native -fopenmp -o poynting_multipole poynting_multipole.c -lzstd -lm
 */

#define SFA_IMPLEMENTATION
#include "../../sfa/format/sfa.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <omp.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ---- f16 to f64 conversion ---- */
static inline double f16_to_f64(uint16_t h) {
    uint16_t sign = h & 0x8000;
    int exp = (h >> 10) & 0x1F;
    uint16_t mant = h & 0x3FF;
    if (exp == 0) return 0.0;
    if (exp == 31) return sign ? -1e30 : 1e30;
    float fv; uint32_t x = ((uint32_t)sign << 16) | ((uint32_t)(exp-15+127) << 23) | ((uint32_t)mant << 13);
    memcpy(&fv, &x, 4); return (double)fv;
}

/* ---- Associated Legendre polynomials P_l^m(x), m >= 0 ---- */
static double assoc_legendre(int l, int m, double x) {
    if (m < 0 || m > l) return 0.0;
    double pmm = 1.0;
    if (m > 0) {
        double sqx = sqrt(1.0 - x*x);
        double fact = 1.0;
        for (int i = 1; i <= m; i++) {
            pmm *= -fact * sqx;
            fact += 2.0;
        }
    }
    if (l == m) return pmm;
    double pm1m = x * (2*m + 1) * pmm;
    if (l == m + 1) return pm1m;
    double pll = 0.0;
    for (int ll = m + 2; ll <= l; ll++) {
        pll = (x * (2*ll - 1) * pm1m - (ll + m - 1) * pmm) / (ll - m);
        pmm = pm1m;
        pm1m = pll;
    }
    return pll;
}

static double factorial(int n) {
    double f = 1.0;
    for (int i = 2; i <= n; i++) f *= i;
    return f;
}

/* Real spherical harmonic Y_l^m(theta, phi)
 * m > 0: cos(m*phi), m < 0: sin(|m|*phi), m = 0: standard
 * Normalized: integral |Y_lm|^2 dOmega = 1
 */
static double real_ylm(int l, int m, double theta, double phi_ang) {
    int am = abs(m);
    double norm = sqrt((2.0*l + 1.0) / (4.0 * M_PI) * factorial(l - am) / factorial(l + am));
    double plm = assoc_legendre(l, am, cos(theta));
    if (m > 0)      return norm * sqrt(2.0) * plm * cos(m * phi_ang);
    else if (m < 0) return norm * sqrt(2.0) * plm * sin(am * phi_ang);
    else            return norm * plm;
}

/* ---- Trilinear interpolation ---- */
static double interp3d(const double *field, int N, double L, double dx,
                       double x, double y, double z) {
    double gx = (x + L) / dx;
    double gy = (y + L) / dx;
    double gz = (z + L) / dx;
    int ix = (int)floor(gx), iy = (int)floor(gy), iz = (int)floor(gz);
    if (ix < 0 || ix >= N-1 || iy < 0 || iy >= N-1 || iz < 0 || iz >= N-1) return 0.0;
    double fx = gx - ix, fy = gy - iy, fz = gz - iz;
    int NN = N * N;
    double c000 = field[ix*NN + iy*N + iz];
    double c001 = field[ix*NN + iy*N + (iz+1)];
    double c010 = field[ix*NN + (iy+1)*N + iz];
    double c011 = field[ix*NN + (iy+1)*N + (iz+1)];
    double c100 = field[(ix+1)*NN + iy*N + iz];
    double c101 = field[(ix+1)*NN + iy*N + (iz+1)];
    double c110 = field[(ix+1)*NN + (iy+1)*N + iz];
    double c111 = field[(ix+1)*NN + (iy+1)*N + (iz+1)];
    return c000*(1-fx)*(1-fy)*(1-fz) + c001*(1-fx)*(1-fy)*fz
         + c010*(1-fx)*fy*(1-fz) + c011*(1-fx)*fy*fz
         + c100*fx*(1-fy)*(1-fz) + c101*fx*(1-fy)*fz
         + c110*fx*fy*(1-fz) + c111*fx*fy*fz;
}

#define LMAX 6
#define NTHETA 90
#define NPHI 180

/* Extract field data from a frame buffer, handling different dtypes */
static void extract_fields(SFA *sfa, void *buf, long N3,
                           double *phi[3], double *theta[3]) {
    uint64_t off = 0;
    for (int c = 0; c < (int)sfa->n_columns; c++) {
        int dtype = sfa->columns[c].dtype;
        int sem = sfa->columns[c].semantic;
        int comp = sfa->columns[c].component;
        int es = sfa_dtype_size[dtype];
        uint8_t *src = (uint8_t *)buf + off;

        double *target = NULL;
        if (sem == SFA_POSITION && comp < 3) target = phi[comp];
        else if (sem == SFA_ANGLE && comp < 3) target = theta[comp];

        if (target) {
            if (dtype == SFA_F64)
                for (long i = 0; i < N3; i++) target[i] = ((double*)src)[i];
            else if (dtype == SFA_F32)
                for (long i = 0; i < N3; i++) target[i] = (double)((float*)src)[i];
            else if (dtype == SFA_F16)
                for (long i = 0; i < N3; i++) target[i] = f16_to_f64(((uint16_t*)src)[i]);
        }
        off += (uint64_t)N3 * es;
    }
}

int main(int argc, char **argv) {
    if (argc < 2) {
        fprintf(stderr, "Usage: %s input.sfa\n", argv[0]);
        return 1;
    }

    SFA *sfa = sfa_open(argv[1]);
    if (!sfa) { fprintf(stderr, "Cannot open %s\n", argv[1]); return 1; }

    int N = sfa->Nx;
    double L = sfa->Lx;
    long N3 = (long)N * N * N;
    int NN = N * N;
    double dx = 2.0 * L / (N - 1);
    int nframes = sfa->total_frames;

    printf("# Poynting flux multipole analysis: %s\n", argv[1]);
    printf("# Grid: N=%d, L=%.1f, dx=%.4f, frames=%d\n", N, L, dx, nframes);

    if (nframes < 2) {
        fprintf(stderr, "Need at least 2 frames for time derivatives\n");
        return 1;
    }

    /* Shell radii */
    double R_list[] = {3.0, 5.0, 8.0, 10.0, 12.0};
    int n_R = 5;

    /* Filter out radii that exceed domain */
    int valid_R[5];
    int n_valid = 0;
    for (int ri = 0; ri < n_R; ri++) {
        if (R_list[ri] <= L - 2*dx) {
            valid_R[n_valid++] = ri;
        } else {
            printf("# R=%.1f too large for domain L=%.1f, skipping\n", R_list[ri], L);
        }
    }

    /* Allocate time-averaged |S| on angular grids for each shell */
    /* avg_S[ri][it * NPHI + ip] = accumulated |S| at that angle */
    double **avg_S = (double **)malloc(n_valid * sizeof(double *));
    for (int i = 0; i < n_valid; i++) {
        avg_S[i] = (double *)calloc(NTHETA * NPHI, sizeof(double));
    }

    /* Also track total flux per shell for 1/R^2 check */
    double *total_flux = (double *)calloc(n_valid, sizeof(double));

    /* Allocate frame buffers */
    void *buf_cur = malloc(sfa->frame_bytes);
    void *buf_nxt = malloc(sfa->frame_bytes);
    if (!buf_cur || !buf_nxt) {
        fprintf(stderr, "Out of memory for frame buffers (%lu bytes each)\n", sfa->frame_bytes);
        return 1;
    }

    /* Field arrays (reused per frame pair) */
    double *phi_cur[3], *theta_cur[3], *theta_nxt[3], *phi_nxt[3];
    for (int a = 0; a < 3; a++) {
        phi_cur[a]   = (double *)malloc(N3 * sizeof(double));
        theta_cur[a] = (double *)malloc(N3 * sizeof(double));
        phi_nxt[a]   = (double *)malloc(N3 * sizeof(double));
        theta_nxt[a] = (double *)malloc(N3 * sizeof(double));
    }

    /* |S| field on the grid (reused per frame pair) */
    double *S_mag = (double *)malloc(N3 * sizeof(double));

    /* Centroid accumulator (time-averaged, |P|-weighted) */
    double cm_avg[3] = {0, 0, 0};
    int n_pairs = 0;

    /* Read first frame */
    printf("# Reading frame 0 ...\n");
    sfa_read_frame(sfa, 0, buf_cur);
    extract_fields(sfa, buf_cur, N3, phi_cur, theta_cur);
    double t_cur = sfa_frame_time(sfa, 0);

    /* Process consecutive frame pairs */
    for (int fn = 0; fn < nframes - 1; fn++) {
        /* Read next frame */
        sfa_read_frame(sfa, fn + 1, buf_nxt);
        extract_fields(sfa, buf_nxt, N3, phi_nxt, theta_nxt);
        double t_nxt = sfa_frame_time(sfa, fn + 1);
        double dt_frame = t_nxt - t_cur;

        if (dt_frame <= 0) {
            fprintf(stderr, "# WARNING: non-positive dt at frame %d (t=%.4f -> %.4f), skipping\n",
                    fn, t_cur, t_nxt);
            /* Swap buffers */
            goto swap;
        }

        if (fn % 50 == 0)
            printf("# Processing frame pair %d/%d (t=%.2f, dt=%.4f) ...\n",
                   fn, nframes - 1, t_cur, dt_frame);

        /* Compute centroid from phi fields (average of cur and nxt) */
        double cm[3] = {0, 0, 0}, wt = 0;
        for (long idx = 0; idx < N3; idx++) {
            double P = fabs(phi_cur[0][idx] * phi_cur[1][idx] * phi_cur[2][idx]);
            int i = (int)(idx / NN), j = (int)((idx / N) % N), k = (int)(idx % N);
            double x = -L + i * dx, y = -L + j * dx, z = -L + k * dx;
            cm[0] += x * P; cm[1] += y * P; cm[2] += z * P; wt += P;
        }
        if (wt > 0) { cm[0] /= wt; cm[1] /= wt; cm[2] /= wt; }
        cm_avg[0] += cm[0]; cm_avg[1] += cm[1]; cm_avg[2] += cm[2];

        /* Compute E = -dtheta/dt and B = curl(theta) at each grid point,
         * then |S| = |E x B| */
        #pragma omp parallel for schedule(static)
        for (long idx = 0; idx < N3; idx++) {
            int i = (int)(idx / NN), j = (int)((idx / N) % N), k = (int)(idx % N);

            /* E = -dtheta/dt */
            double E[3];
            for (int a = 0; a < 3; a++)
                E[a] = -(theta_nxt[a][idx] - theta_cur[a][idx]) / dt_frame;

            /* B = curl(theta) at current frame using centered differences */
            /* With periodic boundary handling */
            int ip = (i < N-1) ? i+1 : i;
            int im = (i > 0) ? i-1 : i;
            int jp = (j < N-1) ? j+1 : j;
            int jm = (j > 0) ? j-1 : j;
            int kp = (k < N-1) ? k+1 : k;
            int km = (k > 0) ? k-1 : k;

            /* Denominators for finite differences (handle boundaries) */
            double dxi = (ip - im) * dx;
            double dyi = (jp - jm) * dx;
            double dzi = (kp - km) * dx;
            if (dxi < 1e-15) dxi = dx;
            if (dyi < 1e-15) dyi = dx;
            if (dzi < 1e-15) dzi = dx;

            /* d(theta_a)/dx = (theta_a[ip,...] - theta_a[im,...]) / dxi */
            /* B_x = dtheta_z/dy - dtheta_y/dz */
            double dtz_dy = (theta_cur[2][ip*NN + jp*N + k] - theta_cur[2][ip*NN + jm*N + k]);
            /* Wait, need proper indexing. Let me fix: */
            /* theta_cur[a][i*NN + j*N + k] */

            double dthetaz_dy = (theta_cur[2][i*NN + jp*N + k] - theta_cur[2][i*NN + jm*N + k]) / dyi;
            double dthetay_dz = (theta_cur[1][i*NN + j*N + kp] - theta_cur[1][i*NN + j*N + km]) / dzi;

            double dthetax_dz = (theta_cur[0][i*NN + j*N + kp] - theta_cur[0][i*NN + j*N + km]) / dzi;
            double dthetaz_dx = (theta_cur[2][ip*NN + j*N + k] - theta_cur[2][im*NN + j*N + k]) / dxi;

            double dthetay_dx = (theta_cur[1][ip*NN + j*N + k] - theta_cur[1][im*NN + j*N + k]) / dxi;
            double dthetax_dy = (theta_cur[0][i*NN + jp*N + k] - theta_cur[0][i*NN + jm*N + k]) / dyi;

            double B[3];
            B[0] = dthetaz_dy - dthetay_dz;
            B[1] = dthetax_dz - dthetaz_dx;
            B[2] = dthetay_dx - dthetax_dy;

            /* S = E x B */
            double Sx = E[1]*B[2] - E[2]*B[1];
            double Sy = E[2]*B[0] - E[0]*B[2];
            double Sz = E[0]*B[1] - E[1]*B[0];

            S_mag[idx] = sqrt(Sx*Sx + Sy*Sy + Sz*Sz);
        }

        /* Sample |S| on angular grids at each shell radius */
        for (int vi = 0; vi < n_valid; vi++) {
            int ri = valid_R[vi];
            double R = R_list[ri];
            double dtheta = M_PI / NTHETA;
            double dphi = 2.0 * M_PI / NPHI;

            for (int it = 0; it < NTHETA; it++) {
                double th = (it + 0.5) * dtheta;
                double sth = sin(th);
                for (int ip = 0; ip < NPHI; ip++) {
                    double ph = (ip + 0.5) * dphi;
                    double xs = cm[0] + R * sth * cos(ph);
                    double ys = cm[1] + R * sth * sin(ph);
                    double zs = cm[2] + R * cos(th);
                    double val = interp3d(S_mag, N, L, dx, xs, ys, zs);
                    avg_S[vi][it * NPHI + ip] += val;
                }
            }
        }

        n_pairs++;

    swap:
        /* Swap cur <-> nxt */
        t_cur = t_nxt;
        for (int a = 0; a < 3; a++) {
            double *tmp;
            tmp = phi_cur[a]; phi_cur[a] = phi_nxt[a]; phi_nxt[a] = tmp;
            tmp = theta_cur[a]; theta_cur[a] = theta_nxt[a]; theta_nxt[a] = tmp;
        }
        /* buf_nxt becomes buf_cur for next iteration */
        void *tmp_buf = buf_cur; buf_cur = buf_nxt; buf_nxt = tmp_buf;
    }

    printf("# Processed %d frame pairs\n", n_pairs);

    if (n_pairs == 0) {
        fprintf(stderr, "No valid frame pairs processed\n");
        return 1;
    }

    /* Normalize time average */
    cm_avg[0] /= n_pairs; cm_avg[1] /= n_pairs; cm_avg[2] /= n_pairs;
    for (int vi = 0; vi < n_valid; vi++)
        for (int k = 0; k < NTHETA * NPHI; k++)
            avg_S[vi][k] /= n_pairs;

    printf("# Time-averaged centroid: (%.3f, %.3f, %.3f)\n",
           cm_avg[0], cm_avg[1], cm_avg[2]);

    /* ---- Spherical harmonic decomposition ---- */

    /* Open output files */
    FILE *ftsv = fopen("poynting_multipole.tsv", "w");
    fprintf(ftsv, "R\tl\tC_l\tfraction\ttotal_power\ttotal_flux\n");

    FILE *fmd = fopen("poynting_multipole_results.md", "w");
    fprintf(fmd, "# Poynting Flux Multipole Decomposition — Single Braid (Null Model)\n\n");
    fprintf(fmd, "**Source**: `%s`\n", argv[1]);
    fprintf(fmd, "**Grid**: N=%d, L=%.1f, dx=%.4f, frames=%d\n", N, L, dx, nframes);
    fprintf(fmd, "**Frame pairs processed**: %d\n", n_pairs);
    fprintf(fmd, "**Time-averaged centroid**: (%.3f, %.3f, %.3f)\n\n", cm_avg[0], cm_avg[1], cm_avg[2]);

    fprintf(fmd, "## Method\n\n");
    fprintf(fmd, "For each consecutive frame pair (n, n+1):\n");
    fprintf(fmd, "1. E = -dtheta/dt via finite difference\n");
    fprintf(fmd, "2. B = curl(theta) via centered spatial differences\n");
    fprintf(fmd, "3. S = E x B (Cosserat Poynting vector)\n");
    fprintf(fmd, "4. |S| computed at each grid point\n\n");
    fprintf(fmd, "Time-averaged <|S|> accumulated over %d frame pairs.\n", n_pairs);
    fprintf(fmd, "Decomposed into real spherical harmonics Y_l^m for l=0..%d.\n\n", LMAX);
    fprintf(fmd, "Angular grid: %d x %d (theta x phi), trapezoidal integration.\n\n", NTHETA, NPHI);

    fprintf(fmd, "## Power Spectrum\n\n");
    fprintf(fmd, "C_l = (1/(2l+1)) sum_m |c_lm|^2\n\n");

    printf("\n# Multipole power spectrum: C_l = (1/(2l+1)) sum_m |c_lm|^2\n");
    printf("# R      ");
    for (int l = 0; l <= LMAX; l++) printf("  l=%-6d", l);
    printf("  f(l=0)   f(l=1)   f(l=2)   total_flux\n");

    fprintf(fmd, "| R | l=0 | l=1 | l=2 | l=3 | l=4 | l=5 | l=6 | f(l=0) | f(l=1) | f(l=2) | flux | dominant |\n");
    fprintf(fmd, "|---|-----|-----|-----|-----|-----|-----|-----|--------|--------|--------|------|----------|\n");

    double prev_flux = -1;

    for (int vi = 0; vi < n_valid; vi++) {
        int ri = valid_R[vi];
        double R = R_list[ri];
        double dtheta = M_PI / NTHETA;
        double dphi = 2.0 * M_PI / NPHI;

        /* Compute c_lm coefficients */
        int n_coeff = (LMAX + 1) * (LMAX + 1);
        double *c_lm = (double *)calloc(n_coeff, sizeof(double));

        /* Also compute total flux integral <|S|> dOmega */
        double flux_integral = 0;

        for (int it = 0; it < NTHETA; it++) {
            double th = (it + 0.5) * dtheta;
            double sth = sin(th);
            for (int ip = 0; ip < NPHI; ip++) {
                double ph = (ip + 0.5) * dphi;
                double val = avg_S[vi][it * NPHI + ip];
                double dA = sth * dtheta * dphi;  /* solid angle element */

                flux_integral += val * dA;

                for (int l = 0; l <= LMAX; l++) {
                    for (int m = -l; m <= l; m++) {
                        double ylm = real_ylm(l, m, th, ph);
                        c_lm[l*l + l + m] += val * ylm * dA;
                    }
                }
            }
        }

        /* Compute power spectrum C_l = (1/(2l+1)) sum_m |c_lm|^2 */
        double C[LMAX + 1];
        double C_total = 0;
        for (int l = 0; l <= LMAX; l++) {
            C[l] = 0;
            for (int m = -l; m <= l; m++) {
                double c = c_lm[l*l + l + m];
                C[l] += c * c;
            }
            C[l] /= (2*l + 1);
            C_total += C[l];
        }

        /* Also multiply flux by R^2 to check 1/R^2 scaling */
        double flux_R2 = flux_integral * R * R;

        /* Print results */
        printf("  %4.1f  ", R);
        for (int l = 0; l <= LMAX; l++) printf("%9.4e ", C[l]);
        printf("%8.4f %8.4f %8.4f %10.4e",
               C_total > 0 ? C[0]/C_total : 0,
               C_total > 0 ? C[1]/C_total : 0,
               C_total > 0 ? C[2]/C_total : 0,
               flux_integral);
        if (prev_flux > 0) {
            double ratio = flux_integral / prev_flux;
            printf("  (ratio=%.3f)", ratio);
        }
        printf("\n");

        /* Determine dominant mode */
        int dom = 0;
        for (int l = 1; l <= LMAX; l++) if (C[l] > C[dom]) dom = l;
        const char *names[] = {"monopole","dipole","quadrupole","octupole",
                               "hexadecapole","l=5","l=6"};

        /* TSV output */
        for (int l = 0; l <= LMAX; l++) {
            fprintf(ftsv, "%.1f\t%d\t%.6e\t%.6f\t%.6e\t%.6e\n",
                    R, l, C[l], C_total > 0 ? C[l]/C_total : 0, C_total, flux_integral);
        }

        /* Markdown */
        fprintf(fmd, "| %.1f ", R);
        for (int l = 0; l <= LMAX; l++)
            fprintf(fmd, "| %.4e ", C[l]);
        fprintf(fmd, "| %.1f%% | %.1f%% | %.1f%% ",
                C_total > 0 ? 100.0*C[0]/C_total : 0,
                C_total > 0 ? 100.0*C[1]/C_total : 0,
                C_total > 0 ? 100.0*C[2]/C_total : 0);
        fprintf(fmd, "| %.4e | %s |\n", flux_integral, names[dom]);

        /* Detailed a_lm breakdown */
        printf("#   c_lm detail for R=%.1f:\n", R);
        for (int l = 0; l <= LMAX; l++) {
            printf("#     l=%d:", l);
            for (int m = -l; m <= l; m++) {
                printf(" m=%+d:%.3e", m, c_lm[l*l + l + m]);
            }
            printf("\n");
        }

        prev_flux = flux_integral;
        free(c_lm);
    }

    /* 1/R^2 scaling analysis */
    fprintf(fmd, "\n## Flux vs Radius (1/R^2 Test)\n\n");
    fprintf(fmd, "If radiative, total flux F(R) should scale as 1/R^2.\n");
    fprintf(fmd, "Check: F(R) * R^2 should be constant.\n\n");
    fprintf(fmd, "| R | F(R) | F(R)*R^2 |\n");
    fprintf(fmd, "|---|------|----------|\n");

    for (int vi = 0; vi < n_valid; vi++) {
        int ri = valid_R[vi];
        double R = R_list[ri];
        double dtheta = M_PI / NTHETA;
        double dphi = 2.0 * M_PI / NPHI;
        double flux = 0;
        for (int it = 0; it < NTHETA; it++) {
            double th = (it + 0.5) * dtheta;
            for (int ip = 0; ip < NPHI; ip++) {
                flux += avg_S[vi][it * NPHI + ip] * sin(th) * dtheta * dphi;
            }
        }
        fprintf(fmd, "| %.1f | %.4e | %.4e |\n", R, flux, flux * R * R);
    }

    /* Anisotropy ratio */
    fprintf(fmd, "\n## Anisotropy\n\n");
    fprintf(fmd, "| R | max(<|S|>) | min(<|S|>) | max/min |\n");
    fprintf(fmd, "|---|------------|------------|--------|\n");
    for (int vi = 0; vi < n_valid; vi++) {
        int ri = valid_R[vi];
        double R = R_list[ri];
        double vmax = 0, vmin = 1e30;
        for (int k = 0; k < NTHETA * NPHI; k++) {
            double v = avg_S[vi][k];
            if (v > vmax) vmax = v;
            if (v < vmin) vmin = v;
        }
        fprintf(fmd, "| %.1f | %.4e | %.4e | %.2f |\n", R, vmax, vmin,
                vmin > 0 ? vmax/vmin : -1);
    }

    /* Summary */
    fprintf(fmd, "\n## Interpretation\n\n");
    fprintf(fmd, "This is the **null model**: a single z-aligned braid.\n");
    fprintf(fmd, "A z-aligned magnetic dipole radiates predominantly in the l=1 mode\n");
    fprintf(fmd, "(dipolar Poynting pattern). The monopole fraction f(l=0) should be\n");
    fprintf(fmd, "small (< 0.2) for a single braid, establishing the baseline.\n\n");
    fprintf(fmd, "For the UUD three-braid composite (future analysis), the three\n");
    fprintf(fmd, "orthogonal dipoles should partially cancel, potentially producing\n");
    fprintf(fmd, "a more isotropic radiation pattern with higher f(l=0).\n");

    fclose(ftsv);
    fclose(fmd);

    printf("\n# Results written to poynting_multipole_results.md and poynting_multipole.tsv\n");

    /* Cleanup */
    for (int a = 0; a < 3; a++) {
        free(phi_cur[a]); free(phi_nxt[a]);
        free(theta_cur[a]); free(theta_nxt[a]);
    }
    free(S_mag);
    free(buf_cur); free(buf_nxt);
    for (int i = 0; i < n_valid; i++) free(avg_S[i]);
    free(avg_S);
    free(total_flux);
    sfa_close(sfa);
    return 0;
}
