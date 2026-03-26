/*  multipole.c — Spherical harmonic decomposition of theta field around UUD composite
 *
 *  Reads the last frame of a V41 UUD SFA file, finds the |P|-weighted centroid,
 *  samples theta_rms on spherical shells at radii R=5,8,10,12,15,20, and decomposes
 *  into real spherical harmonics up to l=4.
 *
 *  Key question: Does l=0 (monopole) dominate? If so, the UUD composite radiates
 *  theta like an effective point charge despite being 3 magnetic dipoles.
 *
 *  Build: gcc -O3 -march=native -fopenmp -o multipole multipole.c -lzstd -lm
 */

#define SFA_IMPLEMENTATION
#include "../../sfa/format/sfa.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ---- f16 → f64 conversion ---- */
static inline double f16_to_f64(uint16_t h) {
    uint16_t sign = h & 0x8000;
    int exp = (h >> 10) & 0x1F;
    uint16_t mant = h & 0x3FF;
    if (exp == 0) return 0.0;
    if (exp == 31) return sign ? -1e30 : 1e30;
    float fv; uint32_t x = ((uint32_t)sign << 16) | ((uint32_t)(exp-15+127) << 23) | ((uint32_t)mant << 13);
    memcpy(&fv, &x, 4); return (double)fv;
}

/* ---- Associated Legendre polynomials P_l^m(x) (unnormalized) ---- */
/* Returns P_l^m(x) for m >= 0. Uses recurrence. */
static double assoc_legendre(int l, int m, double x) {
    if (m < 0 || m > l) return 0.0;
    /* P_m^m */
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
    /* P_{m+1}^m */
    double pm1m = x * (2*m + 1) * pmm;
    if (l == m + 1) return pm1m;
    /* Recurrence */
    double pll = 0.0;
    for (int ll = m + 2; ll <= l; ll++) {
        pll = (x * (2*ll - 1) * pm1m - (ll + m - 1) * pmm) / (ll - m);
        pmm = pm1m;
        pm1m = pll;
    }
    return pll;
}

/* Factorial helper */
static double factorial(int n) {
    double f = 1.0;
    for (int i = 2; i <= n; i++) f *= i;
    return f;
}

/* Real spherical harmonic Y_l^m(theta, phi)
 * Convention: m > 0 → cos(m*phi), m < 0 → sin(|m|*phi), m = 0 → standard
 * Normalized so that ∫|Y_lm|² dΩ = 1
 */
static double real_ylm(int l, int m, double theta, double phi_ang) {
    int am = abs(m);
    /* Normalization factor */
    double norm = sqrt((2.0*l + 1.0) / (4.0 * M_PI) * factorial(l - am) / factorial(l + am));
    double plm = assoc_legendre(l, am, cos(theta));
    if (m > 0)      return norm * sqrt(2.0) * plm * cos(m * phi_ang);
    else if (m < 0) return norm * sqrt(2.0) * plm * sin(am * phi_ang);
    else            return norm * plm;
}

/* ---- Trilinear interpolation ---- */
static double interp3d(const double *field, int N, double L, double dx,
                       double x, double y, double z) {
    /* Map to grid coords: grid runs from 0 to N-1, coord from -L to +L */
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

#define LMAX 4
#define NTHETA 72
#define NPHI 144

int main(int argc, char **argv) {
    if (argc < 2) { fprintf(stderr, "Usage: %s input.sfa [frame_idx]\n", argv[0]); return 1; }

    int target_frame = -1;
    if (argc > 2) target_frame = atoi(argv[2]);

    SFA *sfa = sfa_open(argv[1]);
    if (!sfa) { fprintf(stderr, "Cannot open %s\n", argv[1]); return 1; }

    int N = sfa->Nx; double L = sfa->Lx;
    long N3 = (long)N*N*N; int NN = N*N;
    double dx = 2.0*L / (N - 1);

    printf("# Multipole analysis: %s (N=%d L=%.1f %u frames)\n", argv[1], N, L, sfa->total_frames);

    if (target_frame < 0) target_frame = sfa->total_frames - 1;
    if (target_frame >= (int)sfa->total_frames) target_frame = sfa->total_frames - 1;

    printf("# Reading frame %d ...\n", target_frame);

    void *buf = malloc(sfa->frame_bytes);
    if (!buf) { fprintf(stderr, "Out of memory (%lu bytes)\n", sfa->frame_bytes); return 1; }
    sfa_read_frame(sfa, target_frame, buf);
    double t = sfa_frame_time(sfa, target_frame);

    /* Extract phi and theta fields */
    double *phi[3], *theta[3];
    for (int a = 0; a < 3; a++) {
        phi[a] = (double *)calloc(N3, sizeof(double));
        theta[a] = (double *)calloc(N3, sizeof(double));
    }

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
            if (dtype == SFA_F64) for (long i = 0; i < N3; i++) target[i] = ((double*)src)[i];
            else if (dtype == SFA_F32) for (long i = 0; i < N3; i++) target[i] = (double)((float*)src)[i];
            else if (dtype == SFA_F16) for (long i = 0; i < N3; i++) target[i] = f16_to_f64(((uint16_t*)src)[i]);
        }
        off += (uint64_t)N3 * es;
    }
    free(buf);

    /* Find centroid weighted by |P| = |phi_x * phi_y * phi_z| */
    double cm[3] = {0}, wt = 0;
    for (long idx = 0; idx < N3; idx++) {
        double P = fabs(phi[0][idx] * phi[1][idx] * phi[2][idx]);
        int i = (int)(idx / NN), j = (int)((idx / N) % N), k = (int)(idx % N);
        double x = -L + i * dx, y = -L + j * dx, z = -L + k * dx;
        cm[0] += x * P; cm[1] += y * P; cm[2] += z * P; wt += P;
    }
    if (wt > 0) { cm[0] /= wt; cm[1] /= wt; cm[2] /= wt; }
    printf("# Frame %d, t=%.2f, centroid=(%.3f, %.3f, %.3f), weight=%.4f\n\n",
           target_frame, t, cm[0], cm[1], cm[2], wt);

    /* Also compute theta_rms field for trilinear interpolation */
    double *theta_rms_field = (double *)calloc(N3, sizeof(double));
    for (long idx = 0; idx < N3; idx++) {
        theta_rms_field[idx] = sqrt(theta[0][idx]*theta[0][idx]
                                  + theta[1][idx]*theta[1][idx]
                                  + theta[2][idx]*theta[2][idx]);
    }

    /* Also compute per-component fields for vector decomposition */
    /* theta_r (radial component) for checking if field is radial or tangential */

    /* Radii to analyze */
    double R_list[] = {3.0, 5.0, 8.0, 10.0, 12.0, 15.0};
    int n_R = sizeof(R_list) / sizeof(R_list[0]);

    /* Open TSV output */
    FILE *ftsv = fopen("multipole.tsv", "w");
    fprintf(ftsv, "R\tl\tC_l\tfraction\tsum_check\n");

    /* Open markdown output */
    FILE *fmd = fopen("multipole_results.md", "w");
    fprintf(fmd, "# Multipole Expansion of Theta Field — UUD Composite (V41)\n\n");
    fprintf(fmd, "**Source**: `%s` frame %d (t=%.2f)\n", argv[1], target_frame, t);
    fprintf(fmd, "**Grid**: N=%d, L=%.1f, dx=%.4f\n", N, L, dx);
    fprintf(fmd, "**Centroid**: (%.3f, %.3f, %.3f)\n\n", cm[0], cm[1], cm[2]);
    fprintf(fmd, "## Method\n\n");
    fprintf(fmd, "- Sample theta_rms = sqrt(theta_x^2 + theta_y^2 + theta_z^2) on spherical shells\n");
    fprintf(fmd, "- Decompose into real spherical harmonics Y_l^m up to l=%d\n", LMAX);
    fprintf(fmd, "- Integration grid: %d x %d (theta x phi) points\n", NTHETA, NPHI);
    fprintf(fmd, "- Trilinear interpolation from 3D grid to shell points\n\n");

    printf("# Multipole power spectrum: C_l = sum_m |a_lm|^2\n");
    printf("# R        l=0       l=1       l=2       l=3       l=4       f(l=0)    f(l=1)    f(l=2)\n");

    fprintf(fmd, "## Results\n\n");
    fprintf(fmd, "| R | l=0 frac | l=1 frac | l=2 frac | l=3 frac | l=4 frac | dominant |\n");
    fprintf(fmd, "|---|---------|---------|---------|---------|---------|----------|\n");

    for (int ri = 0; ri < n_R; ri++) {
        double R = R_list[ri];

        /* Check if R is within domain */
        if (R > L - 2*dx) {
            printf("# R=%.1f too large for domain L=%.1f, skipping\n", R, L);
            continue;
        }

        /* Compute a_lm coefficients by numerical integration over (theta, phi) */
        double dtheta = M_PI / NTHETA;
        double dphi = 2.0 * M_PI / NPHI;

        /* a_lm storage: l from 0 to LMAX, m from -l to +l */
        double a_lm[LMAX+1][2*LMAX+1];
        memset(a_lm, 0, sizeof(a_lm));

        double total_signal = 0;
        int n_nonzero = 0;

        for (int it = 0; it < NTHETA; it++) {
            double th = (it + 0.5) * dtheta;  /* midpoint rule */
            double sth = sin(th);

            for (int ip = 0; ip < NPHI; ip++) {
                double ph = (ip + 0.5) * dphi;

                /* 3D coordinates on shell */
                double xs = cm[0] + R * sth * cos(ph);
                double ys = cm[1] + R * sth * sin(ph);
                double zs = cm[2] + R * cos(th);

                /* Interpolate theta_rms */
                double val = interp3d(theta_rms_field, N, L, dx, xs, ys, zs);
                if (val > 0) n_nonzero++;
                total_signal += val * val * sth * dtheta * dphi;

                /* Accumulate a_lm */
                for (int l = 0; l <= LMAX; l++) {
                    for (int m = -l; m <= l; m++) {
                        double ylm = real_ylm(l, m, th, ph);
                        a_lm[l][m + LMAX] += val * ylm * sth * dtheta * dphi;
                    }
                }
            }
        }

        /* Compute power per multipole C_l = sum_m |a_lm|^2 */
        double C[LMAX+1];
        double C_total = 0;
        for (int l = 0; l <= LMAX; l++) {
            C[l] = 0;
            for (int m = -l; m <= l; m++) {
                C[l] += a_lm[l][m + LMAX] * a_lm[l][m + LMAX];
            }
            C_total += C[l];
        }

        /* Print results */
        printf("  %4.1f  ", R);
        for (int l = 0; l <= LMAX; l++) printf("%9.4e ", C[l]);
        for (int l = 0; l <= 2; l++) printf("%8.4f ", C_total > 0 ? C[l]/C_total : 0);
        printf("  [%d/%d nonzero]\n", n_nonzero, NTHETA*NPHI);

        /* Determine dominant mode */
        int dom = 0;
        for (int l = 1; l <= LMAX; l++) if (C[l] > C[dom]) dom = l;
        const char *names[] = {"monopole", "dipole", "quadrupole", "octupole", "hexadecapole"};

        /* TSV */
        for (int l = 0; l <= LMAX; l++) {
            fprintf(ftsv, "%.1f\t%d\t%.6e\t%.6f\t%.6e\n",
                    R, l, C[l], C_total > 0 ? C[l]/C_total : 0, C_total);
        }

        /* Markdown */
        fprintf(fmd, "| %.1f ", R);
        for (int l = 0; l <= LMAX; l++)
            fprintf(fmd, "| %.1f%% ", C_total > 0 ? 100.0*C[l]/C_total : 0);
        fprintf(fmd, "| %s |\n", names[dom]);

        /* Detailed per-m breakdown */
        printf("#   a_lm detail for R=%.1f:\n", R);
        for (int l = 0; l <= LMAX; l++) {
            printf("#     l=%d:", l);
            for (int m = -l; m <= l; m++) {
                printf(" m=%+d:%.3e", m, a_lm[l][m + LMAX]);
            }
            printf("\n");
        }
    }

    /* Also do per-component analysis at R=10 */
    fprintf(fmd, "\n## Per-Component Analysis at R=10\n\n");
    fprintf(fmd, "Decompose each theta component separately:\n\n");

    double R_detail = 10.0;
    if (R_detail <= L - 2*dx) {
        double dtheta = M_PI / NTHETA;
        double dphi = 2.0 * M_PI / NPHI;

        const char *comp_names[] = {"theta_x", "theta_y", "theta_z"};
        fprintf(fmd, "| component | l=0 frac | l=1 frac | l=2 frac | l=3 frac | l=4 frac |\n");
        fprintf(fmd, "|-----------|---------|---------|---------|---------|----------|\n");

        for (int comp = 0; comp < 3; comp++) {
            double a_lm[LMAX+1][2*LMAX+1];
            memset(a_lm, 0, sizeof(a_lm));

            for (int it = 0; it < NTHETA; it++) {
                double th = (it + 0.5) * dtheta;
                double sth = sin(th);
                for (int ip = 0; ip < NPHI; ip++) {
                    double ph = (ip + 0.5) * dphi;
                    double xs = cm[0] + R_detail * sth * cos(ph);
                    double ys = cm[1] + R_detail * sth * sin(ph);
                    double zs = cm[2] + R_detail * cos(th);

                    double val = interp3d(theta[comp], N, L, dx, xs, ys, zs);

                    for (int l = 0; l <= LMAX; l++)
                        for (int m = -l; m <= l; m++)
                            a_lm[l][m + LMAX] += val * real_ylm(l, m, th, ph) * sth * dtheta * dphi;
                }
            }

            double C[LMAX+1], C_total = 0;
            for (int l = 0; l <= LMAX; l++) {
                C[l] = 0;
                for (int m = -l; m <= l; m++)
                    C[l] += a_lm[l][m + LMAX] * a_lm[l][m + LMAX];
                C_total += C[l];
            }

            fprintf(fmd, "| %s ", comp_names[comp]);
            for (int l = 0; l <= LMAX; l++)
                fprintf(fmd, "| %.1f%% ", C_total > 0 ? 100.0*C[l]/C_total : 0);
            fprintf(fmd, "|\n");

            printf("# Per-component R=%.0f %s: ", R_detail, comp_names[comp]);
            for (int l = 0; l <= LMAX; l++)
                printf("l=%d:%.1f%% ", l, C_total > 0 ? 100.0*C[l]/C_total : 0);
            printf("\n");
        }
    }

    /* Summary */
    fprintf(fmd, "\n## Interpretation\n\n");
    fprintf(fmd, "If l=0 dominates at large R, the UUD composite acts as an effective\n");
    fprintf(fmd, "point source for the theta field, despite being composed of three\n");
    fprintf(fmd, "magnetic dipoles. This would be the mechanism by which Coulomb's law\n");
    fprintf(fmd, "emerges from the composite structure.\n\n");
    fprintf(fmd, "If l=1 dominates, the composite still looks like a dipole from afar.\n");
    fprintf(fmd, "If l=2 dominates, the composite has quadrupole character.\n");

    fclose(ftsv);
    fclose(fmd);

    printf("\n# Results written to multipole_results.md and multipole.tsv\n");

    /* Cleanup */
    for (int a = 0; a < 3; a++) { free(phi[a]); free(theta[a]); }
    free(theta_rms_field);
    sfa_close(sfa);
    return 0;
}
