/*  multipole_analysis.c — Angular power spectrum of theta field around baryon
 *
 *  Computes spherical harmonic decomposition of theta_rms on shells at
 *  various radii from the baryon centroid. Tests whether a UUD composite
 *  acts as an effective electric monopole (l=0 dominant) despite being
 *  made of magnetic dipoles.
 *
 *  Method: On each shell at radius R, sample theta_rms at N_samp points
 *  uniformly distributed on the sphere (Fibonacci lattice). Decompose
 *  into real spherical harmonics up to l_max using least-squares projection.
 *  Report power per l: C_l = sum_m |a_lm|^2 / (2l+1).
 *
 *  Build: gcc -O3 -fopenmp -o multipole_analysis multipole_analysis.c -lzstd -lm
 *  Usage: ./multipole_analysis input.sfa [frame_idx]
 */

#define SFA_IMPLEMENTATION
#include "../sfa/format/sfa.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static const double MU = -41.345, KAPPA = 50.0, MASS2 = 2.25;
static const double GOLDEN_RATIO = 1.6180339887498949;

static inline double f16_to_f64(uint16_t h) {
    uint16_t sign = h & 0x8000;
    int exp = (h >> 10) & 0x1F;
    uint16_t mant = h & 0x3FF;
    if (exp == 0) return 0.0;
    if (exp == 31) return sign ? -1e30 : 1e30;
    float fv; uint32_t x = ((uint32_t)sign << 16) | ((uint32_t)(exp-15+127) << 23) | ((uint32_t)mant << 13);
    memcpy(&fv, &x, 4); return (double)fv;
}

/* Associated Legendre polynomials P_l^m(x) (unnormalized, Condon-Shortley) */
static double assoc_legendre(int l, int m, double x) {
    if (m < 0) {
        m = -m;
        double sign = (m % 2 == 0) ? 1.0 : -1.0;
        double fac = 1.0;
        for (int i = l-m+1; i <= l+m; i++) fac *= i;
        return sign * assoc_legendre(l, m, x) / fac;
    }
    /* P_m^m */
    double pmm = 1.0;
    if (m > 0) {
        double sq = sqrt(1.0 - x*x);
        double fact = 1.0;
        for (int i = 1; i <= m; i++) {
            pmm *= -fact * sq;
            fact += 2.0;
        }
    }
    if (l == m) return pmm;
    /* P_{m+1}^m */
    double pmmp1 = x * (2*m + 1) * pmm;
    if (l == m+1) return pmmp1;
    /* Recurrence */
    double pll = 0.0;
    for (int ll = m+2; ll <= l; ll++) {
        pll = (x * (2*ll - 1) * pmmp1 - (ll + m - 1) * pmm) / (ll - m);
        pmm = pmmp1;
        pmmp1 = pll;
    }
    return pll;
}

/* Real spherical harmonic Y_lm(theta, phi) — fully normalized */
static double real_sph_harm(int l, int m, double theta, double phi) {
    int am = abs(m);
    /* Normalization factor */
    double num = (2*l + 1);
    double fac = 1.0;
    for (int i = l - am + 1; i <= l + am; i++) fac *= i;
    double norm = sqrt(num / (4.0 * M_PI * fac));
    double plm = assoc_legendre(l, am, cos(theta));
    if (m > 0) return norm * plm * sqrt(2.0) * cos(m * phi);
    if (m < 0) return norm * plm * sqrt(2.0) * sin(am * phi);
    return norm * plm;  /* m == 0 */
}

/* Trilinear interpolation */
static double interp3d(const double *field, int N, double fi, double fj, double fk) {
    int i0 = (int)floor(fi), j0 = (int)floor(fj), k0 = (int)floor(fk);
    if (i0 < 0 || i0 >= N-1 || j0 < 0 || j0 >= N-1 || k0 < 0 || k0 >= N-1)
        return 0.0;
    double di = fi - i0, dj = fj - j0, dk = fk - k0;
    int NN = N * N;
    double c000 = field[i0*NN + j0*N + k0];
    double c001 = field[i0*NN + j0*N + k0+1];
    double c010 = field[i0*NN + (j0+1)*N + k0];
    double c011 = field[i0*NN + (j0+1)*N + k0+1];
    double c100 = field[(i0+1)*NN + j0*N + k0];
    double c101 = field[(i0+1)*NN + j0*N + k0+1];
    double c110 = field[(i0+1)*NN + (j0+1)*N + k0];
    double c111 = field[(i0+1)*NN + (j0+1)*N + k0+1];
    return c000*(1-di)*(1-dj)*(1-dk) + c001*(1-di)*(1-dj)*dk
         + c010*(1-di)*dj*(1-dk)     + c011*(1-di)*dj*dk
         + c100*di*(1-dj)*(1-dk)     + c101*di*(1-dj)*dk
         + c110*di*dj*(1-dk)         + c111*di*dj*dk;
}

int main(int argc, char **argv) {
    if (argc < 2) {
        fprintf(stderr, "Usage: %s input.sfa [frame_idx]\n", argv[0]);
        return 1;
    }
    int target_frame = -1;
    if (argc > 2) target_frame = atoi(argv[2]);

    SFA *sfa = sfa_open(argv[1]);
    if (!sfa) { fprintf(stderr, "Cannot open %s\n", argv[1]); return 1; }
    int N = sfa->Nx; double L = sfa->Lx;
    long N3 = (long)N*N*N; int NN = N*N;
    double dx = 2.0*L/(N-1);

    printf("# Multipole analysis: %s (N=%d L=%.1f %u frames)\n", argv[1], N, L, sfa->total_frames);

    if (target_frame < 0) target_frame = sfa->total_frames - 1;
    if (target_frame >= (int)sfa->total_frames) target_frame = sfa->total_frames - 1;

    void *buf = malloc(sfa->frame_bytes);
    sfa_read_frame(sfa, target_frame, buf);
    double t = sfa_frame_time(sfa, target_frame);

    /* Extract phi and theta fields */
    double *phi[3], *theta[3];
    for (int a = 0; a < 3; a++) {
        phi[a] = (double*)calloc(N3, sizeof(double));
        theta[a] = (double*)calloc(N3, sizeof(double));
    }

    uint64_t off = 0;
    for (int c = 0; c < sfa->n_columns; c++) {
        int dtype = sfa->columns[c].dtype;
        int sem = sfa->columns[c].semantic;
        int comp = sfa->columns[c].component;
        int es = sfa_dtype_size[dtype];
        uint8_t *src = (uint8_t*)buf + off;
        double *target = NULL;
        if (sem == SFA_POSITION && comp < 3) target = phi[comp];
        else if (sem == SFA_ANGLE && comp < 3) target = theta[comp];
        if (target) {
            if (dtype == SFA_F64) for(long i=0;i<N3;i++) target[i] = ((double*)src)[i];
            else if (dtype == SFA_F32) for(long i=0;i<N3;i++) target[i] = (double)((float*)src)[i];
            else if (dtype == SFA_F16) for(long i=0;i<N3;i++) target[i] = f16_to_f64(((uint16_t*)src)[i]);
        }
        off += (uint64_t)N3 * es;
    }
    free(buf);

    /* Find centroid weighted by |P| */
    double cm[3] = {0}, wt = 0;
    for (long idx = 0; idx < N3; idx++) {
        double P = fabs(phi[0][idx] * phi[1][idx] * phi[2][idx]);
        int ii = idx / NN, jj = (idx / N) % N, kk = idx % N;
        double x = -L + ii * dx, y = -L + jj * dx, z = -L + kk * dx;
        cm[0] += P * x; cm[1] += P * y; cm[2] += P * z;
        wt += P;
    }
    cm[0] /= wt; cm[1] /= wt; cm[2] /= wt;
    printf("# Centroid: (%.3f, %.3f, %.3f)  frame=%d  t=%.2f\n", cm[0], cm[1], cm[2], target_frame, t);

    /* Compute theta_rms field (magnitude at each voxel) */
    double *theta_mag = (double*)calloc(N3, sizeof(double));
    for (long i = 0; i < N3; i++)
        theta_mag[i] = sqrt(theta[0][i]*theta[0][i] + theta[1][i]*theta[1][i] + theta[2][i]*theta[2][i]);

    /* Also compute each theta component field for vector decomposition */
    /* We'll decompose theta_mag (scalar) AND each theta component */

    /* Multipole analysis parameters */
    int L_MAX = 6;       /* up to hexadecapole */
    int N_SAMP = 2000;   /* Fibonacci sphere samples */
    double radii[] = {2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 15.0, 20.0};
    int n_radii = sizeof(radii) / sizeof(radii[0]);

    int n_coeffs = (L_MAX + 1) * (L_MAX + 1);  /* total (l,m) pairs */

    printf("\n# Multipole power spectrum C_l = sum_m |a_lm|^2 / (2l+1)\n");
    printf("# Each row: radius, then C_0 through C_%d, then fractional powers\n", L_MAX);
    printf("#\n");

    /* Header */
    printf("# %6s", "r");
    for (int l = 0; l <= L_MAX; l++) printf("  %10s%d", "C_", l);
    printf("  ");
    for (int l = 0; l <= L_MAX; l++) printf("  %9s%d", "frac_", l);
    printf("  %10s\n", "theta_rms");

    for (int ri = 0; ri < n_radii; ri++) {
        double R = radii[ri];
        if (R > L * 0.8) continue;  /* skip if too close to boundary */

        /* Generate Fibonacci sphere sample points */
        double *sph_theta = (double*)malloc(N_SAMP * sizeof(double));
        double *sph_phi = (double*)malloc(N_SAMP * sizeof(double));
        double *f_vals = (double*)malloc(N_SAMP * sizeof(double));

        for (int s = 0; s < N_SAMP; s++) {
            /* Fibonacci lattice on sphere */
            double y_fib = 1.0 - (2.0 * s + 1.0) / N_SAMP;
            sph_theta[s] = acos(y_fib);
            sph_phi[s] = fmod(2.0 * M_PI * s / GOLDEN_RATIO, 2.0 * M_PI);

            /* World coordinates */
            double st = sin(sph_theta[s]), ct = cos(sph_theta[s]);
            double sp = sin(sph_phi[s]), cp = cos(sph_phi[s]);
            double wx = cm[0] + R * st * cp;
            double wy = cm[1] + R * st * sp;
            double wz = cm[2] + R * ct;

            /* Grid fractional indices */
            double fi = (wx + L) / dx;
            double fj = (wy + L) / dx;
            double fk = (wz + L) / dx;

            f_vals[s] = interp3d(theta_mag, N, fi, fj, fk);
        }

        /* Project onto spherical harmonics: a_lm = (4pi/N) sum_s f(s) * Y_lm(s) */
        double *a_lm = (double*)calloc(n_coeffs, sizeof(double));
        double total_power = 0;
        double mean_theta = 0;

        for (int s = 0; s < N_SAMP; s++)
            mean_theta += f_vals[s];
        mean_theta /= N_SAMP;

        for (int l = 0; l <= L_MAX; l++) {
            for (int m = -l; m <= l; m++) {
                double sum = 0;
                for (int s = 0; s < N_SAMP; s++) {
                    double ylm = real_sph_harm(l, m, sph_theta[s], sph_phi[s]);
                    sum += f_vals[s] * ylm;
                }
                int idx = l * l + l + m;
                a_lm[idx] = sum * 4.0 * M_PI / N_SAMP;
            }
        }

        /* Compute power per l */
        double C_l[L_MAX + 1];
        double total = 0;
        for (int l = 0; l <= L_MAX; l++) {
            C_l[l] = 0;
            for (int m = -l; m <= l; m++) {
                int idx = l * l + l + m;
                C_l[l] += a_lm[idx] * a_lm[idx];
            }
            C_l[l] /= (2*l + 1);
            total += C_l[l] * (2*l + 1);
        }

        /* Print */
        printf("  %6.1f", R);
        for (int l = 0; l <= L_MAX; l++) printf("  %11.6f", C_l[l]);
        printf("  ");
        for (int l = 0; l <= L_MAX; l++) {
            double frac = total > 0 ? C_l[l] * (2*l + 1) / total : 0;
            printf("  %9.4f", frac);
        }
        printf("  %10.6f\n", mean_theta);

        free(sph_theta); free(sph_phi); free(f_vals); free(a_lm);
    }

    /* Also decompose each theta COMPONENT separately at r=10 and r=20 */
    printf("\n# Per-component decomposition at r=10 and r=20\n");
    printf("# component  radius  C_0       C_1       C_2       C_3       frac_0    frac_1    frac_2\n");

    double comp_radii[] = {10.0, 20.0};
    const char *comp_names[] = {"theta_x", "theta_y", "theta_z"};

    for (int ci = 0; ci < 2; ci++) {
        double R = comp_radii[ci];
        if (R > L * 0.8) continue;

        double *sph_theta_s = (double*)malloc(N_SAMP * sizeof(double));
        double *sph_phi_s = (double*)malloc(N_SAMP * sizeof(double));

        for (int s = 0; s < N_SAMP; s++) {
            double y_fib = 1.0 - (2.0 * s + 1.0) / N_SAMP;
            sph_theta_s[s] = acos(y_fib);
            sph_phi_s[s] = fmod(2.0 * M_PI * s / GOLDEN_RATIO, 2.0 * M_PI);
        }

        for (int comp = 0; comp < 3; comp++) {
            double *f_c = (double*)malloc(N_SAMP * sizeof(double));
            for (int s = 0; s < N_SAMP; s++) {
                double st = sin(sph_theta_s[s]), ct = cos(sph_theta_s[s]);
                double sp = sin(sph_phi_s[s]), cp = cos(sph_phi_s[s]);
                double wx = cm[0] + R * st * cp;
                double wy = cm[1] + R * st * sp;
                double wz = cm[2] + R * ct;
                double fi = (wx + L) / dx, fj = (wy + L) / dx, fk = (wz + L) / dx;
                f_c[s] = interp3d(theta[comp], N, fi, fj, fk);
            }

            /* Project */
            int lmax_c = 3;
            double C_l_c[4] = {0};
            double total_c = 0;
            for (int l = 0; l <= lmax_c; l++) {
                for (int m = -l; m <= l; m++) {
                    double sum = 0;
                    for (int s = 0; s < N_SAMP; s++) {
                        double ylm = real_sph_harm(l, m, sph_theta_s[s], sph_phi_s[s]);
                        sum += f_c[s] * ylm;
                    }
                    double alm = sum * 4.0 * M_PI / N_SAMP;
                    C_l_c[l] += alm * alm;
                }
                C_l_c[l] /= (2*l + 1);
                total_c += C_l_c[l] * (2*l + 1);
            }

            printf("  %-8s  %5.1f", comp_names[comp], R);
            for (int l = 0; l <= lmax_c; l++) printf("  %9.6f", C_l_c[l]);
            for (int l = 0; l <= 2; l++) {
                double frac = total_c > 0 ? C_l_c[l] * (2*l + 1) / total_c : 0;
                printf("  %9.4f", frac);
            }
            printf("\n");
            free(f_c);
        }
        free(sph_theta_s); free(sph_phi_s);
    }

    /* Cleanup */
    for (int a = 0; a < 3; a++) { free(phi[a]); free(theta[a]); }
    free(theta_mag);
    sfa_close(sfa);
    return 0;
}
