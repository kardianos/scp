/*  virial_soliton.c — Compute per-voxel virial integrals over soliton region
 *
 *  Reads an SFA file, extracts φ and θ fields, computes:
 *  - P = φ₀φ₁φ₂ (triple product)
 *  - ∇φ, ∇×φ, ∇×θ (finite differences)
 *  - Energy densities: gradient, mass, potential, coupling
 *  - Virial integrals restricted to soliton region (|P| > threshold)
 *  - η₁ from the self-consistency equation
 *
 *  Build: gcc -O3 -fopenmp -o virial_soliton virial_soliton.c -lzstd -lm
 *  Usage: ./virial_soliton input.sfa [P_threshold]
 */

#define SFA_IMPLEMENTATION
#include "../../sfa/format/sfa.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <sys/stat.h>

static const double MU = -41.345, KAPPA = 50.0, MASS2 = 2.25, ETA0 = 0.5;

static inline double f16_to_f64(uint16_t h) {
    uint16_t sign = h & 0x8000;
    int exp = (h >> 10) & 0x1F;
    uint16_t mant = h & 0x3FF;
    if (exp == 0) return 0.0;
    if (exp == 31) return sign ? -1e30 : 1e30;
    float fv; uint32_t x = ((uint32_t)sign << 16) | ((uint32_t)(exp-15+127) << 23) | ((uint32_t)mant << 13);
    memcpy(&fv, &x, 4); return (double)fv;
}

/* V(P) and V'(P) */
static inline double V_of_P(double P) {
    return (MU/2.0) * P*P / (1.0 + KAPPA*P*P);
}

int main(int argc, char **argv) {
    if (argc < 2) {
        fprintf(stderr, "Usage: %s input.sfa [P_threshold]\n", argv[0]);
        return 1;
    }
    double P_thresh = 0.001;  /* default: include voxels with |P| > 0.001 */
    if (argc > 2) P_thresh = atof(argv[2]);

    SFA *sfa = sfa_open(argv[1]);
    if (!sfa) { fprintf(stderr, "Cannot open %s\n", argv[1]); return 1; }
    int N = sfa->Nx; double L = sfa->Lx;
    long N3 = (long)N*N*N; int NN = N*N;
    double dx = 2.0*L/(N-1);
    double dV = dx*dx*dx;
    double idx1 = 0.5/dx;  /* 1/(2dx) for central differences */

    /* Read last frame */
    int frame = sfa->total_frames - 1;
    if (frame < 0) {
        /* Try fixup for streaming files */
        struct stat st;
        fstat(fileno(sfa->fp), &st);
        int valid = sfa_count_valid_frames(sfa, st.st_size);
        frame = valid - 1;
    }
    if (frame < 0) { fprintf(stderr, "No valid frames\n"); return 1; }

    void *buf = malloc(sfa->frame_bytes);
    sfa_read_frame(sfa, frame, buf);
    double t = sfa_frame_time(sfa, frame);

    printf("# Virial soliton analysis: %s\n", argv[1]);
    printf("# N=%d L=%.1f dx=%.4f frame=%d t=%.2f P_thresh=%.4f\n",
           N, L, dx, frame, t, P_thresh);

    /* Extract φ and θ fields */
    double *phi[3], *theta[3];
    for (int a = 0; a < 3; a++) {
        phi[a] = (double*)calloc(N3, sizeof(double));
        theta[a] = (double*)calloc(N3, sizeof(double));
    }

    uint64_t off = 0;
    for (int c = 0; c < sfa->n_columns && c < 6; c++) {
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

    /* Compute virial integrals over soliton region */
    double I_grad_sol = 0;      /* ∫|∇φ|² over soliton */
    double I_tgrad_sol = 0;     /* ∫|∇θ|² over soliton */
    double I_mass_sol = 0;      /* ∫m²|φ|² over soliton */
    double I_pot_sol = 0;       /* ∫V(P) over soliton */
    double I_C0_sol = 0;        /* ∫φ·(∇×θ) over soliton */
    double I_C1_sol = 0;        /* ∫|P|·φ·(∇×θ) over soliton */
    double I_P2_sol = 0;        /* ∫P² over soliton */
    long n_soliton = 0;         /* number of soliton voxels */

    /* Also compute background estimates for comparison */
    double I_grad_bg = 0, I_mass_bg = 0, I_C0_bg = 0;
    long n_bg = 0;

    #pragma omp parallel for schedule(dynamic,4) reduction(+:I_grad_sol,I_tgrad_sol,I_mass_sol,I_pot_sol,I_C0_sol,I_C1_sol,I_P2_sol,n_soliton,I_grad_bg,I_mass_bg,I_C0_bg,n_bg)
    for (int ii = 1; ii < N-1; ii++) {
        for (int jj = 1; jj < N-1; jj++) {
            for (int kk = 1; kk < N-1; kk++) {
                long idx = (long)ii*NN + jj*N + kk;

                /* Triple product */
                double P = phi[0][idx] * phi[1][idx] * phi[2][idx];
                double absP = fabs(P);

                /* Gradients (central differences) */
                double grad_phi2 = 0, grad_theta2 = 0;
                double curl_phi[3], curl_theta[3];

                for (int a = 0; a < 3; a++) {
                    double dpx = (phi[a][(ii+1)*NN+jj*N+kk] - phi[a][(ii-1)*NN+jj*N+kk]) * idx1;
                    double dpy = (phi[a][ii*NN+(jj+1)*N+kk] - phi[a][ii*NN+(jj-1)*N+kk]) * idx1;
                    double dpz = (phi[a][ii*NN+jj*N+kk+1] - phi[a][ii*NN+jj*N+kk-1]) * idx1;
                    grad_phi2 += dpx*dpx + dpy*dpy + dpz*dpz;

                    double dtx = (theta[a][(ii+1)*NN+jj*N+kk] - theta[a][(ii-1)*NN+jj*N+kk]) * idx1;
                    double dty = (theta[a][ii*NN+(jj+1)*N+kk] - theta[a][ii*NN+(jj-1)*N+kk]) * idx1;
                    double dtz = (theta[a][ii*NN+jj*N+kk+1] - theta[a][ii*NN+jj*N+kk-1]) * idx1;
                    grad_theta2 += dtx*dtx + dty*dty + dtz*dtz;

                    /* Store for curl computation */
                    if (a == 0) {
                        /* Will be used for curl components */
                    }
                }

                /* Curl of φ: (∇×φ)_a = ε_abc ∂φ_c/∂x_b */
                /* (∇×φ)_0 = ∂φ_2/∂y - ∂φ_1/∂z */
                curl_phi[0] = (phi[2][ii*NN+(jj+1)*N+kk] - phi[2][ii*NN+(jj-1)*N+kk]) * idx1
                            - (phi[1][ii*NN+jj*N+kk+1] - phi[1][ii*NN+jj*N+kk-1]) * idx1;
                /* (∇×φ)_1 = ∂φ_0/∂z - ∂φ_2/∂x */
                curl_phi[1] = (phi[0][ii*NN+jj*N+kk+1] - phi[0][ii*NN+jj*N+kk-1]) * idx1
                            - (phi[2][(ii+1)*NN+jj*N+kk] - phi[2][(ii-1)*NN+jj*N+kk]) * idx1;
                /* (∇×φ)_2 = ∂φ_1/∂x - ∂φ_0/∂y */
                curl_phi[2] = (phi[1][(ii+1)*NN+jj*N+kk] - phi[1][(ii-1)*NN+jj*N+kk]) * idx1
                            - (phi[0][ii*NN+(jj+1)*N+kk] - phi[0][ii*NN+(jj-1)*N+kk]) * idx1;

                /* Curl of θ */
                curl_theta[0] = (theta[2][ii*NN+(jj+1)*N+kk] - theta[2][ii*NN+(jj-1)*N+kk]) * idx1
                              - (theta[1][ii*NN+jj*N+kk+1] - theta[1][ii*NN+jj*N+kk-1]) * idx1;
                curl_theta[1] = (theta[0][ii*NN+jj*N+kk+1] - theta[0][ii*NN+jj*N+kk-1]) * idx1
                              - (theta[2][(ii+1)*NN+jj*N+kk] - theta[2][(ii-1)*NN+jj*N+kk]) * idx1;
                curl_theta[2] = (theta[1][(ii+1)*NN+jj*N+kk] - theta[1][(ii-1)*NN+jj*N+kk]) * idx1
                              - (theta[0][ii*NN+(jj+1)*N+kk] - theta[0][ii*NN+(jj-1)*N+kk]) * idx1;

                /* φ · (∇×θ) */
                double phi_dot_curl_theta = 0;
                for (int a = 0; a < 3; a++)
                    phi_dot_curl_theta += phi[a][idx] * curl_theta[a];

                /* Mass energy density */
                double phi2 = 0;
                for (int a = 0; a < 3; a++)
                    phi2 += phi[a][idx] * phi[a][idx];
                double e_mass = MASS2 * phi2;

                /* Classify: soliton or background */
                if (absP > P_thresh) {
                    /* SOLITON region */
                    I_grad_sol += 0.5 * grad_phi2 * dV;
                    I_tgrad_sol += 0.5 * grad_theta2 * dV;
                    I_mass_sol += 0.5 * e_mass * dV;
                    I_pot_sol += V_of_P(P) * dV;
                    I_C0_sol += phi_dot_curl_theta * dV;
                    I_C1_sol += absP * phi_dot_curl_theta * dV;
                    I_P2_sol += P*P * dV;
                    n_soliton++;
                } else {
                    /* BACKGROUND region */
                    I_grad_bg += 0.5 * grad_phi2 * dV;
                    I_mass_bg += 0.5 * e_mass * dV;
                    I_C0_bg += phi_dot_curl_theta * dV;
                    n_bg++;
                }
            }
        }
    }

    printf("\n# Soliton region (|P| > %.4f): %ld voxels (%.1f%%)\n",
           P_thresh, n_soliton, 100.0*n_soliton/(N3));
    printf("# Background region: %ld voxels\n\n", n_bg);

    printf("# SOLITON-ONLY ENERGY INTEGRALS:\n");
    printf("#   I_grad   = %.6f  (½∫|∇φ|² over soliton)\n", I_grad_sol);
    printf("#   I_tgrad  = %.6f  (½∫|∇θ|² over soliton)\n", I_tgrad_sol);
    printf("#   I_mass   = %.6f  (½∫m²φ² over soliton)\n", I_mass_sol);
    printf("#   I_pot    = %.6f  (∫V(P) over soliton)\n", I_pot_sol);
    printf("#   I_C0     = %.6f  (∫φ·(∇×θ) over soliton)\n", I_C0_sol);
    printf("#   I_C1     = %.6f  (∫|P|·φ·(∇×θ) over soliton)\n", I_C1_sol);
    printf("#   I_P2     = %.6f  (∫P² over soliton)\n", I_P2_sol);

    printf("\n# BACKGROUND INTEGRALS (for comparison):\n");
    printf("#   I_grad_bg = %.1f\n", I_grad_bg);
    printf("#   I_mass_bg = %.1f\n", I_mass_bg);
    printf("#   I_C0_bg   = %.6f\n", I_C0_bg);

    /* Virial identity: E₁ + 2E₂ + 3E₃ = 0
       E₁ = I_grad + I_tgrad  (λ¹)
       E₂ = η₀ I_C0           (λ²)
       E₃ = I_mass + I_pot    (λ³)
    */
    double E1 = I_grad_sol + I_tgrad_sol;
    double E2 = ETA0 * I_C0_sol;
    double E3 = I_mass_sol + I_pot_sol;
    double R = E1 + 2*E2 + 3*E3;

    printf("\n# VIRIAL ANALYSIS (soliton only):\n");
    printf("#   E₁ (λ¹) = I_grad + I_tgrad = %.6f\n", E1);
    printf("#   E₂ (λ²) = η₀·I_C0 = %.6f\n", E2);
    printf("#   E₃ (λ³) = I_mass + I_pot = %.6f\n", E3);
    printf("#   R = E₁ + 2E₂ + 3E₃ = %.6f  (virial residual)\n", R);

    if (R > 0)
        printf("#   R > 0 → soliton wants to EXPAND (too compressed)\n");
    else
        printf("#   R < 0 → soliton wants to CONTRACT (under-compressed)\n");

    /* η₁ from virial balance:
       η₁ = -R / (2 I_C1)
       (with I_C1 already containing the η₀ factor from θ sourcing) */
    double eta1_virial = 0;
    if (fabs(I_C1_sol) > 1e-20) {
        eta1_virial = -R / (2.0 * I_C1_sol);
        printf("\n# η₁ FROM VIRIAL SELF-CONSISTENCY:\n");
        printf("#   η₁ = -R / (2·I_C1) = %.4f\n", eta1_virial);
        printf("#   η_eff at |P|=0.1: %.4f\n", ETA0 + eta1_virial * 0.1);
        printf("#   η_eff at |P|=0.08: %.4f\n", ETA0 + eta1_virial * 0.08);
    } else {
        printf("\n# I_C1 too small — cannot determine η₁\n");
    }

    /* Cross-check with |μ|/η₀ */
    printf("\n# CROSS-CHECKS:\n");
    printf("#   |μ|/η₀ = %.4f\n", fabs(MU)/ETA0);
    printf("#   Ratio η₁_virial / (|μ|/η₀) = %.4f\n",
           eta1_virial / (fabs(MU)/ETA0));

    /* Also compute at multiple thresholds */
    printf("\n# SENSITIVITY TO P_THRESHOLD:\n");
    printf("# (rerun with different threshold to check stability)\n");

    for (int a = 0; a < 3; a++) { free(phi[a]); free(theta[a]); }
    sfa_close(sfa);
    return 0;
}
