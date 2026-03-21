/*  cosserat_cuda.cu — 6-field Cosserat simulation on NVIDIA GPU
 *
 *  CUDA port of the V34 Cosserat equation:
 *    ∂²φ_a/∂t² = ∇²φ_a - m²φ_a - V'(P) + η×curl(θ)_a
 *    ∂²θ_a/∂t² = ∇²θ_a - m_θ²θ_a       + η×curl(φ)_a
 *
 *  One thread per grid point. All 18 field arrays in GPU global memory.
 *  SFA output with async CPU compression.
 *
 *  Build: nvcc -O3 -arch=sm_70 -o sim cosserat_cuda.cu -lzstd -lm -lpthread
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <sys/stat.h>

/* SFA library (CPU side) */
#define SFA_IMPLEMENTATION
#include "../../sfa/format/sfa.h"

#define NFIELDS 3
#define PI 3.14159265358979323846

/* ================================================================
   Physics parameters (constant memory for fast GPU access)
   ================================================================ */

__constant__ double d_MU;
__constant__ double d_KAPPA;
__constant__ double d_MASS2;
__constant__ double d_MTHETA2;
__constant__ double d_ETA;
__constant__ double d_idx2;  /* 1/(dx²) */
__constant__ double d_idx1;  /* 1/(2dx) */
__constant__ int d_N;
__constant__ int d_NN;       /* N² */
__constant__ long d_N3;      /* N³ */

/* Host-side parameters */
static double MU = -41.345, KAPPA = 50.0, MASS2 = 2.25;
static double MTHETA2 = 0.0, A_BG = 0.1, ETA = 0.5;

/* ================================================================
   GPU Kernel: compute forces for all 6 fields
   ================================================================ */

__global__ void compute_forces_kernel(
    double *phi0, double *phi1, double *phi2,
    double *theta0, double *theta1, double *theta2,
    double *acc_phi0, double *acc_phi1, double *acc_phi2,
    double *acc_theta0, double *acc_theta1, double *acc_theta2)
{
    long idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= d_N3) return;

    int N = d_N, NN = d_NN;
    int i = (int)(idx / NN);
    int j = (int)((idx / N) % N);
    int k = (int)(idx % N);

    /* Periodic neighbors */
    int ip = (i+1)%N, im = (i-1+N)%N;
    int jp = (j+1)%N, jm = (j-1+N)%N;
    int kp = (k+1)%N, km = (k-1+N)%N;

    long n_ip = (long)ip*NN + j*N + k;
    long n_im = (long)im*NN + j*N + k;
    long n_jp = (long)i*NN + jp*N + k;
    long n_jm = (long)i*NN + jm*N + k;
    long n_kp = (long)i*NN + j*N + kp;
    long n_km = (long)i*NN + j*N + km;

    /* Load all 6 field values at center */
    double p0 = phi0[idx], p1 = phi1[idx], p2 = phi2[idx];

    /* Triple product + potential derivative (computed once) */
    double P = p0 * p1 * p2;
    double den = 1.0 + d_KAPPA * P * P;
    double mPd2 = d_MU * P / (den * den);

    /* === φ₀ force (fully unrolled) === */
    double lap0 = (phi0[n_ip]+phi0[n_im]+phi0[n_jp]+phi0[n_jm]+phi0[n_kp]+phi0[n_km]-6.0*p0) * d_idx2;
    double curl_t0 = (theta2[n_jp]-theta2[n_jm] - theta1[n_kp]+theta1[n_km]) * d_idx1;
    acc_phi0[idx] = lap0 - d_MASS2*p0 - mPd2*(p1*p2) + d_ETA*curl_t0;

    /* === φ₁ force === */
    double lap1 = (phi1[n_ip]+phi1[n_im]+phi1[n_jp]+phi1[n_jm]+phi1[n_kp]+phi1[n_km]-6.0*p1) * d_idx2;
    double curl_t1 = (theta0[n_kp]-theta0[n_km] - theta2[n_ip]+theta2[n_im]) * d_idx1;
    acc_phi1[idx] = lap1 - d_MASS2*p1 - mPd2*(p0*p2) + d_ETA*curl_t1;

    /* === φ₂ force === */
    double lap2 = (phi2[n_ip]+phi2[n_im]+phi2[n_jp]+phi2[n_jm]+phi2[n_kp]+phi2[n_km]-6.0*p2) * d_idx2;
    double curl_t2 = (theta1[n_ip]-theta1[n_im] - theta0[n_jp]+theta0[n_jm]) * d_idx1;
    acc_phi2[idx] = lap2 - d_MASS2*p2 - mPd2*(p0*p1) + d_ETA*curl_t2;

    /* === θ₀ force === */
    double lapt0 = (theta0[n_ip]+theta0[n_im]+theta0[n_jp]+theta0[n_jm]+theta0[n_kp]+theta0[n_km]-6.0*theta0[idx]) * d_idx2;
    double curl_p0 = (phi2[n_jp]-phi2[n_jm] - phi1[n_kp]+phi1[n_km]) * d_idx1;
    acc_theta0[idx] = lapt0 - d_MTHETA2*theta0[idx] + d_ETA*curl_p0;

    /* === θ₁ force === */
    double lapt1 = (theta1[n_ip]+theta1[n_im]+theta1[n_jp]+theta1[n_jm]+theta1[n_kp]+theta1[n_km]-6.0*theta1[idx]) * d_idx2;
    double curl_p1 = (phi0[n_kp]-phi0[n_km] - phi2[n_ip]+phi2[n_im]) * d_idx1;
    acc_theta1[idx] = lapt1 - d_MTHETA2*theta1[idx] + d_ETA*curl_p1;

    /* === θ₂ force === */
    double lapt2 = (theta2[n_ip]+theta2[n_im]+theta2[n_jp]+theta2[n_jm]+theta2[n_kp]+theta2[n_km]-6.0*theta2[idx]) * d_idx2;
    double curl_p2 = (phi1[n_ip]-phi1[n_im] - phi0[n_jp]+phi0[n_jm]) * d_idx1;
    acc_theta2[idx] = lapt2 - d_MTHETA2*theta2[idx] + d_ETA*curl_p2;
}

/* ================================================================
   GPU Kernel: Verlet half-kick + drift
   ================================================================ */

__global__ void verlet_halfkick_kernel(
    double *vel0, double *vel1, double *vel2,
    double *vel3, double *vel4, double *vel5,
    double *acc0, double *acc1, double *acc2,
    double *acc3, double *acc4, double *acc5,
    double hdt)
{
    long idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= d_N3) return;
    vel0[idx] += hdt * acc0[idx];
    vel1[idx] += hdt * acc1[idx];
    vel2[idx] += hdt * acc2[idx];
    vel3[idx] += hdt * acc3[idx];
    vel4[idx] += hdt * acc4[idx];
    vel5[idx] += hdt * acc5[idx];
}

__global__ void verlet_drift_kernel(
    double *phi0, double *phi1, double *phi2,
    double *theta0, double *theta1, double *theta2,
    double *vel0, double *vel1, double *vel2,
    double *vel3, double *vel4, double *vel5,
    double dt)
{
    long idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= d_N3) return;
    phi0[idx] += dt * vel0[idx];
    phi1[idx] += dt * vel1[idx];
    phi2[idx] += dt * vel2[idx];
    theta0[idx] += dt * vel3[idx];
    theta1[idx] += dt * vel4[idx];
    theta2[idx] += dt * vel5[idx];
}

/* ================================================================
   Host Code
   ================================================================ */

/* GPU device pointers */
static double *d_phi[3], *d_vel_phi[3], *d_acc_phi[3];
static double *d_theta[3], *d_vel_theta[3], *d_acc_theta[3];

static void gpu_alloc(int N) {
    long N3 = (long)N * N * N;
    size_t bytes = N3 * sizeof(double);
    for (int a = 0; a < 3; a++) {
        cudaMalloc(&d_phi[a], bytes);
        cudaMalloc(&d_vel_phi[a], bytes);
        cudaMalloc(&d_acc_phi[a], bytes);
        cudaMalloc(&d_theta[a], bytes);
        cudaMalloc(&d_vel_theta[a], bytes);
        cudaMalloc(&d_acc_theta[a], bytes);
    }
    printf("GPU: allocated %.2f GB for %d³ × 18 arrays\n",
           18.0 * bytes / 1e9, N);
}

static void gpu_upload(double *phi[3], double *vel_phi[3],
                       double *theta[3], double *vel_theta[3], long N3) {
    size_t bytes = N3 * sizeof(double);
    for (int a = 0; a < 3; a++) {
        cudaMemcpy(d_phi[a], phi[a], bytes, cudaMemcpyHostToDevice);
        cudaMemcpy(d_vel_phi[a], vel_phi[a], bytes, cudaMemcpyHostToDevice);
        cudaMemcpy(d_theta[a], theta[a], bytes, cudaMemcpyHostToDevice);
        cudaMemcpy(d_vel_theta[a], vel_theta[a], bytes, cudaMemcpyHostToDevice);
        /* Zero accelerations */
        cudaMemset(d_acc_phi[a], 0, bytes);
        cudaMemset(d_acc_theta[a], 0, bytes);
    }
}

static void gpu_download(double *phi[3], double *theta[3], long N3) {
    size_t bytes = N3 * sizeof(double);
    for (int a = 0; a < 3; a++) {
        cudaMemcpy(phi[a], d_phi[a], bytes, cudaMemcpyDeviceToHost);
        cudaMemcpy(theta[a], d_theta[a], bytes, cudaMemcpyDeviceToHost);
    }
}

static void gpu_set_constants(int N, double dx) {
    int NN = N * N;
    long N3 = (long)N * N * N;
    double idx2 = 1.0 / (dx * dx);
    double idx1 = 1.0 / (2.0 * dx);

    cudaMemcpyToSymbol(d_MU, &MU, sizeof(double));
    cudaMemcpyToSymbol(d_KAPPA, &KAPPA, sizeof(double));
    cudaMemcpyToSymbol(d_MASS2, &MASS2, sizeof(double));
    cudaMemcpyToSymbol(d_MTHETA2, &MTHETA2, sizeof(double));
    cudaMemcpyToSymbol(d_ETA, &ETA, sizeof(double));
    cudaMemcpyToSymbol(d_idx2, &idx2, sizeof(double));
    cudaMemcpyToSymbol(d_idx1, &idx1, sizeof(double));
    cudaMemcpyToSymbol(d_N, &N, sizeof(int));
    cudaMemcpyToSymbol(d_NN, &NN, sizeof(int));
    cudaMemcpyToSymbol(d_N3, &N3, sizeof(long));
}

static void gpu_step(int N, double dt) {
    long N3 = (long)N * N * N;
    int threads = 256;
    int blocks = (N3 + threads - 1) / threads;
    double hdt = 0.5 * dt;

    /* Half-kick */
    verlet_halfkick_kernel<<<blocks, threads>>>(
        d_vel_phi[0], d_vel_phi[1], d_vel_phi[2],
        d_vel_theta[0], d_vel_theta[1], d_vel_theta[2],
        d_acc_phi[0], d_acc_phi[1], d_acc_phi[2],
        d_acc_theta[0], d_acc_theta[1], d_acc_theta[2], hdt);

    /* Drift */
    verlet_drift_kernel<<<blocks, threads>>>(
        d_phi[0], d_phi[1], d_phi[2],
        d_theta[0], d_theta[1], d_theta[2],
        d_vel_phi[0], d_vel_phi[1], d_vel_phi[2],
        d_vel_theta[0], d_vel_theta[1], d_vel_theta[2], dt);

    /* Compute forces */
    compute_forces_kernel<<<blocks, threads>>>(
        d_phi[0], d_phi[1], d_phi[2],
        d_theta[0], d_theta[1], d_theta[2],
        d_acc_phi[0], d_acc_phi[1], d_acc_phi[2],
        d_acc_theta[0], d_acc_theta[1], d_acc_theta[2]);

    /* Half-kick */
    verlet_halfkick_kernel<<<blocks, threads>>>(
        d_vel_phi[0], d_vel_phi[1], d_vel_phi[2],
        d_vel_theta[0], d_vel_theta[1], d_vel_theta[2],
        d_acc_phi[0], d_acc_phi[1], d_acc_phi[2],
        d_acc_theta[0], d_acc_theta[1], d_acc_theta[2], hdt);
}

/* ================================================================
   Initialization (CPU side, then upload)
   ================================================================ */

static void init_braid(double *phi[3], double *vel[3], int N, double L) {
    double dx = 2.0 * L / (N - 1);
    int NN = N * N;
    double A[3] = {0.8, 0.8, 0.8};
    double delta[3] = {0, 3.0005, 4.4325};
    double R_tube = 3.0, ellip = 0.3325;
    double kw = PI / L, omega = sqrt(kw*kw + MASS2);
    double sx = 1+ellip, sy = 1-ellip;
    double inv2R2 = 1.0 / (2*R_tube*R_tube);
    double k_bg = PI/L, omega_bg = sqrt(k_bg*k_bg + MASS2);

    for (int i = 0; i < N; i++) {
        double x = -L + i*dx;
        for (int j = 0; j < N; j++) {
            double y = -L + j*dx;
            for (int kk = 0; kk < N; kk++) {
                double z = -L + kk*dx;
                long idx = (long)i*NN + j*N + kk;
                double r2e = x*x/(sx*sx) + y*y/(sy*sy);
                double env = exp(-r2e * inv2R2);
                for (int a = 0; a < 3; a++) {
                    double ph = kw*z + delta[a];
                    double ph_bg = k_bg*z + 2*PI*a/3.0;
                    phi[a][idx] = A[a]*env*cos(ph) + A_BG*cos(ph_bg);
                    vel[a][idx] = omega*A[a]*env*sin(ph) + omega_bg*A_BG*sin(ph_bg);
                }
            }
        }
    }
}

/* ================================================================
   Main
   ================================================================ */

int main(int argc, char **argv) {
    int N = 128;
    double L = 20.0, T = 300.0;
    double dt_factor = 1.0;
    double diag_dt = 5.0, snap_dt = 50.0;
    char outdir[256] = "output";

    for (int i = 1; i < argc; i++) {
        if      (!strcmp(argv[i],"-N"))    N = atoi(argv[++i]);
        else if (!strcmp(argv[i],"-L"))    L = atof(argv[++i]);
        else if (!strcmp(argv[i],"-T"))    T = atof(argv[++i]);
        else if (!strcmp(argv[i],"-eta"))  ETA = atof(argv[++i]);
        else if (!strcmp(argv[i],"-mt"))   { double m=atof(argv[++i]); MTHETA2=m*m; }
        else if (!strcmp(argv[i],"-m"))    { double m=atof(argv[++i]); MASS2=m*m; }
        else if (!strcmp(argv[i],"-bg"))   A_BG = atof(argv[++i]);
        else if (!strcmp(argv[i],"-diag")) diag_dt = atof(argv[++i]);
        else if (!strcmp(argv[i],"-snap")) snap_dt = atof(argv[++i]);
        else if (!strcmp(argv[i],"-dtf")) dt_factor = atof(argv[++i]);
        else if (!strcmp(argv[i],"-o"))    strncpy(outdir, argv[++i], 255);
    }

    long N3 = (long)N * N * N;
    double dx = 2.0 * L / (N - 1);
    double dt = 0.10 * dx * dt_factor;

    printf("=== V36 Cosserat CUDA — 6-field GPU simulation ===\n");
    printf("N=%d L=%.0f T=%.0f dx=%.4f dt=%.5f\n", N, L, T, dx, dt);
    printf("m²=%.4f m_θ²=%.4f η=%.3f\n\n", MASS2, MTHETA2, ETA);

    /* ---- CPU arrays for init + download ---- */
    double *h_phi[3], *h_vel_phi[3], *h_theta[3], *h_vel_theta[3];
    for (int a = 0; a < 3; a++) {
        h_phi[a] = (double*)calloc(N3, sizeof(double));
        h_vel_phi[a] = (double*)calloc(N3, sizeof(double));
        h_theta[a] = (double*)calloc(N3, sizeof(double));
        h_vel_theta[a] = (double*)calloc(N3, sizeof(double));
    }

    /* ---- Init on CPU ---- */
    init_braid(h_phi, h_vel_phi, N, L);

    /* ---- GPU setup ---- */
    gpu_alloc(N);
    gpu_set_constants(N, dx);
    gpu_upload(h_phi, h_vel_phi, h_theta, h_vel_theta, N3);

    /* Initial force computation */
    int threads = 256;
    int blocks = (N3 + threads - 1) / threads;
    compute_forces_kernel<<<blocks, threads>>>(
        d_phi[0], d_phi[1], d_phi[2],
        d_theta[0], d_theta[1], d_theta[2],
        d_acc_phi[0], d_acc_phi[1], d_acc_phi[2],
        d_acc_theta[0], d_acc_theta[1], d_acc_theta[2]);
    cudaDeviceSynchronize();

    /* ---- SFA output ---- */
    mkdir(outdir, 0755);
    char sfapath[512];
    snprintf(sfapath, sizeof(sfapath), "%s/sim.sfa", outdir);
    SFA *sfa = sfa_create(sfapath, N, N, N, L, L, L, dt);
    sfa_add_column(sfa, "phi_x",   SFA_F64, SFA_POSITION, 0);
    sfa_add_column(sfa, "phi_y",   SFA_F64, SFA_POSITION, 1);
    sfa_add_column(sfa, "phi_z",   SFA_F64, SFA_POSITION, 2);
    sfa_add_column(sfa, "theta_x", SFA_F64, SFA_ANGLE,    0);
    sfa_add_column(sfa, "theta_y", SFA_F64, SFA_ANGLE,    1);
    sfa_add_column(sfa, "theta_z", SFA_F64, SFA_ANGLE,    2);
    sfa_finalize_header(sfa);

    int n_steps = (int)(T / dt);
    int diag_every = (int)(diag_dt / dt); if (diag_every < 1) diag_every = 1;
    int snap_every = (int)(snap_dt / dt); if (snap_every < 1) snap_every = 1;

    printf("Steps=%d diag=%d snap=%d\n", n_steps, diag_every, snap_every);
    printf("SFA: %s\n\n", sfapath);

    /* ---- Main loop ---- */
    cudaEvent_t start, stop;
    cudaEventCreate(&start); cudaEventCreate(&stop);
    cudaEventRecord(start);

    /* Pre-allocate diagnostic buffer (avoid per-step malloc) */
    double *h_vel_buf = (double*)malloc(N3 * sizeof(double));

    for (int step = 0; step <= n_steps; step++) {
        if (step > 0) gpu_step(N, dt);
        double t = step * dt;

        int need_diag = (step % diag_every == 0);
        int need_snap = (step % snap_every == 0);

        /* Only download when needed, and combine diag+snap downloads */
        if (need_diag || need_snap) {
            /* Single download for both diag and snap */
            gpu_download(h_phi, h_theta, N3);

            if (need_snap) {
                void *cols[] = {h_phi[0], h_phi[1], h_phi[2],
                                h_theta[0], h_theta[1], h_theta[2]};
                sfa_write_frame(sfa, t, cols);
            }

            if (need_diag) {
                /* Energy: download only vel_phi[0] for quick KE estimate */
                cudaMemcpy(h_vel_buf, d_vel_phi[0], N3*sizeof(double), cudaMemcpyDeviceToHost);
                double ek = 0, ep = 0;
                for (long idx = 0; idx < N3; idx++) ek += 0.5*h_vel_buf[idx]*h_vel_buf[idx];
                for (long idx = 0; idx < N3; idx++) {
                    double P = h_phi[0][idx]*h_phi[1][idx]*h_phi[2][idx];
                    ep += (MU/2.0)*P*P/(1.0+KAPPA*P*P);
                }
                double dV = dx*dx*dx;
                ek *= dV; ep *= dV;

                double trms = 0;
                for (long idx = 0; idx < N3; idx++)
                    for (int a = 0; a < 3; a++)
                        trms += h_theta[a][idx]*h_theta[a][idx];
                trms = sqrt(trms / (3*N3));

                printf("t=%7.1f E_kin=%.4e E_pot=%.1f θ_rms=%.3e [%.0f%%]\n",
                       t, ek, ep, trms, 100.0*step/n_steps);
            }
        }
    }
    free(h_vel_buf);

    cudaEventRecord(stop);
    cudaEventSynchronize(stop);
    float ms;
    cudaEventElapsedTime(&ms, start, stop);

    printf("\n=== Complete: %.1f seconds (%.1f ms/step) ===\n",
           ms/1000.0, ms/n_steps);

    sfa_close(sfa);

    /* Cleanup */
    for (int a = 0; a < 3; a++) {
        cudaFree(d_phi[a]); cudaFree(d_vel_phi[a]); cudaFree(d_acc_phi[a]);
        cudaFree(d_theta[a]); cudaFree(d_vel_theta[a]); cudaFree(d_acc_theta[a]);
        free(h_phi[a]); free(h_vel_phi[a]);
        free(h_theta[a]); free(h_vel_theta[a]);
    }

    return 0;
}
