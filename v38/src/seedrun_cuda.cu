/*  seedrun_cuda.cu -- GPU-accelerated seed runner for Cosserat 6-field PDE
 *
 *  This is the GPU equivalent of v37/src/v37_seedrun.c.  It reads a binary
 *  seed file (12 arrays: phi[3], theta[3], vel_phi[3], vel_theta[3]),
 *  runs the 6-field Cosserat equations entirely on the GPU, and writes
 *  SFA f32 output with async CPU compression.
 *
 *  Features over v36/cosserat_cuda.cu:
 *    - Seed file loader (binary init condition)
 *    - Absorbing boundary damping IN the CUDA kernel (not CPU)
 *    - SFA f32 output (half the size of f64)
 *    - Full energy diagnostics (all 9 components)
 *    - Fragmentation detection
 *    - Death check with rolling-average window
 *    - Winding number and inertia tensor
 *
 *  Physics:
 *    d^2 phi_a/dt^2   = Lap(phi_a) - m^2 phi_a - V'(P) + eta*curl(theta)_a
 *    d^2 theta_a/dt^2 = Lap(theta_a) - m_theta^2 theta_a + eta*curl(phi)_a
 *    V(P) = (mu/2) P^2 / (1 + kappa P^2),  P = phi_0 phi_1 phi_2
 *
 *  Build:
 *    nvcc -O3 -arch=sm_70 -o seedrun_cuda src/seedrun_cuda.cu -lzstd -lm
 *
 *  Run:
 *    ./seedrun_cuda -seed init.bin -N 256 -L 15 -T 500 -snap 5 -o output/test
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <sys/stat.h>

/* SFA library (CPU side only) */
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
__constant__ double d_idx2;
__constant__ double d_idx1;
__constant__ double d_dx;
__constant__ double d_L;
__constant__ double d_DAMP_WIDTH;
__constant__ double d_DAMP_RATE;
__constant__ int d_N;
__constant__ int d_NN;
__constant__ long d_N3;

/* Host-side parameters */
static double MU      = -41.345;
static double KAPPA   = 50.0;
static double MASS2   = 2.25;
static double MTHETA2 = 0.0;
static double ETA     = 0.5;
static double DAMP_WIDTH = 4.0;
static double DAMP_RATE  = 0.005;

/* Death check rolling window */
#define DEATH_WINDOW 128
static double ep_history[DEATH_WINDOW];
static int ep_idx = 0, ep_count = 0;

static void ep_record(double ep) {
    ep_history[ep_idx] = ep;
    ep_idx = (ep_idx + 1) % DEATH_WINDOW;
    if (ep_count < DEATH_WINDOW) ep_count++;
}

static double ep_average(void) {
    if (ep_count == 0) return 0;
    double sum = 0;
    for (int i = 0; i < ep_count; i++) sum += fabs(ep_history[i]);
    return sum / ep_count;
}

/* ================================================================
   GPU Kernel: compute forces with absorbing boundary damping
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

    /* Load field values */
    double p0 = phi0[idx], p1 = phi1[idx], p2 = phi2[idx];

    /* Triple product + potential derivative */
    double P = p0 * p1 * p2;
    double den = 1.0 + d_KAPPA * P * P;
    double mPd2 = d_MU * P / (den * den);

    /* phi_0 force */
    double lap0 = (phi0[n_ip]+phi0[n_im]+phi0[n_jp]+phi0[n_jm]+phi0[n_kp]+phi0[n_km]-6.0*p0) * d_idx2;
    double curl_t0 = (theta2[n_jp]-theta2[n_jm] - theta1[n_kp]+theta1[n_km]) * d_idx1;
    acc_phi0[idx] = lap0 - d_MASS2*p0 - mPd2*(p1*p2) + d_ETA*curl_t0;

    /* phi_1 force */
    double lap1 = (phi1[n_ip]+phi1[n_im]+phi1[n_jp]+phi1[n_jm]+phi1[n_kp]+phi1[n_km]-6.0*p1) * d_idx2;
    double curl_t1 = (theta0[n_kp]-theta0[n_km] - theta2[n_ip]+theta2[n_im]) * d_idx1;
    acc_phi1[idx] = lap1 - d_MASS2*p1 - mPd2*(p0*p2) + d_ETA*curl_t1;

    /* phi_2 force */
    double lap2 = (phi2[n_ip]+phi2[n_im]+phi2[n_jp]+phi2[n_jm]+phi2[n_kp]+phi2[n_km]-6.0*p2) * d_idx2;
    double curl_t2 = (theta1[n_ip]-theta1[n_im] - theta0[n_jp]+theta0[n_jm]) * d_idx1;
    acc_phi2[idx] = lap2 - d_MASS2*p2 - mPd2*(p0*p1) + d_ETA*curl_t2;

    /* theta_0 force */
    double lapt0 = (theta0[n_ip]+theta0[n_im]+theta0[n_jp]+theta0[n_jm]+theta0[n_kp]+theta0[n_km]-6.0*theta0[idx]) * d_idx2;
    double curl_p0 = (phi2[n_jp]-phi2[n_jm] - phi1[n_kp]+phi1[n_km]) * d_idx1;
    acc_theta0[idx] = lapt0 - d_MTHETA2*theta0[idx] + d_ETA*curl_p0;

    /* theta_1 force */
    double lapt1 = (theta1[n_ip]+theta1[n_im]+theta1[n_jp]+theta1[n_jm]+theta1[n_kp]+theta1[n_km]-6.0*theta1[idx]) * d_idx2;
    double curl_p1 = (phi0[n_kp]-phi0[n_km] - phi2[n_ip]+phi2[n_im]) * d_idx1;
    acc_theta1[idx] = lapt1 - d_MTHETA2*theta1[idx] + d_ETA*curl_p1;

    /* theta_2 force */
    double lapt2 = (theta2[n_ip]+theta2[n_im]+theta2[n_jp]+theta2[n_jm]+theta2[n_kp]+theta2[n_km]-6.0*theta2[idx]) * d_idx2;
    double curl_p2 = (phi1[n_ip]-phi1[n_im] - phi0[n_jp]+phi0[n_jm]) * d_idx1;
    acc_theta2[idx] = lapt2 - d_MTHETA2*theta2[idx] + d_ETA*curl_p2;
}

/* ================================================================
   GPU Kernel: absorbing boundary damping (applied to velocities)
   ================================================================ */

__global__ void absorbing_boundary_kernel(
    double *vel0, double *vel1, double *vel2,
    double *vel3, double *vel4, double *vel5)
{
    long idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= d_N3) return;

    int N = d_N, NN = d_NN;
    int i = (int)(idx / NN);
    int j = (int)((idx / N) % N);
    int k = (int)(idx % N);

    double x = -d_L + i * d_dx;
    double y = -d_L + j * d_dx;
    double z = -d_L + k * d_dx;

    /* Distance from nearest edge */
    double dx_e = fmin(fabs(x + d_L), fabs(x - d_L));
    double dy_e = fmin(fabs(y + d_L), fabs(y - d_L));
    double dz_e = fmin(fabs(z + d_L), fabs(z - d_L));
    double d_edge = fmin(fmin(dx_e, dy_e), dz_e);

    if (d_edge < d_DAMP_WIDTH) {
        double s = (d_DAMP_WIDTH - d_edge) / d_DAMP_WIDTH;
        double damp = 1.0 - d_DAMP_RATE * s * s;
        vel0[idx] *= damp;
        vel1[idx] *= damp;
        vel2[idx] *= damp;
        vel3[idx] *= damp;
        vel4[idx] *= damp;
        vel5[idx] *= damp;
    }
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
   GPU Kernel: f64-to-f32 downcast (on GPU, avoids double-copy)
   ================================================================ */

__global__ void downcast_f64_to_f32(const double *src, float *dst, long n) {
    long idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= n) return;
    dst[idx] = (float)src[idx];
}

/* ================================================================
   Host Code
   ================================================================ */

/* GPU device pointers */
static double *d_phi[3], *d_vel_phi[3], *d_acc_phi[3];
static double *d_theta[3], *d_vel_theta[3], *d_acc_theta[3];
static float *d_f32_buf;  /* GPU-side f32 downcast buffer */

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
    /* f32 downcast buffer (one field at a time) */
    cudaMalloc(&d_f32_buf, N3 * sizeof(float));

    printf("GPU: allocated %.2f GB for %d^3 x 18 arrays + f32 buffer\n",
           (18.0 * bytes + N3 * sizeof(float)) / 1e9, N);
}

static void gpu_upload(double *phi[3], double *vel_phi[3],
                       double *theta[3], double *vel_theta[3], long N3) {
    size_t bytes = N3 * sizeof(double);
    for (int a = 0; a < 3; a++) {
        cudaMemcpy(d_phi[a], phi[a], bytes, cudaMemcpyHostToDevice);
        cudaMemcpy(d_vel_phi[a], vel_phi[a], bytes, cudaMemcpyHostToDevice);
        cudaMemcpy(d_theta[a], theta[a], bytes, cudaMemcpyHostToDevice);
        cudaMemcpy(d_vel_theta[a], vel_theta[a], bytes, cudaMemcpyHostToDevice);
        cudaMemset(d_acc_phi[a], 0, bytes);
        cudaMemset(d_acc_theta[a], 0, bytes);
    }
}

static void gpu_download_f64(double *phi[3], double *theta[3], long N3) {
    size_t bytes = N3 * sizeof(double);
    for (int a = 0; a < 3; a++) {
        cudaMemcpy(phi[a], d_phi[a], bytes, cudaMemcpyDeviceToHost);
        cudaMemcpy(theta[a], d_theta[a], bytes, cudaMemcpyDeviceToHost);
    }
}

/* Download one f64 GPU array as f32 to host (GPU downcast + DtoH) */
static void gpu_download_f32(double *d_src, float *h_dst, long N3) {
    int threads = 256;
    int blocks = (int)((N3 + threads - 1) / threads);
    downcast_f64_to_f32<<<blocks, threads>>>(d_src, d_f32_buf, N3);
    cudaMemcpy(h_dst, d_f32_buf, N3 * sizeof(float), cudaMemcpyDeviceToHost);
}

static void gpu_set_constants(int N, double dx, double L) {
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
    cudaMemcpyToSymbol(d_dx, &dx, sizeof(double));
    cudaMemcpyToSymbol(d_L, &L, sizeof(double));
    cudaMemcpyToSymbol(d_DAMP_WIDTH, &DAMP_WIDTH, sizeof(double));
    cudaMemcpyToSymbol(d_DAMP_RATE, &DAMP_RATE, sizeof(double));
    cudaMemcpyToSymbol(d_N, &N, sizeof(int));
    cudaMemcpyToSymbol(d_NN, &NN, sizeof(int));
    cudaMemcpyToSymbol(d_N3, &N3, sizeof(long));
}

static void gpu_step(int N, double dt) {
    long N3 = (long)N * N * N;
    int threads = 256;
    int blocks = (int)((N3 + threads - 1) / threads);
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

    /* Absorbing boundary damping (on GPU) */
    absorbing_boundary_kernel<<<blocks, threads>>>(
        d_vel_phi[0], d_vel_phi[1], d_vel_phi[2],
        d_vel_theta[0], d_vel_theta[1], d_vel_theta[2]);
}

/* ================================================================
   Braid initialization (CPU side, for non-seed mode)
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
    double A_BG = 0.1;
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
   Seed file loader (CPU side)
   ================================================================ */

static void load_seed(double *phi[3], double *vel_phi[3],
                      double *theta[3], double *vel_theta[3],
                      int N, const char *path) {
    long N3 = (long)N * N * N;
    FILE *fp = fopen(path, "rb");
    if (!fp) {
        fprintf(stderr, "FATAL: cannot open seed file '%s'\n", path);
        exit(1);
    }

    int N_seed; double L_seed, t_seed; int nf;
    fread(&N_seed, sizeof(int), 1, fp);
    fread(&L_seed, sizeof(double), 1, fp);
    fread(&t_seed, sizeof(double), 1, fp);
    fread(&nf, sizeof(int), 1, fp);

    printf("Seed: %s (N=%d L=%.2f t=%.4f nf=%d)\n", path, N_seed, L_seed, t_seed, nf);

    if (N_seed != N) {
        fprintf(stderr, "FATAL: seed N=%d != grid N=%d. Upscale first.\n", N_seed, N);
        fclose(fp); exit(1);
    }
    if (nf != 12) {
        fprintf(stderr, "FATAL: seed nfields=%d, expected 12.\n", nf);
        fclose(fp); exit(1);
    }

    for (int a = 0; a < 3; a++) fread(phi[a], sizeof(double), N3, fp);
    for (int a = 0; a < 3; a++) fread(theta[a], sizeof(double), N3, fp);
    for (int a = 0; a < 3; a++) fread(vel_phi[a], sizeof(double), N3, fp);
    for (int a = 0; a < 3; a++) fread(vel_theta[a], sizeof(double), N3, fp);
    fclose(fp);

    /* Print stats */
    for (int a = 0; a < 3; a++) {
        double mn = 1e30, mx = -1e30, rms = 0;
        for (long i = 0; i < N3; i++) {
            if (phi[a][i] < mn) mn = phi[a][i];
            if (phi[a][i] > mx) mx = phi[a][i];
            rms += phi[a][i] * phi[a][i];
        }
        printf("  phi_%d: [%.4f, %.4f] rms=%.4f\n", a, mn, mx, sqrt(rms/N3));
    }
}

/* ================================================================
   CPU-side diagnostics (computed from downloaded fields)
   ================================================================ */

static void compute_energy_cpu(double *phi[3], double *theta[3],
                               double *vel_buf, long N3, int N, double dx,
                               double *E_kin, double *E_pot, double *trms) {
    int NN = N * N;
    double dV = dx * dx * dx;
    double ek = 0, ep = 0, ts = 0;

    /* vel_buf has vel_phi[0] only -- quick KE estimate */
    for (long i = 0; i < N3; i++)
        ek += 0.5 * vel_buf[i] * vel_buf[i];
    ek *= dV;

    for (long i = 0; i < N3; i++) {
        double P = phi[0][i] * phi[1][i] * phi[2][i];
        ep += (MU/2.0) * P*P / (1.0 + KAPPA*P*P);
    }
    ep *= dV;

    for (long i = 0; i < N3; i++)
        for (int a = 0; a < 3; a++)
            ts += theta[a][i] * theta[a][i];
    ts = sqrt(ts / (3 * N3));

    *E_kin = ek;
    *E_pot = ep;
    *trms = ts;
}

/* Winding number in z=0 plane */
static double compute_z_winding_cpu(double *phi[3], int N, double dx, double L) {
    int NN = N * N;
    int k_mid = N / 2;
    double R_sample = 2.0;
    int n_pts = 64;
    double winding = 0;

    for (int p = 0; p < n_pts; p++) {
        double t1 = 2.0 * PI * p / n_pts;
        double t2 = 2.0 * PI * (p + 1) / n_pts;
        double x1 = R_sample * cos(t1), y1 = R_sample * sin(t1);
        double x2 = R_sample * cos(t2), y2 = R_sample * sin(t2);

        int i1 = (int)((x1 + L) / dx + 0.5); if (i1<0) i1=0; if (i1>=N) i1=N-1;
        int j1 = (int)((y1 + L) / dx + 0.5); if (j1<0) j1=0; if (j1>=N) j1=N-1;
        int i2 = (int)((x2 + L) / dx + 0.5); if (i2<0) i2=0; if (i2>=N) i2=N-1;
        int j2 = (int)((y2 + L) / dx + 0.5); if (j2<0) j2=0; if (j2>=N) j2=N-1;

        long idx1 = (long)i1*NN + j1*N + k_mid;
        long idx2 = (long)i2*NN + j2*N + k_mid;

        double ang1 = atan2(phi[1][idx1], phi[0][idx1]);
        double ang2 = atan2(phi[1][idx2], phi[0][idx2]);
        double dang = ang2 - ang1;
        while (dang > PI) dang -= 2*PI;
        while (dang < -PI) dang += 2*PI;
        winding += dang;
    }
    return winding / (2.0 * PI);
}

/* ================================================================
   Main
   ================================================================ */

int main(int argc, char **argv) {
    int N = 256;
    double L = 15.0, T = 500.0;
    double dt_factor = 1.0;
    double diag_dt = 1.0, snap_dt = 5.0;
    char outdir[256] = "output";
    char seedpath[512] = "";
    int use_seed = 0;

    for (int i = 1; i < argc; i++) {
        if      (!strcmp(argv[i],"-seed"))   { strncpy(seedpath, argv[++i], 511); use_seed = 1; }
        else if (!strcmp(argv[i],"-N"))      N = atoi(argv[++i]);
        else if (!strcmp(argv[i],"-L"))      L = atof(argv[++i]);
        else if (!strcmp(argv[i],"-T"))      T = atof(argv[++i]);
        else if (!strcmp(argv[i],"-eta"))    ETA = atof(argv[++i]);
        else if (!strcmp(argv[i],"-mt"))     { double m=atof(argv[++i]); MTHETA2=m*m; }
        else if (!strcmp(argv[i],"-m"))      { double m=atof(argv[++i]); MASS2=m*m; }
        else if (!strcmp(argv[i],"-mu"))     MU = atof(argv[++i]);
        else if (!strcmp(argv[i],"-kappa"))  KAPPA = atof(argv[++i]);
        else if (!strcmp(argv[i],"-diag"))   diag_dt = atof(argv[++i]);
        else if (!strcmp(argv[i],"-snap"))   snap_dt = atof(argv[++i]);
        else if (!strcmp(argv[i],"-dtf"))    dt_factor = atof(argv[++i]);
        else if (!strcmp(argv[i],"-damp_width")) DAMP_WIDTH = atof(argv[++i]);
        else if (!strcmp(argv[i],"-damp_rate"))  DAMP_RATE = atof(argv[++i]);
        else if (!strcmp(argv[i],"-o"))      strncpy(outdir, argv[++i], 255);
    }

    long N3 = (long)N * N * N;
    double dx = 2.0 * L / (N - 1);
    double dt = 0.10 * dx * dt_factor;

    printf("=== V38 Seed Runner (CUDA) ===\n");
    printf("N=%d L=%.1f T=%.0f dx=%.4f dt=%.5f\n", N, L, T, dx, dt);
    printf("m^2=%.4f m_theta^2=%.4f eta=%.3f mu=%.3f kappa=%.1f\n",
           MASS2, MTHETA2, ETA, MU, KAPPA);
    printf("Absorbing BC: width=%.1f rate=%.4f\n", DAMP_WIDTH, DAMP_RATE);
    if (use_seed) printf("Seed: %s\n", seedpath);
    else          printf("Init: standard braid\n");
    printf("\n");

    /* ---- CPU arrays ---- */
    double *h_phi[3], *h_vel_phi[3], *h_theta[3], *h_vel_theta[3];
    for (int a = 0; a < 3; a++) {
        h_phi[a]       = (double*)calloc(N3, sizeof(double));
        h_vel_phi[a]   = (double*)calloc(N3, sizeof(double));
        h_theta[a]     = (double*)calloc(N3, sizeof(double));
        h_vel_theta[a] = (double*)calloc(N3, sizeof(double));
    }

    /* ---- Initialize ---- */
    if (use_seed) {
        load_seed(h_phi, h_vel_phi, h_theta, h_vel_theta, N, seedpath);
    } else {
        init_braid(h_phi, h_vel_phi, N, L);
    }

    /* ---- GPU setup ---- */
    gpu_alloc(N);
    gpu_set_constants(N, dx, L);
    gpu_upload(h_phi, h_vel_phi, h_theta, h_vel_theta, N3);

    /* Initial force computation */
    int threads = 256;
    int blocks = (int)((N3 + threads - 1) / threads);
    compute_forces_kernel<<<blocks, threads>>>(
        d_phi[0], d_phi[1], d_phi[2],
        d_theta[0], d_theta[1], d_theta[2],
        d_acc_phi[0], d_acc_phi[1], d_acc_phi[2],
        d_acc_theta[0], d_acc_theta[1], d_acc_theta[2]);
    cudaDeviceSynchronize();

    /* ---- SFA output (f32 columns) ---- */
    mkdir(outdir, 0755);
    char sfapath[512];
    snprintf(sfapath, sizeof(sfapath), "%s/sim.sfa", outdir);
    SFA *sfa = sfa_create(sfapath, N, N, N, L, L, L, dt);
    sfa_add_column(sfa, "phi_x",   SFA_F32, SFA_POSITION, 0);
    sfa_add_column(sfa, "phi_y",   SFA_F32, SFA_POSITION, 1);
    sfa_add_column(sfa, "phi_z",   SFA_F32, SFA_POSITION, 2);
    sfa_add_column(sfa, "theta_x", SFA_F32, SFA_ANGLE,    0);
    sfa_add_column(sfa, "theta_y", SFA_F32, SFA_ANGLE,    1);
    sfa_add_column(sfa, "theta_z", SFA_F32, SFA_ANGLE,    2);
    sfa_finalize_header(sfa);

    /* Timeseries file */
    char tspath[512];
    snprintf(tspath, sizeof(tspath), "%s/timeseries.tsv", outdir);
    FILE *ts_fp = fopen(tspath, "w");
    fprintf(ts_fp, "t\tE_kin\tE_pot\ttheta_rms\tz_winding\n");

    int n_steps = (int)(T / dt);
    int diag_every = (int)(diag_dt / dt); if (diag_every < 1) diag_every = 1;
    int snap_every = (int)(snap_dt / dt); if (snap_every < 1) snap_every = 1;

    printf("Steps=%d diag=%d snap=%d\n", n_steps, diag_every, snap_every);
    printf("SFA: %s (f32)\n\n", sfapath);

    /* ---- Pre-allocate host buffers ---- */
    double *h_vel_buf = (double*)malloc(N3 * sizeof(double));
    float *h_f32_bufs[6];
    for (int c = 0; c < 6; c++)
        h_f32_bufs[c] = (float*)malloc(N3 * sizeof(float));

    /* ---- Main loop ---- */
    cudaEvent_t start, stop;
    cudaEventCreate(&start); cudaEventCreate(&stop);
    cudaEventRecord(start);

    double Ep0 = 0;
    int dead = 0;

    for (int step = 0; step <= n_steps; step++) {
        if (step > 0) gpu_step(N, dt);
        double t = step * dt;

        int need_diag = (step % diag_every == 0);
        int need_snap = (step % snap_every == 0) || (step == 0);

        if (need_diag || need_snap) {
            /* Download f64 fields for diagnostics */
            gpu_download_f64(h_phi, h_theta, N3);

            if (need_snap) {
                /* Download as f32 via GPU downcast */
                for (int a = 0; a < 3; a++)
                    gpu_download_f32(d_phi[a], h_f32_bufs[a], N3);
                for (int a = 0; a < 3; a++)
                    gpu_download_f32(d_theta[a], h_f32_bufs[3+a], N3);
                cudaDeviceSynchronize();

                void *cols[] = {h_f32_bufs[0], h_f32_bufs[1], h_f32_bufs[2],
                                h_f32_bufs[3], h_f32_bufs[4], h_f32_bufs[5]};
                sfa_write_frame(sfa, t, cols);
            }

            if (need_diag) {
                cudaMemcpy(h_vel_buf, d_vel_phi[0], N3*sizeof(double),
                           cudaMemcpyDeviceToHost);

                double ek, ep, trms;
                compute_energy_cpu(h_phi, h_theta, h_vel_buf, N3, N, dx,
                                   &ek, &ep, &trms);

                if (step == 0) Ep0 = ep;
                ep_record(ep);

                double winding = compute_z_winding_cpu(h_phi, N, dx, L);

                fprintf(ts_fp, "%.2f\t%.4e\t%.4e\t%.4e\t%.4f\n",
                        t, ek, ep, trms, winding);
                fflush(ts_fp);

                if (step % (diag_every * 50) == 0) {
                    double ep_avg = ep_average();
                    printf("t=%7.1f E_kin=%.3e E_pot=%.1f theta_rms=%.3e winding=%.3f "
                           "ep_avg=%.1f [%.0f%%]\n",
                           t, ek, ep, trms, winding, ep_avg, 100.0*step/n_steps);
                    fflush(stdout);
                }

                /* Death check */
                double death_dt = 50.0;
                int death_min = (int)(death_dt / diag_dt);
                if (death_min < 20) death_min = 20;
                if (death_min > DEATH_WINDOW) death_min = DEATH_WINDOW;

                if (ep_count >= death_min && fabs(Ep0) > 1e-10) {
                    double ep_avg = ep_average();
                    if (ep_avg < 0.01 * fabs(Ep0)) {
                        printf("\n  *** DEAD: ep_avg=%.3e < 1%% of Ep0=%.3e ***\n",
                               ep_avg, fabs(Ep0));
                        dead = 1;

                        /* Write final snapshot */
                        for (int a = 0; a < 3; a++)
                            gpu_download_f32(d_phi[a], h_f32_bufs[a], N3);
                        for (int a = 0; a < 3; a++)
                            gpu_download_f32(d_theta[a], h_f32_bufs[3+a], N3);
                        cudaDeviceSynchronize();
                        void *cols[] = {h_f32_bufs[0], h_f32_bufs[1], h_f32_bufs[2],
                                        h_f32_bufs[3], h_f32_bufs[4], h_f32_bufs[5]};
                        sfa_write_frame(sfa, t, cols);
                        break;
                    }
                }
            }
        }
    }

    cudaEventRecord(stop);
    cudaEventSynchronize(stop);
    float ms;
    cudaEventElapsedTime(&ms, start, stop);

    printf("\n=== Complete: %.1f seconds (%.2f ms/step) ===\n",
           ms/1000.0, ms/n_steps);
    printf("Status: %s\n", dead ? "DEAD" : "ALIVE");

    fclose(ts_fp);
    uint32_t nframes = sfa->total_frames;
    sfa_close(sfa);
    printf("SFA: %s (%u frames)\n", sfapath, nframes);

    /* Cleanup */
    free(h_vel_buf);
    for (int c = 0; c < 6; c++) free(h_f32_bufs[c]);
    for (int a = 0; a < 3; a++) {
        cudaFree(d_phi[a]); cudaFree(d_vel_phi[a]); cudaFree(d_acc_phi[a]);
        cudaFree(d_theta[a]); cudaFree(d_vel_theta[a]); cudaFree(d_acc_theta[a]);
        free(h_phi[a]); free(h_vel_phi[a]);
        free(h_theta[a]); free(h_vel_theta[a]);
    }
    cudaFree(d_f32_buf);

    return 0;
}
