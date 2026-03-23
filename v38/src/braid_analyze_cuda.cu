/*  braid_analyze_cuda.cu -- GPU braid characterization tool
 *
 *  Initializes the standard V33/V34 braid, runs at high resolution
 *  (N=256 or N=512), and every ANALYSIS_EVERY timesteps computes:
 *
 *    1. Radial theta profile at z=0 (azimuthally averaged)
 *    2. Comparison with 1/r Biot-Savart prediction
 *    3. Winding number of phi field in xy-plane
 *    4. Energy decomposition (kin, grad, mass, pot, coupling)
 *    5. Braid centroid velocity along z
 *
 *  Outputs SFA for field visualization and TSV for analysis.
 *
 *  Physics (6-field Cosserat, same as seedrun_cuda.cu):
 *    d^2 phi_a/dt^2   = Lap(phi_a) - m^2 phi_a - V'(P) + eta*curl(theta)_a
 *    d^2 theta_a/dt^2 = Lap(theta_a) - m_theta^2 theta_a + eta*curl(phi)_a
 *
 *  Build:
 *    nvcc -O3 -arch=sm_70 -o braid_analyze src/braid_analyze_cuda.cu -lzstd -lm
 *
 *  Run:
 *    ./braid_analyze -N 256 -L 20 -T 200 -eta 0.5 -o output/braid_n256
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <sys/stat.h>

#define SFA_IMPLEMENTATION
#include "../../sfa/format/sfa.h"

#define NFIELDS 3
#define PI 3.14159265358979323846
#define N_RADIAL_BINS 100
#define N_ANGLE_SAMPLE 64

/* ================================================================
   Physics parameters (constant memory)
   ================================================================ */

__constant__ double d_MU;
__constant__ double d_KAPPA;
__constant__ double d_MASS2;
__constant__ double d_MTHETA2;
__constant__ double d_ETA;
__constant__ double d_idx2;
__constant__ double d_idx1;
__constant__ int d_N;
__constant__ int d_NN;
__constant__ long d_N3;

static double MU = -41.345, KAPPA = 50.0, MASS2 = 2.25;
static double MTHETA2 = 0.0, ETA = 0.5, A_BG = 0.1;

/* ================================================================
   GPU Kernels (same as seedrun_cuda.cu)
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

    int ip = (i+1)%N, im = (i-1+N)%N;
    int jp = (j+1)%N, jm = (j-1+N)%N;
    int kp = (k+1)%N, km = (k-1+N)%N;

    long n_ip = (long)ip*NN + j*N + k;
    long n_im = (long)im*NN + j*N + k;
    long n_jp = (long)i*NN + jp*N + k;
    long n_jm = (long)i*NN + jm*N + k;
    long n_kp = (long)i*NN + j*N + kp;
    long n_km = (long)i*NN + j*N + km;

    double p0 = phi0[idx], p1 = phi1[idx], p2 = phi2[idx];
    double P = p0 * p1 * p2;
    double den = 1.0 + d_KAPPA * P * P;
    double mPd2 = d_MU * P / (den * den);

    double lap0 = (phi0[n_ip]+phi0[n_im]+phi0[n_jp]+phi0[n_jm]+phi0[n_kp]+phi0[n_km]-6.0*p0) * d_idx2;
    double curl_t0 = (theta2[n_jp]-theta2[n_jm] - theta1[n_kp]+theta1[n_km]) * d_idx1;
    acc_phi0[idx] = lap0 - d_MASS2*p0 - mPd2*(p1*p2) + d_ETA*curl_t0;

    double lap1 = (phi1[n_ip]+phi1[n_im]+phi1[n_jp]+phi1[n_jm]+phi1[n_kp]+phi1[n_km]-6.0*p1) * d_idx2;
    double curl_t1 = (theta0[n_kp]-theta0[n_km] - theta2[n_ip]+theta2[n_im]) * d_idx1;
    acc_phi1[idx] = lap1 - d_MASS2*p1 - mPd2*(p0*p2) + d_ETA*curl_t1;

    double lap2 = (phi2[n_ip]+phi2[n_im]+phi2[n_jp]+phi2[n_jm]+phi2[n_kp]+phi2[n_km]-6.0*p2) * d_idx2;
    double curl_t2 = (theta1[n_ip]-theta1[n_im] - theta0[n_jp]+theta0[n_jm]) * d_idx1;
    acc_phi2[idx] = lap2 - d_MASS2*p2 - mPd2*(p0*p1) + d_ETA*curl_t2;

    double lapt0 = (theta0[n_ip]+theta0[n_im]+theta0[n_jp]+theta0[n_jm]+theta0[n_kp]+theta0[n_km]-6.0*theta0[idx]) * d_idx2;
    double curl_p0 = (phi2[n_jp]-phi2[n_jm] - phi1[n_kp]+phi1[n_km]) * d_idx1;
    acc_theta0[idx] = lapt0 - d_MTHETA2*theta0[idx] + d_ETA*curl_p0;

    double lapt1 = (theta1[n_ip]+theta1[n_im]+theta1[n_jp]+theta1[n_jm]+theta1[n_kp]+theta1[n_km]-6.0*theta1[idx]) * d_idx2;
    double curl_p1 = (phi0[n_kp]-phi0[n_km] - phi2[n_ip]+phi2[n_im]) * d_idx1;
    acc_theta1[idx] = lapt1 - d_MTHETA2*theta1[idx] + d_ETA*curl_p1;

    double lapt2 = (theta2[n_ip]+theta2[n_im]+theta2[n_jp]+theta2[n_jm]+theta2[n_kp]+theta2[n_km]-6.0*theta2[idx]) * d_idx2;
    double curl_p2 = (phi1[n_ip]-phi1[n_im] - phi0[n_jp]+phi0[n_jm]) * d_idx1;
    acc_theta2[idx] = lapt2 - d_MTHETA2*theta2[idx] + d_ETA*curl_p2;
}

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

__global__ void downcast_f64_to_f32(const double *src, float *dst, long n) {
    long idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= n) return;
    dst[idx] = (float)src[idx];
}

/* ================================================================
   Host GPU management
   ================================================================ */

static double *d_phi[3], *d_vel_phi[3], *d_acc_phi[3];
static double *d_theta[3], *d_vel_theta[3], *d_acc_theta[3];
static float *d_f32_buf;

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
    cudaMalloc(&d_f32_buf, N3 * sizeof(float));
    printf("GPU: allocated %.2f GB for %d^3 grid\n",
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

static void gpu_download(double *phi[3], double *theta[3], long N3) {
    size_t bytes = N3 * sizeof(double);
    for (int a = 0; a < 3; a++) {
        cudaMemcpy(phi[a], d_phi[a], bytes, cudaMemcpyDeviceToHost);
        cudaMemcpy(theta[a], d_theta[a], bytes, cudaMemcpyDeviceToHost);
    }
}

static void gpu_download_vel(double *vel_phi[3], double *vel_theta[3], long N3) {
    size_t bytes = N3 * sizeof(double);
    for (int a = 0; a < 3; a++) {
        cudaMemcpy(vel_phi[a], d_vel_phi[a], bytes, cudaMemcpyDeviceToHost);
        cudaMemcpy(vel_theta[a], d_vel_theta[a], bytes, cudaMemcpyDeviceToHost);
    }
}

static void gpu_download_f32(double *d_src, float *h_dst, long N3) {
    int threads = 256;
    int blocks = (int)((N3 + threads - 1) / threads);
    downcast_f64_to_f32<<<blocks, threads>>>(d_src, d_f32_buf, N3);
    cudaMemcpy(h_dst, d_f32_buf, N3 * sizeof(float), cudaMemcpyDeviceToHost);
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
    int blocks = (int)((N3 + threads - 1) / threads);
    double hdt = 0.5 * dt;

    verlet_halfkick_kernel<<<blocks, threads>>>(
        d_vel_phi[0], d_vel_phi[1], d_vel_phi[2],
        d_vel_theta[0], d_vel_theta[1], d_vel_theta[2],
        d_acc_phi[0], d_acc_phi[1], d_acc_phi[2],
        d_acc_theta[0], d_acc_theta[1], d_acc_theta[2], hdt);

    verlet_drift_kernel<<<blocks, threads>>>(
        d_phi[0], d_phi[1], d_phi[2],
        d_theta[0], d_theta[1], d_theta[2],
        d_vel_phi[0], d_vel_phi[1], d_vel_phi[2],
        d_vel_theta[0], d_vel_theta[1], d_vel_theta[2], dt);

    compute_forces_kernel<<<blocks, threads>>>(
        d_phi[0], d_phi[1], d_phi[2],
        d_theta[0], d_theta[1], d_theta[2],
        d_acc_phi[0], d_acc_phi[1], d_acc_phi[2],
        d_acc_theta[0], d_acc_theta[1], d_acc_theta[2]);

    verlet_halfkick_kernel<<<blocks, threads>>>(
        d_vel_phi[0], d_vel_phi[1], d_vel_phi[2],
        d_vel_theta[0], d_vel_theta[1], d_vel_theta[2],
        d_acc_phi[0], d_acc_phi[1], d_acc_phi[2],
        d_acc_theta[0], d_acc_theta[1], d_acc_theta[2], hdt);
}

/* ================================================================
   Braid initialization (same as v36/v33)
   ================================================================ */

static void init_braid(double *phi[3], double *vel[3], int N, double L,
                       double x_cen, double y_cen) {
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
                double xc = x - x_cen, yc = y - y_cen;
                double r2e = xc*xc/(sx*sx) + yc*yc/(sy*sy);
                double env = exp(-r2e * inv2R2);
                for (int a = 0; a < 3; a++) {
                    double ph = kw*z + delta[a];
                    double ph_bg = k_bg*z + 2*PI*a/3.0;
                    phi[a][idx] += A[a]*env*cos(ph) + A_BG*cos(ph_bg);
                    vel[a][idx] += omega*A[a]*env*sin(ph) + omega_bg*A_BG*sin(ph_bg);
                }
            }
        }
    }
}

/* ================================================================
   Analysis: azimuthally-averaged radial theta profile at z=z_slice
   ================================================================ */

static void compute_radial_theta_profile(
    double *theta[3], int N, double dx, double L,
    double z_slice, int n_bins, double R_max,
    double *r_centers, double *theta_rms_profile, double *theta_mag_profile)
{
    int NN = N * N;
    double dr = R_max / n_bins;

    /* Find the k index closest to z_slice */
    int k_slice = (int)((z_slice + L) / dx + 0.5);
    if (k_slice < 0) k_slice = 0;
    if (k_slice >= N) k_slice = N - 1;

    double *sum2 = (double*)calloc(n_bins, sizeof(double));
    double *sum_mag = (double*)calloc(n_bins, sizeof(double));
    int *count = (int*)calloc(n_bins, sizeof(int));

    for (int i = 0; i < N; i++) {
        double x = -L + i * dx;
        for (int j = 0; j < N; j++) {
            double y = -L + j * dx;
            double r = sqrt(x*x + y*y);
            int bin = (int)(r / dr);
            if (bin >= n_bins) continue;

            long idx = (long)i * NN + j * N + k_slice;
            double t2 = 0;
            for (int a = 0; a < 3; a++)
                t2 += theta[a][idx] * theta[a][idx];
            sum2[bin] += t2;
            sum_mag[bin] += sqrt(t2);
            count[bin]++;
        }
    }

    for (int b = 0; b < n_bins; b++) {
        r_centers[b] = (b + 0.5) * dr;
        if (count[b] > 0) {
            theta_rms_profile[b] = sqrt(sum2[b] / count[b]);
            theta_mag_profile[b] = sum_mag[b] / count[b];
        } else {
            theta_rms_profile[b] = 0;
            theta_mag_profile[b] = 0;
        }
    }

    free(sum2); free(sum_mag); free(count);
}

/* ================================================================
   Analysis: full energy decomposition
   ================================================================ */

static void compute_full_energy(
    double *phi[3], double *theta[3],
    double *vel_phi[3], double *vel_theta[3],
    int N, double dx,
    double *E_phi_kin, double *E_theta_kin,
    double *E_grad, double *E_mass, double *E_pot,
    double *E_tgrad, double *E_tmass, double *E_coupling)
{
    int NN = N * N;
    long N3 = (long)N * N * N;
    double dV = dx * dx * dx;
    double idx1 = 1.0 / (2.0 * dx);
    double epk=0, etk=0, eg=0, em=0, ep=0, etg=0, etm=0, ec=0;

    for (long idx = 0; idx < N3; idx++) {
        int i=(int)(idx/NN), j=(int)((idx/N)%N), k=(int)(idx%N);
        int ip=(i+1)%N, im=(i-1+N)%N;
        int jp=(j+1)%N, jm=(j-1+N)%N;
        int kp=(k+1)%N, km=(k-1+N)%N;
        long n_ip=(long)ip*NN+j*N+k, n_im=(long)im*NN+j*N+k;
        long n_jp=(long)i*NN+jp*N+k, n_jm=(long)i*NN+jm*N+k;
        long n_kp=(long)i*NN+j*N+kp, n_km=(long)i*NN+j*N+km;

        for (int a = 0; a < 3; a++) {
            epk += 0.5*vel_phi[a][idx]*vel_phi[a][idx]*dV;
            etk += 0.5*vel_theta[a][idx]*vel_theta[a][idx]*dV;
            double gx=(phi[a][n_ip]-phi[a][n_im])*idx1;
            double gy=(phi[a][n_jp]-phi[a][n_jm])*idx1;
            double gz=(phi[a][n_kp]-phi[a][n_km])*idx1;
            eg += 0.5*(gx*gx+gy*gy+gz*gz)*dV;
            em += 0.5*MASS2*phi[a][idx]*phi[a][idx]*dV;
            double tgx=(theta[a][n_ip]-theta[a][n_im])*idx1;
            double tgy=(theta[a][n_jp]-theta[a][n_jm])*idx1;
            double tgz=(theta[a][n_kp]-theta[a][n_km])*idx1;
            etg += 0.5*(tgx*tgx+tgy*tgy+tgz*tgz)*dV;
            etm += 0.5*MTHETA2*theta[a][idx]*theta[a][idx]*dV;
        }
        double P = phi[0][idx]*phi[1][idx]*phi[2][idx];
        ep += (MU/2.0)*P*P/(1.0+KAPPA*P*P)*dV;

        /* Coupling: -eta * phi_a * curl(theta)_a */
        double ct0 = (theta[2][n_jp]-theta[2][n_jm] - theta[1][n_kp]+theta[1][n_km]) * idx1;
        double ct1 = (theta[0][n_kp]-theta[0][n_km] - theta[2][n_ip]+theta[2][n_im]) * idx1;
        double ct2 = (theta[1][n_ip]-theta[1][n_im] - theta[0][n_jp]+theta[0][n_jm]) * idx1;
        ec -= ETA * (phi[0][idx]*ct0 + phi[1][idx]*ct1 + phi[2][idx]*ct2) * dV;
    }
    *E_phi_kin=epk; *E_theta_kin=etk; *E_grad=eg; *E_mass=em;
    *E_pot=ep; *E_tgrad=etg; *E_tmass=etm; *E_coupling=ec;
}

/* ================================================================
   Analysis: winding number in z=0 plane
   ================================================================ */

static double compute_winding(double *phi[3], int N, double dx, double L) {
    int NN = N * N;
    int k_mid = N / 2;
    double R_sample = 2.0;
    double winding = 0;

    for (int p = 0; p < N_ANGLE_SAMPLE; p++) {
        double t1 = 2.0 * PI * p / N_ANGLE_SAMPLE;
        double t2 = 2.0 * PI * (p + 1) / N_ANGLE_SAMPLE;
        double x1 = R_sample * cos(t1), y1 = R_sample * sin(t1);
        double x2 = R_sample * cos(t2), y2 = R_sample * sin(t2);

        int i1 = (int)((x1+L)/dx+0.5); if(i1<0)i1=0; if(i1>=N)i1=N-1;
        int j1 = (int)((y1+L)/dx+0.5); if(j1<0)j1=0; if(j1>=N)j1=N-1;
        int i2 = (int)((x2+L)/dx+0.5); if(i2<0)i2=0; if(i2>=N)i2=N-1;
        int j2 = (int)((y2+L)/dx+0.5); if(j2<0)j2=0; if(j2>=N)j2=N-1;

        long idx1 = (long)i1*NN+j1*N+k_mid;
        long idx2 = (long)i2*NN+j2*N+k_mid;

        double ang1 = atan2(phi[1][idx1], phi[0][idx1]);
        double ang2 = atan2(phi[1][idx2], phi[0][idx2]);
        double dang = ang2 - ang1;
        while(dang > PI) dang -= 2*PI;
        while(dang < -PI) dang += 2*PI;
        winding += dang;
    }
    return winding / (2.0 * PI);
}

/* Braid z-centroid for velocity tracking */
static double braid_z_centroid(double *phi[3], int N, double dx, double L) {
    int NN = N * N;
    long N3 = (long)N * N * N;
    double dV = dx * dx * dx;

    double avg = 0;
    for (long i = 0; i < N3; i++) {
        double p2 = 0;
        for (int a = 0; a < 3; a++) p2 += phi[a][i]*phi[a][i];
        avg += p2;
    }
    avg /= N3;

    double wz = 0, wt = 0;
    for (long idx = 0; idx < N3; idx++) {
        double p2 = 0;
        for (int a = 0; a < 3; a++) p2 += phi[a][idx]*phi[a][idx];
        if (p2 < 5.0*avg) continue;
        int k = (int)(idx % N);
        double z = -L + k*dx;
        wz += z*p2; wt += p2;
    }
    return (wt > 0) ? wz/wt : 0;
}

/* ================================================================
   Main
   ================================================================ */

int main(int argc, char **argv) {
    int N = 256;
    double L = 20.0, T = 200.0;
    double dt_factor = 1.0;
    double snap_dt = 10.0;
    int analysis_every = 10;  /* timesteps between analyses */
    int n_braids = 1;
    double D = 20.0;
    char outdir[256] = "output";

    for (int i = 1; i < argc; i++) {
        if      (!strcmp(argv[i],"-N"))      N = atoi(argv[++i]);
        else if (!strcmp(argv[i],"-L"))      L = atof(argv[++i]);
        else if (!strcmp(argv[i],"-T"))      T = atof(argv[++i]);
        else if (!strcmp(argv[i],"-eta"))    ETA = atof(argv[++i]);
        else if (!strcmp(argv[i],"-mt"))     { double m=atof(argv[++i]); MTHETA2=m*m; }
        else if (!strcmp(argv[i],"-m"))      { double m=atof(argv[++i]); MASS2=m*m; }
        else if (!strcmp(argv[i],"-bg"))     A_BG = atof(argv[++i]);
        else if (!strcmp(argv[i],"-dtf"))    dt_factor = atof(argv[++i]);
        else if (!strcmp(argv[i],"-snap"))   snap_dt = atof(argv[++i]);
        else if (!strcmp(argv[i],"-aevery")) analysis_every = atoi(argv[++i]);
        else if (!strcmp(argv[i],"-braids")) n_braids = atoi(argv[++i]);
        else if (!strcmp(argv[i],"-D"))      D = atof(argv[++i]);
        else if (!strcmp(argv[i],"-o"))      strncpy(outdir, argv[++i], 255);
    }

    long N3 = (long)N * N * N;
    double dx = 2.0 * L / (N - 1);
    double dt = 0.10 * dx * dt_factor;

    printf("=== V38 Braid Analyzer (CUDA) ===\n");
    printf("N=%d L=%.0f T=%.0f dx=%.4f dt=%.5f\n", N, L, T, dx, dt);
    printf("m^2=%.4f m_theta^2=%.4f eta=%.3f braids=%d\n\n", MASS2, MTHETA2, ETA, n_braids);

    /* ---- CPU arrays ---- */
    double *h_phi[3], *h_vel_phi[3], *h_theta[3], *h_vel_theta[3];
    for (int a = 0; a < 3; a++) {
        h_phi[a]       = (double*)calloc(N3, sizeof(double));
        h_vel_phi[a]   = (double*)calloc(N3, sizeof(double));
        h_theta[a]     = (double*)calloc(N3, sizeof(double));
        h_vel_theta[a] = (double*)calloc(N3, sizeof(double));
    }

    /* Initialize braid(s) */
    if (n_braids == 1) {
        init_braid(h_phi, h_vel_phi, N, L, 0, 0);
    } else {
        for (int b = 0; b < n_braids; b++) {
            double bx = (n_braids==2) ? (b==0?-D/2:D/2) : D/2*cos(2*PI*b/n_braids);
            double by = (n_braids==2) ? 0 : D/2*sin(2*PI*b/n_braids);
            init_braid(h_phi, h_vel_phi, N, L, bx, by);
        }
    }

    /* ---- GPU setup ---- */
    gpu_alloc(N);
    gpu_set_constants(N, dx);
    gpu_upload(h_phi, h_vel_phi, h_theta, h_vel_theta, N3);

    int threads = 256;
    int blocks = (int)((N3 + threads - 1) / threads);
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
    sfa_add_column(sfa, "phi_x",   SFA_F32, SFA_POSITION, 0);
    sfa_add_column(sfa, "phi_y",   SFA_F32, SFA_POSITION, 1);
    sfa_add_column(sfa, "phi_z",   SFA_F32, SFA_POSITION, 2);
    sfa_add_column(sfa, "theta_x", SFA_F32, SFA_ANGLE,    0);
    sfa_add_column(sfa, "theta_y", SFA_F32, SFA_ANGLE,    1);
    sfa_add_column(sfa, "theta_z", SFA_F32, SFA_ANGLE,    2);
    sfa_finalize_header(sfa);

    /* ---- Output files ---- */
    char tspath[512], profpath[512];
    snprintf(tspath, sizeof(tspath), "%s/timeseries.tsv", outdir);
    snprintf(profpath, sizeof(profpath), "%s/radial_profiles.tsv", outdir);
    FILE *ts_fp = fopen(tspath, "w");
    FILE *prof_fp = fopen(profpath, "w");
    fprintf(ts_fp, "t\tE_phi_kin\tE_theta_kin\tE_grad\tE_mass\tE_pot\t"
                   "E_tgrad\tE_tmass\tE_coupling\tE_total\twinding\tz_centroid\n");
    fprintf(prof_fp, "# Radial theta profile: t, r, theta_rms, theta_mag, biot_savart_1_over_r\n");

    float *h_f32_bufs[6];
    for (int c = 0; c < 6; c++)
        h_f32_bufs[c] = (float*)malloc(N3 * sizeof(float));

    int n_steps = (int)(T / dt);
    int snap_every = (int)(snap_dt / dt); if (snap_every < 1) snap_every = 1;
    double R_max_profile = L * 0.8;

    printf("Steps=%d snap=%d analysis_every=%d steps\n", n_steps, snap_every, analysis_every);
    printf("SFA: %s\n\n", sfapath);

    /* ---- Main loop ---- */
    cudaEvent_t start, stop;
    cudaEventCreate(&start); cudaEventCreate(&stop);
    cudaEventRecord(start);

    double prev_z_cen = 0;
    double prev_t = 0;

    for (int step = 0; step <= n_steps; step++) {
        if (step > 0) gpu_step(N, dt);
        double t = step * dt;

        int need_analysis = (step % analysis_every == 0);
        int need_snap = (step % snap_every == 0) || (step == 0);

        if (need_analysis || need_snap) {
            /* Download fields */
            gpu_download(h_phi, h_theta, N3);

            if (need_snap) {
                for (int a = 0; a < 3; a++)
                    gpu_download_f32(d_phi[a], h_f32_bufs[a], N3);
                for (int a = 0; a < 3; a++)
                    gpu_download_f32(d_theta[a], h_f32_bufs[3+a], N3);
                cudaDeviceSynchronize();
                void *cols[] = {h_f32_bufs[0], h_f32_bufs[1], h_f32_bufs[2],
                                h_f32_bufs[3], h_f32_bufs[4], h_f32_bufs[5]};
                sfa_write_frame(sfa, t, cols);
            }

            if (need_analysis) {
                /* Full energy decomposition (needs velocities) */
                gpu_download_vel(h_vel_phi, h_vel_theta, N3);

                double epk, etk, eg, em, ep, etg, etm, ec;
                compute_full_energy(h_phi, h_theta, h_vel_phi, h_vel_theta, N, dx,
                                    &epk, &etk, &eg, &em, &ep, &etg, &etm, &ec);
                double E_total = epk+etk+eg+em+ep+etg+etm+ec;

                double winding = compute_winding(h_phi, N, dx, L);
                double z_cen = braid_z_centroid(h_phi, N, dx, L);

                fprintf(ts_fp, "%.4f\t%.4e\t%.4e\t%.4e\t%.4e\t%.4e\t%.4e\t%.4e\t%.4e\t%.4e\t%.4f\t%.4f\n",
                        t, epk, etk, eg, em, ep, etg, etm, ec, E_total, winding, z_cen);
                fflush(ts_fp);

                /* Radial theta profile at z=0 */
                double r_cen[N_RADIAL_BINS], trms_prof[N_RADIAL_BINS], tmag_prof[N_RADIAL_BINS];
                compute_radial_theta_profile(h_theta, N, dx, L, 0.0,
                                             N_RADIAL_BINS, R_max_profile,
                                             r_cen, trms_prof, tmag_prof);

                /* Find theta_rms at r=1 for Biot-Savart normalization */
                double theta_at_r1 = 0;
                for (int b = 0; b < N_RADIAL_BINS; b++) {
                    if (r_cen[b] > 1.0) { theta_at_r1 = tmag_prof[b]; break; }
                }

                for (int b = 0; b < N_RADIAL_BINS; b++) {
                    double biot_savart = (r_cen[b] > 0.1) ? theta_at_r1 * 1.0 / r_cen[b] : 0;
                    fprintf(prof_fp, "%.4f\t%.4f\t%.6e\t%.6e\t%.6e\n",
                            t, r_cen[b], trms_prof[b], tmag_prof[b], biot_savart);
                }
                fflush(prof_fp);

                /* Console output at major intervals */
                if (step % (analysis_every * 100) == 0) {
                    double v_braid = (t > prev_t) ? (z_cen - prev_z_cen)/(t - prev_t) : 0;
                    printf("t=%6.1f E=%.3e Ep=%.1f winding=%.3f z_cen=%.2f v_z=%.4f [%.0f%%]\n",
                           t, E_total, ep, winding, z_cen, v_braid, 100.0*step/n_steps);
                    fflush(stdout);
                    prev_z_cen = z_cen;
                    prev_t = t;
                }
            }
        }
    }

    cudaEventRecord(stop);
    cudaEventSynchronize(stop);
    float ms;
    cudaEventElapsedTime(&ms, start, stop);
    printf("\n=== Complete: %.1f seconds (%.2f ms/step) ===\n", ms/1000.0, ms/n_steps);

    fclose(ts_fp);
    fclose(prof_fp);
    uint32_t nframes = sfa->total_frames;
    sfa_close(sfa);
    printf("SFA: %s (%u frames)\n", sfapath, nframes);
    printf("Profiles: %s\n", profpath);

    /* Cleanup */
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
