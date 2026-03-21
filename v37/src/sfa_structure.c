/*  sfa_structure.c — Deep structural analysis of SFA field archives
 *
 *  Reads SFA files and computes:
 *    1. Radial profile: field amplitude, P, energy density vs r
 *    2. Phase coherence: are the three fields phase-locked?
 *    3. Binding concentration: how localized is V(P)?
 *    4. Mode analysis: what fraction of energy is in the core vs radiation?
 *    5. Gradient energy: ratio of gradient to mass energy (stability indicator)
 *    6. Time evolution of all metrics
 *
 *  Usage: ./sfa_structure <file.sfa> [-every N] [-o output.tsv]
 *
 *  Build: gcc -O3 -march=native -fopenmp -o sfa_structure src/sfa_structure.c -lzstd -lm
 */

#define SFA_IMPLEMENTATION
#include "../../sfa/format/sfa.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <omp.h>

#define PI 3.14159265358979323846
#define MU    -41.345
#define KAPPA  50.0

#define N_SHELLS 40
#define N_PHASE_BINS 64

typedef struct {
    /* Per-shell (radial profile) */
    double rho_shell[N_SHELLS];     /* field amplitude^2 */
    double P_shell[N_SHELLS];       /* triple product */
    double Vpot_shell[N_SHELLS];    /* V(P) potential energy density */
    double grad_shell[N_SHELLS];    /* gradient energy density */
    double vol_shell[N_SHELLS];     /* shell volume (for normalization) */

    /* Global metrics */
    double E_pot_total;
    double E_grad_total;
    double E_mass_total;
    double rho_total;               /* total field amplitude */
    double P_int;                   /* integrated |P| */

    /* Core vs radiation */
    double E_pot_core;              /* E_pot within R_core */
    double E_pot_outer;             /* E_pot outside R_core */
    double rho_core;                /* amplitude in core */
    double rho_outer;               /* amplitude outside */
    double R_core;                  /* core radius (half-energy radius) */

    /* Phase coherence */
    double phase_coherence;         /* 1 = perfectly locked, 0 = random */
    double phase_01;                /* mean phase difference phi0-phi1 */
    double phase_12;                /* mean phase difference phi1-phi2 */

    /* Compactness */
    double R_rms;                   /* RMS radius of field energy */
    double aspect;                  /* max/min inertia eigenvalue */
    double concentration;           /* E_pot_core / E_pot_total */
} FrameAnalysis;

static void analyze_frame(double *phi[3], int N, double dx, double L,
                          FrameAnalysis *out) {
    memset(out, 0, sizeof(FrameAnalysis));
    const int NN = N * N;
    const long N3 = (long)N * N * N;
    const double dV = dx * dx * dx;
    const double idx1 = 1.0 / (2.0 * dx);
    const double R_max = L;
    const double dr = R_max / N_SHELLS;
    const double MASS2 = 2.25;

    /* Pass 1: radial profiles, energies, phase coherence */
    double epot = 0, egrad = 0, emass = 0, rho_tot = 0, p_int = 0;
    double phase_cos_01 = 0, phase_sin_01 = 0;
    double phase_cos_12 = 0, phase_sin_12 = 0;
    double phase_weight = 0;
    double cm_x = 0, cm_y = 0, cm_z = 0, cm_w = 0;

    double rho_shell[N_SHELLS], P_shell[N_SHELLS], Vpot_shell[N_SHELLS];
    double grad_shell[N_SHELLS], vol_shell[N_SHELLS];
    memset(rho_shell, 0, sizeof(rho_shell));
    memset(P_shell, 0, sizeof(P_shell));
    memset(Vpot_shell, 0, sizeof(Vpot_shell));
    memset(grad_shell, 0, sizeof(grad_shell));
    memset(vol_shell, 0, sizeof(vol_shell));

    #pragma omp parallel for reduction(+:epot,egrad,emass,rho_tot,p_int, \
        phase_cos_01,phase_sin_01,phase_cos_12,phase_sin_12,phase_weight, \
        cm_x,cm_y,cm_z,cm_w) schedule(static)
    for (long idx = 0; idx < N3; idx++) {
        int i = (int)(idx / NN), j = (int)((idx / N) % N), k = (int)(idx % N);
        double x = -L + i * dx, y = -L + j * dx, z = -L + k * dx;
        double r = sqrt(x*x + y*y + z*z);

        int ip = (i+1)%N, im = (i-1+N)%N;
        int jp = (j+1)%N, jm = (j-1+N)%N;
        int kp = (k+1)%N, km = (k-1+N)%N;
        long n_ip = (long)ip*NN+j*N+k, n_im = (long)im*NN+j*N+k;
        long n_jp = (long)i*NN+jp*N+k, n_jm = (long)i*NN+jm*N+k;
        long n_kp = (long)i*NN+j*N+kp, n_km = (long)i*NN+j*N+km;

        double p0 = phi[0][idx], p1 = phi[1][idx], p2 = phi[2][idx];
        double rho = p0*p0 + p1*p1 + p2*p2;
        double P = p0 * p1 * p2;
        double V = (MU/2.0)*P*P/(1.0+KAPPA*P*P);

        /* Gradient energy */
        double eg = 0;
        for (int a = 0; a < 3; a++) {
            double gx = (phi[a][n_ip] - phi[a][n_im]) * idx1;
            double gy = (phi[a][n_jp] - phi[a][n_jm]) * idx1;
            double gz = (phi[a][n_kp] - phi[a][n_km]) * idx1;
            eg += 0.5 * (gx*gx + gy*gy + gz*gz);
        }

        epot += V * dV;
        egrad += eg * dV;
        emass += 0.5 * MASS2 * rho * dV;
        rho_tot += rho * dV;
        p_int += fabs(P) * dV;

        /* Centroid */
        cm_x += x * rho * dV;
        cm_y += y * rho * dV;
        cm_z += z * rho * dV;
        cm_w += rho * dV;

        /* Phase coherence: measure angle between field pairs */
        double amp01 = sqrt(p0*p0 + p1*p1);
        if (amp01 > 0.01) {
            double ang01 = atan2(p1, p0);
            phase_cos_01 += cos(ang01) * amp01 * dV;
            phase_sin_01 += sin(ang01) * amp01 * dV;
            phase_weight += amp01 * dV;
        }
        double amp12 = sqrt(p1*p1 + p2*p2);
        if (amp12 > 0.01) {
            double ang12 = atan2(p2, p1);
            phase_cos_12 += cos(ang12) * amp12 * dV;
            phase_sin_12 += sin(ang12) * amp12 * dV;
        }

        /* Radial shell */
        int shell = (int)(r / dr);
        if (shell < N_SHELLS) {
            #pragma omp atomic
            rho_shell[shell] += rho * dV;
            #pragma omp atomic
            P_shell[shell] += P * dV;
            #pragma omp atomic
            Vpot_shell[shell] += V * dV;
            #pragma omp atomic
            grad_shell[shell] += eg * dV;
            #pragma omp atomic
            vol_shell[shell] += dV;
        }
    }

    /* Store results */
    out->E_pot_total = epot;
    out->E_grad_total = egrad;
    out->E_mass_total = emass;
    out->rho_total = rho_tot;
    out->P_int = p_int;
    memcpy(out->rho_shell, rho_shell, sizeof(rho_shell));
    memcpy(out->P_shell, P_shell, sizeof(P_shell));
    memcpy(out->Vpot_shell, Vpot_shell, sizeof(Vpot_shell));
    memcpy(out->grad_shell, grad_shell, sizeof(grad_shell));
    memcpy(out->vol_shell, vol_shell, sizeof(vol_shell));

    /* Phase coherence */
    if (phase_weight > 1e-10) {
        double pc01 = sqrt(phase_cos_01*phase_cos_01 + phase_sin_01*phase_sin_01) / phase_weight;
        double pc12 = sqrt(phase_cos_12*phase_cos_12 + phase_sin_12*phase_sin_12) / phase_weight;
        out->phase_coherence = (pc01 + pc12) / 2.0;
        out->phase_01 = atan2(phase_sin_01, phase_cos_01);
        out->phase_12 = atan2(phase_sin_12, phase_cos_12);
    }

    /* RMS radius */
    double cx = (cm_w > 0) ? cm_x/cm_w : 0;
    double cy = (cm_w > 0) ? cm_y/cm_w : 0;
    double cz = (cm_w > 0) ? cm_z/cm_w : 0;
    double r2_sum = 0, w_sum = 0;
    #pragma omp parallel for reduction(+:r2_sum,w_sum)
    for (long idx = 0; idx < N3; idx++) {
        int i = (int)(idx/NN), j = (int)((idx/N)%N), k = (int)(idx%N);
        double x = -L+i*dx-cx, y = -L+j*dx-cy, z = -L+k*dx-cz;
        double rho = phi[0][idx]*phi[0][idx] + phi[1][idx]*phi[1][idx] + phi[2][idx]*phi[2][idx];
        r2_sum += (x*x+y*y+z*z) * rho * dV;
        w_sum += rho * dV;
    }
    out->R_rms = (w_sum > 0) ? sqrt(r2_sum / w_sum) : 0;

    /* Core radius: find R where 50% of |V(P)| is enclosed */
    double vpot_cumul = 0, vpot_half = fabs(epot) * 0.5;
    out->R_core = R_max;
    for (int s = 0; s < N_SHELLS; s++) {
        vpot_cumul += fabs(Vpot_shell[s]);
        if (vpot_cumul >= vpot_half && out->R_core >= R_max) {
            out->R_core = (s + 1) * dr;
        }
    }

    /* Core vs outer energy */
    out->E_pot_core = 0;
    out->rho_core = 0;
    int core_shells = (int)(out->R_core / dr);
    for (int s = 0; s < N_SHELLS; s++) {
        if (s < core_shells) {
            out->E_pot_core += Vpot_shell[s];
            out->rho_core += rho_shell[s];
        } else {
            out->E_pot_outer += Vpot_shell[s];
            out->rho_outer += rho_shell[s];
        }
    }
    out->concentration = (fabs(epot) > 1e-10) ?
        fabs(out->E_pot_core) / fabs(epot) : 0;
}

int main(int argc, char **argv) {
    if (argc < 2) {
        fprintf(stderr, "Usage: sfa_structure <file.sfa> [-every N] [-o file.tsv]\n");
        return 1;
    }

    const char *path = argv[1];
    int every = 5;
    const char *outpath = NULL;

    for (int i = 2; i < argc; i++) {
        if (!strcmp(argv[i], "-every")) every = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-o")) outpath = argv[++i];
    }

    SFA *sfa = sfa_open(path);
    if (!sfa) { fprintf(stderr, "Failed to open %s\n", path); return 1; }

    int N = sfa->Nx;
    double L = sfa->Lx;
    double dx = 2.0 * L / (N - 1);
    long N3 = (long)N * N * N;

    fprintf(stderr, "SFA: %s\n", path);
    fprintf(stderr, "  Grid: %d³, L=%.1f, dx=%.4f\n", N, L, dx);
    fprintf(stderr, "  Frames: %d, every=%d\n", sfa->total_frames, every);

    /* Determine dtype size */
    int elem_size = sfa_dtype_size[sfa->columns[0].dtype];
    int is_f32 = (sfa->columns[0].dtype == SFA_F32);
    fprintf(stderr, "  Column dtype: %s\n", is_f32 ? "f32" : "f64");

    /* Allocate frame buffer and field arrays */
    uint8_t *frame_buf = (uint8_t*)malloc(sfa->frame_bytes);
    double *phi[3];
    for (int a = 0; a < 3; a++)
        phi[a] = (double*)malloc(N3 * sizeof(double));

    /* Output */
    FILE *fp = outpath ? fopen(outpath, "w") : stdout;
    fprintf(fp, "frame\ttime\tE_pot\tE_grad\tE_mass\tP_int\tR_rms\tR_core\t"
                "phase_coh\tphase_01\tphase_12\tconcentration\t"
                "grad_mass_ratio\trho_core_frac\n");

    /* Radial profile header (to stderr) */
    fprintf(stderr, "\n=== Radial profiles ===\n");

    int n_analyzed = 0;
    for (uint32_t f = 0; f < sfa->total_frames; f += every) {
        double t = sfa_frame_time(sfa, f);
        int rc = sfa_read_frame(sfa, f, frame_buf);
        if (rc < 0) {
            fprintf(stderr, "Warning: failed to read frame %u (rc=%d), skipping\n", f, rc);
            continue;
        }

        /* Extract phi_0, phi_1, phi_2 from frame buffer */
        uint64_t col_bytes = N3 * elem_size;
        for (int a = 0; a < 3; a++) {
            if (is_f32) {
                float *src = (float*)(frame_buf + a * col_bytes);
                for (long i = 0; i < N3; i++)
                    phi[a][i] = (double)src[i];
            } else {
                double *src = (double*)(frame_buf + a * col_bytes);
                memcpy(phi[a], src, N3 * sizeof(double));
            }
        }

        FrameAnalysis fa;
        analyze_frame(phi, N, dx, L, &fa);

        double grad_mass = (fa.E_mass_total > 1e-10) ?
            fa.E_grad_total / fa.E_mass_total : 0;
        double rho_core_frac = (fa.rho_total > 1e-10) ?
            fa.rho_core / fa.rho_total : 0;

        fprintf(fp, "%u\t%.2f\t%.4e\t%.4e\t%.4e\t%.4e\t%.3f\t%.3f\t"
                    "%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n",
                f, t, fa.E_pot_total, fa.E_grad_total, fa.E_mass_total,
                fa.P_int, fa.R_rms, fa.R_core, fa.phase_coherence,
                fa.phase_01, fa.phase_12, fa.concentration,
                grad_mass, rho_core_frac);

        /* Print radial profile for first and last frames */
        if (n_analyzed == 0 || f + every >= sfa->total_frames) {
            double dr = L / N_SHELLS;
            fprintf(stderr, "\nFrame %u (t=%.1f):\n", f, t);
            fprintf(stderr, "  r_mid    rho_dens    P_dens    Vpot_dens   grad_dens\n");
            for (int s = 0; s < N_SHELLS; s++) {
                if (fa.vol_shell[s] < 1e-20) continue;
                double rm = (s + 0.5) * dr;
                fprintf(stderr, "  %5.2f  %10.4e  %10.4e  %10.4e  %10.4e\n",
                    rm, fa.rho_shell[s]/fa.vol_shell[s],
                    fa.P_shell[s]/fa.vol_shell[s],
                    fa.Vpot_shell[s]/fa.vol_shell[s],
                    fa.grad_shell[s]/fa.vol_shell[s]);
            }
            fprintf(stderr, "  R_rms=%.3f  R_core=%.3f  phase_coh=%.4f  "
                    "concentration=%.4f  grad/mass=%.4f\n",
                    fa.R_rms, fa.R_core, fa.phase_coherence,
                    fa.concentration, grad_mass);
        }

        n_analyzed++;
    }

    if (outpath) fclose(fp);

    /* Summary comparison metrics */
    fprintf(stderr, "\n=== Structure Summary ===\n");
    fprintf(stderr, "  Frames analyzed: %d\n", n_analyzed);

    for (int a = 0; a < 3; a++) free(phi[a]);
    free(frame_buf);
    sfa_close(sfa);
    return 0;
}
