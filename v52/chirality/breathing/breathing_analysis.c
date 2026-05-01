/*  breathing_analysis.c — Track intra-particle dynamics during breathing cycle.
 *
 *  For each frame: find centroid, compute radial shell diagnostics.
 *  Tracks φ and θ independently to see if θ persists when φ disperses.
 *
 *  Outputs two TSV files:
 *    <base>_global.tsv — per-frame global metrics
 *    <base>_radial.tsv — per-frame per-shell radial profiles
 *
 *  Build: gcc -O3 -fopenmp -o breathing_analysis breathing_analysis.c -lzstd -lm
 *  Usage: breathing_analysis input.sfa output_base [t_start] [t_end] [frame_stride]
 */
#define SFA_IMPLEMENTATION
#include "../../sfa/format/sfa.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define MAX_SHELLS 20
#define SHELL_DR 1.0  /* code units per shell */

int main(int argc, char **argv) {
    if (argc < 3) {
        fprintf(stderr, "Usage: %s input.sfa output_base [t_start] [t_end] [stride]\n", argv[0]);
        return 1;
    }
    const char *inpath = argv[1];
    const char *base = argv[2];
    double t_start = argc > 3 ? atof(argv[3]) : 0;
    double t_end = argc > 4 ? atof(argv[4]) : 1e9;
    int stride = argc > 5 ? atoi(argv[5]) : 1;

    SFA *s = sfa_open(inpath);
    if (!s) { fprintf(stderr, "Cannot open %s\n", inpath); return 1; }

    /* Fix total_frames from JTOP chain */
    {
        uint64_t jt = s->first_jtop_offset;
        uint32_t total = 0;
        while (jt) {
            fseek(s->fp, (long)jt + 12, SEEK_SET);
            uint32_t mx, cur; uint64_t nxt;
            fread(&mx, 4, 1, s->fp); fread(&cur, 4, 1, s->fp); fread(&nxt, 8, 1, s->fp);
            for (uint32_t i = 0; i < cur; i++) {
                uint64_t jo; uint32_t ff, fc;
                fread(&jo, 8, 1, s->fp); fread(&ff, 4, 1, s->fp); fread(&fc, 4, 1, s->fp);
                total += fc;
            }
            jt = nxt;
        }
        if (total > s->total_frames) s->total_frames = total;
    }

    int N = s->Nx;
    long N3 = (long)N*N*N;
    int NN = N*N;
    double L = s->Lx;
    double dx = 2.0 * L / (N - 1);
    double dV = dx*dx*dx;
    double idx1 = 1.0 / (2.0 * dx);
    double idx2 = 1.0 / (dx * dx);
    double mu = -41.345, kappa = 50.0, m2 = 2.25, eta = 0.5;

    printf("Breathing analysis: %s\n", inpath);
    printf("  Grid: %d³, L=%.1f, %d frames, stride=%d\n", N, L, s->total_frames, stride);

    /* Open output files */
    char gpath[512], rpath[512];
    snprintf(gpath, 512, "%s_global.tsv", base);
    snprintf(rpath, 512, "%s_radial.tsv", base);
    FILE *gf = fopen(gpath, "w");
    FILE *rf = fopen(rpath, "w");

    fprintf(gf, "frame\tt\tphi_max\tP_int\tP_peak\tphi_rms_core\ttheta_rms_core\t"
                "curl_phi_rms_core\tmismatch_rms_core\tE_kin_core\tE_pot_core\t"
                "phi_rms_all\ttheta_rms_all\tcx\tcy\tcz\n");
    fprintf(rf, "frame\tt\tshell\tr_inner\tr_outer\t"
                "phi_rms\ttheta_rms\tP_mean\tcurl_rms\tmismatch_rms\t"
                "E_kin\tE_mass\tE_pot\tnvox\n");

    void *buf = malloc(s->frame_bytes);

    /* Allocate field arrays */
    double *phi[3], *theta[3], *vel[3];
    for (int a = 0; a < 3; a++) {
        phi[a] = (double*)malloc(N3 * sizeof(double));
        theta[a] = (double*)malloc(N3 * sizeof(double));
        vel[a] = (double*)malloc(N3 * sizeof(double));
    }

    int n_analyzed = 0;
    for (uint32_t fi = 0; fi < s->total_frames; fi += stride) {
        SFA_L2Entry entry;
        if (sfa_find_frame(s, fi, &entry) < 0) continue;
        if (entry.time < t_start || entry.time > t_end) continue;

        if (sfa_read_frame(s, fi, buf) < 0) continue;

        /* Extract fields from frame buffer */
        uint64_t col_off[12] = {0};
        for (int c = 1; c < (int)s->n_columns && c < 12; c++)
            col_off[c] = col_off[c-1] + N3 * sfa_dtype_size[s->columns[c-1].dtype];

        for (long i = 0; i < N3; i++) {
            for (int a = 0; a < 3; a++) {
                int es = sfa_dtype_size[s->columns[a].dtype];
                if (es == 2) {
                    uint16_t h; memcpy(&h, (uint8_t*)buf + col_off[a] + i*2, 2);
                    phi[a][i] = sfa_f16_to_f32(h);
                } else {
                    float f; memcpy(&f, (uint8_t*)buf + col_off[a] + i*4, 4);
                    phi[a][i] = f;
                }
                if (a+3 < (int)s->n_columns) {
                    if (es == 2) {
                        uint16_t h; memcpy(&h, (uint8_t*)buf + col_off[a+3] + i*2, 2);
                        theta[a][i] = sfa_f16_to_f32(h);
                    } else {
                        float f; memcpy(&f, (uint8_t*)buf + col_off[a+3] + i*4, 4);
                        theta[a][i] = f;
                    }
                }
                if (a+6 < (int)s->n_columns) {
                    if (es == 2) {
                        uint16_t h; memcpy(&h, (uint8_t*)buf + col_off[a+6] + i*2, 2);
                        vel[a][i] = sfa_f16_to_f32(h);
                    } else {
                        float f; memcpy(&f, (uint8_t*)buf + col_off[a+6] + i*4, 4);
                        vel[a][i] = f;
                    }
                }
            }
        }

        /* Find centroid (P-weighted) */
        double cx=0, cy=0, cz=0, wsum=0;
        double P_peak=0, phi_max=0, P_int=0;
        for (long idx = 0; idx < N3; idx++) {
            int i = (int)(idx/NN), j = (int)((idx/N)%N), k = (int)(idx%N);
            double P = fabs(phi[0][idx]*phi[1][idx]*phi[2][idx]);
            double ps = phi[0][idx]*phi[0][idx]+phi[1][idx]*phi[1][idx]+phi[2][idx]*phi[2][idx];
            if (sqrt(ps) > phi_max) phi_max = sqrt(ps);
            if (P > P_peak) P_peak = P;
            P_int += P * dV;
            double x = -L + i*dx, y = -L + j*dx, z = -L + k*dx;
            cx += P*x; cy += P*y; cz += P*z; wsum += P;
        }
        if (wsum > 0) { cx /= wsum; cy /= wsum; cz /= wsum; }

        /* Radial shell analysis */
        double shell_phi_sq[MAX_SHELLS] = {0};
        double shell_theta_sq[MAX_SHELLS] = {0};
        double shell_P[MAX_SHELLS] = {0};
        double shell_curl_sq[MAX_SHELLS] = {0};
        double shell_mismatch_sq[MAX_SHELLS] = {0};
        double shell_Ekin[MAX_SHELLS] = {0};
        double shell_Emass[MAX_SHELLS] = {0};
        double shell_Epot[MAX_SHELLS] = {0};
        int shell_nvox[MAX_SHELLS] = {0};

        /* Also track core (r < 3) and global */
        double phi_rms_core=0, theta_rms_core=0, curl_core=0, mis_core=0;
        double Ekin_core=0, Epot_core=0;
        int nvox_core=0;
        double phi_rms_all=0, theta_rms_all=0;

        for (long idx = 0; idx < N3; idx++) {
            int i = (int)(idx/NN), j = (int)((idx/N)%N), k = (int)(idx%N);
            double x = -L + i*dx, y = -L + j*dx, z = -L + k*dx;
            double r = sqrt((x-cx)*(x-cx)+(y-cy)*(y-cy)+(z-cz)*(z-cz));
            int shell = (int)(r / SHELL_DR);
            if (shell >= MAX_SHELLS) shell = MAX_SHELLS - 1;

            double ps = phi[0][idx]*phi[0][idx]+phi[1][idx]*phi[1][idx]+phi[2][idx]*phi[2][idx];
            double ts = theta[0][idx]*theta[0][idx]+theta[1][idx]*theta[1][idx]+theta[2][idx]*theta[2][idx];
            double P = fabs(phi[0][idx]*phi[1][idx]*phi[2][idx]);
            double Vpot = (mu/2.0)*P*P/(1.0+kappa*P*P);
            double vs = vel[0][idx]*vel[0][idx]+vel[1][idx]*vel[1][idx]+vel[2][idx]*vel[2][idx];

            /* Curl(phi) — periodic BC */
            int ip=(i+1)%N, im=(i-1+N)%N, jp=(j+1)%N, jm=(j-1+N)%N, kp=(k+1)%N, km=(k-1+N)%N;
            double curlx = (phi[2][(long)i*NN+jp*N+k]-phi[2][(long)i*NN+jm*N+k])*idx1
                         - (phi[1][(long)i*NN+j*N+kp]-phi[1][(long)i*NN+j*N+km])*idx1;
            double curly = (phi[0][(long)i*NN+j*N+kp]-phi[0][(long)i*NN+j*N+km])*idx1
                         - (phi[2][(long)ip*NN+j*N+k]-phi[2][(long)im*NN+j*N+k])*idx1;
            double curlz = (phi[1][(long)ip*NN+j*N+k]-phi[1][(long)im*NN+j*N+k])*idx1
                         - (phi[0][(long)i*NN+jp*N+k]-phi[0][(long)i*NN+jm*N+k])*idx1;
            double curl_sq = curlx*curlx + curly*curly + curlz*curlz;

            /* Mismatch: M = curl(phi)/2 - theta */
            double mx = curlx/2 - theta[0][idx];
            double my = curly/2 - theta[1][idx];
            double mz = curlz/2 - theta[2][idx];
            double mis_sq = mx*mx + my*my + mz*mz;

            shell_phi_sq[shell] += ps * dV;
            shell_theta_sq[shell] += ts * dV;
            shell_P[shell] += P * dV;
            shell_curl_sq[shell] += curl_sq * dV;
            shell_mismatch_sq[shell] += mis_sq * dV;
            shell_Ekin[shell] += 0.5 * vs * dV;
            shell_Emass[shell] += 0.5 * m2 * ps * dV;
            shell_Epot[shell] += Vpot * dV;
            shell_nvox[shell]++;

            phi_rms_all += ps;
            theta_rms_all += ts;

            if (r < 3.0) {
                phi_rms_core += ps; theta_rms_core += ts;
                curl_core += curl_sq; mis_core += mis_sq;
                Ekin_core += 0.5 * vs * dV; Epot_core += Vpot * dV;
                nvox_core++;
            }
        }

        phi_rms_all = sqrt(phi_rms_all / N3);
        theta_rms_all = sqrt(theta_rms_all / N3);
        if (nvox_core > 0) {
            phi_rms_core = sqrt(phi_rms_core / nvox_core);
            theta_rms_core = sqrt(theta_rms_core / nvox_core);
            curl_core = sqrt(curl_core / nvox_core);
            mis_core = sqrt(mis_core / nvox_core);
        }

        /* Write global */
        fprintf(gf, "%d\t%.4f\t%.6f\t%.4f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.4f\t%.4f\t%.6f\t%.6f\t%.3f\t%.3f\t%.3f\n",
                fi, entry.time, phi_max, P_int, P_peak,
                phi_rms_core, theta_rms_core, curl_core, mis_core,
                Ekin_core, Epot_core, phi_rms_all, theta_rms_all,
                cx, cy, cz);

        /* Write radial */
        for (int sh = 0; sh < MAX_SHELLS; sh++) {
            if (shell_nvox[sh] == 0) continue;
            double r_in = sh * SHELL_DR, r_out = (sh+1) * SHELL_DR;
            int nv = shell_nvox[sh];
            fprintf(rf, "%d\t%.4f\t%d\t%.1f\t%.1f\t%.6f\t%.6f\t%.6e\t%.6f\t%.6f\t%.4f\t%.4f\t%.4f\t%d\n",
                    fi, entry.time, sh, r_in, r_out,
                    sqrt(shell_phi_sq[sh] / (nv * dV)),
                    sqrt(shell_theta_sq[sh] / (nv * dV)),
                    shell_P[sh] / (nv * dV),
                    sqrt(shell_curl_sq[sh] / (nv * dV)),
                    sqrt(shell_mismatch_sq[sh] / (nv * dV)),
                    shell_Ekin[sh], shell_Emass[sh], shell_Epot[sh], nv);
        }

        n_analyzed++;
        if (n_analyzed % 50 == 0)
            printf("  frame %d/%d (t=%.1f) phi_max=%.3f P_int=%.1f θ_core=%.4f\n",
                   fi, s->total_frames, entry.time, phi_max, P_int, theta_rms_core);
    }

    fclose(gf); fclose(rf);
    printf("\nDone: %d frames analyzed\n", n_analyzed);
    printf("  %s\n  %s\n", gpath, rpath);

    for (int a = 0; a < 3; a++) { free(phi[a]); free(theta[a]); free(vel[a]); }
    free(buf);
    sfa_close(s);
    return 0;
}
