/*  phase_binding.c — Characterize binding structure, phi phases, density distribution
 *
 *  For each frame of the proton steep gradient test SFA, computes:
 *    1. Radial phi_rms, phi_max, E_phi, dominant phi component
 *    2. Binding density |P| = |phi_x*phi_y*phi_z| profiles and volumes
 *    3. Phase structure (sign patterns at centroid and r=5)
 *    4. Density asymmetry (left/right split for gravity diagnostic)
 *    5. Theta characterization
 *    6. Velocity/momentum profiles
 *
 *  Build: gcc -O3 -fopenmp -o phase_binding phase_binding.c -lzstd -lm
 */

#define SFA_IMPLEMENTATION
#include "../../sfa/format/sfa.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <omp.h>

#define PI 3.14159265358979323846
#define NSHELLS 30
#define DR 0.5
#define R_MAX 15.0

/* --- f16 conversion --- */
static inline double f16_to_f64(uint16_t h) {
    uint16_t sign = h & 0x8000;
    int exp = (h >> 10) & 0x1F;
    uint16_t mant = h & 0x3FF;
    if (exp == 0) return 0.0;
    if (exp == 31) return sign ? -1e30 : 1e30;
    float fv;
    uint32_t x = ((uint32_t)sign << 16) | ((uint32_t)(exp - 15 + 127) << 23) | ((uint32_t)mant << 13);
    memcpy(&fv, &x, 4);
    return (double)fv;
}

/* Extract a column by semantic+component from frame buffer */
static double *extract_col(void *buf, SFA *sfa, int sem, int comp, long N3) {
    double *arr = (double *)calloc(N3, sizeof(double));
    uint64_t off = 0;
    for (uint32_t c = 0; c < sfa->n_columns; c++) {
        int dt = sfa->columns[c].dtype;
        int s = sfa->columns[c].semantic;
        int cp = sfa->columns[c].component;
        int es = sfa_dtype_size[dt];
        if (s == sem && cp == comp) {
            uint8_t *src = (uint8_t *)buf + off;
            if (dt == SFA_F64) memcpy(arr, src, N3 * 8);
            else if (dt == SFA_F32) for (long i = 0; i < N3; i++) arr[i] = (double)((float *)src)[i];
            else if (dt == SFA_F16) for (long i = 0; i < N3; i++) arr[i] = f16_to_f64(((uint16_t *)src)[i]);
            return arr;
        }
        off += (uint64_t)N3 * es;
    }
    return arr; /* zeros if not found */
}

/* --- Per-shell data --- */
typedef struct {
    double phi_rms;
    double P_mean;      /* |phi_x * phi_y * phi_z| */
    double theta_rms;
    double phi_v_rms;
    double v_radial;    /* mean radial velocity (+ = expanding) */
    long count;
} ShellData;

/* --- Per-frame summary --- */
typedef struct {
    double t;
    double cx, cy, cz;         /* |P|-weighted centroid */
    double phi_rms_core;       /* r<2 */
    double phi_rms_r5;         /* 4.5<r<5.5 */
    double phi_max;
    double phi_max_rx, phi_max_ry, phi_max_rz;
    double P_total, P_max, P_max_r;
    long binding_vol_001, binding_vol_01, binding_vol_1;
    long binding_area;         /* surface voxels crossing 0.001 */
    double E_phi, E_theta, E_kin;
    double asym_phi, asym_P;   /* (left-right)/(left+right) */
    double theta_rms_core, theta_rms_r5, theta_total;
    double v_radial_core, v_radial_r5;
    /* Phase sign at centroid */
    int sign_cx, sign_cy, sign_cz;
    /* Phase sign counts (8 octants) */
    long phase_counts[8];
    /* Dominant component counts at core (r<5) */
    long dom_x, dom_y, dom_z;

    ShellData shells[NSHELLS];
} FrameResult;

int main(int argc, char **argv) {
    const char *sfa_path = "/space/scp/v43/gradient_test/results/proton_steep_output.sfa";
    if (argc >= 2) sfa_path = argv[1];

    SFA *sfa = sfa_open(sfa_path);
    if (!sfa) { fprintf(stderr, "Cannot open %s\n", sfa_path); return 1; }

    int N = sfa->Nx;
    double L = sfa->Lx;
    long N3 = (long)N * N * N;
    int NN = N * N;
    double dx = 2.0 * L / (N - 1);
    double dV = dx * dx * dx;
    int nf = sfa->total_frames;

    printf("=== Phase Binding Analysis ===\n");
    printf("Grid: %d^3, L=%.1f, dx=%.4f, frames=%d\n\n", N, L, dx, nf);

    /* Allocate frame buffer */
    void *buf = malloc(sfa->frame_bytes);
    if (!buf) { fprintf(stderr, "Cannot alloc %.1f GB\n", sfa->frame_bytes / 1e9); return 1; }

    FrameResult *results = (FrameResult *)calloc(nf, sizeof(FrameResult));

    /* Open output files */
    FILE *fp_summary = fopen("phase_binding.tsv", "w");
    FILE *fp_profiles = fopen("phase_profiles.tsv", "w");

    fprintf(fp_summary, "t\tcx\tcy\tcz\tphi_rms_core\tphi_rms_r5\tphi_max\t"
            "P_total\tP_max\tP_max_r\tbinding_vol_001\tbinding_vol_01\tbinding_vol_1\t"
            "binding_area\tE_phi\tE_theta\tE_kin\tasym_phi\tasym_P\t"
            "theta_rms_core\ttheta_rms_r5\ttheta_total\tv_radial_core\tv_radial_r5\n");

    fprintf(fp_profiles, "frame\tr\tphi_rms\tP_mean\ttheta_rms\tphi_v_rms\tv_radial\n");

    for (int fi = 0; fi < nf; fi++) {
        double t0 = omp_get_wtime();
        double time = sfa_frame_time(sfa, fi);

        printf("Frame %d (t=%.1f) ... ", fi, time);
        fflush(stdout);

        sfa_read_frame(sfa, fi, buf);
        double t_read = omp_get_wtime();
        printf("read %.1fs, ", t_read - t0);
        fflush(stdout);

        /* Extract all 12 columns */
        double *phi_x  = extract_col(buf, sfa, SFA_POSITION, 0, N3);
        double *phi_y  = extract_col(buf, sfa, SFA_POSITION, 1, N3);
        double *phi_z  = extract_col(buf, sfa, SFA_POSITION, 2, N3);
        double *th_x   = extract_col(buf, sfa, SFA_ANGLE, 0, N3);
        double *th_y   = extract_col(buf, sfa, SFA_ANGLE, 1, N3);
        double *th_z   = extract_col(buf, sfa, SFA_ANGLE, 2, N3);
        double *vx     = extract_col(buf, sfa, SFA_VELOCITY, 0, N3);
        double *vy     = extract_col(buf, sfa, SFA_VELOCITY, 1, N3);
        double *vz     = extract_col(buf, sfa, SFA_VELOCITY, 2, N3);
        double *vth_x  = extract_col(buf, sfa, SFA_VELOCITY, 3, N3);
        double *vth_y  = extract_col(buf, sfa, SFA_VELOCITY, 4, N3);
        double *vth_z  = extract_col(buf, sfa, SFA_VELOCITY, 5, N3);

        FrameResult *fr = &results[fi];
        fr->t = time;

        /* --- Pass 1: compute |P|-weighted centroid --- */
        double wcx = 0, wcy = 0, wcz = 0, wtot = 0;
        #pragma omp parallel for reduction(+:wcx,wcy,wcz,wtot) schedule(static)
        for (long idx = 0; idx < N3; idx++) {
            double P = fabs(phi_x[idx] * phi_y[idx] * phi_z[idx]);
            if (P < 1e-6) continue;
            int iz = idx / NN;
            int iy = (idx - (long)iz * NN) / N;
            int ix = idx - (long)iz * NN - iy * N;
            double x = -L + ix * dx;
            double y = -L + iy * dx;
            double z = -L + iz * dx;
            wcx += P * x;
            wcy += P * y;
            wcz += P * z;
            wtot += P;
        }
        if (wtot > 0) {
            fr->cx = wcx / wtot;
            fr->cy = wcy / wtot;
            fr->cz = wcz / wtot;
        }

        /* --- Pass 2: main analysis in proton region (r < R_MAX from centroid) --- */
        /* Thread-local accumulators */
        int nthreads;
        #pragma omp parallel
        { nthreads = omp_get_num_threads(); }

        /* Per-thread shell data */
        ShellData (*tshells)[NSHELLS] = calloc(nthreads, sizeof(ShellData[NSHELLS]));

        /* Per-thread scalars */
        double *t_E_phi    = calloc(nthreads, sizeof(double));
        double *t_E_theta  = calloc(nthreads, sizeof(double));
        double *t_E_kin    = calloc(nthreads, sizeof(double));
        double *t_P_total  = calloc(nthreads, sizeof(double));
        double *t_phi_max  = calloc(nthreads, sizeof(double));
        double *t_phi_max_rx = calloc(nthreads, sizeof(double));
        double *t_phi_max_ry = calloc(nthreads, sizeof(double));
        double *t_phi_max_rz = calloc(nthreads, sizeof(double));
        double *t_P_max    = calloc(nthreads, sizeof(double));
        double *t_P_max_r  = calloc(nthreads, sizeof(double));
        double *t_E_phi_left  = calloc(nthreads, sizeof(double));
        double *t_E_phi_right = calloc(nthreads, sizeof(double));
        double *t_P_left   = calloc(nthreads, sizeof(double));
        double *t_P_right  = calloc(nthreads, sizeof(double));
        double *t_theta_total = calloc(nthreads, sizeof(double));
        long *t_bv001 = calloc(nthreads, sizeof(long));
        long *t_bv01  = calloc(nthreads, sizeof(long));
        long *t_bv1   = calloc(nthreads, sizeof(long));
        long (*t_phase)[8] = calloc(nthreads, sizeof(long[8]));
        long *t_dom_x = calloc(nthreads, sizeof(long));
        long *t_dom_y = calloc(nthreads, sizeof(long));
        long *t_dom_z = calloc(nthreads, sizeof(long));

        #pragma omp parallel
        {
            int tid = omp_get_thread_num();

            #pragma omp for schedule(static)
            for (long idx = 0; idx < N3; idx++) {
                int iz = idx / NN;
                int iy = (idx - (long)iz * NN) / N;
                int ix = idx - (long)iz * NN - iy * N;
                double x = -L + ix * dx;
                double y = -L + iy * dx;
                double z = -L + iz * dx;

                double rx = x - fr->cx;
                double ry = y - fr->cy;
                double rz = z - fr->cz;
                double r = sqrt(rx * rx + ry * ry + rz * rz);
                if (r > R_MAX) continue;

                double px = phi_x[idx], py = phi_y[idx], pz = phi_z[idx];
                double phi2 = px * px + py * py + pz * pz;
                double P = fabs(px * py * pz);
                double th2 = th_x[idx] * th_x[idx] + th_y[idx] * th_y[idx] + th_z[idx] * th_z[idx];
                double v2 = vx[idx] * vx[idx] + vy[idx] * vy[idx] + vz[idx] * vz[idx];

                /* Shell index */
                int si = (int)(r / DR);
                if (si >= NSHELLS) si = NSHELLS - 1;

                tshells[tid][si].phi_rms += phi2;
                tshells[tid][si].P_mean += P;
                tshells[tid][si].theta_rms += th2;
                tshells[tid][si].phi_v_rms += v2;
                /* Radial velocity */
                if (r > 0.1) {
                    double vr = (vx[idx] * rx + vy[idx] * ry + vz[idx] * rz) / r;
                    tshells[tid][si].v_radial += vr;
                }
                tshells[tid][si].count++;

                /* Scalars */
                t_E_phi[tid] += phi2 * dV;
                t_E_theta[tid] += th2 * dV;
                t_E_kin[tid] += v2 * dV;
                t_P_total[tid] += P * dV;
                t_theta_total[tid] += th2 * dV;

                /* phi_max */
                double phi_amp = sqrt(phi2);
                if (phi_amp > t_phi_max[tid]) {
                    t_phi_max[tid] = phi_amp;
                    t_phi_max_rx[tid] = rx;
                    t_phi_max_ry[tid] = ry;
                    t_phi_max_rz[tid] = rz;
                }

                /* P_max */
                if (P > t_P_max[tid]) {
                    t_P_max[tid] = P;
                    t_P_max_r[tid] = r;
                }

                /* Binding volumes */
                if (P > 0.001) t_bv001[tid]++;
                if (P > 0.01)  t_bv01[tid]++;
                if (P > 0.1)   t_bv1[tid]++;

                /* Asymmetry: left (x < cx) vs right (x > cx) */
                if (x < fr->cx) {
                    t_E_phi_left[tid] += phi2 * dV;
                    t_P_left[tid] += P * dV;
                } else {
                    t_E_phi_right[tid] += phi2 * dV;
                    t_P_right[tid] += P * dV;
                }

                /* Phase sign pattern (within r<10) */
                if (r < 10.0) {
                    int si_x = (px >= 0) ? 1 : 0;
                    int si_y = (py >= 0) ? 1 : 0;
                    int si_z = (pz >= 0) ? 1 : 0;
                    int oct = (si_x << 2) | (si_y << 1) | si_z;
                    t_phase[tid][oct]++;
                }

                /* Dominant component at core (r<5) */
                if (r < 5.0) {
                    double ax = fabs(px), ay = fabs(py), az = fabs(pz);
                    if (ax >= ay && ax >= az) t_dom_x[tid]++;
                    else if (ay >= az) t_dom_y[tid]++;
                    else t_dom_z[tid]++;
                }
            }
        } /* end parallel */

        /* Reduce across threads */
        memset(fr->shells, 0, sizeof(fr->shells));
        for (int t = 0; t < nthreads; t++) {
            for (int s = 0; s < NSHELLS; s++) {
                fr->shells[s].phi_rms   += tshells[t][s].phi_rms;
                fr->shells[s].P_mean    += tshells[t][s].P_mean;
                fr->shells[s].theta_rms += tshells[t][s].theta_rms;
                fr->shells[s].phi_v_rms += tshells[t][s].phi_v_rms;
                fr->shells[s].v_radial  += tshells[t][s].v_radial;
                fr->shells[s].count     += tshells[t][s].count;
            }
            fr->E_phi      += t_E_phi[t];
            fr->E_theta    += t_E_theta[t];
            fr->E_kin      += t_E_kin[t];
            fr->P_total    += t_P_total[t];
            fr->theta_total += t_theta_total[t];
            fr->binding_vol_001 += t_bv001[t];
            fr->binding_vol_01  += t_bv01[t];
            fr->binding_vol_1   += t_bv1[t];

            /* phi_max: take global max */
            if (t_phi_max[t] > fr->phi_max) {
                fr->phi_max = t_phi_max[t];
                fr->phi_max_rx = t_phi_max_rx[t];
                fr->phi_max_ry = t_phi_max_ry[t];
                fr->phi_max_rz = t_phi_max_rz[t];
            }
            if (t_P_max[t] > fr->P_max) {
                fr->P_max = t_P_max[t];
                fr->P_max_r = t_P_max_r[t];
            }

            for (int o = 0; o < 8; o++) fr->phase_counts[o] += t_phase[t][o];
            fr->dom_x += t_dom_x[t];
            fr->dom_y += t_dom_y[t];
            fr->dom_z += t_dom_z[t];
        }

        /* Compute asymmetry */
        double E_phi_left = 0, E_phi_right = 0, P_left = 0, P_right = 0;
        for (int t = 0; t < nthreads; t++) {
            E_phi_left  += t_E_phi_left[t];
            E_phi_right += t_E_phi_right[t];
            P_left  += t_P_left[t];
            P_right += t_P_right[t];
        }
        double E_phi_sum = E_phi_left + E_phi_right;
        double P_sum = P_left + P_right;
        fr->asym_phi = (E_phi_sum > 0) ? (E_phi_left - E_phi_right) / E_phi_sum : 0;
        fr->asym_P   = (P_sum > 0) ? (P_left - P_right) / P_sum : 0;

        /* Finalize shells: convert sums to rms/mean */
        for (int s = 0; s < NSHELLS; s++) {
            long c = fr->shells[s].count;
            if (c > 0) {
                fr->shells[s].phi_rms   = sqrt(fr->shells[s].phi_rms / c);
                fr->shells[s].P_mean    /= c;
                fr->shells[s].theta_rms = sqrt(fr->shells[s].theta_rms / c);
                fr->shells[s].phi_v_rms = sqrt(fr->shells[s].phi_v_rms / c);
                fr->shells[s].v_radial  /= c;
            }
        }

        /* Core and r=5 values */
        /* core: shells 0-3 (r < 2) */
        {
            double sum_phi2 = 0, sum_th2 = 0;
            long cnt = 0;
            for (int s = 0; s < 4 && s < NSHELLS; s++) {
                /* We need to reconstruct from rms: rms^2 * count */
                double p2 = fr->shells[s].phi_rms * fr->shells[s].phi_rms * fr->shells[s].count;
                double t2 = fr->shells[s].theta_rms * fr->shells[s].theta_rms * fr->shells[s].count;
                sum_phi2 += p2;
                sum_th2  += t2;
                cnt += fr->shells[s].count;
            }
            fr->phi_rms_core = (cnt > 0) ? sqrt(sum_phi2 / cnt) : 0;
            fr->theta_rms_core = (cnt > 0) ? sqrt(sum_th2 / cnt) : 0;
            /* v_radial at core: shells 0-3 weighted avg */
            double vr_sum = 0;
            long vr_cnt = 0;
            for (int s = 0; s < 4 && s < NSHELLS; s++) {
                vr_sum += fr->shells[s].v_radial * fr->shells[s].count;
                vr_cnt += fr->shells[s].count;
            }
            fr->v_radial_core = (vr_cnt > 0) ? vr_sum / vr_cnt : 0;
        }
        /* r=5: shells 9-10 (r = 4.5 to 5.5) */
        {
            double sum_phi2 = 0, sum_th2 = 0;
            long cnt = 0;
            for (int s = 9; s <= 10 && s < NSHELLS; s++) {
                double p2 = fr->shells[s].phi_rms * fr->shells[s].phi_rms * fr->shells[s].count;
                double t2 = fr->shells[s].theta_rms * fr->shells[s].theta_rms * fr->shells[s].count;
                sum_phi2 += p2;
                sum_th2  += t2;
                cnt += fr->shells[s].count;
            }
            fr->phi_rms_r5 = (cnt > 0) ? sqrt(sum_phi2 / cnt) : 0;
            fr->theta_rms_r5 = (cnt > 0) ? sqrt(sum_th2 / cnt) : 0;
            double vr_sum = 0;
            long vr_cnt = 0;
            for (int s = 9; s <= 10 && s < NSHELLS; s++) {
                vr_sum += fr->shells[s].v_radial * fr->shells[s].count;
                vr_cnt += fr->shells[s].count;
            }
            fr->v_radial_r5 = (vr_cnt > 0) ? vr_sum / vr_cnt : 0;
        }

        /* Binding surface area: count voxels where |P| crosses 0.001 threshold
         * (i.e., |P| > 0.001 but at least one face neighbor has |P| < 0.001) */
        long binding_area = 0;
        #pragma omp parallel for reduction(+:binding_area) schedule(static)
        for (long idx = 0; idx < N3; idx++) {
            double P = fabs(phi_x[idx] * phi_y[idx] * phi_z[idx]);
            if (P <= 0.001) continue;

            int iz = idx / NN;
            int iy = (idx - (long)iz * NN) / N;
            int ix = idx - (long)iz * NN - iy * N;
            double x = -L + ix * dx;
            double y = -L + iy * dx;
            double z = -L + iz * dx;
            double r = sqrt((x - fr->cx) * (x - fr->cx) + (y - fr->cy) * (y - fr->cy) + (z - fr->cz) * (z - fr->cz));
            if (r > R_MAX) continue;

            /* Check 6 face neighbors */
            int is_surface = 0;
            int nbr[6][3] = {
                {ix-1,iy,iz}, {ix+1,iy,iz},
                {ix,iy-1,iz}, {ix,iy+1,iz},
                {ix,iy,iz-1}, {ix,iy,iz+1}
            };
            for (int n = 0; n < 6 && !is_surface; n++) {
                int nx = nbr[n][0], ny = nbr[n][1], nz = nbr[n][2];
                if (nx < 0 || nx >= N || ny < 0 || ny >= N || nz < 0 || nz >= N) {
                    is_surface = 1;
                } else {
                    long nidx = (long)nz * NN + ny * N + nx;
                    double Pn = fabs(phi_x[nidx] * phi_y[nidx] * phi_z[nidx]);
                    if (Pn <= 0.001) is_surface = 1;
                }
            }
            if (is_surface) binding_area++;
        }
        fr->binding_area = binding_area;

        /* Phase sign at centroid voxel */
        {
            int cix = (int)((fr->cx + L) / dx + 0.5);
            int ciy = (int)((fr->cy + L) / dx + 0.5);
            int ciz = (int)((fr->cz + L) / dx + 0.5);
            if (cix >= 0 && cix < N && ciy >= 0 && ciy < N && ciz >= 0 && ciz < N) {
                long cidx = (long)ciz * NN + ciy * N + cix;
                fr->sign_cx = (phi_x[cidx] >= 0) ? 1 : -1;
                fr->sign_cy = (phi_y[cidx] >= 0) ? 1 : -1;
                fr->sign_cz = (phi_z[cidx] >= 0) ? 1 : -1;
            }
        }

        /* --- Write outputs --- */

        /* Summary row */
        fprintf(fp_summary, "%.2f\t%.4f\t%.4f\t%.4f\t%.6f\t%.6f\t%.6f\t"
                "%.6f\t%.6f\t%.3f\t%ld\t%ld\t%ld\t"
                "%ld\t%.4f\t%.4f\t%.4f\t%.6f\t%.6f\t"
                "%.6f\t%.6f\t%.4f\t%.6f\t%.6f\n",
                fr->t, fr->cx, fr->cy, fr->cz,
                fr->phi_rms_core, fr->phi_rms_r5, fr->phi_max,
                fr->P_total, fr->P_max, fr->P_max_r,
                fr->binding_vol_001, fr->binding_vol_01, fr->binding_vol_1,
                fr->binding_area,
                fr->E_phi, fr->E_theta, fr->E_kin,
                fr->asym_phi, fr->asym_P,
                fr->theta_rms_core, fr->theta_rms_r5, fr->theta_total,
                fr->v_radial_core, fr->v_radial_r5);

        /* Profile rows */
        for (int s = 0; s < NSHELLS; s++) {
            double r = (s + 0.5) * DR;
            fprintf(fp_profiles, "%d\t%.2f\t%.6f\t%.6e\t%.6f\t%.6f\t%.6f\n",
                    fi, r,
                    fr->shells[s].phi_rms,
                    fr->shells[s].P_mean,
                    fr->shells[s].theta_rms,
                    fr->shells[s].phi_v_rms,
                    fr->shells[s].v_radial);
        }

        /* Free per-thread arrays */
        free(tshells);
        free(t_E_phi); free(t_E_theta); free(t_E_kin);
        free(t_P_total); free(t_phi_max);
        free(t_phi_max_rx); free(t_phi_max_ry); free(t_phi_max_rz);
        free(t_P_max); free(t_P_max_r);
        free(t_E_phi_left); free(t_E_phi_right);
        free(t_P_left); free(t_P_right);
        free(t_theta_total);
        free(t_bv001); free(t_bv01); free(t_bv1);
        free(t_phase);
        free(t_dom_x); free(t_dom_y); free(t_dom_z);

        /* Free field arrays */
        free(phi_x); free(phi_y); free(phi_z);
        free(th_x);  free(th_y);  free(th_z);
        free(vx);    free(vy);    free(vz);
        free(vth_x); free(vth_y); free(vth_z);

        double t_end = omp_get_wtime();
        printf("done %.1fs\n", t_end - t0);
    }

    fclose(fp_summary);
    fclose(fp_profiles);
    sfa_close(sfa);
    free(buf);

    /* === Print summary to stdout === */
    printf("\n=== FRAME-BY-FRAME SUMMARY ===\n");
    printf("%-6s  %-8s  %-8s  %-8s  %-10s  %-10s  %-10s  %-10s  %-10s  %-10s  %-8s  %-8s  %-6s %-6s %-6s  %-7s\n",
           "Frame", "t", "cx", "cy", "phi_rms_c", "phi_max", "P_total",
           "P_max", "E_phi", "E_kin", "asym_phi", "asym_P",
           "s_x", "s_y", "s_z", "v_r_c");
    printf("------  --------  --------  --------  ----------  ----------  ----------  ----------  ----------  ----------  --------  --------  ------  ------  ------  -------\n");

    double P_total_max = 0, P_total_min = 1e30;
    int P_peak_frame = 0, P_trough_frame = 0;

    for (int fi = 0; fi < nf; fi++) {
        FrameResult *fr = &results[fi];
        printf("%-6d  %8.1f  %8.3f  %8.3f  %10.5f  %10.5f  %10.5f  %10.5f  %10.3f  %10.3f  %8.5f  %8.5f  %+4d   %+4d   %+4d   %+7.4f\n",
               fi, fr->t, fr->cx, fr->cy,
               fr->phi_rms_core, fr->phi_max, fr->P_total,
               fr->P_max, fr->E_phi, fr->E_kin,
               fr->asym_phi, fr->asym_P,
               fr->sign_cx, fr->sign_cy, fr->sign_cz,
               fr->v_radial_core);
        if (fr->P_total > P_total_max) { P_total_max = fr->P_total; P_peak_frame = fi; }
        if (fr->P_total < P_total_min) { P_total_min = fr->P_total; P_trough_frame = fi; }
    }

    printf("\n=== BINDING PEAKS AND TROUGHS ===\n");
    printf("Peak binding:   frame %d (t=%.1f), P_total=%.5f\n",
           P_peak_frame, results[P_peak_frame].t, P_total_max);
    printf("Trough binding: frame %d (t=%.1f), P_total=%.5f\n",
           P_trough_frame, results[P_trough_frame].t, P_total_min);
    printf("Binding range:  %.2f%%\n", 100.0 * (P_total_max - P_total_min) / ((P_total_max + P_total_min) / 2.0));

    printf("\n=== ASYMMETRY ANALYSIS ===\n");
    printf("%-6s  %10s  %10s  %10s  %10s  %s\n",
           "Frame", "asym_phi", "asym_P", "v_rad_core", "centroid_x", "interpretation");
    for (int fi = 0; fi < nf; fi++) {
        FrameResult *fr = &results[fi];
        const char *interp = "";
        if (fr->asym_P > 0.01) interp = "LEFT-heavy (binding pulls left)";
        else if (fr->asym_P < -0.01) interp = "RIGHT-heavy (binding pulls right)";
        else interp = "symmetric";

        printf("%-6d  %+10.5f  %+10.5f  %+10.5f  %10.4f  %s\n",
               fi, fr->asym_phi, fr->asym_P, fr->v_radial_core, fr->cx, interp);
    }

    printf("\n=== PHASE SIGN DISTRIBUTION (r<10) ===\n");
    printf("%-6s  %8s  %8s  %8s  %8s  %8s  %8s  %8s  %8s\n",
           "Frame", "(+,+,+)", "(+,+,-)", "(+,-,+)", "(+,-,-)",
           "(-,+,+)", "(-,+,-)", "(-,-,+)", "(-,-,-)");
    for (int fi = 0; fi < nf; fi++) {
        FrameResult *fr = &results[fi];
        long total = 0;
        for (int o = 0; o < 8; o++) total += fr->phase_counts[o];
        printf("%-6d", fi);
        for (int o = 7; o >= 0; o--) {
            printf("  %7.2f%%", 100.0 * fr->phase_counts[o] / (total > 0 ? total : 1));
        }
        printf("\n");
    }

    printf("\n=== DOMINANT PHI COMPONENT AT CORE (r<5) ===\n");
    for (int fi = 0; fi < nf; fi++) {
        FrameResult *fr = &results[fi];
        long total = fr->dom_x + fr->dom_y + fr->dom_z;
        printf("Frame %d: phi_x=%5.1f%%  phi_y=%5.1f%%  phi_z=%5.1f%%\n",
               fi,
               100.0 * fr->dom_x / (total > 0 ? total : 1),
               100.0 * fr->dom_y / (total > 0 ? total : 1),
               100.0 * fr->dom_z / (total > 0 ? total : 1));
    }

    printf("\n=== BREATHING-DRIFT CORRELATION ===\n");
    printf("Checking if binding asymmetry correlates with radial velocity phase...\n");
    /* Compare v_radial_core sign with asym_P change */
    for (int fi = 1; fi < nf; fi++) {
        double d_cx = results[fi].cx - results[fi - 1].cx;
        double d_P = results[fi].P_total - results[fi - 1].P_total;
        double asym = results[fi].asym_P;
        printf("  Frame %d→%d: Δcx=%+.4f, ΔP=%+.5f, asym_P=%+.5f, v_rad=%+.5f  ",
               fi - 1, fi, d_cx, d_P, asym, results[fi].v_radial_core);
        if ((asym > 0 && d_cx > 0) || (asym < 0 && d_cx < 0))
            printf("[CORRELATED: binding pulls in drift direction]\n");
        else if (fabs(asym) < 0.005)
            printf("[SYMMETRIC]\n");
        else
            printf("[ANTI-CORRELATED]\n");
    }

    printf("\n=== THETA STABILITY ===\n");
    for (int fi = 0; fi < nf; fi++) {
        printf("Frame %d: theta_total=%.4f  theta_rms_core=%.6f  theta_rms_r5=%.6f\n",
               fi, results[fi].theta_total, results[fi].theta_rms_core, results[fi].theta_rms_r5);
    }

    printf("\nWrote: phase_binding.tsv, phase_profiles.tsv\n");

    free(results);
    return 0;
}
