/* slice_analysis.c — Planar cross-section analysis of SFA frame
 *
 * Finds the two strongest P-maxima, then outputs 2D slice data
 * through each core in all three planes (XY, XZ, YZ).
 *
 * Build: gcc -O3 -o slice_analysis slice_analysis.c -lzstd -lm
 * Usage: ./slice_analysis input.sfa [frame_idx]
 *
 * Outputs TSV files:
 *   slice_xy_z{N}.tsv  — XY plane at z=N (through core)
 *   slice_xz_y{N}.tsv  — XZ plane at y=N
 *   slice_yz_x{N}.tsv  — YZ plane at x=N
 * Each row: x/y/z  coord1  coord2  phi_rms  |P|  theta_rms  curl_rms  harden
 */

#define SFA_IMPLEMENTATION
#include "../../sfa/format/sfa.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

int main(int argc, char **argv) {
    if (argc < 2) {
        fprintf(stderr, "Usage: %s input.sfa [frame_idx] [output_dir]\n", argv[0]);
        return 1;
    }

    SFA *sfa = sfa_open(argv[1]);
    if (!sfa) { fprintf(stderr, "Cannot open %s\n", argv[1]); return 1; }

    int frame = (argc > 2) ? atoi(argv[2]) : sfa->total_frames - 1;
    const char *outdir = (argc > 3) ? argv[3] : ".";
    int N = sfa->Nx;
    double L = sfa->Lx;
    double dx = 2.0 * L / (N - 1);
    long N3 = (long)N * N * N;
    int NN = N * N;

    printf("File: %s  N=%d L=%.1f dx=%.4f  Frame %d\n", argv[1], N, L, dx, frame);

    /* Read frame */
    size_t frame_bytes = 0;
    for (uint32_t c = 0; c < sfa->n_columns; c++)
        frame_bytes += N3 * sfa_dtype_size[sfa->columns[c].dtype];
    void *buf = malloc(frame_bytes);
    if (sfa_read_frame(sfa, frame, buf) < 0) {
        fprintf(stderr, "Cannot read frame %d\n", frame);
        return 1;
    }

    double t = sfa_frame_time(sfa, frame);
    printf("Time: %.4f\n", t);

    float *phi[3], *theta[3];
    size_t col_bytes = N3 * sizeof(float);
    for (int a = 0; a < 3; a++) {
        phi[a] = (float*)((char*)buf + a * col_bytes);
        theta[a] = (float*)((char*)buf + (3 + a) * col_bytes);
    }

    /* Find two strongest P-maxima */
    typedef struct { int i, j, k; double P; } Peak;
    Peak pk1 = {0,0,0,0}, pk2 = {0,0,0,0};
    int min_sep = N / 6;  /* minimum separation between peaks */

    for (long idx = 0; idx < N3; idx++) {
        double P = fabs((double)phi[0][idx] * phi[1][idx] * phi[2][idx]);
        if (P > pk1.P) {
            pk1.P = P; pk1.i = idx/NN; pk1.j = (idx/N)%N; pk1.k = idx%N;
        }
    }
    for (long idx = 0; idx < N3; idx++) {
        int i = idx/NN, j = (idx/N)%N, k = idx%N;
        int di = abs(i - pk1.i), dj = abs(j - pk1.j), dk = abs(k - pk1.k);
        /* Handle periodic wrapping */
        if (di > N/2) di = N - di;
        if (dj > N/2) dj = N - dj;
        if (dk > N/2) dk = N - dk;
        if (di + dj + dk < min_sep) continue;
        double P = fabs((double)phi[0][idx] * phi[1][idx] * phi[2][idx]);
        if (P > pk2.P) {
            pk2.P = P; pk2.i = i; pk2.j = j; pk2.k = k;
        }
    }

    printf("\nPeak 1: grid(%d,%d,%d) = world(%.2f,%.2f,%.2f)  |P|=%.5f\n",
           pk1.i, pk1.j, pk1.k,
           -L + pk1.i*dx, -L + pk1.j*dx, -L + pk1.k*dx, pk1.P);
    printf("Peak 2: grid(%d,%d,%d) = world(%.2f,%.2f,%.2f)  |P|=%.5f\n",
           pk2.i, pk2.j, pk2.k,
           -L + pk2.i*dx, -L + pk2.j*dx, -L + pk2.k*dx, pk2.P);
    double sep = sqrt(pow((pk1.i-pk2.i)*dx,2)+pow((pk1.j-pk2.j)*dx,2)+pow((pk1.k-pk2.k)*dx,2));
    printf("Separation: %.2f code units (%d grid cells)\n\n", sep,
           abs(pk1.i-pk2.i)+abs(pk1.j-pk2.j)+abs(pk1.k-pk2.k));

    /* Midpoint between the two peaks */
    int mi = (pk1.i + pk2.i) / 2;
    int mj = (pk1.j + pk2.j) / 2;
    int mk = (pk1.k + pk2.k) / 2;

    double idx1 = 1.0 / (2.0 * dx);

    /* Helper: compute fields at a voxel */
    #define COMPUTE_FIELDS(i,j,k) do { \
        long _idx = (long)(i)*NN + (j)*N + (k); \
        double p0=phi[0][_idx], p1=phi[1][_idx], p2=phi[2][_idx]; \
        double t0=theta[0][_idx], t1=theta[1][_idx], t2=theta[2][_idx]; \
        f_phi_rms = sqrt(p0*p0 + p1*p1 + p2*p2); \
        f_P = fabs(p0*p1*p2); \
        f_theta_rms = sqrt(t0*t0 + t1*t1 + t2*t2); \
        int _ip=((i)+1)%N, _im=((i)-1+N)%N; \
        int _jp=((j)+1)%N, _jm=((j)-1+N)%N; \
        int _kp=((k)+1)%N, _km=((k)-1+N)%N; \
        long _nip=(long)_ip*NN+(j)*N+(k), _nim=(long)_im*NN+(j)*N+(k); \
        long _njp=(long)(i)*NN+_jp*N+(k), _njm=(long)(i)*NN+_jm*N+(k); \
        long _nkp=(long)(i)*NN+(j)*N+_kp, _nkm=(long)(i)*NN+(j)*N+_km; \
        double c0=(phi[2][_njp]-phi[2][_njm]-phi[1][_nkp]+phi[1][_nkm])*idx1; \
        double c1=(phi[0][_nkp]-phi[0][_nkm]-phi[2][_nip]+phi[2][_nim])*idx1; \
        double c2=(phi[1][_nip]-phi[1][_nim]-phi[0][_njp]+phi[0][_njm])*idx1; \
        f_curl_rms = sqrt(c0*c0 + c1*c1 + c2*c2); \
        f_harden = (t0*t0+t1*t1+t2*t2) * (c0*c0+c1*c1+c2*c2); \
    } while(0)

    double f_phi_rms, f_P, f_theta_rms, f_curl_rms, f_harden;

    /* Write slices through the midpoint */
    char path[1024];

    /* XY slice at z = midpoint_k */
    snprintf(path, 1024, "%s/slice_xy_z%d.tsv", outdir, mk);
    printf("Writing %s (XY plane at z=%d, world z=%.2f)\n", path, mk, -L+mk*dx);
    FILE *fp = fopen(path, "w");
    fprintf(fp, "i\tj\tx\ty\tphi_rms\tP\ttheta_rms\tcurl_rms\tharden\n");
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            COMPUTE_FIELDS(i, j, mk);
            fprintf(fp, "%d\t%d\t%.3f\t%.3f\t%.5f\t%.6f\t%.5f\t%.5f\t%.8f\n",
                    i, j, -L+i*dx, -L+j*dx, f_phi_rms, f_P, f_theta_rms, f_curl_rms, f_harden);
        }
    }
    fclose(fp);

    /* XZ slice at y = midpoint_j */
    snprintf(path, 1024, "%s/slice_xz_y%d.tsv", outdir, mj);
    printf("Writing %s (XZ plane at y=%d, world y=%.2f)\n", path, mj, -L+mj*dx);
    fp = fopen(path, "w");
    fprintf(fp, "i\tk\tx\tz\tphi_rms\tP\ttheta_rms\tcurl_rms\tharden\n");
    for (int i = 0; i < N; i++) {
        for (int k = 0; k < N; k++) {
            COMPUTE_FIELDS(i, mj, k);
            fprintf(fp, "%d\t%d\t%.3f\t%.3f\t%.5f\t%.6f\t%.5f\t%.5f\t%.8f\n",
                    i, k, -L+i*dx, -L+k*dx, f_phi_rms, f_P, f_theta_rms, f_curl_rms, f_harden);
        }
    }
    fclose(fp);

    /* YZ slice at x = midpoint_i */
    snprintf(path, 1024, "%s/slice_yz_x%d.tsv", outdir, mi);
    printf("Writing %s (YZ plane at x=%d, world x=%.2f)\n", path, mi, -L+mi*dx);
    fp = fopen(path, "w");
    fprintf(fp, "j\tk\ty\tz\tphi_rms\tP\ttheta_rms\tcurl_rms\tharden\n");
    for (int j = 0; j < N; j++) {
        for (int k = 0; k < N; k++) {
            COMPUTE_FIELDS(mi, j, k);
            fprintf(fp, "%d\t%d\t%.3f\t%.3f\t%.5f\t%.6f\t%.5f\t%.5f\t%.8f\n",
                    j, k, -L+j*dx, -L+k*dx, f_phi_rms, f_P, f_theta_rms, f_curl_rms, f_harden);
        }
    }
    fclose(fp);

    /* Also write slices through each peak */
    for (int p = 0; p < 2; p++) {
        Peak *pk = (p == 0) ? &pk1 : &pk2;
        snprintf(path, 1024, "%s/slice_xy_peak%d_z%d.tsv", outdir, p+1, pk->k);
        printf("Writing %s (XY through peak %d at z=%d)\n", path, p+1, pk->k);
        fp = fopen(path, "w");
        fprintf(fp, "i\tj\tx\ty\tphi_rms\tP\ttheta_rms\tcurl_rms\tharden\n");
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                COMPUTE_FIELDS(i, j, pk->k);
                fprintf(fp, "%d\t%d\t%.3f\t%.3f\t%.5f\t%.6f\t%.5f\t%.5f\t%.8f\n",
                        i, j, -L+i*dx, -L+j*dx, f_phi_rms, f_P, f_theta_rms, f_curl_rms, f_harden);
            }
        }
        fclose(fp);
    }

    printf("\nDone. Slice files in %s/\n", outdir);
    free(buf);
    sfa_close(sfa);
    return 0;
}
