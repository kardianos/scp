/*  sfa_qcomp.c — per-component ("quark-level") bookkeeping for complex SFA runs
 *
 *  The complexified theory has three component fields Phi_a (a = x,y,z in
 *  column naming; a = 0,1,2 here). The potential binds only where ALL THREE
 *  overlap (s = prod_a |Phi_a|^2), making the components quark-like: this tool
 *  reports, per frame and per component a:
 *
 *    Q_a      = int (u_a vdot_a - v_a udot_a) dV     (component Noether charge;
 *               sums to the diagonal-U(1) charge of the phi sector)
 *    mass_a   = int (u_a^2 + v_a^2) dV               (component density norm)
 *    centroid (rho2_a-weighted, min-image)           (component position)
 *    rms_a    (rho2_a-weighted rms radius about centroid)
 *
 *  Build: gcc -O3 -fopenmp -o sfa_qcomp sfa_qcomp.c -lzstd -lm
 *  Usage: sfa_qcomp input.sfa [--tsv out.tsv]
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define SFA_IMPLEMENTATION
#include "../format/sfa.h"

static int find_column(const SFA *sfa, const char *name) {
    for (uint32_t c = 0; c < sfa->n_columns; c++)
        if (strncmp(sfa->columns[c].name, name, sizeof(sfa->columns[c].name)) == 0)
            return (int)c;
    return -1;
}

static float *extract_column_f32(const void *buf, const SFA *sfa, int col_idx) {
    uint64_t N3 = sfa->N_total, off = 0;
    for (int c = 0; c < col_idx; c++)
        off += N3 * sfa_dtype_size[sfa->columns[c].dtype];
    int dtype = sfa->columns[col_idx].dtype;
    const uint8_t *src = (const uint8_t *)buf + off;
    float *arr = (float *)malloc(N3 * sizeof(float));
    if (!arr) { fprintf(stderr, "error: malloc\n"); exit(1); }
    if (dtype == SFA_F32) memcpy(arr, src, N3 * sizeof(float));
    else if (dtype == SFA_F64) {
        const double *d = (const double *)src;
        #pragma omp parallel for
        for (uint64_t i = 0; i < N3; i++) arr[i] = (float)d[i];
    } else if (dtype == SFA_F16) {
        const uint16_t *h = (const uint16_t *)src;
        #pragma omp parallel for
        for (uint64_t i = 0; i < N3; i++) arr[i] = sfa_f16_to_f32(h[i]);
    } else { fprintf(stderr, "error: dtype %d\n", dtype); exit(1); }
    return arr;
}

int main(int argc, char **argv) {
    const char *in_path = NULL, *tsv_path = NULL;
    for (int i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "--tsv") && i + 1 < argc) tsv_path = argv[++i];
        else in_path = argv[i];
    }
    if (!in_path) {
        fprintf(stderr, "usage: %s input.sfa [--tsv out.tsv]\n", argv[0]);
        return 1;
    }
    SFA *sfa = sfa_open(in_path);
    if (!sfa) { fprintf(stderr, "error: cannot open %s\n", in_path); return 1; }

    static const char *u_n[3]  = {"phi_x","phi_y","phi_z"};
    static const char *v_n[3]  = {"phiim_x","phiim_y","phiim_z"};
    static const char *ud_n[3] = {"phi_vx","phi_vy","phi_vz"};
    static const char *vd_n[3] = {"phiim_vx","phiim_vy","phiim_vz"};
    int cu[3], cv[3], cud[3], cvd[3];
    for (int a = 0; a < 3; a++) {
        cu[a] = find_column(sfa, u_n[a]);   cv[a] = find_column(sfa, v_n[a]);
        cud[a] = find_column(sfa, ud_n[a]); cvd[a] = find_column(sfa, vd_n[a]);
        if (cu[a] < 0 || cv[a] < 0 || cud[a] < 0 || cvd[a] < 0) {
            fprintf(stderr, "error: %s lacks complex columns\n", in_path);
            return 1;
        }
    }

    int N = (int)sfa->Nx;
    double L = sfa->Lx, dx = 2.0 * L / (N - 1), dV = dx * dx * dx;
    long N3 = (long)N * N * N, NN = (long)N * N;

    FILE *tsv = NULL;
    if (tsv_path) {
        tsv = fopen(tsv_path, "w");
        if (!tsv) { fprintf(stderr, "error: cannot write %s\n", tsv_path); return 1; }
        fprintf(tsv, "frame\tt\tcomp\tQ\tmass\tcx\tcy\tcz\trms\n");
    }

    void *buf = malloc(sfa->frame_bytes);
    printf("# %s: N=%d L=%g frames=%u\n", in_path, N, L, sfa->total_frames);
    printf("%5s %8s %4s %12s %12s %8s %8s %8s %7s\n",
           "frame", "t", "a", "Q_a", "mass_a", "cx", "cy", "cz", "rms");

    for (uint32_t fi = 0; fi < sfa->total_frames; fi++) {
        double t = sfa_frame_time(sfa, fi);
        if (sfa_read_frame(sfa, fi, buf) != 0) break;
        double Qtot = 0;
        for (int a = 0; a < 3; a++) {
            float *u = extract_column_f32(buf, sfa, cu[a]);
            float *v = extract_column_f32(buf, sfa, cv[a]);
            float *du = extract_column_f32(buf, sfa, cud[a]);
            float *dv = extract_column_f32(buf, sfa, cvd[a]);
            double Q = 0, mass = 0, sx = 0, sy = 0, sz = 0, srr = 0;
            /* centroid via angular mean (periodic-safe) */
            double cxr = 0, cxi = 0, cyr = 0, cyi = 0, czr = 0, czi = 0;
            #pragma omp parallel for reduction(+:Q,mass,cxr,cxi,cyr,cyi,czr,czi)
            for (long i = 0; i < N3; i++) {
                double r2 = (double)u[i]*u[i] + (double)v[i]*v[i];
                Q += (double)u[i]*dv[i] - (double)v[i]*du[i];
                mass += r2;
                int ix = (int)(i / NN), iy = (int)((i / N) % N), iz = (int)(i % N);
                double ax = 2*M_PI*ix/N, ay = 2*M_PI*iy/N, az = 2*M_PI*iz/N;
                cxr += r2*cos(ax); cxi += r2*sin(ax);
                cyr += r2*cos(ay); cyi += r2*sin(ay);
                czr += r2*cos(az); czi += r2*sin(az);
            }
            double ax_ = atan2(cxi, cxr); if (ax_ < 0) ax_ += 2*M_PI;
            double ay_ = atan2(cyi, cyr); if (ay_ < 0) ay_ += 2*M_PI;
            double az_ = atan2(czi, czr); if (az_ < 0) az_ += 2*M_PI;
            double cx = -L + ax_/(2*M_PI)*2*L;
            double cy = -L + ay_/(2*M_PI)*2*L;
            double cz = -L + az_/(2*M_PI)*2*L;
            #pragma omp parallel for reduction(+:srr)
            for (long i = 0; i < N3; i++) {
                double r2 = (double)u[i]*u[i] + (double)v[i]*v[i];
                int ix = (int)(i / NN), iy = (int)((i / N) % N), iz = (int)(i % N);
                double px = -L + ix*dx - cx, py = -L + iy*dx - cy, pz = -L + iz*dx - cz;
                if (px > L) px -= 2*L; if (px < -L) px += 2*L;
                if (py > L) py -= 2*L; if (py < -L) py += 2*L;
                if (pz > L) pz -= 2*L; if (pz < -L) pz += 2*L;
                srr += r2*(px*px + py*py + pz*pz);
            }
            double rms = mass > 0 ? sqrt(srr/mass) : 0;
            Q *= dV; mass *= dV; Qtot += Q;
            printf("%5u %8.2f %4d %12.5f %12.5f %8.3f %8.3f %8.3f %7.3f\n",
                   fi, t, a, Q, mass, cx, cy, cz, rms);
            if (tsv) fprintf(tsv, "%u\t%.4f\t%d\t%.8e\t%.8e\t%.4f\t%.4f\t%.4f\t%.4f\n",
                             fi, t, a, Q, mass, cx, cy, cz, rms);
            free(u); free(v); free(du); free(dv);
        }
        printf("%5s %8s %4s %12.5f   (Q_phi total)\n", "", "", "sum", Qtot);
    }
    if (tsv) fclose(tsv);
    sfa_close(sfa);
    return 0;
}
