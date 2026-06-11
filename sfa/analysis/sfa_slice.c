/*  sfa_slice.c — 2D slice / 1D lineout extractor for SFA simulation output
 *
 *  Complex- and gauge-aware: derives the physically meaningful densities
 *  from whatever columns are present (12-col real, 24-col complex,
 *  30-col gauged) and dumps planes / lines for offline visual inspection
 *  (render with sfa/analysis/render_slices.py).
 *
 *  Derived quantities (auto-selected by available columns):
 *    rho2  = sum_a (u_a^2 + v_a^2)            field density (phase-invariant)
 *    s     = prod_a (u_a^2 + v_a^2)           potential source
 *    rhoQ  = sum_a (u v' - v u' + tu tv' - tv tu')   Noether charge density
 *    th2   = sum_a (tu_a^2 + tv_a^2)          theta density
 *    Emag  = |E|                              gauge electric field magnitude
 *
 *  Slice output: <prefix>.meta.txt (text header) + <prefix>.f32
 *  (raw float32, layout [frame][quantity][n2][n1], n1 fastest).
 *  Lineout output: <prefix>.line.tsv — all raw columns + rho2 along an
 *  axis-aligned line through the per-frame rho2 maximum.
 *
 *  Build: gcc -O3 -fopenmp -o sfa_slice sfa_slice.c -lzstd -lm
 *
 *  Usage:
 *    sfa_slice input.sfa [--axis x|y|z] [--pos f] [--frames 0,4,-1]
 *              [--out prefix] [--line x|y|z]
 *
 *    --axis   slice normal (default z; slice plane passes through --pos)
 *    --pos    slice position in physical coords (default 0.0)
 *    --frames comma list of frame indices, negative = from end (default all)
 *    --line   ALSO write a 1D lineout along this axis through the
 *             per-frame rho2-max voxel (phase-tilt analysis input)
 *    --out    output prefix (default: input basename)
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
    if (!arr) { fprintf(stderr, "error: malloc failed\n"); exit(1); }
    if (dtype == SFA_F32) {
        memcpy(arr, src, N3 * sizeof(float));
    } else if (dtype == SFA_F64) {
        const double *d = (const double *)src;
        #pragma omp parallel for
        for (uint64_t i = 0; i < N3; i++) arr[i] = (float)d[i];
    } else if (dtype == SFA_F16) {
        const uint16_t *h = (const uint16_t *)src;
        #pragma omp parallel for
        for (uint64_t i = 0; i < N3; i++) arr[i] = sfa_f16_to_f32(h[i]);
    } else {
        fprintf(stderr, "error: unsupported dtype %d\n", dtype);
        exit(1);
    }
    return arr;
}

static void usage(const char *prog) {
    fprintf(stderr,
        "usage: %s input.sfa [--axis x|y|z] [--pos f] [--frames 0,4,-1] "
        "[--out prefix] [--line x|y|z]\n", prog);
}

int main(int argc, char **argv) {
    const char *in_path = NULL, *out_prefix = NULL;
    char axis = 'z', line_axis = 0;
    double pos = 0.0;
    const char *frames_arg = NULL;

    for (int i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "--axis") && i + 1 < argc)        axis = argv[++i][0];
        else if (!strcmp(argv[i], "--pos") && i + 1 < argc)    pos = atof(argv[++i]);
        else if (!strcmp(argv[i], "--frames") && i + 1 < argc) frames_arg = argv[++i];
        else if (!strcmp(argv[i], "--out") && i + 1 < argc)    out_prefix = argv[++i];
        else if (!strcmp(argv[i], "--line") && i + 1 < argc)   line_axis = argv[++i][0];
        else if (!in_path)                                     in_path = argv[i];
        else { usage(argv[0]); return 1; }
    }
    if (!in_path || (axis != 'x' && axis != 'y' && axis != 'z')) { usage(argv[0]); return 1; }

    char prefix[1024];
    if (out_prefix) {
        snprintf(prefix, sizeof(prefix), "%s", out_prefix);
    } else {
        const char *base = strrchr(in_path, '/');
        base = base ? base + 1 : in_path;
        snprintf(prefix, sizeof(prefix), "%s", base);
        char *dot = strrchr(prefix, '.');
        if (dot && !strcmp(dot, ".sfa")) *dot = '\0';
    }

    SFA *sfa = sfa_open(in_path);
    if (!sfa) { fprintf(stderr, "error: cannot open %s\n", in_path); return 1; }

    int N = (int)sfa->Nx;
    double L = sfa->Lx;
    double dx = 2.0 * L / (N - 1);
    long NN = (long)N * N, N3 = (long)N * N * N;
    uint32_t total = sfa->total_frames;

    /* ---- column groups ---- */
    static const char *u_n[3]   = {"phi_x","phi_y","phi_z"};
    static const char *v_n[3]   = {"phiim_x","phiim_y","phiim_z"};
    static const char *ud_n[3]  = {"phi_vx","phi_vy","phi_vz"};
    static const char *vd_n[3]  = {"phiim_vx","phiim_vy","phiim_vz"};
    static const char *tu_n[3]  = {"theta_x","theta_y","theta_z"};
    static const char *tv_n[3]  = {"thetaim_x","thetaim_y","thetaim_z"};
    static const char *tud_n[3] = {"theta_vx","theta_vy","theta_vz"};
    static const char *tvd_n[3] = {"thetaim_vx","thetaim_vy","thetaim_vz"};
    static const char *E_n[3]   = {"E_x","E_y","E_z"};
    int cu[3],cv[3],cud[3],cvd[3],ctu[3],ctv[3],ctud[3],ctvd[3],cE[3];
    int has_im = 1, has_th = 1, has_thim = 1, has_E = 1, has_vel = 1;
    for (int a = 0; a < 3; a++) {
        cu[a]=find_column(sfa,u_n[a]);   cv[a]=find_column(sfa,v_n[a]);
        cud[a]=find_column(sfa,ud_n[a]); cvd[a]=find_column(sfa,vd_n[a]);
        ctu[a]=find_column(sfa,tu_n[a]); ctv[a]=find_column(sfa,tv_n[a]);
        ctud[a]=find_column(sfa,tud_n[a]); ctvd[a]=find_column(sfa,tvd_n[a]);
        cE[a]=find_column(sfa,E_n[a]);
        if (cu[a]<0) { fprintf(stderr,"error: missing %s\n",u_n[a]); return 1; }
        if (cv[a]<0||cvd[a]<0) has_im = 0;
        if (cud[a]<0) has_vel = 0;
        if (ctu[a]<0) has_th = 0;
        if (ctv[a]<0||ctvd[a]<0) has_thim = 0;
        if (cE[a]<0) has_E = 0;
    }

    /* ---- quantity list ---- */
    const char *qnames[12]; int nq = 0;
    int q_rho2=-1, q_s=-1, q_rhoQ=-1, q_th2=-1, q_E=-1, q_comp=-1;
    q_rho2 = nq; qnames[nq++] = "rho2";
    q_s    = nq; qnames[nq++] = "s";
    if (has_im) {           /* per-component densities ("quark" view) */
        q_comp = nq;
        qnames[nq++] = "rho2_0"; qnames[nq++] = "rho2_1"; qnames[nq++] = "rho2_2";
    }
    if (has_im && has_vel) { q_rhoQ = nq; qnames[nq++] = "rhoQ"; }
    if (has_th)            { q_th2  = nq; qnames[nq++] = "th2";  }
    if (has_E)             { q_E    = nq; qnames[nq++] = "Emag"; }

    /* ---- frame selection ---- */
    int *fsel = (int *)malloc(total * sizeof(int)), nf = 0;
    if (frames_arg) {
        char tmp[512]; snprintf(tmp, sizeof(tmp), "%s", frames_arg);
        for (char *tok = strtok(tmp, ","); tok; tok = strtok(NULL, ",")) {
            long fi = atol(tok);
            if (fi < 0) fi += total;
            if (fi < 0 || fi >= (long)total) {
                fprintf(stderr, "error: frame %s out of range (0..%u)\n", tok, total - 1);
                return 1;
            }
            fsel[nf++] = (int)fi;
        }
    } else {
        for (uint32_t i = 0; i < total; i++) fsel[nf++] = (int)i;
    }

    /* ---- slice geometry: axis index + in-plane dims ----
     * storage: idx = ix*N*N + iy*N + iz  (x slowest, z fastest)
     * coords:  c = -L + i*dx                                          */
    int k_slice = (int)lround((pos + L) / dx);
    if (k_slice < 0) k_slice = 0;
    if (k_slice >= N) k_slice = N - 1;

    char meta_path[1100], data_path[1100], line_path[1100];
    snprintf(meta_path, sizeof(meta_path), "%s.meta.txt", prefix);
    snprintf(data_path, sizeof(data_path), "%s.f32", prefix);
    snprintf(line_path, sizeof(line_path), "%s.line.tsv", prefix);

    FILE *fd = fopen(data_path, "wb");
    if (!fd) { fprintf(stderr, "error: cannot write %s\n", data_path); return 1; }
    FILE *fl = NULL;
    if (line_axis) {
        fl = fopen(line_path, "w");
        if (!fl) { fprintf(stderr, "error: cannot write %s\n", line_path); return 1; }
        fprintf(fl, "frame\tt\tcoord");
        for (uint32_t c = 0; c < sfa->n_columns; c++)
            fprintf(fl, "\t%s", sfa->columns[c].name);
        fprintf(fl, "\trho2\n");
    }

    void *buf = malloc(sfa->frame_bytes);
    float *plane = (float *)malloc((long)N * N * sizeof(float));
    double *times = (double *)malloc(nf * sizeof(double));
    if (!buf || !plane || !times) { fprintf(stderr, "error: alloc\n"); return 1; }

    /* full-volume scratch for the derived fields we need at once */
    float *rho2 = (float *)malloc(N3 * sizeof(float));
    float *sP   = (float *)malloc(N3 * sizeof(float));
    float *rhoQ = q_rhoQ >= 0 ? (float *)malloc(N3 * sizeof(float)) : NULL;
    float *th2  = q_th2  >= 0 ? (float *)malloc(N3 * sizeof(float)) : NULL;
    float *Emag = q_E    >= 0 ? (float *)malloc(N3 * sizeof(float)) : NULL;
    float *rca[3] = {NULL, NULL, NULL};
    if (q_comp >= 0)
        for (int a = 0; a < 3; a++) rca[a] = (float *)malloc(N3 * sizeof(float));
    if (!rho2 || !sP) { fprintf(stderr, "error: alloc\n"); return 1; }

    for (int si = 0; si < nf; si++) {
        int fi = fsel[si];
        times[si] = sfa_frame_time(sfa, fi);
        if (sfa_read_frame(sfa, fi, buf) != 0) {
            fprintf(stderr, "error: read frame %d failed\n", fi);
            return 1;
        }

        #pragma omp parallel for
        for (long i = 0; i < N3; i++) { rho2[i] = 0.f; sP[i] = 1.f; }
        if (rhoQ) memset(rhoQ, 0, N3 * sizeof(float));
        if (th2)  memset(th2,  0, N3 * sizeof(float));

        for (int a = 0; a < 3; a++) {
            float *u = extract_column_f32(buf, sfa, cu[a]);
            float *v = has_im ? extract_column_f32(buf, sfa, cv[a]) : NULL;
            float *du = (rhoQ) ? extract_column_f32(buf, sfa, cud[a]) : NULL;
            float *dv = (rhoQ) ? extract_column_f32(buf, sfa, cvd[a]) : NULL;
            #pragma omp parallel for
            for (long i = 0; i < N3; i++) {
                float w = u[i]*u[i] + (v ? v[i]*v[i] : 0.f);
                rho2[i] += w;
                sP[i]   *= w;
                if (rca[0]) rca[a][i] = w;
                if (rhoQ) rhoQ[i] += u[i]*dv[i] - v[i]*du[i];
            }
            free(u); if (v) free(v); if (du) free(du); if (dv) free(dv);

            if (th2) {
                float *tu = extract_column_f32(buf, sfa, ctu[a]);
                float *tv = has_thim ? extract_column_f32(buf, sfa, ctv[a]) : NULL;
                float *dtu = (rhoQ && has_thim) ? extract_column_f32(buf, sfa, ctud[a]) : NULL;
                float *dtv = (rhoQ && has_thim) ? extract_column_f32(buf, sfa, ctvd[a]) : NULL;
                #pragma omp parallel for
                for (long i = 0; i < N3; i++) {
                    th2[i] += tu[i]*tu[i] + (tv ? tv[i]*tv[i] : 0.f);
                    if (dtu) rhoQ[i] += tu[i]*dtv[i] - tv[i]*dtu[i];
                }
                free(tu); if (tv) free(tv); if (dtu) free(dtu); if (dtv) free(dtv);
            }
        }
        if (Emag) {
            memset(Emag, 0, N3 * sizeof(float));
            for (int a = 0; a < 3; a++) {
                float *E = extract_column_f32(buf, sfa, cE[a]);
                #pragma omp parallel for
                for (long i = 0; i < N3; i++) Emag[i] += E[i]*E[i];
                free(E);
            }
            #pragma omp parallel for
            for (long i = 0; i < N3; i++) Emag[i] = sqrtf(Emag[i]);
        }

        /* ---- write slice planes ---- */
        float *qsrc[12];
        qsrc[q_rho2] = rho2; qsrc[q_s] = sP;
        if (q_comp >= 0)
            for (int a = 0; a < 3; a++) qsrc[q_comp + a] = rca[a];
        if (q_rhoQ >= 0) qsrc[q_rhoQ] = rhoQ;
        if (q_th2  >= 0) qsrc[q_th2]  = th2;
        if (q_E    >= 0) qsrc[q_E]    = Emag;

        for (int q = 0; q < nq; q++) {
            const float *src = qsrc[q];
            if (axis == 'z') {        /* plane (x,y): n1=x, n2=y */
                #pragma omp parallel for
                for (int iy = 0; iy < N; iy++)
                    for (int ix = 0; ix < N; ix++)
                        plane[(long)iy*N + ix] = src[(long)ix*NN + (long)iy*N + k_slice];
            } else if (axis == 'y') { /* plane (x,z): n1=x, n2=z */
                #pragma omp parallel for
                for (int iz = 0; iz < N; iz++)
                    for (int ix = 0; ix < N; ix++)
                        plane[(long)iz*N + ix] = src[(long)ix*NN + (long)k_slice*N + iz];
            } else {                  /* plane (y,z): n1=y, n2=z */
                #pragma omp parallel for
                for (int iz = 0; iz < N; iz++)
                    for (int iy = 0; iy < N; iy++)
                        plane[(long)iz*N + iy] = src[(long)k_slice*NN + (long)iy*N + iz];
            }
            fwrite(plane, sizeof(float), (long)N * N, fd);
        }

        /* ---- lineout through the rho2 max ---- */
        if (fl) {
            long imax = 0; float rmax = rho2[0];
            for (long i = 1; i < N3; i++)
                if (rho2[i] > rmax) { rmax = rho2[i]; imax = i; }
            int mx = (int)(imax / NN), my = (int)((imax / N) % N), mz = (int)(imax % N);
            float **cols = (float **)malloc(sfa->n_columns * sizeof(float *));
            for (uint32_t c = 0; c < sfa->n_columns; c++)
                cols[c] = extract_column_f32(buf, sfa, (int)c);
            for (int i = 0; i < N; i++) {
                long idx; double coord = -L + i * dx;
                if (line_axis == 'x')      idx = (long)i*NN + (long)my*N + mz;
                else if (line_axis == 'y') idx = (long)mx*NN + (long)i*N + mz;
                else                       idx = (long)mx*NN + (long)my*N + i;
                fprintf(fl, "%d\t%.6f\t%.6f", fi, times[si], coord);
                for (uint32_t c = 0; c < sfa->n_columns; c++)
                    fprintf(fl, "\t%.7e", cols[c][idx]);
                fprintf(fl, "\t%.7e\n", rho2[idx]);
            }
            for (uint32_t c = 0; c < sfa->n_columns; c++) free(cols[c]);
            free(cols);
            fprintf(stderr, "frame %d t=%.2f: rho2_max=%.4g at (%.2f,%.2f,%.2f)\n",
                    fi, times[si], rmax,
                    -L + mx*dx, -L + my*dx, -L + mz*dx);
        }
    }
    fclose(fd);
    if (fl) fclose(fl);

    FILE *fm = fopen(meta_path, "w");
    if (!fm) { fprintf(stderr, "error: cannot write %s\n", meta_path); return 1; }
    fprintf(fm, "input=%s\nN=%d\nL=%.10g\ndx=%.10g\naxis=%c\nslice_index=%d\n"
                "slice_pos=%.10g\nnq=%d\n", in_path, N, L, dx, axis, k_slice,
                -L + k_slice * dx, nq);
    fprintf(fm, "quantities=");
    for (int q = 0; q < nq; q++) fprintf(fm, "%s%s", q ? "," : "", qnames[q]);
    fprintf(fm, "\nnframes=%d\nframes=", nf);
    for (int i = 0; i < nf; i++) fprintf(fm, "%s%d", i ? "," : "", fsel[i]);
    fprintf(fm, "\ntimes=");
    for (int i = 0; i < nf; i++) fprintf(fm, "%s%.6f", i ? "," : "", times[i]);
    fprintf(fm, "\n");
    fclose(fm);

    printf("sfa_slice: wrote %s + %s (%d frames x %d quantities, %dx%d planes)%s%s\n",
           meta_path, data_path, nf, nq, N, N,
           line_axis ? " + " : "", line_axis ? line_path : "");
    sfa_close(sfa);
    return 0;
}
