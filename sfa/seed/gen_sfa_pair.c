/*  gen_sfa_pair.c — superpose two copies of a 24-column seed SFA at ±D/2.
 *
 *  Purpose (v72): build a two-ball seed from a RELAXED stationary η-soliton
 *  (eta_qflow output) so that pair experiments (clock-interference force,
 *  QFI witness) run on clean eigenstates instead of Θ=0 transients. Each
 *  ball is locally the exact stationary state; the residual dynamics is the
 *  genuine two-body interaction.
 *
 *  field(x) = f(x − s x̂_ax) + R(Δφ)[f](x + s x̂_ax),  s = round(D/2/dx)·dx
 *
 *  R(Δφ) rotates ball 2's global U(1) phase: every complex pair
 *  (re, im) → (re cosΔφ − im sinΔφ, re sinΔφ + im cosΔφ). Pairs are matched
 *  by semantic/component (POSITION a↔a+3 = Φ, ANGLE a↔a+3 = Θ,
 *  VELOCITY a↔a+6 = Φ̇, VELOCITY a+3↔a+9 = Θ̇). The shift is an integer
 *  number of cells (no interpolation); out-of-range cells are zero.
 *
 *  Build: gcc -O3 -fopenmp -o gen_sfa_pair gen_sfa_pair.c -lzstd -lm
 *  Usage: gen_sfa_pair in.sfa out.sfa D [dphi_deg=0] [axis=x]
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#define SFA_IMPLEMENTATION
#include "../format/sfa.h"

int main(int argc, char **argv) {
    if (argc < 4) { fprintf(stderr, "usage: %s in.sfa out.sfa D [dphi_deg=0] [axis=x]\n", argv[0]); return 1; }
    const char *inp = argv[1], *outp = argv[2];
    double D = atof(argv[3]);
    double dphi = (argc > 4) ? atof(argv[4]) * M_PI / 180.0 : 0.0;
    char axis = (argc > 5) ? argv[5][0] : 'x';

    SFA *s = sfa_open(inp);
    if (!s) { fprintf(stderr, "cannot open %s\n", inp); return 1; }
    int N = s->Nx;
    long N3 = (long)N * N * N, NN = (long)N * N;
    double L = s->Lx, dx = 2.0 * L / (N - 1);
    int nc = s->n_columns;
    int sh = (int)lround(0.5 * D / dx);
    printf("gen_sfa_pair: N=%d L=%g dx=%.4f  D=%.3f -> shift %d cells each way (D_actual=%.3f)  dphi=%.1f deg  axis=%c\n",
           N, L, dx, D, sh, 2.0 * sh * dx, dphi * 180.0 / M_PI, axis);

    /* read frame 0, upconvert all columns to double */
    void *buf = malloc(s->frame_bytes);
    if (sfa_read_frame(s, 0, buf) != 0) { fprintf(stderr, "read_frame failed\n"); return 1; }
    double **col = malloc(nc * sizeof(double *));
    uint64_t off = 0;
    for (int c = 0; c < nc; c++) {
        int dtype = s->columns[c].dtype;
        int es = sfa_dtype_size[dtype];
        col[c] = malloc(N3 * sizeof(double));
        uint8_t *src = (uint8_t *)buf + off;
        if (dtype == SFA_F64) memcpy(col[c], src, N3 * 8);
        else if (dtype == SFA_F32) for (long i = 0; i < N3; i++) col[c][i] = (double)((float *)src)[i];
        else { fprintf(stderr, "unsupported dtype in column %d\n", c); return 1; }
        off += (uint64_t)N3 * es;
    }
    free(buf);

    /* find the imaginary partner of each column (or -1) */
    int *partner = malloc(nc * sizeof(int));
    for (int c = 0; c < nc; c++) partner[c] = -1;
    for (int c = 0; c < nc; c++) {
        int sem = s->columns[c].semantic, cmp = s->columns[c].component;
        int want = -1, wc = -1;
        if ((sem == SFA_POSITION || sem == SFA_ANGLE) && cmp < 3) { want = sem; wc = cmp + 3; }
        if (sem == SFA_VELOCITY && cmp < 3)                        { want = sem; wc = cmp + 6; }
        if (sem == SFA_VELOCITY && cmp >= 3 && cmp < 6)            { want = sem; wc = cmp + 6; }
        if (want < 0) continue;
        for (int c2 = 0; c2 < nc; c2++)
            if (s->columns[c2].semantic == want && s->columns[c2].component == wc) { partner[c] = c2; break; }
    }

    /* superpose: out = shift(col, -sh) + R(dphi) shift(col, +sh) along axis */
    float **out = malloc(nc * sizeof(float *));
    for (int c = 0; c < nc; c++) out[c] = calloc((size_t)N3, sizeof(float));
    double cph = cos(dphi), sph = sin(dphi);
    int done = 0;
    for (int c = 0; c < nc; c++) {
        int im = partner[c];
        if (im < 0) continue;               /* handled as someone's partner */
        #pragma omp parallel for
        for (long ii = 0; ii < N3; ii++) {
            int i = (int)(ii / NN), j = (int)((ii / N) % N), k = (int)(ii % N);
            int i1 = i, j1 = j, k1 = k, i2 = i, j2 = j, k2 = k;
            if (axis == 'x') { i1 = i - sh; i2 = i + sh; }
            else if (axis == 'y') { j1 = j - sh; j2 = j + sh; }
            else { k1 = k - sh; k2 = k + sh; }
            double re1 = 0, im1 = 0, re2 = 0, im2 = 0;
            if (i1 >= 0 && i1 < N && j1 >= 0 && j1 < N && k1 >= 0 && k1 < N) {
                long q = ((long)i1 * N + j1) * N + k1;
                re1 = col[c][q]; im1 = col[im][q];
            }
            if (i2 >= 0 && i2 < N && j2 >= 0 && j2 < N && k2 >= 0 && k2 < N) {
                long q = ((long)i2 * N + j2) * N + k2;
                re2 = col[c][q]; im2 = col[im][q];
            }
            out[c][ii]  = (float)(re1 + (re2 * cph - im2 * sph));
            out[im][ii] = (float)(im1 + (re2 * sph + im2 * cph));
        }
        done += 2;
    }
    if (done != nc)
        fprintf(stderr, "WARNING: %d of %d columns paired; unpaired columns left zero\n", done, nc);

    /* write with the same schema */
    SFA *o = sfa_create(outp, N, N, N, L, L, L, s->dt);
    for (int c = 0; c < nc; c++)
        sfa_add_column(o, s->columns[c].name, SFA_F32, s->columns[c].semantic, s->columns[c].component);
    sfa_finalize_header(o);
    void **fc = malloc(nc * sizeof(void *));
    for (int c = 0; c < nc; c++) fc[c] = out[c];
    sfa_write_frame(o, 0.0, fc);
    sfa_close(o);
    sfa_close(s);
    printf("gen_sfa_pair: wrote %s\n", outp);
    return 0;
}
