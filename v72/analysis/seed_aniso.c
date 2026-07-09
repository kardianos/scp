/*  seed_aniso.c — quadrupole/anisotropy diagnostics of a 24-column seed SFA.
 *
 *  The stationary η-coupled ball is predicted to be axisymmetric about
 *  û = (1,1,1)/√3 (the Θ texture g(r)(r̂×û) is toroidal about û and its
 *  back-reaction deforms Φ along û). This tool quantifies that deformation:
 *
 *    weights:  s(x) = Π_a |Φ_a|²  (binding density)  and  ρ_Θ = Σ_a |Θ_a|²
 *    Q20[w]   = ⟨3cos²ϑ − 1⟩_w,  cosϑ = (x·û)/r     (weighted quadrupole)
 *    r_par    = sqrt(⟨(x·û)²⟩_w),  r_perp = sqrt(⟨|x−(x·û)û|²⟩_w / 2)
 *    aspect   = r_par / r_perp    (1 = spherical)
 *
 *  The same numbers with û replaced by ẑ act as a control: for a body
 *  axisymmetric about û, Q20 measured about an axis at angle α to the true
 *  symmetry axis scales by P2(cos α), and cos(û,ẑ) = 1/√3 is exactly the
 *  P2 = 0 magic angle — so a genuinely û-locked quadrupole gives
 *  Q20(ẑ) = 0 identically (measured: 0.00000), not merely smaller.
 *
 *  Build: gcc -O3 -o seed_aniso seed_aniso.c -lzstd -lm
 *  Usage: seed_aniso seed.sfa
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#define SFA_IMPLEMENTATION
#include "../../sfa/format/sfa.h"

int main(int argc, char **argv) {
    if (argc < 2) { fprintf(stderr, "usage: %s seed.sfa\n", argv[0]); return 1; }
    SFA *s = sfa_open(argv[1]);
    if (!s) { fprintf(stderr, "cannot open %s\n", argv[1]); return 1; }
    int N = s->Nx;
    long N3 = (long)N * N * N;
    double L = s->Lx, dx = 2.0 * L / (N - 1);

    void *buf = malloc(s->frame_bytes);
    if (sfa_read_frame(s, 0, buf) != 0) { fprintf(stderr, "read_frame failed\n"); return 1; }

    /* extract by semantic (same pattern as analyze_sfa.c) */
    double *u[3] = {0}, *v[3] = {0}, *tu[3] = {0}, *tv[3] = {0};
    uint64_t off = 0;
    for (uint32_t c = 0; c < s->n_columns; c++) {
        int dtype = s->columns[c].dtype, sem = s->columns[c].semantic, comp = s->columns[c].component;
        int es = sfa_dtype_size[dtype];
        double *arr = malloc(N3 * sizeof(double));
        uint8_t *src = (uint8_t *)buf + off;
        if (dtype == SFA_F64) memcpy(arr, src, N3 * 8);
        else if (dtype == SFA_F32) for (long i = 0; i < N3; i++) arr[i] = (double)((float *)src)[i];
        else { free(arr); arr = NULL; }
        if      (arr && sem == SFA_POSITION && comp < 3)              u[comp] = arr;
        else if (arr && sem == SFA_POSITION && comp >= 3 && comp < 6) v[comp-3] = arr;
        else if (arr && sem == SFA_ANGLE && comp < 3)                 tu[comp] = arr;
        else if (arr && sem == SFA_ANGLE && comp >= 3 && comp < 6)    tv[comp-3] = arr;
        else free(arr);
        off += (uint64_t)N3 * es;
    }
    for (int a = 0; a < 3; a++)
        if (!u[a] || !v[a] || !tu[a] || !tv[a]) {
            fprintf(stderr, "missing field columns (24-col complex seed required)\n");
            return 1;
        }

    const double inv3 = 1.0 / sqrt(3.0);
    double axes[2][3] = {{inv3, inv3, inv3}, {0, 0, 1}};
    const char *axname[2] = {"u=(1,1,1)/sqrt3", "z (control)"};

    for (int w = 0; w < 2; w++) {            /* weight: 0 = s(x), 1 = |Theta|^2 */
        for (int ax = 0; ax < 2; ax++) {
            double W = 0, q20 = 0, par2 = 0, perp2 = 0;
            for (long ii = 0; ii < N3; ii++) {
                int i = (int)(ii / ((long)N * N)), j = (int)((ii / N) % N), k = (int)(ii % N);
                double x = -L + i * dx, y = -L + j * dx, z = -L + k * dx;
                double r2 = x * x + y * y + z * z;
                if (r2 < 1e-12) continue;
                double wt;
                if (w == 0) {
                    double r0 = u[0][ii]*u[0][ii] + v[0][ii]*v[0][ii];
                    double r1 = u[1][ii]*u[1][ii] + v[1][ii]*v[1][ii];
                    double rr = u[2][ii]*u[2][ii] + v[2][ii]*v[2][ii];
                    wt = r0 * r1 * rr;
                } else {
                    wt = 0;
                    for (int a = 0; a < 3; a++)
                        wt += tu[a][ii]*tu[a][ii] + tv[a][ii]*tv[a][ii];
                }
                if (wt <= 0) continue;
                double dpar = x * axes[ax][0] + y * axes[ax][1] + z * axes[ax][2];
                double c2 = dpar * dpar / r2;
                W += wt;
                q20 += wt * (3.0 * c2 - 1.0);
                par2 += wt * dpar * dpar;
                perp2 += wt * (r2 - dpar * dpar);
            }
            if (W > 0) {
                double rpar = sqrt(par2 / W), rperp = sqrt(perp2 / (2.0 * W));
                printf("%-8s axis %-18s  Q20=%+9.5f  r_par=%.4f  r_perp=%.4f  aspect=%.5f\n",
                       w == 0 ? "s(x)" : "|Theta|2", axname[ax], q20 / W, rpar, rperp, rpar / rperp);
            }
        }
    }
    sfa_close(s);
    return 0;
}
