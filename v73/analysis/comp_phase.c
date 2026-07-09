/*  comp_phase.c — per-component amplitude and phase structure of a 24-col
 *  SFA frame: at the amplitude peak and averaged over the core (|Phi|^2 >
 *  half max), print |Phi_a|, arg(Phi_a), and the relative phases — the
 *  internal-polarization (spin) state of a ball. Also the spin prediction
 *  l_spin = int 2[u_0 v_1 - u_1 v_0] + cyc? No: Jhat spin part about z acts
 *  on (0,1) only: l_spin = int 2(u_0 v_1 - v_0 u_1)  [checked vs eta_qflow].
 *
 *  Build: gcc -O3 -fopenmp -o comp_phase comp_phase.c -lzstd -lm
 *  Usage: comp_phase file.sfa [frame=-1]
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define SFA_IMPLEMENTATION
#include "../../sfa/format/sfa.h"

int main(int argc, char **argv) {
    if (argc < 2) { fprintf(stderr, "usage: %s file.sfa [frame]\n", argv[0]); return 1; }
    SFA *s = sfa_open(argv[1]);
    if (!s) { fprintf(stderr, "cannot open\n"); return 1; }
    int N = (int)s->Nx; long N3 = (long)N*N*N;
    double L = s->Lx, dx = 2.0*L/(N-1), dV = dx*dx*dx;
    int fr = (argc > 2) ? atoi(argv[2]) : (int)s->total_frames - 1;
    if (fr < 0) fr = (int)s->total_frames - 1;
    int cu[3], cv[3];
    for (uint32_t c = 0; c < s->n_columns; c++) {
        int sem = s->columns[c].semantic, cmp = s->columns[c].component;
        if (sem == SFA_POSITION && cmp < 3) cu[cmp] = c;
        else if (sem == SFA_POSITION && cmp < 6) cv[cmp-3] = c;
    }
    void *buf = malloc(s->frame_bytes);
    if (sfa_read_frame(s, fr, buf) != 0) { fprintf(stderr, "read fail\n"); return 1; }
    uint64_t *coff = malloc(s->n_columns * sizeof(uint64_t));
    { uint64_t off = 0; for (uint32_t c = 0; c < s->n_columns; c++) {
        coff[c] = off; off += (uint64_t)N3 * sfa_dtype_size[s->columns[c].dtype]; } }
    #define F(c, ii) ((float*)((uint8_t*)buf + coff[(c)]))[(ii)]

    /* peak */
    long pk = 0; double apk = 0;
    for (long ii = 0; ii < N3; ii++) {
        double a2 = 0; for (int a=0;a<3;a++) a2 += F(cu[a],ii)*F(cu[a],ii)+F(cv[a],ii)*F(cv[a],ii);
        if (a2 > apk) { apk = a2; pk = ii; }
    }
    printf("frame %d (t=%.2f): peak |Phi|^2=%.4e\n", fr, sfa_frame_time(s, fr), apk);
    printf("  at peak:   ");
    double ph[3];
    for (int a=0;a<3;a++) { ph[a]=atan2(F(cv[a],pk),F(cu[a],pk));
        printf("|P%d|=%.4f arg=%+7.2f  ", a, sqrt(F(cu[a],pk)*F(cu[a],pk)+F(cv[a],pk)*F(cv[a],pk)), ph[a]*180/M_PI); }
    printf("\n  rel phases: a10=%+7.2f  a20=%+7.2f  a21=%+7.2f deg\n",
        (ph[1]-ph[0])*180/M_PI, (ph[2]-ph[0])*180/M_PI, (ph[2]-ph[1])*180/M_PI);

    /* core-averaged relative-phase coherences + spin integrals */
    double c10r=0,c10i=0,c20r=0,c20i=0, lz01=0, lz12=0, lz20=0, Nn=0;
    #pragma omp parallel for reduction(+:c10r,c10i,c20r,c20i,lz01,lz12,lz20,Nn)
    for (long ii = 0; ii < N3; ii++) {
        double u0=F(cu[0],ii),v0=F(cv[0],ii),u1=F(cu[1],ii),v1=F(cv[1],ii),u2=F(cu[2],ii),v2=F(cv[2],ii);
        /* coherence <Phi_1 conj(Phi_0)> etc. */
        c10r += u1*u0+v1*v0;  c10i += v1*u0-u1*v0;
        c20r += u2*u0+v2*v0;  c20i += v2*u0-u2*v0;
        lz01 += 2.0*(u0*v1 - v0*u1);   /* spin about z: (0,1) pair */
        lz12 += 2.0*(u1*v2 - v1*u2);   /* spin about x */
        lz20 += 2.0*(u2*v0 - v2*u0);   /* spin about y */
        Nn += u0*u0+v0*v0+u1*u1+v1*v1+u2*u2+v2*v2;
    }
    printf("  core coherences: arg<P1 P0*>=%+7.2f deg  arg<P2 P0*>=%+7.2f deg\n",
        atan2(c10i,c10r)*180/M_PI, atan2(c20i,c20r)*180/M_PI);
    printf("  spin integrals (x dV): l_z(01)=%.4f  l_x(12)=%.4f  l_y(20)=%.4f   N=%.4f\n",
        lz01*dV, lz12*dV, lz20*dV, Nn*dV);
    sfa_close(s);
    return 0;
}
