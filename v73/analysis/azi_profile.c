/*  azi_profile.c — azimuthal structure of a 24-col SFA frame at the peak
 *  circle (z=0): |Phi|^2(phi), arg Phi_0(phi), and the rotating/standing
 *  decomposition of the n=1 azimuthal mode:
 *      u_a ~ A_u cos(phi+a_u),  v_a ~ A_v cos(phi+a_v)
 *  A rotating ring has A_u = A_v with a_v = a_u - pi/2; a standing wave
 *  has A_v ≈ 0 (or phases aligned). Prints the mode amplitudes.
 *
 *  Build: gcc -O3 -o azi_profile azi_profile.c -lzstd -lm
 *  Usage: azi_profile file.sfa [frame=-1] [rho=auto]
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define SFA_IMPLEMENTATION
#include "../../sfa/format/sfa.h"

int main(int argc, char **argv) {
    if (argc < 2) { fprintf(stderr, "usage: %s file.sfa [frame] [rho]\n", argv[0]); return 1; }
    SFA *s = sfa_open(argv[1]);
    if (!s) { fprintf(stderr, "cannot open\n"); return 1; }
    int N = (int)s->Nx; long N3 = (long)N*N*N;
    double L = s->Lx, dx = 2.0*L/(N-1);
    int fr = (argc > 2) ? atoi(argv[2]) : (int)s->total_frames - 1;
    if (fr < 0) fr = (int)s->total_frames - 1;
    int cu[3],cv[3];
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
    #define IX(i,j,k) (((long)(i)*N+(j))*N+(k))
    int kz = N/2;
    double rho_pk = 0, amp_pk = 0;
    for (int i = 1; i < N-1; i++) for (int j = 1; j < N-1; j++) {
        double x = -L+i*dx, y = -L+j*dx; long ii = IX(i,j,kz);
        double a2 = 0; for (int a=0;a<3;a++) a2 += F(cu[a],ii)*F(cu[a],ii)+F(cv[a],ii)*F(cv[a],ii);
        if (a2 > amp_pk) { amp_pk = a2; rho_pk = sqrt(x*x+y*y); }
    }
    if (argc > 3) rho_pk = atof(argv[3]);
    printf("# frame %d t=%.2f rho=%.3f  (phi_deg |Phi|^2 argPhi0 u0 v0)\n",
           fr, sfa_frame_time(s, fr), rho_pk);
    int M = 32;
    double uc=0,us=0,vc=0,vs=0, u0m=0,v0m=0;
    for (int t = 0; t < M; t++) {
        double ph = 2.0*M_PI*t/M, x = rho_pk*cos(ph), y = rho_pk*sin(ph);
        int i = (int)lround((x+L)/dx), j = (int)lround((y+L)/dx);
        if (i<1||i>=N-1||j<1||j>=N-1) continue;
        long ii = IX(i,j,kz);
        double a2 = 0; for (int a=0;a<3;a++) a2 += F(cu[a],ii)*F(cu[a],ii)+F(cv[a],ii)*F(cv[a],ii);
        double u0 = F(cu[0],ii), v0 = F(cv[0],ii);
        printf("%7.1f  %.5e  %+8.4f  %+.5e  %+.5e\n", ph*180/M_PI, a2, atan2(v0,u0), u0, v0);
        uc += u0*cos(ph)*2.0/M; us += u0*sin(ph)*2.0/M;
        vc += v0*cos(ph)*2.0/M; vs += v0*sin(ph)*2.0/M;
        u0m += fabs(u0)/M; v0m += fabs(v0)/M;
    }
    double Au = sqrt(uc*uc+us*us), Av = sqrt(vc*vc+vs*vs);
    printf("# n=1 mode: A_u=%.5e (phase %.1f deg)  A_v=%.5e (phase %.1f deg)\n",
           Au, atan2(-us,uc)*180/M_PI, Av, atan2(-vs,vc)*180/M_PI);
    printf("# rotating ring: A_v/A_u=1, dphase=90; standing: A_v/A_u->0 or dphase->0\n");
    printf("# A_v/A_u = %.4f   mean|u0|=%.4e mean|v0|=%.4e\n", Av/(Au>0?Au:1), u0m, v0m);
    sfa_close(s);
    return 0;
}
