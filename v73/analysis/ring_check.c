/*  ring_check.c — winding number, L_z, and azimuthal charge current of a
 *  24-column complex SFA frame (kernel output or seed). Verifies that the
 *  Q-ring's real-space circulation survives true dynamics.
 *
 *    L_z  = ∫ Σ_a [ u̇_a ∂_φ u_a + v̇_a ∂_φ v_a ] ... computed as
 *           −∫ Σ_a ( u̇ ∂_φu + v̇ ∂_φv )  with  ∂_φ = x∂_y − y∂_x
 *           (angular momentum of the scalar sector; sign fixed so the
 *           relaxer's rotating ring gives L_z > 0)
 *    winding = (1/2π) ∮ d(arg Φ_0) around the circle of peak |Φ|² (z=0)
 *    J_φ(ρ_pk) = azimuthal charge current at the peak circle
 *
 *  Build: gcc -O3 -fopenmp -o ring_check ring_check.c -lzstd -lm
 *  Usage: ring_check file.sfa [frame_idx=-1 (last)]
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#define SFA_IMPLEMENTATION
#include "../../sfa/format/sfa.h"

int main(int argc, char **argv) {
    if (argc < 2) { fprintf(stderr, "usage: %s file.sfa [frame]\n", argv[0]); return 1; }
    SFA *s = sfa_open(argv[1]);
    if (!s) { fprintf(stderr, "cannot open\n"); return 1; }
    int N = (int)s->Nx; long N3 = (long)N*N*N, NN = (long)N*N;
    double L = s->Lx, dx = 2.0*L/(N-1), dV = dx*dx*dx;
    int fr = (argc > 2) ? atoi(argv[2]) : (int)s->total_frames - 1;
    if (fr < 0) fr = (int)s->total_frames - 1;

    int cu[3],cv[3],cud[3],cvd[3];
    for (int a=0;a<3;a++) cu[a]=cv[a]=cud[a]=cvd[a]=-1;
    for (uint32_t c = 0; c < s->n_columns; c++) {
        int sem = s->columns[c].semantic, cmp = s->columns[c].component;
        if (sem == SFA_POSITION && cmp < 3) cu[cmp] = c;
        else if (sem == SFA_POSITION && cmp < 6) cv[cmp-3] = c;
        else if (sem == SFA_VELOCITY && cmp < 3) cud[cmp] = c;
        else if (sem == SFA_VELOCITY && cmp >= 6 && cmp < 9) cvd[cmp-6] = c;
    }
    void *buf = malloc(s->frame_bytes);
    if (sfa_read_frame(s, fr, buf) != 0) { fprintf(stderr, "read fail\n"); return 1; }
    uint64_t *coff = malloc(s->n_columns * sizeof(uint64_t));
    { uint64_t off = 0; for (uint32_t c = 0; c < s->n_columns; c++) {
        coff[c] = off; off += (uint64_t)N3 * sfa_dtype_size[s->columns[c].dtype]; } }
    #define F(c, ii) ((float*)((uint8_t*)buf + coff[(c)]))[(ii)]
    #define IX(i,j,k) (((long)(i)*N+(j))*N+(k))

    /* peak circle in z=0 plane */
    int kz = N/2; double rho_pk = 0, amp_pk = 0;
    for (int i = 1; i < N-1; i++) for (int j = 1; j < N-1; j++) {
        double x = -L+i*dx, y = -L+j*dx; long ii = IX(i,j,kz);
        double a2 = 0; for (int a=0;a<3;a++) a2 += F(cu[a],ii)*F(cu[a],ii)+F(cv[a],ii)*F(cv[a],ii);
        if (a2 > amp_pk) { amp_pk = a2; rho_pk = sqrt(x*x+y*y); }
    }

    /* L_z */
    double Lz = 0;
    #pragma omp parallel for reduction(+:Lz)
    for (long ii = 0; ii < N3; ii++) {
        int i=(int)(ii/NN), j=(int)((ii/N)%N), k=(int)(ii%N);
        if (i<1||i>=N-1||j<1||j>=N-1||k<1||k>=N-1) continue;
        double x=-L+i*dx, y=-L+j*dx;
        long Xp=IX(i+1,j,k),Xm=IX(i-1,j,k),Yp=IX(i,j+1,k),Ym=IX(i,j-1,k);
        for (int a=0;a<3;a++){
            double dfu = x*(F(cu[a],Yp)-F(cu[a],Ym))/(2*dx) - y*(F(cu[a],Xp)-F(cu[a],Xm))/(2*dx);
            double dfv = x*(F(cv[a],Yp)-F(cv[a],Ym))/(2*dx) - y*(F(cv[a],Xp)-F(cv[a],Xm))/(2*dx);
            /* same orientation convention as eta_qflow's ring report:
             * +n omega A^2 for Phi = A e^{i(omega t + n phi)} */
            Lz += F(cud[a],ii)*dfu + F(cvd[a],ii)*dfv;
        }
    }
    Lz *= dV;

    /* charge and azimuthal current at peak circle */
    double Q = 0;
    #pragma omp parallel for reduction(+:Q)
    for (long ii = 0; ii < N3; ii++)
        for (int a=0;a<3;a++) Q += F(cu[a],ii)*F(cvd[a],ii) - F(cv[a],ii)*F(cud[a],ii);
    Q *= dV;

    int M = 256; double wsum = 0, prev = 0, jphi = 0; int jn = 0;
    for (int t = 0; t <= M; t++) {
        double ph = 2.0*M_PI*t/M, x = rho_pk*cos(ph), y = rho_pk*sin(ph);
        int i = (int)lround((x+L)/dx), j = (int)lround((y+L)/dx);
        if (i<1||i>=N-1||j<1||j>=N-1) continue;
        long ii = IX(i,j,kz);
        double arg = atan2(F(cv[0],ii), F(cu[0],ii));
        if (t > 0) { double d = arg-prev;
            while (d > M_PI) d -= 2*M_PI; while (d < -M_PI) d += 2*M_PI; wsum += d; }
        prev = arg;
        /* J_phi = phihat . J, J_i = sum (v d_i u - u d_i v) */
        long Xp=IX(i+1,j,kz),Xm=IX(i-1,j,kz),Yp=IX(i,j+1,kz),Ym=IX(i,j-1,kz);
        double Jx=0, Jy=0;
        for (int a=0;a<3;a++){
            Jx += F(cv[a],ii)*(F(cu[a],Xp)-F(cu[a],Xm))/(2*dx) - F(cu[a],ii)*(F(cv[a],Xp)-F(cv[a],Xm))/(2*dx);
            Jy += F(cv[a],ii)*(F(cu[a],Yp)-F(cu[a],Ym))/(2*dx) - F(cu[a],ii)*(F(cv[a],Yp)-F(cv[a],Ym))/(2*dx);
        }
        jphi += -sin(ph)*Jx + cos(ph)*Jy; jn++;
    }
    /* poloidal winding: tube center R_c = amp^2-weighted mean rho in the
     * y=0 half-plane (x>0); poloidal circle radius = distance of the peak
     * from (R_c, 0); sample chi around it and sum wrapped phase steps. */
    int jy = N/2;
    double Rc_n = 0, Rc_d = 0, ap2 = 0, prho = 0, pz = 0;
    for (int i = N/2+1; i < N-1; i++) for (int k = 1; k < N-1; k++) {
        double x = -L+i*dx, z = -L+k*dx; long ii = IX(i,jy,k);
        double a2 = 0; for (int a=0;a<3;a++) a2 += F(cu[a],ii)*F(cu[a],ii)+F(cv[a],ii)*F(cv[a],ii);
        Rc_n += a2*x; Rc_d += a2;
        if (a2 > ap2) { ap2 = a2; prho = x; pz = z; }
    }
    double Rc = (Rc_d>0)? Rc_n/Rc_d : 0;
    double ds = sqrt((prho-Rc)*(prho-Rc) + pz*pz);
    double pw = 0; int first = 1; double prevp = 0;
    for (int t = 0; t <= M; t++) {
        double ch = 2.0*M_PI*t/M;
        double x = Rc + ds*cos(ch), z = ds*sin(ch);
        int i = (int)lround((x+L)/dx), k = (int)lround((z+L)/dx);
        if (i<1||i>=N-1||k<1||k>=N-1) continue;
        long ii = IX(i,jy,k);
        double arg = atan2(F(cv[0],ii), F(cu[0],ii));
        if (!first) { double dd = arg-prevp;
            while (dd > M_PI) dd -= 2*M_PI; while (dd < -M_PI) dd += 2*M_PI; pw += dd; }
        prevp = arg; first = 0;
    }
    printf("frame %d (t=%.2f): Q=%.4f  L_z=%.4f  L_z/Q=%.4f  rho_peak=%.3f  winding=%.3f  J_phi(rho_pk)=%.4e\n",
           fr, sfa_frame_time(s, fr), Q, Lz, Lz/Q, rho_pk, wsum/(2*M_PI), jphi/(jn>0?jn:1));
    printf("  poloidal: R_c=%.3f  d_pk=%.3f  twist_winding=%.3f\n", Rc, ds, pw/(2*M_PI));
    sfa_close(s);
    return 0;
}
