/*  magneton.c — magnetic moment and g-factor of a gauged SCP object from a
 *  30-column kernel SFA frame (matter 24 cols + links th_* ANGLE 6-8 +
 *  electric field E_* VELOCITY 12-14; v69 SPEC conventions:
 *  theta_i = g*dx*A_i on links, B = Theta_plaq/(g*dx^2), E stored physical).
 *
 *  Outputs per frame:
 *    Q            = int sum_a (u vdot - v udot)                (matter charge)
 *    L_z          = int sum_a (udot dphi_u + vdot dphi_v)      (orbital ang mom)
 *    mu_z         = (1/2) int (x J_y - y J_x),                 (magnetic moment
 *                   J_i = sum_a [u d_i v - v d_i u] - (th_i/dx)|Phi|^2
 *                   — physical Im[Phibar D_i Phi] convention;   in units of g)
 *    g_factor     = mu_z / ( Q * L_z / (2M) )   (M = total energy, CLI arg;
 *                   classical rigid rotor with rho_Q ∝ rho_E gives 1)
 *    Phi_B(hole)  = int B_z dA over z=0 disk rho < rho_in      (trapped flux)
 *    B_z(0)       = field at ring center;  B_z profile on the z=0 axis line
 *    E_B, E_E     = compact magnetic / electric field energies (sum = kernel
 *                   E_em diagnostic — calibration check)
 *
 *  Build: gcc -O3 -fopenmp -o magneton magneton.c -lzstd -lm
 *  Usage: magneton file.sfa g_gauge M [frame=-1] [rho_in=2.0]
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define SFA_IMPLEMENTATION
#include "../../sfa/format/sfa.h"

int main(int argc, char **argv) {
    if (argc < 4) { fprintf(stderr, "usage: %s file.sfa g_gauge M [frame] [rho_in]\n", argv[0]); return 1; }
    SFA *s = sfa_open(argv[1]);
    if (!s) { fprintf(stderr, "cannot open\n"); return 1; }
    double gg = atof(argv[2]), M = atof(argv[3]);
    int N = (int)s->Nx; long N3 = (long)N*N*N, NN = (long)N*N;
    double L = s->Lx, dx = 2.0*L/(N-1), dV = dx*dx*dx;
    int fr = (argc > 4) ? atoi(argv[4]) : (int)s->total_frames - 1;
    if (fr < 0) fr = (int)s->total_frames - 1;
    double rho_in = (argc > 5) ? atof(argv[5]) : 2.0;

    int cu[3],cv[3],cud[3],cvd[3],cth[3],cE[3];
    for (int a=0;a<3;a++) cu[a]=cv[a]=cud[a]=cvd[a]=cth[a]=cE[a]=-1;
    for (uint32_t c = 0; c < s->n_columns; c++) {
        int sem = s->columns[c].semantic, cmp = s->columns[c].component;
        if      (sem == SFA_POSITION && cmp < 3)  cu[cmp] = c;
        else if (sem == SFA_POSITION && cmp < 6)  cv[cmp-3] = c;
        else if (sem == SFA_VELOCITY && cmp < 3)  cud[cmp] = c;
        else if (sem == SFA_VELOCITY && cmp >= 6 && cmp < 9)  cvd[cmp-6] = c;
        else if (sem == SFA_ANGLE    && cmp >= 6 && cmp < 9)  cth[cmp-6] = c;
        else if (sem == SFA_VELOCITY && cmp >= 12 && cmp < 15) cE[cmp-12] = c;
    }
    if (cth[0] < 0) { fprintf(stderr, "no gauge columns (need complex_gauge output)\n"); return 1; }
    void *buf = malloc(s->frame_bytes);
    if (sfa_read_frame(s, fr, buf) != 0) { fprintf(stderr, "read fail\n"); return 1; }
    uint64_t *coff = malloc(s->n_columns * sizeof(uint64_t));
    { uint64_t off = 0; for (uint32_t c = 0; c < s->n_columns; c++) {
        coff[c] = off; off += (uint64_t)N3 * sfa_dtype_size[s->columns[c].dtype]; } }
    #define F(c, ii) ((float*)((uint8_t*)buf + coff[(c)]))[(ii)]
    #define IX(i,j,k) (((long)(i)*N+(j))*N+(k))
    double inv2dx = 1.0/(2.0*dx), bpre = 1.0/(gg*gg*dx*dx*dx*dx), binv = 1.0/(gg*dx*dx);

    double Q=0, Lz=0, muz=0, EB=0, EE=0;
    #pragma omp parallel for reduction(+:Q,Lz,muz,EB,EE)
    for (long ii = 0; ii < N3; ii++) {
        int i=(int)(ii/NN), j=(int)((ii/N)%N), k=(int)(ii%N);
        for (int a=0;a<3;a++) Q += F(cu[a],ii)*F(cvd[a],ii) - F(cv[a],ii)*F(cud[a],ii);
        EE += 0.5*(F(cE[0],ii)*F(cE[0],ii)+F(cE[1],ii)*F(cE[1],ii)+F(cE[2],ii)*F(cE[2],ii));
        if (i<1||i>=N-1||j<1||j>=N-1||k<1||k>=N-1) continue;
        double x=-L+i*dx, y=-L+j*dx;
        long Xp=IX(i+1,j,k),Xm=IX(i-1,j,k),Yp=IX(i,j+1,k),Ym=IX(i,j-1,k);
        /* compact magnetic energy: plaquettes (a,b) = (x,y),(y,z),(z,x) */
        long np[3] = {Xp, Yp, IX(i,j,k+1)};
        const int pa[3]={0,1,2}, pb[3]={1,2,0};
        for (int p=0;p<3;p++) {
            int a=pa[p], b=pb[p];
            double ang = F(cth[a],ii) + F(cth[b],np[a]) - F(cth[a],np[b]) - F(cth[b],ii);
            EB += bpre*(1.0 - cos(ang));
        }
        double amp2 = 0, Jx = 0, Jy = 0;
        for (int a=0;a<3;a++){
            double u=F(cu[a],ii), v=F(cv[a],ii);
            amp2 += u*u + v*v;
            Jx += u*(F(cv[a],Xp)-F(cv[a],Xm))*inv2dx - v*(F(cu[a],Xp)-F(cu[a],Xm))*inv2dx;
            Jy += u*(F(cv[a],Yp)-F(cv[a],Ym))*inv2dx - v*(F(cu[a],Yp)-F(cu[a],Ym))*inv2dx;
            double dfu = x*(F(cu[a],Yp)-F(cu[a],Ym))*inv2dx - y*(F(cu[a],Xp)-F(cu[a],Xm))*inv2dx;
            double dfv = x*(F(cv[a],Yp)-F(cv[a],Ym))*inv2dx - y*(F(cv[a],Xp)-F(cv[a],Xm))*inv2dx;
            Lz += F(cud[a],ii)*dfu + F(cvd[a],ii)*dfv;
        }
        Jx -= (F(cth[0],ii)/dx)*amp2;
        Jy -= (F(cth[1],ii)/dx)*amp2;
        muz += 0.5*(x*Jy - y*Jx);
    }
    Q*=dV; Lz*=dV; muz*=dV; EB*=dV; EE*=dV;

    /* B_z on the z=0 plane: flux through the hole and radial profile */
    int kz = N/2;
    double fluxin = 0, fluxall = 0, B0 = 0;
    printf("frame %d (t=%.2f):\n", fr, sfa_frame_time(s, fr));
    for (int i = 1; i < N-1; i++) for (int j = 1; j < N-1; j++) {
        double x=-L+i*dx, y=-L+j*dx, rho=sqrt(x*x+y*y);
        long ii = IX(i,j,kz);
        double ang = F(cth[0],ii) + F(cth[1],IX(i+1,j,kz)) - F(cth[0],IX(i,j+1,kz)) - F(cth[1],ii);
        double Bz = ang*binv;
        fluxall += Bz*dx*dx;
        if (rho < rho_in) fluxin += Bz*dx*dx;
        if (fabs(x) < 0.5*dx && fabs(y) < 0.5*dx) B0 = Bz;
    }
    printf("  Q=%.4f  L_z=%.4f  M(arg)=%.2f\n", Q, Lz, M);
    printf("  mu_z=%.6f (units of g)   mu*=Q·L_z/2M=%.6f   g_factor=%.4f\n",
           muz, Q*Lz/(2.0*M), muz/(Q*Lz/(2.0*M)));
    printf("  B_z(0)=%.4e  Flux(rho<%.1f)=%.6f  Flux(z=0 plane)=%.4e   [units of g]\n",
           B0*(1.0/1.0), rho_in, fluxin, fluxall);
    printf("  E_B=%.4f  E_E=%.4f  E_B+E_E=%.4f  (kernel E_em diag check)\n", EB, EE, EB+EE);
    /* radial B_z profile at z=0 along +x */
    printf("  B_z(rho) z=0: ");
    for (double rr = 0; rr <= 8.01; rr += 1.0) {
        int i = (int)lround((rr+L)/dx), j = N/2;
        if (i<1||i>=N-1) continue;
        long ii = IX(i,j,kz);
        double ang = F(cth[0],ii) + F(cth[1],IX(i+1,j,kz)) - F(cth[0],IX(i,j+1,kz)) - F(cth[1],ii);
        printf("%.1f:%+.3e ", rr, ang*binv);
    }
    printf("\n");
    sfa_close(s);
    return 0;
}
