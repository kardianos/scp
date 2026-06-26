/*  gen_blob_field.c — free-form gauge-connected seed for fitness-driven search
 *
 *  Inverted methodology (FUTURE.md / 2026-06-26): do NOT start from a Lagrangian
 *  ansatz. Represent a candidate particle as an arbitrary superposition of
 *  rotating complex "blobs" — any shape, knotted or not — and let the search
 *  (fitness = shape-retention + QFI) decide what is viable. The blobs are tied
 *  together through the gauge field automatically (Gauss law builds E from the
 *  joint charge density when the kernel runs complex_gauge=1).
 *
 *  Each component field:  Phi_a(x) = sum_b c_{a,b} exp(-|x-r_b|^2 / (2 sigma_b^2))
 *  Rotating at the common clock omega (charge carrier):
 *    Phi_a(t,x) = A_a(x) e^{i omega t}
 *    u_a = Re A_a,  v_a = Im A_a,  udot_a = -omega v_a,  vdot_a = +omega u_a
 *  Theta sector seeded to zero (kernel relaxes / Gauss builds E).
 *
 *  Genome file (text):
 *    omega <w>
 *    blob <cx> <cy> <cz> <sigma> <a0r> <a0i> <a1r> <a1i> <a2r> <a2i>
 *    blob ...
 *
 *  Build: gcc -O3 -fopenmp -o gen_blob_field gen_blob_field.c -lzstd -lm
 *  Usage: gen_blob_field N L genome.txt out.sfa [curl]
 *    optional 5th arg "curl": seed the torsion partner consistently,
 *    Theta_a = (1/2)(curl Phi)_a (complex, rotating at omega), instead of 0.
 *    This makes the configuration satisfy the eta-coupled (Cosserat) equation
 *    at t=0, removing the Theta-from-zero transient that otherwise drives a
 *    ~2% energy drift over long T (FUTURE.md drift-floor finding).
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#define SFA_IMPLEMENTATION
#include "../format/sfa.h"

#define NCOLS 24
#define MAXBLOB 64

int main(int argc, char **argv) {
    if (argc < 5) { fprintf(stderr, "usage: %s N L genome.txt out.sfa\n", argv[0]); return 1; }
    int N = atoi(argv[1]); double L = atof(argv[2]);
    const char *genome = argv[3], *out = argv[4];
    int theta_curl = (argc > 5 && !strcmp(argv[5], "curl"));

    double omega = 1.45;
    double bc[MAXBLOB][3], bsig[MAXBLOB], bamp[MAXBLOB][6]; int nb = 0;
    FILE *fp = fopen(genome, "r"); if (!fp) { fprintf(stderr, "cannot open %s\n", genome); return 1; }
    char line[512];
    while (fgets(line, sizeof(line), fp)) {
        if (line[0]=='#'||line[0]=='\n') continue;
        if (!strncmp(line, "omega", 5)) sscanf(line+5, "%lf", &omega);
        else if (!strncmp(line, "blob", 4)) {
            if (nb >= MAXBLOB) continue;
            double *c=bc[nb], *a=bamp[nb];
            if (sscanf(line+4, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
                &c[0],&c[1],&c[2],&bsig[nb],&a[0],&a[1],&a[2],&a[3],&a[4],&a[5]) == 10) nb++;
        }
    }
    fclose(fp);
    if (nb == 0) { fprintf(stderr, "genome has no blobs\n"); return 1; }

    long N3 = (long)N*N*N; double dx = 2.0*L/(N-1);
    float *cols[NCOLS]; for (int c=0;c<NCOLS;c++) cols[c]=calloc((size_t)N3,sizeof(float));

    #pragma omp parallel for
    for (long ii=0; ii<N3; ii++) {
        int i=(int)(ii/((long)N*N)), j=(int)((ii/N)%N), k=(int)(ii%N);
        double x=-L+i*dx, y=-L+j*dx, z=-L+k*dx;
        double ur[3]={0,0,0}, ui[3]={0,0,0};
        for (int b=0;b<nb;b++) {
            double dxx=x-bc[b][0], dyy=y-bc[b][1], dzz=z-bc[b][2];
            double g=exp(-(dxx*dxx+dyy*dyy+dzz*dzz)/(2.0*bsig[b]*bsig[b]));
            for (int a=0;a<3;a++){ ur[a]+=bamp[b][2*a]*g; ui[a]+=bamp[b][2*a+1]*g; }
        }
        for (int a=0;a<3;a++){
            double u=ur[a], v=ui[a];
            cols[0 + a][ii] = (float)u;            /* phi_x..z   (u_a)   */
            cols[12+ a][ii] = (float)v;            /* phiim_x..z (v_a)   */
            cols[6 + a][ii] = (float)(-omega*v);   /* phi_v   (udot_a)   */
            cols[18+ a][ii] = (float)( omega*u);   /* phiim_v (vdot_a)   */
        }
    }

    if (theta_curl) {
        /* Theta_a = (1/2)(curl Phi)_a via central differences (periodic), complex,
         * rotating at omega: theta_v = i*omega*theta. Cols: theta re=3+a, im=15+a,
         * theta_v re=9+a, im=21+a. Operates on phi u=cols[0..2], v=cols[12..14]. */
        double inv2dx = 1.0/(2.0*dx);
        #define IDX(i,j,k) ((((long)(i)*N+(j))*N)+(k))
        #pragma omp parallel for
        for (long ii=0; ii<N3; ii++) {
            int i=(int)(ii/((long)N*N)), j=(int)((ii/N)%N), k=(int)(ii%N);
            int ip=(i+1)%N, im=(i-1+N)%N, jp=(j+1)%N, jm=(j-1+N)%N, kp=(k+1)%N, km=(k-1+N)%N;
            for (int part=0; part<2; part++) {
                float *F0 = cols[(part?12:0)+0], *F1 = cols[(part?12:0)+1], *F2 = cols[(part?12:0)+2];
                double dxF1=(F1[IDX(ip,j,k)]-F1[IDX(im,j,k)])*inv2dx;
                double dxF2=(F2[IDX(ip,j,k)]-F2[IDX(im,j,k)])*inv2dx;
                double dyF0=(F0[IDX(i,jp,k)]-F0[IDX(i,jm,k)])*inv2dx;
                double dyF2=(F2[IDX(i,jp,k)]-F2[IDX(i,jm,k)])*inv2dx;
                double dzF0=(F0[IDX(i,j,kp)]-F0[IDX(i,j,km)])*inv2dx;
                double dzF1=(F1[IDX(i,j,kp)]-F1[IDX(i,j,km)])*inv2dx;
                double c0=0.5*(dyF2-dzF1), c1=0.5*(dzF0-dxF2), c2=0.5*(dxF1-dyF0);
                int re = part?15:3;       /* theta re/im base col */
                cols[re+0][ii]=(float)c0; cols[re+1][ii]=(float)c1; cols[re+2][ii]=(float)c2;
            }
        }
        /* theta velocity = i*omega*Theta  (vel_re=-omega*Im, vel_im=+omega*Re) */
        #pragma omp parallel for
        for (long ii=0; ii<N3; ii++)
            for (int a=0;a<3;a++) {
                double tre=cols[3+a][ii], tim=cols[15+a][ii];
                cols[9 +a][ii]=(float)(-omega*tim);
                cols[21+a][ii]=(float)( omega*tre);
            }
        #undef IDX
    }

    double dt = 0.025 * dx;
    SFA *sfa = sfa_create(out, N, N, N, L, L, L, dt);
    sfa_add_column(sfa,"phi_x",SFA_F32,SFA_POSITION,0);  sfa_add_column(sfa,"phi_y",SFA_F32,SFA_POSITION,1);  sfa_add_column(sfa,"phi_z",SFA_F32,SFA_POSITION,2);
    sfa_add_column(sfa,"theta_x",SFA_F32,SFA_ANGLE,0);   sfa_add_column(sfa,"theta_y",SFA_F32,SFA_ANGLE,1);   sfa_add_column(sfa,"theta_z",SFA_F32,SFA_ANGLE,2);
    sfa_add_column(sfa,"phi_vx",SFA_F32,SFA_VELOCITY,0); sfa_add_column(sfa,"phi_vy",SFA_F32,SFA_VELOCITY,1); sfa_add_column(sfa,"phi_vz",SFA_F32,SFA_VELOCITY,2);
    sfa_add_column(sfa,"theta_vx",SFA_F32,SFA_VELOCITY,3);sfa_add_column(sfa,"theta_vy",SFA_F32,SFA_VELOCITY,4);sfa_add_column(sfa,"theta_vz",SFA_F32,SFA_VELOCITY,5);
    sfa_add_column(sfa,"phiim_x",SFA_F32,SFA_POSITION,3);sfa_add_column(sfa,"phiim_y",SFA_F32,SFA_POSITION,4);sfa_add_column(sfa,"phiim_z",SFA_F32,SFA_POSITION,5);
    sfa_add_column(sfa,"thetaim_x",SFA_F32,SFA_ANGLE,3); sfa_add_column(sfa,"thetaim_y",SFA_F32,SFA_ANGLE,4); sfa_add_column(sfa,"thetaim_z",SFA_F32,SFA_ANGLE,5);
    sfa_add_column(sfa,"phiim_vx",SFA_F32,SFA_VELOCITY,6);sfa_add_column(sfa,"phiim_vy",SFA_F32,SFA_VELOCITY,7);sfa_add_column(sfa,"phiim_vz",SFA_F32,SFA_VELOCITY,8);
    sfa_add_column(sfa,"thetaim_vx",SFA_F32,SFA_VELOCITY,9);sfa_add_column(sfa,"thetaim_vy",SFA_F32,SFA_VELOCITY,10);sfa_add_column(sfa,"thetaim_vz",SFA_F32,SFA_VELOCITY,11);
    sfa_finalize_header(sfa);
    void *fc[NCOLS]; for(int c=0;c<NCOLS;c++) fc[c]=cols[c];
    sfa_write_frame(sfa, 0.0, fc);
    sfa_close(sfa);

    /* report total charge proxy and core triple-product */
    double Q=0, smax=0;
    for (long ii=0; ii<N3; ii++){ double s=1;
        for(int a=0;a<3;a++){ double u=cols[0+a][ii],v=cols[12+a][ii]; double ud=cols[6+a][ii],vd=cols[18+a][ii];
            Q += u*vd - v*ud; s *= (u*u+v*v); }
        if(s>smax)smax=s; }
    double dV=dx*dx*dx;
    printf("gen_blob_field: %s  N=%d L=%g blobs=%d omega=%.3f  Q~%.2f s_max=%.3e\n",
           out, N, L, nb, omega, Q*dV, smax);
    for (int c=0;c<NCOLS;c++) free(cols[c]);
    return 0;
}
