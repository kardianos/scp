/*  gen_3braid.c — Generate 3 perpendicular braids (xyz) as a single-frame SFA
 *
 *  Three braids along x, y, z axes, each with configurable chirality.
 *  For UUD: two same-chirality + one flipped. For UDD: one same + two flipped.
 *
 *  Build: gcc -O3 -o gen_3braid gen_3braid.c -lzstd -lm
 *  Usage: ./gen_3braid -N 128 -L 15 -chirality UUD -A 0.4 -o seed.sfa
 */

#define SFA_IMPLEMENTATION
#include "../../sfa/format/sfa.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define PI 3.14159265358979323846
#define NFIELDS 3

static inline uint16_t f64_to_f16(double v) {
    float f=(float)v; uint32_t x; memcpy(&x,&f,4);
    uint16_t sign=(x>>16)&0x8000; int exp=((x>>23)&0xFF)-127+15; uint16_t mant=(x>>13)&0x3FF;
    if(exp<=0)return sign; if(exp>=31)return sign|0x7C00; return sign|(exp<<10)|mant;
}

int main(int argc, char **argv) {
    int N = 128;
    double L = 15.0, A = 0.4, A_bg = 0.1;
    double R_tube = 3.0, ellip = 0.3325;
    double delta[3] = {0.0, 3.0005, 4.4325};
    double m2 = 2.25;
    double dt_factor = 0.025;
    char chirality[8] = "UUU";  /* U=up, D=down per axis */
    char outpath[512] = "3braid_seed.sfa";
    int precision = 1;

    for (int i=1; i<argc-1; i+=2) {
        if      (!strcmp(argv[i],"-N"))         N = atoi(argv[i+1]);
        else if (!strcmp(argv[i],"-L"))         L = atof(argv[i+1]);
        else if (!strcmp(argv[i],"-A"))         A = atof(argv[i+1]);
        else if (!strcmp(argv[i],"-A_bg"))      A_bg = atof(argv[i+1]);
        else if (!strcmp(argv[i],"-R"))         R_tube = atof(argv[i+1]);
        else if (!strcmp(argv[i],"-ellip"))     ellip = atof(argv[i+1]);
        else if (!strcmp(argv[i],"-chirality")) strncpy(chirality, argv[i+1], 7);
        else if (!strcmp(argv[i],"-o"))         strncpy(outpath, argv[i+1], 511);
        else if (!strcmp(argv[i],"-precision")) {
            if(!strcmp(argv[i+1],"f16")) precision=0;
            else if(!strcmp(argv[i+1],"f32")) precision=1;
            else if(!strcmp(argv[i+1],"f64")) precision=2;
        }
    }

    long N3 = (long)N*N*N;
    double dx = 2.0*L/(N-1);
    double dt = dt_factor * dx;
    int NN = N*N;

    /* Chirality: +1 = up (standard k direction), -1 = down (flipped) */
    double chiral[3] = {1, 1, 1};
    for (int b=0; b<3 && chirality[b]; b++)
        chiral[b] = (chirality[b]=='D' || chirality[b]=='d') ? -1.0 : 1.0;

    fprintf(stderr, "gen_3braid: N=%d L=%.1f A=%.2f R=%.1f chirality=%s\n",
            N, L, A, R_tube, chirality);

    double *phi[3], *phi_vel[3], *theta[3], *theta_vel[3];
    for (int a=0; a<3; a++) {
        phi[a]       = (double*)calloc(N3, sizeof(double));
        phi_vel[a]   = (double*)calloc(N3, sizeof(double));
        theta[a]     = (double*)calloc(N3, sizeof(double));
        theta_vel[a] = (double*)calloc(N3, sizeof(double));
    }

    double kw = PI/L;
    double sx = 1+ellip, sy = 1-ellip;
    double inv2R2 = 1.0/(2*R_tube*R_tube);

    /* Background field (along z) */
    double k_bg = PI/L, omega_bg = sqrt(k_bg*k_bg + m2);

    for (int i=0; i<N; i++) { double x=-L+i*dx;
    for (int j=0; j<N; j++) { double y=-L+j*dx;
    for (int k=0; k<N; k++) { double z=-L+k*dx;
        long idx = (long)i*NN + j*N + k;

        /* Background */
        for (int a=0; a<NFIELDS; a++) {
            double ph_bg = k_bg*z + 2*PI*a/3.0;
            phi[a][idx] = A_bg*cos(ph_bg);
            phi_vel[a][idx] = omega_bg*A_bg*sin(ph_bg);
        }

        /* Braid 1: along z-axis (standard orientation) */
        {
            double r2e = x*x/(sx*sx) + y*y/(sy*sy);
            double env = exp(-r2e * inv2R2);
            double omega = sqrt(kw*kw + m2);
            for (int a=0; a<NFIELDS; a++) {
                double ph = chiral[0] * kw * z + delta[a];
                phi[a][idx]     += A*env*cos(ph);
                phi_vel[a][idx] += chiral[0]*omega*A*env*sin(ph);
            }
        }

        /* Braid 2: along x-axis (rotate z→x via y-axis 90° rotation) */
        {
            /* In rotated frame: zr=x, xr=-z, yr=y */
            double zr = x, xr = -z, yr = y;
            double r2e = xr*xr/(sx*sx) + yr*yr/(sy*sy);
            double env = exp(-r2e * inv2R2);
            double omega = sqrt(kw*kw + m2);
            for (int a=0; a<NFIELDS; a++) {
                double ph = chiral[1] * kw * zr + delta[a];
                phi[a][idx]     += A*env*cos(ph);
                phi_vel[a][idx] += chiral[1]*omega*A*env*sin(ph);
            }
        }

        /* Braid 3: along y-axis (rotate z→y via x-axis 90° rotation) */
        {
            /* In rotated frame: zr=y, yr=-z, xr=x */
            double zr = y, yr = -z, xr = x;
            double r2e = xr*xr/(sx*sx) + yr*yr/(sy*sy);
            double env = exp(-r2e * inv2R2);
            double omega = sqrt(kw*kw + m2);
            for (int a=0; a<NFIELDS; a++) {
                double ph = chiral[2] * kw * zr + delta[a];
                phi[a][idx]     += A*env*cos(ph);
                phi_vel[a][idx] += chiral[2]*omega*A*env*sin(ph);
            }
        }
    }}}

    /* Write SFA */
    uint8_t sfa_dtype = (precision==0)?SFA_F16:(precision==1)?SFA_F32:SFA_F64;
    SFA *sfa = sfa_create(outpath, N, N, N, L, L, L, dt);

    char vN[32],vL[32],vA[32],vR[32],vel[32],vch[32];
    snprintf(vN,32,"%d",N); snprintf(vL,32,"%.6f",L); snprintf(vA,32,"%.6f",A);
    snprintf(vR,32,"%.6f",R_tube); snprintf(vel,32,"%.6f",ellip); snprintf(vch,32,"%s",chirality);
    const char *keys[]={"N","L","A","R_tube","ellip","chirality","type"};
    const char *vals[]={vN,vL,vA,vR,vel,vch,"3braid_xyz"};
    sfa_add_kvmd(sfa, 0, 0xFFFFFFFF, 0xFFFFFFFF, keys, vals, 7);

    sfa_add_column(sfa,"phi_x",sfa_dtype,SFA_POSITION,0);
    sfa_add_column(sfa,"phi_y",sfa_dtype,SFA_POSITION,1);
    sfa_add_column(sfa,"phi_z",sfa_dtype,SFA_POSITION,2);
    sfa_add_column(sfa,"theta_x",sfa_dtype,SFA_ANGLE,0);
    sfa_add_column(sfa,"theta_y",sfa_dtype,SFA_ANGLE,1);
    sfa_add_column(sfa,"theta_z",sfa_dtype,SFA_ANGLE,2);
    sfa_add_column(sfa,"phi_vx",sfa_dtype,SFA_VELOCITY,0);
    sfa_add_column(sfa,"phi_vy",sfa_dtype,SFA_VELOCITY,1);
    sfa_add_column(sfa,"phi_vz",sfa_dtype,SFA_VELOCITY,2);
    sfa_add_column(sfa,"theta_vx",sfa_dtype,SFA_VELOCITY,3);
    sfa_add_column(sfa,"theta_vy",sfa_dtype,SFA_VELOCITY,4);
    sfa_add_column(sfa,"theta_vz",sfa_dtype,SFA_VELOCITY,5);
    sfa_finalize_header(sfa);

    if (precision == 2) {
        void *cols[12] = {phi[0],phi[1],phi[2],theta[0],theta[1],theta[2],
                          phi_vel[0],phi_vel[1],phi_vel[2],theta_vel[0],theta_vel[1],theta_vel[2]};
        sfa_write_frame(sfa, 0.0, cols);
    } else {
        void *cols[12];
        int es = (precision==0)?2:4;
        for (int c=0;c<12;c++) {
            double *src=(c<3)?phi[c]:(c<6)?theta[c-3]:(c<9)?phi_vel[c-6]:theta_vel[c-9];
            cols[c]=malloc(N3*es);
            if(precision==1){float *p=(float*)cols[c];for(long i=0;i<N3;i++) p[i]=(float)src[i];}
            else{uint16_t *p=(uint16_t*)cols[c];for(long i=0;i<N3;i++) p[i]=f64_to_f16(src[i]);}
        }
        sfa_write_frame(sfa, 0.0, cols);
        for(int c=0;c<12;c++) free(cols[c]);
    }

    sfa_close(sfa);
    fprintf(stderr, "gen_3braid: wrote %s (%d^3, chirality=%s)\n", outpath, N, chirality);
    for(int a=0;a<3;a++){free(phi[a]);free(phi_vel[a]);free(theta[a]);free(theta_vel[a]);}
    return 0;
}
