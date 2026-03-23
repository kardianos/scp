/*  gen_oscillon.c — Generate a Gaussian oscillon seed as a single-frame SFA
 *
 *  Build: gcc -O3 -o gen_oscillon gen_oscillon.c -lzstd -lm
 *  Usage: ./gen_oscillon -o seed.sfa [-N 128] [-L 10] [-A 0.8] [-sigma 3.0]
 */

#define SFA_IMPLEMENTATION
#include "../format/sfa.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define NFIELDS 3
#define PI 3.14159265358979323846

int main(int argc, char **argv) {
    int N = 128;
    double L = 10.0, A = 0.8, sigma = 3.0;
    double delta[3] = {0.0, 3.0005, 4.4325};
    double m2 = 2.25, eta = 0.5, mu = -41.345, kappa = 50.0;
    double dt_factor = 0.025;
    double cx = 0, cy = 0, cz = 0;
    char outpath[512] = "oscillon_seed.sfa";
    int precision = 1;

    for (int i = 1; i < argc - 1; i += 2) {
        const char *k = argv[i], *v = argv[i+1];
        if      (!strcmp(k,"-N"))      N = atoi(v);
        else if (!strcmp(k,"-L"))      L = atof(v);
        else if (!strcmp(k,"-A"))      A = atof(v);
        else if (!strcmp(k,"-sigma"))  sigma = atof(v);
        else if (!strcmp(k,"-delta"))  sscanf(v, "%lf,%lf,%lf", &delta[0], &delta[1], &delta[2]);
        else if (!strcmp(k,"-cx"))     cx = atof(v);
        else if (!strcmp(k,"-cy"))     cy = atof(v);
        else if (!strcmp(k,"-cz"))     cz = atof(v);
        else if (!strcmp(k,"-o"))      strncpy(outpath, v, 511);
        else if (!strcmp(k,"-precision")) {
            if (!strcmp(v,"f16")) precision = 0;
            else if (!strcmp(v,"f32")) precision = 1;
            else if (!strcmp(v,"f64")) precision = 2;
        }
    }

    long N3 = (long)N * N * N;
    double dx = 2.0 * L / (N - 1);
    double dt = dt_factor * dx;
    int NN = N * N;

    fprintf(stderr, "gen_oscillon: N=%d L=%.1f A=%.2f sigma=%.1f center=(%.1f,%.1f,%.1f)\n",
            N, L, A, sigma, cx, cy, cz);

    double *phi[3], *theta[3], *phi_vel[3], *theta_vel[3];
    for (int a = 0; a < 3; a++) {
        phi[a]       = (double*)calloc(N3, sizeof(double));
        phi_vel[a]   = (double*)calloc(N3, sizeof(double));
        theta[a]     = (double*)calloc(N3, sizeof(double));
        theta_vel[a] = (double*)calloc(N3, sizeof(double));
    }

    for (int i = 0; i < N; i++) { double x = -L + i*dx;
    for (int j = 0; j < N; j++) { double y = -L + j*dx;
    for (int k = 0; k < N; k++) { double z = -L + k*dx;
        long idx = (long)i*NN + j*N + k;
        double r2 = (x-cx)*(x-cx) + (y-cy)*(y-cy) + (z-cz)*(z-cz);
        double env = A * exp(-r2 / (2.0 * sigma * sigma));
        for (int a = 0; a < NFIELDS; a++)
            phi[a][idx] = env * cos(delta[a]);
    }}}

    uint8_t sfa_dtype = (precision == 0) ? SFA_F16 : (precision == 1) ? SFA_F32 : SFA_F64;
    SFA *sfa = sfa_create(outpath, N, N, N, L, L, L, dt);

    char vN[32],vL[32],vA[32],vs[32],vdelta[64],vcx[32],vcy[32],vcz[32];
    snprintf(vN,32,"%d",N); snprintf(vL,32,"%.6f",L); snprintf(vA,32,"%.6f",A);
    snprintf(vs,32,"%.6f",sigma);
    snprintf(vdelta,64,"%.6f,%.6f,%.6f",delta[0],delta[1],delta[2]);
    snprintf(vcx,32,"%.6f",cx); snprintf(vcy,32,"%.6f",cy); snprintf(vcz,32,"%.6f",cz);
    const char *keys[] = {"N","L","A","sigma","delta","cx","cy","cz"};
    const char *vals[] = {vN,vL,vA,vs,vdelta,vcx,vcy,vcz};
    sfa_add_kvmd(sfa, 0, 0xFFFFFFFF, 0xFFFFFFFF, keys, vals, 8);

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
        void *cols[12] = { phi[0],phi[1],phi[2], theta[0],theta[1],theta[2],
                           phi_vel[0],phi_vel[1],phi_vel[2], theta_vel[0],theta_vel[1],theta_vel[2] };
        sfa_write_frame(sfa, 0.0, cols);
    } else {
        void *cols[12];
        int es = (precision == 0) ? 2 : 4;
        for (int c = 0; c < 12; c++) {
            double *src = (c<3)?phi[c]:(c<6)?theta[c-3]:(c<9)?phi_vel[c-6]:theta_vel[c-9];
            cols[c] = malloc(N3 * es);
            if (precision == 1) { float *p=(float*)cols[c]; for(long i=0;i<N3;i++) p[i]=(float)src[i]; }
            else { uint16_t *p=(uint16_t*)cols[c]; for(long i=0;i<N3;i++) {
                float f=(float)src[i]; uint32_t x; memcpy(&x,&f,4);
                uint16_t sign=(x>>16)&0x8000; int exp=((x>>23)&0xFF)-127+15; uint16_t mant=(x>>13)&0x3FF;
                p[i] = (exp<=0)?sign:(exp>=31)?(sign|0x7C00):(sign|(exp<<10)|mant);
            }}
        }
        sfa_write_frame(sfa, 0.0, cols);
        for (int c = 0; c < 12; c++) free(cols[c]);
    }

    sfa_close(sfa);
    fprintf(stderr, "gen_oscillon: wrote %s (%d^3, %s)\n", outpath, N,
            (const char*[]){"f16","f32","f64"}[precision]);

    for (int a = 0; a < 3; a++) { free(phi[a]); free(phi_vel[a]); free(theta[a]); free(theta_vel[a]); }
    return 0;
}
