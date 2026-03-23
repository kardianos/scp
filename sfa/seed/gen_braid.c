/*  gen_braid.c — Generate a braid seed as a single-frame 12-column SFA
 *
 *  Writes a single SFA frame containing the standard V34 braid initialization
 *  (3 phi + 3 theta fields + 6 velocities). Theta and theta velocities are zero.
 *  Physics parameters are embedded as KVMD metadata.
 *
 *  Build: gcc -O3 -o gen_braid gen_braid.c -lzstd -lm
 *
 *  Usage:
 *    ./gen_braid -o seed.sfa                          # default braid
 *    ./gen_braid -N 256 -L 20 -A 0.9 -o big_seed.sfa # custom
 *    ./gen_braid -o -                                  # stdout (for scp_sim exec init)
 */

#define SFA_IMPLEMENTATION
#include "../format/sfa.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <unistd.h>

#define NFIELDS 3
#define PI 3.14159265358979323846

int main(int argc, char **argv) {
    /* Defaults */
    int N = 128;
    double L = 15.0;
    double A = 0.8, A_bg = 0.1;
    double R_tube = 3.0, ellip = 0.3325;
    double delta[3] = {0.0, 3.0005, 4.4325};
    double m2 = 2.25, mtheta2 = 0.0, eta = 0.5;
    double mu = -41.345, kappa = 50.0;
    double dt_factor = 0.025;
    double cx = 0, cy = 0, cz = 0;  /* braid center in world coords */
    char outpath[512] = "seed.sfa";
    int precision = 1;  /* f32 */

    for (int i = 1; i < argc - 1; i += 2) {
        const char *k = argv[i], *v = argv[i+1];
        if      (!strcmp(k,"-N"))      N = atoi(v);
        else if (!strcmp(k,"-L"))      L = atof(v);
        else if (!strcmp(k,"-A"))      A = atof(v);
        else if (!strcmp(k,"-A_bg"))   A_bg = atof(v);
        else if (!strcmp(k,"-R"))      R_tube = atof(v);
        else if (!strcmp(k,"-ellip"))  ellip = atof(v);
        else if (!strcmp(k,"-delta"))  sscanf(v, "%lf,%lf,%lf", &delta[0], &delta[1], &delta[2]);
        else if (!strcmp(k,"-m"))    { double m = atof(v); m2 = m*m; }
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

    fprintf(stderr, "gen_braid: N=%d L=%.1f A=%.2f R=%.1f ellip=%.4f center=(%.1f,%.1f,%.1f)\n",
            N, L, A, R_tube, ellip, cx, cy, cz);

    /* Allocate field arrays */
    double *phi[3], *phi_vel[3], *theta[3], *theta_vel[3];
    for (int a = 0; a < 3; a++) {
        phi[a]       = (double*)calloc(N3, sizeof(double));
        phi_vel[a]   = (double*)calloc(N3, sizeof(double));
        theta[a]     = (double*)calloc(N3, sizeof(double));
        theta_vel[a] = (double*)calloc(N3, sizeof(double));
    }

    /* Generate braid */
    double kw = PI / L, omega = sqrt(kw*kw + m2);
    double sx = 1 + ellip, sy = 1 - ellip;
    double inv2R2 = 1.0 / (2 * R_tube * R_tube);
    double k_bg = PI / L, omega_bg = sqrt(k_bg*k_bg + m2);

    for (int i = 0; i < N; i++) { double x = -L + i*dx;
    for (int j = 0; j < N; j++) { double y = -L + j*dx;
    for (int k = 0; k < N; k++) { double z = -L + k*dx;
        long idx = (long)i*NN + j*N + k;
        double xc = x - cx, yc = y - cy;
        double r2e = xc*xc/(sx*sx) + yc*yc/(sy*sy);
        double env = exp(-r2e * inv2R2);
        for (int a = 0; a < NFIELDS; a++) {
            double ph = kw*z + delta[a];
            double ph_bg = k_bg*z + 2*PI*a/3.0;
            phi[a][idx] = A*env*cos(ph) + A_bg*cos(ph_bg);
            phi_vel[a][idx] = omega*A*env*sin(ph) + omega_bg*A_bg*sin(ph_bg);
        }
    }}}

    /* Write SFA */
    uint8_t sfa_dtype = (precision == 0) ? SFA_F16 : (precision == 1) ? SFA_F32 : SFA_F64;
    SFA *sfa = sfa_create(outpath, N, N, N, L, L, L, dt);

    /* Embed parameters as KVMD */
    char vN[32],vL[32],vm[32],veta[32],vmu[32],vkappa[32],vdelta[64];
    char vA[32],vAbg[32],vR[32],vellip[32],vcx[32],vcy[32],vcz[32];
    snprintf(vN,32,"%d",N); snprintf(vL,32,"%.6f",L);
    snprintf(vm,32,"%.6f",sqrt(m2)); snprintf(veta,32,"%.6f",eta);
    snprintf(vmu,32,"%.6f",mu); snprintf(vkappa,32,"%.6f",kappa);
    snprintf(vdelta,64,"%.6f,%.6f,%.6f",delta[0],delta[1],delta[2]);
    snprintf(vA,32,"%.6f",A); snprintf(vAbg,32,"%.6f",A_bg);
    snprintf(vR,32,"%.6f",R_tube); snprintf(vellip,32,"%.6f",ellip);
    snprintf(vcx,32,"%.6f",cx); snprintf(vcy,32,"%.6f",cy); snprintf(vcz,32,"%.6f",cz);
    const char *keys[] = {"N","L","m","eta","mu","kappa","delta",
                          "A","A_bg","R_tube","ellip","cx","cy","cz"};
    const char *vals[] = {vN,vL,vm,veta,vmu,vkappa,vdelta,
                          vA,vAbg,vR,vellip,vcx,vcy,vcz};
    sfa_add_kvmd(sfa, 0, 0xFFFFFFFF, 0xFFFFFFFF, keys, vals, 14);

    sfa_add_column(sfa, "phi_x",    sfa_dtype, SFA_POSITION, 0);
    sfa_add_column(sfa, "phi_y",    sfa_dtype, SFA_POSITION, 1);
    sfa_add_column(sfa, "phi_z",    sfa_dtype, SFA_POSITION, 2);
    sfa_add_column(sfa, "theta_x",  sfa_dtype, SFA_ANGLE,    0);
    sfa_add_column(sfa, "theta_y",  sfa_dtype, SFA_ANGLE,    1);
    sfa_add_column(sfa, "theta_z",  sfa_dtype, SFA_ANGLE,    2);
    sfa_add_column(sfa, "phi_vx",   sfa_dtype, SFA_VELOCITY, 0);
    sfa_add_column(sfa, "phi_vy",   sfa_dtype, SFA_VELOCITY, 1);
    sfa_add_column(sfa, "phi_vz",   sfa_dtype, SFA_VELOCITY, 2);
    sfa_add_column(sfa, "theta_vx", sfa_dtype, SFA_VELOCITY, 3);
    sfa_add_column(sfa, "theta_vy", sfa_dtype, SFA_VELOCITY, 4);
    sfa_add_column(sfa, "theta_vz", sfa_dtype, SFA_VELOCITY, 5);
    sfa_finalize_header(sfa);

    /* Write single frame */
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
    fprintf(stderr, "gen_braid: wrote %s (%d^3, %s)\n", outpath, N,
            (const char*[]){"f16","f32","f64"}[precision]);

    for (int a = 0; a < 3; a++) { free(phi[a]); free(phi_vel[a]); free(theta[a]); free(theta_vel[a]); }
    return 0;
}
