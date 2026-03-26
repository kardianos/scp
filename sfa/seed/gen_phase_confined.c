/*  gen_phase_confined.c — Generate phase-confined 3-braid seed
 *
 *  Three braids along x/y/z with carrier phase offsets 0, 2π/3, 4π/3.
 *  The phase offset prevents merging (P → 0 at triple overlap).
 *  Pre-loads θ and contracting velocity profile per V41 stability signatures.
 *
 *  Build: gcc -O3 -fopenmp -o gen_phase_confined gen_phase_confined.c \
 *         -I/home/d/code/scp -lzstd -lm
 *
 *  Usage: ./gen_phase_confined -chirality UDD -o neutron_seed.sfa
 *         ./gen_phase_confined -chirality UUD -o proton_seed.sfa
 */

#define SFA_IMPLEMENTATION
#include "../format/sfa.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <omp.h>

#define NFIELDS 3
#define PI 3.14159265358979323846

static inline uint16_t f64_to_f16(double v) {
    float f=(float)v; uint32_t x; memcpy(&x,&f,4);
    uint16_t sign=(x>>16)&0x8000; int exp=((x>>23)&0xFF)-127+15; uint16_t mant=(x>>13)&0x3FF;
    if(exp<=0)return sign; if(exp>=31)return sign|0x7C00; return sign|(exp<<10)|mant;
}

int main(int argc, char **argv) {
    int N = 192;
    double L = 25.0, A = 0.3, A_bg = 0.1;
    double R_tube = 4.5, ellip = 0.3325;
    double delta[3] = {0.0, 3.0005, 4.4325};
    double m2 = 2.25, eta = 0.5;
    double dt_factor = 0.025;
    double theta_init_frac = 0.5;
    char chirality[8] = "UDD";
    char outpath[512] = "phase_confined_seed.sfa";
    int precision = 1;

    for (int i=1; i<argc-1; i+=2) {
        if      (!strcmp(argv[i],"-N"))          N = atoi(argv[i+1]);
        else if (!strcmp(argv[i],"-L"))          L = atof(argv[i+1]);
        else if (!strcmp(argv[i],"-A"))          A = atof(argv[i+1]);
        else if (!strcmp(argv[i],"-R"))          R_tube = atof(argv[i+1]);
        else if (!strcmp(argv[i],"-chirality"))  strncpy(chirality, argv[i+1], 7);
        else if (!strcmp(argv[i],"-theta_init")) theta_init_frac = atof(argv[i+1]);
        else if (!strcmp(argv[i],"-o"))          strncpy(outpath, argv[i+1], 511);
        else if (!strcmp(argv[i],"-precision"))  {
            if(!strcmp(argv[i+1],"f16")) precision=0;
            else if(!strcmp(argv[i+1],"f32")) precision=1;
            else if(!strcmp(argv[i+1],"f64")) precision=2;
        }
    }

    long N3 = (long)N*N*N;
    int NN = N*N;
    double dx = 2.0*L/(N-1);
    double dt = dt_factor * dx;
    double kw = PI/L;
    double omega = sqrt(kw*kw + m2);
    double k_bg = PI/L, omega_bg = sqrt(k_bg*k_bg + m2);
    double sx = 1+ellip, sy = 1-ellip;
    double inv2R2 = 1.0/(2*R_tube*R_tube);

    /* Chirality per braid */
    double chiral[3] = {1, 1, 1};
    for (int b=0; b<3 && chirality[b]; b++)
        chiral[b] = (chirality[b]=='D'||chirality[b]=='d') ? -1.0 : 1.0;

    /* CARRIER PHASE OFFSETS — the color charge */
    double carrier_phase[3] = {0.0, 2.0*PI/3.0, 4.0*PI/3.0};

    fprintf(stderr, "gen_phase_confined: chirality=%s phases={0, 2π/3, 4π/3}\n", chirality);
    fprintf(stderr, "  N=%d L=%.0f A=%.2f R=%.1f theta_init=%.1f\n", N, L, A, R_tube, theta_init_frac);

    /* Verify P cancellation at origin */
    double test_P = 0;
    for (int b=0; b<3; b++) test_P += cos(carrier_phase[b]);
    fprintf(stderr, "  Phase sum test: cos(0)+cos(2π/3)+cos(4π/3) = %.6f (should be 0)\n", test_P);

    /* Allocate */
    double *phi[3], *phi_vel[3], *theta[3], *theta_vel[3];
    for (int a=0; a<3; a++) {
        phi[a]       = (double*)calloc(N3, sizeof(double));
        phi_vel[a]   = (double*)calloc(N3, sizeof(double));
        theta[a]     = (double*)calloc(N3, sizeof(double));
        theta_vel[a] = (double*)calloc(N3, sizeof(double));
    }

    #pragma omp parallel for schedule(static)
    for (long idx = 0; idx < N3; idx++) {
        int i=(int)(idx/NN), j=(int)((idx/N)%N), k=(int)(idx%N);
        double x = -L + i*dx, y = -L + j*dx, z = -L + k*dx;
        double r_center = sqrt(x*x + y*y + z*z);

        /* Background */
        for (int a=0; a<NFIELDS; a++) {
            double ph_bg = k_bg*z + 2*PI*a/3.0;
            phi[a][idx] = A_bg*cos(ph_bg);
            phi_vel[a][idx] = omega_bg*A_bg*sin(ph_bg);
        }

        /* Three braids with carrier phase offsets */

        /* Braid 0: along z-axis, carrier phase = 0 */
        {
            double r2e = x*x/(sx*sx) + y*y/(sy*sy);
            double env = A * exp(-r2e * inv2R2);
            for (int a=0; a<NFIELDS; a++) {
                double ph = chiral[0] * kw * z + delta[a] + carrier_phase[0];
                phi[a][idx] += env*cos(ph);
                phi_vel[a][idx] += chiral[0]*omega*env*sin(ph);
            }
        }

        /* Braid 1: along x-axis, carrier phase = 2π/3 */
        {
            double zr = x, xr = -z, yr = y;
            double r2e = xr*xr/(sx*sx) + yr*yr/(sy*sy);
            double env = A * exp(-r2e * inv2R2);
            for (int a=0; a<NFIELDS; a++) {
                double ph = chiral[1] * kw * zr + delta[a] + carrier_phase[1];
                phi[a][idx] += env*cos(ph);
                phi_vel[a][idx] += chiral[1]*omega*env*sin(ph);
            }
        }

        /* Braid 2: along y-axis, carrier phase = 4π/3 */
        {
            double zr = y, yr = -z, xr = x;
            double r2e = xr*xr/(sx*sx) + yr*yr/(sy*sy);
            double env = A * exp(-r2e * inv2R2);
            for (int a=0; a<NFIELDS; a++) {
                double ph = chiral[2] * kw * zr + delta[a] + carrier_phase[2];
                phi[a][idx] += env*cos(ph);
                phi_vel[a][idx] += chiral[2]*omega*env*sin(ph);
            }
        }

        /* Pre-initialize θ with per-braid phase-shifted pattern */
        if (theta_init_frac > 0) {
            for (int b = 0; b < 3; b++) {
                /* Each braid's theta contribution, phase-shifted by carrier_phase[b] */
                double braid_env;
                if (b == 0) {
                    double r2e = x*x/(sx*sx) + y*y/(sy*sy);
                    braid_env = A * exp(-r2e * inv2R2);
                } else if (b == 1) {
                    double xr = -z, yr = y;
                    double r2e = xr*xr/(sx*sx) + yr*yr/(sy*sy);
                    braid_env = A * exp(-r2e * inv2R2);
                } else {
                    double xr = x, yr = -z;
                    double r2e = xr*xr/(sx*sx) + yr*yr/(sy*sy);
                    braid_env = A * exp(-r2e * inv2R2);
                }
                double t_env = theta_init_frac * eta * braid_env * 0.1;
                /* Theta with braid's carrier phase — ensures cancellation for UDD */
                theta[0][idx] += chiral[b] * t_env * cos(kw*z + carrier_phase[b]);
                theta[1][idx] += chiral[b] * t_env * sin(kw*z + carrier_phase[b]);
                theta[2][idx] += chiral[b] * t_env * cos(kw*z + carrier_phase[b] + PI/4);
            }
        }

        /* Contracting velocity profile */
        if (r_center > 0.1) {
            double v_radial = -0.05 * exp(-r_center*r_center / (4*R_tube*R_tube));
            double rhat[3] = {x/r_center, y/r_center, z/r_center};
            for (int a=0; a<NFIELDS; a++)
                phi_vel[a][idx] += v_radial * rhat[a%3] * phi[a][idx];
        }
    }

    /* Measure net theta (should be ~0 for UDD, nonzero for UUD) */
    double net_theta = 0;
    for (long i = 0; i < N3; i++)
        net_theta += theta[0][i]*theta[0][i] + theta[1][i]*theta[1][i] + theta[2][i]*theta[2][i];
    net_theta = sqrt(net_theta / (3.0*N3));
    fprintf(stderr, "  Net theta_rms = %.6f\n", net_theta);

    /* Write SFA */
    uint8_t sfa_dtype = (precision==0)?SFA_F16:(precision==1)?SFA_F32:SFA_F64;
    SFA *sfa = sfa_create(outpath, N, N, N, L, L, L, dt);

    char vch[32], vph[64];
    snprintf(vch, 32, "%s", chirality);
    snprintf(vph, 64, "0,%.6f,%.6f", carrier_phase[1], carrier_phase[2]);
    const char *keys[] = {"chirality", "carrier_phases", "type", "A", "R_tube", "theta_init"};
    char vA[32], vR[32], vti[32];
    snprintf(vA,32,"%.6f",A); snprintf(vR,32,"%.6f",R_tube); snprintf(vti,32,"%.6f",theta_init_frac);
    const char *vals[] = {vch, vph, "phase_confined_3braid", vA, vR, vti};
    sfa_add_kvmd(sfa, 0, 0xFFFFFFFF, 0xFFFFFFFF, keys, vals, 6);

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
        void *cols[12]; int es = (precision==0)?2:4;
        for (int c=0;c<12;c++) {
            double *src=(c<3)?phi[c]:(c<6)?theta[c-3]:(c<9)?phi_vel[c-6]:theta_vel[c-9];
            cols[c]=malloc(N3*es);
            if(precision==1){float *p=(float*)cols[c];for(long i=0;i<N3;i++)p[i]=(float)src[i];}
            else{uint16_t *p=(uint16_t*)cols[c];for(long i=0;i<N3;i++)p[i]=f64_to_f16(src[i]);}
        }
        sfa_write_frame(sfa, 0.0, cols);
        for(int c=0;c<12;c++) free(cols[c]);
    }
    sfa_close(sfa);
    fprintf(stderr, "  Wrote %s\n", outpath);

    for(int a=0;a<3;a++){free(phi[a]);free(phi_vel[a]);free(theta[a]);free(theta_vel[a]);}
    return 0;
}
