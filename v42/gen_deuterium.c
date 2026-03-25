/*  gen_deuterium.c — Generate deuterium seed: UUD proton + UDD neutron
 *
 *  Two phase-confined baryons separated along x-axis.
 *  Each baryon has 3 braids with carrier phase offsets {0, 2π/3, 4π/3}.
 *
 *  Build: gcc -O3 -fopenmp -o gen_deuterium gen_deuterium.c \
 *         -I/home/d/code/scp -lzstd -lm
 *
 *  Usage: ./gen_deuterium -N 512 -L 100 -sep 40 -o deuterium_seed.sfa
 */

#define SFA_IMPLEMENTATION
#include "sfa/format/sfa.h"

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

/* Add one baryon (3 phase-confined braids) centered at (cx,cy,cz) */
static void add_baryon(double *phi[3], double *phi_vel[3], double *theta[3],
    int N, double L, double dx, double cx, double cy, double cz,
    double A, double R_tube, double ellip, double delta[3],
    double chiral[3], double carrier_phase[3],
    double theta_init_frac, double m2, double eta) {

    int NN = N*N;
    long N3 = (long)N*N*N;
    double kw = PI/L;
    double omega = sqrt(kw*kw + m2);
    double sx = 1+ellip, sy = 1-ellip;
    double inv2R2 = 1.0/(2*R_tube*R_tube);

    #pragma omp parallel for schedule(static)
    for (long idx = 0; idx < N3; idx++) {
        int i=(int)(idx/NN), j=(int)((idx/N)%N), k=(int)(idx%N);
        double x = -L + i*dx - cx;
        double y = -L + j*dx - cy;
        double z = -L + k*dx - cz;
        double r_center = sqrt(x*x + y*y + z*z);

        /* Only contribute within ~3*R_tube of baryon center */
        if (r_center > 5.0 * R_tube) continue;

        /* Braid 0: along z-axis */
        {
            double r2e = x*x/(sx*sx) + y*y/(sy*sy);
            double env = A * exp(-r2e * inv2R2);
            for (int a=0; a<NFIELDS; a++) {
                double ph = chiral[0] * kw * z + delta[a] + carrier_phase[0];
                phi[a][idx] += env*cos(ph);
                phi_vel[a][idx] += chiral[0]*omega*env*sin(ph);
            }
        }

        /* Braid 1: along x-axis */
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

        /* Braid 2: along y-axis */
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

        /* Pre-load θ per braid */
        if (theta_init_frac > 0) {
            for (int b = 0; b < 3; b++) {
                double braid_env;
                if (b == 0) { double r2e = x*x/(sx*sx)+y*y/(sy*sy); braid_env = A*exp(-r2e*inv2R2); }
                else if (b == 1) { double xr=-z,yr=y; double r2e=xr*xr/(sx*sx)+yr*yr/(sy*sy); braid_env=A*exp(-r2e*inv2R2); }
                else { double xr=x,yr=-z; double r2e=xr*xr/(sx*sx)+yr*yr/(sy*sy); braid_env=A*exp(-r2e*inv2R2); }
                double t_env = theta_init_frac * eta * braid_env * 0.1;
                theta[0][idx] += chiral[b] * t_env * cos(kw*z + carrier_phase[b]);
                theta[1][idx] += chiral[b] * t_env * sin(kw*z + carrier_phase[b]);
                theta[2][idx] += chiral[b] * t_env * cos(kw*z + carrier_phase[b] + PI/4);
            }
        }

        /* Contracting velocity */
        if (r_center > 0.1) {
            double v_radial = -0.05 * exp(-r_center*r_center / (4*R_tube*R_tube));
            double rhat[3] = {x/r_center, y/r_center, z/r_center};
            for (int a=0; a<NFIELDS; a++)
                phi_vel[a][idx] += v_radial * rhat[a%3] * phi[a][idx];
        }
    }
}

int main(int argc, char **argv) {
    int N = 512;
    double L = 100.0;
    double A = 0.3, A_bg = 0.1;
    double R_tube = 4.5, ellip = 0.3325;
    double separation = 40.0;
    double theta_init = 0.8;
    double m2 = 2.25, eta = 0.5;
    char outpath[512] = "deuterium_seed.sfa";
    int precision = 1;

    for (int i=1; i<argc-1; i+=2) {
        if      (!strcmp(argv[i],"-N"))     N = atoi(argv[i+1]);
        else if (!strcmp(argv[i],"-L"))     L = atof(argv[i+1]);
        else if (!strcmp(argv[i],"-A"))     A = atof(argv[i+1]);
        else if (!strcmp(argv[i],"-R"))     R_tube = atof(argv[i+1]);
        else if (!strcmp(argv[i],"-sep"))   separation = atof(argv[i+1]);
        else if (!strcmp(argv[i],"-theta")) theta_init = atof(argv[i+1]);
        else if (!strcmp(argv[i],"-o"))     strncpy(outpath, argv[i+1], 511);
        else if (!strcmp(argv[i],"-precision")) {
            if(!strcmp(argv[i+1],"f16")) precision=0;
            else if(!strcmp(argv[i+1],"f32")) precision=1;
            else if(!strcmp(argv[i+1],"f64")) precision=2;
        }
    }

    long N3 = (long)N*N*N;
    int NN = N*N;
    double dx = 2.0*L/(N-1);
    double dt = 0.025 * dx;

    fprintf(stderr, "gen_deuterium: N=%d L=%.0f sep=%.0f A=%.2f R=%.1f θ_init=%.1f\n",
            N, L, separation, A, R_tube, theta_init);
    fprintf(stderr, "  Memory: %.1f GB for field arrays\n", 12.0*N3*8/1e9);

    double *phi[3], *phi_vel[3], *theta[3], *theta_vel[3];
    for (int a=0; a<3; a++) {
        phi[a]       = (double*)calloc(N3, sizeof(double));
        phi_vel[a]   = (double*)calloc(N3, sizeof(double));
        theta[a]     = (double*)calloc(N3, sizeof(double));
        theta_vel[a] = (double*)calloc(N3, sizeof(double));
    }

    /* Background */
    double k_bg = PI/L, omega_bg = sqrt(k_bg*k_bg + m2);
    fprintf(stderr, "  Initializing background...\n");
    #pragma omp parallel for schedule(static)
    for (long idx = 0; idx < N3; idx++) {
        int i=(int)(idx/NN), j=(int)((idx/N)%N), k=(int)(idx%N);
        double z = -L + k*dx;
        for (int a=0; a<NFIELDS; a++) {
            double ph_bg = k_bg*z + 2*PI*a/3.0;
            phi[a][idx] = A_bg*cos(ph_bg);
            phi_vel[a][idx] = omega_bg*A_bg*sin(ph_bg);
        }
    }

    double delta[3] = {0.0, 3.0005, 4.4325};
    double carrier[3] = {0.0, 2.0*PI/3.0, 4.0*PI/3.0};

    /* Proton (UUD) at x = -sep/2 */
    double chiral_UUD[3] = {1.0, 1.0, -1.0};
    fprintf(stderr, "  Adding UUD proton at x=%.0f...\n", -separation/2);
    add_baryon(phi, phi_vel, theta, N, L, dx,
               -separation/2, 0, 0,
               A, R_tube, ellip, delta, chiral_UUD, carrier,
               theta_init, m2, eta);

    /* Neutron (UDD) at x = +sep/2 */
    double chiral_UDD[3] = {1.0, -1.0, -1.0};
    fprintf(stderr, "  Adding UDD neutron at x=+%.0f...\n", separation/2);
    add_baryon(phi, phi_vel, theta, N, L, dx,
               separation/2, 0, 0,
               A, R_tube, ellip, delta, chiral_UDD, carrier,
               theta_init, m2, eta);

    /* Write SFA */
    uint8_t sfa_dtype = (precision==0)?SFA_F16:(precision==1)?SFA_F32:SFA_F64;
    SFA *sfa = sfa_create(outpath, N, N, N, L, L, L, dt);

    char vsep[32]; snprintf(vsep, 32, "%.1f", separation);
    const char *keys[] = {"type", "separation", "proton_pos", "neutron_pos"};
    const char *vals[] = {"deuterium_2H", vsep, "-20,0,0", "+20,0,0"};
    sfa_add_kvmd(sfa, 0, 0xFFFFFFFF, 0xFFFFFFFF, keys, vals, 4);

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

    fprintf(stderr, "  Writing SFA...\n");
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
