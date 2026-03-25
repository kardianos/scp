/*  construct_seed.c — First-principles seed construction with stability prediction
 *
 *  Builds a 3-braid xyz seed following the three stability signatures:
 *    1. θ confinement (pre-loaded, not zero)
 *    2. Velocity structure (breathing or contracting profile)
 *    3. |P| concentration (tuned amplitude for P_peak ≈ 0.082)
 *
 *  After construction, computes 7 stability metrics and outputs a prediction score.
 *
 *  Build: gcc -O3 -fopenmp -o construct_seed construct_seed.c -I/home/d/code/scp -lzstd -lm
 *  Usage: ./construct_seed -A 0.5 -R 3.5 -theta_init 0.5 -v_profile contracting \
 *             -chirality UDD -N 192 -L 25 -o seed.sfa
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

/* Physics parameters */
static const double MU = -41.345, KAPPA = 50.0, MASS2 = 2.25;
static const double MTHETA2 = 0.0, ETA = 0.5;

/* Optimal P for maximum binding force: d/dP[V'(P)] = 0 */
/* V'(P) = mu*P/(1+kappa*P^2)^2, maximum |V'| at P_opt = 1/sqrt(3*kappa) */
static const double P_OPT = 0.0816;  /* 1/sqrt(3*50) = 0.0816 */

/* ================================================================
   Velocity profile types
   ================================================================ */

enum { VP_ZERO=0, VP_BREATHING=1, VP_CONTRACTING=2, VP_MIXED=3 };

static const char *vp_names[] = {"zero", "breathing", "contracting", "mixed"};

static int parse_vp(const char *s) {
    if (!strcmp(s,"zero")) return VP_ZERO;
    if (!strcmp(s,"breathing")) return VP_BREATHING;
    if (!strcmp(s,"contracting")) return VP_CONTRACTING;
    if (!strcmp(s,"mixed")) return VP_MIXED;
    return VP_ZERO;
}

/* ================================================================
   f16 helper
   ================================================================ */

static inline uint16_t f64_to_f16(double v) {
    float f=(float)v; uint32_t x; memcpy(&x,&f,4);
    uint16_t sign=(x>>16)&0x8000; int exp=((x>>23)&0xFF)-127+15; uint16_t mant=(x>>13)&0x3FF;
    if(exp<=0)return sign; if(exp>=31)return sign|0x7C00; return sign|(exp<<10)|mant;
}

/* ================================================================
   Construct the seed
   ================================================================ */

int main(int argc, char **argv) {
    int N = 192;
    double L = 25.0, A = 0.5, R_tube = 3.5;
    double theta_init_frac = 0.5;  /* fraction of equilibrium θ to pre-load */
    int v_profile = VP_CONTRACTING;
    char chirality[8] = "UDD";
    char outpath[512] = "seed.sfa";
    double delta[3] = {0.0, 3.0005, 4.4325};
    double ellip = 0.3325;
    double A_bg = 0.1;
    int precision = 1;  /* f32 */
    int quiet = 0;

    for (int i=1; i<argc-1; i+=2) {
        if      (!strcmp(argv[i],"-N"))          N = atoi(argv[i+1]);
        else if (!strcmp(argv[i],"-L"))          L = atof(argv[i+1]);
        else if (!strcmp(argv[i],"-A"))          A = atof(argv[i+1]);
        else if (!strcmp(argv[i],"-R"))          R_tube = atof(argv[i+1]);
        else if (!strcmp(argv[i],"-theta_init")) theta_init_frac = atof(argv[i+1]);
        else if (!strcmp(argv[i],"-v_profile"))  v_profile = parse_vp(argv[i+1]);
        else if (!strcmp(argv[i],"-chirality"))  strncpy(chirality, argv[i+1], 7);
        else if (!strcmp(argv[i],"-ellip"))      ellip = atof(argv[i+1]);
        else if (!strcmp(argv[i],"-A_bg"))       A_bg = atof(argv[i+1]);
        else if (!strcmp(argv[i],"-o"))          strncpy(outpath, argv[i+1], 511);
        else if (!strcmp(argv[i],"-precision"))  {
            if(!strcmp(argv[i+1],"f16")) precision=0;
            else if(!strcmp(argv[i+1],"f32")) precision=1;
            else if(!strcmp(argv[i+1],"f64")) precision=2;
        }
        else if (!strcmp(argv[i],"-quiet"))      { quiet=1; i--; }
    }

    long N3 = (long)N*N*N;
    int NN = N*N;
    double dx = 2.0*L/(N-1);
    double dt = 0.025 * dx;
    double kw = PI/L;
    double omega = sqrt(kw*kw + MASS2);
    double k_bg = PI/L, omega_bg = sqrt(k_bg*k_bg + MASS2);
    double sx = 1+ellip, sy = 1-ellip;
    double inv2R2 = 1.0/(2*R_tube*R_tube);

    /* Chirality: +1 = up, -1 = down */
    double chiral[3] = {1, 1, 1};
    for (int b=0; b<3 && chirality[b]; b++)
        chiral[b] = (chirality[b]=='D'||chirality[b]=='d') ? -1.0 : 1.0;

    /* Allocate */
    double *phi[3], *phi_vel[3], *theta[3], *theta_vel[3];
    for (int a=0; a<3; a++) {
        phi[a]       = (double*)calloc(N3, sizeof(double));
        phi_vel[a]   = (double*)calloc(N3, sizeof(double));
        theta[a]     = (double*)calloc(N3, sizeof(double));
        theta_vel[a] = (double*)calloc(N3, sizeof(double));
    }

    /* ================================================================
       Construct 3-braid field
       ================================================================ */

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

        /* Braid 1: along z-axis */
        {
            double r2e = x*x/(sx*sx) + y*y/(sy*sy);
            double env = A * exp(-r2e * inv2R2);
            for (int a=0; a<NFIELDS; a++) {
                double ph = chiral[0] * kw * z + delta[a];
                phi[a][idx] += env*cos(ph);
                phi_vel[a][idx] += chiral[0]*omega*env*sin(ph);
            }
        }

        /* Braid 2: along x-axis (rotate z→x) */
        {
            double zr = x, xr = -z, yr = y;
            double r2e = xr*xr/(sx*sx) + yr*yr/(sy*sy);
            double env = A * exp(-r2e * inv2R2);
            for (int a=0; a<NFIELDS; a++) {
                double ph = chiral[1] * kw * zr + delta[a];
                phi[a][idx] += env*cos(ph);
                phi_vel[a][idx] += chiral[1]*omega*env*sin(ph);
            }
        }

        /* Braid 3: along y-axis (rotate z→y) */
        {
            double zr = y, yr = -z, xr = x;
            double r2e = xr*xr/(sx*sx) + yr*yr/(sy*sy);
            double env = A * exp(-r2e * inv2R2);
            for (int a=0; a<NFIELDS; a++) {
                double ph = chiral[2] * kw * zr + delta[a];
                phi[a][idx] += env*cos(ph);
                phi_vel[a][idx] += chiral[2]*omega*env*sin(ph);
            }
        }

        /* Pre-initialize θ based on curl(φ) structure */
        if (theta_init_frac > 0) {
            /* θ equilibrium ≈ η*curl(φ)/m_θ² (for massive θ) or time-integrated
               For massless θ, equilibrium is where ∂²θ/∂t² = 0:
               θ_eq such that ∇²θ + η*curl(φ) = 0
               Approximate: θ ~ -η * integral(curl(φ)) ~ η * φ_envelope / k
               Use a fraction of the phi envelope as theta init */
            double env_total = 0;
            for (int a=0; a<3; a++) env_total += phi[a][idx]*phi[a][idx];
            env_total = sqrt(env_total / 3.0);
            double theta_env = theta_init_frac * ETA * env_total * 0.1;  /* scale factor */
            /* Give θ a rotational pattern matching the braid helicity */
            double angle = atan2(y, x);
            theta[0][idx] = theta_env * cos(angle + kw*z);
            theta[1][idx] = theta_env * sin(angle + kw*z);
            theta[2][idx] = theta_env * cos(kw*z + PI/4);
        }

        /* Velocity profile modification */
        if (v_profile != VP_ZERO && r_center > 0.1) {
            double r_hat[3] = {x/r_center, y/r_center, z/r_center};
            double v_radial = 0;

            switch (v_profile) {
            case VP_BREATHING:
                /* |v| increases outward: v_radial > 0 in outer shell */
                v_radial = 0.1 * (r_center / (R_tube * 2.0) - 0.5);
                break;
            case VP_CONTRACTING:
                /* |v| decreases outward: v_radial < 0 everywhere */
                v_radial = -0.05 * exp(-r_center*r_center / (4*R_tube*R_tube));
                break;
            case VP_MIXED:
                /* Inner contracting, outer breathing */
                if (r_center < R_tube * 1.5)
                    v_radial = -0.03 * exp(-r_center*r_center / (R_tube*R_tube));
                else
                    v_radial = 0.05 * (r_center / (R_tube * 3.0) - 0.5);
                break;
            }

            for (int a=0; a<NFIELDS; a++)
                phi_vel[a][idx] += v_radial * r_hat[a%3] * phi[a][idx];
        }
    }

    /* ================================================================
       Compute stability metrics BEFORE writing
       ================================================================ */

    double rho_inner=0, rho_outer=0, P_inner=0, P_outer=0;
    double theta_inner=0, theta_outer=0, v_inner=0, v_outer=0;
    double E_pot_total=0, P_int_total=0, P_int_core=0;
    double force_sum=0, force_max_sum=0;
    long n_inner=0, n_outer=0;

    double R_core = R_tube * 1.5;  /* inner region radius */
    double R_outer = R_tube * 4.0; /* outer region radius */

    #pragma omp parallel for reduction(+:rho_inner,rho_outer,P_inner,P_outer, \
        theta_inner,theta_outer,v_inner,v_outer,E_pot_total,P_int_total,P_int_core, \
        force_sum,force_max_sum,n_inner,n_outer) schedule(static)
    for (long idx = 0; idx < N3; idx++) {
        int i=(int)(idx/NN), j=(int)((idx/N)%N), k=(int)(idx%N);
        double x=-L+i*dx, y=-L+j*dx, z=-L+k*dx;
        double r = sqrt(x*x + y*y + z*z);
        double dV = dx*dx*dx;

        double p0=phi[0][idx], p1=phi[1][idx], p2=phi[2][idx];
        double rho = p0*p0 + p1*p1 + p2*p2;
        double P = p0*p1*p2;
        double P2 = P*P;
        double Pabs = fabs(P);
        double thr = theta[0][idx]*theta[0][idx] + theta[1][idx]*theta[1][idx] + theta[2][idx]*theta[2][idx];
        double vmag = 0;
        for (int a=0; a<3; a++) vmag += phi_vel[a][idx]*phi_vel[a][idx];
        vmag = sqrt(vmag);

        /* E_pot */
        E_pot_total += (MU/2.0)*P2/(1.0+KAPPA*P2)*dV;
        P_int_total += Pabs * dV;
        if (r < R_core) P_int_core += Pabs * dV;

        /* Force computation at sample points */
        if (i>0 && i<N-1 && j>0 && j<N-1 && k>0 && k<N-1 && (idx % 100 == 0)) {
            double idx2 = 1.0/(dx*dx);
            for (int a=0; a<3; a++) {
                long nip=(long)(i+1)*NN+j*N+k, nim=(long)(i-1)*NN+j*N+k;
                long njp=(long)i*NN+(j+1)*N+k, njm=(long)i*NN+(j-1)*N+k;
                long nkp=(long)i*NN+j*N+(k+1), nkm=(long)i*NN+j*N+(k-1);
                double lap = (phi[a][nip]+phi[a][nim]+phi[a][njp]+phi[a][njm]+phi[a][nkp]+phi[a][nkm]-6*phi[a][idx])*idx2;
                double dPda = (a==0)?p1*p2:(a==1)?p0*p2:p0*p1;
                double dVdP = MU*P/((1+KAPPA*P2)*(1+KAPPA*P2));
                double F_lap = fabs(lap), F_mass = fabs(MASS2*phi[a][idx]);
                double F_pot = fabs(dVdP*dPda);
                double F_net = fabs(lap - MASS2*phi[a][idx] - dVdP*dPda);
                double F_max = F_lap; if(F_mass>F_max) F_max=F_mass; if(F_pot>F_max) F_max=F_pot;
                force_sum += F_net;
                force_max_sum += F_max;
            }
        }

        /* Inner/outer accumulation */
        if (r < R_core) {
            rho_inner += rho; P_inner += Pabs; theta_inner += sqrt(thr); v_inner += vmag;
            n_inner++;
        } else if (r > R_outer) {
            rho_outer += rho; P_outer += Pabs; theta_outer += sqrt(thr); v_outer += vmag;
            n_outer++;
        }
    }

    /* Normalize */
    if (n_inner > 0) { rho_inner/=n_inner; P_inner/=n_inner; theta_inner/=n_inner; v_inner/=n_inner; }
    if (n_outer > 0) { rho_outer/=n_outer; P_outer/=n_outer; theta_outer/=n_outer; v_outer/=n_outer; }

    /* Metrics */
    double m_P_peak = P_inner;  /* average |P| in core ≈ P_peak proxy */
    double m_theta_ratio = (theta_inner > 1e-10) ? theta_outer / theta_inner : 0;
    double m_v_ratio = (v_inner > 1e-10) ? v_outer / v_inner : 0;
    double m_rho_ratio = (rho_outer > 1e-10) ? rho_inner / rho_outer : 0;
    double m_P_conc = (P_int_total > 1e-10) ? P_int_core / P_int_total : 0;
    double m_force_bal = (force_max_sum > 1e-10) ? force_sum / force_max_sum : 1;

    /* Scoring functions — map each metric to [0, 1] */
    double s_P = exp(-pow((m_P_peak - P_OPT) / 0.05, 2));  /* Gaussian around P_opt */
    double s_theta = (m_theta_ratio < 0.7) ? 1.0 : (m_theta_ratio < 1.5) ? (1.5 - m_theta_ratio)/0.8 : 0;
    double s_v = 0;
    if (v_profile == VP_CONTRACTING || v_profile == VP_MIXED)
        s_v = (m_v_ratio < 0.7) ? 1.0 : (m_v_ratio < 1.0) ? (1.0 - m_v_ratio)/0.3 : 0;
    else
        s_v = (m_v_ratio > 1.3) ? 1.0 : (m_v_ratio > 0.8) ? (m_v_ratio - 0.8)/0.5 : 0;
    double s_rho = (m_rho_ratio > 10) ? 1.0 : (m_rho_ratio > 2) ? (m_rho_ratio - 2)/8.0 : 0;
    double s_Epot = fmin(1.0, fabs(E_pot_total) / 200.0);
    double s_Pconc = m_P_conc;  /* already 0-1 */
    double s_force = 1.0 - fmin(1.0, m_force_bal);  /* lower = better balanced */

    /* Weights */
    double S_pred = 0.20*s_P + 0.20*s_theta + 0.15*s_v + 0.10*s_rho +
                    0.15*s_Epot + 0.10*s_Pconc + 0.10*s_force;

    /* Output metrics */
    if (!quiet) {
        fprintf(stderr, "construct_seed: A=%.2f R=%.1f theta_init=%.1f v=%s chirality=%s\n",
                A, R_tube, theta_init_frac, vp_names[v_profile], chirality);
        fprintf(stderr, "  P_peak(core)=%.4f (opt=%.4f) theta_out/in=%.3f v_out/in=%.3f\n",
                m_P_peak, P_OPT, m_theta_ratio, m_v_ratio);
        fprintf(stderr, "  rho_in/out=%.1f E_pot=%.1f P_conc=%.3f force_bal=%.3f\n",
                m_rho_ratio, E_pot_total, m_P_conc, m_force_bal);
        fprintf(stderr, "  S_pred=%.4f (components: P=%.2f θ=%.2f v=%.2f ρ=%.2f Ep=%.2f Pc=%.2f F=%.2f)\n",
                S_pred, s_P, s_theta, s_v, s_rho, s_Epot, s_Pconc, s_force);
    }

    /* TSV output to stdout: one line per seed for easy collection */
    printf("%s\t%.2f\t%.1f\t%.1f\t%s\t%.4f\t%.4f\t%.3f\t%.3f\t%.1f\t%.1f\t%.3f\t%.3f\t%.4f\t%s\n",
           chirality, A, R_tube, theta_init_frac, vp_names[v_profile],
           m_P_peak, P_OPT, m_theta_ratio, m_v_ratio, m_rho_ratio,
           E_pot_total, m_P_conc, m_force_bal, S_pred, outpath);

    /* ================================================================
       Write SFA
       ================================================================ */

    uint8_t sfa_dtype = (precision==0)?SFA_F16:(precision==1)?SFA_F32:SFA_F64;
    SFA *sfa = sfa_create(outpath, N, N, N, L, L, L, dt);

    char vA[32],vR[32],vti[32],vvp[32],vch[32];
    snprintf(vA,32,"%.6f",A); snprintf(vR,32,"%.6f",R_tube);
    snprintf(vti,32,"%.6f",theta_init_frac); snprintf(vvp,32,"%s",vp_names[v_profile]);
    snprintf(vch,32,"%s",chirality);
    const char *keys[]={"A","R_tube","theta_init","v_profile","chirality","S_pred"};
    char vsp[32]; snprintf(vsp,32,"%.6f",S_pred);
    const char *vals[]={vA,vR,vti,vvp,vch,vsp};
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
        void *cols[12];
        int es = (precision==0)?2:4;
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

    for (int a=0;a<3;a++){free(phi[a]);free(phi_vel[a]);free(theta[a]);free(theta_vel[a]);}
    return 0;
}
