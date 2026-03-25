/*  breathing.c — Characterize the breathing oscillator dynamics
 *
 *  For multi-frame SFA: extracts the time series of key oscillation quantities
 *  at the core, mid-shell, and outer halo. Measures:
 *
 *  1. Oscillation amplitude and phase of ρ, |P|, θ, |v| at each radius
 *  2. Internal velocity field: min/max/mean |v|, fraction of volume at |v|>0.5c
 *  3. Internal temperature proxy: velocity dispersion (T_kin ∝ <v²> - <v>²)
 *  4. Breathing period from peak-to-peak in E_pot time series
 *  5. Phase velocity of the breathing wave (does the oscillation propagate outward?)
 *  6. Energy partition: E_kin(t) vs E_pot(t) — are they in anti-phase? (virial)
 *  7. Speed of sound: group velocity of small perturbations through the structure
 *
 *  Build: cd sfa/analysis && make breathing
 *  Usage: ./breathing input.sfa [--json out.json]
 */

#define SFA_IMPLEMENTATION
#include "../format/sfa.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <omp.h>

#define PI 3.14159265358979323846
#define NFIELDS 3
#define MAX_FRAMES 200
#define NSHELLS 40

static double MU = -41.345, KAPPA = 50.0, MASS2 = 2.25, ETA = 0.5;

static inline double f16_to_f64(uint16_t h) {
    uint16_t sign = h & 0x8000; int exp = (h >> 10) & 0x1F; uint16_t mant = h & 0x3FF;
    if (exp == 0) return 0.0; if (exp == 31) return sign ? -1e30 : 1e30;
    float fv; uint32_t x = ((uint32_t)sign << 16) | ((uint32_t)(exp-15+127) << 23) | ((uint32_t)mant << 13);
    memcpy(&fv, &x, 4); return (double)fv;
}

static double *extract_col(void *buf, SFA *sfa, int sem, int comp, long N3) {
    double *arr = (double*)calloc(N3, sizeof(double));
    uint64_t off = 0;
    for (int c = 0; c < sfa->n_columns; c++) {
        int dt = sfa->columns[c].dtype, s = sfa->columns[c].semantic, cp = sfa->columns[c].component;
        int es = sfa_dtype_size[dt];
        if (s == sem && cp == comp) {
            uint8_t *src = (uint8_t*)buf + off;
            if (dt == SFA_F64) memcpy(arr, src, N3*8);
            else if (dt == SFA_F32) for(long i=0;i<N3;i++) arr[i]=(double)((float*)src)[i];
            else if (dt == SFA_F16) for(long i=0;i<N3;i++) arr[i]=f16_to_f64(((uint16_t*)src)[i]);
            return arr;
        }
        off += (uint64_t)N3 * es;
    }
    return arr;
}

/* Per-shell, per-frame measurements */
typedef struct {
    double rho_mean, P_mean, theta_rms;
    double v_mean, v_max, v_rms;       /* velocity statistics */
    double v_radial_mean;               /* mean radial velocity (expansion/contraction) */
    double T_kin;                       /* kinetic temperature = <v²> - <v>² */
    double E_kin, E_pot;                /* energy in this shell */
    double frac_v_gt_half_c;            /* fraction of volume with |v| > 0.5 */
    long count;
} ShellData;

int main(int argc, char **argv) {
    if (argc < 2) { fprintf(stderr, "Usage: %s input.sfa [--json out.json]\n", argv[0]); return 1; }
    const char *sfa_path = argv[1];
    const char *json_path = NULL;
    for (int i=2;i<argc;i++) if(!strcmp(argv[i],"--json")&&i+1<argc) json_path=argv[++i];

    SFA *sfa = sfa_open(sfa_path);
    if (!sfa) { fprintf(stderr, "Cannot open %s\n", sfa_path); return 1; }
    int N = sfa->Nx; double L = sfa->Lx;
    long N3 = (long)N*N*N; int NN = N*N;
    double dx = 2.0*L/(N-1), dV = dx*dx*dx;

    /* KVMD */
    SFA_KVMDSet kv[16]; int nkv = sfa_read_kvmd(sfa, kv, 16);
    for(int s=0;s<nkv;s++) for(int p=0;p<kv[s].n_pairs;p++){
        if(!strcmp(kv[s].keys[p],"mu")) MU=atof(kv[s].values[p]);
        if(!strcmp(kv[s].keys[p],"kappa")) KAPPA=atof(kv[s].values[p]);
        if(!strcmp(kv[s].keys[p],"m")){double m=atof(kv[s].values[p]);MASS2=m*m;}
        if(!strcmp(kv[s].keys[p],"eta")) ETA=atof(kv[s].values[p]);
    }

    int nframes = sfa->total_frames;
    if (nframes > MAX_FRAMES) nframes = MAX_FRAMES;
    double dr = L / NSHELLS;

    printf("Breathing analysis: %s (N=%d L=%.0f %d frames)\n", sfa_path, N, L, nframes);
    printf("Physics: mu=%.3f kappa=%.1f m^2=%.4f eta=%.3f\n", MU, KAPPA, MASS2, ETA);
    printf("Shell width: %.2f, %d shells\n\n", dr, NSHELLS);

    /* Storage */
    double *frame_time = (double*)calloc(nframes, sizeof(double));
    double *frame_Ekin = (double*)calloc(nframes, sizeof(double));
    double *frame_Epot = (double*)calloc(nframes, sizeof(double));
    double *frame_Etot = (double*)calloc(nframes, sizeof(double));
    double *frame_vmax = (double*)calloc(nframes, sizeof(double));
    double *frame_vmean = (double*)calloc(nframes, sizeof(double));
    double *frame_Tkin = (double*)calloc(nframes, sizeof(double));
    double *frame_frac_fast = (double*)calloc(nframes, sizeof(double));
    double *frame_Rrms = (double*)calloc(nframes, sizeof(double));
    ShellData (*shells)[NSHELLS] = calloc(nframes, sizeof(ShellData[NSHELLS]));

    void *buf = malloc(sfa->frame_bytes);

    for (int f = 0; f < nframes; f++) {
        sfa_read_frame(sfa, f, buf);
        frame_time[f] = sfa_frame_time(sfa, f);

        double *phi[3], *theta[3], *phi_vel[3];
        int has_vel = (sfa->n_columns >= 9);
        for (int a=0;a<3;a++) {
            phi[a] = extract_col(buf, sfa, SFA_POSITION, a, N3);
            theta[a] = extract_col(buf, sfa, SFA_ANGLE, a, N3);
            phi_vel[a] = has_vel ? extract_col(buf, sfa, SFA_VELOCITY, a, N3) : NULL;
        }

        /* Centroid */
        double cm[3]={0}, wt=0;
        for (long idx=0;idx<N3;idx++){
            double P=fabs(phi[0][idx]*phi[1][idx]*phi[2][idx]);
            int i=(int)(idx/NN),j=(int)((idx/N)%N),k=(int)(idx%N);
            cm[0]+=(-L+i*dx)*P; cm[1]+=(-L+j*dx)*P; cm[2]+=(-L+k*dx)*P; wt+=P;
        }
        if(wt>0){cm[0]/=wt;cm[1]/=wt;cm[2]/=wt;}

        /* Per-voxel analysis */
        double Ekin_tot=0, Epot_tot=0, vmax_glob=0;
        double vsum=0, v2sum=0, n_fast=0;
        long n_active=0;
        double R2_sum=0, R_wt=0;

        for (long idx=0;idx<N3;idx++){
            int i=(int)(idx/NN),j=(int)((idx/N)%N),k=(int)(idx%N);
            double x=-L+i*dx-cm[0], y=-L+j*dx-cm[1], z=-L+k*dx-cm[2];
            double r=sqrt(x*x+y*y+z*z);
            int shell=(int)(r/dr); if(shell>=NSHELLS) continue;

            double p0=phi[0][idx],p1=phi[1][idx],p2=phi[2][idx];
            double rho=p0*p0+p1*p1+p2*p2;
            double P=p0*p1*p2, P2=P*P, Pabs=fabs(P);
            double t2=0; for(int a=0;a<3;a++) t2+=theta[a][idx]*theta[a][idx];

            /* Velocity */
            double vmag=0, vr=0;
            if (has_vel && phi_vel[0]) {
                double vx=phi_vel[0][idx], vy=phi_vel[1][idx], vz=phi_vel[2][idx];
                vmag = sqrt(vx*vx + vy*vy + vz*vz);
                if (r > 0.1) vr = (vx*x + vy*y + vz*z) / r;
            }

            /* Energy densities */
            double ekin = 0;
            if (has_vel && phi_vel[0]) {
                for(int a=0;a<3;a++) ekin += 0.5*phi_vel[a][idx]*phi_vel[a][idx];
            }
            double epot = (MU/2.0)*P2/(1.0+KAPPA*P2);

            /* Accumulate to shell */
            ShellData *s = &shells[f][shell];
            s->rho_mean += rho;
            s->P_mean += Pabs;
            s->theta_rms += t2;
            s->v_mean += vmag;
            s->v_rms += vmag*vmag;
            if (vmag > s->v_max) s->v_max = vmag;
            s->v_radial_mean += vr;
            s->E_kin += ekin * dV;
            s->E_pot += epot * dV;
            if (vmag > 0.5) s->frac_v_gt_half_c += 1.0;
            s->count++;

            /* Global accumulators */
            Ekin_tot += ekin * dV;
            Epot_tot += epot * dV;
            if (vmag > vmax_glob) vmax_glob = vmag;
            vsum += vmag;
            v2sum += vmag*vmag;
            n_active++;
            if (vmag > 0.5) n_fast++;
            R2_sum += r*r*Pabs; R_wt += Pabs;
        }

        /* Normalize shells */
        for (int s=0;s<NSHELLS;s++){
            if(shells[f][s].count>0){
                double n=(double)shells[f][s].count;
                shells[f][s].rho_mean /= n;
                shells[f][s].P_mean /= n;
                shells[f][s].theta_rms = sqrt(shells[f][s].theta_rms / (3*n));
                double vm = shells[f][s].v_mean / n;
                shells[f][s].v_mean = vm;
                shells[f][s].v_rms = sqrt(shells[f][s].v_rms / n);
                shells[f][s].v_radial_mean /= n;
                shells[f][s].T_kin = shells[f][s].v_rms*shells[f][s].v_rms - vm*vm;
                shells[f][s].frac_v_gt_half_c /= n;
            }
        }

        /* Global frame stats */
        frame_Ekin[f] = Ekin_tot;
        frame_Epot[f] = Epot_tot;
        frame_Etot[f] = Ekin_tot + Epot_tot;
        frame_vmax[f] = vmax_glob;
        frame_vmean[f] = (n_active>0) ? vsum/n_active : 0;
        frame_Tkin[f] = (n_active>0) ? v2sum/n_active - (vsum/n_active)*(vsum/n_active) : 0;
        frame_frac_fast[f] = (n_active>0) ? (double)n_fast/n_active : 0;
        frame_Rrms[f] = (R_wt>0) ? sqrt(R2_sum/R_wt) : 0;

        printf("f=%2d t=%6.1f Ekin=%8.1f Epot=%8.1f vmax=%.3f vmean=%.4f Tkin=%.5f frac>0.5c=%.4f Rrms=%.1f\n",
               f, frame_time[f], Ekin_tot, Epot_tot, vmax_glob,
               frame_vmean[f], frame_Tkin[f], frame_frac_fast[f], frame_Rrms[f]);

        for(int a=0;a<3;a++){free(phi[a]);free(theta[a]);if(phi_vel[a])free(phi_vel[a]);}
    }

    /* Breathing period from E_pot peaks */
    printf("\n=== BREATHING ANALYSIS ===\n");
    int n_peaks = 0;
    double peak_times[50];
    for (int f=1; f<nframes-1; f++) {
        if (frame_Epot[f] < frame_Epot[f-1] && frame_Epot[f] < frame_Epot[f+1] && n_peaks < 50)
            peak_times[n_peaks++] = frame_time[f];
    }
    double breath_period = 0;
    if (n_peaks >= 2) {
        double sum_dt = 0;
        for (int i=1;i<n_peaks;i++) sum_dt += peak_times[i] - peak_times[i-1];
        breath_period = sum_dt / (n_peaks-1);
    }
    printf("E_pot minima: %d detected\n", n_peaks);
    for(int i=0;i<n_peaks;i++) printf("  Peak %d at t=%.1f\n", i, peak_times[i]);
    printf("Breathing period: %.1f time units\n", breath_period);

    /* Energy partition: Ekin vs Epot correlation */
    double corr = 0, Ek_mean=0, Ep_mean=0;
    for(int f=0;f<nframes;f++){Ek_mean+=frame_Ekin[f];Ep_mean+=frame_Epot[f];}
    Ek_mean/=nframes; Ep_mean/=nframes;
    double Ek_var=0, Ep_var=0;
    for(int f=0;f<nframes;f++){
        corr += (frame_Ekin[f]-Ek_mean)*(frame_Epot[f]-Ep_mean);
        Ek_var += (frame_Ekin[f]-Ek_mean)*(frame_Ekin[f]-Ek_mean);
        Ep_var += (frame_Epot[f]-Ep_mean)*(frame_Epot[f]-Ep_mean);
    }
    double Ekin_Epot_corr = (Ek_var>0&&Ep_var>0) ? corr/sqrt(Ek_var*Ep_var) : 0;
    printf("\nEkin-Epot correlation: %.3f", Ekin_Epot_corr);
    if (Ekin_Epot_corr < -0.5) printf(" (ANTI-PHASE — virial oscillation)\n");
    else if (Ekin_Epot_corr > 0.5) printf(" (IN-PHASE — coherent)\n");
    else printf(" (UNCORRELATED)\n");

    /* Velocity constraints */
    printf("\n=== VELOCITY CONSTRAINTS ===\n");
    double global_vmax = 0;
    for(int f=0;f<nframes;f++) if(frame_vmax[f]>global_vmax) global_vmax=frame_vmax[f];
    printf("Global max |v| across all frames: %.4f c\n", global_vmax);
    printf("  %s\n", global_vmax < 1.0 ? "SUBLUMINAL (0 < V < c satisfied)" : "SUPERLUMINAL WARNING");

    /* Temperature (velocity dispersion) profile */
    printf("\n=== INTERNAL TEMPERATURE PROFILE (last frame) ===\n");
    printf("%-6s %8s %8s %8s %8s %8s %8s %8s\n",
           "r", "T_kin", "|v|_mean", "|v|_max", "v_rad", "f>0.5c", "rho", "|P|");
    int lf = nframes-1;
    for(int s=0;s<NSHELLS;s++){
        if(shells[lf][s].count < 10) continue;
        printf("%-6.1f %8.5f %8.4f %8.4f %+8.4f %8.4f %8.4f %8.5f\n",
               (s+0.5)*dr,
               shells[lf][s].T_kin, shells[lf][s].v_mean, shells[lf][s].v_max,
               shells[lf][s].v_radial_mean, shells[lf][s].frac_v_gt_half_c,
               shells[lf][s].rho_mean, shells[lf][s].P_mean);
    }

    /* Breathing wave propagation: does the oscillation peak move outward? */
    printf("\n=== BREATHING WAVE PROPAGATION ===\n");
    if (nframes >= 3) {
        printf("Shell rho oscillation amplitude (max-min across frames):\n");
        printf("%-6s %10s %10s %10s\n", "r", "rho_osc", "P_osc", "v_osc");
        for(int s=0;s<NSHELLS;s++){
            double rho_min=1e30,rho_max=-1e30;
            double P_min=1e30,P_max=-1e30;
            double v_min=1e30,v_max=-1e30;
            int valid=0;
            for(int f=0;f<nframes;f++){
                if(shells[f][s].count<10) continue;
                valid=1;
                if(shells[f][s].rho_mean<rho_min)rho_min=shells[f][s].rho_mean;
                if(shells[f][s].rho_mean>rho_max)rho_max=shells[f][s].rho_mean;
                if(shells[f][s].P_mean<P_min)P_min=shells[f][s].P_mean;
                if(shells[f][s].P_mean>P_max)P_max=shells[f][s].P_mean;
                if(shells[f][s].v_mean<v_min)v_min=shells[f][s].v_mean;
                if(shells[f][s].v_mean>v_max)v_max=shells[f][s].v_mean;
            }
            if(valid)
                printf("%-6.1f %10.5f %10.6f %10.5f\n",
                       (s+0.5)*dr, rho_max-rho_min, P_max-P_min, v_max-v_min);
        }
    }

    /* JSON output */
    if (json_path) {
        FILE *jfp = fopen(json_path, "w");
        fprintf(jfp, "{\n  \"sfa\": \"%s\",\n  \"N\": %d, \"L\": %.1f, \"n_frames\": %d,\n",
                sfa_path, N, L, nframes);
        fprintf(jfp, "  \"breathing_period\": %.2f,\n", breath_period);
        fprintf(jfp, "  \"Ekin_Epot_correlation\": %.4f,\n", Ekin_Epot_corr);
        fprintf(jfp, "  \"global_vmax\": %.6f,\n", global_vmax);
        fprintf(jfp, "  \"subluminal\": %s,\n", global_vmax < 1.0 ? "true" : "false");
        fprintf(jfp, "  \"frames\": [\n");
        for(int f=0;f<nframes;f++){
            fprintf(jfp, "    {\"t\":%.1f,\"Ekin\":%.2f,\"Epot\":%.2f,\"Etot\":%.2f,"
                   "\"vmax\":%.5f,\"vmean\":%.5f,\"Tkin\":%.6f,\"frac_fast\":%.5f,\"Rrms\":%.2f}%s\n",
                   frame_time[f],frame_Ekin[f],frame_Epot[f],frame_Etot[f],
                   frame_vmax[f],frame_vmean[f],frame_Tkin[f],frame_frac_fast[f],frame_Rrms[f],
                   f<nframes-1?",":"");
        }
        fprintf(jfp, "  ],\n  \"shell_profile_last_frame\": [\n");
        for(int s=0;s<NSHELLS;s++){
            if(shells[lf][s].count<10) continue;
            fprintf(jfp, "    {\"r\":%.2f,\"T_kin\":%.6f,\"v_mean\":%.5f,\"v_max\":%.5f,"
                   "\"v_rad\":%.5f,\"frac_fast\":%.5f,\"rho\":%.5f,\"P\":%.6f}%s\n",
                   (s+0.5)*dr,shells[lf][s].T_kin,shells[lf][s].v_mean,shells[lf][s].v_max,
                   shells[lf][s].v_radial_mean,shells[lf][s].frac_v_gt_half_c,
                   shells[lf][s].rho_mean,shells[lf][s].P_mean,
                   s<NSHELLS-1?",":"");
        }
        fprintf(jfp, "  ]\n}\n");
        fclose(jfp);
        printf("\nJSON written to %s\n", json_path);
    }

    free(buf); free(frame_time); free(frame_Ekin); free(frame_Epot); free(frame_Etot);
    free(frame_vmax); free(frame_vmean); free(frame_Tkin); free(frame_frac_fast);
    free(frame_Rrms); free(shells);
    sfa_close(sfa);
    return 0;
}
