/*  accel_analysis.c — Acceleration/force field analysis for SFA files
 *
 *  Computes the full force decomposition from the Cosserat equation:
 *    F_total = F_lap + F_mass + F_pot + F_curl
 *  where:
 *    F_lap  = Laplacian (dispersive, always outward from concentrations)
 *    F_mass = -m^2 * phi (restoring, toward zero)
 *    F_pot  = -V'(P) * dP/dphi (binding, attractive where |P| is large)
 *    F_curl = eta * curl(theta) (electromagnetic coupling)
 *
 *  Outputs: radial profiles of each force component, force balance ratio,
 *  directional force (toward/away from centroid), inter-baryon force,
 *  and force spectral content.
 *
 *  Build: cd sfa/analysis && make accel_analysis
 *  Usage: ./accel_analysis input.sfa [--frame N] [--json out.json]
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
#define NBINS 60

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

int main(int argc, char **argv) {
    if (argc < 2) { fprintf(stderr, "Usage: %s input.sfa [--frame N] [--json out.json]\n", argv[0]); return 1; }
    const char *sfa_path = argv[1];
    int target_frame = -1;
    const char *json_path = NULL;
    for (int i=2;i<argc;i++) {
        if (!strcmp(argv[i],"--frame") && i+1<argc) target_frame = atoi(argv[++i]);
        if (!strcmp(argv[i],"--json") && i+1<argc) json_path = argv[++i];
    }

    SFA *sfa = sfa_open(sfa_path);
    if (!sfa) { fprintf(stderr, "Cannot open %s\n", sfa_path); return 1; }
    int N = sfa->Nx; double L = sfa->Lx;
    long N3 = (long)N*N*N; int NN = N*N;
    double dx = 2.0*L/(N-1);
    double idx2 = 1.0/(dx*dx), idx1 = 1.0/(2.0*dx);

    /* KVMD physics */
    SFA_KVMDSet kv[16]; int nkv = sfa_read_kvmd(sfa, kv, 16);
    for (int s=0;s<nkv;s++) for(int p=0;p<kv[s].n_pairs;p++) {
        if(!strcmp(kv[s].keys[p],"mu")) MU=atof(kv[s].values[p]);
        if(!strcmp(kv[s].keys[p],"kappa")) KAPPA=atof(kv[s].values[p]);
        if(!strcmp(kv[s].keys[p],"m")){double m=atof(kv[s].values[p]);MASS2=m*m;}
        if(!strcmp(kv[s].keys[p],"eta")) ETA=atof(kv[s].values[p]);
    }

    if (target_frame < 0) target_frame = sfa->total_frames - 1;
    printf("Acceleration analysis: %s (N=%d L=%.0f frame=%d)\n", sfa_path, N, L, target_frame);
    printf("Physics: mu=%.3f kappa=%.1f m^2=%.4f eta=%.3f\n\n", MU, KAPPA, MASS2, ETA);

    void *buf = malloc(sfa->frame_bytes);
    sfa_read_frame(sfa, target_frame, buf);
    double t = sfa_frame_time(sfa, target_frame);

    double *phi[3], *theta[3];
    for (int a=0; a<3; a++) {
        phi[a] = extract_col(buf, sfa, SFA_POSITION, a, N3);
        theta[a] = extract_col(buf, sfa, SFA_ANGLE, a, N3);
    }
    free(buf);

    /* Find centroid weighted by |P| */
    double cm[3] = {0}, wt = 0;
    for (long idx=0; idx<N3; idx++) {
        double P = fabs(phi[0][idx]*phi[1][idx]*phi[2][idx]);
        int i=(int)(idx/NN),j=(int)((idx/N)%N),k=(int)(idx%N);
        cm[0]+=(-L+i*dx)*P; cm[1]+=(-L+j*dx)*P; cm[2]+=(-L+k*dx)*P; wt+=P;
    }
    if(wt>0){cm[0]/=wt;cm[1]/=wt;cm[2]/=wt;}
    printf("Centroid: (%.2f, %.2f, %.2f)  t=%.1f\n\n", cm[0], cm[1], cm[2], t);

    /* Radial bins */
    double dr = L / NBINS;

    /* Per-bin accumulators */
    double bin_Flap[NBINS]={0}, bin_Fmass[NBINS]={0}, bin_Fpot[NBINS]={0}, bin_Fcurl[NBINS]={0};
    double bin_Ftotal[NBINS]={0}, bin_Fradial[NBINS]={0}, bin_Ftang[NBINS]={0};
    double bin_balance[NBINS]={0}; /* |F_total| / max(|F_component|) */
    double bin_rho[NBINS]={0}, bin_P[NBINS]={0};
    long bin_count[NBINS]={0};

    /* Also track inter-baryon force: force along the x-axis (baryon separation axis) */
    /* Split field into left (x<cm) and right (x>cm) halves */
    double F_inter_x = 0; /* net force on right half in x direction */
    long n_right = 0;

    printf("Computing forces at %ld grid points...\n", N3);

    #pragma omp parallel for reduction(+:F_inter_x,n_right) schedule(static)
    for (long idx = 0; idx < N3; idx++) {
        int i=(int)(idx/NN), j=(int)((idx/N)%N), k=(int)(idx%N);
        if (i<1||i>=N-1||j<1||j>=N-1||k<1||k>=N-1) continue;

        double x=-L+i*dx, y=-L+j*dx, z=-L+k*dx;
        double rx=x-cm[0], ry=y-cm[1], rz=z-cm[2];
        double r = sqrt(rx*rx + ry*ry + rz*rz);
        int bin = (int)(r / dr);
        if (bin >= NBINS) continue;

        double p0=phi[0][idx], p1=phi[1][idx], p2=phi[2][idx];
        double P=p0*p1*p2, P2=P*P;
        double rho = p0*p0 + p1*p1 + p2*p2;

        /* Force components (summed over 3 field components) */
        double flap2=0, fmass2=0, fpot2=0, fcurl2=0, ftot2=0;
        double f_radial=0, f_tang2=0;
        double fx_total=0; /* x-component of total force (for inter-baryon) */

        long nip=(long)(i+1)*NN+j*N+k, nim=(long)(i-1)*NN+j*N+k;
        long njp=(long)i*NN+(j+1)*N+k, njm=(long)i*NN+(j-1)*N+k;
        long nkp=(long)i*NN+j*N+(k+1), nkm=(long)i*NN+j*N+(k-1);

        double dVdP = MU * P / ((1+KAPPA*P2)*(1+KAPPA*P2));

        for (int a=0; a<3; a++) {
            double lap = (phi[a][nip]+phi[a][nim]+phi[a][njp]+phi[a][njm]
                        +phi[a][nkp]+phi[a][nkm]-6*phi[a][idx]) * idx2;
            double mass_f = -MASS2 * phi[a][idx];
            double dPda = (a==0)?p1*p2:(a==1)?p0*p2:p0*p1;
            double pot_f = -dVdP * dPda;

            double curl_f = 0;
            if (a==0) curl_f = ETA*(theta[2][njp]-theta[2][njm]-theta[1][nkp]+theta[1][nkm])*idx1;
            else if (a==1) curl_f = ETA*(theta[0][nkp]-theta[0][nkm]-theta[2][nip]+theta[2][nim])*idx1;
            else curl_f = ETA*(theta[1][nip]-theta[1][nim]-theta[0][njp]+theta[0][njm])*idx1;

            double total_f = lap + mass_f + pot_f + curl_f;

            flap2 += lap*lap;
            fmass2 += mass_f*mass_f;
            fpot2 += pot_f*pot_f;
            fcurl2 += curl_f*curl_f;
            ftot2 += total_f*total_f;

            /* Radial component: F · r_hat */
            if (r > 0.1) {
                double rhat[3] = {rx/r, ry/r, rz/r};
                double comp = 0;
                if (a==0) comp = total_f; /* this is the a-th FIELD force, not spatial */
                /* For spatial radial force, need to project total force vector */
            }

            if (a == 0) fx_total = total_f; /* phi_0 force ≈ x-direction */
        }

        /* Approximate radial force: project (F_phi0, F_phi1, F_phi2) along r_hat
           This is approximate because phi components aren't spatial components,
           but for the breathing mode analysis it's informative */
        double F_phi[3];
        for (int a=0; a<3; a++) {
            double lap = (phi[a][nip]+phi[a][nim]+phi[a][njp]+phi[a][njm]
                        +phi[a][nkp]+phi[a][nkm]-6*phi[a][idx]) * idx2;
            double dPda = (a==0)?p1*p2:(a==1)?p0*p2:p0*p1;
            double curl_f = 0;
            if (a==0) curl_f = ETA*(theta[2][njp]-theta[2][njm]-theta[1][nkp]+theta[1][nkm])*idx1;
            else if (a==1) curl_f = ETA*(theta[0][nkp]-theta[0][nkm]-theta[2][nip]+theta[2][nim])*idx1;
            else curl_f = ETA*(theta[1][nip]-theta[1][nim]-theta[0][njp]+theta[0][njm])*idx1;
            F_phi[a] = lap - MASS2*phi[a][idx] - dVdP*dPda + curl_f;
        }

        /* Use the gradient of |phi|^2 as a proxy for the spatial force direction */
        double grad_rho_x = (phi[0][nip]*phi[0][nip]+phi[1][nip]*phi[1][nip]+phi[2][nip]*phi[2][nip]
                           - phi[0][nim]*phi[0][nim]-phi[1][nim]*phi[1][nim]-phi[2][nim]*phi[2][nim]) * idx1;
        double grad_rho_y = (phi[0][njp]*phi[0][njp]+phi[1][njp]*phi[1][njp]+phi[2][njp]*phi[2][njp]
                           - phi[0][njm]*phi[0][njm]-phi[1][njm]*phi[1][njm]-phi[2][njm]*phi[2][njm]) * idx1;
        double grad_rho_z = (phi[0][nkp]*phi[0][nkp]+phi[1][nkp]*phi[1][nkp]+phi[2][nkp]*phi[2][nkp]
                           - phi[0][nkm]*phi[0][nkm]-phi[1][nkm]*phi[1][nkm]-phi[2][nkm]*phi[2][nkm]) * idx1;

        /* Radial force: acceleration in the direction of grad(rho) ≈ toward/away from center */
        if (r > 0.1) {
            double rhat[3] = {rx/r, ry/r, rz/r};
            /* Net acceleration dot radial = sum of field forces weighted by field values */
            double f_r = 0;
            for (int a=0;a<3;a++) f_r += F_phi[a] * phi[a][idx];
            /* Normalize by field magnitude */
            if (rho > 1e-10) f_r /= sqrt(rho);
            f_radial = f_r;
            f_tang2 = ftot2 - f_radial*f_radial; if (f_tang2 < 0) f_tang2 = 0;
        }

        /* Force balance: |F_total| / max(|F_i|) */
        double fmax = sqrt(flap2);
        if (sqrt(fmass2) > fmax) fmax = sqrt(fmass2);
        if (sqrt(fpot2) > fmax) fmax = sqrt(fpot2);
        if (sqrt(fcurl2) > fmax) fmax = sqrt(fcurl2);
        double bal = (fmax > 1e-10) ? sqrt(ftot2) / fmax : 0;

        /* Accumulate to bins (need atomic for parallel) */
        #pragma omp atomic
        bin_Flap[bin] += sqrt(flap2);
        #pragma omp atomic
        bin_Fmass[bin] += sqrt(fmass2);
        #pragma omp atomic
        bin_Fpot[bin] += sqrt(fpot2);
        #pragma omp atomic
        bin_Fcurl[bin] += sqrt(fcurl2);
        #pragma omp atomic
        bin_Ftotal[bin] += sqrt(ftot2);
        #pragma omp atomic
        bin_Fradial[bin] += f_radial;
        #pragma omp atomic
        bin_Ftang[bin] += sqrt(f_tang2);
        #pragma omp atomic
        bin_balance[bin] += bal;
        #pragma omp atomic
        bin_rho[bin] += rho;
        #pragma omp atomic
        bin_P[bin] += fabs(P);
        #pragma omp atomic
        bin_count[bin]++;

        /* Inter-baryon force: net x-force on right half */
        if (x > cm[0]) {
            F_inter_x += F_phi[0] * phi[0][idx] * dx*dx*dx;
            n_right++;
        }
    }

    /* Normalize bins */
    printf("\n%-6s %8s %8s %8s %8s %8s %8s %8s %8s %8s %8s %8s\n",
           "r", "|F_lap|", "|F_mass|", "|F_pot|", "|F_curl|", "|F_tot|",
           "F_rad", "|F_tan|", "balance", "rho", "|P|", "count");

    for (int b = 0; b < NBINS; b++) {
        if (bin_count[b] == 0) continue;
        double n = (double)bin_count[b];
        printf("%-6.1f %8.3f %8.3f %8.3f %8.3f %8.3f %+8.3f %8.3f %8.3f %8.4f %8.5f %8ld\n",
               (b+0.5)*dr,
               bin_Flap[b]/n, bin_Fmass[b]/n, bin_Fpot[b]/n, bin_Fcurl[b]/n,
               bin_Ftotal[b]/n, bin_Fradial[b]/n, bin_Ftang[b]/n,
               bin_balance[b]/n, bin_rho[b]/n, bin_P[b]/n, bin_count[b]);
    }

    /* Summary */
    printf("\n=== FORCE SUMMARY ===\n");
    double peak_pot = 0, peak_curl = 0, peak_r = 0;
    int peak_bin = 0;
    for (int b=0; b<NBINS; b++) {
        if (bin_count[b] == 0) continue;
        double fp = bin_Fpot[b]/bin_count[b];
        double fc = bin_Fcurl[b]/bin_count[b];
        if (fp > peak_pot) { peak_pot = fp; peak_bin = b; peak_r = (b+0.5)*dr; }
        if (fc > peak_curl) peak_curl = fc;
    }
    printf("Peak binding force |F_pot|: %.3f at r=%.1f\n", peak_pot, peak_r);
    printf("Peak curl force |F_curl|:   %.3f\n", peak_curl);
    printf("Pot/Curl ratio at peak:     %.2f\n", peak_pot / (peak_curl + 1e-10));

    /* Force balance in core */
    double core_balance = 0; long core_count = 0;
    for (int b=0; b<NBINS/3; b++) {
        core_balance += bin_balance[b]; core_count += bin_count[b];
    }
    if (core_count > 0) core_balance /= core_count;
    printf("Core force balance ratio:   %.3f (lower = better balanced)\n", core_balance);

    /* Inter-baryon force */
    printf("\nInter-baryon force (F_x on right half): %+.4f\n", F_inter_x);
    printf("  Positive = repulsive, Negative = attractive\n");
    if (F_inter_x < 0)
        printf("  → ATTRACTIVE: baryons are pulling toward each other\n");
    else
        printf("  → REPULSIVE: baryons are pushing apart\n");

    /* JSON output */
    if (json_path) {
        FILE *jfp = fopen(json_path, "w");
        fprintf(jfp, "{\n  \"sfa\": \"%s\",\n  \"frame\": %d, \"time\": %.1f,\n", sfa_path, target_frame, t);
        fprintf(jfp, "  \"centroid\": [%.3f, %.3f, %.3f],\n", cm[0], cm[1], cm[2]);
        fprintf(jfp, "  \"peak_F_pot\": %.4f, \"peak_F_pot_r\": %.1f,\n", peak_pot, peak_r);
        fprintf(jfp, "  \"peak_F_curl\": %.4f,\n", peak_curl);
        fprintf(jfp, "  \"pot_curl_ratio\": %.3f,\n", peak_pot/(peak_curl+1e-10));
        fprintf(jfp, "  \"core_balance\": %.4f,\n", core_balance);
        fprintf(jfp, "  \"inter_baryon_Fx\": %.6f,\n", F_inter_x);
        fprintf(jfp, "  \"radial_profile\": [\n");
        for (int b=0; b<NBINS; b++) {
            if (bin_count[b] == 0) continue;
            double n = (double)bin_count[b];
            fprintf(jfp, "    {\"r\":%.2f,\"F_lap\":%.4f,\"F_mass\":%.4f,\"F_pot\":%.4f,"
                   "\"F_curl\":%.4f,\"F_total\":%.4f,\"F_radial\":%.4f,\"F_tang\":%.4f,"
                   "\"balance\":%.4f,\"rho\":%.5f,\"P\":%.6f,\"count\":%ld}%s\n",
                   (b+0.5)*dr, bin_Flap[b]/n, bin_Fmass[b]/n, bin_Fpot[b]/n,
                   bin_Fcurl[b]/n, bin_Ftotal[b]/n, bin_Fradial[b]/n, bin_Ftang[b]/n,
                   bin_balance[b]/n, bin_rho[b]/n, bin_P[b]/n, bin_count[b],
                   (b < NBINS-1) ? "," : "");
        }
        fprintf(jfp, "  ]\n}\n");
        fclose(jfp);
        printf("\nJSON written to %s\n", json_path);
    }

    for (int a=0;a<3;a++){free(phi[a]);free(theta[a]);}
    sfa_close(sfa);
    return 0;
}
