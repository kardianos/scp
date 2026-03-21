/*  wedge_standalone.c — Run the electron orbital independently
 *
 *  The braid is FROZEN (static V_eff). Only the 1D Schrödinger wedge evolves.
 *  This allows running for T=10^12 to see orbital dynamics.
 *
 *  Build: gcc -O3 -o wedge_standalone src/wedge_standalone.c -lm
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <sys/stat.h>

#define PI 3.14159265358979323846

/* Parameters */
static double HBAR    = 22727.0;
static double M_EFF   = 1535.0;
static double V_DEPTH = 1.27;
static double N_POWER = 1.189;
static double R_BRAID = 5.0;
static double R_MIN   = 100.0;  /* start far from core — no orbital physics below this */
static double R_MAX   = 500000.0;

typedef struct {
    int Nr;
    double *r, *dr, *V;
    double *psi_re, *psi_im;
    double dt;
    /* CN work arrays */
    double *u_re, *u_im, *Hu_re, *Hu_im;
    double *rhs_re, *rhs_im;
    double *a_coeff, *b_coeff, *c_coeff;
    double *d_re, *d_im, *f_re, *f_im;
} Wedge;

static Wedge *wedge_alloc(int Nr) {
    Wedge *w = calloc(1, sizeof(Wedge));
    w->Nr = Nr;
    w->r = malloc(Nr * sizeof(double));
    w->dr = malloc(Nr * sizeof(double));
    w->V = malloc(Nr * sizeof(double));
    w->psi_re = calloc(Nr, sizeof(double));
    w->psi_im = calloc(Nr, sizeof(double));
    w->u_re = malloc(Nr * sizeof(double));
    w->u_im = malloc(Nr * sizeof(double));
    w->Hu_re = malloc(Nr * sizeof(double));
    w->Hu_im = malloc(Nr * sizeof(double));
    w->rhs_re = malloc(Nr * sizeof(double));
    w->rhs_im = malloc(Nr * sizeof(double));
    w->a_coeff = malloc(Nr * sizeof(double));
    w->b_coeff = malloc(Nr * sizeof(double));
    w->c_coeff = malloc(Nr * sizeof(double));
    w->d_re = malloc(Nr * sizeof(double));
    w->d_im = malloc(Nr * sizeof(double));
    w->f_re = malloc(Nr * sizeof(double));
    w->f_im = malloc(Nr * sizeof(double));

    /* Logarithmic r spacing */
    double log_min = log(R_MIN), log_max = log(R_MAX);
    for (int i = 0; i < Nr; i++) {
        double frac = (double)i / (Nr - 1);
        w->r[i] = exp(log_min + frac * (log_max - log_min));
    }
    for (int i = 0; i < Nr - 1; i++)
        w->dr[i] = w->r[i+1] - w->r[i];
    w->dr[Nr-1] = w->dr[Nr-2];

    /* Potential */
    for (int i = 0; i < Nr; i++) {
        if (w->r[i] < 2.0)
            w->V[i] = V_DEPTH * 10.0 * (2.0/w->r[i] - 1.0);
        else
            w->V[i] = -V_DEPTH / pow(w->r[i]/R_BRAID, N_POWER);
    }

    /* CN tridiagonal coefficients: H = -hbar^2/(2m) d^2/dr^2 + V
     * Using u=r*psi, the radial equation becomes:
     * -hbar^2/(2m) u'' + V(r)*u = E*u
     * Tridiag: a_i u_{i-1} + b_i u_i + c_i u_{i+1} = E u_i */
    double alpha = HBAR*HBAR / (2.0 * M_EFF);
    for (int i = 0; i < Nr; i++) {
        double dr_avg = (i > 0 && i < Nr-1) ?
            0.5*(w->r[i+1] - w->r[i-1]) : w->dr[i < Nr-1 ? i : Nr-2];
        double dr_m = (i > 0) ? (w->r[i] - w->r[i-1]) : w->dr[0];
        double dr_p = (i < Nr-1) ? (w->r[i+1] - w->r[i]) : w->dr[Nr-2];

        w->a_coeff[i] = (i > 0) ? -alpha / (dr_avg * dr_m) : 0;
        w->c_coeff[i] = (i < Nr-1) ? -alpha / (dr_avg * dr_p) : 0;
        w->b_coeff[i] = -w->a_coeff[i] - w->c_coeff[i] + w->V[i];
    }

    return w;
}

static void wedge_init_gaussian(Wedge *w, double r0, double sigma) {
    double norm = 0;
    for (int i = 0; i < w->Nr; i++) {
        double dr = w->r[i] - r0;
        w->psi_re[i] = exp(-dr*dr / (2*sigma*sigma));
        w->psi_im[i] = 0;
        norm += (w->psi_re[i]*w->psi_re[i]) * w->r[i]*w->r[i] * w->dr[i];
    }
    double inv = 1.0 / sqrt(norm);
    for (int i = 0; i < w->Nr; i++) {
        w->psi_re[i] *= inv;
        w->psi_im[i] = 0;
    }
}

static void wedge_step(Wedge *w) {
    const int Nr = w->Nr;
    const double hdt_h = w->dt / (2.0 * HBAR);
    double *u_re=w->u_re, *u_im=w->u_im;
    double *Hu_re=w->Hu_re, *Hu_im=w->Hu_im;
    double *rhs_re=w->rhs_re, *rhs_im=w->rhs_im;
    double *d_re=w->d_re, *d_im=w->d_im;
    double *f_re=w->f_re, *f_im=w->f_im;

    for (int i = 0; i < Nr; i++) {
        u_re[i] = w->r[i] * w->psi_re[i];
        u_im[i] = w->r[i] * w->psi_im[i];
    }

    /* H*u */
    Hu_re[0] = w->b_coeff[0]*u_re[0]; Hu_im[0] = w->b_coeff[0]*u_im[0];
    Hu_re[Nr-1] = w->b_coeff[Nr-1]*u_re[Nr-1]; Hu_im[Nr-1] = w->b_coeff[Nr-1]*u_im[Nr-1];
    for (int i = 1; i < Nr-1; i++) {
        Hu_re[i] = w->a_coeff[i]*u_re[i-1] + w->b_coeff[i]*u_re[i] + w->c_coeff[i]*u_re[i+1];
        Hu_im[i] = w->a_coeff[i]*u_im[i-1] + w->b_coeff[i]*u_im[i] + w->c_coeff[i]*u_im[i+1];
    }

    /* RHS = u - i*(dt/2ℏ)*H*u */
    for (int i = 0; i < Nr; i++) {
        rhs_re[i] = u_re[i] + hdt_h * Hu_im[i];
        rhs_im[i] = u_im[i] - hdt_h * Hu_re[i];
    }

    /* Complex Thomas algorithm */
    memcpy(f_re, rhs_re, Nr*sizeof(double));
    memcpy(f_im, rhs_im, Nr*sizeof(double));
    d_re[0] = 1.0; d_im[0] = hdt_h * w->b_coeff[0];
    for (int i = 1; i < Nr; i++) {
        double L_im = hdt_h * w->a_coeff[i];
        double U_im = hdt_h * w->c_coeff[i-1];
        double den = d_re[i-1]*d_re[i-1] + d_im[i-1]*d_im[i-1];
        double w_re = L_im * d_im[i-1] / den;
        double w_im = L_im * d_re[i-1] / den;
        d_re[i] = 1.0 + w_im * U_im;
        d_im[i] = hdt_h * w->b_coeff[i] - w_re * U_im;
        double fr = f_re[i] - (w_re*f_re[i-1] - w_im*f_im[i-1]);
        double fi = f_im[i] - (w_re*f_im[i-1] + w_im*f_re[i-1]);
        f_re[i] = fr; f_im[i] = fi;
    }
    { int i=Nr-1; double den=d_re[i]*d_re[i]+d_im[i]*d_im[i];
      u_re[i]=(f_re[i]*d_re[i]+f_im[i]*d_im[i])/den;
      u_im[i]=(f_im[i]*d_re[i]-f_re[i]*d_im[i])/den; }
    for (int i=Nr-2; i>=0; i--) {
        double U_im = hdt_h * w->c_coeff[i];
        double t_re = f_re[i] + U_im * u_im[i+1];
        double t_im = f_im[i] - U_im * u_re[i+1];
        double den = d_re[i]*d_re[i] + d_im[i]*d_im[i];
        u_re[i] = (t_re*d_re[i] + t_im*d_im[i]) / den;
        u_im[i] = (t_im*d_re[i] - t_re*d_im[i]) / den;
    }

    for (int i = 0; i < Nr; i++) {
        w->psi_re[i] = u_re[i] / w->r[i];
        w->psi_im[i] = u_im[i] / w->r[i];
    }
    w->psi_re[0]=0; w->psi_im[0]=0;
    w->psi_re[Nr-1]=0; w->psi_im[Nr-1]=0;
}

static void measure(Wedge *w, double *mean_r, double *rms_r, double *norm,
                    double *peak_r, double *energy) {
    double mr=0, mr2=0, n=0, pk=0, pk_r=0;
    for (int i = 0; i < w->Nr; i++) {
        double prob = (w->psi_re[i]*w->psi_re[i] + w->psi_im[i]*w->psi_im[i])
                      * w->r[i]*w->r[i] * w->dr[i];
        mr += w->r[i] * prob;
        mr2 += w->r[i]*w->r[i] * prob;
        n += prob;
        if (prob/w->dr[i] > pk) { pk = prob/w->dr[i]; pk_r = w->r[i]; }
    }
    *mean_r = mr/n; *rms_r = sqrt(mr2/n); *norm = n; *peak_r = pk_r;

    /* Energy via <ψ|H|ψ> */
    double ekin=0, epot=0;
    double alpha = HBAR*HBAR / (2.0 * M_EFF);
    for (int i = 1; i < w->Nr-1; i++) {
        double ri = w->r[i];
        double u_re = ri * w->psi_re[i], u_im = ri * w->psi_im[i];
        double up_re = w->r[i+1]*w->psi_re[i+1], um_re = w->r[i-1]*w->psi_re[i-1];
        double up_im = w->r[i+1]*w->psi_im[i+1], um_im = w->r[i-1]*w->psi_im[i-1];
        double dr_avg = 0.5*(w->r[i+1]-w->r[i-1]);
        double dr_m = w->r[i]-w->r[i-1], dr_p = w->r[i+1]-w->r[i];
        double d2u_re = (up_re/dr_p - u_re*(1/dr_m+1/dr_p) + um_re/dr_m) / dr_avg;
        double d2u_im = (up_im/dr_p - u_im*(1/dr_m+1/dr_p) + um_im/dr_m) / dr_avg;
        double Tpsi_re = -alpha * d2u_re / ri;
        double Tpsi_im = -alpha * d2u_im / ri;
        ekin += (w->psi_re[i]*Tpsi_re + w->psi_im[i]*Tpsi_im) * ri*ri * w->dr[i];
        epot += w->V[i] * (w->psi_re[i]*w->psi_re[i]+w->psi_im[i]*w->psi_im[i]) * ri*ri * w->dr[i];
    }
    *energy = (ekin + epot) / n;
}

int main(int argc, char **argv) {
    int Nr = 4000;
    double T = 1e12;
    double r_frac = 0.3;
    double sigma_frac = 0.2;
    double dt = 1e8;
    int diag_every = 100;
    char outdir[256] = "data/wedge_standalone";

    for (int i=1; i<argc; i++) {
        if (!strcmp(argv[i],"-Nr")) Nr = atoi(argv[++i]);
        else if (!strcmp(argv[i],"-T")) T = atof(argv[++i]);
        else if (!strcmp(argv[i],"-dt")) dt = atof(argv[++i]);
        else if (!strcmp(argv[i],"-rfrac")) r_frac = atof(argv[++i]);
        else if (!strcmp(argv[i],"-sfrac")) sigma_frac = atof(argv[++i]);
        else if (!strcmp(argv[i],"-hbar")) HBAR = atof(argv[++i]);
        else if (!strcmp(argv[i],"-meff")) M_EFF = atof(argv[++i]);
        else if (!strcmp(argv[i],"-Vdepth")) V_DEPTH = atof(argv[++i]);
        else if (!strcmp(argv[i],"-npower")) N_POWER = atof(argv[++i]);
        else if (!strcmp(argv[i],"-diag")) diag_every = atoi(argv[++i]);
        else if (!strcmp(argv[i],"-o")) strncpy(outdir, argv[++i], 255);
    }

    double BOHR = HBAR*HBAR / (M_EFF * V_DEPTH * R_BRAID);
    double T_orbital = 2*PI * BOHR * M_EFF / HBAR;
    long n_steps = (long)(T / dt);

    printf("=== Standalone Wedge: Electron in Frozen Braid Potential ===\n\n");
    printf("ℏ=%.0f  m_eff=%.0f  V_depth=%.3f  n_power=%.3f\n", HBAR, M_EFF, V_DEPTH, N_POWER);
    printf("Bohr radius = %.0f  (%.0f × r_braid)\n", BOHR, BOHR/R_BRAID);
    printf("Orbital period = %.3e\n", T_orbital);
    printf("Nr=%d  dt=%.2e  T=%.2e  steps=%ld\n", Nr, dt, T, n_steps);
    printf("T/T_orbital = %.1f orbits\n", T/T_orbital);
    printf("Packet: r=%.2f×Bohr=%.0f  σ=%.2f×Bohr=%.0f\n\n",
           r_frac, r_frac*BOHR, sigma_frac, sigma_frac*BOHR);

    mkdir("data", 0755); mkdir(outdir, 0755);

    Wedge *w = wedge_alloc(Nr);
    w->dt = dt;
    wedge_init_gaussian(w, r_frac * BOHR, sigma_frac * BOHR);

    char tspath[512];
    snprintf(tspath, sizeof(tspath), "%s/timeseries.tsv", outdir);
    FILE *fp = fopen(tspath, "w");
    fprintf(fp, "step\tt\tmean_r\trms_r\tnorm\tpeak_r\tenergy\tmean_r_over_bohr\n");

    printf("step         t              <r>        rms_r      norm     peak_r     energy     <r>/a₀\n");
    printf("--------------------------------------------------------------------------------------------\n");

    for (long step = 0; step <= n_steps; step++) {
        if (step > 0) wedge_step(w);

        if (step % diag_every == 0 || step == n_steps) {
            double mr, rr, nm, pk, E;
            measure(w, &mr, &rr, &nm, &pk, &E);
            double t = step * dt;

            fprintf(fp, "%ld\t%.6e\t%.2f\t%.2f\t%.8f\t%.2f\t%.6e\t%.4f\n",
                    step, t, mr, rr, nm, pk, E, mr/BOHR);
            fflush(fp);

            if (step % (diag_every * 100) == 0 || step == 0 || step == n_steps) {
                printf("%10ld  %13.4e  %10.0f  %10.0f  %8.6f  %10.0f  %10.3e  %6.3f\n",
                       step, t, mr, rr, nm, pk, E, mr/BOHR);
                fflush(stdout);
            }
        }
    }

    fclose(fp);
    printf("\nComplete. %ld steps.\n", n_steps);

    /* Save final wavefunction */
    char wfpath[512];
    snprintf(wfpath, sizeof(wfpath), "%s/wavefunction_final.tsv", outdir);
    FILE *wf = fopen(wfpath, "w");
    fprintf(wf, "r\tpsi_re\tpsi_im\tprob\tV\n");
    for (int i = 0; i < Nr; i++) {
        double prob = (w->psi_re[i]*w->psi_re[i]+w->psi_im[i]*w->psi_im[i]) * w->r[i]*w->r[i];
        fprintf(wf, "%.2f\t%.8e\t%.8e\t%.8e\t%.8e\n",
                w->r[i], w->psi_re[i], w->psi_im[i], prob, w->V[i]);
    }
    fclose(wf);

    free(w->r); free(w->dr); free(w->V);
    free(w->psi_re); free(w->psi_im);
    free(w->u_re); free(w->u_im);
    free(w->Hu_re); free(w->Hu_im);
    free(w->rhs_re); free(w->rhs_im);
    free(w->a_coeff); free(w->b_coeff); free(w->c_coeff);
    free(w->d_re); free(w->d_im); free(w->f_re); free(w->f_im);
    free(w);
    return 0;
}
