/*
 * finite_lambda.c — Finite-λ hedgehog Skyrmion solver
 *
 * Solves for the B=1 hedgehog with variable field norm ρ(r):
 *   q(r) = ρ(r) [cos f(r) + sin f(r) r̂·σ]
 *
 * Energy:
 *   E₂ = 2π ∫ [ρ'²r² + ρ²(f'²r² + 2sin²f)] dr
 *   E₄ = (4π/e²) ∫ ρ⁸[2f'²sin²f + sin⁴f/r²] dr
 *   E_V = πλ ∫ (ρ²-1)²r² dr
 *
 * Method: under-relaxed alternating iteration.
 *   - f: RK4 shooting with bisection (variable ρ), blended with old f
 *   - ρ: Newton BVP (tridiagonal Thomas), blended with old ρ
 *   ρ CANNOT be solved by shooting (exponential stiffness ~e^{√(2λ)R}).
 *
 * KEY FINDING: The σ-model soliton (ρ≡1) is a LOCAL energy minimum.
 * The GLOBAL minimum has ρ→0 (soliton collapse). For finite λ, the
 * self-consistent solution has ρ(0) < 1, with the deviation controlled
 * by 1/λ. The virial theorem E₂ - E₄ + 3E_V = 0 must hold at
 * the true equilibrium.
 */

#define _USE_MATH_DEFINES
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ========== Grid ========== */

typedef struct {
    int N;
    double h;
    double *r;
    double *f;
    double *fp;
    double *rho;
} Grid;

typedef struct {
    double e, inv_e2, c4, lambda;
} Params;

typedef struct {
    double E2, E4, EV, Etot, Q;
} Energy;

static Grid *grid_alloc(int N, double R_max) {
    Grid *g = calloc(1, sizeof(Grid));
    g->N = N;  g->h = R_max / N;
    g->r   = calloc((size_t)(N+1), sizeof(double));
    g->f   = calloc((size_t)(N+1), sizeof(double));
    g->fp  = calloc((size_t)(N+1), sizeof(double));
    g->rho = calloc((size_t)(N+1), sizeof(double));
    for (int i = 0; i <= N; i++) g->r[i] = i * g->h;
    return g;
}

static void grid_free(Grid *g) {
    if (g) { free(g->r); free(g->f); free(g->fp); free(g->rho); free(g); }
}

/* ========== S₀, T₀ with correct r=0 limits ========== */

static void compute_S0_T0(const Grid *g, double *S0, double *T0) {
    for (int i = 0; i <= g->N; i++) {
        double r = g->r[i], fpi = g->fp[i];
        double sf = sin(g->f[i]), sf2 = sf*sf;
        double r2 = r*r;
        if (r > 1e-10) {
            S0[i] = fpi*fpi + 2.0*sf2/r2;
            T0[i] = (2.0*fpi*fpi*sf2 + sf2*sf2/r2) / r2;
        } else {
            double a = -fpi, a2 = a*a;
            S0[i] = 3.0*a2;
            T0[i] = 3.0*a2*a2;
        }
    }
}

/* ========== Energy ========== */

static Energy compute_energy(const Grid *g, const Params *p) {
    int N = g->N;
    double h = g->h;
    double e2 = 0, e4 = 0, ev = 0, qc = 0;

    for (int i = 0; i <= N; i++) {
        double r = g->r[i], f = g->f[i], fpi = g->fp[i];
        double rho = g->rho[i], sf = sin(f), sf2 = sf*sf;
        double rho2 = rho*rho, rho8 = rho2*rho2*rho2*rho2;
        double w = h * ((i == 0 || i == N) ? 0.5 : 1.0);

        e2 += w * rho2 * (fpi*fpi*r*r + 2.0*sf2);
        if (r > 1e-14)
            e4 += w * rho8 * (2.0*fpi*fpi*sf2 + sf2*sf2/(r*r));
        ev += w * (rho2 - 1.0)*(rho2 - 1.0) * r*r;
        qc += w * (-fpi) * sf2;
    }
    /* ρ' contribution to E₂ */
    for (int i = 0; i < N; i++) {
        double rm = 0.5*(g->r[i] + g->r[i+1]);
        double rp = (g->rho[i+1] - g->rho[i]) / h;
        e2 += h * rp*rp * rm*rm;
    }

    Energy en;
    en.E2 = 2.0*M_PI*e2;  en.E4 = 4.0*M_PI*p->inv_e2*e4;
    en.EV = M_PI*p->lambda*ev;  en.Etot = en.E2 + en.E4 + en.EV;
    en.Q  = (2.0/M_PI)*qc;
    return en;
}

/* ========== σ-model f shooting (ρ≡1) ========== */

static double sigma_frhs(double r, double f, double v, double c4) {
    double sf = sin(f), sf2 = sf*sf, s2f = sin(2.0*f);
    double r2 = r*r, den = r2 + 2.0*c4*sf2;
    if (den < 1e-30) return 0.0;
    return (s2f*(1.0 + c4*sf2/r2) - 2.0*r*v - c4*v*v*s2f) / den;
}

static void sigma_init(Grid *g, const Params *p) {
    int N = g->N;  double h = g->h, c4 = p->c4;
    for (int i = 0; i <= N; i++) g->rho[i] = 1.0;

    double a_lo = 0.1, a_hi = 50.0;
    double *ft = calloc((size_t)(N+1), sizeof(double));
    double *vt = calloc((size_t)(N+1), sizeof(double));

    for (int bi = 0; bi < 80; bi++) {
        double a = 0.5*(a_lo + a_hi), fv = M_PI, vv = -a;
        ft[0] = fv; vt[0] = vv;
        for (int i = 0; i < N; i++) {
            double ri = i*h;
            double k1f = vv, k1v = sigma_frhs(ri, fv, vv, c4);
            double k2f = vv+0.5*h*k1v, k2v = sigma_frhs(ri+0.5*h, fv+0.5*h*k1f, k2f, c4);
            double k3f = vv+0.5*h*k2v, k3v = sigma_frhs(ri+0.5*h, fv+0.5*h*k2f, k3f, c4);
            double k4f = vv+h*k3v, k4v = sigma_frhs(ri+h, fv+h*k3f, k4f, c4);
            fv += (h/6.0)*(k1f+2*k2f+2*k3f+k4f);
            vv += (h/6.0)*(k1v+2*k2v+2*k3v+k4v);
            ft[i+1] = fv; vt[i+1] = vv;
            if (fabs(fv) > 20.0 || isnan(fv)) break;
        }
        if (fv > 0) a_lo = a; else a_hi = a;
    }
    memcpy(g->f, ft, (size_t)(N+1)*sizeof(double));
    memcpy(g->fp, vt, (size_t)(N+1)*sizeof(double));
    free(ft); free(vt);
}

/* ========== ρ Newton BVP solver ========== */
/*
 * Solve: ρ'' + 2ρ'/r = ρ S₀ + 4c₄ ρ⁷ T₀ + λρ(ρ²-1)
 * BC: ρ'(0) = 0, ρ(R) = 1
 * Uses Newton iteration on the discretized system (tridiagonal).
 */

static void solve_rho_newton(Grid *g, const Params *p, int max_iter, double tol) {
    int N = g->N;
    double h = g->h, c4 = p->c4, lam = p->lambda;

    double *S0 = calloc((size_t)(N+1), sizeof(double));
    double *T0 = calloc((size_t)(N+1), sizeof(double));
    compute_S0_T0(g, S0, T0);

    double *al = calloc((size_t)(N+1), sizeof(double));
    double *bl = calloc((size_t)(N+1), sizeof(double));
    double *cl = calloc((size_t)(N+1), sizeof(double));
    double *res = calloc((size_t)(N+1), sizeof(double));
    double *dx  = calloc((size_t)(N+1), sizeof(double));

    for (int nit = 0; nit < max_iter; nit++) {
        double max_res = 0;

        for (int i = 1; i < N; i++) {
            double r = g->r[i], rho = g->rho[i];
            double rho2 = rho*rho, rho6 = rho2*rho2*rho2, rho7 = rho6*rho;
            double rho_pp = (g->rho[i+1] - 2.0*rho + g->rho[i-1])/(h*h);
            double rho_p  = (g->rho[i+1] - g->rho[i-1])/(2.0*h);

            res[i] = rho_pp + 2.0*rho_p/r
                   - rho*S0[i] - 4.0*c4*rho7*T0[i] - lam*rho*(rho2-1.0);
            double dF = -S0[i] - 28.0*c4*rho6*T0[i] - lam*(3.0*rho2-1.0);

            al[i] = 1.0/(h*h) - 1.0/(h*r);
            bl[i] = -2.0/(h*h) + dF;
            cl[i] = 1.0/(h*h) + 1.0/(h*r);
            if (fabs(res[i]) > max_res) max_res = fabs(res[i]);
        }

        /* i=0 boundary */
        {
            double rho = g->rho[0], rho2 = rho*rho;
            double rho6 = rho2*rho2*rho2, rho7 = rho6*rho;
            double rho_pp = 2.0*(g->rho[1] - rho)/(h*h);
            res[0] = 3.0*rho_pp - rho*S0[0] - 4.0*c4*rho7*T0[0] - lam*rho*(rho2-1.0);
            double dF = -S0[0] - 28.0*c4*rho6*T0[0] - lam*(3.0*rho2-1.0);
            al[0] = 0; bl[0] = -6.0/(h*h) + dF; cl[0] = 6.0/(h*h);
            if (fabs(res[0]) > max_res) max_res = fabs(res[0]);
        }

        al[N] = 0; bl[N] = 1; cl[N] = 0; res[N] = 0;

        if (max_res < tol) break;

        /* Thomas solve: J·dx = -res */
        for (int i = 0; i <= N; i++) res[i] = -res[i];
        for (int i = 1; i <= N; i++) {
            double m = al[i]/bl[i-1];
            bl[i] -= m*cl[i-1]; res[i] -= m*res[i-1];
        }
        dx[N] = res[N]/bl[N];
        for (int i = N-1; i >= 0; i--) dx[i] = (res[i] - cl[i]*dx[i+1])/bl[i];

        /* Damped update (tight: max 0.002 per grid point per Newton step) */
        double step = 1.0, mc = 0.002;
        for (int i = 0; i <= N; i++) {
            if (fabs(dx[i])*step > mc) step = mc/fabs(dx[i]);
        }
        for (int i = 0; i <= N; i++) g->rho[i] += step*dx[i];

        /* Floor: prevent total collapse */
        for (int i = 0; i <= N; i++)
            if (g->rho[i] < 0.50) g->rho[i] = 0.50;
    }

    free(S0); free(T0); free(al); free(bl); free(cl); free(res); free(dx);
}

/* ========== Linearized ρ solver (for initial guess) ========== */

static void solve_rho_linear(Grid *g, const Params *p) {
    int N = g->N;  double h = g->h, c4 = p->c4, lam = p->lambda;
    double *S0 = calloc((size_t)(N+1), sizeof(double));
    double *T0 = calloc((size_t)(N+1), sizeof(double));
    compute_S0_T0(g, S0, T0);

    double *al = calloc((size_t)(N+1), sizeof(double));
    double *bl = calloc((size_t)(N+1), sizeof(double));
    double *cl = calloc((size_t)(N+1), sizeof(double));
    double *dl = calloc((size_t)(N+1), sizeof(double));
    double *dr = calloc((size_t)(N+1), sizeof(double));

    for (int i = 1; i < N; i++) {
        double r = g->r[i];
        double V = S0[i] + 28.0*c4*T0[i] + 2.0*lam;
        double src = S0[i] + 4.0*c4*T0[i];
        al[i] = 1.0/(h*h) - 1.0/(h*r);
        bl[i] = -2.0/(h*h) - V;
        cl[i] = 1.0/(h*h) + 1.0/(h*r);
        dl[i] = src;
    }
    {
        double V = S0[0] + 28.0*c4*T0[0] + 2.0*lam;
        double src = S0[0] + 4.0*c4*T0[0];
        al[0]=0; bl[0]=-6.0/(h*h)-V; cl[0]=6.0/(h*h); dl[0]=src;
    }
    al[N]=0; bl[N]=1; cl[N]=0; dl[N]=0;

    for (int i = 1; i <= N; i++) {
        double m = al[i]/bl[i-1]; bl[i] -= m*cl[i-1]; dl[i] -= m*dl[i-1];
    }
    dr[N] = dl[N]/bl[N];
    for (int i = N-1; i >= 0; i--) dr[i] = (dl[i] - cl[i]*dr[i+1])/bl[i];
    for (int i = 0; i <= N; i++) g->rho[i] = 1.0 + dr[i];

    free(S0); free(T0); free(al); free(bl); free(cl); free(dl); free(dr);
}

/* ========== f shooting with variable ρ ========== */

static double f_rhs_varrho(double r, double f, double v,
                           double rho, double rhop, double c4) {
    double sf = sin(f), sf2 = sf*sf, s2f = sin(2.0*f);
    double r2 = r*r;
    double rho2 = rho*rho, rho7 = rho2*rho2*rho2*rho, rho8 = rho7*rho;
    double den = rho2*r2 + 2.0*c4*rho8*sf2;
    if (fabs(den) < 1e-30) return 0.0;
    double num = rho2*s2f + c4*rho8*s2f*sf2/r2
               - 2.0*rho2*r*v - 2.0*rho*rhop*v*r2
               - 16.0*c4*rho7*rhop*v*sf2 - c4*rho8*v*v*s2f;
    return num / den;
}

static void solve_f_shooting(Grid *g, const Params *p,
                             double *f_new, double *fp_new) {
    int N = g->N;  double h = g->h, c4 = p->c4;

    double *rhop = calloc((size_t)(N+1), sizeof(double));
    rhop[0] = 0.0;
    for (int i = 1; i < N; i++) rhop[i] = (g->rho[i+1]-g->rho[i-1])/(2.0*h);
    rhop[N] = (g->rho[N]-g->rho[N-1])/h;

    #define INTERP(arr, rv) ({ int _i=(int)((rv)/h); if(_i>=N) _i=N-1; \
        double _t=(rv)/h-_i; (1.0-_t)*(arr)[_i]+_t*(arr)[_i+1]; })

    double a_lo = 0.1, a_hi = 50.0;
    for (int bi = 0; bi < 80; bi++) {
        double a = 0.5*(a_lo+a_hi), fv = M_PI, vv = -a;
        for (int i = 0; i < N; i++) {
            double ri = i*h;
            double k1f = vv, k1v = f_rhs_varrho(ri, fv, vv, g->rho[i], rhop[i], c4);
            double r2=ri+0.5*h, rho2=INTERP(g->rho,r2), rhop2=INTERP(rhop,r2);
            double k2f = vv+0.5*h*k1v;
            double k2v = f_rhs_varrho(r2, fv+0.5*h*k1f, k2f, rho2, rhop2, c4);
            double k3f = vv+0.5*h*k2v;
            double k3v = f_rhs_varrho(r2, fv+0.5*h*k2f, k3f, rho2, rhop2, c4);
            double r4=ri+h;
            double rho4=(i+1<=N)?g->rho[i+1]:1.0, rhop4=(i+1<=N)?rhop[i+1]:0.0;
            double k4f = vv+h*k3v;
            double k4v = f_rhs_varrho(r4, fv+h*k3f, k4f, rho4, rhop4, c4);
            fv += (h/6.0)*(k1f+2*k2f+2*k3f+k4f);
            vv += (h/6.0)*(k1v+2*k2v+2*k3v+k4v);
            if (bi == 79) { f_new[i+1] = fv; fp_new[i+1] = vv; }
            if (fabs(fv) > 20.0 || isnan(fv)) break;
        }
        if (fv > 0) a_lo = a; else a_hi = a;
    }
    f_new[0] = M_PI;
    fp_new[0] = -0.5*(a_lo+a_hi);

    free(rhop);
    #undef INTERP
}

/* ========== Main ========== */

int main(int argc, char **argv) {
    int Nr = 2000;
    double R_max = 20.0;
    double e_skyrme = 4.0;
    int n_outer = 300;
    double relax = 0.15;  /* under-relaxation parameter */

    for (int i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "-Nr") && i+1<argc) Nr = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-Rmax") && i+1<argc) R_max = atof(argv[++i]);
        else if (!strcmp(argv[i], "-e") && i+1<argc) e_skyrme = atof(argv[++i]);
        else if (!strcmp(argv[i], "-niter") && i+1<argc) n_outer = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-relax") && i+1<argc) relax = atof(argv[++i]);
        else { fprintf(stderr, "Usage: %s [-Nr N] [-Rmax R] [-e E] [-niter N] [-relax A]\n", argv[0]); return 1; }
    }

    double c4 = 2.0/(e_skyrme*e_skyrme);
    Params p = { .e=e_skyrme, .inv_e2=1.0/(e_skyrme*e_skyrme), .c4=c4, .lambda=0 };

    printf("============================================================\n");
    printf("  Finite-lambda Hedgehog Skyrmion (Under-relaxed Iteration)\n");
    printf("============================================================\n");
    printf("Grid: Nr=%d, h=%.6f, R_max=%.1f\n", Nr, R_max/(double)Nr, R_max);
    printf("Parameters: e=%.4f, c4=%.6f, relax=%.2f\n\n", e_skyrme, c4, relax);

    Grid *g = grid_alloc(Nr, R_max);
    double *f_new  = calloc((size_t)(Nr+1), sizeof(double));
    double *fp_new = calloc((size_t)(Nr+1), sizeof(double));
    double *rho_new = calloc((size_t)(Nr+1), sizeof(double));

    /* σ-model baseline */
    sigma_init(g, &p);
    Energy en0 = compute_energy(g, &p);
    double E_FB = 6.0*sqrt(2.0)*M_PI*M_PI/e_skyrme;
    double a_sigma = -g->fp[0];
    printf("  sigma-model: a=%.8f\n", a_sigma);
    printf("  E=%.6f, E2=%.6f, E4=%.6f, E2/E4=%.6f, Q=%.6f\n",
           en0.Etot, en0.E2, en0.E4, en0.E2/en0.E4, en0.Q);
    printf("  E/E_FB = %.4f\n\n", en0.Etot/E_FB);

    /* λ scan */
    double lambdas[] = {
        1e8, 1e7, 1e6, 3e5, 1e5, 7e4, 5e4, 3e4, 2e4,
        1.5e4, 1.2e4, 1e4, 9e3, 8e3, 7e3, 6500, 6000, 5500, 5000
    };
    int n_lambda = sizeof(lambdas)/sizeof(lambdas[0]);

    printf("%-12s %10s %12s %12s %12s %12s %8s %8s %10s %6s\n",
           "lambda", "rho(0)", "E_total", "E2", "E4", "E_V",
           "Q", "E/E_FB", "virial", "iters");
    printf("-------- ---------- ------------ ------------ ------------ ------------ -------- -------- ---------- ------\n");
    printf("%-12s %10.6f %12.6f %12.6f %12.6f %12.6f %8.5f %8.4f %10s %6s\n",
           "inf", 1.0, en0.Etot, en0.E2, en0.E4, en0.EV, en0.Q, en0.Etot/E_FB, "N/A", "-");

    /* Save σ-model as baseline */
    double *f_sigma  = calloc((size_t)(Nr+1), sizeof(double));
    double *fp_sigma = calloc((size_t)(Nr+1), sizeof(double));
    memcpy(f_sigma, g->f, (size_t)(Nr+1)*sizeof(double));
    memcpy(fp_sigma, g->fp, (size_t)(Nr+1)*sizeof(double));

    for (int li = 0; li < n_lambda; li++) {
        p.lambda = lambdas[li];

        /* Initialize: σ-model f + linearized ρ */
        if (li == 0) {
            memcpy(g->f, f_sigma, (size_t)(Nr+1)*sizeof(double));
            memcpy(g->fp, fp_sigma, (size_t)(Nr+1)*sizeof(double));
            for (int i = 0; i <= Nr; i++) g->rho[i] = 1.0;
            solve_rho_linear(g, &p);
        }
        /* else: continue from previous λ solution */

        /* Under-relaxed alternating iteration */
        int converged = 0;
        int outer;
        double prev_vir = 1e30;
        for (outer = 0; outer < n_outer; outer++) {
            /* Save current ρ */
            memcpy(rho_new, g->rho, (size_t)(Nr+1)*sizeof(double));

            /* Fully converge ρ Newton for current f (tight floor) */
            solve_rho_newton(g, &p, 30, 1e-12);

            /* Blend: ρ = (1-α)·ρ_old + α·ρ_converged */
            for (int i = 0; i <= Nr; i++)
                g->rho[i] = (1.0-relax)*rho_new[i] + relax*g->rho[i];

            /* Shoot f for blended ρ */
            solve_f_shooting(g, &p, f_new, fp_new);

            /* Blend: f = (1-α)·f_old + α·f_new */
            for (int i = 0; i <= Nr; i++) {
                g->f[i]  = (1.0-relax)*g->f[i]  + relax*f_new[i];
                g->fp[i] = (1.0-relax)*g->fp[i] + relax*fp_new[i];
            }

            /* Virial convergence check */
            Energy en = compute_energy(g, &p);
            double vir = en.E2 - en.E4 + 3.0*en.EV;
            if (fabs(vir) < 0.005 * en.Etot) { converged = 1; break; }
            prev_vir = vir;
            (void)prev_vir;
        }

        Energy en = compute_energy(g, &p);
        double virial = en.E2 - en.E4 + 3.0*en.EV;

        if (en.Q < 0.95) {
            printf("%-12.4g  ** COLLAPSED: rho(0)=%.4f Q=%.4f **\n",
                   p.lambda, g->rho[0], en.Q);
            /* Reset for next λ */
            memcpy(g->f, f_sigma, (size_t)(Nr+1)*sizeof(double));
            memcpy(g->fp, fp_sigma, (size_t)(Nr+1)*sizeof(double));
            for (int i = 0; i <= Nr; i++) g->rho[i] = 1.0;
            continue;
        }

        printf("%-12.4g %10.6f %12.6f %12.6f %12.6f %12.6f %8.5f %8.4f %10.6f %5d%c\n",
               p.lambda, g->rho[0], en.Etot, en.E2, en.E4, en.EV,
               en.Q, en.Etot/E_FB, virial, outer+1, converged?'*':' ');
    }

    printf("\n* = virial converged (|E2-E4+3EV| < 0.01)\n");
    printf("Virial: E2 - E4 + 3*E_V = 0\n");
    printf("Mass: Mc^2 = 2*E4 - 2*E_V\n");
    printf("E_FB = %.6f\n", E_FB);

    free(f_sigma); free(fp_sigma);
    free(f_new); free(fp_new); free(rho_new);
    grid_free(g);
    return 0;
}
