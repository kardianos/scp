/*  T7: Two-Braid Interaction at the Bimodal Sweet Spot
 *
 *  Two braids separated by D=30 along x-axis.
 *  Three configs: same twist, opposite twist, single (control).
 *  Tracks separation D(t) via center-of-mass of |phi|^2 in each half.
 *  Uses BIMODAL params with mass^2=2.25 (mass2_override=-1).
 *
 *  Build: cd T7_twobraid && gcc -O3 -fopenmp -o t7 src/t7.c -lm
 */

#include "../../src/braid_core.h"

/* ================================================================
   Two-braid initialization
   Places braid 1 at (x1, 0, 0) and braid 2 at (x2, 0, 0).
   If flip_second=1, negate phi[1] for the second braid (opposite twist).
   ================================================================ */

static void init_two_braids(Grid *g, const double *phys,
                            double cx1, double cx2, int flip_second)
{
    double A[3]     = {phys[0], phys[1], phys[2]};
    double delta[3] = {0.0, phys[3], phys[4]};
    double R_tube   = phys[5];
    double ellip    = phys[6];
    double ell_ang  = phys[7];
    double k_fac    = phys[8];
    double A_bg     = phys[9];
    double R_disp   = phys[10];
    double ell_rot  = phys[11];
    double m_init   = phys[14];

    int N = g->N, NN = N * N;
    double dx = g->dx, L = g->L;
    double k = k_fac * PI / L;
    double omega = sqrt(k * k + m_init * m_init);
    double inv_2R2 = 1.0 / (2.0 * R_tube * R_tube);
    double sx = 1.0 + ellip, sy = 1.0 - ellip;

    /* Strand offsets in (y,z) plane for 3-strand braid */
    double scx[3], scy[3], ea[3];
    for (int a = 0; a < 3; a++) {
        double ang = 2.0 * PI * a / 3.0;
        scx[a] = R_disp * cos(ang);
        scy[a] = R_disp * sin(ang);
        ea[a] = (ell_rot > 0.5) ? ell_ang + 2.0*PI*a/3.0 : ell_ang;
    }

    grid_zero(g);

    for (int i = 0; i < N; i++) {
        double x = -L + i * dx;
        for (int j = 0; j < N; j++) {
            double y = -L + j * dx;
            for (int kk = 0; kk < N; kk++) {
                double z = -L + kk * dx;
                int idx = i * NN + j * N + kk;

                /* Braid 1 at (cx1, 0, 0) */
                for (int a = 0; a < NFIELDS; a++) {
                    double xc = (x - cx1) - scx[a];
                    double yc = y - scy[a];
                    double ca = cos(ea[a]), sa = sin(ea[a]);
                    double xr = xc*ca + yc*sa;
                    double yr = -xc*sa + yc*ca;
                    double r2e = xr*xr/(sx*sx) + yr*yr/(sy*sy);
                    double env = exp(-r2e * inv_2R2);
                    double ph = k * z + delta[a];
                    double amp = A[a] * env + A_bg;
                    g->phi[a][idx] += amp * cos(ph);
                    g->vel[a][idx] += omega * amp * sin(ph);
                }

                /* Braid 2 at (cx2, 0, 0) */
                for (int a = 0; a < NFIELDS; a++) {
                    double xc = (x - cx2) - scx[a];
                    double yc = y - scy[a];
                    double ca = cos(ea[a]), sa = sin(ea[a]);
                    double xr = xc*ca + yc*sa;
                    double yr = -xc*sa + yc*ca;
                    double r2e = xr*xr/(sx*sx) + yr*yr/(sy*sy);
                    double env = exp(-r2e * inv_2R2);
                    double ph = k * z + delta[a];
                    double amp = A[a] * env + A_bg;
                    double sign = (flip_second && a == 1) ? -1.0 : 1.0;
                    g->phi[a][idx] += sign * amp * cos(ph);
                    g->vel[a][idx] += sign * omega * amp * sin(ph);
                }
            }
        }
    }
}

/* ================================================================
   Center-of-mass of |phi|^2 in x-half (left: x<0, right: x>=0)
   Returns x_cm.  Also returns fc (fraction of |phi|^2 in core R<8
   around the half's center).
   ================================================================ */

static double com_half(Grid *g, int right, double *fc_out) {
    int N = g->N, NN = N * N;
    double dx = g->dx, L = g->L;
    double sum_rho = 0, sum_xrho = 0;
    double sum_core = 0;
    double R_core = 8.0;
    double Rc2 = R_core * R_core;

    for (int i = 1; i < N-1; i++) {
        double x = -L + i * dx;
        if (right && x < 0) continue;
        if (!right && x >= 0) continue;
        for (int j = 1; j < N-1; j++) {
            double y = -L + j * dx;
            for (int kk = 0; kk < N; kk++) {
                int idx = i * NN + j * N + kk;
                double rho = 0;
                for (int a = 0; a < NFIELDS; a++)
                    rho += g->phi[a][idx] * g->phi[a][idx];
                sum_rho += rho;
                sum_xrho += rho * x;
                /* core: r_perp < R_core from half center (±15) */
                double cx_half = right ? 15.0 : -15.0;
                double dr2 = (x - cx_half)*(x - cx_half) + y*y;
                if (dr2 < Rc2) sum_core += rho;
            }
        }
    }
    if (fc_out) *fc_out = sum_core / (sum_rho + 1e-30);
    return sum_xrho / (sum_rho + 1e-30);
}

/* ================================================================
   Total energy (same as in braid_core diagnostics)
   ================================================================ */

static double total_energy(Grid *g, const double *phys) {
    int N = g->N, NN = N*N;
    double dx = g->dx, L = g->L, dV = dx*dx*dx;
    double mass2 = phys[14]*phys[14];
    double mu = phys[12], kappa = phys[13], lpw = phys[15];

    double E = 0;
    #pragma omp parallel for reduction(+:E) schedule(static)
    for (int i = 1; i < N-1; i++) {
        double x = -L + i*dx;
        for (int j = 1; j < N-1; j++) {
            double y = -L + j*dx;
            for (int k = 0; k < N; k++) {
                int idx = i*NN + j*N + k;
                double p0 = g->phi[0][idx], p1 = g->phi[1][idx], p2 = g->phi[2][idx];
                double rho = p0*p0 + p1*p1 + p2*p2;

                double ek = 0, eg = 0;
                for (int a = 0; a < NFIELDS; a++) {
                    ek += 0.5 * g->vel[a][idx] * g->vel[a][idx];
                    int kp = (k+1)%N, km = (k-1+N)%N;
                    double gx = (g->phi[a][idx+NN]-g->phi[a][idx-NN])/(2*dx);
                    double gy = (g->phi[a][idx+N] -g->phi[a][idx-N]) /(2*dx);
                    double gz = (g->phi[a][i*NN+j*N+kp]-g->phi[a][i*NN+j*N+km])/(2*dx);
                    eg += 0.5*(gx*gx + gy*gy + gz*gz);
                }
                double em = 0.5 * mass2 * rho;
                double P = p0*p1*p2;
                double ep = (mu/2.0)*P*P/(1.0+kappa*P*P);
                double epw = lpw*(p0*p1 + p1*p2 + p2*p0);
                E += (ek + eg + em + ep + epw) * dV;
            }
        }
    }
    return E;
}

/* ================================================================ */

#define MAX_SNAP 100

static void run_config(const char *label, double cx1, double cx2,
                       int flip_second, int single_braid,
                       double T_max, double T_diag,
                       double *t_out, double *D_out, double *E_out,
                       double *fcL_out, double *fcR_out, int *n_snap_out)
{
    printf("\n============================================================\n");
    printf("  Config: %s\n", label);
    printf("============================================================\n\n");

    int N = 192;
    double L = 40.0;
    Grid *g = grid_alloc(N, L);
    printf("  Grid: N=%d, L=%.1f, dx=%.4f, dt=%.6f\n", N, L, g->dx, g->dt);

    if (single_braid) {
        init_braid(g, BIMODAL, -1);
        printf("  Single braid at center (control)\n");
    } else {
        init_two_braids(g, BIMODAL, cx1, cx2, flip_second);
        printf("  Braid 1 at x=%.1f, Braid 2 at x=%.1f, flip=%d\n", cx1, cx2, flip_second);
    }

    compute_forces(g, BIMODAL, -1);

    int steps_per_diag = (int)(T_diag / g->dt + 0.5);
    int total_steps = (int)(T_max / g->dt + 0.5);
    int n_diag = total_steps / steps_per_diag;
    if (n_diag > MAX_SNAP - 1) n_diag = MAX_SNAP - 1;

    printf("  Total steps: %d, diag every %d steps\n\n", total_steps, steps_per_diag);

    /* Header */
    printf("  %8s %10s %12s %8s %8s\n", "T", "D", "E", "fc_L", "fc_R");
    printf("  -------- ---------- ------------ -------- --------\n");

    /* Initial snapshot */
    int ns = 0;
    {
        double fcL, fcR;
        double xL = com_half(g, 0, &fcL);
        double xR = com_half(g, 1, &fcR);
        double D = xR - xL;
        double E = total_energy(g, BIMODAL);
        t_out[ns] = 0; D_out[ns] = D; E_out[ns] = E;
        fcL_out[ns] = fcL; fcR_out[ns] = fcR;
        ns++;
        printf("  %8.1f %10.4f %12.2f %8.4f %8.4f\n", 0.0, D, E, fcL, fcR);
    }

    int step = 0;
    int blown = 0;
    double wall0 = omp_get_wtime();

    for (int diag = 1; diag <= n_diag && !blown; diag++) {
        int target = diag * steps_per_diag;
        if (target > total_steps) target = total_steps;

        for (; step < target; step++) {
            verlet_full_step(g, BIMODAL, -1);
            apply_damping_xy(g);
        }

        if (check_blowup(g)) {
            blown = 1;
            printf("  BLOWUP at step %d\n", step);
            break;
        }

        double T = step * g->dt;
        double fcL, fcR;
        double xL = com_half(g, 0, &fcL);
        double xR = com_half(g, 1, &fcR);
        double D = xR - xL;
        double E = total_energy(g, BIMODAL);

        if (ns < MAX_SNAP) {
            t_out[ns] = T; D_out[ns] = D; E_out[ns] = E;
            fcL_out[ns] = fcL; fcR_out[ns] = fcR;
            ns++;
        }

        printf("  %8.1f %10.4f %12.2f %8.4f %8.4f\n", T, D, E, fcL, fcR);
        fflush(stdout);

        if (diag == 1 || diag % 5 == 0) {
            double wall = omp_get_wtime() - wall0;
            double frac = (double)step / total_steps;
            printf("    [wall %.1fs, %.1f%% done, ETA %.0fs]\n",
                   wall, 100*frac, wall*(1-frac)/(frac+1e-30));
        }
    }

    double wall_total = omp_get_wtime() - wall0;
    printf("\n  Finished %d steps in %.1f s (%.4f s/step)\n", step, wall_total,
           wall_total/(step+1));

    *n_snap_out = ns;
    grid_free(g);
}

/* ================================================================ */

int main(void) {
    printf("=== T7: Two-Braid Interaction ===\n\n");

    bimodal_init_params();
    printf("BIMODAL params:\n");
    for (int d = 0; d < NDIM; d++)
        printf("  %8s = %8.4f\n", PNAME[d], BIMODAL[d]);
    printf("\nmass = %.4f => mass^2 = %.4f\n", BIMODAL[14], BIMODAL[14]*BIMODAL[14]);
    printf("mass2_override = -1 (use mass^2 from params)\n\n");

    double T_max = 500.0, T_diag = 25.0;

    /* Arrays for each config */
    double t1[MAX_SNAP], D1[MAX_SNAP], E1[MAX_SNAP], fcL1[MAX_SNAP], fcR1[MAX_SNAP];
    double t2[MAX_SNAP], D2[MAX_SNAP], E2[MAX_SNAP], fcL2[MAX_SNAP], fcR2[MAX_SNAP];
    double t3[MAX_SNAP], D3[MAX_SNAP], E3[MAX_SNAP], fcL3[MAX_SNAP], fcR3[MAX_SNAP];
    int ns1, ns2, ns3;

    /* Config A: Same twist */
    run_config("Same twist (both W=-1)", -15.0, +15.0, 0, 0,
               T_max, T_diag, t1, D1, E1, fcL1, fcR1, &ns1);

    /* Config B: Opposite twist */
    run_config("Opposite twist (W=-1, W=+1)", -15.0, +15.0, 1, 0,
               T_max, T_diag, t2, D2, E2, fcL2, fcR2, &ns2);

    /* Config C: Single braid control */
    run_config("Single braid (control)", 0.0, 0.0, 0, 1,
               T_max, T_diag, t3, D3, E3, fcL3, fcR3, &ns3);

    /* ================================================================
       Summary
       ================================================================ */
    printf("\n============================================================\n");
    printf("  SUMMARY\n");
    printf("============================================================\n\n");

    if (ns1 >= 2) {
        double dD = D1[ns1-1] - D1[0];
        printf("Same twist:     D(0)=%8.4f  D(end)=%8.4f  deltaD=%+.4f  %s\n",
               D1[0], D1[ns1-1], dD, dD < 0 ? "ATTRACTION" : (dD > 0 ? "REPULSION" : "NEUTRAL"));
        printf("                E(0)=%10.2f  E(end)=%10.2f  dE/E=%.4f\n",
               E1[0], E1[ns1-1], (E1[ns1-1]-E1[0])/(E1[0]+1e-30));
    }
    if (ns2 >= 2) {
        double dD = D2[ns2-1] - D2[0];
        printf("Opposite twist: D(0)=%8.4f  D(end)=%8.4f  deltaD=%+.4f  %s\n",
               D2[0], D2[ns2-1], dD, dD < 0 ? "ATTRACTION" : (dD > 0 ? "REPULSION" : "NEUTRAL"));
        printf("                E(0)=%10.2f  E(end)=%10.2f  dE/E=%.4f\n",
               E2[0], E2[ns2-1], (E2[ns2-1]-E2[0])/(E2[0]+1e-30));
    }
    if (ns3 >= 2) {
        printf("Control:        E(0)=%10.2f  E(end)=%10.2f  dE/E=%.4f\n",
               E3[0], E3[ns3-1], (E3[ns3-1]-E3[0])/(E3[0]+1e-30));
    }

    /* Interpretation */
    if (ns1 >= 2 && ns2 >= 2) {
        double dD1 = D1[ns1-1] - D1[0];
        double dD2 = D2[ns2-1] - D2[0];
        printf("\n--- Interpretation ---\n");
        if (dD1 < -0.1 && dD2 < -0.1 && fabs(dD1 - dD2) < 0.5 * fabs(dD1)) {
            printf("Both attract similarly => UNIVERSAL (gravity-like)\n");
        } else if (dD1 < -0.1 && dD2 > 0.1) {
            printf("Same attract, opposite repel => CHARGE-DEPENDENT (EM-like)\n");
        } else if (dD1 > 0.1 && dD2 < -0.1) {
            printf("Same repel, opposite attract => CHARGE-DEPENDENT (EM-like, reversed)\n");
        } else if (fabs(dD1) < 0.1 && fabs(dD2) < 0.1) {
            printf("No significant motion => NEGLIGIBLE interaction at D=30\n");
        } else {
            printf("Mixed result: dD_same=%.4f, dD_opp=%.4f\n", dD1, dD2);
        }

        /* Acceleration estimate (quadratic fit: D = D0 + v0*t + 0.5*a*t^2) */
        if (ns1 >= 5) {
            double sum_t2 = 0, sum_t4 = 0, sum_Dt2 = 0;
            for (int i = 0; i < ns1; i++) {
                double ti = t1[i], dDi = D1[i] - D1[0];
                sum_t2 += ti*ti; sum_t4 += ti*ti*ti*ti;
                sum_Dt2 += dDi*ti*ti;
            }
            double a_fit = 2.0 * sum_Dt2 / (sum_t4 + 1e-30);
            printf("\nQuadratic fit (same twist): a = %.6e (D units/T^2)\n", a_fit);
        }
        if (ns2 >= 5) {
            double sum_t2 = 0, sum_t4 = 0, sum_Dt2 = 0;
            for (int i = 0; i < ns2; i++) {
                double ti = t2[i], dDi = D2[i] - D2[0];
                sum_t2 += ti*ti; sum_t4 += ti*ti*ti*ti;
                sum_Dt2 += dDi*ti*ti;
            }
            double a_fit = 2.0 * sum_Dt2 / (sum_t4 + 1e-30);
            printf("Quadratic fit (opp twist):  a = %.6e (D units/T^2)\n", a_fit);
        }
    }

    /* Write TSV */
    FILE *fp = fopen("data/t7_results.tsv", "w");
    if (fp) {
        fprintf(fp, "config\tT\tD\tE\tfc_left\tfc_right\n");
        for (int i = 0; i < ns1; i++)
            fprintf(fp, "same\t%.2f\t%.6f\t%.6f\t%.6f\t%.6f\n",
                    t1[i], D1[i], E1[i], fcL1[i], fcR1[i]);
        for (int i = 0; i < ns2; i++)
            fprintf(fp, "opposite\t%.2f\t%.6f\t%.6f\t%.6f\t%.6f\n",
                    t2[i], D2[i], E2[i], fcL2[i], fcR2[i]);
        for (int i = 0; i < ns3; i++)
            fprintf(fp, "control\t%.2f\t%.6f\t%.6f\t%.6f\t%.6f\n",
                    t3[i], D3[i], E3[i], fcL3[i], fcR3[i]);
        fclose(fp);
        printf("\nWrote data/t7_results.tsv\n");
    }

    printf("\n=== T7 Complete ===\n");
    return 0;
}
