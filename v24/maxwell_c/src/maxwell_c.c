/*
 * maxwell_c.c — 2D three-field oscillon with topological charge measurement
 *
 * V24-MC: Topological Current Coupling
 *
 * Lagrangian (2D):
 *   L = sum_a [ (1/2)(dt phi_a)^2 - (1/2)|grad phi_a|^2 - (m^2/2)phi_a^2 ]
 *     - (mu/2) P^2 / (1 + kappa P^2)
 *   P = phi_1 phi_2 phi_3
 *
 * Topological charge density (2D, three fields):
 *   rho_top = eps_{ij} eps_{abc} phi_c (d_i phi_a)(d_j phi_b)
 *           = eps_{abc} phi_c [ (d_x phi_a)(d_y phi_b) - (d_y phi_a)(d_x phi_b) ]
 *
 * This is NOT the baby-Skyrmion charge (which uses 2 fields mapping to S^1).
 * For 3 fields in 2D, this density measures whether the field configuration
 * has a nontrivial winding.
 *
 * Initialize: spherically symmetric Gaussian oscillon (phi_1=phi_2=phi_3).
 * Key question: does rho_top = 0 for this symmetric configuration?
 *
 * Compile: gcc -O3 -Wall -o maxwell_c src/maxwell_c.c -lm
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Parameters */
static double mu     = -20.0;
static double kappa  = 20.0;
static double mass   = 1.0;
static double A_init = 0.8;
static double sigma  = 3.0;
static int    Nx     = 256;
static int    Ny     = 256;
static double L      = 20.0;   /* half-extent: domain is [-L, L] */
static double tfinal = 2000.0;
static char   outdir[512] = "data";

static void parse_args(int argc, char **argv)
{
    for (int i = 1; i < argc - 1; i += 2) {
        if      (!strcmp(argv[i], "-mu"))     mu     = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-kappa"))  kappa  = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-mass"))   mass   = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-A"))      A_init = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-sigma"))  sigma  = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-Nx"))     Nx     = atoi(argv[i+1]);
        else if (!strcmp(argv[i], "-Ny"))     Ny     = atoi(argv[i+1]);
        else if (!strcmp(argv[i], "-L"))      L      = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-tfinal")) tfinal = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-o"))      strncpy(outdir, argv[i+1], sizeof(outdir)-1);
        else { fprintf(stderr, "Unknown flag: %s\n", argv[i]); exit(1); }
    }
}

/* Index into flat 2D array */
#define IDX(i,j) ((i)*Ny + (j))

/* -dV/dphi_a where V = (mu/2)P^2/(1+kappa P^2), P = phi1*phi2*phi3 */
static inline double force_pot(double p1, double p2, double p3, int a)
{
    double P  = p1 * p2 * p3;
    double P2 = P * P;
    double denom2 = (1.0 + kappa * P2) * (1.0 + kappa * P2);
    double dP;
    switch (a) {
        case 0: dP = p2 * p3; break;
        case 1: dP = p1 * p3; break;
        case 2: dP = p1 * p2; break;
        default: dP = 0.0;
    }
    return -mu * P * dP / denom2;
}

/*
 * Compute topological charge density:
 *   rho_top = eps_{abc} phi_c [ (d_x phi_a)(d_y phi_b) - (d_y phi_a)(d_x phi_b) ]
 *
 * Expand eps_{abc} (6 terms, 3 even + 3 odd permutations):
 *   (a,b,c) = (1,2,3): +phi3 [ dx_phi1 dy_phi2 - dy_phi1 dx_phi2 ]
 *   (a,b,c) = (2,3,1): +phi1 [ dx_phi2 dy_phi3 - dy_phi2 dx_phi3 ]
 *   (a,b,c) = (3,1,2): +phi2 [ dx_phi3 dy_phi1 - dy_phi3 dx_phi1 ]
 *   (a,b,c) = (2,1,3): -phi3 [ dx_phi2 dy_phi1 - dy_phi2 dx_phi1 ]
 *   (a,b,c) = (1,3,2): -phi2 [ dx_phi1 dy_phi3 - dy_phi1 dx_phi3 ]
 *   (a,b,c) = (3,2,1): -phi1 [ dx_phi3 dy_phi2 - dy_phi3 dx_phi2 ]
 *
 * Note: even perms give same cross-product as odd perms (just with sign flip),
 * so rho_top = 2 * { phi3[dx1 dy2 - dy1 dx2] + phi1[dx2 dy3 - dy2 dx3]
 *                   + phi2[dx3 dy1 - dy3 dx1] }
 *
 * Actually let's just compute it carefully. Define:
 *   J_{ab} = (d_x phi_a)(d_y phi_b) - (d_y phi_a)(d_x phi_b)
 * Then rho_top = eps_{abc} phi_c J_{ab}
 *              = phi_3 J_{12} + phi_1 J_{23} + phi_2 J_{31}
 *              - phi_3 J_{21} - phi_2 J_{13} - phi_1 J_{32}
 * Since J_{ab} = -J_{ba}:
 *   rho_top = 2 [ phi_3 J_{12} + phi_1 J_{23} + phi_2 J_{31} ]
 * where J_{ab} = (d_x phi_a)(d_y phi_b) - (d_y phi_a)(d_x phi_b)
 */
static double compute_topo_charge(double **phi, double dx, double dy,
                                   double *rho_peak, double *rho_rms)
{
    double Q = 0;
    double peak = 0, sum2 = 0;
    int count = 0;

    for (int i = 1; i < Nx-1; i++)
        for (int j = 1; j < Ny-1; j++) {
            /* Gradients */
            double dx1 = (phi[0][IDX(i+1,j)] - phi[0][IDX(i-1,j)]) / (2.0*dx);
            double dy1 = (phi[0][IDX(i,j+1)] - phi[0][IDX(i,j-1)]) / (2.0*dy);
            double dx2 = (phi[1][IDX(i+1,j)] - phi[1][IDX(i-1,j)]) / (2.0*dx);
            double dy2 = (phi[1][IDX(i,j+1)] - phi[1][IDX(i,j-1)]) / (2.0*dy);
            double dx3 = (phi[2][IDX(i+1,j)] - phi[2][IDX(i-1,j)]) / (2.0*dx);
            double dy3 = (phi[2][IDX(i,j+1)] - phi[2][IDX(i,j-1)]) / (2.0*dy);

            /* Cross products J_{ab} = dx_a * dy_b - dy_a * dx_b */
            double J12 = dx1 * dy2 - dy1 * dx2;
            double J23 = dx2 * dy3 - dy2 * dx3;
            double J31 = dx3 * dy1 - dy3 * dx1;

            double rho = 2.0 * (phi[2][IDX(i,j)] * J12
                              + phi[0][IDX(i,j)] * J23
                              + phi[1][IDX(i,j)] * J31);

            Q += rho * dx * dy;
            if (fabs(rho) > peak) peak = fabs(rho);
            sum2 += rho * rho;
            count++;
        }

    *rho_peak = peak;
    *rho_rms  = (count > 0) ? sqrt(sum2 / count) : 0.0;
    return Q;
}

/*
 * Also compute the topological current spatial components:
 *   j^x = -eps_{abc} (dt phi_a)(dy phi_b) phi_c  (times 2 for antisym)
 *   j^y = +eps_{abc} (dt phi_a)(dx phi_b) phi_c  (times 2 for antisym)
 * These are analogous to the 2+1D reduction; useful for checking conservation.
 */

int main(int argc, char **argv)
{
    parse_args(argc, argv);

    double dx = 2.0 * L / (Nx - 1);
    double dy = 2.0 * L / (Ny - 1);
    double dx2 = dx * dx;
    double dy2 = dy * dy;
    double m2 = mass * mass;

    /* CFL for 2D */
    double cfl = 1.0 / sqrt(1.0/dx2 + 1.0/dy2 + m2);
    double dt = 0.4 * cfl;
    int Nt = (int)(tfinal / dt) + 1;

    int N = Nx * Ny;

    printf("maxwell_c: 2D oscillon topological charge\n");
    printf("  mu=%.3f kappa=%.4f mass=%.4f A=%.3f sigma=%.3f\n",
           mu, kappa, mass, A_init, sigma);
    printf("  Nx=%d Ny=%d L=%.1f dx=%.5f dy=%.5f\n", Nx, Ny, L, dx, dy);
    printf("  dt=%.6f Nt=%d tfinal=%.0f\n", dt, Nt, tfinal);
    printf("  Grid points: %d  Memory: %.1f MB\n",
           N, 10.0 * N * sizeof(double) / (1024.0*1024.0));

    /* Allocate: 3 fields x (phi, vel, acc) + damp */
    double *phi[3], *vel[3], *acc[3];
    for (int a = 0; a < 3; a++) {
        phi[a] = calloc(N, sizeof(double));
        vel[a] = calloc(N, sizeof(double));
        acc[a] = calloc(N, sizeof(double));
        if (!phi[a] || !vel[a] || !acc[a]) {
            fprintf(stderr, "Allocation failed\n");
            return 1;
        }
    }

    /* Absorbing boundary: damp in outer 25% */
    double *damp = malloc(N * sizeof(double));
    double abs_start = L * 0.75;
    for (int i = 0; i < Nx; i++) {
        double x = -L + i * dx;
        double ax = fabs(x);
        double fx = (ax > abs_start) ? (ax - abs_start) / (L - abs_start) : 0.0;
        for (int j = 0; j < Ny; j++) {
            double y = -L + j * dy;
            double ay = fabs(y);
            double fy = (ay > abs_start) ? (ay - abs_start) / (L - abs_start) : 0.0;
            double f = (fx > fy) ? fx : fy;
            damp[IDX(i,j)] = 1.0 - 0.98 * f * f;
        }
    }

    /* Initialize: spherically symmetric Gaussian oscillon
     * All three fields identical: phi_1 = phi_2 = phi_3 = A * exp(-r^2/(2*sigma^2))
     * This is the KEY configuration: symmetric so rho_top should vanish. */
    for (int i = 0; i < Nx; i++) {
        double x = -L + i * dx;
        for (int j = 0; j < Ny; j++) {
            double y = -L + j * dy;
            double r2 = x*x + y*y;
            double g = A_init * exp(-r2 / (2.0 * sigma * sigma));
            phi[0][IDX(i,j)] = g;
            phi[1][IDX(i,j)] = g;
            phi[2][IDX(i,j)] = g;
        }
    }

    /* Compute initial topological charge */
    double rho_peak, rho_rms;
    double Q0 = compute_topo_charge(phi, dx, dy, &rho_peak, &rho_rms);
    printf("\n  Initial topological charge Q = %.6e\n", Q0);
    printf("  Initial |rho_top| peak = %.6e, rms = %.6e\n", rho_peak, rho_rms);

    /* Analytical argument printout */
    printf("\n  --- Analytical check ---\n");
    printf("  When phi_1=phi_2=phi_3=f(r), all gradients are parallel:\n");
    printf("    d_x phi_a = f'(r) * x/r  for all a\n");
    printf("    d_y phi_a = f'(r) * y/r  for all a\n");
    printf("  So J_{ab} = (d_x phi_a)(d_y phi_b) - (d_y phi_a)(d_x phi_b)\n");
    printf("            = f'^2 (x/r)(y/r) - f'^2 (y/r)(x/r) = 0  for all a,b\n");
    printf("  Therefore rho_top = 0 identically (no topology in symmetric config)\n");
    printf("  ---\n\n");

    /* Compute acceleration macro */
    #define COMPUTE_ACC() do { \
        for (int a = 0; a < 3; a++) { \
            for (int ii = 0; ii < Nx; ii++) { \
                acc[a][IDX(ii,0)] = 0; \
                acc[a][IDX(ii,Ny-1)] = 0; \
            } \
            for (int jj = 0; jj < Ny; jj++) { \
                acc[a][IDX(0,jj)] = 0; \
                acc[a][IDX(Nx-1,jj)] = 0; \
            } \
            for (int ii = 1; ii < Nx-1; ii++) \
                for (int jj = 1; jj < Ny-1; jj++) { \
                    double lapl = (phi[a][IDX(ii+1,jj)] + phi[a][IDX(ii-1,jj)] \
                                   - 2.0*phi[a][IDX(ii,jj)]) / dx2 \
                                + (phi[a][IDX(ii,jj+1)] + phi[a][IDX(ii,jj-1)] \
                                   - 2.0*phi[a][IDX(ii,jj)]) / dy2; \
                    double fp = force_pot(phi[0][IDX(ii,jj)], phi[1][IDX(ii,jj)], \
                                          phi[2][IDX(ii,jj)], a); \
                    acc[a][IDX(ii,jj)] = lapl - m2*phi[a][IDX(ii,jj)] + fp; \
                } \
        } \
    } while(0)

    COMPUTE_ACC();

    /* Open time series file */
    char tspath[600];
    snprintf(tspath, sizeof(tspath), "%s/maxwell_c_ts.tsv", outdir);
    FILE *fts = fopen(tspath, "w");
    if (!fts) {
        fprintf(stderr, "Cannot open %s — creating directory\n", tspath);
        char cmd[700];
        snprintf(cmd, sizeof(cmd), "mkdir -p %s", outdir);
        system(cmd);
        fts = fopen(tspath, "w");
        if (!fts) { fprintf(stderr, "Still can't open %s\n", tspath); return 1; }
    }
    fprintf(fts, "time\tphi1_center\tpeak_phi\tE_kin\tE_grad\tE_mass\tE_pot\tE_total\t"
                 "Q_top\trho_top_peak\trho_top_rms\tf_core\n");

    int rec_every  = Nt / 10000;
    if (rec_every < 1) rec_every = 1;
    int print_every = Nt / 40;
    if (print_every < 1) print_every = 1;

    /* Snapshot times */
    double snap_times[] = {0, 10, 50, 100, 200, 500, 1000, 2000};
    int n_snaps = sizeof(snap_times) / sizeof(snap_times[0]);
    int next_snap = 0;

    /* DFT storage for spectrum */
    int max_dft = 50000;
    double *phi0_hist = malloc(max_dft * sizeof(double));
    double *t_hist    = malloc(max_dft * sizeof(double));
    int n_dft = 0;
    int dft_every = Nt / max_dft;
    if (dft_every < 1) dft_every = 1;

    int ic = Nx / 2, jc = Ny / 2;  /* center indices */
    double core_r2 = 9.0 * sigma * sigma;  /* 3*sigma */

    /* Track max |Q_top| and max rho_rms over the run */
    double max_Q = 0, max_rho_rms = 0, max_rho_peak_all = 0;

    for (int n = 0; n <= Nt; n++) {
        double t = n * dt;

        /* DFT sampling */
        if (n % dft_every == 0 && n_dft < max_dft) {
            phi0_hist[n_dft] = phi[0][IDX(ic,jc)];
            t_hist[n_dft] = t;
            n_dft++;
        }

        /* Snapshots — write rho_top profile */
        if (next_snap < n_snaps && t >= snap_times[next_snap] - 0.5*dt) {
            char snappath[600];
            snprintf(snappath, sizeof(snappath),
                     "%s/maxwell_c_profile_t%04d.tsv", outdir, (int)snap_times[next_snap]);
            FILE *fsnap = fopen(snappath, "w");
            if (fsnap) {
                fprintf(fsnap, "x\ty\tphi1\tphi2\tphi3\trho_top\n");
                int step = 4;
                for (int i = 1; i < Nx-1; i += step)
                    for (int j = 1; j < Ny-1; j += step) {
                        double x = -L + i * dx;
                        double y = -L + j * dy;
                        /* Compute rho_top at this point */
                        double dx1 = (phi[0][IDX(i+1,j)] - phi[0][IDX(i-1,j)]) / (2.0*dx);
                        double dy1 = (phi[0][IDX(i,j+1)] - phi[0][IDX(i,j-1)]) / (2.0*dy);
                        double dx2 = (phi[1][IDX(i+1,j)] - phi[1][IDX(i-1,j)]) / (2.0*dx);
                        double dy2 = (phi[1][IDX(i,j+1)] - phi[1][IDX(i,j-1)]) / (2.0*dy);
                        double dx3 = (phi[2][IDX(i+1,j)] - phi[2][IDX(i-1,j)]) / (2.0*dx);
                        double dy3 = (phi[2][IDX(i,j+1)] - phi[2][IDX(i,j-1)]) / (2.0*dy);
                        double J12 = dx1 * dy2 - dy1 * dx2;
                        double J23 = dx2 * dy3 - dy2 * dx3;
                        double J31 = dx3 * dy1 - dy3 * dx1;
                        double rho = 2.0 * (phi[2][IDX(i,j)] * J12
                                          + phi[0][IDX(i,j)] * J23
                                          + phi[1][IDX(i,j)] * J31);
                        fprintf(fsnap, "%.4f\t%.4f\t%.6e\t%.6e\t%.6e\t%.6e\n",
                                x, y, phi[0][IDX(i,j)], phi[1][IDX(i,j)],
                                phi[2][IDX(i,j)], rho);
                    }
                fclose(fsnap);
                printf("  Snapshot t=%g written: %s\n", snap_times[next_snap], snappath);
            }
            next_snap++;
        }

        int do_rec   = (n % rec_every == 0);
        int do_print = (n % print_every == 0);

        if (do_rec || do_print) {
            double Ek = 0, Eg = 0, Em = 0, Ep = 0;
            double peak_phi = 0;
            double Ecore = 0, Eall = 0;

            for (int i = 1; i < Nx-1; i++)
                for (int j = 1; j < Ny-1; j++) {
                    int k = IDX(i,j);
                    double x = -L + i * dx;
                    double y = -L + j * dy;
                    double r2_here = x*x + y*y;

                    double e_local = 0;
                    for (int a = 0; a < 3; a++) {
                        Ek += 0.5 * vel[a][k] * vel[a][k] * dx * dy;
                        double dpx = (phi[a][IDX(i+1,j)] - phi[a][IDX(i-1,j)]) / (2.0*dx);
                        double dpy = (phi[a][IDX(i,j+1)] - phi[a][IDX(i,j-1)]) / (2.0*dy);
                        Eg += 0.5 * (dpx*dpx + dpy*dpy) * dx * dy;
                        Em += 0.5 * m2 * phi[a][k] * phi[a][k] * dx * dy;
                        if (fabs(phi[a][k]) > peak_phi) peak_phi = fabs(phi[a][k]);

                        e_local += 0.5*vel[a][k]*vel[a][k] + 0.5*(dpx*dpx+dpy*dpy)
                                 + 0.5*m2*phi[a][k]*phi[a][k];
                    }
                    double P = phi[0][k] * phi[1][k] * phi[2][k];
                    double P2 = P * P;
                    double V = 0.5 * mu * P2 / (1.0 + kappa * P2);
                    Ep += V * dx * dy;
                    e_local += V;
                    Eall += e_local * dx * dy;
                    if (r2_here < core_r2) Ecore += e_local * dx * dy;
                }

            double Et = Ek + Eg + Em + Ep;
            double fc = (Eall > 1e-20) ? Ecore / Eall : 0.0;

            /* Topological charge */
            double rp, rrms;
            double Qtop = compute_topo_charge(phi, dx, dy, &rp, &rrms);
            if (fabs(Qtop) > max_Q) max_Q = fabs(Qtop);
            if (rrms > max_rho_rms) max_rho_rms = rrms;
            if (rp > max_rho_peak_all) max_rho_peak_all = rp;

            if (do_rec)
                fprintf(fts, "%.4f\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t"
                             "%.6e\t%.6e\t%.6e\t%.6f\n",
                        t, phi[0][IDX(ic,jc)], peak_phi,
                        Ek, Eg, Em, Ep, Et,
                        Qtop, rp, rrms, fc);

            if (do_print)
                printf("  t=%7.1f  phi_c=%.4f  pk=%.4f  E=%.3f(k=%.3f g=%.3f m=%.3f p=%.3f)"
                       "  Q_top=%.2e  rho_pk=%.2e  fc=%.3f\n",
                       t, phi[0][IDX(ic,jc)], peak_phi,
                       Et, Ek, Eg, Em, Ep,
                       Qtop, rp, fc);
        }

        if (n == Nt) break;

        /* Velocity Verlet */
        for (int a = 0; a < 3; a++)
            for (int k = 0; k < N; k++)
                vel[a][k] += 0.5 * dt * acc[a][k];
        for (int a = 0; a < 3; a++)
            for (int k = 0; k < N; k++)
                phi[a][k] += dt * vel[a][k];
        COMPUTE_ACC();
        for (int a = 0; a < 3; a++)
            for (int k = 0; k < N; k++)
                vel[a][k] += 0.5 * dt * acc[a][k];

        /* Absorbing boundary */
        for (int a = 0; a < 3; a++)
            for (int k = 0; k < N; k++) {
                vel[a][k] *= damp[k];
                phi[a][k] *= damp[k];
            }
    }

    fclose(fts);

    /* Summary */
    printf("\n===== TOPOLOGICAL CHARGE SUMMARY =====\n");
    printf("  Max |Q_top| over run: %.6e\n", max_Q);
    printf("  Max rho_top_rms:      %.6e\n", max_rho_rms);
    printf("  Max rho_top_peak:     %.6e\n", max_rho_peak_all);
    printf("\n  RESULT: ");
    if (max_Q < 1e-10 && max_rho_peak_all < 1e-10)
        printf("Q_top = 0 to machine precision. Oscillon has NO topological charge.\n");
    else if (max_Q < 1e-6)
        printf("Q_top ~ %.1e (numerical noise only). No physical topological charge.\n", max_Q);
    else
        printf("Q_top ~ %.3e — NONZERO topological charge detected!\n", max_Q);

    printf("\n  EXPLANATION:\n");
    printf("  For phi_1=phi_2=phi_3 (symmetric init), the gradients are identical:\n");
    printf("    d_i phi_a = d_i phi_b for all a,b\n");
    printf("  So J_{ab} = (d_x phi_a)(d_y phi_b) - (d_y phi_a)(d_x phi_b) = 0\n");
    printf("  because the cross product of parallel vectors vanishes.\n");
    printf("  rho_top = eps_{abc} phi_c * J_{ab} = 0 identically.\n");
    printf("  Even after time evolution, the symmetry phi_1=phi_2=phi_3 is preserved\n");
    printf("  because the EOM and initial data treat all three fields identically.\n");
    printf("  Topological charge requires DISTINCT field configurations (e.g., a baby\n");
    printf("  Skyrmion where the fields map S^2 -> S^2 with winding number).\n");
    printf("======================================\n");

    /* DFT of phi_1(center, t) — second half for spectrum */
    int dft_start = n_dft / 2;
    if (n_dft - dft_start > 100) {
        char dftpath[600];
        snprintf(dftpath, sizeof(dftpath), "%s/maxwell_c_spectrum.tsv", outdir);
        FILE *fdft = fopen(dftpath, "w");
        fprintf(fdft, "omega\tpower\n");
        double T = t_hist[n_dft-1] - t_hist[dft_start];
        int nf = 500;
        double peak_pow = 0, peak_om = 0;
        for (int k = 0; k < nf; k++) {
            double omega = 3.0 * mass * k / nf;
            double re = 0, im = 0;
            for (int jj = dft_start; jj < n_dft; jj++) {
                double dtj = (jj > dft_start) ?
                    (t_hist[jj]-t_hist[jj-1]) : (t_hist[dft_start+1]-t_hist[dft_start]);
                re += phi0_hist[jj] * cos(omega * t_hist[jj]) * dtj;
                im += phi0_hist[jj] * sin(omega * t_hist[jj]) * dtj;
            }
            double pw = (re*re + im*im) / (T*T);
            fprintf(fdft, "%.6f\t%.6e\n", omega, pw);
            if (pw > peak_pow) { peak_pow = pw; peak_om = omega; }
        }
        fclose(fdft);
        printf("\nSpectrum: peak omega = %.4f (mass gap = %.4f)\n", peak_om, mass);
        printf("Oscillon (omega < m)? %s\n",
               (peak_om > 0.01 && peak_om < mass) ? "YES" : "NO");
        printf("Output: %s\n", dftpath);
    }

    printf("Time series: %s\n", tspath);

    /* Cleanup */
    for (int a = 0; a < 3; a++) { free(phi[a]); free(vel[a]); free(acc[a]); }
    free(damp); free(phi0_hist); free(t_hist);
    return 0;
}
