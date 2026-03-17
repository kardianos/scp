/*
 * dynAC.c — V26-DynAC: Propagating Massless Braid
 *
 * Helical wave propagating at c along z-axis (m=0, omega=k).
 * Triple-product potential only; periodic BC in z, absorbing in x,y.
 *
 * Lagrangian:
 *   L = Sum_a [1/2(dt phi_a)^2 - 1/2|grad phi_a|^2]
 *     - V(P)     where P = phi_1 phi_2 phi_3
 *     V(P) = (1/2) mu P^2 / (1 + kappa P^2)
 *
 * Init:
 *   phi_a(x,0) = A(r_perp) cos(kz + 2 pi a/3)
 *   vel_a(x,0) = k A(r_perp) sin(kz + 2 pi a/3)   [omega=k for m=0]
 *
 * Compile: gcc -O3 -fopenmp -Wall -o dynAC src/dynAC.c -lm
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>

/* ─── Parameters ─── */
static double mu_pot    = -20.0;
static double kappa     = 20.0;
static double mass      = 0.0;
static double A0        = 0.8;
static double R_tube    = 3.0;    /* transverse Gaussian width */

static int    N         = 128;
static double L         = 20.0;
static double tfinal    = 500.0;
static double cfl_frac  = 0.20;
static char   outdir[512] = "data";

static void parse_args(int argc, char **argv)
{
    for (int i = 1; i < argc - 1; i += 2) {
        if      (!strcmp(argv[i], "-mu"))       mu_pot    = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-kappa"))    kappa     = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-mass"))     mass      = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-A"))        A0        = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-N"))        N         = atoi(argv[i+1]);
        else if (!strcmp(argv[i], "-L"))        L         = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-tfinal"))   tfinal    = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-cfl"))      cfl_frac  = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-Rtube"))    R_tube    = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-o"))        strncpy(outdir, argv[i+1], sizeof(outdir)-1);
        else { fprintf(stderr, "Unknown flag: %s\n", argv[i]); exit(1); }
    }
}

/* ─── Index helpers ─── */
#define IDX(i,j,k) ((long)(i)*N*N + (long)(j)*N + (long)(k))

/* ─── Globals ─── */
static double *phi[3], *vel[3], *acc[3];
static double *damp;
static double dx, dx2, m2, dt;
static long Ngrid;
static double Lz;  /* z domain length = 2L */

/* ─── Periodic z wrap ─── */
static inline int wrap_z(int k)
{
    if (k < 0)   return k + N;
    if (k >= N)  return k - N;
    return k;
}

/* ─── Initialization: Propagating helical braid ─── */
static void init_propagating(void)
{
    double k = 2.0 * M_PI / Lz;

    printf("  Init: propagating helical braid, k=%.6f, Lz=%.2f\n", k, Lz);
    printf("  omega=k=%.6f (massless dispersion), v_phase=c=1\n", k);

    #pragma omp parallel for schedule(static)
    for (long idx = 0; idx < Ngrid; idx++) {
        int ii = idx / (N * N);
        int jj = (idx / N) % N;
        int kk = idx % N;
        double x = -L + ii * dx;
        double y = -L + jj * dx;
        double z = -L + kk * dx;

        double rperp2 = x * x + y * y;
        double envelope = A0 * exp(-rperp2 / (2.0 * R_tube * R_tube));

        for (int a = 0; a < 3; a++) {
            double phase = k * z + 2.0 * M_PI * a / 3.0;
            phi[a][idx] = envelope * cos(phase);
            vel[a][idx] = k * envelope * sin(phase);
        }
    }
}

/* ─── Compute acceleration ─── */
static void compute_acc(void)
{
    for (int a = 0; a < 3; a++) {
        #pragma omp parallel for schedule(static)
        for (int i = 2; i < N-2; i++) {
            for (int j = 2; j < N-2; j++) {
                for (int k = 0; k < N; k++) {
                    long idx = IDX(i, j, k);

                    /* 7-point Laplacian with periodic z */
                    int km1 = wrap_z(k - 1);
                    int kp1 = wrap_z(k + 1);

                    double lapl = (phi[a][IDX(i+1,j,k)] + phi[a][IDX(i-1,j,k)]
                                 + phi[a][IDX(i,j+1,k)] + phi[a][IDX(i,j-1,k)]
                                 + phi[a][IDX(i,j,kp1)] + phi[a][IDX(i,j,km1)]
                                 - 6.0 * phi[a][idx]) / dx2;

                    /* Triple-product potential force */
                    double p0 = phi[0][idx], p1 = phi[1][idx], p2 = phi[2][idx];
                    double P = p0 * p1 * p2;
                    double P2 = P * P;
                    double denom2 = (1.0 + kappa * P2);
                    denom2 *= denom2;

                    double dP;
                    switch (a) {
                        case 0: dP = p1 * p2; break;
                        case 1: dP = p0 * p2; break;
                        default: dP = p0 * p1; break;
                    }
                    double dVdphi = mu_pot * P * dP / denom2;

                    acc[a][idx] = lapl - m2 * phi[a][idx] - dVdphi;
                }
            }
        }

        /* Boundary: zero acceleration for non-periodic (x,y) directions */
        #pragma omp parallel for schedule(static)
        for (long idx = 0; idx < Ngrid; idx++) {
            int ii = idx / (N * N);
            int jj = (idx / N) % N;
            if (ii < 2 || ii >= N-2 || jj < 2 || jj >= N-2) {
                acc[a][idx] = 0.0;
            }
        }
    }
}

/* ─── Velocity Verlet step ─── */
static void verlet_step(void)
{
    /* 1. Half-kick */
    for (int a = 0; a < 3; a++) {
        #pragma omp parallel for schedule(static)
        for (long idx = 0; idx < Ngrid; idx++)
            vel[a][idx] += 0.5 * dt * acc[a][idx];
    }

    /* 2. Drift */
    for (int a = 0; a < 3; a++) {
        #pragma omp parallel for schedule(static)
        for (long idx = 0; idx < Ngrid; idx++)
            phi[a][idx] += dt * vel[a][idx];
    }

    /* 3. Recompute acceleration */
    compute_acc();

    /* 4. Half-kick */
    for (int a = 0; a < 3; a++) {
        #pragma omp parallel for schedule(static)
        for (long idx = 0; idx < Ngrid; idx++)
            vel[a][idx] += 0.5 * dt * acc[a][idx];
    }

    /* 5. Absorbing boundary damping (x,y only; z is periodic) */
    for (int a = 0; a < 3; a++) {
        #pragma omp parallel for schedule(static)
        for (long idx = 0; idx < Ngrid; idx++) {
            vel[a][idx] *= damp[idx];
            phi[a][idx] *= damp[idx];
        }
    }
}

/* ─── Diagnostics structure ─── */
typedef struct {
    double Ek, Eg, Em, Ep, Et;
    double peak[3];
    double peak_P;
    double phi_center[3];
    double fc;             /* core fraction */
    double rho_center;     /* energy density at center */
    double z_com;          /* z center-of-energy (momentum tracking) */
} Diag;

/* ─── Compute diagnostics ─── */
static Diag compute_diag(double core_radius)
{
    Diag d;
    memset(&d, 0, sizeof(d));
    long ic = IDX(N/2, N/2, N/2);
    for (int a = 0; a < 3; a++)
        d.phi_center[a] = phi[a][ic];

    double Ek=0, Eg=0, Em=0, Ep=0;
    double Ecore=0, Eall=0;
    double peak[3] = {0,0,0};
    double peak_P = 0;
    double z_num = 0;  /* weighted z for COM */

    /* Energy density at center */
    double rho_ctr = 0;
    {
        int ci = N/2, cj = N/2, ck = N/2;
        for (int a = 0; a < 3; a++) {
            rho_ctr += 0.5 * vel[a][ic] * vel[a][ic];
            double gx = (phi[a][IDX(ci+1,cj,ck)] - phi[a][IDX(ci-1,cj,ck)]) / (2*dx);
            double gy = (phi[a][IDX(ci,cj+1,ck)] - phi[a][IDX(ci,cj-1,ck)]) / (2*dx);
            int ckp = wrap_z(ck+1), ckm = wrap_z(ck-1);
            double gz = (phi[a][IDX(ci,cj,ckp)] - phi[a][IDX(ci,cj,ckm)]) / (2*dx);
            rho_ctr += 0.5 * (gx*gx + gy*gy + gz*gz);
            rho_ctr += 0.5 * m2 * phi[a][ic] * phi[a][ic];
        }
        double p0c = phi[0][ic], p1c = phi[1][ic], p2c = phi[2][ic];
        double Pc = p0c * p1c * p2c;
        double Pc2 = Pc * Pc;
        rho_ctr += 0.5 * mu_pot * Pc2 / (1.0 + kappa * Pc2);
    }

    #pragma omp parallel
    {
        double lEk=0, lEg=0, lEm=0, lEp=0, lEc=0, lEa=0;
        double lpk[3] = {0,0,0};
        double lpkP = 0;
        double lz_num = 0;

        #pragma omp for schedule(static) nowait
        for (int i = 2; i < N-2; i++) {
            double x = -L + i * dx;
            for (int j = 2; j < N-2; j++) {
                double y = -L + j * dx;
                for (int k = 0; k < N; k++) {
                    double z = -L + k * dx;
                    long idx = IDX(i,j,k);
                    double dV = dx * dx * dx;
                    double e_loc = 0;

                    for (int a = 0; a < 3; a++) {
                        double v2 = vel[a][idx] * vel[a][idx];
                        lEk += 0.5 * v2 * dV;
                        e_loc += 0.5 * v2;

                        double gx = (phi[a][IDX(i+1,j,k)] - phi[a][IDX(i-1,j,k)]) / (2*dx);
                        double gy = (phi[a][IDX(i,j+1,k)] - phi[a][IDX(i,j-1,k)]) / (2*dx);
                        int kp1 = wrap_z(k+1), km1 = wrap_z(k-1);
                        double gz = (phi[a][IDX(i,j,kp1)] - phi[a][IDX(i,j,km1)]) / (2*dx);
                        double grad2 = gx*gx + gy*gy + gz*gz;
                        lEg += 0.5 * grad2 * dV;
                        e_loc += 0.5 * grad2;

                        double mass_e = 0.5 * m2 * phi[a][idx] * phi[a][idx];
                        lEm += mass_e * dV;
                        e_loc += mass_e;

                        double ap = fabs(phi[a][idx]);
                        if (ap > lpk[a]) lpk[a] = ap;
                    }

                    double p0 = phi[0][idx], p1 = phi[1][idx], p2 = phi[2][idx];
                    double P = p0 * p1 * p2;
                    double P2 = P * P;
                    double Vloc = 0.5 * mu_pot * P2 / (1.0 + kappa * P2);
                    lEp += Vloc * dV;
                    e_loc += Vloc;

                    double absP = fabs(P);
                    if (absP > lpkP) lpkP = absP;

                    lEa += e_loc * dV;

                    double rperp = sqrt(x*x + y*y);
                    if (rperp < core_radius) lEc += e_loc * dV;

                    /* z-COM weighted by energy density */
                    lz_num += z * e_loc * dV;
                }
            }
        }

        #pragma omp critical
        {
            Ek += lEk; Eg += lEg; Em += lEm; Ep += lEp;
            Ecore += lEc; Eall += lEa;
            z_num += lz_num;
            for (int a = 0; a < 3; a++)
                if (lpk[a] > peak[a]) peak[a] = lpk[a];
            if (lpkP > peak_P) peak_P = lpkP;
        }
    }

    d.Ek = Ek; d.Eg = Eg; d.Em = Em; d.Ep = Ep;
    d.Et = Ek + Eg + Em + Ep;
    for (int a = 0; a < 3; a++) d.peak[a] = peak[a];
    d.fc = (Eall > 1e-20) ? Ecore / Eall : 0.0;
    d.peak_P = peak_P;
    d.rho_center = rho_ctr;
    d.z_com = (Eall > 1e-20) ? z_num / Eall : 0.0;

    return d;
}

/* ─── Fast center density ─── */
static double center_rho(void)
{
    long ic = IDX(N/2, N/2, N/2);
    int ci = N/2, cj = N/2, ck = N/2;
    double rho = 0;
    for (int a = 0; a < 3; a++) {
        rho += 0.5 * vel[a][ic] * vel[a][ic];
        double gx = (phi[a][IDX(ci+1,cj,ck)] - phi[a][IDX(ci-1,cj,ck)]) / (2*dx);
        double gy = (phi[a][IDX(ci,cj+1,ck)] - phi[a][IDX(ci,cj-1,ck)]) / (2*dx);
        int ckp = wrap_z(ck+1), ckm = wrap_z(ck-1);
        double gz = (phi[a][IDX(ci,cj,ckp)] - phi[a][IDX(ci,cj,ckm)]) / (2*dx);
        rho += 0.5 * (gx*gx + gy*gy + gz*gz);
        rho += 0.5 * m2 * phi[a][ic] * phi[a][ic];
    }
    double p0 = phi[0][ic], p1 = phi[1][ic], p2 = phi[2][ic];
    double P = p0 * p1 * p2;
    double P2 = P * P;
    rho += 0.5 * mu_pot * P2 / (1.0 + kappa * P2);
    return rho;
}

/* ─── DFT peak finder ─── */
static double find_peak_omega(double *hist, double *t_hist, int n_pts, int start,
                              double omega_max, double *peak_power_out)
{
    double T = t_hist[n_pts-1] - t_hist[start];
    if (T < 10.0 || n_pts - start < 50) return 0.0;

    double mean = 0;
    for (int j = start; j < n_pts; j++) mean += hist[j];
    mean /= (n_pts - start);

    int nf = 500;
    double peak_pow = 0.0, peak_om = 0.0;
    for (int kk = 1; kk < nf; kk++) {
        double omega = omega_max * kk / nf;
        double re = 0, im = 0;
        for (int j = start; j < n_pts; j++) {
            double dtj = (j > start) ?
                (t_hist[j]-t_hist[j-1]) : (t_hist[start+1]-t_hist[start]);
            double val = hist[j] - mean;
            re += val * cos(omega * t_hist[j]) * dtj;
            im += val * sin(omega * t_hist[j]) * dtj;
        }
        double pw = re*re + im*im;
        if (pw > peak_pow) { peak_pow = pw; peak_om = omega; }
    }
    if (peak_power_out) *peak_power_out = peak_pow;
    return peak_om;
}

/* ─── Measure z-momentum (Pz = -sum_a integral dt_phi_a dz_phi_a) ─── */
static double compute_Pz(void)
{
    double Pz = 0;
    #pragma omp parallel for reduction(+:Pz) schedule(static)
    for (int i = 2; i < N-2; i++) {
        for (int j = 2; j < N-2; j++) {
            for (int k = 0; k < N; k++) {
                long idx = IDX(i,j,k);
                int kp1 = wrap_z(k+1), km1 = wrap_z(k-1);
                double dV = dx * dx * dx;
                for (int a = 0; a < 3; a++) {
                    double gz = (phi[a][IDX(i,j,kp1)] - phi[a][IDX(i,j,km1)]) / (2*dx);
                    Pz += vel[a][idx] * gz * dV;
                }
            }
        }
    }
    return -Pz;  /* field momentum = -integral(dot_phi grad_phi) */
}

/* ─── Compute transverse RMS radius ─── */
static double compute_rperp_rms(void)
{
    double num = 0, den = 0;
    #pragma omp parallel for reduction(+:num,den) schedule(static)
    for (int i = 2; i < N-2; i++) {
        double x = -L + i * dx;
        for (int j = 2; j < N-2; j++) {
            double y = -L + j * dx;
            double rperp2 = x*x + y*y;
            for (int k = 0; k < N; k++) {
                long idx = IDX(i,j,k);
                double dV = dx * dx * dx;
                double e_loc = 0;
                for (int a = 0; a < 3; a++) {
                    e_loc += 0.5 * vel[a][idx] * vel[a][idx];
                    double gx = (phi[a][IDX(i+1,j,k)] - phi[a][IDX(i-1,j,k)]) / (2*dx);
                    double gy = (phi[a][IDX(i,j+1,k)] - phi[a][IDX(i,j-1,k)]) / (2*dx);
                    int kp1 = wrap_z(k+1), km1 = wrap_z(k-1);
                    double gz = (phi[a][IDX(i,j,kp1)] - phi[a][IDX(i,j,km1)]) / (2*dx);
                    e_loc += 0.5 * (gx*gx + gy*gy + gz*gz);
                }
                double p0 = phi[0][idx], p1 = phi[1][idx], p2 = phi[2][idx];
                double P = p0 * p1 * p2;
                double P2 = P * P;
                e_loc += 0.5 * mu_pot * P2 / (1.0 + kappa * P2);

                num += rperp2 * e_loc * dV;
                den += e_loc * dV;
            }
        }
    }
    return (den > 1e-20) ? sqrt(num / den) : 0.0;
}

/* ─── Main ─── */
int main(int argc, char **argv)
{
    parse_args(argc, argv);

    dx  = 2.0 * L / (N - 1);
    dx2 = dx * dx;
    m2  = mass * mass;
    dt  = cfl_frac * dx;
    Ngrid = (long)N * N * N;
    Lz  = 2.0 * L;

    double k_wave = 2.0 * M_PI / Lz;

    printf("=== V26-DynAC: Propagating Massless Braid ===\n");
    printf("Parameters:\n");
    printf("  mu=%.1f  kappa=%.1f  mass=%.3f  A0=%.3f\n", mu_pot, kappa, mass, A0);
    printf("  R_tube=%.1f  k=%.6f  Lz=%.2f\n", R_tube, k_wave, Lz);
    printf("  N=%d  L=%.1f  dx=%.4f  dt=%.5f\n", N, L, dx, dt);
    printf("  Ngrid=%ld (%.1f M)  Memory: %.1f MB\n",
           Ngrid, Ngrid/1e6, Ngrid*8.0*10/1e6);
    printf("  Threads: %d\n", omp_get_max_threads());
    printf("  tfinal=%.0f  CFL=%.2f\n", tfinal, cfl_frac);
    fflush(stdout);

    /* Allocate */
    for (int a = 0; a < 3; a++) {
        phi[a] = calloc(Ngrid, sizeof(double));
        vel[a] = calloc(Ngrid, sizeof(double));
        acc[a] = calloc(Ngrid, sizeof(double));
        if (!phi[a] || !vel[a] || !acc[a]) {
            fprintf(stderr, "Allocation failed\n");
            return 1;
        }
    }
    damp = malloc(Ngrid * sizeof(double));
    if (!damp) { fprintf(stderr, "Allocation failed for damp\n"); return 1; }

    /* Damping: absorb in transverse (x,y) only; z is periodic */
    double R_abs_inner = L * 0.70;
    double R_abs_outer = L * 0.95;
    #pragma omp parallel for schedule(static)
    for (long idx = 0; idx < Ngrid; idx++) {
        int ii = idx / (N * N);
        int jj = (idx / N) % N;
        double x = -L + ii * dx;
        double y = -L + jj * dx;
        double rperp = sqrt(x*x + y*y);
        if (rperp > R_abs_inner) {
            double f = (rperp - R_abs_inner) / (R_abs_outer - R_abs_inner);
            if (f > 1.0) f = 1.0;
            damp[idx] = 1.0 - 0.98 * f * f;
        } else {
            damp[idx] = 1.0;
        }
    }

    /* Initialize propagating helical braid */
    init_propagating();

    /* Compute initial acceleration */
    compute_acc();

    /* Time integration */
    int Nt = (int)(tfinal / dt) + 1;
    double core_radius = 8.0;

    /* DFT storage */
    int max_dft = 50000;
    double *rho_hist  = malloc(max_dft * sizeof(double));
    double *phi0_hist = malloc(max_dft * sizeof(double));
    double *t_hist    = malloc(max_dft * sizeof(double));
    int n_dft = 0;

    int diag_every  = Nt / 5000;  if (diag_every < 1)  diag_every = 1;
    int print_every = Nt / 30;    if (print_every < 1)  print_every = 1;
    int dft_every   = Nt / max_dft; if (dft_every < 1)  dft_every = 1;

    /* Output file */
    char path[600];
    snprintf(path, sizeof(path), "%s/dynAC_phase1.tsv", outdir);
    FILE *f1 = fopen(path, "w");
    if (!f1) { fprintf(stderr, "Cannot open %s\n", path); return 1; }
    fprintf(f1, "time\tE_total\tE_kin\tE_grad\tE_mass\tE_pot\t"
                "fc\tpeak0\tpeak1\tpeak2\tpeak_P\trho_center\tPz\tRperp\n");

    double wall_start = omp_get_wtime();

    /* Initial diagnostics */
    {
        Diag d0 = compute_diag(core_radius);
        double Pz0 = compute_Pz();
        double Rp0 = compute_rperp_rms();
        printf("\n--- Evolution (t=0..%.0f) ---\n", tfinal);
        printf("  t=%7.1f  E=%.2f  Ek=%.2f  Eg=%.2f  Ep=%.2f  fc=%.4f  pk=(%.4f,%.4f,%.4f)  |P|=%.6f  Pz=%.2f  Rp=%.3f\n",
               0.0, d0.Et, d0.Ek, d0.Eg, d0.Ep, d0.fc, d0.peak[0], d0.peak[1], d0.peak[2],
               d0.peak_P, Pz0, Rp0);
        fprintf(f1, "%.4f\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6f\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\n",
                0.0, d0.Et, d0.Ek, d0.Eg, d0.Em, d0.Ep, d0.fc,
                d0.peak[0], d0.peak[1], d0.peak[2], d0.peak_P, d0.rho_center, Pz0, Rp0);
    }

    for (int n = 1; n <= Nt; n++) {
        verlet_step();
        double t = n * dt;

        /* DFT history */
        if (n % dft_every == 0 && n_dft < max_dft) {
            long ic = IDX(N/2, N/2, N/2);
            rho_hist[n_dft]  = center_rho();
            phi0_hist[n_dft] = phi[0][ic];
            t_hist[n_dft]    = t;
            n_dft++;
        }

        int do_diag  = (n % diag_every == 0);
        int do_print = (n % print_every == 0) || (n == Nt);

        if (do_diag || do_print) {
            Diag d = compute_diag(core_radius);
            double Pz = compute_Pz();

            if (do_diag) {
                double Rp = (do_print) ? compute_rperp_rms() : 0.0;
                fprintf(f1, "%.4f\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6f\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\n",
                        t, d.Et, d.Ek, d.Eg, d.Em, d.Ep, d.fc,
                        d.peak[0], d.peak[1], d.peak[2], d.peak_P, d.rho_center, Pz, Rp);
            }

            if (do_print) {
                double Rp = compute_rperp_rms();
                double elapsed = omp_get_wtime() - wall_start;
                double frac = (double)n / Nt;
                double eta_t = (frac > 0.001) ? elapsed * (1.0-frac)/frac : 0;
                printf("  t=%7.1f  E=%.2f  Ek=%.2f  Eg=%.2f  Ep=%.2f  fc=%.4f  pk=(%.4f,%.4f,%.4f)  |P|=%.6f  Pz=%.2f  Rp=%.3f  [%.0fs ETA %.0fs]\n",
                       t, d.Et, d.Ek, d.Eg, d.Ep, d.fc, d.peak[0], d.peak[1], d.peak[2],
                       d.peak_P, Pz, Rp, elapsed, eta_t);
                fflush(stdout);
            }
        }
    }

    fclose(f1);

    double elapsed_total = omp_get_wtime() - wall_start;
    printf("\nEvolution complete (%.1f sec)\n", elapsed_total);

    /* Final diagnostics */
    Diag dfinal = compute_diag(core_radius);
    double Pz_final = compute_Pz();
    double Rp_final = compute_rperp_rms();
    printf("  Final: E=%.2f  fc=%.4f  pk=(%.4f,%.4f,%.4f)  |P|max=%.6f  Pz=%.2f  Rp=%.3f\n",
           dfinal.Et, dfinal.fc, dfinal.peak[0], dfinal.peak[1], dfinal.peak[2],
           dfinal.peak_P, Pz_final, Rp_final);

    int survived = (dfinal.fc > 0.01 && dfinal.peak_P > 1e-6);
    printf("  Survived? %s (fc=%.4f, |P|max=%.6f)\n",
           survived ? "YES" : "NO", dfinal.fc, dfinal.peak_P);

    /* ─── Phase 2: DFT analysis ─── */
    printf("\n--- Phase 2: Non-Breathing Verification ---\n");

    if (n_dft < 100) {
        printf("  Not enough DFT samples (%d), skipping.\n", n_dft);
    } else {
        int start = n_dft / 2;

        double peak_power_rho = 0;
        double omega_rho = find_peak_omega(rho_hist, t_hist, n_dft, start, 5.0, &peak_power_rho);

        double peak_power_phi = 0;
        double omega_phi = find_peak_omega(phi0_hist, t_hist, n_dft, start, 5.0, &peak_power_phi);

        double mean_rho = 0;
        for (int j = start; j < n_dft; j++) mean_rho += rho_hist[j];
        mean_rho /= (n_dft - start);

        double var_rho = 0;
        for (int j = start; j < n_dft; j++) {
            double dr = rho_hist[j] - mean_rho;
            var_rho += dr * dr;
        }
        var_rho /= (n_dft - start);

        printf("  rho(center): mean=%.4e, variance=%.4e, relative var=%.4e\n",
               mean_rho, var_rho,
               (fabs(mean_rho) > 1e-30) ? var_rho / (mean_rho * mean_rho) : 0.0);
        printf("  Peak omega (rho_center): %.4f (power=%.4e)\n", omega_rho, peak_power_rho);
        printf("  Peak omega (phi_0_center): %.4f (power=%.4e)\n", omega_phi, peak_power_phi);

        int is_breathing = (omega_rho > 0.1 && var_rho / (mean_rho * mean_rho + 1e-30) > 0.01);
        printf("  Breathing? %s\n", is_breathing ? "YES" : "NO (static or dispersed)");

        /* Write DFT data */
        snprintf(path, sizeof(path), "%s/dynAC_dft.tsv", outdir);
        FILE *fdft = fopen(path, "w");
        if (fdft) {
            fprintf(fdft, "# DFT of rho(center,t) and phi_0(center,t) — second half\n");
            fprintf(fdft, "omega\tpower_rho\tpower_phi\n");
            double mean_phi0 = 0;
            for (int j = start; j < n_dft; j++) mean_phi0 += phi0_hist[j];
            mean_phi0 /= (n_dft - start);

            int nf = 500;
            for (int kk = 1; kk < nf; kk++) {
                double omega = 5.0 * kk / nf;
                double re_r = 0, im_r = 0, re_p = 0, im_p = 0;
                for (int j = start; j < n_dft; j++) {
                    double dtj = (j > start) ?
                        (t_hist[j]-t_hist[j-1]) : (t_hist[start+1]-t_hist[start]);
                    double c = cos(omega * t_hist[j]);
                    double s = sin(omega * t_hist[j]);
                    re_r += (rho_hist[j] - mean_rho) * c * dtj;
                    im_r += (rho_hist[j] - mean_rho) * s * dtj;
                    re_p += (phi0_hist[j] - mean_phi0) * c * dtj;
                    im_p += (phi0_hist[j] - mean_phi0) * s * dtj;
                }
                fprintf(fdft, "%.6f\t%.6e\t%.6e\n", omega, re_r*re_r+im_r*im_r, re_p*re_p+im_p*im_p);
            }
            fclose(fdft);
        }

        /* Write time histories */
        snprintf(path, sizeof(path), "%s/dynAC_rho_history.tsv", outdir);
        FILE *frho = fopen(path, "w");
        if (frho) {
            fprintf(frho, "time\trho_center\tphi0_center\n");
            for (int j = 0; j < n_dft; j++)
                fprintf(frho, "%.4f\t%.6e\t%.6e\n", t_hist[j], rho_hist[j], phi0_hist[j]);
            fclose(frho);
        }
    }

    /* Cleanup */
    for (int a = 0; a < 3; a++) {
        free(phi[a]); free(vel[a]); free(acc[a]);
    }
    free(damp);
    free(rho_hist); free(phi0_hist); free(t_hist);

    printf("\n=== V26-DynAC Complete ===\n");
    return 0;
}
