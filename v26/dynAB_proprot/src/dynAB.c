/*
 * dynAB.c — V26-DynAB: Propagating + Rotating Braid
 *
 * Combines axial propagation (A) with azimuthal rotation (B).
 * The braid both travels along z AND rotates around z, giving
 * both linear momentum Pz and angular momentum Lz.
 *
 * Initialization:
 *   phase_a = theta + k*z + 2*pi*a/3
 *   phi_a(x,0)   = A(r_perp) * cos(phase_a)
 *   vel_a(x,0)   = (omega + Omega) * A(r_perp) * sin(phase_a)
 *
 * where omega = sqrt(k^2 + m^2), k = 2*pi/Lz, Omega = rotation rate.
 *
 * Lagrangian (same as V26):
 *   L = Sum_a [1/2(dt phi_a)^2 - 1/2|grad phi_a|^2 - 1/2 m^2 phi_a^2]
 *     - V(P)      with V = mu*P^2/(2(1+kappa*P^2))
 *
 * Diagnostics: fc, |P|_max, l=2 multipole, Pz (linear mom), Lz (angular mom)
 *
 * Compile: gcc -O3 -fopenmp -Wall -o dynAB src/dynAB.c -lm
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>

/* ─── Parameters ─── */
static double mu_pot    = -20.0;
static double kappa     = 20.0;
static double mass      = 1.0;
static double A0        = 0.8;
static double Omega_rot = 0.1;     /* azimuthal rotation rate */
static double R_tube    = 3.0;     /* Gaussian envelope width */

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
        else if (!strcmp(argv[i], "-Omega"))    Omega_rot = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-Rtube"))    R_tube    = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-N"))        N         = atoi(argv[i+1]);
        else if (!strcmp(argv[i], "-L"))        L         = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-tfinal"))   tfinal    = atof(argv[i+1]);
        else if (!strcmp(argv[i], "-cfl"))      cfl_frac  = atof(argv[i+1]);
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

/* Periodic z BC */
static inline int wrap_z(int k)
{
    if (k < 0)   return k + N;
    if (k >= N)  return k - N;
    return k;
}

/* ─── Initialization: propagating + rotating braid ─── */
static void init_proprot(void)
{
    double Lz = 2.0 * L;           /* domain length in z */
    double k  = 2.0 * M_PI / Lz;   /* axial wavenumber */
    double omega = sqrt(k * k + mass * mass);  /* dispersion */

    printf("  Init: k=%.4f, omega=%.4f, Omega=%.4f, omega_eff=%.4f\n",
           k, omega, Omega_rot, omega + Omega_rot);

    #pragma omp parallel for schedule(static)
    for (long idx = 0; idx < Ngrid; idx++) {
        int ii = idx / (N * N);
        int jj = (idx / N) % N;
        int kk = idx % N;
        double x = -L + ii * dx;
        double y = -L + jj * dx;
        double z = -L + kk * dx;

        double rperp = sqrt(x * x + y * y);
        double envelope = A0 * exp(-rperp * rperp / (2.0 * R_tube * R_tube));
        double theta = atan2(y, x);

        for (int a = 0; a < 3; a++) {
            double phase = theta + k * z + 2.0 * M_PI * a / 3.0;
            phi[a][idx] = envelope * cos(phase);
            vel[a][idx] = (omega + Omega_rot) * envelope * sin(phase);
        }
    }
}

/* ─── Compute acceleration (triple-product potential only, periodic z) ─── */
static void compute_acc(void)
{
    for (int a = 0; a < 3; a++) {
        #pragma omp parallel for schedule(static)
        for (int i = 2; i < N-2; i++) {
            for (int j = 2; j < N-2; j++) {
                for (int k = 0; k < N; k++) {
                    long idx = IDX(i, j, k);

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

        /* Boundary: zero acceleration for x,y edges */
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
    /* Half-kick */
    for (int a = 0; a < 3; a++) {
        #pragma omp parallel for schedule(static)
        for (long idx = 0; idx < Ngrid; idx++)
            vel[a][idx] += 0.5 * dt * acc[a][idx];
    }

    /* Drift */
    for (int a = 0; a < 3; a++) {
        #pragma omp parallel for schedule(static)
        for (long idx = 0; idx < Ngrid; idx++)
            phi[a][idx] += dt * vel[a][idx];
    }

    /* Recompute acceleration */
    compute_acc();

    /* Half-kick */
    for (int a = 0; a < 3; a++) {
        #pragma omp parallel for schedule(static)
        for (long idx = 0; idx < Ngrid; idx++)
            vel[a][idx] += 0.5 * dt * acc[a][idx];
    }

    /* Absorbing boundary damping (transverse only, z is periodic) */
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
    double fc;
    double peak_P;
    double Pz;       /* linear momentum (z component) */
    double Lz;       /* angular momentum (z component) */
} Diag;

/* ─── Compute diagnostics ─── */
static Diag compute_diag(double core_radius)
{
    Diag d;
    memset(&d, 0, sizeof(d));

    double Ek = 0, Eg = 0, Em = 0, Ep = 0;
    double Ecore = 0, Eall = 0;
    double peak[3] = {0, 0, 0};
    double peak_P = 0;
    double Pz_tot = 0, Lz_tot = 0;

    #pragma omp parallel
    {
        double lEk = 0, lEg = 0, lEm = 0, lEp = 0, lEc = 0, lEa = 0;
        double lpk[3] = {0, 0, 0};
        double lpkP = 0;
        double lPz = 0, lLz = 0;

        #pragma omp for schedule(static) nowait
        for (int i = 2; i < N-2; i++) {
            double x = -L + i * dx;
            for (int j = 2; j < N-2; j++) {
                double y = -L + j * dx;
                for (int k = 0; k < N; k++) {
                    long idx = IDX(i, j, k);
                    double dV = dx * dx * dx;
                    double e_loc = 0;

                    int km1 = wrap_z(k - 1);
                    int kp1 = wrap_z(k + 1);

                    for (int a = 0; a < 3; a++) {
                        double v2 = vel[a][idx] * vel[a][idx];
                        lEk += 0.5 * v2 * dV;
                        e_loc += 0.5 * v2;

                        double gx = (phi[a][IDX(i+1,j,k)] - phi[a][IDX(i-1,j,k)]) / (2*dx);
                        double gy = (phi[a][IDX(i,j+1,k)] - phi[a][IDX(i,j-1,k)]) / (2*dx);
                        double gz = (phi[a][IDX(i,j,kp1)] - phi[a][IDX(i,j,km1)]) / (2*dx);
                        double grad2 = gx*gx + gy*gy + gz*gz;
                        lEg += 0.5 * grad2 * dV;
                        e_loc += 0.5 * grad2;

                        double mass_e = 0.5 * m2 * phi[a][idx] * phi[a][idx];
                        lEm += mass_e * dV;
                        e_loc += mass_e;

                        double ap = fabs(phi[a][idx]);
                        if (ap > lpk[a]) lpk[a] = ap;

                        /* Linear momentum: P_z^a = -dot(phi_a) * d_z(phi_a) */
                        lPz += -vel[a][idx] * gz * dV;

                        /* Angular momentum: L_z = sum_a integral (x * T^{0y}_a - y * T^{0x}_a) dV
                         *   T^{0j}_a = -dot(phi_a) * d_j(phi_a)
                         *   L_z = sum_a integral dot(phi_a) * (y * d_x - x * d_y) phi_a dV */
                        lLz += vel[a][idx] * (y * gx - x * gy) * dV;
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
                    double rperp = sqrt(x * x + y * y);
                    if (rperp < core_radius) lEc += e_loc * dV;
                }
            }
        }

        #pragma omp critical
        {
            Ek += lEk; Eg += lEg; Em += lEm; Ep += lEp;
            Ecore += lEc; Eall += lEa;
            Pz_tot += lPz; Lz_tot += lLz;
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
    d.Pz = Pz_tot;
    d.Lz = Lz_tot;

    return d;
}

/* ─── Legendre polynomials ─── */
static double legendre_P_fn(int l, double x)
{
    switch (l) {
        case 0: return 1.0;
        case 1: return x;
        case 2: return 0.5 * (3.0 * x * x - 1.0);
        case 3: return 0.5 * (5.0 * x * x * x - 3.0 * x);
        case 4: return (35.0*x*x*x*x - 30.0*x*x + 3.0) / 8.0;
        default: return 0.0;
    }
}

/* ─── Compute strain multipoles on cylindrical shell ─── */
static void compute_strain_multipoles(double R_shell, FILE *fmulti, double *c_l_out)
{
    /* For a z-periodic tube, sample on a cylinder r_perp = R_shell */
    int N_theta = 100, N_z_samp = 100;
    double c_l[5] = {0, 0, 0, 0, 0};  /* l = 0..4 */
    double norm = 0;
    int n_valid = 0;

    for (int it = 0; it < N_theta; it++) {
        double th = 2.0 * M_PI * it / N_theta;
        double xs = R_shell * cos(th);
        double ys = R_shell * sin(th);

        int i = (int)((xs + L) / dx + 0.5);
        int j = (int)((ys + L) / dx + 0.5);
        if (i < 2 || i >= N-2 || j < 2 || j >= N-2) continue;

        for (int iz = 0; iz < N_z_samp; iz++) {
            double zs = -L + 2.0 * L * iz / N_z_samp;
            int k = (int)((zs + L) / dx + 0.5);
            k = wrap_z(k);

            long idx = IDX(i, j, k);

            /* Compute energy density at this point */
            int km1 = wrap_z(k - 1);
            int kp1 = wrap_z(k + 1);
            double eloc = 0;
            for (int a = 0; a < 3; a++) {
                eloc += 0.5 * vel[a][idx] * vel[a][idx];
                double gx = (phi[a][IDX(i+1,j,k)] - phi[a][IDX(i-1,j,k)]) / (2*dx);
                double gy = (phi[a][IDX(i,j+1,k)] - phi[a][IDX(i,j-1,k)]) / (2*dx);
                double gz = (phi[a][IDX(i,j,kp1)] - phi[a][IDX(i,j,km1)]) / (2*dx);
                eloc += 0.5 * (gx*gx + gy*gy + gz*gz);
                eloc += 0.5 * m2 * phi[a][idx] * phi[a][idx];
            }
            double p0 = phi[0][idx], p1 = phi[1][idx], p2 = phi[2][idx];
            double P = p0 * p1 * p2;
            double P2 = P * P;
            eloc += 0.5 * mu_pot * P2 / (1.0 + kappa * P2);

            /* cos_angle for multipole: use theta angle on cylinder */
            double cos_th = cos(th);
            for (int l = 0; l <= 4; l++)
                c_l[l] += eloc * legendre_P_fn(l, cos_th);
            norm += eloc;
            n_valid++;
        }
    }

    if (fmulti) {
        fprintf(fmulti, "# Multipole decomposition of energy density on cylinder r=%.1f\n", R_shell);
        fprintf(fmulti, "l\tcoeff\tfraction\n");
        double sum_abs = 0;
        for (int l = 0; l <= 4; l++) sum_abs += fabs(c_l[l]);
        for (int l = 0; l <= 4; l++) {
            double frac = (sum_abs > 1e-30) ? fabs(c_l[l]) / sum_abs : 0.0;
            fprintf(fmulti, "%d\t%.6e\t%.6f\n", l, c_l[l], frac);
            printf("    l=%d: coefficient = %.6e, fraction = %.4f\n", l, c_l[l], frac);
        }
        printf("    Total energy on shell: %.6e (from %d points)\n",
               (n_valid > 0) ? norm / n_valid : 0.0, n_valid);
    }

    if (c_l_out) {
        for (int l = 0; l <= 4; l++) c_l_out[l] = c_l[l];
    }
}

/* ─── DFT to find peak frequency ─── */
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
                (t_hist[j] - t_hist[j-1]) : (t_hist[start+1] - t_hist[start]);
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

/* ─── Main ─── */
int main(int argc, char **argv)
{
    parse_args(argc, argv);

    dx  = 2.0 * L / (N - 1);
    dx2 = dx * dx;
    m2  = mass * mass;
    dt  = cfl_frac * dx;
    Ngrid = (long)N * N * N;

    printf("=== V26-DynAB: Propagating + Rotating Braid ===\n");
    printf("Parameters:\n");
    printf("  mu=%.1f  kappa=%.1f  mass=%.3f  A0=%.3f  Omega=%.3f\n",
           mu_pot, kappa, mass, A0, Omega_rot);
    printf("  R_tube=%.1f\n", R_tube);
    printf("  N=%d  L=%.1f  dx=%.4f  dt=%.5f\n", N, L, dx, dt);
    printf("  Ngrid=%ld (%.1f M)  Memory: %.1f MB\n",
           Ngrid, Ngrid/1e6, Ngrid * 8.0 * 10 / 1e6);
    printf("  Threads: %d\n", omp_get_max_threads());
    printf("  tfinal=%.0f  BC: periodic-z, absorbing-xy\n", tfinal);
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

    /* Damping profile: absorbing in transverse (r_perp), none in z */
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

    /* Initialize fields */
    init_proprot();

    double core_radius = 8.0;

    /* Time stepping */
    int Nt = (int)(tfinal / dt) + 1;

    /* DFT storage — lightweight center-only sampling */
    int max_dft = 5000;
    double *rho_hist = malloc(max_dft * sizeof(double));
    double *phi0_hist = malloc(max_dft * sizeof(double));
    double *t_hist   = malloc(max_dft * sizeof(double));
    int n_dft = 0;

    int diag_every = Nt / 200;
    if (diag_every < 1) diag_every = 1;
    int print_every = Nt / 30;
    if (print_every < 1) print_every = 1;
    int dft_every = Nt / max_dft;
    if (dft_every < 1) dft_every = 1;
    printf("  Nt=%d  diag_every=%d  dft_every=%d\n", Nt, diag_every, dft_every);

    /* Phase 1 output file */
    char path[600];
    snprintf(path, sizeof(path), "%s/dynAB_phase1.tsv", outdir);
    FILE *f1 = fopen(path, "w");
    if (!f1) { fprintf(stderr, "Cannot open %s\n", path); return 1; }
    fprintf(f1, "time\tE_total\tE_kin\tE_grad\tE_mass\tE_pot\t"
                "fc\tpeak0\tpeak1\tpeak2\tpeak_P\tPz\tLz\n");

    compute_acc();

    double wall_start = omp_get_wtime();

    /* Initial diagnostics */
    Diag d0 = compute_diag(core_radius);
    printf("\n--- Evolution (t=0..%.0f) ---\n", tfinal);
    printf("  t=%7.1f  E=%.2f  fc=%.4f  pk=(%.4f,%.4f,%.4f)  |P|=%.6f  Pz=%.4f  Lz=%.4f\n",
           0.0, d0.Et, d0.fc, d0.peak[0], d0.peak[1], d0.peak[2],
           d0.peak_P, d0.Pz, d0.Lz);
    fflush(stdout);

    for (int n = 0; n <= Nt; n++) {
        double t = n * dt;

        /* DFT history — lightweight: center density + phi0 only */
        if (n % dft_every == 0 && n_dft < max_dft) {
            long ic = IDX(N/2, N/2, N/2);
            int ci = N/2, cj = N/2, ck = N/2;
            double rho_c = 0;
            for (int a = 0; a < 3; a++) {
                rho_c += 0.5 * vel[a][ic] * vel[a][ic];
                double gx = (phi[a][IDX(ci+1,cj,ck)] - phi[a][IDX(ci-1,cj,ck)]) / (2*dx);
                double gy = (phi[a][IDX(ci,cj+1,ck)] - phi[a][IDX(ci,cj-1,ck)]) / (2*dx);
                int ckm = wrap_z(ck-1), ckp = wrap_z(ck+1);
                double gz = (phi[a][IDX(ci,cj,ckp)] - phi[a][IDX(ci,cj,ckm)]) / (2*dx);
                rho_c += 0.5 * (gx*gx + gy*gy + gz*gz);
                rho_c += 0.5 * m2 * phi[a][ic] * phi[a][ic];
            }
            rho_hist[n_dft]  = rho_c;
            phi0_hist[n_dft] = phi[0][ic];
            t_hist[n_dft]    = t;
            n_dft++;
        }

        int do_diag  = (n % diag_every == 0);
        int do_print = (n % print_every == 0);

        if (do_diag || do_print) {
            Diag d = compute_diag(core_radius);

            if (do_diag)
                fprintf(f1, "%.4f\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t"
                            "%.6f\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\n",
                        t, d.Et, d.Ek, d.Eg, d.Em, d.Ep,
                        d.fc, d.peak[0], d.peak[1], d.peak[2], d.peak_P, d.Pz, d.Lz);

            if (do_print) {
                double elapsed = omp_get_wtime() - wall_start;
                double frac = (double)n / Nt;
                double eta_t = (frac > 0.001) ? elapsed * (1.0 - frac) / frac : 0;
                printf("  t=%7.1f  E=%.2f  fc=%.4f  pk=(%.4f,%.4f,%.4f)  |P|=%.6f  Pz=%.2f  Lz=%.2f  [%.0fs, ETA %.0fs]\n",
                       t, d.Et, d.fc, d.peak[0], d.peak[1], d.peak[2],
                       d.peak_P, d.Pz, d.Lz, elapsed, eta_t);
                fflush(stdout);
            }
        }

        if (n == Nt) break;
        verlet_step();
    }

    fclose(f1);

    /* Final diagnostics */
    Diag dfinal = compute_diag(core_radius);
    double elapsed1 = omp_get_wtime() - wall_start;

    printf("\nEvolution complete (%.1f sec)\n", elapsed1);
    printf("  Final: E=%.2f  fc=%.4f  pk=(%.4f,%.4f,%.4f)  |P|max=%.6f\n",
           dfinal.Et, dfinal.fc, dfinal.peak[0], dfinal.peak[1], dfinal.peak[2], dfinal.peak_P);
    printf("  Final: Pz=%.4f  Lz=%.4f\n", dfinal.Pz, dfinal.Lz);

    int survived = (dfinal.fc > 0.01 && dfinal.peak_P > 1e-6);
    printf("  Survived? %s (fc=%.4f, |P|max=%.6f)\n",
           survived ? "YES" : "NO", dfinal.fc, dfinal.peak_P);

    /* Phase 2: DFT of center density */
    printf("\n--- DFT Analysis ---\n");

    if (n_dft >= 100) {
        int start = n_dft / 2;

        double pp_rho = 0, pp_phi = 0;
        double omega_rho = find_peak_omega(rho_hist,  t_hist, n_dft, start, 5.0, &pp_rho);
        double omega_phi = find_peak_omega(phi0_hist, t_hist, n_dft, start, 5.0, &pp_phi);

        printf("  Peak omega (rho_center): %.4f (power=%.4e)\n", omega_rho, pp_rho);
        printf("  Peak omega (phi0_center): %.4f (power=%.4e)\n", omega_phi, pp_phi);

        /* Variance of rho */
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

        int is_breathing = (omega_rho > 0.1 && var_rho / (mean_rho * mean_rho + 1e-30) > 0.01);
        printf("  Breathing? %s\n", is_breathing ? "YES" : "NO (static or dispersed)");

        /* Write history */
        snprintf(path, sizeof(path), "%s/dynAB_rho_history.tsv", outdir);
        FILE *fh = fopen(path, "w");
        if (fh) {
            fprintf(fh, "time\trho_center\tphi0_center\n");
            for (int j = 0; j < n_dft; j++)
                fprintf(fh, "%.4f\t%.6e\t%.6e\n",
                        t_hist[j], rho_hist[j], phi0_hist[j]);
            fclose(fh);
        }
    }

    /* Phase 3: Strain multipoles on cylindrical shell */
    printf("\n--- Multipole Analysis ---\n");

    double R_shell = 6.0;
    snprintf(path, sizeof(path), "%s/dynAB_multipoles.tsv", outdir);
    FILE *fmulti = fopen(path, "w");
    double c_l[5] = {0};
    if (fmulti) {
        compute_strain_multipoles(R_shell, fmulti, c_l);
        fclose(fmulti);
    }

    /* Summary */
    printf("\n=== SUMMARY ===\n");
    printf("  Configuration: Propagating + Rotating Braid\n");
    printf("  mu=%.1f  kappa=%.1f  mass=%.3f  Omega=%.3f\n",
           mu_pot, kappa, mass, Omega_rot);
    printf("  N=%d  L=%.1f  tfinal=%.0f\n", N, L, tfinal);
    printf("  Final E=%.4f  fc=%.4f  |P|max=%.6f\n",
           dfinal.Et, dfinal.fc, dfinal.peak_P);
    printf("  Final Pz=%.4f  Lz=%.4f\n", dfinal.Pz, dfinal.Lz);
    printf("  Survived: %s\n", survived ? "YES" : "NO");
    double l2_frac = 0;
    {
        double sum_abs = 0;
        for (int l = 0; l <= 4; l++) sum_abs += fabs(c_l[l]);
        l2_frac = (sum_abs > 1e-30) ? fabs(c_l[2]) / sum_abs : 0.0;
    }
    printf("  l=2 fraction: %.4f\n", l2_frac);
    printf("  Wall time: %.1f sec\n", elapsed1);

    /* Cleanup */
    for (int a = 0; a < 3; a++) {
        free(phi[a]); free(vel[a]); free(acc[a]);
    }
    free(damp);
    free(rho_hist); free(phi0_hist); free(t_hist);

    printf("\n=== V26-DynAB Complete ===\n");
    return 0;
}
