/*  gen_galerkin_seed.c — Relaxation-based equilibrium seed generator
 *
 *  Loads a proton template SFA file and iteratively relaxes the fields
 *  toward the static equilibrium of the Cosserat equations:
 *
 *    ∇²φ - m²φ - V'(P) + η·curl(θ) = 0
 *    ∇²θ - m_θ²θ + η·curl(φ) = 0
 *
 *  The method is overdamped dynamics (gradient descent on the potential
 *  energy): at each step, compute the force (PDE residual) and update
 *  φ += α * force. No velocities are evolved — this is a pure relaxation.
 *
 *  KEY DESIGN: The proton exists as a perturbation on top of a background
 *  carrier wave. A spherical mask limits updates to the proton's core
 *  region, and a soft taper smoothly transitions to the pinned background
 *  at the boundary. This prevents the trivial solution (zero fields) from
 *  dominating the relaxation.
 *
 *  Build: gcc -O3 -march=native -fopenmp -o gen_galerkin_seed gen_galerkin_seed.c -lzstd -lm
 *
 *  Usage:
 *    ./gen_galerkin_seed input.sfa output.sfa [--iterations 1000] [--alpha 0.01]
 *    ./gen_galerkin_seed input.sfa output.sfa --iterations 500 --alpha 0.02 --mask_r 6.0
 */

#define SFA_IMPLEMENTATION
#include "../format/sfa.h"
#include "../sim/scp_config.h"

#include <omp.h>

#define NCOLS 12  /* 3 phi + 3 theta + 3 phi_vel + 3 theta_vel */

/* ================================================================
   Grid: 6 field arrays (phi[3], theta[3]) in double precision
   ================================================================ */

typedef struct {
    double *phi[NFIELDS];
    double *theta[NFIELDS];
    double *phi_vel[NFIELDS];    /* original velocities from template */
    double *theta_vel[NFIELDS];
    double *mask;     /* spatial relaxation mask: 1 at center, 0 at boundary */
    int N;
    long N3;
    double L, dx;
    int has_velocities;
} RelaxGrid;

static RelaxGrid *rgrid_alloc(int N, double L) {
    RelaxGrid *g = calloc(1, sizeof(RelaxGrid));
    g->N = N;
    g->N3 = (long)N * N * N;
    g->L = L;
    g->dx = 2.0 * L / (N - 1);
    for (int a = 0; a < NFIELDS; a++) {
        g->phi[a] = calloc(g->N3, sizeof(double));
        g->theta[a] = calloc(g->N3, sizeof(double));
        g->phi_vel[a] = calloc(g->N3, sizeof(double));
        g->theta_vel[a] = calloc(g->N3, sizeof(double));
    }
    g->mask = calloc(g->N3, sizeof(double));
    g->has_velocities = 0;
    return g;
}

static void rgrid_free(RelaxGrid *g) {
    for (int a = 0; a < NFIELDS; a++) {
        free(g->phi[a]); free(g->theta[a]);
        free(g->phi_vel[a]); free(g->theta_vel[a]);
    }
    free(g->mask);
    free(g);
}

/* Initialize the spatial mask: 1 inside R_inner, tapers to 0 at R_outer.
 * Uses a smooth cos² taper in the transition zone. */
static void init_mask(RelaxGrid *g, double R_inner, double R_taper) {
    const int N = g->N, NN = N * N;
    const double L = g->L, dx = g->dx;
    double R_outer = R_inner + R_taper;

    #pragma omp parallel for schedule(static)
    for (long idx = 0; idx < g->N3; idx++) {
        int i = (int)(idx / NN), j = (int)((idx / N) % N), k = (int)(idx % N);
        double x = -L + i * dx;
        double y = -L + j * dx;
        double z = -L + k * dx;
        double r = sqrt(x*x + y*y + z*z);
        if (r <= R_inner)
            g->mask[idx] = 1.0;
        else if (r >= R_outer)
            g->mask[idx] = 0.0;
        else {
            double s = (r - R_inner) / R_taper;
            g->mask[idx] = 0.5 * (1.0 + cos(PI * s));  /* cos² taper */
        }
    }
}

/* ================================================================
   Load template SFA into RelaxGrid (fields only, ignore velocities)
   ================================================================ */

static int load_template(RelaxGrid *g, const char *path) {
    SFA *sfa = sfa_open(path);
    if (!sfa) { fprintf(stderr, "FATAL: cannot open SFA '%s'\n", path); return -1; }

    if ((int)sfa->Nx != g->N || (int)sfa->Ny != g->N || (int)sfa->Nz != g->N) {
        fprintf(stderr, "FATAL: SFA grid %ux%ux%u != expected %d^3\n",
                sfa->Nx, sfa->Ny, sfa->Nz, g->N);
        sfa_close(sfa);
        return -1;
    }

    int frame = sfa->total_frames - 1;
    printf("Loading '%s': %ux%ux%u, %u frames, reading frame %d\n",
           path, sfa->Nx, sfa->Ny, sfa->Nz, sfa->total_frames, frame);

    void *buf = malloc(sfa->frame_bytes);
    if (!buf) { fprintf(stderr, "FATAL: malloc failed\n"); sfa_close(sfa); return -1; }
    sfa_read_frame(sfa, frame, buf);

    /* Walk columns and load fields */
    int loaded = 0;
    uint64_t off = 0;
    for (uint32_t col = 0; col < sfa->n_columns; col++) {
        int dtype = sfa->columns[col].dtype;
        int sem = sfa->columns[col].semantic;
        int comp = sfa->columns[col].component;
        int es = sfa_dtype_size[dtype];
        uint8_t *src = (uint8_t*)buf + off;

        double *target = NULL;
        if (sem == SFA_POSITION && comp < 3)           target = g->phi[comp];
        else if (sem == SFA_ANGLE && comp < 3)          target = g->theta[comp];
        else if (sem == SFA_VELOCITY && comp < 3)       target = g->phi_vel[comp];
        else if (sem == SFA_VELOCITY && comp >= 3 && comp < 6) target = g->theta_vel[comp - 3];

        if (target) {
            long N3 = g->N3;
            if (dtype == SFA_F64)
                for (long i = 0; i < N3; i++) target[i] = ((double*)src)[i];
            else if (dtype == SFA_F32)
                for (long i = 0; i < N3; i++) target[i] = (double)((float*)src)[i];
            else if (dtype == SFA_F16)
                for (long i = 0; i < N3; i++) target[i] = f16_to_f64(((uint16_t*)src)[i]);
            loaded++;
            if (sem == SFA_VELOCITY) g->has_velocities = 1;
        }
        off += (uint64_t)g->N3 * es;
    }

    printf("  Loaded %d/6 field arrays (phi + theta)\n", loaded);
    free(buf);
    sfa_close(sfa);
    return (loaded >= 3) ? 0 : -1;  /* need at least phi */
}

/* ================================================================
   Compute equilibrium residual (force) at every voxel
   R_phi[a] = ∇²φ_a - m²φ_a - dV/dP · dP/dφ_a + η·curl(θ)_a
   R_theta[a] = ∇²θ_a - m_θ²θ_a + η·curl(φ)_a
   ================================================================ */

static void compute_residual(RelaxGrid *g,
    double m2, double mtheta2, double mu, double kappa, double eta,
    double *R_phi[NFIELDS], double *R_theta[NFIELDS]) {

    const int N = g->N, NN = N * N;
    const long N3 = g->N3;
    const double idx2 = 1.0 / (g->dx * g->dx);
    const double idx1 = 1.0 / (2.0 * g->dx);

    #pragma omp parallel for schedule(static)
    for (long idx = 0; idx < N3; idx++) {
        int i = (int)(idx / NN), j = (int)((idx / N) % N), k = (int)(idx % N);
        /* Periodic boundary conditions */
        int ip = (i+1) % N, im = (i-1+N) % N;
        int jp = (j+1) % N, jm = (j-1+N) % N;
        int kp = (k+1) % N, km = (k-1+N) % N;
        long n_ip = (long)ip*NN + j*N + k;
        long n_im = (long)im*NN + j*N + k;
        long n_jp = (long)i*NN + jp*N + k;
        long n_jm = (long)i*NN + jm*N + k;
        long n_kp = (long)i*NN + j*N + kp;
        long n_km = (long)i*NN + j*N + km;

        double p0 = g->phi[0][idx], p1 = g->phi[1][idx], p2 = g->phi[2][idx];
        double P = p0 * p1 * p2;
        double P2 = P * P;
        double den = 1.0 + kappa * P2;
        double dVdP = mu * P / (den * den);

        /* Phi forces */
        for (int a = 0; a < NFIELDS; a++) {
            double lap = (g->phi[a][n_ip] + g->phi[a][n_im]
                        + g->phi[a][n_jp] + g->phi[a][n_jm]
                        + g->phi[a][n_kp] + g->phi[a][n_km]
                        - 6.0 * g->phi[a][idx]) * idx2;

            double dPda = (a == 0) ? p1*p2 : (a == 1) ? p0*p2 : p0*p1;

            /* curl(θ)_a */
            double ct;
            if (a == 0)
                ct = (g->theta[2][n_jp] - g->theta[2][n_jm]
                    - g->theta[1][n_kp] + g->theta[1][n_km]) * idx1;
            else if (a == 1)
                ct = (g->theta[0][n_kp] - g->theta[0][n_km]
                    - g->theta[2][n_ip] + g->theta[2][n_im]) * idx1;
            else
                ct = (g->theta[1][n_ip] - g->theta[1][n_im]
                    - g->theta[0][n_jp] + g->theta[0][n_jm]) * idx1;

            R_phi[a][idx] = lap - m2 * g->phi[a][idx] - dVdP * dPda + eta * ct;
        }

        /* Theta forces */
        for (int a = 0; a < NFIELDS; a++) {
            double lapt = (g->theta[a][n_ip] + g->theta[a][n_im]
                         + g->theta[a][n_jp] + g->theta[a][n_jm]
                         + g->theta[a][n_kp] + g->theta[a][n_km]
                         - 6.0 * g->theta[a][idx]) * idx2;

            /* curl(φ)_a */
            double cp;
            if (a == 0)
                cp = (g->phi[2][n_jp] - g->phi[2][n_jm]
                    - g->phi[1][n_kp] + g->phi[1][n_km]) * idx1;
            else if (a == 1)
                cp = (g->phi[0][n_kp] - g->phi[0][n_km]
                    - g->phi[2][n_ip] + g->phi[2][n_im]) * idx1;
            else
                cp = (g->phi[1][n_ip] - g->phi[1][n_im]
                    - g->phi[0][n_jp] + g->phi[0][n_jm]) * idx1;

            R_theta[a][idx] = lapt - mtheta2 * g->theta[a][idx] + eta * cp;
        }
    }
}

/* ================================================================
   Compute diagnostics: RMS residual, max|phi|, integrated |P|
   ================================================================ */

typedef struct {
    double rms_phi;
    double rms_theta;
    double phi_max;
    double P_int;
    double theta_rms;
    /* Masked versions (only where mask > 0.5) */
    double rms_phi_core;
    double rms_theta_core;
} DiagInfo;

static DiagInfo compute_diag(RelaxGrid *g,
    double *R_phi[NFIELDS], double *R_theta[NFIELDS]) {

    const long N3 = g->N3;
    const double dV = g->dx * g->dx * g->dx;
    double sum_rp = 0, sum_rt = 0, pm = 0, pi = 0, trms = 0;
    double core_rp = 0, core_rt = 0;
    long core_count = 0;

    #pragma omp parallel for reduction(+:sum_rp,sum_rt,pi,trms,core_rp,core_rt,core_count) reduction(max:pm)
    for (long idx = 0; idx < N3; idx++) {
        int is_core = (g->mask[idx] > 0.5) ? 1 : 0;
        for (int a = 0; a < NFIELDS; a++) {
            double rp2 = R_phi[a][idx] * R_phi[a][idx];
            double rt2 = R_theta[a][idx] * R_theta[a][idx];
            sum_rp += rp2;
            sum_rt += rt2;
            if (is_core) { core_rp += rp2; core_rt += rt2; }
            double ap = fabs(g->phi[a][idx]);
            if (ap > pm) pm = ap;
            trms += g->theta[a][idx] * g->theta[a][idx];
        }
        if (is_core) core_count++;
        double P = fabs(g->phi[0][idx] * g->phi[1][idx] * g->phi[2][idx]);
        pi += P * dV;
    }
    DiagInfo d;
    d.rms_phi = sqrt(sum_rp / (3.0 * N3));
    d.rms_theta = sqrt(sum_rt / (3.0 * N3));
    d.phi_max = pm;
    d.P_int = pi;
    d.theta_rms = sqrt(trms / (3.0 * N3));
    d.rms_phi_core = (core_count > 0) ? sqrt(core_rp / (3.0 * core_count)) : 0;
    d.rms_theta_core = (core_count > 0) ? sqrt(core_rt / (3.0 * core_count)) : 0;
    return d;
}

/* ================================================================
   Compute total potential energy
   ================================================================ */

static double compute_energy(RelaxGrid *g,
    double m2, double mtheta2, double mu, double kappa, double eta) {

    const int N = g->N, NN = N * N;
    const long N3 = g->N3;
    const double dx = g->dx, dV = dx * dx * dx;
    const double idx1 = 1.0 / (2.0 * dx);
    double E = 0;

    #pragma omp parallel for reduction(+:E)
    for (long idx = 0; idx < N3; idx++) {
        int i = (int)(idx / NN), j = (int)((idx / N) % N), k = (int)(idx % N);
        int ip = (i+1) % N, im = (i-1+N) % N;
        int jp = (j+1) % N, jm = (j-1+N) % N;
        int kp = (k+1) % N, km = (k-1+N) % N;
        long n_ip = (long)ip*NN + j*N + k, n_im = (long)im*NN + j*N + k;
        long n_jp = (long)i*NN + jp*N + k, n_jm = (long)i*NN + jm*N + k;
        long n_kp = (long)i*NN + j*N + kp, n_km = (long)i*NN + j*N + km;

        double p0 = g->phi[0][idx], p1 = g->phi[1][idx], p2 = g->phi[2][idx];
        double P = p0 * p1 * p2, P2 = P * P;

        for (int a = 0; a < NFIELDS; a++) {
            double gx = (g->phi[a][n_ip] - g->phi[a][n_im]) * idx1;
            double gy = (g->phi[a][n_jp] - g->phi[a][n_jm]) * idx1;
            double gz = (g->phi[a][n_kp] - g->phi[a][n_km]) * idx1;
            E += 0.5 * (gx*gx + gy*gy + gz*gz) * dV;
            E += 0.5 * m2 * g->phi[a][idx] * g->phi[a][idx] * dV;
            double tgx = (g->theta[a][n_ip] - g->theta[a][n_im]) * idx1;
            double tgy = (g->theta[a][n_jp] - g->theta[a][n_jm]) * idx1;
            double tgz = (g->theta[a][n_kp] - g->theta[a][n_km]) * idx1;
            E += 0.5 * (tgx*tgx + tgy*tgy + tgz*tgz) * dV;
            E += 0.5 * mtheta2 * g->theta[a][idx] * g->theta[a][idx] * dV;
        }
        E += (mu / 2.0) * P2 / (1.0 + kappa * P2) * dV;

        for (int a = 0; a < NFIELDS; a++) {
            double ct;
            if (a == 0)
                ct = (g->theta[2][n_jp] - g->theta[2][n_jm]
                    - g->theta[1][n_kp] + g->theta[1][n_km]) * idx1;
            else if (a == 1)
                ct = (g->theta[0][n_kp] - g->theta[0][n_km]
                    - g->theta[2][n_ip] + g->theta[2][n_im]) * idx1;
            else
                ct = (g->theta[1][n_ip] - g->theta[1][n_im]
                    - g->theta[0][n_jp] + g->theta[0][n_jm]) * idx1;
            E -= eta * g->phi[a][idx] * ct * dV;
        }
    }
    return E;
}

/* ================================================================
   Write output SFA: 12 columns (fields + background velocities)
   The background carrier wave has velocity omega*A_bg*sin(k*z + 2*pi*a/3).
   Without this, cold-starting from zero velocity creates a massive
   transient that destroys the proton.
   ================================================================ */

static void write_seed(RelaxGrid *g, const char *path, int precision,
                       double m2, double A_bg) {
    long N3 = g->N3;
    int N = g->N;
    double L = g->L;
    double dx = g->dx;
    double dt = 0.025 * dx;

    /* Use original template velocities if available, else generate background */
    double *phi_vel[NFIELDS];
    double *theta_vel[NFIELDS];
    int vel_allocated = 0;

    if (g->has_velocities) {
        printf("  Using original template velocities\n");
        for (int a = 0; a < NFIELDS; a++) {
            phi_vel[a] = g->phi_vel[a];
            theta_vel[a] = g->theta_vel[a];
        }
    } else {
        printf("  Generating background velocities (A_bg=%.3f)\n", A_bg);
        double k_bg = PI / L;
        double omega_bg = sqrt(k_bg * k_bg + m2);
        vel_allocated = 1;
        for (int a = 0; a < NFIELDS; a++) {
            phi_vel[a] = calloc(N3, sizeof(double));
            theta_vel[a] = calloc(N3, sizeof(double));
            #pragma omp parallel for schedule(static)
            for (long idx = 0; idx < N3; idx++) {
                int k = (int)(idx % N);
                double z = -L + k * dx;
                double phase = k_bg * z + 2.0 * PI * a / 3.0;
                phi_vel[a][idx] = omega_bg * A_bg * sin(phase);
            }
        }
    }

    uint8_t sfa_dtype = (precision == 0) ? SFA_F16 :
                        (precision == 1) ? SFA_F32 : SFA_F64;
    SFA *sfa = sfa_create(path, N, N, N, L, L, L, dt);

    const char *keys[] = {"generator", "N", "L"};
    char buf_N[32], buf_L[32];
    snprintf(buf_N, 32, "%d", N);
    snprintf(buf_L, 32, "%.6f", L);
    const char *vals[] = {"gen_galerkin_seed", buf_N, buf_L};
    sfa_add_kvmd(sfa, 0, 0xFFFFFFFF, 0xFFFFFFFF, keys, vals, 3);

    sfa_add_column(sfa, "phi_x",    sfa_dtype, SFA_POSITION, 0);
    sfa_add_column(sfa, "phi_y",    sfa_dtype, SFA_POSITION, 1);
    sfa_add_column(sfa, "phi_z",    sfa_dtype, SFA_POSITION, 2);
    sfa_add_column(sfa, "theta_x",  sfa_dtype, SFA_ANGLE,    0);
    sfa_add_column(sfa, "theta_y",  sfa_dtype, SFA_ANGLE,    1);
    sfa_add_column(sfa, "theta_z",  sfa_dtype, SFA_ANGLE,    2);
    sfa_add_column(sfa, "phi_vx",   sfa_dtype, SFA_VELOCITY, 0);
    sfa_add_column(sfa, "phi_vy",   sfa_dtype, SFA_VELOCITY, 1);
    sfa_add_column(sfa, "phi_vz",   sfa_dtype, SFA_VELOCITY, 2);
    sfa_add_column(sfa, "theta_vx", sfa_dtype, SFA_VELOCITY, 3);
    sfa_add_column(sfa, "theta_vy", sfa_dtype, SFA_VELOCITY, 4);
    sfa_add_column(sfa, "theta_vz", sfa_dtype, SFA_VELOCITY, 5);
    sfa_finalize_header(sfa);

    if (precision == 2) {
        void *cols[12] = {
            g->phi[0], g->phi[1], g->phi[2],
            g->theta[0], g->theta[1], g->theta[2],
            phi_vel[0], phi_vel[1], phi_vel[2],
            theta_vel[0], theta_vel[1], theta_vel[2]
        };
        sfa_write_frame(sfa, 0.0, cols);
    } else {
        void *cols[12];
        int es = (precision == 0) ? 2 : 4;
        for (int c = 0; c < 12; c++) {
            double *src;
            if (c < 3) src = g->phi[c];
            else if (c < 6) src = g->theta[c - 3];
            else if (c < 9) src = phi_vel[c - 6];
            else src = theta_vel[c - 9];

            cols[c] = malloc(N3 * es);
            if (precision == 1) {
                float *p = (float*)cols[c];
                for (long i = 0; i < N3; i++) p[i] = (float)src[i];
            } else {
                uint16_t *p = (uint16_t*)cols[c];
                for (long i = 0; i < N3; i++) p[i] = f64_to_f16(src[i]);
            }
        }
        sfa_write_frame(sfa, 0.0, cols);
        for (int c = 0; c < 12; c++) free(cols[c]);
    }
    if (vel_allocated) {
        for (int a = 0; a < NFIELDS; a++) { free(phi_vel[a]); free(theta_vel[a]); }
    }
    sfa_close(sfa);
}

/* ================================================================
   Main
   ================================================================ */

int main(int argc, char **argv) {
    if (argc < 3) {
        fprintf(stderr,
            "Usage: %s input.sfa output.sfa [options]\n\n"
            "Options:\n"
            "  --iterations N   Number of relaxation iterations (default 1000)\n"
            "  --alpha STEP     Step size for gradient descent (default 0.01)\n"
            "  --mask_r RADIUS  Inner radius of relaxation mask (default 5.0)\n"
            "  --mask_taper W   Width of mask taper zone (default 2.0)\n"
            "  --precision P    Output precision: f16, f32, f64 (default f32)\n"
            "  --m2 VAL         Mass squared (default 2.25)\n"
            "  --mu VAL         Potential coupling (default -41.345)\n"
            "  --kappa VAL      Potential saturation (default 50)\n"
            "  --eta VAL        Curl coupling (default 0.5)\n"
            "  --mtheta2 VAL    Theta mass squared (default 0)\n"
            "  --A_bg VAL       Background amplitude for velocities (default 0.1)\n"
            "  --diag FILE      Write diagnostics to TSV file\n"
            "\n", argv[0]);
        return 1;
    }

    const char *input = argv[1];
    const char *output = argv[2];
    int iterations = 1000;
    double alpha = 0.01;
    double mask_r = 5.0;
    double mask_taper = 2.0;
    int precision = 1;  /* f32 */
    double m2 = 2.25, mu = -41.345, kappa = 50.0, eta = 0.5, mtheta2 = 0.0;
    double A_bg = 0.1;
    const char *diag_file = NULL;

    for (int i = 3; i < argc - 1; i++) {
        if (!strcmp(argv[i], "--iterations"))  iterations = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--alpha"))  alpha = atof(argv[++i]);
        else if (!strcmp(argv[i], "--mask_r")) mask_r = atof(argv[++i]);
        else if (!strcmp(argv[i], "--mask_taper")) mask_taper = atof(argv[++i]);
        else if (!strcmp(argv[i], "--m2"))     m2 = atof(argv[++i]);
        else if (!strcmp(argv[i], "--mu"))     mu = atof(argv[++i]);
        else if (!strcmp(argv[i], "--kappa"))  kappa = atof(argv[++i]);
        else if (!strcmp(argv[i], "--eta"))    eta = atof(argv[++i]);
        else if (!strcmp(argv[i], "--mtheta2")) mtheta2 = atof(argv[++i]);
        else if (!strcmp(argv[i], "--A_bg"))   A_bg = atof(argv[++i]);
        else if (!strcmp(argv[i], "--diag"))   diag_file = argv[++i];
        else if (!strcmp(argv[i], "--precision")) {
            i++;
            if (!strcmp(argv[i], "f16")) precision = 0;
            else if (!strcmp(argv[i], "f32")) precision = 1;
            else if (!strcmp(argv[i], "f64")) precision = 2;
        }
        else fprintf(stderr, "WARNING: unknown option '%s'\n", argv[i]);
    }

    /* OpenMP setup */
    int nth = 4;
    char *env_t = getenv("OMP_NUM_THREADS");
    if (env_t) nth = atoi(env_t);
    omp_set_num_threads(nth);

    printf("gen_galerkin_seed: Cosserat equilibrium relaxation\n");
    printf("  Input:  %s\n", input);
    printf("  Output: %s\n", output);
    printf("  Iterations: %d, alpha: %.6f\n", iterations, alpha);
    printf("  Mask: R_inner=%.2f, taper=%.2f (R_outer=%.2f)\n",
           mask_r, mask_taper, mask_r + mask_taper);
    printf("  Physics: m^2=%.4f mu=%.3f kappa=%.1f eta=%.3f m_theta^2=%.4f\n",
           m2, mu, kappa, eta, mtheta2);
    printf("  Threads: %d\n", nth);

    /* Open the input to get grid size */
    SFA *probe = sfa_open(input);
    if (!probe) { fprintf(stderr, "FATAL: cannot open '%s'\n", input); return 1; }
    int N = probe->Nx;
    double L = probe->Lx;
    sfa_close(probe);

    double dx = 2.0 * L / (N - 1);
    printf("  Grid: N=%d, L=%.2f, dx=%.4f\n\n", N, L, dx);

    /* Allocate grid and load template */
    RelaxGrid *g = rgrid_alloc(N, L);
    if (load_template(g, input) < 0) { rgrid_free(g); return 1; }

    /* Initialize spatial mask */
    init_mask(g, mask_r, mask_taper);

    /* Count masked voxels */
    long mask_count = 0;
    for (long i = 0; i < g->N3; i++)
        if (g->mask[i] > 0.01) mask_count++;
    printf("Mask: %ld/%ld voxels active (%.1f%%)\n\n",
           mask_count, g->N3, 100.0 * mask_count / g->N3);

    /* Save the initial (template) state for pinning outside mask */
    double *pin_phi[NFIELDS], *pin_theta[NFIELDS];
    for (int a = 0; a < NFIELDS; a++) {
        pin_phi[a] = malloc(g->N3 * sizeof(double));
        pin_theta[a] = malloc(g->N3 * sizeof(double));
        memcpy(pin_phi[a], g->phi[a], g->N3 * sizeof(double));
        memcpy(pin_theta[a], g->theta[a], g->N3 * sizeof(double));
    }

    /* Residual arrays */
    double *R_phi[NFIELDS], *R_theta[NFIELDS];
    for (int a = 0; a < NFIELDS; a++) {
        R_phi[a] = malloc(g->N3 * sizeof(double));
        R_theta[a] = malloc(g->N3 * sizeof(double));
    }

    /* Open diagnostics file if requested */
    FILE *diag_fp = NULL;
    if (diag_file) {
        diag_fp = fopen(diag_file, "w");
        if (diag_fp)
            fprintf(diag_fp, "iter\trms_phi\trms_theta\trms_phi_core\trms_theta_core\tenergy\tP_int\tphi_max\ttheta_rms\n");
    }

    /* Initial diagnostics */
    compute_residual(g, m2, mtheta2, mu, kappa, eta, R_phi, R_theta);
    DiagInfo d0 = compute_diag(g, R_phi, R_theta);
    double E0 = compute_energy(g, m2, mtheta2, mu, kappa, eta);
    printf("Initial: rms_phi=%.6e rms_theta=%.6e E=%.6f P_int=%.4f phi_max=%.4f\n",
           d0.rms_phi, d0.rms_theta, E0, d0.P_int, d0.phi_max);
    printf("  Core:  rms_phi=%.6e rms_theta=%.6e\n\n", d0.rms_phi_core, d0.rms_theta_core);

    if (diag_fp)
        fprintf(diag_fp, "0\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.6f\t%.6f\t%.6e\n",
                d0.rms_phi, d0.rms_theta, d0.rms_phi_core, d0.rms_theta_core,
                E0, d0.P_int, d0.phi_max, d0.theta_rms);

    /* ================================================================
       Relaxation loop
       ================================================================ */

    double prev_E = E0;

    for (int iter = 1; iter <= iterations; iter++) {
        /* Compute residual */
        compute_residual(g, m2, mtheta2, mu, kappa, eta, R_phi, R_theta);

        /* Update: fields += alpha * mask * residual
         * Then blend with pinned values outside mask */
        for (int a = 0; a < NFIELDS; a++) {
            #pragma omp parallel for schedule(static)
            for (long i = 0; i < g->N3; i++) {
                double m = g->mask[i];
                if (m > 0.0) {
                    /* Apply masked update */
                    g->phi[a][i] += alpha * m * R_phi[a][i];
                    g->theta[a][i] += alpha * m * R_theta[a][i];
                }
                /* Blend back toward pinned state where mask < 1 */
                if (m < 1.0) {
                    g->phi[a][i] = m * g->phi[a][i] + (1.0 - m) * pin_phi[a][i];
                    g->theta[a][i] = m * g->theta[a][i] + (1.0 - m) * pin_theta[a][i];
                }
            }
        }

        /* Periodic diagnostics */
        if (iter % 100 == 0 || iter == iterations || iter <= 10) {
            compute_residual(g, m2, mtheta2, mu, kappa, eta, R_phi, R_theta);
            DiagInfo d = compute_diag(g, R_phi, R_theta);
            double E = compute_energy(g, m2, mtheta2, mu, kappa, eta);
            double dE = (E - prev_E) / (fabs(prev_E) + 1e-30);

            printf("  iter %5d: rms=%.4e rms_core=%.4e E=%.4f (dE/E=%+.2e) P_int=%.4f phi_max=%.4f\n",
                   iter, d.rms_phi, d.rms_phi_core, E, dE, d.P_int, d.phi_max);

            if (diag_fp)
                fprintf(diag_fp, "%d\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.6f\t%.6f\t%.6e\n",
                        iter, d.rms_phi, d.rms_theta, d.rms_phi_core, d.rms_theta_core,
                        E, d.P_int, d.phi_max, d.theta_rms);

            /* Adaptive step: reduce alpha if energy increased significantly */
            if (dE > 0.05 && iter > 10) {
                alpha *= 0.5;
                printf("  >> alpha reduced to %.6f (energy rose)\n", alpha);
            }

            prev_E = E;
        }
    }

    /* Final diagnostics */
    compute_residual(g, m2, mtheta2, mu, kappa, eta, R_phi, R_theta);
    DiagInfo df = compute_diag(g, R_phi, R_theta);
    double Ef = compute_energy(g, m2, mtheta2, mu, kappa, eta);
    printf("\nFinal: rms_phi=%.6e rms_theta=%.6e E=%.6f P_int=%.4f phi_max=%.4f\n",
           df.rms_phi, df.rms_theta, Ef, df.P_int, df.phi_max);
    printf("  Core residual: rms_phi=%.6e rms_theta=%.6e\n",
           df.rms_phi_core, df.rms_theta_core);
    printf("  Core residual reduction: phi %.2fx, theta %.2fx\n",
           d0.rms_phi_core / (df.rms_phi_core + 1e-30),
           d0.rms_theta_core / (df.rms_theta_core + 1e-30));
    printf("  Energy change: %.6f -> %.6f (%.4f%%)\n",
           E0, Ef, 100.0 * (Ef - E0) / fabs(E0));

    if (diag_fp) fclose(diag_fp);

    /* Write output */
    printf("\nWriting seed: %s\n", output);
    write_seed(g, output, precision, m2, A_bg);
    printf("Done.\n");

    /* Cleanup */
    for (int a = 0; a < NFIELDS; a++) {
        free(R_phi[a]); free(R_theta[a]);
        free(pin_phi[a]); free(pin_theta[a]);
    }
    rgrid_free(g);
    return 0;
}
