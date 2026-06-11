/* gen_qball_bath.c — stamp ONE Q-ball (center, delta=0) plus an isotropic
 * massless theta-sector radiation bath into an SFA seed for the v66
 * complexified kernel (complex_phi=1, 12 fields + 12 velocities = 24 columns).
 *
 * Ball (v66/THEORY.md §1,§3), for all a = 0,1,2:
 *     u_a = f(r),  v_a = 0,  udot_a = 0,  vdot_a = omega f(r)
 *
 * Bath: superposition of nmodes random traveling plane waves in the THETA
 * sector only (theta is massless: omega_k = |k|). Per mode:
 *   - random direction on the sphere, |k| uniform in [kmin,kmax]
 *   - wavevector SNAPPED to exact lattice-periodic values
 *         k_i = 2*pi*n_i / (N*dx),   N*dx = periodic box length
 *     (grid is x = -L + i*dx with dx = 2L/(N-1) and bc_type=2 wraps
 *      i=0 <-> i=N-1, so the true period is N*dx ~= 2L, NOT 2L;
 *      non-commensurate modes seam at the boundary)
 *   - random phase phi0, random component a in 0..2
 *
 * Two polarizations (bath_pol):
 *   neutral  (default / legacy): random U(1) angle alpha splits the amplitude
 *     between tu (cos alpha) and tv (sin alpha); both carry the SAME phase:
 *         tu  = A ca cos(k.x - phi0)        tv  = A sa cos(k.x - phi0)
 *         tud = A ca w_k sin(k.x - phi0)    tvd = A sa w_k sin(k.x - phi0)
 *     -> q_theta = tu*tvd - tv*tud = 0 per mode (charge-neutral bath; v67 §4
 *     showed such a bath can only ERODE a ball's charge).
 *   charged  (co-rotating): U(1)-circular polarization rotating in the SAME
 *     sense as a ball with omega>0. v66/THEORY.md §2 convention: Q > 0 for
 *     e^{+i omega t}, q = tu*tvd - tv*tud. Take Theta_a = A e^{+i(w_k t - k.x
 *     + phi0)}; with phase phi = w_k t - k.x + phi0:
 *         tu = A cos(phi), tv = A sin(phi)
 *         tud = -A w_k sin(phi), tvd = +A w_k cos(phi)
 *         q_theta = A^2 w_k (cos^2+sin^2) = +A^2 w_k > 0   (uniform in space)
 *     Evaluated at t=0 with ph = k.x - phi0 (so phi = -ph):
 *         tu  =  A cos(ph)        tv  = -A sin(ph)
 *         tud = +A w_k sin(ph)    tvd = +A w_k cos(ph)
 *     (Note: the naive tv = +A sin(k.x - w_k t + phi0) would give
 *     q_theta = -A^2 w_k < 0, i.e. COUNTER-rotating — sign matters.)
 *
 * Amplitudes are normalized NUMERICALLY (against the kernel's central
 * difference energy stencil, compute_energy_complex) so the box-averaged
 * bath energy density (kinetic + gradient over the 6 theta components)
 * equals bath_e exactly; the achieved value is printed. For charged baths
 * the achieved Noether charge density and total Q_bath are also printed.
 * bath_e = 0 reproduces a plain single-ball seed (control).
 *
 * Output: ONE SFA frame, f32, 24 columns named/typed exactly as the kernel
 * registers them (scp_sim.c sfa_add_column block / v66 SPEC §7.1).
 *
 * Build:
 *   gcc -O3 -march=native -fopenmp -o bin/gen_qball_bath \
 *       sfa/seed/gen_qball_bath.c -lzstd -lm
 *
 * Usage:
 *   gen_qball_bath N L profile omega bath_e bath_pol kmin kmax nmodes rngseed out.sfa
 *   (profile = radial_qball output: '#' headers then "r f";
 *    bath_pol = neutral | charged)
 *   Legacy 10-arg form (no bath_pol) is accepted, defaults to neutral, and is
 *   byte-identical to the pre-bath_pol generator (deprecation note printed).
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>

#define SFA_IMPLEMENTATION
#include "../format/sfa.h"

#define NCOLS 24
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

typedef struct {
    double *r, *f;
    size_t  n;
} Profile;

static void load_profile(const char *path, Profile *p) {
    FILE *fp = fopen(path, "r");
    if (!fp) { fprintf(stderr, "FATAL: cannot open profile '%s'\n", path); exit(1); }
    size_t cap = 1024;
    p->r = (double *)malloc(cap * sizeof(double));
    p->f = (double *)malloc(cap * sizeof(double));
    if (!p->r || !p->f) { fprintf(stderr, "FATAL: alloc\n"); exit(1); }
    p->n = 0;
    char line[512];
    while (fgets(line, sizeof(line), fp)) {
        const char *s = line;
        while (*s == ' ' || *s == '\t') s++;
        if (*s == '#' || *s == '\n' || *s == '\0' || *s == '\r') continue;
        double rv, fv;
        if (sscanf(s, "%lf %lf", &rv, &fv) != 2) {
            fprintf(stderr, "FATAL: parse error in '%s': %s", path, line);
            exit(1);
        }
        if (p->n > 0 && rv <= p->r[p->n - 1]) {
            fprintf(stderr, "FATAL: r not strictly increasing in '%s' at r=%g\n",
                    path, rv);
            exit(1);
        }
        if (p->n == cap) {
            cap *= 2;
            p->r = (double *)realloc(p->r, cap * sizeof(double));
            p->f = (double *)realloc(p->f, cap * sizeof(double));
            if (!p->r || !p->f) { fprintf(stderr, "FATAL: realloc\n"); exit(1); }
        }
        p->r[p->n] = rv;
        p->f[p->n] = fv;
        p->n++;
    }
    fclose(fp);
    if (p->n < 2) {
        fprintf(stderr, "FATAL: profile '%s' has <2 data points\n", path);
        exit(1);
    }
    printf("  profile '%s': %zu points, r=[%.3f,%.3f], f0=%.6f, f_last=%.3g\n",
           path, p->n, p->r[0], p->r[p->n - 1], p->f[0], p->f[p->n - 1]);
    if (fabs(p->f[p->n - 1]) > 1e-8)
        printf("  WARNING: profile '%s' truncated while still sizable "
               "(f_last=%.3g)\n", path, p->f[p->n - 1]);
}

/* piecewise-linear interpolation; flat below first r, zero above last */
static double interp(const Profile *p, double r) {
    if (r <= p->r[0]) return p->f[0];
    if (r >= p->r[p->n - 1]) return 0.0;
    size_t lo = 0, hi = p->n - 1;
    while (hi - lo > 1) {
        size_t mid = (lo + hi) / 2;
        if (p->r[mid] <= r) lo = mid; else hi = mid;
    }
    double t = (r - p->r[lo]) / (p->r[hi] - p->r[lo]);
    return p->f[lo] + t * (p->f[hi] - p->f[lo]);
}

/* ---- splitmix64 RNG (deterministic, seedable) ---- */
static uint64_t rng_state;
static uint64_t rng_u64(void) {
    uint64_t z = (rng_state += 0x9E3779B97F4A7C15ULL);
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
    return z ^ (z >> 31);
}
static double rng_uniform(void) {           /* [0,1) */
    return (double)(rng_u64() >> 11) * (1.0 / 9007199254740992.0);
}

typedef struct {
    double kx, ky, kz, omega;  /* snapped wavevector, omega = |k| */
    double phi0;               /* phase */
    double ctu, stv;           /* amplitude split: A*cos(alpha), A*sin(alpha) */
    int    a;                  /* component 0..2 */
    double amp;                /* base amplitude A */
} Mode;

/* box-averaged theta bath energy density (kinetic + gradient over the 6
 * theta components), central differences with periodic wrap — the SAME
 * stencil the kernel's compute_energy_complex uses. cols layout as main(). */
static double bath_energy_density(float **cols, int N, double dx) {
    const long NN = (long)N * N, N3 = NN * N;
    const double idx1 = 1.0 / (2.0 * dx);
    double sum = 0.0;
    #pragma omp parallel for reduction(+:sum) schedule(static)
    for (long idx = 0; idx < N3; idx++) {
        int i = (int)(idx / NN), j = (int)((idx / N) % N), k = (int)(idx % N);
        int ip = (i + 1) % N, im = (i - 1 + N) % N;
        int jp = (j + 1) % N, jm = (j - 1 + N) % N;
        int kp = (k + 1) % N, km = (k - 1 + N) % N;
        long n_ip = (long)ip * NN + (long)j * N + k;
        long n_im = (long)im * NN + (long)j * N + k;
        long n_jp = (long)i * NN + (long)jp * N + k;
        long n_jm = (long)i * NN + (long)jm * N + k;
        long n_kp = (long)i * NN + (long)j * N + kp;
        long n_km = (long)i * NN + (long)j * N + km;
        for (int a = 0; a < 3; a++) {
            /* theta (tu): col 3+a, vel col 9+a; thetaim (tv): col 15+a, vel 21+a */
            const float *tu = cols[3 + a],  *tv = cols[15 + a];
            double tud = cols[9 + a][idx], tvd = cols[21 + a][idx];
            sum += 0.5 * (tud * tud + tvd * tvd);
            double gx = (tu[n_ip] - tu[n_im]) * idx1;
            double gy = (tu[n_jp] - tu[n_jm]) * idx1;
            double gz = (tu[n_kp] - tu[n_km]) * idx1;
            double hx = (tv[n_ip] - tv[n_im]) * idx1;
            double hy = (tv[n_jp] - tv[n_jm]) * idx1;
            double hz = (tv[n_kp] - tv[n_km]) * idx1;
            sum += 0.5 * (gx * gx + gy * gy + gz * gz +
                          hx * hx + hy * hy + hz * hz);
        }
    }
    return sum / (double)N3;
}

/* box-averaged theta-sector Noether charge density
 *   q_theta = sum_a (tu_a * tvd_a - tv_a * tud_a)
 * (v66/THEORY.md §2 sign convention: Q > 0 for e^{+i omega t}). */
static double bath_charge_density(float **cols, int N) {
    const long N3 = (long)N * N * N;
    double sum = 0.0;
    #pragma omp parallel for reduction(+:sum) schedule(static)
    for (long idx = 0; idx < N3; idx++)
        for (int a = 0; a < 3; a++)
            sum += (double)cols[3 + a][idx] * cols[21 + a][idx]
                 - (double)cols[15 + a][idx] * cols[9 + a][idx];
    return sum / (double)N3;
}

int main(int argc, char **argv) {
    if (argc != 11 && argc != 12) {
        fprintf(stderr,
            "Usage: %s N L profile omega bath_e bath_pol kmin kmax nmodes rngseed out.sfa\n"
            "  Ball: u_a=f(r), v_a=0, udot=0, vdot_a=omega*f(r) at box center.\n"
            "  Bath: nmodes random theta-sector plane waves, |k| in [kmin,kmax]\n"
            "        snapped to periodic box wavevectors; box-averaged energy\n"
            "        density normalized to bath_e. bath_e=0 => plain ball seed.\n"
            "  bath_pol: neutral = random independent tu/tv phases, q_theta=0\n"
            "                      (legacy behavior)\n"
            "            charged = U(1)-circular co-rotating modes, q_theta>0\n"
            "                      (same rotation sense as a ball with omega>0)\n"
            "  Legacy 10-arg form (no bath_pol) = neutral, byte-identical output.\n",
            argv[0]);
        return 1;
    }
    int legacy = (argc == 11);
    int s_ = legacy ? 0 : 1;   /* arg index shift for args after bath_pol */
    int      N       = atoi(argv[1]);
    double   L       = atof(argv[2]);
    const char *ppath = argv[3];
    double   omega   = atof(argv[4]);
    double   bath_e  = atof(argv[5]);
    const char *polstr = legacy ? "neutral" : argv[6];
    double   kmin    = atof(argv[6 + s_]);
    double   kmax    = atof(argv[7 + s_]);
    int      nmodes  = atoi(argv[8 + s_]);
    uint64_t rngseed = (uint64_t)strtoull(argv[9 + s_], NULL, 10);
    const char *outpath = argv[10 + s_];

    int charged;
    if      (!strcmp(polstr, "neutral")) charged = 0;
    else if (!strcmp(polstr, "charged")) charged = 1;
    else {
        fprintf(stderr, "FATAL: bath_pol must be 'neutral' or 'charged' "
                "(got '%s')\n", polstr);
        return 1;
    }
    if (legacy)
        printf("NOTE: legacy 10-arg invocation (no bath_pol) is deprecated; "
               "defaulting to bath_pol=neutral.\n");

    if (N < 2 || L <= 0.0) {
        fprintf(stderr, "FATAL: bad N=%d or L=%g\n", N, L);
        return 1;
    }
    if (bath_e < 0.0) { fprintf(stderr, "FATAL: bath_e < 0\n"); return 1; }
    if (bath_e > 0.0 && (nmodes < 1 || kmin <= 0.0 || kmax < kmin)) {
        fprintf(stderr, "FATAL: bath_e>0 needs nmodes>=1 and 0<kmin<=kmax\n");
        return 1;
    }
    rng_state = rngseed;

    printf("gen_qball_bath: N=%d L=%.4f -> %s\n", N, L, outpath);
    Profile prof;
    load_profile(ppath, &prof);
    printf("  ball: omega=%.6f at center; bath_e=%.6g bath_pol=%s kmin=%.3f "
           "kmax=%.3f nmodes=%d rngseed=%llu\n", omega, bath_e, polstr, kmin,
           kmax, nmodes, (unsigned long long)rngseed);

    long   N3 = (long)N * N * N;
    double dx = 2.0 * L / (N - 1);   /* grid: x = -L + i*dx (scp_init.h) */
    double Pbox = (double)N * dx;    /* TRUE periodic length under bc_type=2 */

    /* column order MUST match kernel registration (scp_sim.c / SPEC §7.1):
     *  0-2  phi_x/y/z       (u_a)        12-14 phiim_x/y/z     (v_a)
     *  3-5  theta_x/y/z     (tu_a)       15-17 thetaim_x/y/z   (tv_a)
     *  6-8  phi_vx/vy/vz    (udot_a)     18-20 phiim_vx/vy/vz  (vdot_a)
     *  9-11 theta_vx/vy/vz  (tudot_a)    21-23 thetaim_vx/...  (tvdot_a)  */
    float *cols[NCOLS];
    for (int c = 0; c < NCOLS; c++) {
        cols[c] = (float *)calloc((size_t)N3, sizeof(float));
        if (!cols[c]) { fprintf(stderr, "FATAL: alloc column %d\n", c); return 1; }
    }

    /* ---- ball: u_a = f(r), vdot_a = omega f(r), all a ---- */
    #pragma omp parallel for schedule(static)
    for (long idx = 0; idx < N3; idx++) {
        long NN = (long)N * N;
        int i = (int)(idx / NN), j = (int)((idx / N) % N), k = (int)(idx % N);
        double x = -L + i * dx, y = -L + j * dx, z = -L + k * dx;
        double f = interp(&prof, sqrt(x * x + y * y + z * z));
        for (int a = 0; a < 3; a++) {
            cols[0 + a][idx]  = (float)f;            /* u_a    */
            cols[18 + a][idx] = (float)(omega * f);  /* vdot_a */
        }
    }

    /* ---- bath: nmodes snapped traveling plane waves in theta sector ---- */
    double achieved = 0.0;
    if (bath_e > 0.0) {
        Mode *modes = (Mode *)malloc((size_t)nmodes * sizeof(Mode));
        if (!modes) { fprintf(stderr, "FATAL: alloc modes\n"); return 1; }
        double ksnap_min = 1e300, ksnap_max = 0.0;
        for (int m = 0; m < nmodes; m++) {
            long nx = 0, ny = 0, nz = 0;
            int tries = 0;
            do {
                if (++tries > 10000) {
                    fprintf(stderr, "FATAL: cannot snap a nonzero in-band mode "
                            "(box too small for kmin=%.3f?)\n", kmin);
                    return 1;
                }
                /* random direction on the sphere */
                double cz_ = 2.0 * rng_uniform() - 1.0;
                double az  = 2.0 * M_PI * rng_uniform();
                double sz_ = sqrt(fmax(0.0, 1.0 - cz_ * cz_));
                double kmag = kmin + (kmax - kmin) * rng_uniform();
                /* snap to exact periodic wavevectors k_i = 2*pi*n_i/(N*dx) */
                nx = lround(kmag * sz_ * cos(az) * Pbox / (2.0 * M_PI));
                ny = lround(kmag * sz_ * sin(az) * Pbox / (2.0 * M_PI));
                nz = lround(kmag * cz_           * Pbox / (2.0 * M_PI));
            } while ((nx == 0 && ny == 0 && nz == 0) ||
                     labs(nx) > N / 2 || labs(ny) > N / 2 || labs(nz) > N / 2);
            Mode *M = &modes[m];
            M->kx = 2.0 * M_PI * (double)nx / Pbox;
            M->ky = 2.0 * M_PI * (double)ny / Pbox;
            M->kz = 2.0 * M_PI * (double)nz / Pbox;
            M->omega = sqrt(M->kx * M->kx + M->ky * M->ky + M->kz * M->kz);
            M->phi0 = 2.0 * M_PI * rng_uniform();
            double alpha = 2.0 * M_PI * rng_uniform();   /* U(1) orientation */
            M->ctu = cos(alpha);
            M->stv = sin(alpha);
            M->a = (int)(3.0 * rng_uniform()); if (M->a > 2) M->a = 2;
            /* equal continuum energy share per mode:
             *   neutral: <e> = A^2 omega^2 / 2 (linear mode, spatial average)
             *   charged: e   = A^2 omega^2     (circular, uniform in space) */
            M->amp = charged
                   ? sqrt(bath_e / (double)nmodes) / M->omega
                   : sqrt(2.0 * bath_e / (double)nmodes) / M->omega;
            if (M->omega < ksnap_min) ksnap_min = M->omega;
            if (M->omega > ksnap_max) ksnap_max = M->omega;
        }
        printf("  bath: %d modes, snapped |k| in [%.4f,%.4f] "
               "(period N*dx=%.4f)\n", nmodes, ksnap_min, ksnap_max, Pbox);

        #pragma omp parallel for schedule(static)
        for (long idx = 0; idx < N3; idx++) {
            long NN = (long)N * N;
            int i = (int)(idx / NN), j = (int)((idx / N) % N), k = (int)(idx % N);
            double x = -L + i * dx, y = -L + j * dx, z = -L + k * dx;
            double tu[3] = {0, 0, 0}, tv[3] = {0, 0, 0};
            double tud[3] = {0, 0, 0}, tvd[3] = {0, 0, 0};
            for (int m = 0; m < nmodes; m++) {
                const Mode *M = &modes[m];
                double ph = M->kx * x + M->ky * y + M->kz * z - M->phi0;
                double cp = cos(ph), sp = sin(ph);
                if (charged) {
                    /* Theta_a = A e^{+i(w t - k.x + phi0)} at t=0 (phase = -ph):
                     * q_theta = tu*tvd - tv*tud = +A^2 w > 0 (co-rotating) */
                    tu[M->a]  +=  M->amp * cp;
                    tv[M->a]  += -M->amp * sp;
                    tud[M->a] +=  M->amp * M->omega * sp;
                    tvd[M->a] +=  M->amp * M->omega * cp;
                } else {
                    tu[M->a]  += M->amp * M->ctu * cp;
                    tud[M->a] += M->amp * M->ctu * M->omega * sp;
                    tv[M->a]  += M->amp * M->stv * cp;
                    tvd[M->a] += M->amp * M->stv * M->omega * sp;
                }
            }
            for (int a = 0; a < 3; a++) {
                cols[3 + a][idx]  = (float)tu[a];
                cols[9 + a][idx]  = (float)tud[a];
                cols[15 + a][idx] = (float)tv[a];
                cols[21 + a][idx] = (float)tvd[a];
            }
        }
        free(modes);

        /* normalize against the kernel's discrete energy stencil */
        double e0 = bath_energy_density(cols, N, dx);
        if (e0 <= 0.0) { fprintf(stderr, "FATAL: zero bath energy\n"); return 1; }
        float s = (float)sqrt(bath_e / e0);
        #pragma omp parallel for schedule(static)
        for (long idx = 0; idx < N3; idx++)
            for (int a = 0; a < 3; a++) {
                cols[3 + a][idx]  *= s;
                cols[9 + a][idx]  *= s;
                cols[15 + a][idx] *= s;
                cols[21 + a][idx] *= s;
            }
        achieved = bath_energy_density(cols, N, dx);
        printf("  bath energy density: raw=%.6e, rescaled by %.6f -> "
               "achieved=%.6e (target %.6e)\n", e0, (double)s, achieved, bath_e);
        printf("  bath total energy (kernel metric, N^3 dx^3): %.4f\n",
               achieved * (double)N3 * dx * dx * dx);
        double qdens = bath_charge_density(cols, N);
        double Qbath = qdens * (double)N3 * dx * dx * dx;
        printf("  bath charge density (q_theta = tu*tvd - tv*tud): %.6e, "
               "total Q_bath = %.6f  [pol=%s]\n", qdens, Qbath, polstr);
        if (charged && qdens <= 0.0)
            printf("  WARNING: charged bath produced q_theta <= 0 — sign bug?\n");
    } else {
        printf("  bath_e=0: theta sector left zero (plain single-ball seed)\n");
    }

    double dt = 0.025 * dx;  /* standard dt factor (metadata only) */
    SFA *sfa = sfa_create(outpath, N, N, N, L, L, L, dt);
    if (!sfa) { fprintf(stderr, "FATAL: cannot create '%s'\n", outpath); return 1; }

    /* metadata */
    char bN[32], bL[32], bo[32], be[32], bach[32], bkmin[32], bkmax[32],
         bnm[32], bseed[32];
    snprintf(bN, sizeof(bN), "%d", N);
    snprintf(bL, sizeof(bL), "%.6f", L);
    snprintf(bo, sizeof(bo), "%.6f", omega);
    snprintf(be, sizeof(be), "%.8g", bath_e);
    snprintf(bach, sizeof(bach), "%.8g", achieved);
    snprintf(bkmin, sizeof(bkmin), "%.6f", kmin);
    snprintf(bkmax, sizeof(bkmax), "%.6f", kmax);
    snprintf(bnm, sizeof(bnm), "%d", nmodes);
    snprintf(bseed, sizeof(bseed), "%llu", (unsigned long long)rngseed);
    const char *keys[] = { "generator", "N", "L", "omega", "bath_e",
                           "bath_e_achieved", "bath_kmin", "bath_kmax",
                           "bath_nmodes", "bath_rngseed", "bath_pol" };
    const char *vals[] = { "gen_qball_bath", bN, bL, bo, be, bach,
                           bkmin, bkmax, bnm, bseed, polstr };
    /* legacy invocation: write exactly the old 10 keys (byte-identical) */
    sfa_add_kvmd(sfa, 0, 0xFFFFFFFF, 0xFFFFFFFF, keys, vals, legacy ? 10 : 11);

    /* exact names/semantics/components as the kernel registers (24 cols) */
    sfa_add_column(sfa, "phi_x",      SFA_F32, SFA_POSITION, 0);
    sfa_add_column(sfa, "phi_y",      SFA_F32, SFA_POSITION, 1);
    sfa_add_column(sfa, "phi_z",      SFA_F32, SFA_POSITION, 2);
    sfa_add_column(sfa, "theta_x",    SFA_F32, SFA_ANGLE,    0);
    sfa_add_column(sfa, "theta_y",    SFA_F32, SFA_ANGLE,    1);
    sfa_add_column(sfa, "theta_z",    SFA_F32, SFA_ANGLE,    2);
    sfa_add_column(sfa, "phi_vx",     SFA_F32, SFA_VELOCITY, 0);
    sfa_add_column(sfa, "phi_vy",     SFA_F32, SFA_VELOCITY, 1);
    sfa_add_column(sfa, "phi_vz",     SFA_F32, SFA_VELOCITY, 2);
    sfa_add_column(sfa, "theta_vx",   SFA_F32, SFA_VELOCITY, 3);
    sfa_add_column(sfa, "theta_vy",   SFA_F32, SFA_VELOCITY, 4);
    sfa_add_column(sfa, "theta_vz",   SFA_F32, SFA_VELOCITY, 5);
    sfa_add_column(sfa, "phiim_x",    SFA_F32, SFA_POSITION, 3);
    sfa_add_column(sfa, "phiim_y",    SFA_F32, SFA_POSITION, 4);
    sfa_add_column(sfa, "phiim_z",    SFA_F32, SFA_POSITION, 5);
    sfa_add_column(sfa, "thetaim_x",  SFA_F32, SFA_ANGLE,    3);
    sfa_add_column(sfa, "thetaim_y",  SFA_F32, SFA_ANGLE,    4);
    sfa_add_column(sfa, "thetaim_z",  SFA_F32, SFA_ANGLE,    5);
    sfa_add_column(sfa, "phiim_vx",   SFA_F32, SFA_VELOCITY, 6);
    sfa_add_column(sfa, "phiim_vy",   SFA_F32, SFA_VELOCITY, 7);
    sfa_add_column(sfa, "phiim_vz",   SFA_F32, SFA_VELOCITY, 8);
    sfa_add_column(sfa, "thetaim_vx", SFA_F32, SFA_VELOCITY, 9);
    sfa_add_column(sfa, "thetaim_vy", SFA_F32, SFA_VELOCITY, 10);
    sfa_add_column(sfa, "thetaim_vz", SFA_F32, SFA_VELOCITY, 11);
    sfa_finalize_header(sfa);

    void *frame_cols[NCOLS];
    for (int c = 0; c < NCOLS; c++) frame_cols[c] = cols[c];
    sfa_write_frame(sfa, 0.0, frame_cols);
    sfa_close(sfa);

    printf("gen_qball_bath: wrote %s (%d^3, f32, 24 columns, 1 frame)\n",
           outpath, N);

    for (int c = 0; c < NCOLS; c++) free(cols[c]);
    free(prof.r); free(prof.f);
    return 0;
}
