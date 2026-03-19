/*  v32_sph_two.c — SPH scalar field simulation with TWO braids
 *
 *  Phase 2: Two braids separated by D=20, measuring gravitational attraction.
 *
 *  Based on v32_sph.c (Phase 1, single braid, confirmed stable T=300).
 *  Key changes:
 *    - Domain L=30 (larger)
 *    - Two braids at x_center = -10 and +10
 *    - Extra particles near EACH braid center
 *    - Braid separation tracking via energy-density-weighted centroids
 *    - ~250K total particles
 *
 *  Build: gcc -O3 -march=native -fopenmp -o v32_two src/v32_sph_two.c -lm
 *  Run:   OMP_NUM_THREADS=16 ./v32_two
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <omp.h>
#include <sys/stat.h>
#include <time.h>

/* ================================================================
   Constants and parameters
   ================================================================ */

#define NFIELDS    3
#define PI         3.14159265358979323846
#define MAX_PARTS  600000
#define MAX_CELL   256

/* Field physics (same as Phase 1) */
static double MASS2    = 2.25;    /* m^2 = 1.5^2 */
static double MU       = -41.3;
static double KAPPA    = 50.0;

/* SPH parameters */
static double ETA      = 1.2;    /* smoothing length factor */
static double CFL      = 0.3;
static double BETA     = 0.01;   /* transport strength */
static double EPS_DIV  = 1e-10;  /* prevent /0 in transport */

/* Braid initialization */
static double A_AMP[3]    = {0.8, 0.8, 0.8};
static double DELTA[3]    = {0.0, 3.0005, 4.4325};
static double R_TUBE      = 3.0;
static double ELLIP       = 0.3325;
static double K_FAC       = 1.0;
static double A_BG        = 0.1;

/* Two-braid geometry */
static double BRAID_X1 = -10.0;   /* left braid center */
static double BRAID_X2 = +10.0;   /* right braid center */
static double R_BRAID  = 5.0;     /* extra-particle radius around each braid */
static double E_THRESH_MULT = 5.0; /* threshold for centroid: E > 5x average */

/* ================================================================
   Particle structure
   ================================================================ */

typedef struct {
    double x, y, z;          /* position */
    double phi[NFIELDS];     /* field values */
    double vel[NFIELDS];     /* field time derivatives */
    double vt[3];            /* transport velocity */
    double h;                /* smoothing length */
    double rho_sph;          /* SPH density */
    double m_sph;            /* SPH mass (fixed) */
} Particle;

/* ================================================================
   Spatial hash grid
   ================================================================ */

typedef struct {
    int *cell_start;
    int *cell_count;
    int *sorted_ids;
    int nx, ny, nz;
    double cell_size;
    double xmin, ymin, zmin;
    double xmax, ymax, zmax;
    int n_cells;
    int use_pbc;
    double Lx, Ly, Lz;
} HashGrid;

static inline int hash_cell(HashGrid *hg, int cx, int cy, int cz) {
    if (hg->use_pbc) {
        cx = ((cx % hg->nx) + hg->nx) % hg->nx;
        cy = ((cy % hg->ny) + hg->ny) % hg->ny;
        cz = ((cz % hg->nz) + hg->nz) % hg->nz;
    } else {
        if (cx < 0 || cx >= hg->nx) return -1;
        if (cy < 0 || cy >= hg->ny) return -1;
        if (cz < 0 || cz >= hg->nz) return -1;
    }
    return cx * hg->ny * hg->nz + cy * hg->nz + cz;
}

static void hash_build(HashGrid *hg, Particle *parts, int np,
                        double cell_size, double L, int pbc) {
    hg->cell_size = cell_size;
    hg->use_pbc = pbc;
    hg->xmin = -L; hg->ymin = -L; hg->zmin = -L;
    hg->xmax =  L; hg->ymax =  L; hg->zmax =  L;
    hg->Lx = 2.0 * L; hg->Ly = 2.0 * L; hg->Lz = 2.0 * L;
    hg->nx = (int)ceil(hg->Lx / cell_size);
    hg->ny = (int)ceil(hg->Ly / cell_size);
    hg->nz = (int)ceil(hg->Lz / cell_size);
    hg->n_cells = hg->nx * hg->ny * hg->nz;

    hg->cell_count = realloc(hg->cell_count, hg->n_cells * sizeof(int));
    hg->cell_start = realloc(hg->cell_start, hg->n_cells * sizeof(int));
    hg->sorted_ids = realloc(hg->sorted_ids, np * sizeof(int));

    memset(hg->cell_count, 0, hg->n_cells * sizeof(int));

    for (int i = 0; i < np; i++) {
        int cx = (int)floor((parts[i].x - hg->xmin) / cell_size);
        int cy = (int)floor((parts[i].y - hg->ymin) / cell_size);
        int cz = (int)floor((parts[i].z - hg->zmin) / cell_size);
        if (cx < 0) cx = 0; if (cx >= hg->nx) cx = hg->nx - 1;
        if (cy < 0) cy = 0; if (cy >= hg->ny) cy = hg->ny - 1;
        if (cz < 0) cz = 0; if (cz >= hg->nz) cz = hg->nz - 1;
        int c = cx * hg->ny * hg->nz + cy * hg->nz + cz;
        hg->cell_count[c]++;
    }

    hg->cell_start[0] = 0;
    for (int c = 1; c < hg->n_cells; c++)
        hg->cell_start[c] = hg->cell_start[c-1] + hg->cell_count[c-1];

    int *fill = calloc(hg->n_cells, sizeof(int));
    for (int i = 0; i < np; i++) {
        int cx = (int)floor((parts[i].x - hg->xmin) / cell_size);
        int cy = (int)floor((parts[i].y - hg->ymin) / cell_size);
        int cz = (int)floor((parts[i].z - hg->zmin) / cell_size);
        if (cx < 0) cx = 0; if (cx >= hg->nx) cx = hg->nx - 1;
        if (cy < 0) cy = 0; if (cy >= hg->ny) cy = hg->ny - 1;
        if (cz < 0) cz = 0; if (cz >= hg->nz) cz = hg->nz - 1;
        int c = cx * hg->ny * hg->nz + cy * hg->nz + cz;
        hg->sorted_ids[hg->cell_start[c] + fill[c]] = i;
        fill[c]++;
    }
    free(fill);
}

static void hash_free(HashGrid *hg) {
    free(hg->cell_count);
    free(hg->cell_start);
    free(hg->sorted_ids);
    hg->cell_count = NULL;
    hg->cell_start = NULL;
    hg->sorted_ids = NULL;
}

/* ================================================================
   SPH Cubic Spline Kernel (3D)
   ================================================================ */

static inline double kernel_W(double r, double h) {
    double q = r / h;
    double norm = 1.0 / (PI * h * h * h);
    if (q <= 1.0) {
        return norm * (1.0 - 1.5*q*q + 0.75*q*q*q);
    } else if (q <= 2.0) {
        double t = 2.0 - q;
        return norm * 0.25 * t * t * t;
    }
    return 0.0;
}

static inline double kernel_dWdr(double r, double h) {
    double q = r / h;
    double norm = 1.0 / (PI * h * h * h * h);
    if (q <= 1.0) {
        return norm * (-3.0*q + 2.25*q*q);
    } else if (q <= 2.0) {
        double t = 2.0 - q;
        return norm * (-0.75 * t * t);
    }
    return 0.0;
}

/* ================================================================
   Periodic distance helper
   ================================================================ */

static inline double pbc_wrap(double dx, double Lbox) {
    if (dx >  0.5 * Lbox) dx -= Lbox;
    if (dx < -0.5 * Lbox) dx += Lbox;
    return dx;
}

static inline void pbc_dr(double *dx, double *dy, double *dz,
                           double Lx, double Ly, double Lz, int pbc) {
    if (pbc) {
        *dx = pbc_wrap(*dx, Lx);
        *dy = pbc_wrap(*dy, Ly);
        *dz = pbc_wrap(*dz, Lz);
    }
}

/* ================================================================
   Wrap particle positions into box
   ================================================================ */

static void wrap_positions(Particle *parts, int np, double L) {
    double Lbox = 2.0 * L;
    for (int i = 0; i < np; i++) {
        while (parts[i].x < -L) parts[i].x += Lbox;
        while (parts[i].x >= L) parts[i].x -= Lbox;
        while (parts[i].y < -L) parts[i].y += Lbox;
        while (parts[i].y >= L) parts[i].y -= Lbox;
        while (parts[i].z < -L) parts[i].z += Lbox;
        while (parts[i].z >= L) parts[i].z -= Lbox;
    }
}

/* ================================================================
   Compute SPH density and adaptive smoothing length
   ================================================================ */

static void compute_density(Particle *parts, int np, HashGrid *hg) {
    double Lx = hg->Lx, Ly = hg->Ly, Lz = hg->Lz;
    int pbc = hg->use_pbc;

    #pragma omp parallel for schedule(dynamic, 256)
    for (int i = 0; i < np; i++) {
        double xi = parts[i].x, yi = parts[i].y, zi = parts[i].z;
        double hi = parts[i].h;
        double rho = 0.0;

        int cxi = (int)floor((xi - hg->xmin) / hg->cell_size);
        int cyi = (int)floor((yi - hg->ymin) / hg->cell_size);
        int czi = (int)floor((zi - hg->zmin) / hg->cell_size);

        for (int dcx = -1; dcx <= 1; dcx++)
        for (int dcy = -1; dcy <= 1; dcy++)
        for (int dcz = -1; dcz <= 1; dcz++) {
            int c = hash_cell(hg, cxi+dcx, cyi+dcy, czi+dcz);
            if (c < 0) continue;
            int start = hg->cell_start[c];
            int count = hg->cell_count[c];
            for (int jj = 0; jj < count; jj++) {
                int j = hg->sorted_ids[start + jj];
                double dx = xi - parts[j].x;
                double dy = yi - parts[j].y;
                double dz = zi - parts[j].z;
                pbc_dr(&dx, &dy, &dz, Lx, Ly, Lz, pbc);
                double r = sqrt(dx*dx + dy*dy + dz*dz);
                if (r < 2.0 * hi)
                    rho += parts[j].m_sph * kernel_W(r, hi);
            }
        }
        parts[i].rho_sph = rho;
    }
}

static void update_smoothing_length(Particle *parts, int np) {
    for (int i = 0; i < np; i++) {
        if (parts[i].rho_sph > 1e-30) {
            parts[i].h = ETA * cbrt(parts[i].m_sph / parts[i].rho_sph);
        }
    }
}

static void converge_density_h(Particle *parts, int np, HashGrid *hg,
                                double L) {
    for (int iter = 0; iter < 3; iter++) {
        double h_max = 0;
        for (int i = 0; i < np; i++)
            if (parts[i].h > h_max) h_max = parts[i].h;
        hash_build(hg, parts, np, 2.0 * h_max, L, 1);
        compute_density(parts, np, hg);
        update_smoothing_length(parts, np);
    }
}

/* ================================================================
   SPH Laplacian (Brookshaw) and field forces
   ================================================================ */

static void compute_forces(Particle *parts, int np, HashGrid *hg,
                            double *acc,
                            double *energy_density) {
    double Lx = hg->Lx, Ly = hg->Ly, Lz = hg->Lz;
    int pbc = hg->use_pbc;

    #pragma omp parallel for schedule(dynamic, 128)
    for (int i = 0; i < np; i++) {
        double xi = parts[i].x, yi = parts[i].y, zi = parts[i].z;
        double hi = parts[i].h;
        double h2_reg = 0.01 * hi * hi;

        double lap[NFIELDS] = {0, 0, 0};

        int cxi = (int)floor((xi - hg->xmin) / hg->cell_size);
        int cyi = (int)floor((yi - hg->ymin) / hg->cell_size);
        int czi = (int)floor((zi - hg->zmin) / hg->cell_size);

        for (int dcx = -1; dcx <= 1; dcx++)
        for (int dcy = -1; dcy <= 1; dcy++)
        for (int dcz = -1; dcz <= 1; dcz++) {
            int c = hash_cell(hg, cxi+dcx, cyi+dcy, czi+dcz);
            if (c < 0) continue;
            int start = hg->cell_start[c];
            int count = hg->cell_count[c];
            for (int jj = 0; jj < count; jj++) {
                int j = hg->sorted_ids[start + jj];
                if (j == i) continue;

                double dx = xi - parts[j].x;
                double dy = yi - parts[j].y;
                double dz = zi - parts[j].z;
                pbc_dr(&dx, &dy, &dz, Lx, Ly, Lz, pbc);

                double r2 = dx*dx + dy*dy + dz*dz;
                double r = sqrt(r2);
                if (r >= 2.0 * hi) continue;

                double dWdr = kernel_dWdr(r, hi);
                double r_dot_gradW = r * dWdr;
                double fac = 2.0 * (parts[j].m_sph / (parts[j].rho_sph + 1e-30))
                           * r_dot_gradW / (r2 + h2_reg);

                for (int a = 0; a < NFIELDS; a++) {
                    lap[a] += (parts[i].phi[a] - parts[j].phi[a]) * fac;
                }
            }
        }

        double p0 = parts[i].phi[0], p1 = parts[i].phi[1], p2 = parts[i].phi[2];
        double P = p0 * p1 * p2;
        double denom = 1.0 + KAPPA * P * P;
        double mu_P_d2 = MU * P / (denom * denom);

        for (int a = 0; a < NFIELDS; a++) {
            double dPda;
            if      (a == 0) dPda = p1 * p2;
            else if (a == 1) dPda = p0 * p2;
            else             dPda = p0 * p1;

            acc[i * NFIELDS + a] = lap[a] - MASS2 * parts[i].phi[a] - mu_P_d2 * dPda;
        }

        if (energy_density) {
            double ek = 0, ep_mass = 0;
            for (int a = 0; a < NFIELDS; a++) {
                ek += 0.5 * parts[i].vel[a] * parts[i].vel[a];
                ep_mass += 0.5 * MASS2 * parts[i].phi[a] * parts[i].phi[a];
            }
            double ep_triple = (MU / 2.0) * P * P / (1.0 + KAPPA * P * P);
            energy_density[i] = ek + ep_mass + ep_triple;
            if (energy_density[i] < 0) energy_density[i] = 0;
        }
    }
}

/* ================================================================
   Transport velocity: v_transport = -beta * grad(E) / (E + eps)
   ================================================================ */

static void compute_transport(Particle *parts, int np, HashGrid *hg,
                               double *energy_density) {
    double Lx = hg->Lx, Ly = hg->Ly, Lz = hg->Lz;
    int pbc = hg->use_pbc;

    #pragma omp parallel for schedule(dynamic, 256)
    for (int i = 0; i < np; i++) {
        double xi = parts[i].x, yi = parts[i].y, zi = parts[i].z;
        double hi = parts[i].h;
        double Ei = energy_density[i];

        double grad_E[3] = {0, 0, 0};

        int cxi = (int)floor((xi - hg->xmin) / hg->cell_size);
        int cyi = (int)floor((yi - hg->ymin) / hg->cell_size);
        int czi = (int)floor((zi - hg->zmin) / hg->cell_size);

        for (int dcx = -1; dcx <= 1; dcx++)
        for (int dcy = -1; dcy <= 1; dcy++)
        for (int dcz = -1; dcz <= 1; dcz++) {
            int c = hash_cell(hg, cxi+dcx, cyi+dcy, czi+dcz);
            if (c < 0) continue;
            int start = hg->cell_start[c];
            int count = hg->cell_count[c];
            for (int jj = 0; jj < count; jj++) {
                int j = hg->sorted_ids[start + jj];
                if (j == i) continue;

                double ddx = xi - parts[j].x;
                double ddy = yi - parts[j].y;
                double ddz = zi - parts[j].z;
                pbc_dr(&ddx, &ddy, &ddz, Lx, Ly, Lz, pbc);

                double r2 = ddx*ddx + ddy*ddy + ddz*ddz;
                double r = sqrt(r2);
                if (r >= 2.0 * hi || r < 1e-15) continue;

                double dWdr = kernel_dWdr(r, hi);
                double fac = (parts[j].m_sph / (parts[j].rho_sph + 1e-30))
                           * (energy_density[j] - Ei)
                           * dWdr / r;

                grad_E[0] += fac * ddx;
                grad_E[1] += fac * ddy;
                grad_E[2] += fac * ddz;
            }
        }

        double inv_E = 1.0 / (Ei + EPS_DIV);
        parts[i].vt[0] = -BETA * grad_E[0] * inv_E;
        parts[i].vt[1] = -BETA * grad_E[1] * inv_E;
        parts[i].vt[2] = -BETA * grad_E[2] * inv_E;
    }
}

/* ================================================================
   Initialization: TWO braids + background
   ================================================================ */

static int init_particles(Particle *parts, int np_target, double L) {
    double Lbox = 2.0 * L;
    double vol = Lbox * Lbox * Lbox;

    /* Extra particles in TWO braid spheres */
    double vol_braid = (4.0/3.0) * PI * R_BRAID * R_BRAID * R_BRAID;
    double frac_braid = vol_braid / vol;
    int np_uniform = np_target;
    int np_extra_per = (int)(2.0 * frac_braid * np_target);  /* 2x extra per braid */
    int np_extra = 2 * np_extra_per;  /* for BOTH braids */
    int np_total = np_uniform + np_extra;
    if (np_total > MAX_PARTS) {
        np_extra_per = (MAX_PARTS - np_uniform) / 2;
        np_extra = 2 * np_extra_per;
        np_total = np_uniform + np_extra;
    }

    printf("  Initializing: %d uniform + %d extra (%d per braid) = %d total\n",
           np_uniform, np_extra, np_extra_per, np_total);
    printf("  Braid 1 at x=%.1f, Braid 2 at x=%.1f, separation D=%.1f\n",
           BRAID_X1, BRAID_X2, BRAID_X2 - BRAID_X1);

    double m_sph = vol / np_uniform;
    double h_init = ETA * cbrt(m_sph);

    srand48(42);

    /* Uniform particles */
    for (int i = 0; i < np_uniform; i++) {
        parts[i].x = -L + drand48() * Lbox;
        parts[i].y = -L + drand48() * Lbox;
        parts[i].z = -L + drand48() * Lbox;
        parts[i].m_sph = m_sph;
        parts[i].h = h_init;
        parts[i].rho_sph = 1.0;
        parts[i].vt[0] = parts[i].vt[1] = parts[i].vt[2] = 0;
    }

    /* Extra particles in LEFT braid region (rejection sampling in sphere r<R_BRAID around BRAID_X1) */
    int n_added = 0;
    while (n_added < np_extra_per) {
        double x = -R_BRAID + drand48() * 2.0 * R_BRAID;
        double y = -R_BRAID + drand48() * 2.0 * R_BRAID;
        double z = -R_BRAID + drand48() * 2.0 * R_BRAID;
        if (x*x + y*y + z*z >= R_BRAID * R_BRAID) continue;
        int idx = np_uniform + n_added;
        parts[idx].x = x + BRAID_X1;  /* shift to braid 1 center */
        parts[idx].y = y;
        parts[idx].z = z;
        parts[idx].m_sph = m_sph;
        parts[idx].h = h_init * 0.7;
        parts[idx].rho_sph = 1.0;
        parts[idx].vt[0] = parts[idx].vt[1] = parts[idx].vt[2] = 0;
        n_added++;
    }

    /* Extra particles in RIGHT braid region */
    int n_added2 = 0;
    while (n_added2 < np_extra_per) {
        double x = -R_BRAID + drand48() * 2.0 * R_BRAID;
        double y = -R_BRAID + drand48() * 2.0 * R_BRAID;
        double z = -R_BRAID + drand48() * 2.0 * R_BRAID;
        if (x*x + y*y + z*z >= R_BRAID * R_BRAID) continue;
        int idx = np_uniform + np_extra_per + n_added2;
        parts[idx].x = x + BRAID_X2;  /* shift to braid 2 center */
        parts[idx].y = y;
        parts[idx].z = z;
        parts[idx].m_sph = m_sph;
        parts[idx].h = h_init * 0.7;
        parts[idx].rho_sph = 1.0;
        parts[idx].vt[0] = parts[idx].vt[1] = parts[idx].vt[2] = 0;
        n_added2++;
    }
    np_total = np_uniform + n_added + n_added2;

    /* Initialize field values: helical braid at EACH center + background */
    double kw = K_FAC * PI / L;
    double omega = sqrt(kw * kw + MASS2);
    double sx = 1.0 + ELLIP, sy = 1.0 - ELLIP;
    double inv_2R2 = 1.0 / (2.0 * R_TUBE * R_TUBE);
    double k_bg = PI / L;
    double om_bg = sqrt(k_bg * k_bg + MASS2);

    for (int i = 0; i < np_total; i++) {
        double x = parts[i].x, y = parts[i].y, z = parts[i].z;
        for (int a = 0; a < NFIELDS; a++) {
            /* Background */
            double ph_bg = k_bg * z + 2.0 * PI * a / 3.0;
            double val = A_BG * cos(ph_bg);
            double dval = om_bg * A_BG * sin(ph_bg);

            /* Braid 1: centered at (BRAID_X1, 0, 0) */
            {
                double dx1 = x - BRAID_X1;
                double xr = dx1, yr = y;
                double r2e = xr*xr / (sx*sx) + yr*yr / (sy*sy);
                double env = exp(-r2e * inv_2R2);
                double ph = kw * z + DELTA[a];
                val += A_AMP[a] * env * cos(ph);
                dval += omega * A_AMP[a] * env * sin(ph);
            }

            /* Braid 2: centered at (BRAID_X2, 0, 0) */
            {
                double dx2 = x - BRAID_X2;
                double xr = dx2, yr = y;
                double r2e = xr*xr / (sx*sx) + yr*yr / (sy*sy);
                double env = exp(-r2e * inv_2R2);
                double ph = kw * z + DELTA[a] + PI;  /* phase offset for variety */
                val += A_AMP[a] * env * cos(ph);
                dval += omega * A_AMP[a] * env * sin(ph);
            }

            parts[i].phi[a] = val;
            parts[i].vel[a] = dval;
        }
    }

    printf("  m_sph=%.6f, h_init=%.4f, kw=%.4f, omega=%.4f\n",
           m_sph, h_init, kw, omega);

    return np_total;
}

/* ================================================================
   Braid separation measurement
   ================================================================ */

typedef struct {
    double x_left, y_left, z_left;
    double x_right, y_right, z_right;
    double D;           /* separation: x_right - x_left */
    int n_left, n_right;
    double E_left, E_right;
    double E_threshold;
} SeparationResult;

static SeparationResult measure_separation(Particle *parts, int np,
                                             double *energy_density) {
    SeparationResult res;
    memset(&res, 0, sizeof(res));

    /* Compute average energy density for threshold */
    double E_sum = 0;
    for (int i = 0; i < np; i++)
        E_sum += energy_density[i];
    double E_avg = E_sum / np;
    res.E_threshold = E_THRESH_MULT * E_avg;

    /* Energy-density-weighted centroids, split by x<0 and x>0 */
    /* Only include particles above threshold */
    double wx_left = 0, wy_left = 0, wz_left = 0, wt_left = 0;
    double wx_right = 0, wy_right = 0, wz_right = 0, wt_right = 0;
    int n_left = 0, n_right = 0;

    for (int i = 0; i < np; i++) {
        double E_i = energy_density[i];
        if (E_i < res.E_threshold) continue;

        if (parts[i].x < 0) {
            wx_left += E_i * parts[i].x;
            wy_left += E_i * parts[i].y;
            wz_left += E_i * parts[i].z;
            wt_left += E_i;
            n_left++;
        } else {
            wx_right += E_i * parts[i].x;
            wy_right += E_i * parts[i].y;
            wz_right += E_i * parts[i].z;
            wt_right += E_i;
            n_right++;
        }
    }

    if (wt_left > 0) {
        res.x_left = wx_left / wt_left;
        res.y_left = wy_left / wt_left;
        res.z_left = wz_left / wt_left;
    }
    if (wt_right > 0) {
        res.x_right = wx_right / wt_right;
        res.y_right = wy_right / wt_right;
        res.z_right = wz_right / wt_right;
    }

    res.D = res.x_right - res.x_left;
    res.n_left = n_left;
    res.n_right = n_right;
    res.E_left = wt_left;
    res.E_right = wt_right;

    return res;
}

/* ================================================================
   Diagnostics — two-braid version
   ================================================================ */

static void diagnostics(Particle *parts, int np, double *energy_density,
                         double t, int step, FILE *log_fp) {
    double E_total = 0;
    double max_phi2 = 0;
    double h_sum = 0;
    double vt_max = 0;

    for (int i = 0; i < np; i++) {
        double vol_i = parts[i].m_sph / (parts[i].rho_sph + 1e-30);
        E_total += energy_density[i] * vol_i;

        double phi2 = 0;
        for (int a = 0; a < NFIELDS; a++)
            phi2 += parts[i].phi[a] * parts[i].phi[a];
        if (phi2 > max_phi2) max_phi2 = phi2;

        h_sum += parts[i].h;

        double vt2 = parts[i].vt[0]*parts[i].vt[0]
                    + parts[i].vt[1]*parts[i].vt[1]
                    + parts[i].vt[2]*parts[i].vt[2];
        if (vt2 > vt_max) vt_max = vt2;
    }

    /* Measure braid separation */
    SeparationResult sep = measure_separation(parts, np, energy_density);

    printf("  t=%7.2f  step=%6d  D=%.4f  E=%.4e  max|phi|2=%.4f"
           "  n_L=%d  n_R=%d  <h>=%.4f  |vt|=%.4f\n",
           t, step, sep.D, E_total, max_phi2,
           sep.n_left, sep.n_right, h_sum/np, sqrt(vt_max));

    /* Also log separation data to file */
    if (log_fp) {
        fprintf(log_fp, "%8.2f  %12.6f  %12.6f  %12.6f  %12.6f  %12.6f  %12.6f  "
                "%12.4e  %8d  %8d  %12.4e  %12.4e  %8.4f\n",
                t, sep.D, sep.x_left, sep.y_left, sep.x_right, sep.y_right,
                sep.E_threshold, E_total, sep.n_left, sep.n_right,
                sep.E_left, sep.E_right, max_phi2);
        fflush(log_fp);
    }
}

/* ================================================================
   Save snapshot (binary)
   ================================================================ */

static void save_snapshot(Particle *parts, int np, double t, const char *dir) {
    char fn[512];
    snprintf(fn, sizeof(fn), "%s/sph2_t%04d.bin", dir, (int)(t + 0.5));
    FILE *fp = fopen(fn, "wb");
    if (!fp) { printf("  ERROR: cannot open %s\n", fn); return; }

    fwrite(&np, sizeof(int), 1, fp);
    fwrite(&t, sizeof(double), 1, fp);

    for (int i = 0; i < np; i++) {
        fwrite(&parts[i].x, sizeof(double), 1, fp);
        fwrite(&parts[i].y, sizeof(double), 1, fp);
        fwrite(&parts[i].z, sizeof(double), 1, fp);
        fwrite(parts[i].phi, sizeof(double), NFIELDS, fp);
        fwrite(parts[i].vel, sizeof(double), NFIELDS, fp);
        fwrite(&parts[i].h, sizeof(double), 1, fp);
        fwrite(&parts[i].rho_sph, sizeof(double), 1, fp);
    }
    fclose(fp);

    double mb = (double)np * (3 + NFIELDS + NFIELDS + 1 + 1) * 8.0 / 1e6;
    printf("  Saved: %s (%.1f MB, %d particles)\n", fn, mb, np);
}

/* ================================================================
   Main simulation loop
   ================================================================ */

int main(int argc, char **argv) {
    /* Default parameters */
    int np_target = 240000;
    double T_end = 300.0;
    double L = 30.0;
    int nthreads = 16;
    char outdir[256] = "data/sph2";
    double diag_interval = 5.0;
    double snap_interval = 50.0;

    /* Parse arguments */
    for (int i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "-np") && i+1 < argc) np_target = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-T") && i+1 < argc) T_end = atof(argv[++i]);
        else if (!strcmp(argv[i], "-L") && i+1 < argc) L = atof(argv[++i]);
        else if (!strcmp(argv[i], "-beta") && i+1 < argc) BETA = atof(argv[++i]);
        else if (!strcmp(argv[i], "-threads") && i+1 < argc) nthreads = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-o") && i+1 < argc) strncpy(outdir, argv[++i], 255);
        else if (!strcmp(argv[i], "-diag") && i+1 < argc) diag_interval = atof(argv[++i]);
        else if (!strcmp(argv[i], "-snap") && i+1 < argc) snap_interval = atof(argv[++i]);
        else if (!strcmp(argv[i], "-D") && i+1 < argc) {
            double D = atof(argv[++i]);
            BRAID_X1 = -D/2.0;
            BRAID_X2 = +D/2.0;
        }
    }

    omp_set_num_threads(nthreads);
    mkdir(outdir, 0755);

    printf("V32 SPH Two-Braid Gravity Test\n");
    printf("================================\n");
    printf("  N_particles (target): %d\n", np_target);
    printf("  T_end: %.1f\n", T_end);
    printf("  L: %.1f (domain [-L,L]^3)\n", L);
    printf("  beta (transport): %.4f\n", BETA);
    printf("  m^2=%.2f, mu=%.1f, kappa=%.0f\n", MASS2, MU, KAPPA);
    printf("  Braid 1: x=%.1f, Braid 2: x=%.1f, D=%.1f\n",
           BRAID_X1, BRAID_X2, BRAID_X2-BRAID_X1);
    printf("  E threshold: %.1fx average\n", E_THRESH_MULT);
    printf("  OMP threads: %d\n", nthreads);
    printf("  Output: %s/\n\n", outdir);

    /* Allocate particles */
    Particle *parts = calloc(MAX_PARTS, sizeof(Particle));
    if (!parts) { printf("ERROR: cannot allocate %d particles\n", MAX_PARTS); return 1; }

    /* Initialize */
    printf("Initialization:\n");
    int np = init_particles(parts, np_target, L);

    /* Work arrays */
    double *acc = calloc(MAX_PARTS * (size_t)NFIELDS, sizeof(double));
    double *energy_density = calloc(MAX_PARTS, sizeof(double));

    /* Hash grid */
    HashGrid hg;
    memset(&hg, 0, sizeof(hg));

    /* Initial density convergence */
    printf("\nConverging density and smoothing lengths...\n");
    converge_density_h(parts, np, &hg, L);

    /* Find min h for CFL */
    double h_min = 1e30;
    for (int i = 0; i < np; i++)
        if (parts[i].h < h_min) h_min = parts[i].h;
    double dt = CFL * h_min;
    printf("  h_min=%.4f, dt=%.6f\n", h_min, dt);

    /* Initial forces */
    printf("\nComputing initial forces...\n");
    compute_forces(parts, np, &hg, acc, energy_density);

    /* Open log file for separation data */
    char log_fn[512];
    snprintf(log_fn, sizeof(log_fn), "%s/separation.dat", outdir);
    FILE *log_fp = fopen(log_fn, "w");
    if (log_fp) {
        fprintf(log_fp, "# V32 Two-Braid Separation Log\n");
        fprintf(log_fp, "# t  D  x_left  y_left  x_right  y_right  "
                "E_thresh  E_total  n_left  n_right  E_left  E_right  max_phi2\n");
    }

    /* Save initial snapshot */
    save_snapshot(parts, np, 0.0, outdir);
    diagnostics(parts, np, energy_density, 0.0, 0, log_fp);

    /* ==================== Main time loop ==================== */
    printf("\nStarting time integration...\n");
    double t = 0.0;
    int step = 0;
    double next_diag = diag_interval;
    double next_snap = snap_interval;
    double wall_start = omp_get_wtime();
    double last_report = wall_start;

    int nsteps = (int)(T_end / dt) + 1;
    printf("  Total steps: ~%d (dt=%.6f)\n\n", nsteps, dt);

    while (t < T_end) {
        /* === Velocity Verlet: kick-drift-kick === */

        /* Half kick */
        #pragma omp parallel for schedule(static)
        for (int i = 0; i < np; i++) {
            for (int a = 0; a < NFIELDS; a++) {
                parts[i].vel[a] += 0.5 * dt * acc[i * NFIELDS + a];
            }
        }

        /* Drift: phi += dt*vel */
        #pragma omp parallel for schedule(static)
        for (int i = 0; i < np; i++) {
            for (int a = 0; a < NFIELDS; a++) {
                parts[i].phi[a] += dt * parts[i].vel[a];
            }
        }

        /* Particle transport: move positions */
        #pragma omp parallel for schedule(static)
        for (int i = 0; i < np; i++) {
            parts[i].x += dt * parts[i].vt[0];
            parts[i].y += dt * parts[i].vt[1];
            parts[i].z += dt * parts[i].vt[2];
        }

        /* Wrap positions (periodic BC) */
        wrap_positions(parts, np, L);

        /* Rebuild hash, recompute density + h */
        if (step % 10 == 0) {
            converge_density_h(parts, np, &hg, L);
        } else {
            double h_max = 0;
            for (int i = 0; i < np; i++)
                if (parts[i].h > h_max) h_max = parts[i].h;
            hash_build(&hg, parts, np, 2.0 * h_max, L, 1);
            compute_density(parts, np, &hg);
            update_smoothing_length(parts, np);
        }

        /* Recompute forces */
        compute_forces(parts, np, &hg, acc, energy_density);

        /* Second half kick */
        #pragma omp parallel for schedule(static)
        for (int i = 0; i < np; i++) {
            for (int a = 0; a < NFIELDS; a++) {
                parts[i].vel[a] += 0.5 * dt * acc[i * NFIELDS + a];
            }
        }

        /* Compute transport velocity for next step */
        compute_transport(parts, np, &hg, energy_density);

        t += dt;
        step++;

        /* Update dt based on current h_min */
        if (step % 50 == 0) {
            h_min = 1e30;
            for (int i = 0; i < np; i++)
                if (parts[i].h < h_min) h_min = parts[i].h;
            dt = CFL * h_min;
        }

        /* Diagnostics */
        if (t >= next_diag) {
            diagnostics(parts, np, energy_density, t, step, log_fp);
            next_diag = t + diag_interval;
        }

        /* Snapshot */
        if (t >= next_snap) {
            save_snapshot(parts, np, t, outdir);
            next_snap = t + snap_interval;
        }

        /* Wall-time progress report every 30s */
        double wall_now = omp_get_wtime();
        if (wall_now - last_report > 30.0) {
            double frac = t / T_end;
            double elapsed = wall_now - wall_start;
            double eta = (frac > 0.01) ? elapsed / frac * (1.0 - frac) : 0;
            printf("  [PROGRESS] t=%.1f/%.1f (%.1f%%)  step=%d  elapsed=%.0fs  ETA=%.0fs  dt=%.6f\n",
                   t, T_end, 100*frac, step, elapsed, eta, dt);
            fflush(stdout);
            last_report = wall_now;
        }

        /* Blowup check */
        if (step % 100 == 0) {
            double max_v = 0;
            for (int i = 0; i < np; i += 37) {
                for (int a = 0; a < NFIELDS; a++) {
                    double v = fabs(parts[i].phi[a]);
                    if (v > max_v) max_v = v;
                }
            }
            if (max_v > 100.0) {
                printf("\n  *** BLOWUP at t=%.2f step=%d (max|phi|=%.1f) ***\n", t, step, max_v);
                save_snapshot(parts, np, t, outdir);
                diagnostics(parts, np, energy_density, t, step, log_fp);
                break;
            }
        }
    }

    double wall_end = omp_get_wtime();
    printf("\nSimulation complete: t=%.2f, steps=%d, wall=%.1fs\n",
           t, step, wall_end - wall_start);

    /* Final snapshot + diagnostics */
    save_snapshot(parts, np, t, outdir);
    diagnostics(parts, np, energy_density, t, step, log_fp);

    /* ==================== Post-run analysis ==================== */
    printf("\n=== Post-run Analysis ===\n");

    /* Radial density profile around EACH braid */
    int n_bins = 30;
    double r_max = 15.0;
    double dr = r_max / n_bins;

    for (int braid = 0; braid < 2; braid++) {
        double cx = (braid == 0) ? BRAID_X1 : BRAID_X2;
        printf("\n  Radial profile around Braid %d (center x=%.1f):\n", braid+1, cx);
        printf("  %8s  %8s  %8s  %10s\n", "r", "n_part", "rho_sph", "E_density");

        double *rho_prof = calloc(n_bins, sizeof(double));
        double *E_prof = calloc(n_bins, sizeof(double));
        int *count_prof = calloc(n_bins, sizeof(int));

        for (int i = 0; i < np; i++) {
            double dx = parts[i].x - cx;
            double dy = parts[i].y;
            double dz = parts[i].z;
            /* Use PBC for distance to center */
            double Lbox = 2.0 * L;
            dx = pbc_wrap(dx, Lbox);
            dy = pbc_wrap(dy, Lbox);
            dz = pbc_wrap(dz, Lbox);
            double r = sqrt(dx*dx + dy*dy + dz*dz);
            int bin = (int)(r / dr);
            if (bin >= 0 && bin < n_bins) {
                rho_prof[bin] += parts[i].rho_sph;
                E_prof[bin] += energy_density[i];
                count_prof[bin]++;
            }
        }

        for (int b = 0; b < n_bins; b++) {
            double r = (b + 0.5) * dr;
            double rho_avg = (count_prof[b] > 0) ? rho_prof[b] / count_prof[b] : 0;
            double E_avg = (count_prof[b] > 0) ? E_prof[b] / count_prof[b] : 0;
            printf("  %8.2f  %8d  %8.4f  %10.4e\n",
                   r, count_prof[b], rho_avg, E_avg);
        }

        free(rho_prof); free(E_prof); free(count_prof);
    }

    /* Print separation history summary */
    printf("\n  === Separation Summary ===\n");
    printf("  (Full data in %s)\n", log_fn);

    /* Close log */
    if (log_fp) fclose(log_fp);

    /* Cleanup */
    hash_free(&hg);
    free(acc);
    free(energy_density);
    free(parts);

    return 0;
}
