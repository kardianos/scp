/*  foam_sim.c — V56 stage A: 8-component multivector field on a Voronoi
 *  cell complex.
 *
 *  Field: M[i] for i = 0..7, components of an 8-component subalgebra of
 *  the projective spacetime geometric algebra ℝ(3,1,1) (Lengyel's
 *  "relativistic quaternion"). The metric tensor for the bulk-norm is
 *
 *      g = diag(−,−,−,−, +,+,+,+)
 *
 *  Components 0..3 (rotor part: e₄₁₀, e₄₂₀, e₄₃₀, 𝟙) carry timelike
 *  signs; components 4..7 (spatial: e₁, e₂, e₃, e₃₂₁) carry spacelike.
 *
 *  Equation of motion (per cell c, per component i):
 *
 *      □ M[i] = − λ (|M|²_bulk(c) − v²) g_{ii} M[i](c)
 *
 *  where □ M[i] is the d'Alembertian (spatial Laplacian under
 *  velocity-Verlet), |M|²_bulk = Σ g_{ii} M[i]² is the bulk norm.
 *  This is a Higgs-like quartic potential with vacuum on the
 *  hyperboloid |M|²_bulk = v². See v56/DERIVATION.md.
 *
 *  Build:
 *      gcc -O3 -march=native -fopenmp -o foam_sim foam_sim.c -lzstd -lm
 *
 *  Usage:
 *      ./foam_sim foam_mesh.bin run.cfg
 */

#define _GNU_SOURCE
#define SFA_IMPLEMENTATION
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>
#include <omp.h>
#include <time.h>
#include "../../sfa/format/sfa.h"

/* Field plug-in: defines NCOMP, g_metric, field_names/semantic/component,
 * and FIELD_NAME. The implementations of field_init() and field_forces()
 * are pulled in below, after Sim/Mesh structs are defined.
 *
 * Build switch: -DUSE_FIELD_QBALL → 6-component SO(3)×SO(3) Q-ball plug-in.
 * Default: 8-component PGA relativistic-quaternion field. */
#ifdef USE_FIELD_QBALL
#  include "field_qball.h"
#else
#  include "field_pga.h"
#endif

#define PI 3.14159265358979323846

/* =====================================================================
 *  Mesh structures (unchanged from v55)
 * ===================================================================== */

typedef struct {
    uint32_t a, b;
    double area;
    double dx, dy, dz;
    double dist;
    double nx, ny, nz;
} Face;

typedef struct {
    double L;
    uint32_t N_cells;
    uint32_t N_faces;
    double *cell_x, *cell_y, *cell_z;
    double *cell_vol;
    Face *faces;
    uint32_t *cell_face_off;
    uint32_t *cell_face_idx;
    double dx_min;
} Mesh;

/* Morton reorder helpers — unchanged from v55 */
static inline uint64_t spread_bits_3d(uint32_t v) {
    uint64_t r = v & 0x1FFFFFULL;
    r = (r | (r << 32)) & 0x001F00000000FFFFULL;
    r = (r | (r << 16)) & 0x001F0000FF0000FFULL;
    r = (r | (r << 8))  & 0x100F00F00F00F00FULL;
    r = (r | (r << 4))  & 0x10C30C30C30C30C3ULL;
    r = (r | (r << 2))  & 0x1249249249249249ULL;
    return r;
}

static inline uint64_t morton_code_3d(uint32_t x, uint32_t y, uint32_t z) {
    return spread_bits_3d(x) | (spread_bits_3d(y) << 1) | (spread_bits_3d(z) << 2);
}

typedef struct { uint64_t key; uint32_t idx; } SortPair;

static int compare_sortpair(const void *a, const void *b) {
    const SortPair *p = a, *q = b;
    if (p->key < q->key) return -1;
    if (p->key > q->key) return  1;
    return 0;
}

static void reorder_cells_morton(Mesh *m) {
    uint32_t N = m->N_cells;
    double L = m->L;
    double scale = (double)(1U << 20) / (2.0 * L);

    SortPair *pairs = malloc(sizeof(SortPair) * N);
    for (uint32_t i = 0; i < N; i++) {
        uint32_t x = (uint32_t)((m->cell_x[i] + L) * scale);
        uint32_t y = (uint32_t)((m->cell_y[i] + L) * scale);
        uint32_t z = (uint32_t)((m->cell_z[i] + L) * scale);
        if (x >= (1U<<20)) x = (1U<<20)-1;
        if (y >= (1U<<20)) y = (1U<<20)-1;
        if (z >= (1U<<20)) z = (1U<<20)-1;
        pairs[i].key = morton_code_3d(x, y, z);
        pairs[i].idx = i;
    }
    qsort(pairs, N, sizeof(SortPair), compare_sortpair);

    uint32_t *inv = malloc(sizeof(uint32_t) * N);
    double *tx = malloc(sizeof(double) * N);
    double *ty = malloc(sizeof(double) * N);
    double *tz = malloc(sizeof(double) * N);
    double *tv = malloc(sizeof(double) * N);
    for (uint32_t new_id = 0; new_id < N; new_id++) {
        uint32_t old = pairs[new_id].idx;
        inv[old] = new_id;
        tx[new_id] = m->cell_x[old];
        ty[new_id] = m->cell_y[old];
        tz[new_id] = m->cell_z[old];
        tv[new_id] = m->cell_vol[old];
    }
    memcpy(m->cell_x, tx, sizeof(double) * N);
    memcpy(m->cell_y, ty, sizeof(double) * N);
    memcpy(m->cell_z, tz, sizeof(double) * N);
    memcpy(m->cell_vol, tv, sizeof(double) * N);
    free(tx); free(ty); free(tz); free(tv);
    free(pairs);

    #pragma omp parallel for schedule(static)
    for (uint32_t f = 0; f < m->N_faces; f++) {
        m->faces[f].a = inv[m->faces[f].a];
        m->faces[f].b = inv[m->faces[f].b];
    }
    free(inv);

    SortPair *fpairs = malloc(sizeof(SortPair) * m->N_faces);
    for (uint32_t f = 0; f < m->N_faces; f++) {
        uint32_t a = m->faces[f].a, b = m->faces[f].b;
        uint32_t mn = a < b ? a : b, mx = a < b ? b : a;
        fpairs[f].key = ((uint64_t)mn << 32) | (uint64_t)mx;
        fpairs[f].idx = f;
    }
    qsort(fpairs, m->N_faces, sizeof(SortPair), compare_sortpair);
    Face *new_faces = malloc(sizeof(Face) * m->N_faces);
    for (uint32_t f = 0; f < m->N_faces; f++)
        new_faces[f] = m->faces[fpairs[f].idx];
    free(m->faces);
    m->faces = new_faces;
    free(fpairs);

    uint32_t *new_off = calloc(N + 1, sizeof(uint32_t));
    for (uint32_t f = 0; f < m->N_faces; f++) {
        new_off[m->faces[f].a + 1]++;
        new_off[m->faces[f].b + 1]++;
    }
    for (uint32_t i = 1; i <= N; i++) new_off[i] += new_off[i-1];
    uint32_t total = new_off[N];
    uint32_t *new_idx = malloc(sizeof(uint32_t) * total);
    uint32_t *cursor = malloc(sizeof(uint32_t) * N);
    memcpy(cursor, new_off, sizeof(uint32_t) * N);
    for (uint32_t f = 0; f < m->N_faces; f++) {
        new_idx[cursor[m->faces[f].a]++] = f;
        new_idx[cursor[m->faces[f].b]++] = f;
    }
    free(cursor);
    free(m->cell_face_off); free(m->cell_face_idx);
    m->cell_face_off = new_off;
    m->cell_face_idx = new_idx;

    printf("[mesh] Morton-reordered %u cells + %u faces\n", N, m->N_faces);
}

static Mesh *mesh_read(const char *path) {
    FILE *fp = fopen(path, "rb");
    if (!fp) { fprintf(stderr, "FATAL: cannot open mesh '%s'\n", path); exit(1); }
    char magic[4];
    if (fread(magic, 1, 4, fp) != 4 || memcmp(magic, "FOAM", 4) != 0) {
        fprintf(stderr, "FATAL: bad magic\n"); exit(1);
    }
    uint32_t version;
    fread(&version, 4, 1, fp);
    if (version != 1) { fprintf(stderr, "FATAL: unsupported mesh version %u\n", version); exit(1); }

    Mesh *m = calloc(1, sizeof(Mesh));
    fread(&m->L, 8, 1, fp);
    fread(&m->N_cells, 4, 1, fp);
    fread(&m->N_faces, 4, 1, fp);
    uint64_t reserved; fread(&reserved, 8, 1, fp);

    printf("[mesh] L=%.3f N_cells=%u N_faces=%u\n",
           m->L, m->N_cells, m->N_faces);

    m->cell_x   = malloc(sizeof(double) * m->N_cells);
    m->cell_y   = malloc(sizeof(double) * m->N_cells);
    m->cell_z   = malloc(sizeof(double) * m->N_cells);
    m->cell_vol = malloc(sizeof(double) * m->N_cells);

    for (uint32_t i = 0; i < m->N_cells; i++) {
        double rec[4];
        fread(rec, 8, 4, fp);
        m->cell_x[i] = rec[0];
        m->cell_y[i] = rec[1];
        m->cell_z[i] = rec[2];
        m->cell_vol[i] = rec[3];
    }

    m->faces = malloc(sizeof(Face) * m->N_faces);
    m->dx_min = 1e9;
    for (uint32_t i = 0; i < m->N_faces; i++) {
        Face *f = &m->faces[i];
        fread(&f->a, 4, 1, fp);
        fread(&f->b, 4, 1, fp);
        fread(&f->area, 8, 1, fp);
        fread(&f->dx, 8, 1, fp);
        fread(&f->dy, 8, 1, fp);
        fread(&f->dz, 8, 1, fp);
        double d = sqrt(f->dx*f->dx + f->dy*f->dy + f->dz*f->dz);
        f->dist = d;
        f->nx = f->dx / d;
        f->ny = f->dy / d;
        f->nz = f->dz / d;
        if (d < m->dx_min) m->dx_min = d;
    }

    m->cell_face_off = malloc(sizeof(uint32_t) * (m->N_cells + 1));
    fread(m->cell_face_off, 4, m->N_cells + 1, fp);
    uint32_t total_idx = m->cell_face_off[m->N_cells];
    m->cell_face_idx = malloc(sizeof(uint32_t) * total_idx);
    fread(m->cell_face_idx, 4, total_idx, fp);

    fclose(fp);

    double vmin = m->cell_vol[0], vmax = m->cell_vol[0], vsum = 0;
    for (uint32_t i = 0; i < m->N_cells; i++) {
        if (m->cell_vol[i] < vmin) vmin = m->cell_vol[i];
        if (m->cell_vol[i] > vmax) vmax = m->cell_vol[i];
        vsum += m->cell_vol[i];
    }
    double V_box = pow(2.0 * m->L, 3.0);
    double avg_deg = (double)(2 * m->N_faces) / m->N_cells;
    printf("[mesh] dx_min=%.3f V/cell mean=%.3f range=[%.3f,%.3f] V_total=%.1f (box=%.1f) avg_deg=%.1f\n",
           m->dx_min, vsum / m->N_cells, vmin, vmax, vsum, V_box, avg_deg);

    reorder_cells_morton(m);
    return m;
}

static void mesh_free(Mesh *m) {
    free(m->cell_x); free(m->cell_y); free(m->cell_z); free(m->cell_vol);
    free(m->faces);
    free(m->cell_face_off); free(m->cell_face_idx);
    free(m);
}

/* =====================================================================
 *  Sim state — V56: 8-component multivector field per cell
 * ===================================================================== */

typedef struct {
    /* Physics */
    double mass2;       /* squared mass scale (used in linearised regime) */
    double lambda;      /* quartic coupling (Higgs term: λ(|M|²−v²)²) */
    double v2;          /* vacuum norm-squared, |M|²_bulk = v² */
    double sigma_e2;    /* rotor σ-constraint strength: (e2/4)(|q|²−1)² */
    double skyrme_c4;   /* Skyrme L_4 coupling: c4 = 1/(2 e²) (standard
                         * Skyrme convention has L_4 = (1/32 e²) Tr([L_μ,L_ν]²)
                         * which expands to (c4/2) [(Tr g)² − ‖g‖²_F]) */
    int skyrme_project; /* 1 → hard-project rotor onto |q|=1 each step
                         * (with tangent velocity projection). 0 → off
                         * (rely on soft σ-constraint instead). */
    double damping;     /* velocity damping rate γ (units of 1/time).
                         * After each Verlet step, M_vel *= (1 − γ·dt).
                         * Use to absorb Derrick relaxation kinetic energy
                         * for an oversized seed without dynamic overshoot.
                         * Set to 0 for pure energy-conserving dynamics. */
    double damping_T;   /* time after which damping is disabled (s).
                         * Useful for a damped initial relaxation followed
                         * by undamped continuation. 0 → damping always on. */
    int    damping_sector_size; /* if > 0, tangent-only damping per sector
                         * of this size. Decomposes velocity into radial
                         * (along M_sector) and tangent (⊥ M_sector), damps
                         * ONLY radial — preserves Noether charge of a
                         * rotating soliton (Q-ball). Use 3 for the 6-comp
                         * Q-ball (two 3-vec sectors). 0 → full damping. */
    double m_pi2;       /* pion mass squared. V_π = m_π² (1 − q₀)
                         * with q₀ = M[3]. Sets a finite preferred
                         * Skyrmion size (Compton wavelength of the
                         * pion) — without this term the soliton can
                         * shrink past lattice resolution and unwind. */

    /* Seed knobs (interpretation depends on init type) */
    double A;           /* perturbation amplitude */
    double R_seed;      /* localisation radius */
    double omega_seed;  /* oscillation frequency for time-dependent seeds */
    double omega_theta; /* (Q-ball) rotation frequency for the θ sector.
                         * 0 → θ stays static (single-charge Q-ball).
                         * >0 → θ rotates in (θ_y, θ_z) plane at this rate,
                         *      giving a two-charge Q-ball with both
                         *      SO(3)_φ and SO(3)_θ Noether currents active. */

    /* Braid-seed knobs (mirror v55) */
    double A_bg;        /* background carrier amplitude */
    double R_tube;      /* envelope radius (Gaussian) */
    double ellip;       /* cross-section ellipticity */
    double delta[3];    /* phase offsets for the 3 carriers */

    /* Time stepping */
    double T, dt_factor, dt;

    /* I/O */
    char init[64];
    char output[256];
    char diag_file[256];
    double snap_dt, diag_dt;
    int voxel_N;
    int sfa_output;
    int cell_native;
    int cell_iframe_interval;
    int cell_model_interval;
    double cell_delta_tol;
    double cell_omega;
} Config;

typedef struct {
    Mesh *m;
    Config c;
    double *M[NCOMP];           /* M[i][cell] — multivector components */
    double *M_vel[NCOMP];
    double *M_acc[NCOMP];
    double *grad_M[NCOMP][3];   /* ∂_b M[i] per cell */
    double *bulk_norm;          /* |M|²_bulk per cell, cached for diagnostics */
    /* Skyrme L_4 scratch: J^{a,b}(c) on rotor only. Allocated lazily. */
    double *Jx[NCOMP];
    double *Jy[NCOMP];
    double *Jz[NCOMP];
    double *winding;            /* topological charge density per cell (diag) */
} Sim;

static double parse_double(const char *v) { return strtod(v, NULL); }

static void config_default(Config *c) {
    c->mass2 = 0.0;             /* default: pure Higgs dynamics (no bare mass) */
    c->lambda = 5.0;
    c->v2 = 0.25;               /* v = 0.5 → |M|²_bulk = 0.25 at vacuum */
    c->sigma_e2 = 0.0;          /* default: no σ-constraint */
    c->skyrme_c4 = 0.0;         /* default: no Skyrme L_4 term */
    c->skyrme_project = 0;      /* default: no hard S³ projection */
    c->damping = 0.0;           /* default: no velocity damping */
    c->damping_T = 0.0;
    c->damping_sector_size = 0; /* default: full damping (no sector split) */
    c->m_pi2 = 0.0;             /* default: massless pion (no V_π term) */

    c->A = 0.3;
    c->R_seed = 4.0;
    c->omega_seed = 1.0;
    c->omega_theta = 0.0;       /* default: θ sector static */

    /* v55-mirroring defaults */
    c->A_bg = 0.1;
    c->R_tube = 3.0;
    c->ellip = 0.3325;
    c->delta[0] = 0.0;
    c->delta[1] = 3.0005;
    c->delta[2] = 4.4325;

    c->T = 50.0;
    c->dt_factor = 0.05;
    c->snap_dt = 1.0;
    c->diag_dt = 1.0;

    c->voxel_N = 128;
    c->sfa_output = 0;
    c->cell_native = 1;         /* v56 default: cell-native (bigger fields) */
    c->cell_iframe_interval = 0;
    c->cell_model_interval = 0;
    c->cell_delta_tol = 0.01;
    c->cell_omega = 1.5;

    strcpy(c->init, "vacuum");
    strcpy(c->output, "v56_run.sfa");
    strcpy(c->diag_file, "v56_run_diag.tsv");
}

static void config_load(Config *c, const char *path) {
    FILE *fp = fopen(path, "r");
    if (!fp) { fprintf(stderr, "FATAL: cannot open config '%s'\n", path); exit(1); }
    char line[1024];
    while (fgets(line, sizeof(line), fp)) {
        char *p = line;
        while (*p == ' ' || *p == '\t') p++;
        if (*p == '#' || *p == '\n' || *p == '\0') continue;
        char key[64], val[256];
        if (sscanf(p, "%63[^=]=%255[^\n]", key, val) != 2) continue;
        char *e = key + strlen(key) - 1;
        while (e > key && (*e == ' ' || *e == '\t')) *e-- = '\0';
        char *vs = val;
        while (*vs == ' ' || *vs == '\t') vs++;
        e = vs + strlen(vs) - 1;
        while (e > vs && (*e == ' ' || *e == '\t' || *e == '\n')) *e-- = '\0';

        if      (!strcmp(key, "m"))       { double m = parse_double(vs); c->mass2 = m*m; }
        else if (!strcmp(key, "mass2"))     c->mass2 = parse_double(vs);
        else if (!strcmp(key, "lambda"))    c->lambda = parse_double(vs);
        else if (!strcmp(key, "v"))       { double v = parse_double(vs); c->v2 = v*v; }
        else if (!strcmp(key, "v2"))        c->v2 = parse_double(vs);
        else if (!strcmp(key, "sigma_e2"))  c->sigma_e2 = parse_double(vs);
        else if (!strcmp(key, "skyrme_e2")) c->sigma_e2 = parse_double(vs);   /* legacy alias */
        else if (!strcmp(key, "skyrme_c4")) c->skyrme_c4 = parse_double(vs);
        else if (!strcmp(key, "skyrme_project")) c->skyrme_project = atoi(vs);
        else if (!strcmp(key, "damping"))   c->damping = parse_double(vs);
        else if (!strcmp(key, "damping_T")) c->damping_T = parse_double(vs);
        else if (!strcmp(key, "damping_sector_size")) c->damping_sector_size = atoi(vs);
        else if (!strcmp(key, "m_pi"))      { double m = parse_double(vs); c->m_pi2 = m*m; }
        else if (!strcmp(key, "m_pi2"))      c->m_pi2 = parse_double(vs);
        else if (!strcmp(key, "A"))         c->A = parse_double(vs);
        else if (!strcmp(key, "R_seed"))    c->R_seed = parse_double(vs);
        else if (!strcmp(key, "omega_seed")) c->omega_seed = parse_double(vs);
        else if (!strcmp(key, "omega_theta")) c->omega_theta = parse_double(vs);
        else if (!strcmp(key, "A_bg"))      c->A_bg = parse_double(vs);
        else if (!strcmp(key, "R_tube"))    c->R_tube = parse_double(vs);
        else if (!strcmp(key, "ellip"))     c->ellip = parse_double(vs);
        else if (!strcmp(key, "delta")) {
            sscanf(vs, "%lf,%lf,%lf", &c->delta[0], &c->delta[1], &c->delta[2]);
        }
        else if (!strcmp(key, "T"))         c->T = parse_double(vs);
        else if (!strcmp(key, "dt_factor")) c->dt_factor = parse_double(vs);
        else if (!strcmp(key, "snap_dt"))   c->snap_dt = parse_double(vs);
        else if (!strcmp(key, "diag_dt"))   c->diag_dt = parse_double(vs);
        else if (!strcmp(key, "voxel_N"))   c->voxel_N = atoi(vs);
        else if (!strcmp(key, "sfa_output")) c->sfa_output = atoi(vs);
        else if (!strcmp(key, "cell_native")) c->cell_native = atoi(vs);
        else if (!strcmp(key, "cell_iframe_interval")) c->cell_iframe_interval = atoi(vs);
        else if (!strcmp(key, "cell_model_interval")) c->cell_model_interval = atoi(vs);
        else if (!strcmp(key, "cell_delta_tol")) c->cell_delta_tol = parse_double(vs);
        else if (!strcmp(key, "cell_omega")) c->cell_omega = parse_double(vs);
        else if (!strcmp(key, "init"))      strncpy(c->init, vs, 63);
        else if (!strcmp(key, "output"))    strncpy(c->output, vs, 255);
        else if (!strcmp(key, "diag_file")) strncpy(c->diag_file, vs, 255);
    }
    fclose(fp);
}

static void sim_alloc(Sim *s) {
    uint32_t N = s->m->N_cells;
    for (int k = 0; k < NCOMP; k++) {
        s->M[k]     = calloc(N, sizeof(double));
        s->M_vel[k] = calloc(N, sizeof(double));
        s->M_acc[k] = calloc(N, sizeof(double));
        for (int b = 0; b < 3; b++)
            s->grad_M[k][b] = calloc(N, sizeof(double));
    }
    s->bulk_norm = calloc(N, sizeof(double));
    /* Skyrme J scratch is allocated lazily inside compute_forces only
     * when skyrme_c4 > 0 — saves 4*3*8 = 96 bytes/cell on non-Skyrme runs. */
    for (int k = 0; k < NCOMP; k++) {
        s->Jx[k] = NULL; s->Jy[k] = NULL; s->Jz[k] = NULL;
    }
    s->winding = calloc(N, sizeof(double));
}

static void sim_free(Sim *s) {
    for (int k = 0; k < NCOMP; k++) {
        free(s->M[k]); free(s->M_vel[k]); free(s->M_acc[k]);
        for (int b = 0; b < 3; b++) free(s->grad_M[k][b]);
        if (s->Jx[k]) { free(s->Jx[k]); free(s->Jy[k]); free(s->Jz[k]); }
    }
    free(s->bulk_norm);
    free(s->winding);
}

/* =====================================================================
 *  Field-specific code (init seeds + self-coupling) is in field_pga.h.
 *  Pull in its implementation here, now that Sim/Mesh/Config are
 *  defined. To swap in a different field theory, replace field_pga.h
 *  with another header conforming to the same API contract.
 * ===================================================================== */
#define FIELD_IMPL
#ifdef USE_FIELD_QBALL
#  include "field_qball.h"
#else
#  include "field_pga.h"
#endif

/* =====================================================================
 *  Finite-volume operators
 * ===================================================================== */

/* Combined gradient + Laplacian for all NCOMP fields in a single cell pass.
 * Output layout: grad_x[k][cell], grad_y[k][cell], grad_z[k][cell], lap[k][cell]. */
static void compute_grads_and_lap_all(Sim *s,
        double *grad_x[NCOMP], double *grad_y[NCOMP], double *grad_z[NCOMP],
        double *lap[NCOMP]) {
    Mesh *m = s->m;
    uint32_t N = m->N_cells;
    double *fields[NCOMP];
    for (int k = 0; k < NCOMP; k++) fields[k] = s->M[k];

    #pragma omp parallel for schedule(static)
    for (uint32_t c = 0; c < N; c++) {
        double inv_V = 1.0 / m->cell_vol[c];
        double agx[NCOMP] = {0}, agy[NCOMP] = {0}, agz[NCOMP] = {0}, alp[NCOMP] = {0};
        double fc[NCOMP];
        for (int k = 0; k < NCOMP; k++) fc[k] = fields[k][c];

        uint32_t off0 = m->cell_face_off[c];
        uint32_t off1 = m->cell_face_off[c+1];
        for (uint32_t kk = off0; kk < off1; kk++) {
            uint32_t fid = m->cell_face_idx[kk];
            const Face *F = &m->faces[fid];
            uint32_t nb;
            double sign;
            if (F->a == c) { nb = F->b; sign = +1.0; }
            else           { nb = F->a; sign = -1.0; }
            double area = F->area;
            double d    = F->dist;
            double nx   = sign * F->nx;
            double ny   = sign * F->ny;
            double nz   = sign * F->nz;
            double inv_d = 1.0 / d;

            for (int k = 0; k < NCOMP; k++) {
                double diff = fields[k][nb] - fc[k];
                agx[k] += area * nx * diff;
                agy[k] += area * ny * diff;
                agz[k] += area * nz * diff;
                alp[k] += area * diff * inv_d;
            }
        }
        for (int k = 0; k < NCOMP; k++) {
            grad_x[k][c] = 0.5 * agx[k] * inv_V;
            grad_y[k][c] = 0.5 * agy[k] * inv_V;
            grad_z[k][c] = 0.5 * agz[k] * inv_V;
            lap[k][c]    = alp[k] * inv_V;
        }
    }
}

/* =====================================================================
 *  Generic force step: kernel computes Laplacian (always the same on
 *  the foam mesh), then delegates the field-specific potential / self-
 *  coupling to field_forces() defined in the field plug-in. If the
 *  plug-in declares a Skyrme L_4 flux (field_skyrme_J returns 1), the
 *  kernel applies its discrete divergence via a second face pass.
 * ===================================================================== */
static void compute_forces(Sim *s) {
    Mesh *m = s->m;
    uint32_t N = m->N_cells;

    /* Compute Laplacian of all NCOMP components in one pass. */
    static double *lap_M[NCOMP] = {NULL};
    if (lap_M[0] == NULL) {
        for (int k = 0; k < NCOMP; k++) lap_M[k] = malloc(sizeof(double) * N);
    }
    double *gx_all[NCOMP], *gy_all[NCOMP], *gz_all[NCOMP], *lap_all[NCOMP];
    for (int k = 0; k < NCOMP; k++) {
        gx_all[k] = s->grad_M[k][0];
        gy_all[k] = s->grad_M[k][1];
        gz_all[k] = s->grad_M[k][2];
        lap_all[k] = lap_M[k];
    }
    compute_grads_and_lap_all(s, gx_all, gy_all, gz_all, lap_all);

    /* Local part of forces (lap + potential + σ-constraint). */
    double *acc[NCOMP];
    for (int k = 0; k < NCOMP; k++) acc[k] = s->M_acc[k];
    field_forces(s, (double *const *)lap_all, acc);

    /* Skyrme L_4 contribution: compute J^{a,b}(c) on cells, then add
     * ∂_b J^{a,b} to acc[a] via face divergence. */
    if (s->c.skyrme_c4 > 0.0) {
        if (s->Jx[0] == NULL) {
            for (int a = 0; a < N_SKYRME; a++) {
                s->Jx[a] = calloc(N, sizeof(double));
                s->Jy[a] = calloc(N, sizeof(double));
                s->Jz[a] = calloc(N, sizeof(double));
            }
        }
        int active = field_skyrme_J(s, gx_all, gy_all, gz_all,
                                    s->Jx, s->Jy, s->Jz);
        if (active) {
            /* Discrete divergence of J on faces:
             *   acc[a](c) += (1/V_c) Σ_{f∈c} A_f n̂_f^{out} · (J(a)+J(nb))/2
             * Face stores normal pointing from a→b, so for cell c we
             * pick sign = +1 if c==a, −1 if c==b. */
            #pragma omp parallel for schedule(static)
            for (uint32_t c = 0; c < N; c++) {
                double inv_V = 1.0 / m->cell_vol[c];
                double div_acc[N_SKYRME] = {0};
                uint32_t off0 = m->cell_face_off[c];
                uint32_t off1 = m->cell_face_off[c+1];
                for (uint32_t kk = off0; kk < off1; kk++) {
                    uint32_t fid = m->cell_face_idx[kk];
                    const Face *F = &m->faces[fid];
                    uint32_t nb;
                    double sign;
                    if (F->a == c) { nb = F->b; sign = +1.0; }
                    else           { nb = F->a; sign = -1.0; }
                    double area = F->area;
                    double nx = sign * F->nx;
                    double ny = sign * F->ny;
                    double nz = sign * F->nz;
                    for (int a = 0; a < N_SKYRME; a++) {
                        double Jxf = 0.5 * (s->Jx[a][c] + s->Jx[a][nb]);
                        double Jyf = 0.5 * (s->Jy[a][c] + s->Jy[a][nb]);
                        double Jzf = 0.5 * (s->Jz[a][c] + s->Jz[a][nb]);
                        div_acc[a] += area * (nx * Jxf + ny * Jyf + nz * Jzf);
                    }
                }
                for (int a = 0; a < N_SKYRME; a++)
                    acc[a][c] += div_acc[a] * inv_V;
            }
        }
    }
}

/* =====================================================================
 *  Hard S³ projection (option 1 from SCAN_RESULTS.md)
 *
 *  Constrains M[0..N_SKYRME-1] to |q|=1 and tangent-projects the
 *  velocity so that ⟨q, q̇⟩ = 0. Applied after each Verlet step.
 *  Drift envelope is bounded by O(dt²·|a|) per step — same order as
 *  the integrator itself. Eliminates the bulk-Higgs scalar branch
 *  competition that caused the soft σ-constraint to fail in the
 *  parameter scans (see v56/SCAN_RESULTS.md).
 * ===================================================================== */
static void project_rotor_S3(Sim *s) {
    uint32_t N = s->m->N_cells;
    #pragma omp parallel for schedule(static)
    for (uint32_t i = 0; i < N; i++) {
        double q2 = 0;
        for (int k = 0; k < N_SKYRME; k++) q2 += s->M[k][i] * s->M[k][i];
        if (q2 < 1e-30) continue;
        double inv_qn = 1.0 / sqrt(q2);
        for (int k = 0; k < N_SKYRME; k++) s->M[k][i] *= inv_qn;
        /* Tangent project velocity: v ← v − ⟨q, v⟩ q  (with the new |q|=1 q). */
        double qdotv = 0;
        for (int k = 0; k < N_SKYRME; k++)
            qdotv += s->M[k][i] * s->M_vel[k][i];
        for (int k = 0; k < N_SKYRME; k++)
            s->M_vel[k][i] -= qdotv * s->M[k][i];
    }
}

/* =====================================================================
 *  Velocity Verlet (generic — operates on M[k] components)
 * ===================================================================== */
static void verlet_step(Sim *s, double t_end) {
    uint32_t N = s->m->N_cells;
    double hdt = 0.5 * s->c.dt, dt = s->c.dt;

    #pragma omp parallel for schedule(static)
    for (uint32_t i = 0; i < N; i++) {
        for (int k = 0; k < NCOMP; k++)
            s->M_vel[k][i] += hdt * s->M_acc[k][i];
    }
    #pragma omp parallel for schedule(static)
    for (uint32_t i = 0; i < N; i++) {
        for (int k = 0; k < NCOMP; k++)
            s->M[k][i] += dt * s->M_vel[k][i];
    }
    if (s->c.skyrme_project) project_rotor_S3(s);
    compute_forces(s);
    #pragma omp parallel for schedule(static)
    for (uint32_t i = 0; i < N; i++) {
        for (int k = 0; k < NCOMP; k++)
            s->M_vel[k][i] += hdt * s->M_acc[k][i];
    }
    /* Velocity damping. Applied after the second half-kick so projection's
     * tangent correction still sees the damped velocity.
     *
     * Two modes:
     *   damping_sector_size = 0 → full damping: M_vel *= (1−γ·dt) for all k.
     *   damping_sector_size = S > 0 → tangent-only damping per sector of
     *     size S. For each cell and each sector of S consecutive components,
     *     decompose v into v_radial = (v·M̂)M̂ (along M_sector) and v_tan =
     *     v − v_radial. Then damp only v_radial:
     *       v ← v_tan + (1−γ·dt)·v_radial
     *     This preserves the Noether charge of a rotating soliton (|M_sec|=
     *     const, motion purely tangent), while absorbing the radial
     *     breathing mode excited by lattice-discretization mismatch. */
    if (s->c.damping > 0.0 &&
        (s->c.damping_T <= 0.0 || t_end < s->c.damping_T)) {
        double damp = 1.0 - s->c.damping * dt;
        if (damp < 0.0) damp = 0.0;
        int S = s->c.damping_sector_size;
        if (S <= 0) {
            #pragma omp parallel for schedule(static)
            for (uint32_t i = 0; i < N; i++) {
                for (int k = 0; k < NCOMP; k++)
                    s->M_vel[k][i] *= damp;
            }
        } else {
            double gamma_dt = s->c.damping * dt;
            #pragma omp parallel for schedule(static)
            for (uint32_t i = 0; i < N; i++) {
                for (int sec = 0; sec + S <= NCOMP; sec += S) {
                    double mag2 = 0;
                    for (int k = 0; k < S; k++)
                        mag2 += s->M[sec+k][i] * s->M[sec+k][i];
                    if (mag2 < 1e-30) continue;
                    double v_dot_m = 0;
                    for (int k = 0; k < S; k++)
                        v_dot_m += s->M_vel[sec+k][i] * s->M[sec+k][i];
                    double damp_amt = gamma_dt * v_dot_m / mag2;
                    /* Subtract gamma_dt · (radial component of v) */
                    for (int k = 0; k < S; k++)
                        s->M_vel[sec+k][i] -= damp_amt * s->M[sec+k][i];
                }
            }
        }
    }
    if (s->c.skyrme_project) project_rotor_S3(s);
}

/* =====================================================================
 *  Diagnostics — V56: bulk-norm stats + energy decomposition
 * ===================================================================== */

typedef struct {
    double E_kin;       /* (1/2) Σ V_c (M_vel)² */
    double E_grad;      /* (1/2) Σ_face A_f/d_f (ΔM)²  — discrete face form */
    double E_mass;      /* (1/2) m² Σ V_c M[k]² */
    double E_quartic;   /* (λ/4) Σ V_c (|M|²_bulk − v²)² */
    double E_sigma;     /* (sigma_e2/4) Σ V_c (|q|² − 1)² */
    double E_skyrme;    /* (skyrme_c4/2) Σ V_c [(Tr g)² − ‖g‖²_F]   (real L_4) */
    double E_pion;      /* m_π² Σ V_c (1 − M[3])   (pion-mass term) */
    double E_total;
    double bulk_min, bulk_mean, bulk_max;
    double M_max[NCOMP];
    double rotor_rms, spatial_rms;   /* RMS of M[0..3] vs M[4..7] */
    double B_winding;   /* topological charge ∫ ρ_B d³x = (1/2π²)∫det4(q,∂q,∂q,∂q) */
    double q_norm_min, q_norm_mean, q_norm_max;  /* |q|² stats */
} Diag;

static void compute_diag(Sim *s, Diag *d) {
    Mesh *m = s->m;
    Config *c = &s->c;
    uint32_t N = m->N_cells;

    double E_kin = 0, E_mass = 0, E_quartic = 0, E_sigma = 0, E_skyrme = 0, E_pion = 0;
    double bulk_sum = 0, bulk_min = +1e30, bulk_max = -1e30;
    double rotor_sum = 0, spatial_sum = 0;
    double V_total = 0;
    double M_max_local[NCOMP] = {0};
    double B_winding = 0;
    double q_sum = 0, q_min = +1e30, q_max = -1e30;

    /* Cell-summed pieces. Skyrme L_4 and winding both use grad_M which
     * was computed during the most recent verlet force step. */
    double *gx[NCOMP], *gy[NCOMP], *gz[NCOMP];
    for (int k = 0; k < NCOMP; k++) {
        gx[k] = s->grad_M[k][0];
        gy[k] = s->grad_M[k][1];
        gz[k] = s->grad_M[k][2];
    }

    #pragma omp parallel for reduction(+:E_kin,E_mass,E_quartic,E_sigma,E_skyrme,E_pion,bulk_sum,rotor_sum,spatial_sum,V_total,B_winding,q_sum) \
                              reduction(min:bulk_min,q_min) reduction(max:bulk_max,q_max) \
                              schedule(static)
    for (uint32_t i = 0; i < N; i++) {
        double V_c = m->cell_vol[i];
        V_total += V_c;
        double Mb = 0;
        for (int k = 0; k < NCOMP; k++) {
            double v = s->M[k][i], vv = s->M_vel[k][i];
            E_kin  += 0.5 * vv * vv * V_c;
            E_mass += 0.5 * c->mass2 * v * v * V_c;
            Mb += g_metric[k] * v * v;
            if (k < N_SKYRME) rotor_sum   += v * v;
            else              spatial_sum += v * v;
        }
        bulk_sum += Mb * V_c;
        if (Mb < bulk_min) bulk_min = Mb;
        if (Mb > bulk_max) bulk_max = Mb;
        double dMb = Mb - c->v2;
        E_quartic += 0.25 * c->lambda * dMb * dMb * V_c;
        double q2 = 0.0;
        for (int j = 0; j < N_SKYRME; j++) q2 += s->M[j][i] * s->M[j][i];
        double dq = q2 - 1.0;
        E_sigma += 0.25 * c->sigma_e2 * dq * dq * V_c;
        q_sum += q2 * V_c;
        if (q2 < q_min) q_min = q2;
        if (q2 > q_max) q_max = q2;

        /* Skyrme L_4 energy density:
         *   e_4 = (c4/2) [(Tr g)² − ‖g‖²_F]
         * with g_ij = Σ_a (∂_i q^a)(∂_j q^a) over rotor a = 0..N_SKYRME-1. */
        double gxx=0,gyy=0,gzz=0,gxy=0,gxz=0,gyz=0;
        for (int a = 0; a < N_SKYRME; a++) {
            double dxa = gx[a][i], dya = gy[a][i], dza = gz[a][i];
            gxx += dxa*dxa; gyy += dya*dya; gzz += dza*dza;
            gxy += dxa*dya; gxz += dxa*dza; gyz += dya*dza;
        }
        double trg = gxx + gyy + gzz;
        double frob2 = gxx*gxx + gyy*gyy + gzz*gzz
                     + 2.0*(gxy*gxy + gxz*gxz + gyz*gyz);
        double tr2 = trg * trg;
        double e4_density = 0.5 * c->skyrme_c4 * (tr2 - frob2);
        E_skyrme += e4_density * V_c;
        E_pion += c->m_pi2 * (1.0 - s->M[3][i]) * V_c;

        /* Winding-number density: ρ_B = (1/2π²) det4(q, ∂_x q, ∂_y q, ∂_z q)
         * over rotor 4-vector q = (M[0], M[1], M[2], M[3]). */
        double q0 = s->M[0][i], q1 = s->M[1][i], q2c = s->M[2][i], q3 = s->M[3][i];
        double dx0=gx[0][i], dx1=gx[1][i], dx2=gx[2][i], dx3=gx[3][i];
        double dy0=gy[0][i], dy1=gy[1][i], dy2=gy[2][i], dy3=gy[3][i];
        double dz0=gz[0][i], dz1=gz[1][i], dz2=gz[2][i], dz3=gz[3][i];
        /* 4×4 determinant via cofactor expansion along first row. */
        double m00=q0,  m01=q1,  m02=q2c, m03=q3;
        double m10=dx0, m11=dx1, m12=dx2, m13=dx3;
        double m20=dy0, m21=dy1, m22=dy2, m23=dy3;
        double m30=dz0, m31=dz1, m32=dz2, m33=dz3;
        double c0 =  m11*(m22*m33 - m23*m32) - m12*(m21*m33 - m23*m31) + m13*(m21*m32 - m22*m31);
        double c1 = -m10*(m22*m33 - m23*m32) + m12*(m20*m33 - m23*m30) - m13*(m20*m32 - m22*m30);
        double c2 =  m10*(m21*m33 - m23*m31) - m11*(m20*m33 - m23*m30) + m13*(m20*m31 - m21*m30);
        double c3 = -m10*(m21*m32 - m22*m31) + m11*(m20*m32 - m22*m30) - m12*(m20*m31 - m21*m30);
        double det4 = m00*c0 + m01*c1 + m02*c2 + m03*c3;
        double rho_B = det4 / (2.0 * M_PI * M_PI);
        s->winding[i] = rho_B;
        B_winding += rho_B * V_c;
    }

    /* Per-component max — separate parallel reduction to keep things simple */
    for (int k = 0; k < NCOMP; k++) {
        double mk = 0;
        #pragma omp parallel for reduction(max:mk) schedule(static)
        for (uint32_t i = 0; i < N; i++) {
            double a = fabs(s->M[k][i]);
            if (a > mk) mk = a;
        }
        M_max_local[k] = mk;
    }

    /* Face-summed gradient energy (discretely conserved form) */
    double E_grad = 0;
    #pragma omp parallel for reduction(+:E_grad) schedule(static)
    for (uint32_t fid = 0; fid < m->N_faces; fid++) {
        const Face *F = &m->faces[fid];
        uint32_t a = F->a, b = F->b;
        double w = 0.5 * F->area / F->dist;
        for (int k = 0; k < NCOMP; k++) {
            double dM = s->M[k][b] - s->M[k][a];
            E_grad += w * dM * dM;
        }
    }

    d->E_kin = E_kin;
    d->E_mass = E_mass;
    d->E_quartic = E_quartic;
    d->E_sigma = E_sigma;
    d->E_skyrme = E_skyrme;
    d->E_pion = E_pion;
    d->E_grad = E_grad;
    d->E_total = E_kin + E_mass + E_quartic + E_sigma + E_skyrme + E_pion + E_grad;
    d->bulk_min = bulk_min;
    d->bulk_max = bulk_max;
    d->bulk_mean = bulk_sum / V_total;
    d->rotor_rms = sqrt(rotor_sum / N);
    d->spatial_rms = sqrt(spatial_sum / N);
    d->B_winding = B_winding;
    d->q_norm_min = q_min;
    d->q_norm_max = q_max;
    d->q_norm_mean = q_sum / V_total;
    for (int k = 0; k < NCOMP; k++) d->M_max[k] = M_max_local[k];
}

/* =====================================================================
 *  Foam → voxel resampling (unchanged from v55)
 * ===================================================================== */

typedef struct {
    int Ng;
    double bin_size;
    int *count;
    int *offset;
    int *cells;
    int  voxel_N;
    int *voxel_to_cell;
} ResampleIndex;

static ResampleIndex *resample_index_build(const Mesh *m, int voxel_N) {
    ResampleIndex *r = calloc(1, sizeof(ResampleIndex));
    double box_vol = pow(2.0 * m->L, 3.0);
    double cell_vol = box_vol / m->N_cells;
    r->bin_size = pow(cell_vol, 1.0/3.0) * 1.5;
    r->Ng = (int)ceil((2.0 * m->L) / r->bin_size);
    r->bin_size = (2.0 * m->L) / r->Ng;
    long Ng3 = (long)r->Ng * r->Ng * r->Ng;
    r->count  = calloc(Ng3, sizeof(int));
    r->offset = calloc(Ng3, sizeof(int));

    for (uint32_t i = 0; i < m->N_cells; i++) {
        int bx = (int)((m->cell_x[i] + m->L) / r->bin_size);
        int by = (int)((m->cell_y[i] + m->L) / r->bin_size);
        int bz = (int)((m->cell_z[i] + m->L) / r->bin_size);
        if (bx < 0) bx = 0; if (bx >= r->Ng) bx = r->Ng - 1;
        if (by < 0) by = 0; if (by >= r->Ng) by = r->Ng - 1;
        if (bz < 0) bz = 0; if (bz >= r->Ng) bz = r->Ng - 1;
        r->count[(long)bx * r->Ng * r->Ng + by * r->Ng + bz]++;
    }
    int total = 0;
    for (long b = 0; b < Ng3; b++) { r->offset[b] = total; total += r->count[b]; }
    r->cells = malloc(sizeof(int) * total);
    int *cursor = malloc(sizeof(int) * Ng3);
    memcpy(cursor, r->offset, sizeof(int) * Ng3);
    for (uint32_t i = 0; i < m->N_cells; i++) {
        int bx = (int)((m->cell_x[i] + m->L) / r->bin_size);
        int by = (int)((m->cell_y[i] + m->L) / r->bin_size);
        int bz = (int)((m->cell_z[i] + m->L) / r->bin_size);
        if (bx < 0) bx = 0; if (bx >= r->Ng) bx = r->Ng - 1;
        if (by < 0) by = 0; if (by >= r->Ng) by = r->Ng - 1;
        if (bz < 0) bz = 0; if (bz >= r->Ng) bz = r->Ng - 1;
        r->cells[cursor[(long)bx * r->Ng * r->Ng + by * r->Ng + bz]++] = (int)i;
    }
    free(cursor);

    r->voxel_N = voxel_N;
    long N3 = (long)voxel_N * voxel_N * voxel_N;
    r->voxel_to_cell = malloc(sizeof(int) * N3);
    double L = m->L;
    double dx = 2.0 * L / (voxel_N - 1);
    double S = 2.0 * L;
    int Ng = r->Ng;

    printf("[resample] building voxel→cell index: %d^3 voxels, bin_size=%.3f, Ng=%d\n",
           voxel_N, r->bin_size, Ng);

    #pragma omp parallel for schedule(dynamic, 32)
    for (int i = 0; i < voxel_N; i++) {
        double x = -L + i * dx;
        int bx0 = (int)((x + L) / r->bin_size);
        bx0 = ((bx0 % Ng) + Ng) % Ng;
        for (int j = 0; j < voxel_N; j++) {
            double y = -L + j * dx;
            int by0 = (int)((y + L) / r->bin_size);
            by0 = ((by0 % Ng) + Ng) % Ng;
            for (int k = 0; k < voxel_N; k++) {
                double z = -L + k * dx;
                int bz0 = (int)((z + L) / r->bin_size);
                bz0 = ((bz0 % Ng) + Ng) % Ng;
                int best = -1;
                double best_d2 = 1e30;
                for (int di = -1; di <= 1; di++) {
                    int gi = ((bx0 + di) % Ng + Ng) % Ng;
                    for (int dj = -1; dj <= 1; dj++) {
                        int gj = ((by0 + dj) % Ng + Ng) % Ng;
                        for (int dk = -1; dk <= 1; dk++) {
                            int gk = ((bz0 + dk) % Ng + Ng) % Ng;
                            long b = (long)gi * Ng * Ng + gj * Ng + gk;
                            int n = r->count[b], o = r->offset[b];
                            for (int q = 0; q < n; q++) {
                                int c = r->cells[o + q];
                                double cdx = m->cell_x[c] - x;
                                double cdy = m->cell_y[c] - y;
                                double cdz = m->cell_z[c] - z;
                                if (cdx >  0.5*S) cdx -= S; if (cdx < -0.5*S) cdx += S;
                                if (cdy >  0.5*S) cdy -= S; if (cdy < -0.5*S) cdy += S;
                                if (cdz >  0.5*S) cdz -= S; if (cdz < -0.5*S) cdz += S;
                                double d2 = cdx*cdx + cdy*cdy + cdz*cdz;
                                if (d2 < best_d2) { best_d2 = d2; best = c; }
                            }
                        }
                    }
                }
                r->voxel_to_cell[(long)i * voxel_N * voxel_N + j * voxel_N + k] = best;
            }
        }
    }
    return r;
}

static void resample_index_free(ResampleIndex *r) {
    free(r->count); free(r->offset); free(r->cells); free(r->voxel_to_cell);
    free(r);
}

/* =====================================================================
 *  SFA voxel snapshot — V56: writes 8 columns
 * ===================================================================== */

static void sfa_snap(SFA *sfa, Sim *s, ResampleIndex *r, double t,
                     float **vox_buf) {
    int N = r->voxel_N;
    long N3 = (long)N * N * N;
    for (int k = 0; k < NCOMP; k++) {
        double *src = s->M[k];
        float *dst = vox_buf[k];
        #pragma omp parallel for schedule(static)
        for (long idx = 0; idx < N3; idx++) {
            int c = r->voxel_to_cell[idx];
            dst[idx] = (c >= 0) ? (float)src[c] : 0.0f;
        }
    }
    void *cols[NCOMP];
    for (int k = 0; k < NCOMP; k++) cols[k] = vox_buf[k];
    sfa_write_frame(sfa, t, cols);
}

/* =====================================================================
 *  Cell-native I/P frame writers — V56: NCOMP=8 columns
 * ===================================================================== */

typedef struct {
    uint32_t N_cells;
    int snap_idx;
    int accum_count;
    float omega;
    float *mean, *amp, *phase;
    float *sum_cos, *sum_sin, *sum_mean;
    uint32_t *cell_ids;
    float    *delta_buf;
} CellTempModel;

static CellTempModel *temp_model_alloc(uint32_t N_cells, double omega) {
    CellTempModel *m = calloc(1, sizeof(CellTempModel));
    m->N_cells = N_cells;
    m->omega = (float)omega;
    long n = (long)N_cells * NCOMP;
    m->mean      = calloc(n, sizeof(float));
    m->amp       = calloc(n, sizeof(float));
    m->phase     = calloc(n, sizeof(float));
    m->sum_cos   = calloc(n, sizeof(float));
    m->sum_sin   = calloc(n, sizeof(float));
    m->sum_mean  = calloc(n, sizeof(float));
    m->cell_ids  = malloc(sizeof(uint32_t) * N_cells);
    m->delta_buf = malloc(sizeof(float) * n);
    return m;
}

static void temp_model_free(CellTempModel *m) {
    if (!m) return;
    free(m->mean); free(m->amp); free(m->phase);
    free(m->sum_cos); free(m->sum_sin); free(m->sum_mean);
    free(m->cell_ids); free(m->delta_buf);
    free(m);
}

static void temp_model_accumulate(CellTempModel *m, double t, float * const *fbuf) {
    long n = (long)m->N_cells * NCOMP;
    float cw = (float)cos(m->omega * t);
    float sw = (float)sin(m->omega * t);
    #pragma omp parallel for schedule(static)
    for (long i = 0; i < n; i++) {
        uint32_t cell = (uint32_t)(i / NCOMP);
        int col = (int)(i % NCOMP);
        float v = fbuf[col][cell];
        m->sum_mean[i] += v;
        m->sum_cos[i]  += v * cw;
        m->sum_sin[i]  += v * sw;
    }
    m->accum_count++;
}

static void temp_model_refit(CellTempModel *m, double t) {
    (void)t;
    if (m->accum_count < 3) return;
    long n = (long)m->N_cells * NCOMP;
    float inv_n = 1.0f / (float)m->accum_count;
    #pragma omp parallel for schedule(static)
    for (long i = 0; i < n; i++) {
        m->mean[i] = m->sum_mean[i] * inv_n;
        float sc = m->sum_cos[i] * inv_n;
        float ss = m->sum_sin[i] * inv_n;
        m->amp[i]   = 2.0f * sqrtf(sc*sc + ss*ss);
        m->phase[i] = atan2f(-ss, sc);
        m->sum_mean[i] = 0; m->sum_cos[i] = 0; m->sum_sin[i] = 0;
    }
    m->accum_count = 0;
}

static void temp_model_bootstrap(CellTempModel *m, float * const *fbuf) {
    long n = (long)m->N_cells * NCOMP;
    #pragma omp parallel for schedule(static)
    for (long i = 0; i < n; i++) {
        uint32_t cell = (uint32_t)(i / NCOMP);
        int col = (int)(i % NCOMP);
        m->mean[i]  = fbuf[col][cell];
        m->amp[i]   = 0;
        m->phase[i] = 0;
    }
}

static void cell_doubles_to_floats(Sim *s, float **fbuf) {
    uint32_t N = s->m->N_cells;
    for (int k = 0; k < NCOMP; k++) {
        double *src = s->M[k];
        float *dst = fbuf[k];
        #pragma omp parallel for schedule(static)
        for (uint32_t i = 0; i < N; i++) dst[i] = (float)src[i];
    }
}

static void sfa_snap_cell(SFA *sfa, Sim *s, double t, float **fbuf) {
    uint32_t N = s->m->N_cells;
    cell_doubles_to_floats(s, fbuf);
    void *cols[NCOMP];
    for (int k = 0; k < NCOMP; k++) cols[k] = fbuf[k];
    sfa_write_cell_frame(sfa, t, N, NCOMP, SFA_F32, cols);
}

static void sfa_snap_cell_iframe(SFA *sfa, Sim *s, CellTempModel *m,
                                  double t, float **fbuf) {
    uint32_t N = s->m->N_cells;
    cell_doubles_to_floats(s, fbuf);
    void *cols[NCOMP];
    for (int k = 0; k < NCOMP; k++) cols[k] = fbuf[k];
    sfa_write_cell_iframe_temporal(sfa, t, N, NCOMP, SFA_F32, cols,
                                    m->omega, m->mean, m->amp, m->phase);
}

static void sfa_snap_cell_iframe_lite(SFA *sfa, Sim *s,
                                       double t, float **fbuf) {
    uint32_t N = s->m->N_cells;
    cell_doubles_to_floats(s, fbuf);
    void *cols[NCOMP];
    for (int k = 0; k < NCOMP; k++) cols[k] = fbuf[k];
    sfa_write_cell_frame(sfa, t, N, NCOMP, SFA_F32, cols);
}

static void sfa_snap_cell_model(SFA *sfa, CellTempModel *m, double t) {
    sfa_write_temporal_model_frame(sfa, t, m->N_cells, NCOMP,
                                    m->omega, m->mean, m->amp, m->phase);
}

static void sfa_snap_cell_pframe(SFA *sfa, Sim *s, CellTempModel *m,
                                  double t, float **fbuf, double tol) {
    uint32_t N = s->m->N_cells;
    cell_doubles_to_floats(s, fbuf);
    float ftol = (float)tol;
    uint32_t n_changed = 0;
    for (uint32_t c = 0; c < N; c++) {
        long base = (long)c * NCOMP;
        float worst = 0;
        float row[NCOMP];
        for (int k = 0; k < NCOMP; k++) {
            float pred = m->mean[base + k]
                       + m->amp[base + k]
                       * cosf(m->omega * (float)t + m->phase[base + k]);
            float d = fbuf[k][c] - pred;
            row[k] = d;
            float ad = d < 0 ? -d : d;
            if (ad > worst) worst = ad;
        }
        if (worst > ftol) {
            m->cell_ids[n_changed] = c;
            for (int k = 0; k < NCOMP; k++)
                m->delta_buf[(long)n_changed * NCOMP + k] = row[k];
            n_changed++;
        }
    }
    sfa_write_cell_pframe(sfa, t, N, NCOMP, n_changed,
                           m->cell_ids, m->delta_buf);
    if ((m->snap_idx % 10) == 0) {
        printf("  [P-frame] t=%.3f: %u/%u cells changed (%.2f%% — tol %.4f)\n",
               t, n_changed, N, 100.0 * n_changed / N, tol);
    }
}

/* =====================================================================
 *  Main
 * ===================================================================== */

int main(int argc, char **argv) {
    if (argc < 3) {
        fprintf(stderr, "Usage: %s <foam_mesh.bin> <config.cfg>\n", argv[0]);
        return 1;
    }

    int nthreads = 8;
    char *env = getenv("OMP_NUM_THREADS");
    if (env) nthreads = atoi(env);
    omp_set_num_threads(nthreads);

    Mesh *m = mesh_read(argv[1]);

    Sim s = {0};
    s.m = m;
    config_default(&s.c);
    config_load(&s.c, argv[2]);

    /* Auto-enable hard S³ projection for Skyrmion seeds without a σ-penalty.
     * SCAN_RESULTS.md / PROJECTION_RESULT.md established that the soft
     * σ-constraint alone cannot hold |q|=1 against the bulk-Higgs scalar
     * branch; the unwinding is intrinsic, not a coupling-strength issue. */
    if (!strcmp(s.c.init, "skyrme") && s.c.skyrme_c4 > 0.0 &&
        !s.c.skyrme_project && s.c.sigma_e2 == 0.0) {
        printf("[sim] auto-enable skyrme_project=1 (Skyrme seed without σ-penalty needs hard projection)\n");
        s.c.skyrme_project = 1;
    }

    s.c.dt = s.c.dt_factor * m->dx_min;

    /* Tighten dt for the Skyrme L_4 dispersion. The quartic-gradient term
     * has ω² ~ c4·k⁴·⟨|∇q|²⟩; in a Skyrmion core ⟨|∇q|²⟩ ~ O(1) and
     * k_max ~ π/dx, giving ω_L4 ~ √c4·(π/dx)². Stability requires
     * dt ≲ 2/ω_L4 = (2/π²)·dx²/√c4 ≈ 0.2·dx²/√c4. We use the more
     * conservative 0.15 prefactor with the user's dt_factor as an upper bound. */
    if (s.c.skyrme_c4 > 0.0) {
        double dt_L4 = 0.15 * m->dx_min * m->dx_min / sqrt(s.c.skyrme_c4);
        if (dt_L4 < s.c.dt) {
            printf("[sim] CFL: tightening dt %.6f → %.6f for L_4 dispersion (c4=%.3f, dx_min=%.4f)\n",
                   s.c.dt, dt_L4, s.c.skyrme_c4, m->dx_min);
            s.c.dt = dt_L4;
        }
    }

    printf("[sim] field: %s (NCOMP=%d)\n", FIELD_NAME, NCOMP);
    printf("[sim] m²=%.4f λ=%.3f v²=%.4f (v=%.4f) sigma_e2=%.3f skyrme_c4=%.4f m_π²=%.4f project=%d\n",
           s.c.mass2, s.c.lambda, s.c.v2, sqrt(s.c.v2),
           s.c.sigma_e2, s.c.skyrme_c4, s.c.m_pi2, s.c.skyrme_project);
    if (s.c.damping > 0.0)
        printf("[sim] damping γ=%.4f for t<%.2f (then undamped)  mode=%s\n",
               s.c.damping, s.c.damping_T > 0.0 ? s.c.damping_T : INFINITY,
               s.c.damping_sector_size > 0 ? "tangent-only" : "full");
    printf("[sim] dt=%.6f dt_factor=%.4f T=%.1f threads=%d init=%s\n",
           s.c.dt, s.c.dt_factor, s.c.T, nthreads, s.c.init);

    sim_alloc(&s);

    /* SFA setup */
    SFA *sfa = NULL;
    ResampleIndex *resamp = NULL;
    float *vox_buf[NCOMP] = {0};
    float *cell_fbuf[NCOMP] = {0};
    int writing_sfa = (s.c.cell_native || s.c.sfa_output);
    if (writing_sfa) {
        sfa = sfa_create(s.c.output,
                         s.c.cell_native ? m->N_cells : s.c.voxel_N,
                         s.c.cell_native ? 1          : s.c.voxel_N,
                         s.c.cell_native ? 1          : s.c.voxel_N,
                         m->L, m->L, m->L, s.c.snap_dt);
        if (!sfa) { fprintf(stderr, "FATAL: cannot create SFA '%s'\n", s.c.output); return 1; }
        /* SFA column metadata pulled from the field plug-in. */
        for (int k = 0; k < NCOMP; k++) {
            sfa_add_column(sfa, field_names[k], SFA_F32,
                           field_semantic[k], field_component[k]);
        }
        sfa_finalize_header(sfa);
    }

    CellTempModel *temp = NULL;
    if (s.c.cell_native) {
        for (int k = 0; k < NCOMP; k++)
            cell_fbuf[k] = malloc(sizeof(float) * m->N_cells);
        double *pos_xyz = malloc(sizeof(double) * 3 * m->N_cells);
        for (uint32_t i = 0; i < m->N_cells; i++) {
            pos_xyz[3*i+0] = m->cell_x[i];
            pos_xyz[3*i+1] = m->cell_y[i];
            pos_xyz[3*i+2] = m->cell_z[i];
        }
        sfa_write_mesh_frame(sfa, 0.0, m->N_cells, 0, m->L,
                             pos_xyz, m->cell_vol, NULL, NULL, NULL, 0);
        free(pos_xyz);
        if (s.c.cell_iframe_interval > 0) {
            temp = temp_model_alloc(m->N_cells, s.c.cell_omega);
            printf("[sim] SFA cell-native: %s — FMSH + I/P-frame mix\n", s.c.output);
            printf("        I-frame every %d snaps, ω=%.3f, delta tol=%.4f\n",
                   s.c.cell_iframe_interval, s.c.cell_omega, s.c.cell_delta_tol);
        } else {
            printf("[sim] SFA cell-native: %s (FMSH + FCEL only, %u cells, NCOMP=%d)\n",
                   s.c.output, m->N_cells, NCOMP);
        }
    } else if (s.c.sfa_output) {
        resamp = resample_index_build(m, s.c.voxel_N);
        long N3 = (long)s.c.voxel_N * s.c.voxel_N * s.c.voxel_N;
        for (int k = 0; k < NCOMP; k++) vox_buf[k] = malloc(sizeof(float) * N3);
        printf("[sim] SFA voxel output: %s, voxel_N=%d (%d columns)\n",
               s.c.output, s.c.voxel_N, NCOMP);
    }

    field_init(&s);
    compute_forces(&s);

    FILE *dfp = fopen(s.c.diag_file, "w");
    fprintf(dfp, "t\tE_kin\tE_grad\tE_mass\tE_quartic\tE_sigma\tE_skyrme\tE_pion\tE_total"
                 "\tbulk_min\tbulk_mean\tbulk_max\trotor_rms\tspatial_rms"
                 "\tB_winding\tq_min\tq_mean\tq_max"
                 "\tM0_max\tM1_max\tM2_max\tM3_max\tM4_max\tM5_max\tM6_max\tM7_max\n");

    int n_steps = (int)(s.c.T / s.c.dt);
    int diag_every = (int)(s.c.diag_dt / s.c.dt); if (diag_every < 1) diag_every = 1;
    int snap_every = (int)(s.c.snap_dt / s.c.dt); if (snap_every < 1) snap_every = 1;
    printf("[sim] n_steps=%d diag_every=%d snap_every=%d\n",
           n_steps, diag_every, snap_every);

    Diag d0, d;
    compute_diag(&s, &d0);
    printf("[sim] E0=%.4e bulk[%.4f,%.4f,%.4f] rotor=%.4f spatial=%.4f |q|²[%.4f,%.4f] B=%.4f\n",
           d0.E_total, d0.bulk_min, d0.bulk_mean, d0.bulk_max,
           d0.rotor_rms, d0.spatial_rms,
           d0.q_norm_min, d0.q_norm_max, d0.B_winding);
    fprintf(dfp,
         "%.4f\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f"
         "\t%.6f\t%.6f\t%.6f\t%.6f"
         "\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n",
         0.0, d0.E_kin, d0.E_grad, d0.E_mass, d0.E_quartic, d0.E_sigma, d0.E_skyrme, d0.E_pion, d0.E_total,
        d0.bulk_min, d0.bulk_mean, d0.bulk_max, d0.rotor_rms, d0.spatial_rms,
        d0.B_winding, d0.q_norm_min, d0.q_norm_mean, d0.q_norm_max,
        d0.M_max[0], d0.M_max[1], d0.M_max[2], d0.M_max[3],
        d0.M_max[4], d0.M_max[5], d0.M_max[6], d0.M_max[7]);
    fflush(dfp);

    /* Snapshot dispatcher (mirrors v55 cell-native path) */
    #define SNAP_AT(_t) do {                                                    \
        if (s.c.cell_native && temp) {                                          \
            cell_doubles_to_floats(&s, cell_fbuf);                              \
            int is_iframe = (temp->snap_idx % s.c.cell_iframe_interval == 0);   \
            int iframe_idx = temp->snap_idx / s.c.cell_iframe_interval;         \
            int separate_model = (s.c.cell_model_interval > 0);                 \
            int write_model = separate_model && is_iframe                        \
                              && (iframe_idx % s.c.cell_model_interval == 0);   \
            temp_model_accumulate(temp, (_t), cell_fbuf);                       \
            if (is_iframe) {                                                    \
                if (temp->snap_idx == 0) temp_model_bootstrap(temp, cell_fbuf); \
                else                     temp_model_refit(temp, (_t));          \
                if (separate_model) {                                            \
                    if (write_model) sfa_snap_cell_model(sfa, temp, (_t));      \
                    sfa_snap_cell_iframe_lite(sfa, &s, (_t), cell_fbuf);        \
                } else {                                                        \
                    sfa_snap_cell_iframe(sfa, &s, temp, (_t), cell_fbuf);       \
                }                                                               \
            } else {                                                            \
                sfa_snap_cell_pframe(sfa, &s, temp, (_t), cell_fbuf,            \
                                      s.c.cell_delta_tol);                      \
            }                                                                   \
            temp->snap_idx++;                                                   \
        } else if (s.c.cell_native) {                                           \
            sfa_snap_cell(sfa, &s, (_t), cell_fbuf);                            \
        } else if (sfa) {                                                       \
            sfa_snap(sfa, &s, resamp, (_t), vox_buf);                           \
        }                                                                       \
    } while (0)

    SNAP_AT(0.0);
    double wall0 = omp_get_wtime();

    for (int step = 1; step <= n_steps; step++) {
        double t = step * s.c.dt;
        verlet_step(&s, t);
        if (step % diag_every == 0) {
            compute_diag(&s, &d);
            double drift = 100.0 * (d.E_total - d0.E_total) / (fabs(d0.E_total) + 1e-30);
            fprintf(dfp,
                "%.4f\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f"
                "\t%.6f\t%.6f\t%.6f\t%.6f"
                "\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n",
                t, d.E_kin, d.E_grad, d.E_mass, d.E_quartic, d.E_sigma, d.E_skyrme, d.E_pion, d.E_total,
                d.bulk_min, d.bulk_mean, d.bulk_max, d.rotor_rms, d.spatial_rms,
                d.B_winding, d.q_norm_min, d.q_norm_mean, d.q_norm_max,
                d.M_max[0], d.M_max[1], d.M_max[2], d.M_max[3],
                d.M_max[4], d.M_max[5], d.M_max[6], d.M_max[7]);
            fflush(dfp);
            if (step % (diag_every * 10) == 0) {
                double wall = omp_get_wtime() - wall0;
                printf("t=%7.2f E=%.4e (drift %+.3f%%) bulk=[%.3f,%.3f] rotor=%.3f |q|²=[%.3f,%.3f] B=%.3f [%.0f%% %.0fs %.1fms/step]\n",
                       t, d.E_total, drift, d.bulk_min, d.bulk_max,
                       d.rotor_rms, d.q_norm_min, d.q_norm_max, d.B_winding,
                       100.0 * step / n_steps, wall, 1000.0 * wall / step);
                fflush(stdout);
            }
        }
        if (step % snap_every == 0) SNAP_AT(t);
    }

    SNAP_AT(s.c.T);
    fclose(dfp);
    #undef SNAP_AT

    if (temp) temp_model_free(temp);
    if (sfa) sfa_close(sfa);
    if (resamp) resample_index_free(resamp);
    for (int k = 0; k < NCOMP; k++) {
        if (vox_buf[k])   free(vox_buf[k]);
        if (cell_fbuf[k]) free(cell_fbuf[k]);
    }

    double wall = omp_get_wtime() - wall0;
    printf("\n[sim] V56 stage A complete: %.0fs (%.1f min)\n", wall, wall/60);

    sim_free(&s);
    mesh_free(m);
    return 0;
}
