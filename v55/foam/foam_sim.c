/*  foam_sim.c — 6-field Cosserat on a Voronoi cell complex
 *
 *  Reads foam_mesh.bin (produced by gen_foam_mesh) and runs the same
 *  physics as sfa/sim/scp_sim.c but with finite-volume operators on the
 *  unstructured Voronoi mesh.
 *
 *  Equations (per cell c):
 *      ∂²φ_a/∂t² = ∇²φ_a - m²φ_a - V'(P)∂P/∂φ_a + η × curl(θ)_a
 *      ∂²θ_a/∂t² = ∇²θ_a - m_θ²θ_a              + η × curl(φ)_a
 *      V(P) = (μ/2) P²/(1+κP²),  P = φ₀φ₁φ₂
 *
 *  Finite-volume operators on the Voronoi-Delaunay dual mesh:
 *      ∇²f(c) = (1/V_c) Σ_n A_cn (f(n) − f(c)) / |δ_cn|
 *      ∂_b f(c) = (1/V_c) Σ_n A_cn n̂_cn[b] (f(n) − f(c))
 *      curl(F)_a(c) = ε_abc ∂_b F_c(c)
 *
 *  Build:
 *      gcc -O3 -march=native -fopenmp -o foam_sim foam_sim.c -lm
 *
 *  Usage:
 *      ./foam_sim foam_mesh.bin foam_run.cfg
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

#define NFIELDS 3
#define PI 3.14159265358979323846

/* =====================================================================
 *  Mesh structures
 * ===================================================================== */

/* Single face record — Array-of-Structs packing. One cache-line-friendly
 * read fetches all per-face fields together (vs the older SoA layout which
 * cost 10 separate cache misses per face access). */
typedef struct {
    uint32_t a, b;          /* cell indices                   8 bytes */
    double area;            /* shared face area               8 bytes */
    double dx, dy, dz;      /* δ = pos(b) − pos(a) periodic  24 bytes */
    double dist;            /* |δ|                            8 bytes */
    double nx, ny, nz;      /* unit normal δ/|δ|             24 bytes */
} Face;                     /* 72 bytes total */

typedef struct {
    /* Header */
    double L;
    uint32_t N_cells;
    uint32_t N_faces;

    /* Cell data: pos[3], volume */
    double *cell_x, *cell_y, *cell_z;
    double *cell_vol;

    /* Face data (AoS) */
    Face *faces;

    /* CSR cell→face index (each cell appears on |faces(c)| edges) */
    uint32_t *cell_face_off;   /* size N_cells+1 */
    uint32_t *cell_face_idx;   /* size 2*N_faces */

    double dx_min;             /* min face distance, for CFL */
} Mesh;

/* =====================================================================
 *  Morton-code cell reordering
 *
 *  Sorts cells so that spatially-nearby cells have nearby memory IDs,
 *  greatly improving cache behaviour for neighbour lookups in the
 *  Laplacian/curl computations. Also re-sorts faces by (min cell, max
 *  cell) and rebuilds the CSR. Runs once at mesh load.
 * ===================================================================== */

static inline uint64_t spread_bits_3d(uint32_t v) {
    uint64_t r = v & 0x1FFFFFULL;                 /* keep 21 low bits */
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
    double scale = (double)(1U << 20) / (2.0 * L);   /* 20 bits per axis */

    /* (Morton, oldId) pairs, sort by Morton */
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

    /* Build inv[old] = new and reorder cell_x,y,z,vol */
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

    /* Update face cell IDs through inv[] */
    #pragma omp parallel for schedule(static)
    for (uint32_t f = 0; f < m->N_faces; f++) {
        m->faces[f].a = inv[m->faces[f].a];
        m->faces[f].b = inv[m->faces[f].b];
    }
    free(inv);

    /* Sort faces by (min(a,b), max(a,b)) so face IDs cluster by cell */
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

    /* Rebuild CSR cell_face_off / cell_face_idx from new face order */
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

    printf("[mesh] Morton-reordered %u cells + %u faces (cache-friendly layout)\n",
           N, m->N_faces);
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
    if (version != 1) { fprintf(stderr, "FATAL: unsupported version %u\n", version); exit(1); }

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

    /* Cell records 32 bytes each: x,y,z,vol */
    for (uint32_t i = 0; i < m->N_cells; i++) {
        double rec[4];
        fread(rec, 8, 4, fp);
        m->cell_x[i] = rec[0];
        m->cell_y[i] = rec[1];
        m->cell_z[i] = rec[2];
        m->cell_vol[i] = rec[3];
    }

    m->faces = malloc(sizeof(Face) * m->N_faces);

    /* Disk format: a(4) b(4) area(8) dx(8) dy(8) dz(8) = 40 bytes per face.
     * In-memory: pack into Face struct, precompute dist + unit normal. */
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

    /* CSR offsets and face indices */
    m->cell_face_off = malloc(sizeof(uint32_t) * (m->N_cells + 1));
    fread(m->cell_face_off, 4, m->N_cells + 1, fp);
    uint32_t total_idx = m->cell_face_off[m->N_cells];
    m->cell_face_idx = malloc(sizeof(uint32_t) * total_idx);
    fread(m->cell_face_idx, 4, total_idx, fp);

    fclose(fp);

    /* Stats */
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
 *  Sim state
 * ===================================================================== */

typedef struct {
    /* Physics */
    double m2, mtheta2, eta, mu, kappa;

    /* Seed */
    double A, A_bg, R_tube, ellip;
    double delta[3];

    /* Time stepping */
    double T;
    double dt_factor;
    double dt;

    /* I/O */
    char init[64];
    char output[256];
    char diag_file[256];
    double snap_dt, diag_dt;
    int voxel_N;            /* voxel grid resolution for SFA output */
    int sfa_output;         /* 1 = direct SFA voxel writes, 0 = .fsnp dumps */
    int cell_native;        /* 1 = write FMSH+FCEL (cell-native, no resample),
                               overrides sfa_output */
    int cell_iframe_interval;   /* >0 enables temporal model + sparse P-frames;
                                   I-frame every N snaps. 0 = always I-frame. */
    int cell_model_interval;    /* how often to write a standalone FMTL chunk;
                                   0 (default) = embed model in every I-frame
                                   (legacy v2 behaviour); >0 = write FMTL every
                                   N I-frames and write I-frames model-less. */
    double cell_delta_tol;      /* threshold for FCEP entry inclusion */
    double cell_omega;          /* angular frequency for breathing model
                                   (defaults to 2π/2.2 ≈ 2.857) */
} Config;

typedef struct {
    Mesh *m;
    Config c;
    /* Field arrays: phi[a][cell] etc. */
    double *phi[NFIELDS], *phi_vel[NFIELDS], *phi_acc[NFIELDS];
    double *theta[NFIELDS], *theta_vel[NFIELDS], *theta_acc[NFIELDS];
    /* Per-cell scratch: gradient components, curl components */
    double *grad_phi[NFIELDS][3];   /* ∂_b φ_a */
    double *grad_theta[NFIELDS][3];
} Sim;

static double parse_double(const char *v) { return strtod(v, NULL); }

static void config_default(Config *c) {
    c->m2 = 2.25;
    c->mtheta2 = 0.0;
    c->eta = 0.5;
    c->mu = -41.345;
    c->kappa = 50.0;
    c->A = 0.8;
    c->A_bg = 0.1;
    c->R_tube = 3.0;
    c->ellip = 0.3325;
    c->delta[0] = 0.0;
    c->delta[1] = 3.0005;
    c->delta[2] = 4.4325;
    c->T = 50.0;
    c->dt_factor = 0.05;
    c->snap_dt = 10.0;
    c->diag_dt = 2.0;
    c->voxel_N = 128;
    c->sfa_output = 1;     /* default: write SFA voxel frames */
    c->cell_native = 0;    /* default off; set 1 for cell-native FMSH+FCEL */
    c->cell_iframe_interval = 0;          /* default: no P-frame compression */
    c->cell_model_interval = 0;           /* default: embed model (legacy) */
    c->cell_delta_tol = 0.01;
    /* Default to the carrier-wave temporal frequency at L=20 with default
     * params: ω_carrier = sqrt((π/L)² + m²) ≈ sqrt(2.275) ≈ 1.508. This
     * dominates the time-evolution of vacuum cells (which are most of
     * the box); soliton-interior cells with breathing add a subdominant
     * component that doesn't compress well with a single-frequency fit. */
    c->cell_omega = 1.508;
    strcpy(c->init, "braid");
    strcpy(c->output, "foam_run.sfa");
    strcpy(c->diag_file, "foam_run_diag.tsv");
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
        /* Trim key */
        char *e = key + strlen(key) - 1;
        while (e > key && (*e == ' ' || *e == '\t')) *e-- = '\0';
        char *vs = val;
        while (*vs == ' ' || *vs == '\t') vs++;
        e = vs + strlen(vs) - 1;
        while (e > vs && (*e == ' ' || *e == '\t' || *e == '\n')) *e-- = '\0';

        if (!strcmp(key, "m"))         c->m2 = parse_double(vs) * parse_double(vs);
        else if (!strcmp(key, "m_theta")) c->mtheta2 = parse_double(vs) * parse_double(vs);
        else if (!strcmp(key, "eta"))     c->eta = parse_double(vs);
        else if (!strcmp(key, "mu"))      c->mu = parse_double(vs);
        else if (!strcmp(key, "kappa"))   c->kappa = parse_double(vs);
        else if (!strcmp(key, "A"))       c->A = parse_double(vs);
        else if (!strcmp(key, "A_bg"))    c->A_bg = parse_double(vs);
        else if (!strcmp(key, "R_tube"))  c->R_tube = parse_double(vs);
        else if (!strcmp(key, "ellip"))   c->ellip = parse_double(vs);
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
    for (int a = 0; a < NFIELDS; a++) {
        s->phi[a]       = calloc(N, sizeof(double));
        s->phi_vel[a]   = calloc(N, sizeof(double));
        s->phi_acc[a]   = calloc(N, sizeof(double));
        s->theta[a]     = calloc(N, sizeof(double));
        s->theta_vel[a] = calloc(N, sizeof(double));
        s->theta_acc[a] = calloc(N, sizeof(double));
        for (int b = 0; b < 3; b++) {
            s->grad_phi[a][b]   = calloc(N, sizeof(double));
            s->grad_theta[a][b] = calloc(N, sizeof(double));
        }
    }
}

static void sim_free(Sim *s) {
    for (int a = 0; a < NFIELDS; a++) {
        free(s->phi[a]); free(s->phi_vel[a]); free(s->phi_acc[a]);
        free(s->theta[a]); free(s->theta_vel[a]); free(s->theta_acc[a]);
        for (int b = 0; b < 3; b++) {
            free(s->grad_phi[a][b]);
            free(s->grad_theta[a][b]);
        }
    }
}

/* =====================================================================
 *  Initialization
 * ===================================================================== */

static void init_braid(Sim *s) {
    Mesh *m = s->m;
    Config *c = &s->c;
    double L = m->L;
    double kw = PI / L;
    double omega = sqrt(kw * kw + c->m2);
    double omega_bg = omega;
    double sx = 1.0 + c->ellip;
    double sy = 1.0 - c->ellip;
    double inv2R2 = 1.0 / (2.0 * c->R_tube * c->R_tube);

    printf("[init] braid: A=%.3f R=%.2f ellip=%.4f A_bg=%.3f kw=%.4f omega=%.4f\n",
           c->A, c->R_tube, c->ellip, c->A_bg, kw, omega);

    #pragma omp parallel for schedule(static)
    for (uint32_t i = 0; i < m->N_cells; i++) {
        double x = m->cell_x[i];
        double y = m->cell_y[i];
        double z = m->cell_z[i];
        double r2e = x*x/(sx*sx) + y*y/(sy*sy);
        double env = exp(-r2e * inv2R2);
        for (int a = 0; a < NFIELDS; a++) {
            double ph = kw * z + c->delta[a];
            double ph_bg = kw * z + 2.0 * PI * a / 3.0;
            s->phi[a][i] = c->A * env * cos(ph) + c->A_bg * cos(ph_bg);
            s->phi_vel[a][i] = omega * c->A * env * sin(ph)
                             + omega_bg * c->A_bg * sin(ph_bg);
            s->theta[a][i] = 0.0;
            s->theta_vel[a][i] = 0.0;
        }
    }
}

/* =====================================================================
 *  Finite-volume operators
 *
 *  For each cell c, we accumulate over its faces. Each face stores
 *  (a, b, area, delta_ab, n̂_ab). For cell c at one end of face f:
 *      neighbor n = (a==c) ? b : a
 *      n̂_cn   = (a==c) ? +n̂_ab : −n̂_ab     (always points away from c)
 *      (f(n) − f(c)) appears in both Laplacian and gradient.
 * ===================================================================== */

/* Single-pass version: compute gradients + laplacians of ALL 6 fields
 * (φ₀, φ₁, φ₂, θ₀, θ₁, θ₂) in one cell visit. Each face is read once
 * per cell side but reused across all 6 fields, and the per-cell scratch
 * (accum arrays + neighbour fetches) sits in registers / L1.
 *
 * Replaces 12 separate calls to compute_grads_and_lap (which each iterated
 * every cell and every face). 6× fewer outer-loop traversals; face data
 * stays in cache across the 6 inner field iterations. */
static void compute_grads_and_lap_all(Sim *s,
        double *grad_x[6], double *grad_y[6], double *grad_z[6],
        double *lap[6]) {
    Mesh *m = s->m;
    uint32_t N = m->N_cells;
    double *fields[6] = {
        s->phi[0],   s->phi[1],   s->phi[2],
        s->theta[0], s->theta[1], s->theta[2]
    };

    #pragma omp parallel for schedule(static)
    for (uint32_t c = 0; c < N; c++) {
        double inv_V = 1.0 / m->cell_vol[c];
        double agx[6] = {0}, agy[6] = {0}, agz[6] = {0}, alp[6] = {0};
        double fc[6];
        for (int k = 0; k < 6; k++) fc[k] = fields[k][c];

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

            /* All 6 fields share this face's geometry — reuse it */
            for (int k = 0; k < 6; k++) {
                double diff = fields[k][nb] - fc[k];
                agx[k] += area * nx * diff;
                agy[k] += area * ny * diff;
                agz[k] += area * nz * diff;
                alp[k] += area * diff * inv_d;
            }
        }
        for (int k = 0; k < 6; k++) {
            grad_x[k][c] = 0.5 * agx[k] * inv_V;
            grad_y[k][c] = 0.5 * agy[k] * inv_V;
            grad_z[k][c] = 0.5 * agz[k] * inv_V;
            lap[k][c]    = alp[k] * inv_V;
        }
    }
}

static void compute_grads_and_lap(Sim *s, double *f, double *grad_x, double *grad_y, double *grad_z, double *lap) {
    Mesh *m = s->m;
    uint32_t N = m->N_cells;

    /* Zero arrays */
    memset(grad_x, 0, sizeof(double) * N);
    memset(grad_y, 0, sizeof(double) * N);
    memset(grad_z, 0, sizeof(double) * N);
    memset(lap,    0, sizeof(double) * N);

    /* Sum contributions from each face. We do NOT use #pragma omp here
     * because two cells share each face and we'd race on accumulator.
     * Instead, parallelise per-cell and re-iterate faces from both ends.
     */
    #pragma omp parallel for schedule(static)
    for (uint32_t c = 0; c < N; c++) {
        double inv_V = 1.0 / m->cell_vol[c];
        double gx = 0, gy = 0, gz = 0, lp = 0;
        uint32_t off0 = m->cell_face_off[c];
        uint32_t off1 = m->cell_face_off[c+1];
        double fc = f[c];
        for (uint32_t k = off0; k < off1; k++) {
            uint32_t fid = m->cell_face_idx[k];
            const Face *F = &m->faces[fid];
            uint32_t nb;
            double sign;       /* +1 if a==c, -1 if b==c */
            if (F->a == c) { nb = F->b; sign = +1.0; }
            else           { nb = F->a; sign = -1.0; }
            double area = F->area;
            double d    = F->dist;
            double nx   = sign * F->nx;
            double ny   = sign * F->ny;
            double nz   = sign * F->nz;
            double diff = f[nb] - fc;
            /* Gradient via divergence theorem with two-point flux */
            gx += area * nx * diff;
            gy += area * ny * diff;
            gz += area * nz * diff;
            /* Laplacian via two-point flux */
            lp += area * diff / d;
        }
        /* Gradient: divergence theorem with f_face ≈ (f(c)+f(n))/2 includes
         * a factor of 1/2. (Closed-cell sum makes the f(c) part vanish, so
         * the implementation accumulates A·n̂·(f(n)-f(c)), and we divide by 2.)
         * The Laplacian uses 2-point normal flux (f(n)-f(c))/d_f directly,
         * NO 1/2 factor.
         */
        grad_x[c] = 0.5 * gx * inv_V;
        grad_y[c] = 0.5 * gy * inv_V;
        grad_z[c] = 0.5 * gz * inv_V;
        lap[c]    = lp * inv_V;
    }
}

/* =====================================================================
 *  Force computation
 * ===================================================================== */

static void compute_forces(Sim *s) {
    Mesh *m = s->m;
    Config *c = &s->c;
    uint32_t N = m->N_cells;

    /* Step 1: gradients + laplacians of all 6 fields in one pass */
    static double *lap_phi[NFIELDS] = {NULL};
    static double *lap_theta[NFIELDS] = {NULL};
    if (lap_phi[0] == NULL) {
        for (int a = 0; a < NFIELDS; a++) {
            lap_phi[a]   = malloc(sizeof(double) * N);
            lap_theta[a] = malloc(sizeof(double) * N);
        }
    }
    double *gx_all[6] = {
        s->grad_phi[0][0], s->grad_phi[1][0], s->grad_phi[2][0],
        s->grad_theta[0][0], s->grad_theta[1][0], s->grad_theta[2][0]
    };
    double *gy_all[6] = {
        s->grad_phi[0][1], s->grad_phi[1][1], s->grad_phi[2][1],
        s->grad_theta[0][1], s->grad_theta[1][1], s->grad_theta[2][1]
    };
    double *gz_all[6] = {
        s->grad_phi[0][2], s->grad_phi[1][2], s->grad_phi[2][2],
        s->grad_theta[0][2], s->grad_theta[1][2], s->grad_theta[2][2]
    };
    double *lap_all[6] = {
        lap_phi[0], lap_phi[1], lap_phi[2],
        lap_theta[0], lap_theta[1], lap_theta[2]
    };
    compute_grads_and_lap_all(s, gx_all, gy_all, gz_all, lap_all);

    /* Step 2: assemble forces. curl(F)_0 = ∂_y F_z − ∂_z F_y, etc. */
    #pragma omp parallel for schedule(static)
    for (uint32_t i = 0; i < N; i++) {
        double p0 = s->phi[0][i], p1 = s->phi[1][i], p2 = s->phi[2][i];
        double P = p0 * p1 * p2;
        double den = 1.0 + c->kappa * P * P;
        double mPd2 = c->mu * P / (den * den);

        /* curl(theta) and curl(phi) at cell i */
        double cth0 = s->grad_theta[2][1][i] - s->grad_theta[1][2][i];
        double cth1 = s->grad_theta[0][2][i] - s->grad_theta[2][0][i];
        double cth2 = s->grad_theta[1][0][i] - s->grad_theta[0][1][i];
        double cph0 = s->grad_phi[2][1][i] - s->grad_phi[1][2][i];
        double cph1 = s->grad_phi[0][2][i] - s->grad_phi[2][0][i];
        double cph2 = s->grad_phi[1][0][i] - s->grad_phi[0][1][i];

        double dPda[3];
        dPda[0] = p1 * p2;
        dPda[1] = p0 * p2;
        dPda[2] = p0 * p1;

        double cth_arr[3] = {cth0, cth1, cth2};
        double cph_arr[3] = {cph0, cph1, cph2};

        for (int a = 0; a < NFIELDS; a++) {
            s->phi_acc[a][i] = lap_phi[a][i]
                             - c->m2 * s->phi[a][i]
                             - mPd2 * dPda[a]
                             + c->eta * cth_arr[a];
            s->theta_acc[a][i] = lap_theta[a][i]
                               - c->mtheta2 * s->theta[a][i]
                               + c->eta * cph_arr[a];
        }
    }
}

/* =====================================================================
 *  Velocity Verlet step
 * ===================================================================== */

static void verlet_step(Sim *s) {
    uint32_t N = s->m->N_cells;
    double hdt = 0.5 * s->c.dt, dt = s->c.dt;

    /* Half-kick */
    #pragma omp parallel for schedule(static)
    for (uint32_t i = 0; i < N; i++) {
        for (int a = 0; a < NFIELDS; a++) {
            s->phi_vel[a][i]   += hdt * s->phi_acc[a][i];
            s->theta_vel[a][i] += hdt * s->theta_acc[a][i];
        }
    }
    /* Drift */
    #pragma omp parallel for schedule(static)
    for (uint32_t i = 0; i < N; i++) {
        for (int a = 0; a < NFIELDS; a++) {
            s->phi[a][i]   += dt * s->phi_vel[a][i];
            s->theta[a][i] += dt * s->theta_vel[a][i];
        }
    }
    /* Recompute forces */
    compute_forces(s);
    /* Half-kick */
    #pragma omp parallel for schedule(static)
    for (uint32_t i = 0; i < N; i++) {
        for (int a = 0; a < NFIELDS; a++) {
            s->phi_vel[a][i]   += hdt * s->phi_acc[a][i];
            s->theta_vel[a][i] += hdt * s->theta_acc[a][i];
        }
    }
}

/* =====================================================================
 *  Diagnostics
 * ===================================================================== */

typedef struct {
    double E_phi_kin, E_theta_kin;
    double E_grad, E_mass, E_pot;
    double E_tgrad, E_tmass, E_coupling;
    double E_total;
    double phi_max, P_max, P_int, theta_rms;
} Diag;

static void compute_diag(Sim *s, Diag *d) {
    Mesh *m = s->m;
    Config *c = &s->c;
    uint32_t N = m->N_cells;

    double epk=0, etk=0, eg=0, em=0, ep=0, etg=0, etm=0, ec=0;
    double phi_max_sq = 0, P_max_abs = 0, P_int = 0, t_rms = 0;
    double V_total = 0;

    /* Cell-summed pieces: kinetic, mass, V(P), invariants. */
    #pragma omp parallel for reduction(+:epk,etk,em,ep,etm,P_int,t_rms,V_total) \
                             reduction(max:phi_max_sq,P_max_abs) schedule(static)
    for (uint32_t i = 0; i < N; i++) {
        double V_c = m->cell_vol[i];
        V_total += V_c;
        for (int a = 0; a < NFIELDS; a++) {
            double pv = s->phi_vel[a][i], tv = s->theta_vel[a][i];
            epk += 0.5 * pv * pv * V_c;
            etk += 0.5 * tv * tv * V_c;
            em  += 0.5 * c->m2 * s->phi[a][i] * s->phi[a][i] * V_c;
            etm += 0.5 * c->mtheta2 * s->theta[a][i] * s->theta[a][i] * V_c;
            double phi_sq = s->phi[a][i] * s->phi[a][i];
            if (phi_sq > phi_max_sq) phi_max_sq = phi_sq;
            t_rms += s->theta[a][i] * s->theta[a][i];
        }
        double P = s->phi[0][i] * s->phi[1][i] * s->phi[2][i];
        double Pa = fabs(P);
        if (Pa > P_max_abs) P_max_abs = Pa;
        P_int += Pa * V_c;
        ep += (c->mu / 2.0) * P * P / (1.0 + c->kappa * P * P) * V_c;
    }

    /* Face-summed gradient energy (the discretely-conserved form):
     *   H_grad = (1/2) Σ_f (A_f / d_f) (f_b − f_a)²  for each field.
     */
    #pragma omp parallel for reduction(+:eg,etg) schedule(static)
    for (uint32_t fid = 0; fid < m->N_faces; fid++) {
        const Face *F = &m->faces[fid];
        uint32_t a = F->a, b = F->b;
        double w = 0.5 * F->area / F->dist;
        for (int k = 0; k < NFIELDS; k++) {
            double dphi = s->phi[k][b]   - s->phi[k][a];
            double dtht = s->theta[k][b] - s->theta[k][a];
            eg  += w * dphi * dphi;
            etg += w * dtht * dtht;
        }
    }

    /* Face-symmetric coupling energy (matches the discrete H from which
     * the φ and θ EOMs are both derivable):
     *   H_coup = -(η/2) Σ_f A_f ε_abc n̂_b (φ_a(a)+φ_a(b)) (θ_c(b)−θ_c(a))
     */
    #pragma omp parallel for reduction(+:ec) schedule(static)
    for (uint32_t fid = 0; fid < m->N_faces; fid++) {
        const Face *F = &m->faces[fid];
        uint32_t ai = F->a, bi = F->b;
        double A = F->area;
        double nx = F->nx, ny = F->ny, nz = F->nz;
        double psum[3], tdiff[3];
        for (int k = 0; k < NFIELDS; k++) {
            psum[k]  = s->phi[k][ai]   + s->phi[k][bi];
            tdiff[k] = s->theta[k][bi] - s->theta[k][ai];
        }
        /* ε_abc n̂_b ψ_a(sum) θ_c(diff)
         *   = nx (ψ1 θ2 − ψ2 θ1) + ny (ψ2 θ0 − ψ0 θ2) + nz (ψ0 θ1 − ψ1 θ0)
         */
        double term =
            nx * (psum[1] * tdiff[2] - psum[2] * tdiff[1]) +
            ny * (psum[2] * tdiff[0] - psum[0] * tdiff[2]) +
            nz * (psum[0] * tdiff[1] - psum[1] * tdiff[0]);
        ec += -0.5 * c->eta * A * term;
    }

    d->E_phi_kin = epk;
    d->E_theta_kin = etk;
    d->E_grad = eg;
    d->E_mass = em;
    d->E_pot = ep;
    d->E_tgrad = etg;
    d->E_tmass = etm;
    d->E_coupling = ec;
    d->E_total = epk + etk + eg + em + ep + etg + etm + ec;
    d->phi_max = sqrt(phi_max_sq);
    d->P_max = P_max_abs;
    d->P_int = P_int;
    d->theta_rms = sqrt(t_rms / (3.0 * N));
}

/* =====================================================================
 *  Foam → voxel resampling (nearest-neighbor / Voronoi)
 *
 *  Build a uniform spatial bin index over the cell centers.  For each
 *  voxel position, search a 3×3×3 cube of bins around it and pick the
 *  closest cell.  That cell's field values become the voxel value.
 * ===================================================================== */

typedef struct {
    int Ng;
    double bin_size;
    int *count;        /* Ng^3 ints */
    int *offset;       /* Ng^3 ints */
    int *cells;        /* total ints, indices sorted by bin */
    /* Voxel→cell map (size voxel_N^3) */
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

    /* Pre-compute voxel→cell map (one nearest-neighbor lookup per voxel,
     * done once at startup since cells are stationary in Stage 1a). */
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
 *  SFA voxel snapshot
 * ===================================================================== */

static void sfa_snap(SFA *sfa, Sim *s, ResampleIndex *r, double t,
                     float **vox_buf) {
    int N = r->voxel_N;
    long N3 = (long)N * N * N;

    /* Resample 6 fields in parallel */
    for (int f = 0; f < 6; f++) {
        double *src = (f < 3) ? s->phi[f] : s->theta[f - 3];
        float *dst = vox_buf[f];
        #pragma omp parallel for schedule(static)
        for (long idx = 0; idx < N3; idx++) {
            int c = r->voxel_to_cell[idx];
            dst[idx] = (c >= 0) ? (float)src[c] : 0.0f;
        }
    }
    void *cols[6] = { vox_buf[0], vox_buf[1], vox_buf[2],
                      vox_buf[3], vox_buf[4], vox_buf[5] };
    sfa_write_frame(sfa, t, cols);
}

/* ---- Temporal model accumulator + delta encoder ----
 * Layout: per-(cell × column) state. We store [cell × NCOLS + col] with
 * NCOLS=6. Total per array: N_cells × 6 floats. */
typedef struct {
    uint32_t N_cells;
    int snap_idx;          /* count of snaps written so far */
    int accum_count;       /* frames included in current sum window */
    float omega;
    float *mean;           /* N_cells × 6 */
    float *amp;
    float *phase;
    float *sum_cos;
    float *sum_sin;
    float *sum_mean;
    /* Sparse-delta scratch (P-frame writer) */
    uint32_t *cell_ids;    /* up to N_cells */
    float    *delta_buf;   /* up to N_cells × 6 */
} CellTempModel;

static CellTempModel *temp_model_alloc(uint32_t N_cells, double omega) {
    CellTempModel *m = calloc(1, sizeof(CellTempModel));
    m->N_cells = N_cells;
    m->omega = (float)omega;
    long n = (long)N_cells * 6;
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

/* Update sliding accumulators with the current frame at time t. */
static void temp_model_accumulate(CellTempModel *m, double t, float * const *fbuf) {
    long n = (long)m->N_cells * 6;
    float cw = (float)cos(m->omega * t);
    float sw = (float)sin(m->omega * t);
    #pragma omp parallel for schedule(static)
    for (long i = 0; i < n; i++) {
        uint32_t cell = (uint32_t)(i / 6);
        int col = (int)(i % 6);
        float v = fbuf[col][cell];
        m->sum_mean[i] += v;
        m->sum_cos[i]  += v * cw;
        m->sum_sin[i]  += v * sw;
    }
    m->accum_count++;
}

/* Print summary stats of fit for a few sample cells. */
static void temp_model_inspect(CellTempModel *m, const char *tag) {
    double mean_amp = 0, max_amp = 0;
    long n = (long)m->N_cells * 6;
    for (long i = 0; i < n; i++) {
        if (m->amp[i] > max_amp) max_amp = m->amp[i];
        mean_amp += m->amp[i];
    }
    mean_amp /= n;
    /* Sample a few cells: cell 0 column 0, cell N/2 col 0 */
    long c0 = 0, cM = (long)m->N_cells / 2;
    printf("  [model %s] mean_amp=%.4f max_amp=%.4f | "
           "cell0: m=%.3f a=%.3f φ=%.2f | cellN/2: m=%.3f a=%.3f φ=%.2f\n",
           tag, mean_amp, max_amp,
           m->mean[c0*6], m->amp[c0*6], m->phase[c0*6],
           m->mean[cM*6], m->amp[cM*6], m->phase[cM*6]);
}

/* Refit (mean, amp, phase) from the accumulated sums.
 * Fourier projection: for value(t) = mean + A cos(ω·t + φ),
 *   <value × cos(ω·t)> = A/2 cos(φ),  <value × sin(ω·t)> = -A/2 sin(φ),
 * provided the accumulation window covers near-integer multiples of
 * the period (otherwise mean leaks into sc/ss; we ignore that small
 * correction here). */
static void temp_model_refit(CellTempModel *m, double t) {
    (void)t;
    if (m->accum_count < 3) return;
    long n = (long)m->N_cells * 6;
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

/* Bootstrap: on the very first frame, set mean=actual, amp=phase=0.
 * Reader's predicted = mean, so deltas start at zero. */
static void temp_model_bootstrap(CellTempModel *m, float * const *fbuf) {
    long n = (long)m->N_cells * 6;
    #pragma omp parallel for schedule(static)
    for (long i = 0; i < n; i++) {
        uint32_t cell = (uint32_t)(i / 6);
        int col = (int)(i % 6);
        m->mean[i]  = fbuf[col][cell];
        m->amp[i]   = 0;
        m->phase[i] = 0;
    }
}

/* Convert double phi/theta → float buffer (column-major, per cell) */
static void cell_doubles_to_floats(Sim *s, float **fbuf) {
    uint32_t N = s->m->N_cells;
    for (int f = 0; f < 6; f++) {
        double *src = (f < 3) ? s->phi[f] : s->theta[f - 3];
        float *dst = fbuf[f];
        #pragma omp parallel for schedule(static)
        for (uint32_t i = 0; i < N; i++) dst[i] = (float)src[i];
    }
}

/* Write a single FCEL frame: per-cell field values, f32 column-major.
 * Used in non-temporal cell-native mode (no model, no P-frames). */
static void sfa_snap_cell(SFA *sfa, Sim *s, double t, float **fbuf) {
    uint32_t N = s->m->N_cells;
    cell_doubles_to_floats(s, fbuf);
    void *cols[6] = { fbuf[0], fbuf[1], fbuf[2], fbuf[3], fbuf[4], fbuf[5] };
    sfa_write_cell_frame(sfa, t, N, 6, SFA_F32, cols);
}

/* I-frame writer with embedded temporal model (legacy v2 path). */
static void sfa_snap_cell_iframe(SFA *sfa, Sim *s, CellTempModel *m,
                                  double t, float **fbuf) {
    uint32_t N = s->m->N_cells;
    cell_doubles_to_floats(s, fbuf);
    void *cols[6] = { fbuf[0], fbuf[1], fbuf[2], fbuf[3], fbuf[4], fbuf[5] };
    sfa_write_cell_iframe_temporal(sfa, t, N, 6, SFA_F32, cols,
                                    m->omega, m->mean, m->amp, m->phase);
}

/* Lite I-frame: data only, no embedded model. Use when the model is
 * stored separately as an FMTL chunk. Reader takes the model from the
 * most recent FMTL frame in the stream. */
static void sfa_snap_cell_iframe_lite(SFA *sfa, Sim *s,
                                       double t, float **fbuf) {
    uint32_t N = s->m->N_cells;
    cell_doubles_to_floats(s, fbuf);
    void *cols[6] = { fbuf[0], fbuf[1], fbuf[2], fbuf[3], fbuf[4], fbuf[5] };
    sfa_write_cell_frame(sfa, t, N, 6, SFA_F32, cols);
}

/* Standalone FMTL writer: dumps the current temporal model. */
static void sfa_snap_cell_model(SFA *sfa, CellTempModel *m, double t) {
    sfa_write_temporal_model_frame(sfa, t, m->N_cells, 6,
                                    m->omega, m->mean, m->amp, m->phase);
}

/* P-frame writer: compute residuals against model, threshold, write FCEP. */
static void sfa_snap_cell_pframe(SFA *sfa, Sim *s, CellTempModel *m,
                                  double t, float **fbuf, double tol) {
    uint32_t N = s->m->N_cells;
    cell_doubles_to_floats(s, fbuf);
    float cw = (float)cos(m->omega * t);
    float sw = (float)sin(m->omega * t);
    /* Two-pass: collect changed-cell list serially after computing residuals.
     * Per cell, predict + residual for each column; if any column exceeds tol,
     * keep this cell. */
    float ftol = (float)tol;
    uint32_t n_changed = 0;
    /* Allocate per-cell temp residuals (6 floats per cell), reuse delta_buf. */
    /* Simplest correct: serial pass; for very large N we could parallelise
     * with per-thread chunks then merge. Serial keeps the cell ordering. */
    for (uint32_t c = 0; c < N; c++) {
        long base = (long)c * 6;
        float worst = 0;
        float row[6];
        for (int k = 0; k < 6; k++) {
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
            for (int k = 0; k < 6; k++)
                m->delta_buf[(long)n_changed * 6 + k] = row[k];
            n_changed++;
        }
    }
    sfa_write_cell_pframe(sfa, t, N, 6, n_changed,
                           m->cell_ids, m->delta_buf);
    /* Diagnostic: print sparsity at long intervals */
    if ((m->snap_idx % 10) == 0) {
        printf("  [P-frame] t=%.3f: %u/%u cells changed (%.2f%% — tol %.4f)\n",
               t, n_changed, N, 100.0 * n_changed / N, tol);
    }
    /* Suppress unused-var warning when cw/sw aren't used (cosf inline above) */
    (void)cw; (void)sw;
}

/* Legacy .fsnp writer (per-cell binary, kept for back-compat) */
static void snapshot_write(Sim *s, double t, const char *base) {
    char path[512];
    snprintf(path, sizeof(path), "%s_t%05d.fsnp", base, (int)(t + 0.5));
    FILE *fp = fopen(path, "wb");
    if (!fp) { fprintf(stderr, "WARN: cannot open snapshot '%s'\n", path); return; }
    fwrite("FSNP", 1, 4, fp);
    uint32_t version = 1; fwrite(&version, 4, 1, fp);
    fwrite(&s->m->N_cells, 4, 1, fp);
    fwrite(&t, 8, 1, fp);
    uint64_t flags = 0; fwrite(&flags, 8, 1, fp);
    for (uint32_t i = 0; i < s->m->N_cells; i++) {
        double rec[6] = { s->phi[0][i], s->phi[1][i], s->phi[2][i],
                          s->theta[0][i], s->theta[1][i], s->theta[2][i] };
        fwrite(rec, 8, 6, fp);
    }
    fclose(fp);
    printf("  [SNAP] %s (%lu MB)\n",
           path, (unsigned long)(s->m->N_cells * 48ULL / (1024*1024)));
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

    s.c.dt = s.c.dt_factor * m->dx_min;
    printf("[sim] m²=%.4f m_θ²=%.4f η=%.3f μ=%.4f κ=%.2f\n",
           s.c.m2, s.c.mtheta2, s.c.eta, s.c.mu, s.c.kappa);
    printf("[sim] dt=%.6f dt_factor=%.4f T=%.1f threads=%d\n",
           s.c.dt, s.c.dt_factor, s.c.T, nthreads);

    sim_alloc(&s);

    /* Output paths:
     *   cell_native = 1: write FMSH (mesh) + FCEL (cell data) frames directly.
     *                    No voxel resampling. Output is a foam-native SFA.
     *   sfa_output  = 1: voxel SFA via foam→voxel resampling at every snap.
     *   else (legacy):   .fsnp per-cell binary dumps. */
    SFA *sfa = NULL;
    ResampleIndex *resamp = NULL;
    float *vox_buf[6] = {0};
    float *cell_fbuf[6] = {0};
    int writing_sfa = (s.c.cell_native || s.c.sfa_output);
    if (writing_sfa) {
        sfa = sfa_create(s.c.output,
                         s.c.cell_native ? m->N_cells : s.c.voxel_N,
                         s.c.cell_native ? 1          : s.c.voxel_N,
                         s.c.cell_native ? 1          : s.c.voxel_N,
                         m->L, m->L, m->L, s.c.snap_dt);
        if (!sfa) { fprintf(stderr, "FATAL: cannot create SFA '%s'\n", s.c.output); return 1; }
        sfa_add_column(sfa, "phi_x",   SFA_F32, SFA_POSITION, 0);
        sfa_add_column(sfa, "phi_y",   SFA_F32, SFA_POSITION, 1);
        sfa_add_column(sfa, "phi_z",   SFA_F32, SFA_POSITION, 2);
        sfa_add_column(sfa, "theta_x", SFA_F32, SFA_ANGLE,    0);
        sfa_add_column(sfa, "theta_y", SFA_F32, SFA_ANGLE,    1);
        sfa_add_column(sfa, "theta_z", SFA_F32, SFA_ANGLE,    2);
        sfa_finalize_header(sfa);
    }
    CellTempModel *temp = NULL;
    if (s.c.cell_native) {
        /* Cell-native mode: write a single FMSH frame at t=0 (positions
         * + volumes only — no faces/CSR, since visualisation only needs
         * the spatial index). All subsequent frames are FCEL or FCEP. */
        for (int f = 0; f < 6; f++) cell_fbuf[f] = malloc(sizeof(float) * m->N_cells);
        double *pos_xyz = malloc(sizeof(double) * 3 * m->N_cells);
        for (uint32_t i = 0; i < m->N_cells; i++) {
            pos_xyz[3*i+0] = m->cell_x[i];
            pos_xyz[3*i+1] = m->cell_y[i];
            pos_xyz[3*i+2] = m->cell_z[i];
        }
        sfa_write_mesh_frame(sfa, 0.0, m->N_cells, 0, m->L,
                             pos_xyz, m->cell_vol, NULL, NULL, NULL,
                             /* flags: f32 positions, no faces */ 0);
        free(pos_xyz);
        if (s.c.cell_iframe_interval > 0) {
            temp = temp_model_alloc(m->N_cells, s.c.cell_omega);
            printf("[sim] SFA cell-native: %s — FMSH + I/P-frame mix\n", s.c.output);
            printf("        I-frame every %d snaps, ω=%.3f, delta tol=%.4f\n",
                   s.c.cell_iframe_interval, s.c.cell_omega, s.c.cell_delta_tol);
        } else {
            printf("[sim] SFA cell-native: %s (FMSH + FCEL frames, %u cells)\n",
                   s.c.output, m->N_cells);
        }
    } else if (s.c.sfa_output) {
        resamp = resample_index_build(m, s.c.voxel_N);
        long N3 = (long)s.c.voxel_N * s.c.voxel_N * s.c.voxel_N;
        for (int f = 0; f < 6; f++) vox_buf[f] = malloc(sizeof(float) * N3);
        printf("[sim] SFA voxel output: %s, voxel_N=%d, dx=%.3f\n",
               s.c.output, s.c.voxel_N, 2.0 * m->L / (s.c.voxel_N - 1));
    }

    if (!strcmp(s.c.init, "braid")) init_braid(&s);
    else { fprintf(stderr, "FATAL: unknown init '%s'\n", s.c.init); return 1; }

    compute_forces(&s);

    FILE *dfp = fopen(s.c.diag_file, "w");
    fprintf(dfp, "t\tE_phi_kin\tE_theta_kin\tE_grad\tE_mass\tE_pot\tE_tgrad\tE_tmass\tE_coupling\tE_total\tphi_max\tP_max\tP_int\ttheta_rms\n");

    int n_steps = (int)(s.c.T / s.c.dt);
    int diag_every = (int)(s.c.diag_dt / s.c.dt); if (diag_every < 1) diag_every = 1;
    int snap_every = (int)(s.c.snap_dt / s.c.dt); if (snap_every < 1) snap_every = 1;
    printf("[sim] n_steps=%d diag_every=%d snap_every=%d\n",
           n_steps, diag_every, snap_every);

    Diag d0, d;
    compute_diag(&s, &d0);
    printf("[sim] E0=%.4e E_pot=%.1f phi_max=%.4f P_max=%.4f P_int=%.1f θ_rms=%.4f\n",
           d0.E_total, d0.E_pot, d0.phi_max, d0.P_max, d0.P_int, d0.theta_rms);

    /* Helper macro: dispatch a snapshot to the right writer, including the
     * temporal-model accumulator and I/P-frame logic when active. */
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
                if (temp->snap_idx == 0) {                                      \
                    temp_model_bootstrap(temp, cell_fbuf);                      \
                    temp_model_inspect(temp, "bootstrap");                      \
                } else {                                                        \
                    temp_model_refit(temp, (_t));                               \
                    temp_model_inspect(temp, "refit");                          \
                }                                                               \
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
        } else {                                                                \
            snapshot_write(&s, (_t), s.c.output);                               \
        }                                                                       \
    } while (0)

    SNAP_AT(0.0);
    fprintf(dfp, "%.4f\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.4f\t%.4f\t%.6e\t%.6e\n",
            0.0, d0.E_phi_kin, d0.E_theta_kin, d0.E_grad, d0.E_mass, d0.E_pot,
            d0.E_tgrad, d0.E_tmass, d0.E_coupling, d0.E_total,
            d0.phi_max, d0.P_max, d0.P_int, d0.theta_rms);
    fflush(dfp);

    double wall0 = omp_get_wtime();

    for (int step = 1; step <= n_steps; step++) {
        verlet_step(&s);
        double t = step * s.c.dt;
        if (step % diag_every == 0) {
            compute_diag(&s, &d);
            double drift = 100.0 * (d.E_total - d0.E_total) / (fabs(d0.E_total) + 1e-30);
            fprintf(dfp, "%.4f\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.4f\t%.4f\t%.6e\t%.6e\n",
                    t, d.E_phi_kin, d.E_theta_kin, d.E_grad, d.E_mass, d.E_pot,
                    d.E_tgrad, d.E_tmass, d.E_coupling, d.E_total,
                    d.phi_max, d.P_max, d.P_int, d.theta_rms);
            fflush(dfp);
            if (step % (diag_every * 10) == 0) {
                double wall = omp_get_wtime() - wall0;
                printf("t=%7.2f E=%.4e (drift %+.3f%%) E_pot=%.1f phi=%.4f θ_rms=%.4f [%.0f%% %.0fs %.1fms/step]\n",
                       t, d.E_total, drift, d.E_pot, d.phi_max, d.theta_rms,
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
    for (int f = 0; f < 6; f++) {
        if (vox_buf[f])   free(vox_buf[f]);
        if (cell_fbuf[f]) free(cell_fbuf[f]);
    }

    double wall = omp_get_wtime() - wall0;
    printf("\n[sim] Complete: %.0fs (%.1f min)\n", wall, wall/60);

    sim_free(&s);
    mesh_free(m);
    return 0;
}
