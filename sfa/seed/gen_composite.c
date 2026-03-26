/*  gen_composite.c — Place pre-computed baryon templates into a composite seed
 *
 *  Reads template SFA files (e.g. proton_template.sfa) and stamps them into a
 *  larger output grid at specified positions. Handles multiple particles,
 *  background field initialization, and optional gradient backgrounds.
 *
 *  Build: gcc -O3 -fopenmp -o gen_composite gen_composite.c -lzstd -lm
 *
 *  Usage:
 *    # Single proton at origin
 *    ./gen_composite -N 192 -L 30 -template proton_template.sfa -place 0,0,0 -o proton.sfa
 *
 *    # Deuterium: proton + neutron separated by 40
 *    ./gen_composite -N 512 -L 100 -template proton_template.sfa -place -20,0,0 \
 *        -template neutron_template.sfa -place 20,0,0 -o deuterium.sfa
 *
 *    # Proton in gradient background (gravity test, bc_type=1)
 *    ./gen_composite -N 512 -L 100 -template proton_template.sfa -place 0,0,0 \
 *        -gradient_A_high 0.15 -gradient_A_low 0.05 -o gradient_test.sfa
 *
 *    # Custom precision and background
 *    ./gen_composite -N 256 -L 50 -A_bg 0.12 -precision f16 \
 *        -template p.sfa -place -10,0,0 -template p.sfa -place 10,0,0 -o he2.sfa
 */

#define SFA_IMPLEMENTATION
#include "../format/sfa.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define PI 3.14159265358979323846
#define MAX_PARTICLES 8
#define NCOLS 12   /* 3 phi + 3 theta + 3 phi_vel + 3 theta_vel */

/* ---------- Template data loaded from SFA ---------- */
typedef struct {
    char path[512];
    uint32_t Nx, Ny, Nz;
    double Lx, Ly, Lz;      /* half-domain extents */
    double dx, dy, dz;
    float *data;             /* NCOLS * Nx*Ny*Nz floats (column-major per SFA) */
    double cx, cy, cz;      /* placement center in output world coords */
} Template;

static void usage(const char *prog) {
    fprintf(stderr,
        "Usage: %s [options] -template <file.sfa> -place <cx>,<cy>,<cz> [...] -o <output.sfa>\n"
        "\n"
        "Options:\n"
        "  -N <grid_size>         Output grid size (default 192)\n"
        "  -L <domain_size>       Output half-domain (default 30)\n"
        "  -o <output.sfa>        Output file path\n"
        "  -template <file.sfa>   Template SFA file (repeatable, up to %d)\n"
        "  -place <cx>,<cy>,<cz>  Placement center for preceding template\n"
        "  -A_bg <value>          Background amplitude (default 0.1)\n"
        "  -gradient_A_high <v>   Background amplitude at x=+L (enables gradient)\n"
        "  -gradient_A_low <v>    Background amplitude at x=-L (enables gradient)\n"
        "  -m <mass>              Field mass, m^2 used in dispersion (default 1.5)\n"
        "  -precision f16|f32|f64 Output precision (default f32)\n"
        "  -delta <d0>,<d1>,<d2>  Phase offsets (default 0,3.0005,4.4325)\n"
        "\n"
        "Templates and placements are paired: each -template must be followed by -place.\n"
        "Multiple particles: repeat -template ... -place ... pairs.\n"
        "\n"
        "The template background is subtracted before stamping so that only the\n"
        "baryon perturbation is transferred (avoids double-counting the cosine).\n",
        prog, MAX_PARTICLES);
}

/* Load a template SFA, converting all columns to float.
 * SFA read_frame returns interleaved columns in native dtype. */
static int load_template(Template *t) {
    SFA *sfa = sfa_open(t->path);
    if (!sfa) {
        fprintf(stderr, "ERROR: cannot open template '%s'\n", t->path);
        return -1;
    }
    if (sfa->n_columns != NCOLS) {
        fprintf(stderr, "ERROR: template '%s' has %u columns, expected %d\n",
                t->path, sfa->n_columns, NCOLS);
        sfa_close(sfa);
        return -1;
    }

    t->Nx = sfa->Nx; t->Ny = sfa->Ny; t->Nz = sfa->Nz;
    t->Lx = sfa->Lx; t->Ly = sfa->Ly; t->Lz = sfa->Lz;
    t->dx = 2.0 * t->Lx / (t->Nx - 1);
    t->dy = 2.0 * t->Ly / (t->Ny - 1);
    t->dz = 2.0 * t->Lz / (t->Nz - 1);

    uint64_t N3 = (uint64_t)t->Nx * t->Ny * t->Nz;

    /* Allocate raw frame buffer */
    void *raw = malloc(sfa->frame_bytes);
    if (!raw) { fprintf(stderr, "ERROR: malloc failed for template frame\n"); sfa_close(sfa); return -1; }

    /* Read last frame (most converged) */
    uint32_t frame = (sfa->total_frames > 0) ? sfa->total_frames - 1 : 0;
    if (sfa_read_frame(sfa, frame, raw) < 0) {
        fprintf(stderr, "ERROR: failed to read frame %u from '%s'\n", frame, t->path);
        free(raw); sfa_close(sfa); return -1;
    }

    /* Convert all columns to float, stored as col[c][idx] */
    t->data = (float *)malloc(NCOLS * N3 * sizeof(float));
    if (!t->data) { fprintf(stderr, "ERROR: malloc failed\n"); free(raw); sfa_close(sfa); return -1; }

    uint64_t offset = 0;
    for (int c = 0; c < NCOLS; c++) {
        int es = sfa_dtype_size[sfa->columns[c].dtype];
        uint8_t dtype = sfa->columns[c].dtype;
        float *dst = t->data + (uint64_t)c * N3;

        for (uint64_t i = 0; i < N3; i++) {
            if (dtype == SFA_F32) {
                float v; memcpy(&v, (uint8_t*)raw + offset + i * es, 4);
                dst[i] = v;
            } else if (dtype == SFA_F64) {
                double v; memcpy(&v, (uint8_t*)raw + offset + i * es, 8);
                dst[i] = (float)v;
            } else if (dtype == SFA_F16) {
                uint16_t h; memcpy(&h, (uint8_t*)raw + offset + i * es, 2);
                /* f16 -> f32 */
                uint32_t sign = (h & 0x8000) << 16;
                int exp = (h >> 10) & 0x1F;
                uint32_t mant = h & 0x3FF;
                uint32_t f;
                if (exp == 0) {
                    if (mant == 0) f = sign;
                    else { /* denorm */
                        exp = 1;
                        while (!(mant & 0x400)) { mant <<= 1; exp--; }
                        mant &= 0x3FF;
                        f = sign | ((exp + 127 - 15) << 23) | (mant << 13);
                    }
                } else if (exp == 31) {
                    f = sign | 0x7F800000 | (mant << 13);
                } else {
                    f = sign | ((exp + 127 - 15) << 23) | (mant << 13);
                }
                float v; memcpy(&v, &f, 4);
                dst[i] = v;
            } else {
                fprintf(stderr, "ERROR: unsupported dtype %d in template\n", dtype);
                free(raw); free(t->data); sfa_close(sfa); return -1;
            }
        }
        offset += N3 * es;
    }

    free(raw);
    sfa_close(sfa);

    fprintf(stderr, "  Loaded template '%s': %ux%ux%u, L=(%.2f,%.2f,%.2f), dx=%.4f\n",
            t->path, t->Nx, t->Ny, t->Nz, t->Lx, t->Ly, t->Lz, t->dx);
    return 0;
}

/* Trilinear interpolation of template column c at fractional voxel (fi, fj, fk) */
static float interp_template(const Template *t, int c, double fi, double fj, double fk) {
    uint64_t N3 = (uint64_t)t->Nx * t->Ny * t->Nz;
    const float *col = t->data + (uint64_t)c * N3;
    int NN = t->Ny * t->Nz;

    int i0 = (int)floor(fi), j0 = (int)floor(fj), k0 = (int)floor(fk);
    double fx = fi - i0, fy = fj - j0, fz = fk - k0;

    /* Clamp to valid range */
    int i1 = i0 + 1, j1 = j0 + 1, k1 = k0 + 1;
    if (i0 < 0) i0 = 0; if (i1 >= (int)t->Nx) i1 = t->Nx - 1;
    if (j0 < 0) j0 = 0; if (j1 >= (int)t->Ny) j1 = t->Ny - 1;
    if (k0 < 0) k0 = 0; if (k1 >= (int)t->Nz) k1 = t->Nz - 1;
    if (i0 >= (int)t->Nx) i0 = t->Nx - 1;
    if (j0 >= (int)t->Ny) j0 = t->Ny - 1;
    if (k0 >= (int)t->Nz) k0 = t->Nz - 1;

    #define TV(ii,jj,kk) col[(long)(ii)*NN + (long)(jj)*(int)t->Nz + (kk)]
    float c000 = TV(i0,j0,k0), c100 = TV(i1,j0,k0);
    float c010 = TV(i0,j1,k0), c110 = TV(i1,j1,k0);
    float c001 = TV(i0,j0,k1), c101 = TV(i1,j0,k1);
    float c011 = TV(i0,j1,k1), c111 = TV(i1,j1,k1);
    #undef TV

    float c00 = c000*(1-fx) + c100*fx;
    float c01 = c001*(1-fx) + c101*fx;
    float c10 = c010*(1-fx) + c110*fx;
    float c11 = c011*(1-fx) + c111*fx;
    float c0  = c00*(1-fy) + c10*fy;
    float c1  = c01*(1-fy) + c11*fy;
    return c0*(1-fz) + c1*fz;
}

int main(int argc, char **argv) {
    /* Defaults */
    int N = 192;
    double L = 30.0;
    double A_bg = 0.1;
    double m2 = 2.25;
    double delta[3] = {0.0, 3.0005, 4.4325};
    double gradient_A_high = -1, gradient_A_low = -1;  /* <0 = not set */
    char outpath[512] = "composite.sfa";
    int precision = 1;  /* f32 */

    Template particles[MAX_PARTICLES];
    int n_particles = 0;
    int pending_template = 0;  /* 1 if last arg was -template without -place */

    if (argc < 2) { usage(argv[0]); return 1; }

    /* Parse arguments — templates and placements are paired sequentially */
    for (int i = 1; i < argc; i++) {
        const char *k = argv[i];
        if (!strcmp(k, "-h") || !strcmp(k, "--help")) { usage(argv[0]); return 0; }

        if (i + 1 >= argc && strcmp(k, "-h") && strcmp(k, "--help")) {
            /* Check for flags that need a value */
            if (k[0] == '-') { fprintf(stderr, "ERROR: missing value for %s\n", k); return 1; }
            continue;
        }
        const char *v = argv[i+1];

        if (!strcmp(k, "-N"))       { N = atoi(v); i++; }
        else if (!strcmp(k, "-L")) { L = atof(v); i++; }
        else if (!strcmp(k, "-o")) { strncpy(outpath, v, 511); i++; }
        else if (!strcmp(k, "-A_bg")) { A_bg = atof(v); i++; }
        else if (!strcmp(k, "-m"))  { double m = atof(v); m2 = m*m; i++; }
        else if (!strcmp(k, "-delta")) { sscanf(v, "%lf,%lf,%lf", &delta[0], &delta[1], &delta[2]); i++; }
        else if (!strcmp(k, "-gradient_A_high")) { gradient_A_high = atof(v); i++; }
        else if (!strcmp(k, "-gradient_A_low"))  { gradient_A_low = atof(v); i++; }
        else if (!strcmp(k, "-precision")) {
            if (!strcmp(v,"f16")) precision = 0;
            else if (!strcmp(v,"f32")) precision = 1;
            else if (!strcmp(v,"f64")) precision = 2;
            i++;
        }
        else if (!strcmp(k, "-template")) {
            if (n_particles >= MAX_PARTICLES) {
                fprintf(stderr, "ERROR: max %d particles\n", MAX_PARTICLES);
                return 1;
            }
            if (pending_template) {
                fprintf(stderr, "ERROR: -template without -place for '%s'\n",
                        particles[n_particles-1].path);
                return 1;
            }
            memset(&particles[n_particles], 0, sizeof(Template));
            strncpy(particles[n_particles].path, v, 511);
            pending_template = 1;
            n_particles++;
            i++;
        }
        else if (!strcmp(k, "-place")) {
            if (!pending_template) {
                fprintf(stderr, "ERROR: -place without preceding -template\n");
                return 1;
            }
            double cx, cy, cz;
            if (sscanf(v, "%lf,%lf,%lf", &cx, &cy, &cz) != 3) {
                fprintf(stderr, "ERROR: -place requires cx,cy,cz (got '%s')\n", v);
                return 1;
            }
            particles[n_particles-1].cx = cx;
            particles[n_particles-1].cy = cy;
            particles[n_particles-1].cz = cz;
            pending_template = 0;
            i++;
        }
    }

    if (pending_template) {
        fprintf(stderr, "ERROR: -template without -place for '%s'\n",
                particles[n_particles-1].path);
        return 1;
    }
    if (n_particles == 0) {
        fprintf(stderr, "ERROR: no particles specified (use -template ... -place ...)\n");
        return 1;
    }

    /* Validate gradient */
    int use_gradient = 0;
    if (gradient_A_high >= 0 && gradient_A_low >= 0) {
        use_gradient = 1;
        fprintf(stderr, "gen_composite: gradient mode A_low=%.4f A_high=%.4f\n",
                gradient_A_low, gradient_A_high);
    } else if (gradient_A_high >= 0 || gradient_A_low >= 0) {
        fprintf(stderr, "ERROR: both -gradient_A_high and -gradient_A_low must be set\n");
        return 1;
    }

    long N3 = (long)N * N * N;
    double dx = 2.0 * L / (N - 1);
    int NN = N * N;
    double k_bg = PI / L;
    double omega_bg = sqrt(k_bg * k_bg + m2);

    fprintf(stderr, "gen_composite: N=%d L=%.1f dx=%.4f A_bg=%.3f m=%.4f\n",
            N, L, dx, A_bg, sqrt(m2));
    fprintf(stderr, "gen_composite: %d particle(s)\n", n_particles);

    /* Load all templates */
    for (int p = 0; p < n_particles; p++) {
        fprintf(stderr, "  Particle %d: center=(%.2f,%.2f,%.2f)\n",
                p, particles[p].cx, particles[p].cy, particles[p].cz);
        if (load_template(&particles[p]) < 0) return 1;
    }

    /* Allocate output field arrays (double for accumulation) */
    double *fields[NCOLS];
    for (int c = 0; c < NCOLS; c++) {
        fields[c] = (double *)calloc(N3, sizeof(double));
        if (!fields[c]) { fprintf(stderr, "ERROR: malloc failed for output field\n"); return 1; }
    }

    /* Step 1: Initialize background */
    fprintf(stderr, "gen_composite: initializing background...\n");
    #pragma omp parallel for schedule(static)
    for (long idx = 0; idx < N3; idx++) {
        int i = idx / NN;
        int j = (idx / N) % N;
        int k = idx % N;
        double x = -L + i * dx;
        double y = -L + j * dx;
        double z = -L + k * dx;
        (void)y;  /* background is z-dependent only */

        /* Effective background amplitude (uniform or gradient in x) */
        double A_eff = A_bg;
        if (use_gradient) {
            double frac = (x + L) / (2.0 * L);  /* 0 at x=-L, 1 at x=+L */
            A_eff = gradient_A_low + (gradient_A_high - gradient_A_low) * frac;
        }

        for (int a = 0; a < 3; a++) {
            double phase = k_bg * z + 2.0 * PI * a / 3.0;
            fields[a][idx]     = A_eff * cos(phase);      /* phi_a */
            fields[6 + a][idx] = omega_bg * A_eff * sin(phase);  /* phi_vel_a */
        }
        /* theta and theta_vel remain zero from calloc */
    }

    /* Step 2: For each particle, subtract template background and add perturbation */
    for (int p = 0; p < n_particles; p++) {
        Template *t = &particles[p];
        uint64_t tN3 = (uint64_t)t->Nx * t->Ny * t->Nz;
        int tNN = t->Ny * t->Nz;

        /* Background wavenumber for the template's own domain */
        double k_tmpl = PI / t->Lz;
        double omega_tmpl = sqrt(k_tmpl * k_tmpl + m2);

        /* Determine if we need interpolation (dx mismatch > 1%) */
        int need_interp = (fabs(t->dx - dx) / dx > 0.01);
        if (need_interp) {
            fprintf(stderr, "  Particle %d: interpolation mode (template dx=%.4f, output dx=%.4f)\n",
                    p, t->dx, dx);
        }

        /* Template physical extent: [-Lx, +Lx] centered at (cx, cy, cz) */
        /* For each output voxel in the affected region, sample the template */

        /* Compute output voxel range that the template covers */
        double xmin = t->cx - t->Lx, xmax = t->cx + t->Lx;
        double ymin = t->cy - t->Ly, ymax = t->cy + t->Ly;
        double zmin = t->cz - t->Lz, zmax = t->cz + t->Lz;

        int oi_min = (int)floor((xmin + L) / dx);
        int oi_max = (int)ceil((xmax + L) / dx);
        int oj_min = (int)floor((ymin + L) / dx);
        int oj_max = (int)ceil((ymax + L) / dx);
        int ok_min = (int)floor((zmin + L) / dx);
        int ok_max = (int)ceil((zmax + L) / dx);

        /* Clamp to output grid */
        if (oi_min < 0) oi_min = 0; if (oi_max >= N) oi_max = N - 1;
        if (oj_min < 0) oj_min = 0; if (oj_max >= N) oj_max = N - 1;
        if (ok_min < 0) ok_min = 0; if (ok_max >= N) ok_max = N - 1;

        fprintf(stderr, "  Particle %d: output range [%d:%d, %d:%d, %d:%d]\n",
                p, oi_min, oi_max, oj_min, oj_max, ok_min, ok_max);

        #pragma omp parallel for collapse(2) schedule(static)
        for (int oi = oi_min; oi <= oi_max; oi++) {
            for (int oj = oj_min; oj <= oj_max; oj++) {
                for (int ok = ok_min; ok <= ok_max; ok++) {
                    /* Output world coords */
                    double ox = -L + oi * dx;
                    double oy = -L + oj * dx;
                    double oz = -L + ok * dx;

                    /* Map to template fractional voxel */
                    double tx = ox - t->cx;  /* relative to template center */
                    double ty = oy - t->cy;
                    double tz = oz - t->cz;

                    /* Template voxel indices (template grid: [-Lx..+Lx]) */
                    double fi = (tx + t->Lx) / t->dx;
                    double fj = (ty + t->Ly) / t->dy;
                    double fk = (tz + t->Lz) / t->dz;

                    /* Check bounds — skip if outside template */
                    if (fi < -0.5 || fi > t->Nx - 0.5) continue;
                    if (fj < -0.5 || fj > t->Ny - 0.5) continue;
                    if (fk < -0.5 || fk > t->Nz - 0.5) continue;

                    long out_idx = (long)oi * NN + oj * N + ok;

                    /* Template z coordinate (in template's own frame, for background subtraction) */
                    double z_tmpl = tz;

                    for (int c = 0; c < NCOLS; c++) {
                        float tmpl_val;
                        if (need_interp) {
                            tmpl_val = interp_template(t, c, fi, fj, fk);
                        } else {
                            /* Nearest-neighbor */
                            int ti = (int)round(fi);
                            int tj = (int)round(fj);
                            int tk = (int)round(fk);
                            if (ti < 0) ti = 0; if (ti >= (int)t->Nx) ti = t->Nx - 1;
                            if (tj < 0) tj = 0; if (tj >= (int)t->Ny) tj = t->Ny - 1;
                            if (tk < 0) tk = 0; if (tk >= (int)t->Nz) tk = t->Nz - 1;
                            tmpl_val = t->data[(uint64_t)c * tN3 + (long)ti * tNN + tj * t->Nz + tk];
                        }

                        /* Subtract template background to get perturbation only.
                         * The template contains background + baryon. We want just the baryon.
                         * Background in template: A_bg * cos(k_tmpl * z_tmpl + 2*pi*a/3)
                         * for phi fields (c=0,1,2), and omega_tmpl * A_bg * sin(...) for
                         * phi velocities (c=6,7,8). Theta fields (c=3,4,5,9,10,11) have
                         * no simple background to subtract. */
                        double bg = 0;
                        if (c < 3) {
                            /* phi_a field */
                            double phase = k_tmpl * z_tmpl + 2.0 * PI * c / 3.0;
                            bg = A_bg * cos(phase);
                        } else if (c >= 6 && c < 9) {
                            /* phi_vel_a */
                            int a = c - 6;
                            double phase = k_tmpl * z_tmpl + 2.0 * PI * a / 3.0;
                            bg = omega_bg * A_bg * sin(phase);
                        }
                        /* theta and theta_vel: bg=0, no subtraction needed */

                        double perturbation = (double)tmpl_val - bg;
                        fields[c][out_idx] += perturbation;
                    }
                }
            }
        }
    }

    /* Step 3: Write output SFA */
    fprintf(stderr, "gen_composite: writing %s...\n", outpath);

    double dt = 0.025 * dx;  /* standard dt factor */
    uint8_t sfa_dtype = (precision == 0) ? SFA_F16 : (precision == 1) ? SFA_F32 : SFA_F64;
    SFA *sfa = sfa_create(outpath, N, N, N, L, L, L, dt);

    /* Embed metadata */
    char buf_N[32], buf_L[32], buf_m[32], buf_Abg[32], buf_np[32];
    char buf_delta[64], buf_grad_h[32], buf_grad_l[32];
    snprintf(buf_N, 32, "%d", N);
    snprintf(buf_L, 32, "%.6f", L);
    snprintf(buf_m, 32, "%.6f", sqrt(m2));
    snprintf(buf_Abg, 32, "%.6f", A_bg);
    snprintf(buf_np, 32, "%d", n_particles);
    snprintf(buf_delta, 64, "%.6f,%.6f,%.6f", delta[0], delta[1], delta[2]);

    const char *keys[32], *vals[32];
    int nkv = 0;
    keys[nkv] = "generator"; vals[nkv] = "gen_composite"; nkv++;
    keys[nkv] = "N"; vals[nkv] = buf_N; nkv++;
    keys[nkv] = "L"; vals[nkv] = buf_L; nkv++;
    keys[nkv] = "m"; vals[nkv] = buf_m; nkv++;
    keys[nkv] = "A_bg"; vals[nkv] = buf_Abg; nkv++;
    keys[nkv] = "delta"; vals[nkv] = buf_delta; nkv++;
    keys[nkv] = "n_particles"; vals[nkv] = buf_np; nkv++;

    /* Per-particle metadata */
    char pbuf_tmpl[MAX_PARTICLES][512];
    char pbuf_pos[MAX_PARTICLES][128];
    char pkey_tmpl[MAX_PARTICLES][32];
    char pkey_pos[MAX_PARTICLES][32];
    for (int p = 0; p < n_particles; p++) {
        snprintf(pkey_tmpl[p], 32, "template_%d", p);
        snprintf(pkey_pos[p], 32, "position_%d", p);
        /* Just store the filename, not full path */
        const char *basename = strrchr(particles[p].path, '/');
        basename = basename ? basename + 1 : particles[p].path;
        snprintf(pbuf_tmpl[p], 512, "%s", basename);
        snprintf(pbuf_pos[p], 128, "%.4f,%.4f,%.4f",
                 particles[p].cx, particles[p].cy, particles[p].cz);
        keys[nkv] = pkey_tmpl[p]; vals[nkv] = pbuf_tmpl[p]; nkv++;
        keys[nkv] = pkey_pos[p]; vals[nkv] = pbuf_pos[p]; nkv++;
    }

    if (use_gradient) {
        snprintf(buf_grad_h, 32, "%.6f", gradient_A_high);
        snprintf(buf_grad_l, 32, "%.6f", gradient_A_low);
        keys[nkv] = "gradient_A_high"; vals[nkv] = buf_grad_h; nkv++;
        keys[nkv] = "gradient_A_low"; vals[nkv] = buf_grad_l; nkv++;
    }

    sfa_add_kvmd(sfa, 0, 0xFFFFFFFF, 0xFFFFFFFF, keys, vals, nkv);

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

    /* Convert to output precision and write */
    if (precision == 2) {
        /* f64: pass directly */
        void *cols[NCOLS];
        for (int c = 0; c < NCOLS; c++) cols[c] = fields[c];
        sfa_write_frame(sfa, 0.0, cols);
    } else {
        void *cols[NCOLS];
        int es = (precision == 0) ? 2 : 4;
        for (int c = 0; c < NCOLS; c++) {
            double *src = fields[c];
            cols[c] = malloc(N3 * es);
            if (precision == 1) {
                float *p = (float *)cols[c];
                for (long i = 0; i < N3; i++) p[i] = (float)src[i];
            } else {
                /* f16 */
                uint16_t *p = (uint16_t *)cols[c];
                for (long i = 0; i < N3; i++) {
                    float f = (float)src[i];
                    uint32_t x; memcpy(&x, &f, 4);
                    uint16_t sign = (x >> 16) & 0x8000;
                    int exp = ((x >> 23) & 0xFF) - 127 + 15;
                    uint16_t mant = (x >> 13) & 0x3FF;
                    p[i] = (exp <= 0) ? sign : (exp >= 31) ? (sign | 0x7C00) : (sign | (exp << 10) | mant);
                }
            }
        }
        sfa_write_frame(sfa, 0.0, cols);
        for (int c = 0; c < NCOLS; c++) free(cols[c]);
    }

    sfa_close(sfa);
    fprintf(stderr, "gen_composite: wrote %s (%d^3, %s, %d particles)\n",
            outpath, N, (const char*[]){"f16","f32","f64"}[precision], n_particles);

    /* Cleanup */
    for (int c = 0; c < NCOLS; c++) free(fields[c]);
    for (int p = 0; p < n_particles; p++) free(particles[p].data);

    return 0;
}
