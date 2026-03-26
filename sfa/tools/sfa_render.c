/*  sfa_render.c — Render animated WebP from SFA simulation data
 *
 *  Produces a 4x3 grid of cross-section visualizations:
 *    Col 1: Fields        Col 2: Velocity       Col 3: Acceleration
 *    Row 1: phi           d(phi)/dt             d²(phi)/dt²
 *    Row 2: theta         d(theta)/dt           d²(theta)/dt²
 *    Row 3: |P|           |v_phi|               |a_phi|
 *    Row 4: COMPOSITE     COMPOSITE             COMPOSITE
 *
 *  Build: gcc -O3 -o sfa_render sfa_render.c -lzstd -lm -lwebp -lwebpmux -lsharpyuv
 *  Usage: ./sfa_render input.sfa -o output.webp [-fps 15] [-slice z] [-size 256] [-skip 1]
 */

#define SFA_IMPLEMENTATION
#include "../format/sfa.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <webp/encode.h>
#include <webp/mux.h>

/* ---- Bitmap font (8x8, minimal character set) ---- */

static const uint8_t font8x8[][8] = {
    /* space (0x20=32) */ {0,0,0,0,0,0,0,0},
    /* ! */ {0x18,0x18,0x18,0x18,0x18,0x00,0x18,0x00},
    /* " */ {0x6C,0x6C,0x24,0x00,0x00,0x00,0x00,0x00},
    /* # */ {0x00,0x24,0x7E,0x24,0x7E,0x24,0x00,0x00},
    /* $ */ {0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00},
    /* % */ {0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00},
    /* & */ {0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00},
    /* ' */ {0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00},
    /* ( */ {0x0C,0x18,0x30,0x30,0x30,0x18,0x0C,0x00},
    /* ) */ {0x30,0x18,0x0C,0x0C,0x0C,0x18,0x30,0x00},
    /* * */ {0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00},
    /* + */ {0x00,0x18,0x18,0x7E,0x18,0x18,0x00,0x00},
    /* , */ {0x00,0x00,0x00,0x00,0x00,0x18,0x18,0x30},
    /* - */ {0x00,0x00,0x00,0x7E,0x00,0x00,0x00,0x00},
    /* . */ {0x00,0x00,0x00,0x00,0x00,0x18,0x18,0x00},
    /* / */ {0x06,0x0C,0x18,0x30,0x60,0xC0,0x80,0x00},
    /* 0 */ {0x3C,0x66,0x6E,0x76,0x66,0x66,0x3C,0x00},
    /* 1 */ {0x18,0x38,0x18,0x18,0x18,0x18,0x7E,0x00},
    /* 2 */ {0x3C,0x66,0x06,0x0C,0x18,0x30,0x7E,0x00},
    /* 3 */ {0x3C,0x66,0x06,0x1C,0x06,0x66,0x3C,0x00},
    /* 4 */ {0x0C,0x1C,0x3C,0x6C,0x7E,0x0C,0x0C,0x00},
    /* 5 */ {0x7E,0x60,0x7C,0x06,0x06,0x66,0x3C,0x00},
    /* 6 */ {0x1C,0x30,0x60,0x7C,0x66,0x66,0x3C,0x00},
    /* 7 */ {0x7E,0x06,0x0C,0x18,0x30,0x30,0x30,0x00},
    /* 8 */ {0x3C,0x66,0x66,0x3C,0x66,0x66,0x3C,0x00},
    /* 9 */ {0x3C,0x66,0x66,0x3E,0x06,0x0C,0x38,0x00},
    /* : */ {0x00,0x18,0x18,0x00,0x18,0x18,0x00,0x00},
    /* ; */ {0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00},
    /* < */ {0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00},
    /* = */ {0x00,0x00,0x7E,0x00,0x7E,0x00,0x00,0x00},
    /* > */ {0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00},
    /* ? */ {0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00},
    /* @ */ {0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00},
    /* A */ {0x3C,0x66,0x66,0x7E,0x66,0x66,0x66,0x00},
    /* B */ {0x7C,0x66,0x66,0x7C,0x66,0x66,0x7C,0x00},
    /* C */ {0x3C,0x66,0x60,0x60,0x60,0x66,0x3C,0x00},
    /* D */ {0x78,0x6C,0x66,0x66,0x66,0x6C,0x78,0x00},
    /* E */ {0x7E,0x60,0x60,0x7C,0x60,0x60,0x7E,0x00},
    /* F */ {0x7E,0x60,0x60,0x7C,0x60,0x60,0x60,0x00},
    /* G */ {0x3C,0x66,0x60,0x6E,0x66,0x66,0x3E,0x00},
    /* H */ {0x66,0x66,0x66,0x7E,0x66,0x66,0x66,0x00},
    /* I */ {0x3C,0x18,0x18,0x18,0x18,0x18,0x3C,0x00},
    /* J */ {0x1E,0x0C,0x0C,0x0C,0x0C,0x6C,0x38,0x00},
    /* K */ {0x66,0x6C,0x78,0x70,0x78,0x6C,0x66,0x00},
    /* L */ {0x60,0x60,0x60,0x60,0x60,0x60,0x7E,0x00},
    /* M */ {0x63,0x77,0x7F,0x6B,0x63,0x63,0x63,0x00},
    /* N */ {0x66,0x76,0x7E,0x7E,0x6E,0x66,0x66,0x00},
    /* O */ {0x3C,0x66,0x66,0x66,0x66,0x66,0x3C,0x00},
    /* P */ {0x7C,0x66,0x66,0x7C,0x60,0x60,0x60,0x00},
    /* Q */ {0x3C,0x66,0x66,0x66,0x6A,0x6C,0x36,0x00},
    /* R */ {0x7C,0x66,0x66,0x7C,0x6C,0x66,0x66,0x00},
    /* S */ {0x3C,0x66,0x60,0x3C,0x06,0x66,0x3C,0x00},
    /* T */ {0x7E,0x18,0x18,0x18,0x18,0x18,0x18,0x00},
    /* U */ {0x66,0x66,0x66,0x66,0x66,0x66,0x3C,0x00},
    /* V */ {0x66,0x66,0x66,0x66,0x66,0x3C,0x18,0x00},
    /* W */ {0x63,0x63,0x63,0x6B,0x7F,0x77,0x63,0x00},
    /* X */ {0x66,0x66,0x3C,0x18,0x3C,0x66,0x66,0x00},
    /* Y */ {0x66,0x66,0x66,0x3C,0x18,0x18,0x18,0x00},
    /* Z */ {0x7E,0x06,0x0C,0x18,0x30,0x60,0x7E,0x00},
    /* [ */ {0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00},
    /* \ */ {0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00},
    /* ] */ {0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00},
    /* ^ */ {0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00},
    /* _ */ {0x00,0x00,0x00,0x00,0x00,0x00,0x7E,0x00},
};

static const uint8_t *font_glyph(char c) {
    if (c >= ' ' && c <= '_') return font8x8[c - ' '];
    /* lowercase: map to uppercase */
    if (c >= 'a' && c <= 'z') return font8x8[c - 'a' + ('A' - ' ')];
    if (c == '|') return font8x8['I' - ' ']; /* use I glyph for | */
    return font8x8[0]; /* space */
}

/* Draw text at (x0, y0) into ARGB buffer. Color is 0xAARRGGBB. */
static void draw_text(uint32_t *argb, int stride, int img_w, int img_h,
                      int x0, int y0, const char *text, uint32_t color) {
    for (int ci = 0; text[ci]; ci++) {
        const uint8_t *g = font_glyph(text[ci]);
        int cx = x0 + ci * 8;
        for (int row = 0; row < 8; row++) {
            int py = y0 + row;
            if (py < 0 || py >= img_h) continue;
            uint8_t bits = g[row];
            for (int col = 0; col < 8; col++) {
                int px = cx + col;
                if (px < 0 || px >= img_w) continue;
                if (bits & (0x80 >> col)) {
                    argb[py * stride + px] = color;
                } else {
                    /* semi-transparent shadow for readability */
                    uint32_t cur = argb[py * stride + px];
                    uint32_t r = ((cur >> 16) & 0xFF) * 3 / 4;
                    uint32_t gc = ((cur >> 8) & 0xFF) * 3 / 4;
                    uint32_t b = (cur & 0xFF) * 3 / 4;
                    /* only darken in the label region (first 10 rows) */
                    if (row < 10 && col < 8)
                        ; /* don't darken outside glyph pixels */
                }
            }
        }
    }
}

/* Draw text with a dark background band for readability */
static void draw_label(uint32_t *argb, int stride, int img_w, int img_h,
                       int x0, int y0, const char *text) {
    int len = (int)strlen(text);
    int tw = len * 8 + 4;
    /* darken background behind text */
    for (int row = 0; row < 12; row++) {
        int py = y0 + row - 1;
        if (py < 0 || py >= img_h) continue;
        for (int col = 0; col < tw; col++) {
            int px = x0 + col - 2;
            if (px < 0 || px >= img_w) continue;
            uint32_t cur = argb[py * stride + px];
            uint32_t r = ((cur >> 16) & 0xFF) / 3;
            uint32_t g = ((cur >> 8) & 0xFF) / 3;
            uint32_t b = (cur & 0xFF) / 3;
            argb[py * stride + px] = 0xFF000000u | (r << 16) | (g << 8) | b;
        }
    }
    draw_text(argb, stride, img_w, img_h, x0, y0, text, 0xFFFFFFFFu);
}

/* ---- Half-float conversion ---- */

static inline double f16_to_f64(uint16_t h) {
    uint16_t sign = h & 0x8000;
    int exp = (h >> 10) & 0x1F;
    uint16_t mant = h & 0x3FF;
    if (exp == 0) return 0.0;
    if (exp == 31) return sign ? -1e30 : 1e30;
    float fv;
    uint32_t x = ((uint32_t)sign << 16) | ((uint32_t)(exp - 15 + 127) << 23) | ((uint32_t)mant << 13);
    memcpy(&fv, &x, 4);
    return (double)fv;
}

/* ---- Extract a column as float array ---- */

static float *extract_col_f(void *buf, SFA *sfa, int sem, int comp, long N3) {
    float *arr = (float *)calloc(N3, sizeof(float));
    if (!arr) return NULL;
    uint64_t off = 0;
    for (uint32_t c = 0; c < sfa->n_columns; c++) {
        int dt = sfa->columns[c].dtype;
        int s = sfa->columns[c].semantic;
        int cp = sfa->columns[c].component;
        int es = sfa_dtype_size[dt];
        if (s == sem && cp == comp) {
            uint8_t *src = (uint8_t *)buf + off;
            if (dt == SFA_F64)
                for (long i = 0; i < N3; i++) arr[i] = (float)((double *)src)[i];
            else if (dt == SFA_F32)
                memcpy(arr, src, N3 * 4);
            else if (dt == SFA_F16)
                for (long i = 0; i < N3; i++) arr[i] = (float)f16_to_f64(((uint16_t *)src)[i]);
            return arr;
        }
        off += (uint64_t)N3 * es;
    }
    return arr; /* zeros if not found */
}

/* ---- White-hot colormap: black -> blue -> cyan -> white ---- */

static void colormap_whitehot(float val, uint8_t *r, uint8_t *g, uint8_t *b) {
    /* val in [0,1] */
    if (val < 0) val = 0;
    if (val > 1) val = 1;
    if (val < 0.333f) {
        float t = val / 0.333f;
        *r = 0; *g = 0; *b = (uint8_t)(t * 255);
    } else if (val < 0.667f) {
        float t = (val - 0.333f) / 0.334f;
        *r = 0; *g = (uint8_t)(t * 255); *b = 255;
    } else {
        float t = (val - 0.667f) / 0.333f;
        *r = (uint8_t)(t * 255); *g = 255; *b = 255;
    }
}

/* ---- Main ---- */

int main(int argc, char **argv) {
    if (argc < 2) {
        fprintf(stderr, "Usage: %s input.sfa -o output.webp [-fps 15] [-slice z|y|x] [-size 256] [-skip N]\n", argv[0]);
        return 1;
    }

    const char *input_path = argv[1];
    const char *output_path = "output.webp";
    int fps = 15;
    int slice_axis = 2; /* 0=x, 1=y, 2=z */
    int sub_size = 256;
    int skip = 1;

    for (int i = 2; i < argc; i++) {
        if (!strcmp(argv[i], "-o") && i + 1 < argc) output_path = argv[++i];
        else if (!strcmp(argv[i], "-fps") && i + 1 < argc) fps = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-slice") && i + 1 < argc) {
            i++;
            if (argv[i][0] == 'x') slice_axis = 0;
            else if (argv[i][0] == 'y') slice_axis = 1;
            else slice_axis = 2;
        }
        else if (!strcmp(argv[i], "-size") && i + 1 < argc) sub_size = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-skip") && i + 1 < argc) skip = atoi(argv[++i]);
    }
    if (skip < 1) skip = 1;

    /* Open SFA */
    SFA *sfa = sfa_open(input_path);
    if (!sfa) { fprintf(stderr, "Cannot open %s\n", input_path); return 1; }

    int N = (int)sfa->Nx;
    int Ny = (int)sfa->Ny;
    int Nz = (int)sfa->Nz;
    double L = sfa->Lx;
    long N3 = (long)N * Ny * Nz;
    long NN = (long)Ny * Nz;
    double dx = 2.0 * L / (N - 1);
    double idx2 = 1.0 / (dx * dx);
    double idx1 = 1.0 / (2.0 * dx);
    int total_frames = (int)sfa->total_frames;

    /* Physics parameters (defaults, overridden by KVMD) */
    double MASS2 = 2.25, MU = -41.345, KAPPA = 50.0, ETA = 0.5;
    double DAMP_WIDTH = 3.0;

    SFA_KVMDSet kv[16];
    int nkv = sfa_read_kvmd(sfa, kv, 16);
    for (int s = 0; s < nkv; s++)
        for (int p = 0; p < kv[s].n_pairs; p++) {
            if (!strcmp(kv[s].keys[p], "mu")) MU = atof(kv[s].values[p]);
            if (!strcmp(kv[s].keys[p], "kappa")) KAPPA = atof(kv[s].values[p]);
            if (!strcmp(kv[s].keys[p], "m")) { double m = atof(kv[s].values[p]); MASS2 = m * m; }
            if (!strcmp(kv[s].keys[p], "eta")) ETA = atof(kv[s].values[p]);
            if (!strcmp(kv[s].keys[p], "damp_width")) DAMP_WIDTH = atof(kv[s].values[p]);
        }

    /* Boundary margin in grid cells — voxels within this distance of any edge are masked */
    int boundary_margin = (DAMP_WIDTH > 0) ? (int)(DAMP_WIDTH / dx + 0.5) : 0;
    if (boundary_margin < 0) boundary_margin = 0;

    printf("SFA: %s (%dx%dx%d, L=%.1f, %d frames)\n", input_path, N, Ny, Nz, L, total_frames);
    printf("Physics: m^2=%.4f mu=%.3f kappa=%.1f eta=%.3f\n", MASS2, MU, KAPPA, ETA);
    printf("Output: %s (%d fps, skip=%d, slice=%c, sub=%d)\n", output_path, fps, skip,
           "xyz"[slice_axis], sub_size);

    /* Count output frames (need pairs for velocity) */
    int n_input = 0;
    for (int f = 0; f < total_frames; f += skip) n_input++;
    int n_output = n_input - 1; /* velocity needs pairs */
    if (n_output < 1) { fprintf(stderr, "Need at least 2 frames\n"); sfa_close(sfa); return 1; }
    printf("Input frames: %d, output frames: %d\n", n_input, n_output);

    /* Grid dimensions for slicing */
    int slice_dim[3]; /* dimensions of the 2D slice */
    if (slice_axis == 0) { slice_dim[0] = Ny; slice_dim[1] = Nz; }
    else if (slice_axis == 1) { slice_dim[0] = N; slice_dim[1] = Nz; }
    else { slice_dim[0] = N; slice_dim[1] = Ny; }

    /* Image dimensions */
    int grid_cols = 3, grid_rows = 4;
    int img_w = grid_cols * sub_size;
    int img_h = grid_rows * sub_size;

    /* 12 sub-images per output frame:
     * [0]  phi_rms      [1]  dphi_rms/dt     [2]  a_phi_rms
     * [3]  theta_rms    [4]  dtheta_rms/dt   [5]  a_theta_rms
     * [6]  |P|          [7]  |v_phi|         [8]  |a_phi|
     * [9]  composite    [10] composite       [11] composite
     */
    #define N_SUBS 9 /* only 9 unique quantities; row 4 composites from rows 1-3 */

    /* Pass 1: Read all frames, compute all slice data, find global min/max */
    printf("Pass 1: Reading frames and computing slice data...\n");

    int slice_pixels = slice_dim[0] * slice_dim[1];
    /* Store all slice data: n_output frames x 9 quantities x slice_pixels */
    float **slice_data = (float **)malloc(n_output * N_SUBS * sizeof(float *));
    for (int i = 0; i < n_output * N_SUBS; i++)
        slice_data[i] = (float *)calloc(slice_pixels, sizeof(float));

    float gmin[N_SUBS], gmax[N_SUBS];
    for (int q = 0; q < N_SUBS; q++) { gmin[q] = 1e30f; gmax[q] = -1e30f; }

    /* We process pairs of consecutive (skip-stepped) frames */
    void *buf_cur = malloc(sfa->frame_bytes);
    void *buf_next = malloc(sfa->frame_bytes);
    if (!buf_cur || !buf_next) { fprintf(stderr, "Out of memory\n"); return 1; }

    /* Read first frame */
    int frame_indices[n_input];
    { int idx = 0; for (int f = 0; f < total_frames; f += skip) frame_indices[idx++] = f; }

    sfa_read_frame(sfa, frame_indices[0], buf_cur);

    for (int out_f = 0; out_f < n_output; out_f++) {
        int fi_cur = frame_indices[out_f];
        int fi_next = frame_indices[out_f + 1];
        double t_cur = sfa_frame_time(sfa, fi_cur);
        double t_next = sfa_frame_time(sfa, fi_next);
        double dt_pair = t_next - t_cur;
        if (dt_pair <= 0) dt_pair = 1e-10;

        /* Read next frame */
        sfa_read_frame(sfa, fi_next, buf_next);

        /* Extract fields from current frame */
        float *phi[3], *theta[3];
        for (int a = 0; a < 3; a++) {
            phi[a] = extract_col_f(buf_cur, sfa, SFA_POSITION, a, N3);
            theta[a] = extract_col_f(buf_cur, sfa, SFA_ANGLE, a, N3);
        }

        /* Extract fields from next frame (for velocity) */
        float *phi_next[3], *theta_next[3];
        for (int a = 0; a < 3; a++) {
            phi_next[a] = extract_col_f(buf_next, sfa, SFA_POSITION, a, N3);
            theta_next[a] = extract_col_f(buf_next, sfa, SFA_ANGLE, a, N3);
        }

        /* Check if velocity columns exist in SFA */
        int has_vel = 0;
        for (uint32_t c = 0; c < sfa->n_columns; c++)
            if (sfa->columns[c].semantic == SFA_VELOCITY) { has_vel = 1; break; }

        float *vphi[3] = {NULL, NULL, NULL}, *vtheta[3] = {NULL, NULL, NULL};
        if (has_vel) {
            /* Use stored velocities */
            for (int a = 0; a < 3; a++) {
                vphi[a] = extract_col_f(buf_cur, sfa, SFA_VELOCITY, a, N3);
                vtheta[a] = extract_col_f(buf_cur, sfa, SFA_VELOCITY, a + 3, N3);
            }
        }

        /* Find centroid weighted by |P| */
        double cm[3] = {0, 0, 0}, wt = 0;
        for (long idx = 0; idx < N3; idx++) {
            double P = fabs((double)phi[0][idx] * phi[1][idx] * phi[2][idx]);
            long i = idx / NN, j = (idx / Nz) % Ny, k = idx % Nz;
            cm[0] += (-L + i * dx) * P;
            cm[1] += (-L + j * dx) * P;
            cm[2] += (-L + k * dx) * P;
            wt += P;
        }
        if (wt > 0) { cm[0] /= wt; cm[1] /= wt; cm[2] /= wt; }

        /* Determine slice index */
        int slice_idx;
        if (slice_axis == 0) slice_idx = (int)((cm[0] + L) / dx + 0.5);
        else if (slice_axis == 1) slice_idx = (int)((cm[1] + L) / dx + 0.5);
        else slice_idx = (int)((cm[2] + L) / dx + 0.5);

        int max_idx = (slice_axis == 0) ? N : (slice_axis == 1) ? Ny : Nz;
        if (slice_idx < 1) slice_idx = 1;
        if (slice_idx >= max_idx - 1) slice_idx = max_idx - 2;

        /* For each pixel in the slice, compute all 9 quantities */
        float *sd = slice_data[out_f * N_SUBS]; /* base pointer for this frame */

        for (int u = 0; u < slice_dim[0]; u++) {
            for (int v = 0; v < slice_dim[1]; v++) {
                int sp = u * slice_dim[1] + v;

                /* Map (u, v, slice_idx) to 3D index */
                long i3, j3, k3;
                if (slice_axis == 0) { i3 = slice_idx; j3 = u; k3 = v; }
                else if (slice_axis == 1) { i3 = u; j3 = slice_idx; k3 = v; }
                else { i3 = u; j3 = v; k3 = slice_idx; }
                long idx = i3 * NN + j3 * Nz + k3;

                /* Mask absorbing boundary region */
                if (boundary_margin > 0 &&
                    (i3 < boundary_margin || i3 >= N - boundary_margin ||
                     j3 < boundary_margin || j3 >= Ny - boundary_margin ||
                     k3 < boundary_margin || k3 >= Nz - boundary_margin)) {
                    for (int q = 0; q < N_SUBS; q++)
                        slice_data[out_f * N_SUBS + q][sp] = 0;
                    continue;
                }

                /* Q0: phi_rms = sqrt(phi_x^2 + phi_y^2 + phi_z^2) */
                float prms = 0;
                for (int a = 0; a < 3; a++) prms += phi[a][idx] * phi[a][idx];
                prms = sqrtf(prms);
                slice_data[out_f * N_SUBS + 0][sp] = prms;

                /* Q3: theta_rms */
                float trms = 0;
                for (int a = 0; a < 3; a++) trms += theta[a][idx] * theta[a][idx];
                trms = sqrtf(trms);
                slice_data[out_f * N_SUBS + 3][sp] = trms;

                /* Q6: |P| = |phi_x * phi_y * phi_z| */
                float P = fabsf(phi[0][idx] * phi[1][idx] * phi[2][idx]);
                slice_data[out_f * N_SUBS + 6][sp] = P;

                /* Velocity: Q1 = dphi_rms/dt, Q4 = dtheta_rms/dt, Q7 = |v_phi| */
                float dprms = 0, dtrms = 0, vmag2 = 0;
                if (has_vel) {
                    for (int a = 0; a < 3; a++) {
                        float vp = vphi[a][idx];
                        float vt = vtheta[a] ? vtheta[a][idx] : 0;
                        dprms += vp * vp;
                        dtrms += vt * vt;
                        vmag2 += vp * vp;
                    }
                } else {
                    for (int a = 0; a < 3; a++) {
                        float dp = (phi_next[a][idx] - phi[a][idx]) / (float)dt_pair;
                        float dt_f = (theta_next[a][idx] - theta[a][idx]) / (float)dt_pair;
                        dprms += dp * dp;
                        dtrms += dt_f * dt_f;
                        vmag2 += dp * dp;
                    }
                }
                slice_data[out_f * N_SUBS + 1][sp] = sqrtf(dprms);
                slice_data[out_f * N_SUBS + 4][sp] = sqrtf(dtrms);
                slice_data[out_f * N_SUBS + 7][sp] = sqrtf(vmag2);

                /* Acceleration (from EOM): Q2, Q5, Q8 */
                /* Need boundary check */
                if (i3 < 1 || i3 >= N - 1 || j3 < 1 || j3 >= Ny - 1 || k3 < 1 || k3 >= Nz - 1) {
                    slice_data[out_f * N_SUBS + 2][sp] = 0;
                    slice_data[out_f * N_SUBS + 5][sp] = 0;
                    slice_data[out_f * N_SUBS + 8][sp] = 0;
                    continue;
                }

                long nip = (i3 + 1) * NN + j3 * Nz + k3;
                long nim = (i3 - 1) * NN + j3 * Nz + k3;
                long njp = i3 * NN + (j3 + 1) * Nz + k3;
                long njm = i3 * NN + (j3 - 1) * Nz + k3;
                long nkp = i3 * NN + j3 * Nz + (k3 + 1);
                long nkm = i3 * NN + j3 * Nz + (k3 - 1);

                float p0 = phi[0][idx], p1 = phi[1][idx], p2 = phi[2][idx];
                float Pv = p0 * p1 * p2;
                float P2 = Pv * Pv;
                float denom = (1.0f + (float)KAPPA * P2);
                float dVdP = (float)MU * Pv / (denom * denom);

                float a_phi2 = 0, a_theta2 = 0, a_phi_mag2 = 0;
                for (int a = 0; a < 3; a++) {
                    /* Laplacian phi */
                    float lap = (phi[a][nip] + phi[a][nim] + phi[a][njp] + phi[a][njm]
                                + phi[a][nkp] + phi[a][nkm] - 6.0f * phi[a][idx]) * (float)idx2;
                    float mass_f = -(float)MASS2 * phi[a][idx];
                    float dPda = (a == 0) ? p1 * p2 : (a == 1) ? p0 * p2 : p0 * p1;
                    float pot_f = -dVdP * dPda;

                    /* curl(theta) for phi equation */
                    float curl_phi = 0;
                    if (a == 0)
                        curl_phi = (float)ETA * (theta[2][njp] - theta[2][njm] - theta[1][nkp] + theta[1][nkm]) * (float)idx1;
                    else if (a == 1)
                        curl_phi = (float)ETA * (theta[0][nkp] - theta[0][nkm] - theta[2][nip] + theta[2][nim]) * (float)idx1;
                    else
                        curl_phi = (float)ETA * (theta[1][nip] - theta[1][nim] - theta[0][njp] + theta[0][njm]) * (float)idx1;

                    float af = lap + mass_f + pot_f + curl_phi;
                    a_phi2 += af * af;

                    /* Laplacian theta + curl(phi) for theta equation */
                    float lap_t = (theta[a][nip] + theta[a][nim] + theta[a][njp] + theta[a][njm]
                                  + theta[a][nkp] + theta[a][nkm] - 6.0f * theta[a][idx]) * (float)idx2;
                    float curl_th = 0;
                    if (a == 0)
                        curl_th = (float)ETA * (phi[2][njp] - phi[2][njm] - phi[1][nkp] + phi[1][nkm]) * (float)idx1;
                    else if (a == 1)
                        curl_th = (float)ETA * (phi[0][nkp] - phi[0][nkm] - phi[2][nip] + phi[2][nim]) * (float)idx1;
                    else
                        curl_th = (float)ETA * (phi[1][nip] - phi[1][nim] - phi[0][njp] + phi[0][njm]) * (float)idx1;

                    float at = lap_t + curl_th; /* m_theta^2 = 0 */
                    a_theta2 += at * at;
                    a_phi_mag2 += af * af; /* same as a_phi2 */
                }

                slice_data[out_f * N_SUBS + 2][sp] = sqrtf(a_phi2);
                slice_data[out_f * N_SUBS + 5][sp] = sqrtf(a_theta2);
                slice_data[out_f * N_SUBS + 8][sp] = sqrtf(a_phi_mag2);
            }
        }

        /* Update global min/max */
        for (int q = 0; q < N_SUBS; q++) {
            float *sd_q = slice_data[out_f * N_SUBS + q];
            for (int p = 0; p < slice_pixels; p++) {
                if (sd_q[p] < gmin[q]) gmin[q] = sd_q[p];
                if (sd_q[p] > gmax[q]) gmax[q] = sd_q[p];
            }
        }

        /* Free field arrays */
        for (int a = 0; a < 3; a++) {
            free(phi[a]); free(theta[a]);
            free(phi_next[a]); free(theta_next[a]);
            if (has_vel) { free(vphi[a]); if (vtheta[a]) free(vtheta[a]); }
        }

        /* Swap: next becomes current */
        void *tmp = buf_cur; buf_cur = buf_next; buf_next = tmp;

        if ((out_f + 1) % 10 == 0 || out_f == n_output - 1)
            printf("  Frame %d/%d (t=%.2f)\n", out_f + 1, n_output, t_cur);
    }

    free(buf_cur);
    free(buf_next);
    sfa_close(sfa);

    /* Print ranges */
    const char *qnames[] = {"phi_rms", "dphi/dt", "a_phi", "theta_rms", "dtheta/dt", "a_theta", "|P|", "|v_phi|", "|a_phi_mag|"};
    printf("\nGlobal ranges:\n");
    for (int q = 0; q < N_SUBS; q++)
        printf("  %12s: [%.4e, %.4e]\n", qnames[q], gmin[q], gmax[q]);

    /* Use 99.5th percentile for better contrast (avoid outlier domination) */
    float gmax_p[N_SUBS];
    for (int q = 0; q < N_SUBS; q++) {
        /* Collect histogram for percentile */
        int nbins = 1000;
        long hist[1000];
        memset(hist, 0, sizeof(hist));
        float range = gmax[q] - gmin[q];
        if (range < 1e-20f) range = 1e-20f;
        long total = 0;
        for (int f = 0; f < n_output; f++) {
            float *sd_q = slice_data[f * N_SUBS + q];
            for (int p = 0; p < slice_pixels; p++) {
                int b = (int)((sd_q[p] - gmin[q]) / range * (nbins - 1));
                if (b < 0) b = 0; if (b >= nbins) b = nbins - 1;
                hist[b]++;
                total++;
            }
        }
        /* Find 99.5th percentile */
        long cum = 0;
        long target = (long)(total * 0.995);
        gmax_p[q] = gmax[q];
        for (int b = 0; b < nbins; b++) {
            cum += hist[b];
            if (cum >= target) {
                gmax_p[q] = gmin[q] + (b + 1) * range / nbins;
                break;
            }
        }
        if (gmax_p[q] <= gmin[q]) gmax_p[q] = gmax[q];
    }

    printf("\n99.5th percentile max:\n");
    for (int q = 0; q < N_SUBS; q++)
        printf("  %12s: %.4e\n", qnames[q], gmax_p[q]);

    /* Pass 2: Render frames */
    printf("\nPass 2: Rendering %d frames to WebP...\n", n_output);

    int frame_dur_ms = 1000 / fps;

    WebPAnimEncoderOptions enc_options;
    WebPAnimEncoderOptionsInit(&enc_options);
    enc_options.allow_mixed = 1;
    WebPAnimEncoder *enc = WebPAnimEncoderNew(img_w, img_h, &enc_options);
    if (!enc) { fprintf(stderr, "Failed to create WebP encoder\n"); return 1; }

    /* Sub-image labels */
    const char *labels[12] = {
        "PHI", "dPHI/dt", "d2PHI/dt2",
        "THETA", "dTHETA/dt", "d2THETA/dt2",
        "|P|", "|V PHI|", "|A PHI|",
        "COMPOSITE", "COMPOSITE", "COMPOSITE"
    };

    /* Reopen SFA for timestamps */
    SFA *sfa2 = sfa_open(input_path);

    for (int out_f = 0; out_f < n_output; out_f++) {
        /* Allocate ARGB image */
        uint32_t *argb = (uint32_t *)calloc(img_w * img_h, sizeof(uint32_t));

        /* Fill with black + alpha */
        for (int p = 0; p < img_w * img_h; p++)
            argb[p] = 0xFF000000u;

        /* Render each sub-image (rows 0-2, cols 0-2) */
        for (int row = 0; row < 3; row++) {
            for (int col = 0; col < 3; col++) {
                int q = row * 3 + col;
                float *sd_q = slice_data[out_f * N_SUBS + q];
                float qmin = gmin[q];
                float qrange = gmax_p[q] - qmin;
                if (qrange < 1e-20f) qrange = 1e-20f;

                int x0 = col * sub_size;
                int y0 = row * sub_size;

                for (int py = 0; py < sub_size; py++) {
                    for (int px = 0; px < sub_size; px++) {
                        /* Map pixel to slice coordinates (nearest neighbor) */
                        int su = py * slice_dim[0] / sub_size;
                        int sv = px * slice_dim[1] / sub_size;
                        if (su >= slice_dim[0]) su = slice_dim[0] - 1;
                        if (sv >= slice_dim[1]) sv = slice_dim[1] - 1;
                        int sp = su * slice_dim[1] + sv;

                        float val = (sd_q[sp] - qmin) / qrange;
                        uint8_t r, g, b;
                        colormap_whitehot(val, &r, &g, &b);
                        argb[(y0 + py) * img_w + (x0 + px)] = 0xFF000000u | ((uint32_t)r << 16) | ((uint32_t)g << 8) | b;
                    }
                }
            }
        }

        /* Row 4: Composites (RGB from rows 0-2) */
        for (int col = 0; col < 3; col++) {
            int x0 = col * sub_size;
            int y0 = 3 * sub_size;

            float *sd_r = slice_data[out_f * N_SUBS + 0 * 3 + col]; /* phi row */
            float *sd_g = slice_data[out_f * N_SUBS + 1 * 3 + col]; /* theta row */
            float *sd_b = slice_data[out_f * N_SUBS + 2 * 3 + col]; /* |P| row */

            float rmin = gmin[0 * 3 + col], rrange = gmax_p[0 * 3 + col] - rmin;
            float gmin_c = gmin[1 * 3 + col], grange = gmax_p[1 * 3 + col] - gmin_c;
            float bmin = gmin[2 * 3 + col], brange = gmax_p[2 * 3 + col] - bmin;
            if (rrange < 1e-20f) rrange = 1e-20f;
            if (grange < 1e-20f) grange = 1e-20f;
            if (brange < 1e-20f) brange = 1e-20f;

            for (int py = 0; py < sub_size; py++) {
                for (int px = 0; px < sub_size; px++) {
                    int su = py * slice_dim[0] / sub_size;
                    int sv = px * slice_dim[1] / sub_size;
                    if (su >= slice_dim[0]) su = slice_dim[0] - 1;
                    if (sv >= slice_dim[1]) sv = slice_dim[1] - 1;
                    int sp = su * slice_dim[1] + sv;

                    float rv = (sd_r[sp] - rmin) / rrange;
                    float gv = (sd_g[sp] - gmin_c) / grange;
                    float bv = (sd_b[sp] - bmin) / brange;
                    if (rv < 0) rv = 0; if (rv > 1) rv = 1;
                    if (gv < 0) gv = 0; if (gv > 1) gv = 1;
                    if (bv < 0) bv = 0; if (bv > 1) bv = 1;

                    uint8_t r8 = (uint8_t)(rv * 255);
                    uint8_t g8 = (uint8_t)(gv * 255);
                    uint8_t b8 = (uint8_t)(bv * 255);
                    argb[(y0 + py) * img_w + (x0 + px)] = 0xFF000000u | ((uint32_t)r8 << 16) | ((uint32_t)g8 << 8) | b8;
                }
            }
        }

        /* Draw labels */
        for (int row = 0; row < 4; row++)
            for (int col = 0; col < 3; col++)
                draw_label(argb, img_w, img_w, img_h,
                           col * sub_size + 4, row * sub_size + 4,
                           labels[row * 3 + col]);

        /* Draw time label */
        double t = sfa_frame_time(sfa2, frame_indices[out_f]);
        char tbuf[32];
        snprintf(tbuf, sizeof(tbuf), "t=%.2f", t);
        draw_label(argb, img_w, img_w, img_h, 4, img_h - 16, tbuf);

        /* Encode frame */
        WebPConfig config;
        WebPConfigInit(&config);
        config.quality = 80;
        config.method = 4;

        WebPPicture pic;
        WebPPictureInit(&pic);
        pic.width = img_w;
        pic.height = img_h;
        pic.use_argb = 1;
        WebPPictureAlloc(&pic);

        /* Copy ARGB data */
        memcpy(pic.argb, argb, img_w * img_h * sizeof(uint32_t));

        int timestamp_ms = out_f * frame_dur_ms;
        if (!WebPAnimEncoderAdd(enc, &pic, timestamp_ms, &config)) {
            fprintf(stderr, "WebP encode error at frame %d: %s\n", out_f, WebPAnimEncoderGetError(enc));
            WebPPictureFree(&pic);
            free(argb);
            break;
        }
        WebPPictureFree(&pic);
        free(argb);

        if ((out_f + 1) % 10 == 0 || out_f == n_output - 1)
            printf("  Encoded frame %d/%d\n", out_f + 1, n_output);
    }

    /* Signal end */
    int final_ts = n_output * frame_dur_ms;
    WebPAnimEncoderAdd(enc, NULL, final_ts, NULL);

    /* Assemble */
    WebPData webp_data;
    WebPDataInit(&webp_data);
    if (!WebPAnimEncoderAssemble(enc, &webp_data)) {
        fprintf(stderr, "WebP assembly failed: %s\n", WebPAnimEncoderGetError(enc));
        WebPAnimEncoderDelete(enc);
        return 1;
    }

    /* Write to file */
    FILE *fp = fopen(output_path, "wb");
    if (!fp) { fprintf(stderr, "Cannot open %s for writing\n", output_path); return 1; }
    fwrite(webp_data.bytes, 1, webp_data.size, fp);
    fclose(fp);

    printf("\nWrote %s (%.1f MB, %d frames)\n", output_path,
           (double)webp_data.size / (1024 * 1024), n_output);

    WebPDataClear(&webp_data);
    WebPAnimEncoderDelete(enc);
    sfa_close(sfa2);

    /* Free slice data */
    for (int i = 0; i < n_output * N_SUBS; i++) free(slice_data[i]);
    free(slice_data);

    return 0;
}
