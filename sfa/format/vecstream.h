/*  vecstream.h — Streaming vector compression for 3D volumetric fields
 *
 *  Single-header C library. #define VECSTREAM_IMPLEMENTATION in exactly ONE
 *  .c file before including. Requires: libzstd (-lzstd), libm (-lm).
 *
 *  Writer:
 *    VecStream *vs = vecstream_create("out.vecstream", Nx, Ny, Nz,
 *                                     Lx, Ly, Lz, dt, n_fields, block_size);
 *    vecstream_write_iframe(vs, time, field_idx, patches, n_patches);
 *    vecstream_write_pframe(vs, time, field_idx, indices, delta_coeffs,
 *                           n_deltas, coeffs_per_patch);
 *    vecstream_write_kframe(vs, time, field_idx, voxels, 1);  // 1 = f32
 *    vecstream_write_fourier(vs, field_idx, fourier_data, n_patches, cpp);
 *    vecstream_close(vs);
 *
 *  Reader:
 *    VecStream *vs = vecstream_open("out.vecstream");
 *    vecstream_frame_count(vs);
 *    vecstream_frame_time(vs, idx);
 *    vecstream_frame_type(vs, idx);
 *    vecstream_read_iframe(vs, idx, patches);
 *    vecstream_read_pframe(vs, idx, indices, delta_coeffs, &n_deltas);
 *    vecstream_read_kframe(vs, idx, voxels);
 *    vecstream_reconstruct(vs, frame_idx, field_idx, output_voxels);
 *    vecstream_close(vs);
 */

#ifndef VECSTREAM_H
#define VECSTREAM_H

#include <stdint.h>
#include <stdio.h>

/* ---- Constants ---- */

#define VS_VERSION        1
#define VS_MAGIC_FILE     "VECS"
#define VS_MAGIC_FRAME    "FRVS"
#define VS_MAGIC_FOURIER  "FLVS"
#define VS_MAGIC_INDEX    "IXVS"
#define VS_MAGIC_FOOTER   "VSND"

#define VS_FRAME_I  0
#define VS_FRAME_P  1
#define VS_FRAME_K  2

#define VS_FLAG_HAS_FOURIER  (1u << 0)

#define VS_MAX_ORDER   3
#define VS_NCOEFFS     64   /* (VS_MAX_ORDER+1)^3 = 4^3 */

/* ---- Structures ---- */

typedef struct {
    int16_t  origin_x, origin_y, origin_z;
    uint8_t  size_x, size_y, size_z;
    uint8_t  order;
    uint16_t n_coeffs;
    uint16_t max_err_f16;   /* float16-encoded max error */
    uint16_t rms_err_f16;   /* float16-encoded rms error */
    float    coeffs[VS_NCOEFFS];
} VecPatch;

typedef struct {
    double   time;
    uint64_t offset;
    uint8_t  frame_type;
    uint8_t  field_idx;
    uint8_t  _pad[6];
} VecIndexEntry;

typedef struct {
    FILE        *fp;
    int          mode;           /* 0=read, 1=write */

    /* Header fields */
    uint32_t     version;
    uint32_t     Nx, Ny, Nz;
    double       Lx, Ly, Lz, dt;
    uint16_t     n_fields;
    uint16_t     block_size;
    uint32_t     n_frames;
    uint32_t     flags;

    /* Derived */
    uint32_t     blocks_per_axis_x, blocks_per_axis_y, blocks_per_axis_z;
    uint32_t     total_patches;

    /* Write state */
    VecIndexEntry *index;
    uint32_t     index_cap;
    uint64_t     fourier_offset;   /* 0 if no fourier layer */

    /* Read state */
    uint32_t     index_count;

#ifdef _OPENMP
    /* For thread-safe writes */
    void        *write_mutex;      /* omp_lock_t* */
#endif
} VecStream;

/* ---- API ---- */

/* Writer */
VecStream *vecstream_create(const char *path, uint32_t Nx, uint32_t Ny, uint32_t Nz,
                            double Lx, double Ly, double Lz, double dt,
                            uint16_t n_fields, uint16_t block_size);

int vecstream_write_iframe(VecStream *vs, double time, uint8_t field_idx,
                           const VecPatch *patches, uint32_t n_patches);

int vecstream_write_pframe(VecStream *vs, double time, uint8_t field_idx,
                           const uint32_t *indices, const float *delta_coeffs,
                           uint32_t n_deltas, uint16_t coeffs_per_patch);

int vecstream_write_kframe(VecStream *vs, double time, uint8_t field_idx,
                           const void *voxels, uint8_t dtype);

int vecstream_write_fourier(VecStream *vs, uint8_t field_idx,
                            const float *fourier_data, uint32_t n_patches,
                            uint16_t coeffs_per_patch);

void vecstream_close(VecStream *vs);

/* Reader */
VecStream *vecstream_open(const char *path);

uint32_t vecstream_frame_count(const VecStream *vs);
double   vecstream_frame_time(const VecStream *vs, uint32_t idx);
uint8_t  vecstream_frame_type(const VecStream *vs, uint32_t idx);
uint8_t  vecstream_frame_field(const VecStream *vs, uint32_t idx);

int vecstream_read_iframe(VecStream *vs, uint32_t idx,
                          VecPatch *patches, uint32_t *n_patches);

int vecstream_read_pframe(VecStream *vs, uint32_t idx, uint32_t *indices,
                          float *delta_coeffs, uint32_t *n_deltas,
                          uint16_t *coeffs_per_patch);

int vecstream_read_kframe(VecStream *vs, uint32_t idx, float *voxels);

int vecstream_reconstruct(VecStream *vs, uint32_t frame_idx, uint8_t field_idx,
                          float *output_voxels);

/* Utility: fit a patch from voxel data */
void vecstream_fit_patch(const float *field, uint32_t Nx, uint32_t Ny, uint32_t Nz,
                         int ox, int oy, int oz, int bsx, int bsy, int bsz,
                         int order, VecPatch *patch);

/* Utility: evaluate a patch at a local grid point */
float vecstream_eval_patch(const VecPatch *p, int bs_x, int bs_y, int bs_z,
                           int di, int dj, int dk);

/* Utility: reconstruct full voxel grid from patches */
void vecstream_patches_to_voxels(const VecPatch *patches, uint32_t n_patches,
                                 uint32_t Nx, uint32_t Ny, uint32_t Nz,
                                 float *output);

/* Float16 conversion helpers */
static inline uint16_t vs_f32_to_f16(float f);
static inline float    vs_f16_to_f32(uint16_t h);

#endif /* VECSTREAM_H */


/* ================================================================
   IMPLEMENTATION
   ================================================================ */

#ifdef VECSTREAM_IMPLEMENTATION

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <zstd.h>

#ifdef _OPENMP
#include <omp.h>
#endif

/* ---- Float16 conversion ---- */

static inline uint16_t vs_f32_to_f16(float f) {
    uint32_t x; memcpy(&x, &f, 4);
    uint16_t s = (x >> 16) & 0x8000;
    int e = ((x >> 23) & 0xFF) - 127 + 15;
    uint16_t m = (x >> 13) & 0x3FF;
    if (e <= 0) return s;
    if (e >= 31) return s | 0x7C00;
    return s | (uint16_t)(e << 10) | m;
}

static inline float vs_f16_to_f32(uint16_t h) {
    uint16_t s = h & 0x8000;
    int e = (h >> 10) & 0x1F;
    uint16_t m = h & 0x3FF;
    if (e == 0) return 0.0f;
    if (e == 31) return s ? -1e30f : 1e30f;
    float f;
    uint32_t x = ((uint32_t)s << 16) | ((uint32_t)(e - 15 + 127) << 23) | ((uint32_t)m << 13);
    memcpy(&f, &x, 4);
    return f;
}

/* ---- CRC32 (Castagnoli) ---- */

static uint32_t vs_crc32(const void *data, size_t len) {
    uint32_t crc = 0xFFFFFFFF;
    const uint8_t *p = (const uint8_t *)data;
    for (size_t i = 0; i < len; i++) {
        crc ^= p[i];
        for (int j = 0; j < 8; j++)
            crc = (crc >> 1) ^ (0x82F63B78 & -(crc & 1));
    }
    return ~crc;
}

/* ---- Zstd helpers ---- */

static uint8_t *vs_compress(const void *src, size_t src_size, size_t *out_size) {
    size_t bound = ZSTD_compressBound(src_size);
    uint8_t *dst = (uint8_t *)malloc(bound);
    if (!dst) return NULL;
    size_t cs = ZSTD_compress(dst, bound, src, src_size, 3);
    if (ZSTD_isError(cs)) { free(dst); return NULL; }
    *out_size = cs;
    return dst;
}

static uint8_t *vs_decompress(const void *src, size_t src_size,
                               size_t expected_size, size_t *out_size) {
    uint8_t *dst = (uint8_t *)malloc(expected_size);
    if (!dst) return NULL;
    size_t ds = ZSTD_decompress(dst, expected_size, src, src_size);
    if (ZSTD_isError(ds)) { free(dst); return NULL; }
    *out_size = ds;
    return dst;
}

/* ---- Writer ---- */

VecStream *vecstream_create(const char *path, uint32_t Nx, uint32_t Ny, uint32_t Nz,
                            double Lx, double Ly, double Lz, double dt,
                            uint16_t n_fields, uint16_t block_size) {
    VecStream *vs = (VecStream *)calloc(1, sizeof(VecStream));
    if (!vs) return NULL;

    vs->fp = fopen(path, "wb");
    if (!vs->fp) { free(vs); return NULL; }

    vs->mode = 1;
    vs->version = VS_VERSION;
    vs->Nx = Nx; vs->Ny = Ny; vs->Nz = Nz;
    vs->Lx = Lx; vs->Ly = Ly; vs->Lz = Lz;
    vs->dt = dt;
    vs->n_fields = n_fields;
    vs->block_size = block_size;
    vs->n_frames = 0;
    vs->flags = 0;

    vs->blocks_per_axis_x = (Nx + block_size - 1) / block_size;
    vs->blocks_per_axis_y = (Ny + block_size - 1) / block_size;
    vs->blocks_per_axis_z = (Nz + block_size - 1) / block_size;
    vs->total_patches = vs->blocks_per_axis_x * vs->blocks_per_axis_y * vs->blocks_per_axis_z;

    vs->index_cap = 1024;
    vs->index = (VecIndexEntry *)calloc(vs->index_cap, sizeof(VecIndexEntry));

#ifdef _OPENMP
    vs->write_mutex = malloc(sizeof(omp_lock_t));
    omp_init_lock((omp_lock_t *)vs->write_mutex);
#endif

    /* Write file header (64 bytes) */
    uint8_t hdr[64];
    memset(hdr, 0, 64);
    memcpy(hdr + 0, VS_MAGIC_FILE, 4);
    memcpy(hdr + 4, &vs->version, 4);
    memcpy(hdr + 8, &Nx, 4);
    memcpy(hdr + 12, &Ny, 4);
    memcpy(hdr + 16, &Nz, 4);
    memcpy(hdr + 20, &Lx, 8);
    memcpy(hdr + 28, &Ly, 8);
    memcpy(hdr + 36, &Lz, 8);
    memcpy(hdr + 44, &dt, 8);
    memcpy(hdr + 52, &n_fields, 2);
    memcpy(hdr + 54, &block_size, 2);
    memcpy(hdr + 56, &vs->n_frames, 4);
    memcpy(hdr + 60, &vs->flags, 4);
    fwrite(hdr, 1, 64, vs->fp);

    return vs;
}

static void vs_add_index_entry(VecStream *vs, double time, uint64_t offset,
                               uint8_t frame_type, uint8_t field_idx) {
    if (vs->n_frames >= vs->index_cap) {
        vs->index_cap *= 2;
        vs->index = (VecIndexEntry *)realloc(vs->index,
                                              vs->index_cap * sizeof(VecIndexEntry));
    }
    VecIndexEntry *e = &vs->index[vs->n_frames];
    e->time = time;
    e->offset = offset;
    e->frame_type = frame_type;
    e->field_idx = field_idx;
    memset(e->_pad, 0, 6);
    vs->n_frames++;
}

static int vs_write_frame_header(FILE *fp, uint8_t frame_type, uint8_t field_idx,
                                 double time, uint64_t compressed_size,
                                 uint64_t uncompressed_size) {
    uint8_t hdr[32];
    memset(hdr, 0, 32);
    memcpy(hdr + 0, VS_MAGIC_FRAME, 4);
    hdr[4] = frame_type;
    hdr[5] = field_idx;
    memcpy(hdr + 8, &time, 8);
    memcpy(hdr + 16, &compressed_size, 8);
    memcpy(hdr + 24, &uncompressed_size, 8);
    return fwrite(hdr, 1, 32, fp) == 32 ? 0 : -1;
}

int vecstream_write_iframe(VecStream *vs, double time, uint8_t field_idx,
                           const VecPatch *patches, uint32_t n_patches) {
    if (!vs || vs->mode != 1) return -1;

    /* Build uncompressed payload */
    /* Layout: n_patches(u32) + [PatchHeader(16) + coeffs(n_coeffs*4)] * n_patches */
    size_t patch_data_size = 16; /* PatchHeader */
    /* Each patch: 16-byte header + n_coeffs * 4 bytes of coefficients */
    /* But n_coeffs varies per patch, so compute total */
    size_t total_size = 4; /* n_patches */
    for (uint32_t i = 0; i < n_patches; i++)
        total_size += 16 + (size_t)patches[i].n_coeffs * 4;

    uint8_t *payload = (uint8_t *)malloc(total_size);
    if (!payload) return -1;

    uint8_t *p = payload;
    memcpy(p, &n_patches, 4); p += 4;

    for (uint32_t i = 0; i < n_patches; i++) {
        const VecPatch *vp = &patches[i];
        /* PatchHeader: 16 bytes */
        memcpy(p + 0, &vp->origin_x, 2);
        memcpy(p + 2, &vp->origin_y, 2);
        memcpy(p + 4, &vp->origin_z, 2);
        p[6] = vp->size_x;
        p[7] = vp->size_y;
        p[8] = vp->size_z;
        p[9] = vp->order;
        memcpy(p + 10, &vp->n_coeffs, 2);
        memcpy(p + 12, &vp->max_err_f16, 2);
        memcpy(p + 14, &vp->rms_err_f16, 2);
        p += 16;

        memcpy(p, vp->coeffs, vp->n_coeffs * 4);
        p += vp->n_coeffs * 4;
    }

    /* Compress */
    size_t comp_size;
    uint8_t *comp = vs_compress(payload, total_size, &comp_size);
    free(payload);
    if (!comp) return -1;

    /* Write */
#ifdef _OPENMP
    omp_set_lock((omp_lock_t *)vs->write_mutex);
#endif
    uint64_t offset = (uint64_t)ftello(vs->fp);
    vs_write_frame_header(vs->fp, VS_FRAME_I, field_idx, time, comp_size, total_size);
    fwrite(comp, 1, comp_size, vs->fp);
    vs_add_index_entry(vs, time, offset, VS_FRAME_I, field_idx);
#ifdef _OPENMP
    omp_unset_lock((omp_lock_t *)vs->write_mutex);
#endif

    free(comp);
    return 0;
}

int vecstream_write_pframe(VecStream *vs, double time, uint8_t field_idx,
                           const uint32_t *indices, const float *delta_coeffs,
                           uint32_t n_deltas, uint16_t coeffs_per_patch) {
    if (!vs || vs->mode != 1) return -1;

    /* Layout: n_deltas(u32) + coeffs_per_patch(u16) + pad(u16)
     *       + [patch_index(u32) + coeffs(cpp*4)] * n_deltas */
    size_t total_size = 4 + 2 + 2 + (size_t)n_deltas * (4 + (size_t)coeffs_per_patch * 4);

    uint8_t *payload = (uint8_t *)malloc(total_size);
    if (!payload) return -1;

    uint8_t *p = payload;
    memcpy(p, &n_deltas, 4); p += 4;
    memcpy(p, &coeffs_per_patch, 2); p += 2;
    uint16_t pad = 0;
    memcpy(p, &pad, 2); p += 2;

    for (uint32_t i = 0; i < n_deltas; i++) {
        memcpy(p, &indices[i], 4); p += 4;
        memcpy(p, &delta_coeffs[(size_t)i * coeffs_per_patch],
               coeffs_per_patch * 4);
        p += (size_t)coeffs_per_patch * 4;
    }

    size_t comp_size;
    uint8_t *comp = vs_compress(payload, total_size, &comp_size);
    free(payload);
    if (!comp) return -1;

#ifdef _OPENMP
    omp_set_lock((omp_lock_t *)vs->write_mutex);
#endif
    uint64_t offset = (uint64_t)ftello(vs->fp);
    vs_write_frame_header(vs->fp, VS_FRAME_P, field_idx, time, comp_size, total_size);
    fwrite(comp, 1, comp_size, vs->fp);
    vs_add_index_entry(vs, time, offset, VS_FRAME_P, field_idx);
#ifdef _OPENMP
    omp_unset_lock((omp_lock_t *)vs->write_mutex);
#endif

    free(comp);
    return 0;
}

int vecstream_write_kframe(VecStream *vs, double time, uint8_t field_idx,
                           const void *voxels, uint8_t dtype) {
    if (!vs || vs->mode != 1) return -1;

    uint64_t N3 = (uint64_t)vs->Nx * vs->Ny * vs->Nz;
    size_t voxel_size = N3 * 4; /* always stored as f32 */

    /* Convert to f32 if needed */
    float *f32_buf = NULL;
    const float *src;
    if (dtype == 1) { /* already f32 */
        src = (const float *)voxels;
    } else {
        f32_buf = (float *)malloc(voxel_size);
        if (!f32_buf) return -1;
        if (dtype == 0) { /* f16 */
            const uint16_t *h = (const uint16_t *)voxels;
            for (uint64_t i = 0; i < N3; i++)
                f32_buf[i] = vs_f16_to_f32(h[i]);
        } else if (dtype == 2) { /* f64 */
            const double *d = (const double *)voxels;
            for (uint64_t i = 0; i < N3; i++)
                f32_buf[i] = (float)d[i];
        }
        src = f32_buf;
    }

    size_t comp_size;
    uint8_t *comp = vs_compress(src, voxel_size, &comp_size);
    if (f32_buf) free(f32_buf);
    if (!comp) return -1;

#ifdef _OPENMP
    omp_set_lock((omp_lock_t *)vs->write_mutex);
#endif
    uint64_t offset = (uint64_t)ftello(vs->fp);
    vs_write_frame_header(vs->fp, VS_FRAME_K, field_idx, time, comp_size, voxel_size);
    fwrite(comp, 1, comp_size, vs->fp);
    vs_add_index_entry(vs, time, offset, VS_FRAME_K, field_idx);
#ifdef _OPENMP
    omp_unset_lock((omp_lock_t *)vs->write_mutex);
#endif

    free(comp);
    return 0;
}

int vecstream_write_fourier(VecStream *vs, uint8_t field_idx,
                            const float *fourier_data, uint32_t n_patches,
                            uint16_t coeffs_per_patch) {
    if (!vs || vs->mode != 1) return -1;

    vs->flags |= VS_FLAG_HAS_FOURIER;

#ifdef _OPENMP
    omp_set_lock((omp_lock_t *)vs->write_mutex);
#endif

    vs->fourier_offset = (uint64_t)ftello(vs->fp);

    /* Header: 16 bytes */
    uint8_t hdr[16];
    memset(hdr, 0, 16);
    memcpy(hdr + 0, VS_MAGIC_FOURIER, 4);
    hdr[4] = field_idx;
    hdr[5] = 1; /* n_harmonics = 1 for now */
    memcpy(hdr + 8, &n_patches, 4);
    memcpy(hdr + 12, &coeffs_per_patch, 2);
    fwrite(hdr, 1, 16, vs->fp);

    /* Data: 16 bytes per coeff per harmonic (amp, omega, phase, offset) */
    size_t data_size = (size_t)n_patches * coeffs_per_patch * 4 * sizeof(float);
    fwrite(fourier_data, 1, data_size, vs->fp);

#ifdef _OPENMP
    omp_unset_lock((omp_lock_t *)vs->write_mutex);
#endif

    return 0;
}

void vecstream_close(VecStream *vs) {
    if (!vs) return;

    if (vs->mode == 1) {
        /* Write frame index */
        uint64_t index_offset = (uint64_t)ftello(vs->fp);

        /* Index header: 16 bytes */
        uint8_t ix_hdr[16];
        memset(ix_hdr, 0, 16);
        memcpy(ix_hdr + 0, VS_MAGIC_INDEX, 4);
        memcpy(ix_hdr + 4, &vs->n_frames, 4);
        fwrite(ix_hdr, 1, 16, vs->fp);

        /* Index entries: 24 bytes each */
        for (uint32_t i = 0; i < vs->n_frames; i++) {
            uint8_t entry[24];
            memset(entry, 0, 24);
            memcpy(entry + 0, &vs->index[i].time, 8);
            memcpy(entry + 8, &vs->index[i].offset, 8);
            entry[16] = vs->index[i].frame_type;
            entry[17] = vs->index[i].field_idx;
            fwrite(entry, 1, 24, vs->fp);
        }

        /* Footer: 16 bytes */
        /* First, compute CRC32 of the file header */
        uint8_t file_hdr[64];
        fseeko(vs->fp, 0, SEEK_SET);
        if (fread(file_hdr, 1, 64, vs->fp) != 64) {
            /* Can't read header back, use 0 */
            memset(file_hdr, 0, 64);
        }
        fseeko(vs->fp, 0, SEEK_END);

        uint32_t checksum = vs_crc32(file_hdr, 64);
        uint8_t footer[16];
        memcpy(footer + 0, &index_offset, 8);
        memcpy(footer + 8, &checksum, 4);
        memcpy(footer + 12, VS_MAGIC_FOOTER, 4);
        fwrite(footer, 1, 16, vs->fp);

        /* Update n_frames and flags in file header */
        fseeko(vs->fp, 56, SEEK_SET);
        fwrite(&vs->n_frames, 4, 1, vs->fp);
        fwrite(&vs->flags, 4, 1, vs->fp);
    }

    if (vs->fp) fclose(vs->fp);
    if (vs->index) free(vs->index);

#ifdef _OPENMP
    if (vs->write_mutex) {
        omp_destroy_lock((omp_lock_t *)vs->write_mutex);
        free(vs->write_mutex);
    }
#endif

    free(vs);
}

/* ---- Reader ---- */

VecStream *vecstream_open(const char *path) {
    FILE *fp = fopen(path, "rb");
    if (!fp) return NULL;

    /* Read file header */
    uint8_t hdr[64];
    if (fread(hdr, 1, 64, fp) != 64) { fclose(fp); return NULL; }
    if (memcmp(hdr, VS_MAGIC_FILE, 4) != 0) { fclose(fp); return NULL; }

    VecStream *vs = (VecStream *)calloc(1, sizeof(VecStream));
    vs->fp = fp;
    vs->mode = 0;

    memcpy(&vs->version, hdr + 4, 4);
    memcpy(&vs->Nx, hdr + 8, 4);
    memcpy(&vs->Ny, hdr + 12, 4);
    memcpy(&vs->Nz, hdr + 16, 4);
    memcpy(&vs->Lx, hdr + 20, 8);
    memcpy(&vs->Ly, hdr + 28, 8);
    memcpy(&vs->Lz, hdr + 36, 8);
    memcpy(&vs->dt, hdr + 44, 8);
    memcpy(&vs->n_fields, hdr + 52, 2);
    memcpy(&vs->block_size, hdr + 54, 2);
    memcpy(&vs->n_frames, hdr + 56, 4);
    memcpy(&vs->flags, hdr + 60, 4);

    vs->blocks_per_axis_x = (vs->Nx + vs->block_size - 1) / vs->block_size;
    vs->blocks_per_axis_y = (vs->Ny + vs->block_size - 1) / vs->block_size;
    vs->blocks_per_axis_z = (vs->Nz + vs->block_size - 1) / vs->block_size;
    vs->total_patches = vs->blocks_per_axis_x * vs->blocks_per_axis_y * vs->blocks_per_axis_z;

    /* Read footer to find index */
    fseeko(fp, -16, SEEK_END);
    uint8_t footer[16];
    if (fread(footer, 1, 16, fp) != 16) { fclose(fp); free(vs); return NULL; }
    if (memcmp(footer + 12, VS_MAGIC_FOOTER, 4) != 0) { fclose(fp); free(vs); return NULL; }

    uint64_t index_offset;
    memcpy(&index_offset, footer, 8);

    /* Read index header */
    fseeko(fp, index_offset, SEEK_SET);
    uint8_t ix_hdr[16];
    if (fread(ix_hdr, 1, 16, fp) != 16) { fclose(fp); free(vs); return NULL; }
    if (memcmp(ix_hdr, VS_MAGIC_INDEX, 4) != 0) { fclose(fp); free(vs); return NULL; }

    uint32_t n_entries;
    memcpy(&n_entries, ix_hdr + 4, 4);
    vs->index_count = n_entries;

    /* Read index entries */
    vs->index = (VecIndexEntry *)calloc(n_entries, sizeof(VecIndexEntry));
    for (uint32_t i = 0; i < n_entries; i++) {
        uint8_t entry[24];
        if (fread(entry, 1, 24, fp) != 24) break;
        memcpy(&vs->index[i].time, entry, 8);
        memcpy(&vs->index[i].offset, entry + 8, 8);
        vs->index[i].frame_type = entry[16];
        vs->index[i].field_idx = entry[17];
    }

    return vs;
}

uint32_t vecstream_frame_count(const VecStream *vs) {
    return vs ? vs->index_count : 0;
}

double vecstream_frame_time(const VecStream *vs, uint32_t idx) {
    if (!vs || idx >= vs->index_count) return -1.0;
    return vs->index[idx].time;
}

uint8_t vecstream_frame_type(const VecStream *vs, uint32_t idx) {
    if (!vs || idx >= vs->index_count) return 255;
    return vs->index[idx].frame_type;
}

uint8_t vecstream_frame_field(const VecStream *vs, uint32_t idx) {
    if (!vs || idx >= vs->index_count) return 255;
    return vs->index[idx].field_idx;
}

/* Read compressed frame payload and decompress */
static uint8_t *vs_read_frame_payload(VecStream *vs, uint32_t idx,
                                       size_t *out_size) {
    if (!vs || idx >= vs->index_count) return NULL;

    fseeko(vs->fp, vs->index[idx].offset, SEEK_SET);

    uint8_t hdr[32];
    if (fread(hdr, 1, 32, vs->fp) != 32) return NULL;
    if (memcmp(hdr, VS_MAGIC_FRAME, 4) != 0) return NULL;

    uint64_t comp_size, uncomp_size;
    memcpy(&comp_size, hdr + 16, 8);
    memcpy(&uncomp_size, hdr + 24, 8);

    uint8_t *comp = (uint8_t *)malloc((size_t)comp_size);
    if (!comp) return NULL;
    if (fread(comp, 1, (size_t)comp_size, vs->fp) != (size_t)comp_size) {
        free(comp);
        return NULL;
    }

    size_t ds;
    uint8_t *data = vs_decompress(comp, (size_t)comp_size, (size_t)uncomp_size, &ds);
    free(comp);
    if (!data) return NULL;
    *out_size = ds;
    return data;
}

int vecstream_read_iframe(VecStream *vs, uint32_t idx,
                          VecPatch *patches, uint32_t *n_patches) {
    if (!vs || vs->index[idx].frame_type != VS_FRAME_I) return -1;

    size_t data_size;
    uint8_t *data = vs_read_frame_payload(vs, idx, &data_size);
    if (!data) return -1;

    uint8_t *p = data;
    uint32_t np;
    memcpy(&np, p, 4); p += 4;
    if (n_patches) *n_patches = np;

    for (uint32_t i = 0; i < np; i++) {
        VecPatch *vp = &patches[i];
        memset(vp, 0, sizeof(VecPatch));

        memcpy(&vp->origin_x, p + 0, 2);
        memcpy(&vp->origin_y, p + 2, 2);
        memcpy(&vp->origin_z, p + 4, 2);
        vp->size_x = p[6];
        vp->size_y = p[7];
        vp->size_z = p[8];
        vp->order = p[9];
        memcpy(&vp->n_coeffs, p + 10, 2);
        memcpy(&vp->max_err_f16, p + 12, 2);
        memcpy(&vp->rms_err_f16, p + 14, 2);
        p += 16;

        uint16_t nc = vp->n_coeffs;
        if (nc > VS_NCOEFFS) nc = VS_NCOEFFS;
        memcpy(vp->coeffs, p, nc * 4);
        p += vp->n_coeffs * 4;
    }

    free(data);
    return 0;
}

int vecstream_read_pframe(VecStream *vs, uint32_t idx, uint32_t *indices,
                          float *delta_coeffs, uint32_t *n_deltas,
                          uint16_t *coeffs_per_patch) {
    if (!vs || vs->index[idx].frame_type != VS_FRAME_P) return -1;

    size_t data_size;
    uint8_t *data = vs_read_frame_payload(vs, idx, &data_size);
    if (!data) return -1;

    uint8_t *p = data;
    uint32_t nd;
    uint16_t cpp;
    memcpy(&nd, p, 4); p += 4;
    memcpy(&cpp, p, 2); p += 2;
    p += 2; /* padding */

    if (n_deltas) *n_deltas = nd;
    if (coeffs_per_patch) *coeffs_per_patch = cpp;

    for (uint32_t i = 0; i < nd; i++) {
        uint32_t pidx;
        memcpy(&pidx, p, 4); p += 4;
        if (indices) indices[i] = pidx;
        if (delta_coeffs)
            memcpy(&delta_coeffs[(size_t)i * cpp], p, cpp * 4);
        p += (size_t)cpp * 4;
    }

    free(data);
    return 0;
}

int vecstream_read_kframe(VecStream *vs, uint32_t idx, float *voxels) {
    if (!vs || vs->index[idx].frame_type != VS_FRAME_K) return -1;

    size_t data_size;
    uint8_t *data = vs_read_frame_payload(vs, idx, &data_size);
    if (!data) return -1;

    uint64_t N3 = (uint64_t)vs->Nx * vs->Ny * vs->Nz;
    size_t expected = N3 * 4;
    if (data_size < expected) { free(data); return -1; }

    memcpy(voxels, data, expected);
    free(data);
    return 0;
}

/* ---- Patch fitting ---- */

void vecstream_fit_patch(const float *field, uint32_t Nx, uint32_t Ny, uint32_t Nz,
                         int ox, int oy, int oz, int bsx, int bsy, int bsz,
                         int order, VecPatch *patch) {
    memset(patch, 0, sizeof(VecPatch));
    patch->origin_x = (int16_t)ox;
    patch->origin_y = (int16_t)oy;
    patch->origin_z = (int16_t)oz;
    patch->size_x = (uint8_t)bsx;
    patch->size_y = (uint8_t)bsy;
    patch->size_z = (uint8_t)bsz;
    patch->order = (uint8_t)order;

    int O1 = order + 1;
    int nc = O1 * O1 * O1;
    patch->n_coeffs = (uint16_t)nc;

    /* 1D Vandermonde matrices */
    double Vx[32][4], Vy[32][4], Vz[32][4];
    for (int i = 0; i < bsx; i++) {
        double t = (bsx > 1) ? (double)i / (bsx - 1) : 0.0;
        double tk = 1.0;
        for (int a = 0; a < O1; a++) { Vx[i][a] = tk; tk *= t; }
    }
    for (int j = 0; j < bsy; j++) {
        double t = (bsy > 1) ? (double)j / (bsy - 1) : 0.0;
        double tk = 1.0;
        for (int b = 0; b < O1; b++) { Vy[j][b] = tk; tk *= t; }
    }
    for (int k = 0; k < bsz; k++) {
        double t = (bsz > 1) ? (double)k / (bsz - 1) : 0.0;
        double tk = 1.0;
        for (int c = 0; c < O1; c++) { Vz[k][c] = tk; tk *= t; }
    }

    /* Compute VtV^{-1} for each axis separately (4x4 system) */
    double VtVx[4][4], VtVy[4][4], VtVz[4][4];
    double VtVx_inv[4][4], VtVy_inv[4][4], VtVz_inv[4][4];

    /* Build VtV for x */
    for (int a = 0; a < O1; a++)
        for (int b = 0; b < O1; b++) {
            double s = 0;
            for (int i = 0; i < bsx; i++) s += Vx[i][a] * Vx[i][b];
            VtVx[a][b] = s;
        }

    /* Build VtV for y */
    for (int a = 0; a < O1; a++)
        for (int b = 0; b < O1; b++) {
            double s = 0;
            for (int j = 0; j < bsy; j++) s += Vy[j][a] * Vy[j][b];
            VtVy[a][b] = s;
        }

    /* Build VtV for z */
    for (int a = 0; a < O1; a++)
        for (int b = 0; b < O1; b++) {
            double s = 0;
            for (int k = 0; k < bsz; k++) s += Vz[k][a] * Vz[k][b];
            VtVz[a][b] = s;
        }

    /* Invert each VtV via Gauss-Jordan */
    #define VS_INVERT_4x4(VtV, inv) do { \
        double aug[4][8]; \
        for (int a_ = 0; a_ < O1; a_++) \
            for (int b_ = 0; b_ < O1; b_++) { \
                aug[a_][b_] = VtV[a_][b_]; \
                aug[a_][O1+b_] = (a_==b_) ? 1.0 : 0.0; \
            } \
        for (int col_ = 0; col_ < O1; col_++) { \
            int piv_ = col_; \
            for (int row_ = col_+1; row_ < O1; row_++) \
                if (fabs(aug[row_][col_]) > fabs(aug[piv_][col_])) piv_ = row_; \
            if (piv_ != col_) \
                for (int j_ = 0; j_ < 2*O1; j_++) { \
                    double tmp_ = aug[col_][j_]; aug[col_][j_] = aug[piv_][j_]; aug[piv_][j_] = tmp_; \
                } \
            double d_ = aug[col_][col_]; \
            if (fabs(d_) < 1e-15) continue; \
            for (int j_ = 0; j_ < 2*O1; j_++) aug[col_][j_] /= d_; \
            for (int row_ = 0; row_ < O1; row_++) { \
                if (row_ == col_) continue; \
                double f_ = aug[row_][col_]; \
                for (int j_ = 0; j_ < 2*O1; j_++) aug[row_][j_] -= f_ * aug[col_][j_]; \
            } \
        } \
        for (int a_ = 0; a_ < O1; a_++) \
            for (int b_ = 0; b_ < O1; b_++) \
                inv[a_][b_] = aug[a_][O1+b_]; \
    } while(0)

    VS_INVERT_4x4(VtVx, VtVx_inv);
    VS_INVERT_4x4(VtVy, VtVy_inv);
    VS_INVERT_4x4(VtVz, VtVz_inv);
    #undef VS_INVERT_4x4

    /* Projection matrices: P[i][a] = sum_m VtV_inv[a][m] * V[i][m] */
    double Px[32][4], Py[32][4], Pz[32][4];
    for (int i = 0; i < bsx; i++)
        for (int a = 0; a < O1; a++) {
            double s = 0;
            for (int m = 0; m < O1; m++) s += VtVx_inv[a][m] * Vx[i][m];
            Px[i][a] = s;
        }
    for (int j = 0; j < bsy; j++)
        for (int b = 0; b < O1; b++) {
            double s = 0;
            for (int m = 0; m < O1; m++) s += VtVy_inv[b][m] * Vy[j][m];
            Py[j][b] = s;
        }
    for (int k = 0; k < bsz; k++)
        for (int c = 0; c < O1; c++) {
            double s = 0;
            for (int m = 0; m < O1; m++) s += VtVz_inv[c][m] * Vz[k][m];
            Pz[k][c] = s;
        }

    /* Compute tensor-product coefficients */
    for (int a = 0; a < O1; a++)
    for (int b = 0; b < O1; b++)
    for (int c = 0; c < O1; c++) {
        double coeff = 0;
        for (int di = 0; di < bsx; di++) {
            int gi = ox + di;
            if (gi >= (int)Nx) continue;
            for (int dj = 0; dj < bsy; dj++) {
                int gj = oy + dj;
                if (gj >= (int)Ny) continue;
                for (int dk = 0; dk < bsz; dk++) {
                    int gk = oz + dk;
                    if (gk >= (int)Nz) continue;
                    long idx = (long)gi * (long)Ny * Nz + (long)gj * Nz + gk;
                    coeff += Px[di][a] * Py[dj][b] * Pz[dk][c] * (double)field[idx];
                }
            }
        }
        patch->coeffs[a * O1 * O1 + b * O1 + c] = (float)coeff;
    }

    /* Compute error statistics */
    float max_err = 0;
    double se = 0;
    int count = 0;
    for (int di = 0; di < bsx; di++) {
        int gi = ox + di;
        if (gi >= (int)Nx) continue;
        double tx = (bsx > 1) ? (double)di / (bsx - 1) : 0.0;
        for (int dj = 0; dj < bsy; dj++) {
            int gj = oy + dj;
            if (gj >= (int)Ny) continue;
            double ty = (bsy > 1) ? (double)dj / (bsy - 1) : 0.0;
            for (int dk = 0; dk < bsz; dk++) {
                int gk = oz + dk;
                if (gk >= (int)Nz) continue;
                double tz = (bsz > 1) ? (double)dk / (bsz - 1) : 0.0;

                long idx = (long)gi * (long)Ny * Nz + (long)gj * Nz + gk;
                double pred = 0;
                double txa[4] = {1, tx, tx*tx, tx*tx*tx};
                double tya[4] = {1, ty, ty*ty, ty*ty*ty};
                double tza[4] = {1, tz, tz*tz, tz*tz*tz};
                for (int a = 0; a < O1; a++)
                for (int b = 0; b < O1; b++)
                for (int c2 = 0; c2 < O1; c2++)
                    pred += patch->coeffs[a*O1*O1 + b*O1 + c2] * txa[a] * tya[b] * tza[c2];

                float err = fabsf(field[idx] - (float)pred);
                if (err > max_err) max_err = err;
                se += (double)err * err;
                count++;
            }
        }
    }
    float rms_err = (float)sqrt(se / (count > 0 ? count : 1));
    patch->max_err_f16 = vs_f32_to_f16(max_err);
    patch->rms_err_f16 = vs_f32_to_f16(rms_err);
}

/* ---- Patch evaluation ---- */

float vecstream_eval_patch(const VecPatch *p, int bs_x, int bs_y, int bs_z,
                           int di, int dj, int dk) {
    int O1 = p->order + 1;
    double tx = (bs_x > 1) ? (double)di / (bs_x - 1) : 0.0;
    double ty = (bs_y > 1) ? (double)dj / (bs_y - 1) : 0.0;
    double tz = (bs_z > 1) ? (double)dk / (bs_z - 1) : 0.0;

    double txa[4] = {1, tx, tx*tx, tx*tx*tx};
    double tya[4] = {1, ty, ty*ty, ty*ty*ty};
    double tza[4] = {1, tz, tz*tz, tz*tz*tz};

    double val = 0;
    for (int a = 0; a < O1; a++)
    for (int b = 0; b < O1; b++)
    for (int c = 0; c < O1; c++)
        val += p->coeffs[a * O1 * O1 + b * O1 + c] * txa[a] * tya[b] * tza[c];

    return (float)val;
}

/* ---- Patch-to-voxel reconstruction ---- */

void vecstream_patches_to_voxels(const VecPatch *patches, uint32_t n_patches,
                                 uint32_t Nx, uint32_t Ny, uint32_t Nz,
                                 float *output) {
    memset(output, 0, (size_t)Nx * Ny * Nz * sizeof(float));

    for (uint32_t pi = 0; pi < n_patches; pi++) {
        const VecPatch *p = &patches[pi];
        int ox = p->origin_x, oy = p->origin_y, oz = p->origin_z;
        int bsx = p->size_x, bsy = p->size_y, bsz = p->size_z;

        for (int di = 0; di < bsx; di++) {
            int gi = ox + di;
            if (gi < 0 || gi >= (int)Nx) continue;
            for (int dj = 0; dj < bsy; dj++) {
                int gj = oy + dj;
                if (gj < 0 || gj >= (int)Ny) continue;
                for (int dk = 0; dk < bsz; dk++) {
                    int gk = oz + dk;
                    if (gk < 0 || gk >= (int)Nz) continue;
                    long idx = (long)gi * (long)Ny * Nz + (long)gj * Nz + gk;
                    output[idx] = vecstream_eval_patch(p, bsx, bsy, bsz, di, dj, dk);
                }
            }
        }
    }
}

/* ---- Full reconstruction ---- */

int vecstream_reconstruct(VecStream *vs, uint32_t frame_idx, uint8_t field_idx,
                          float *output_voxels) {
    if (!vs || frame_idx >= vs->index_count) return -1;

    uint32_t tp = vs->total_patches;
    uint16_t cpp = VS_NCOEFFS;

    /* Allocate working patch array */
    VecPatch *patches = (VecPatch *)calloc(tp, sizeof(VecPatch));
    if (!patches) return -1;

    /* Find the most recent I-frame or K-frame for this field at or before frame_idx */
    int base_idx = -1;
    for (int i = (int)frame_idx; i >= 0; i--) {
        if (vs->index[i].field_idx != field_idx) continue;
        if (vs->index[i].frame_type == VS_FRAME_I ||
            vs->index[i].frame_type == VS_FRAME_K) {
            base_idx = i;
            break;
        }
    }

    if (base_idx < 0) { free(patches); return -1; }

    if (vs->index[base_idx].frame_type == VS_FRAME_K) {
        /* K-frame: decompress raw voxels */
        uint64_t N3 = (uint64_t)vs->Nx * vs->Ny * vs->Nz;
        int rc = vecstream_read_kframe(vs, (uint32_t)base_idx, output_voxels);
        /* Apply P-frames after the K-frame if needed, but K-frame already gives
           raw voxels. For P-frame application we need patches, so re-vectorize. */
        /* For simplicity: if the target is the K-frame itself, done. Otherwise
           we need to vectorize the K-frame and apply P-frames. */
        if ((uint32_t)base_idx == frame_idx) {
            free(patches);
            return rc;
        }
        /* Re-vectorize K-frame into patches for P-frame application */
        uint32_t BS = vs->block_size;
        uint32_t bpx = vs->blocks_per_axis_x;
        uint32_t bpy = vs->blocks_per_axis_y;
        uint32_t bpz = vs->blocks_per_axis_z;
        uint32_t pi = 0;
        for (uint32_t bi = 0; bi < bpx; bi++)
        for (uint32_t bj = 0; bj < bpy; bj++)
        for (uint32_t bk = 0; bk < bpz; bk++) {
            vecstream_fit_patch(output_voxels, vs->Nx, vs->Ny, vs->Nz,
                               bi * BS, bj * BS, bk * BS,
                               BS, BS, BS, VS_MAX_ORDER, &patches[pi]);
            pi++;
        }
    } else {
        /* I-frame */
        uint32_t np;
        if (vecstream_read_iframe(vs, (uint32_t)base_idx, patches, &np) != 0) {
            free(patches);
            return -1;
        }
    }

    /* Apply P-frames from base_idx+1 to frame_idx */
    for (uint32_t i = (uint32_t)base_idx + 1; i <= frame_idx; i++) {
        if (vs->index[i].field_idx != field_idx) continue;
        if (vs->index[i].frame_type != VS_FRAME_P) continue;

        uint32_t nd;
        uint16_t frame_cpp;
        /* First read to get count */
        uint32_t *d_indices = (uint32_t *)malloc(tp * sizeof(uint32_t));
        float *d_coeffs = (float *)malloc((size_t)tp * cpp * sizeof(float));
        if (!d_indices || !d_coeffs) {
            free(d_indices); free(d_coeffs); free(patches);
            return -1;
        }

        if (vecstream_read_pframe(vs, i, d_indices, d_coeffs, &nd, &frame_cpp) != 0) {
            free(d_indices); free(d_coeffs); free(patches);
            return -1;
        }

        /* Apply deltas */
        for (uint32_t d = 0; d < nd; d++) {
            uint32_t pidx = d_indices[d];
            if (pidx >= tp) continue;
            uint16_t nc = patches[pidx].n_coeffs;
            if (nc > frame_cpp) nc = frame_cpp;
            for (uint16_t c = 0; c < nc; c++)
                patches[pidx].coeffs[c] += d_coeffs[(size_t)d * frame_cpp + c];
        }

        free(d_indices);
        free(d_coeffs);
    }

    /* Reconstruct voxels from patches */
    vecstream_patches_to_voxels(patches, tp, vs->Nx, vs->Ny, vs->Nz, output_voxels);

    free(patches);
    return 0;
}

#endif /* VECSTREAM_IMPLEMENTATION */
