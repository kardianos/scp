/*  sfa.h — SCP Field Archive: single-header C library
 *
 *  #define SFA_IMPLEMENTATION in exactly ONE .c file before including.
 *  Requires: libzstd (-lzstd)
 *
 *  Writer:
 *    SFA *s = sfa_create("out.sfa", Nx, Ny, Nz, Lx, Ly, Lz, dt);
 *    sfa_add_column(s, "phi_x", SFA_F64, SFA_POSITION, 0);
 *    ...
 *    sfa_finalize_header(s);
 *    void *cols[] = {phi0, phi1, ...};
 *    sfa_write_frame(s, time, cols);
 *    sfa_close(s);
 *
 *  Reader:
 *    SFA *s = sfa_open("out.sfa");
 *    void *buf = (uint8_t*)malloc(s->frame_bytes);
 *    sfa_read_frame(s, frame_idx, buf);
 *    sfa_close(s);
 */

#ifndef SFA_H
#define SFA_H

#include <stdint.h>
#include <stdio.h>

/* ---- Constants ---- */

#define SFA_VERSION       3
#define SFA_MAGIC         "SFA\0"
#define SFA_JTOP_DEFAULT  1024
#define SFA_JMPF_DEFAULT  1024

/* SFAH.flags layout:
 *   bits 0-3:  compression codec
 *   bit 4:     AMR mode (multi-grid)
 *   bit 5:     streaming mode (indexes may be incomplete)
 *   bits 6-31: reserved
 */
#define SFA_FLAG_AMR       (1 << 4)
#define SFA_FLAG_STREAMING (1 << 5)

/* Compression codec flags (in SFAH.flags bits 0-3) */
#define SFA_CODEC_RAW     0   /* no compression */
#define SFA_CODEC_ZSTD    1   /* zstd only */
#define SFA_CODEC_BSS     2   /* byte stream split + zstd (default, best for floats) */
#define SFA_CODEC_F32BSS  3   /* f32 downcast + BSS + zstd (lossy, ~5-18x, viz) */
#define SFA_CODEC_F16BSS  4   /* f16 downcast + BSS + zstd (lossy, ~57-89x, preview) */
#define SFA_CODEC_BQ8     5   /* block-quant 8-bit + zstd (lossy, ~119-163x, thumbnail) */
#define SFA_CODEC_IDELTA  6   /* int64 spatial delta + BSS + zstd (lossless, +5-7%) */
#define SFA_CODEC_COLZSTD 7   /* per-column BSS+zstd: each column BSS-encoded then compressed.
                               * FRMD layout: [n_cols(4)] [comp_size(8)]×n_cols [data...]
                               * Enables parallel decompression and lazy column access.
                               * BSS byte reordering within each column gives ~2x compression
                               * vs whole-frame BSS+zstd, with 12-way parallel decompress. */

/* Chunk type codes (as uint32 for comparison) */
#define SFA_CHUNK_SFAH  0x48414653  /* "SFAH" little-endian */
#define SFA_CHUNK_CDEF  0x46454443  /* "CDEF" */
#define SFA_CHUNK_JTOP  0x504F544A  /* "JTOP" */
#define SFA_CHUNK_JMPF  0x46504D4A  /* "JMPF" */
#define SFA_CHUNK_FRMD  0x444D5246  /* "FRMD" */
#define SFA_CHUNK_KVMD  0x444D564B  /* "KVMD" */
#define SFA_CHUNK_GDEF  0x46454447  /* "GDEF" */
#define SFA_CHUNK_FRVD  0x44565246  /* "FRVD" — vector frame data */
#define SFA_CHUNK_FMSH  0x48534D46  /* "FMSH" — Voronoi mesh frame */
#define SFA_CHUNK_FCEL  0x4C454346  /* "FCEL" — cell-data frame    */
#define SFA_CHUNK_FCEP  0x50454346  /* "FCEP" — cell-data P-frame  */
#define SFA_CHUNK_FMTL  0x4C544D46  /* "FMTL" — per-cell temporal model frame */

/* Frame type codes (stored in SFA_L2Entry.frame_type) */
#define SFA_FRAME_VOXEL  0   /* standard voxel data (FRMD chunk) */
#define SFA_FRAME_VEC_I  1   /* vector I-frame: full polynomial patches (FRVD chunk) */
#define SFA_FRAME_VEC_P  2   /* vector P-frame: sparse delta patches (FRVD chunk) */
#define SFA_FRAME_VEC_K  3   /* vector K-frame: raw voxel keyframe for verification (FRMD chunk) */
#define SFA_FRAME_MESH   4   /* Voronoi mesh frame (FMSH chunk) — cell positions for cell-native data */
#define SFA_FRAME_CELL   5   /* per-cell field values (FCEL chunk) — uses most recently-seen FMSH */
#define SFA_FRAME_CELL_P 6   /* sparse cell delta (FCEP chunk) — applied on top of latest temporal model */
#define SFA_FRAME_TEMPORAL_MODEL 7 /* per-cell Fourier model (FMTL chunk) — supersedes earlier model in stream */

/* dtype codes */
enum {
    SFA_F16=0, SFA_F32=1, SFA_F64=2, SFA_F128=3,
    SFA_I8=4,  SFA_I16=5, SFA_I32=6, SFA_I64=7,
    SFA_U8=8,  SFA_U16=9, SFA_U32=10, SFA_U64=11
};

/* semantic codes */
enum {
    SFA_POSITION=0, SFA_ANGLE=1, SFA_VELOCITY=2, SFA_ACCELERATION=3,
    SFA_ENERGY=4, SFA_BINDING=5, SFA_TORSION=6, SFA_METRIC=7,
    SFA_MASK=8, SFA_CUSTOM=255
};

static const int sfa_dtype_size[] = {2,4,8,16, 1,2,4,8, 1,2,4,8};

/* ---- Common helpers (used by analysis tools — avoids copy-paste) ---- */

/* Convert f16 (IEEE 754 half) to f32. Handles zero/inf, no denormals. */
#ifndef SFA_F16_TO_F32_DEFINED
#define SFA_F16_TO_F32_DEFINED
static inline float sfa_f16_to_f32(uint16_t h) {
    uint16_t s = h & 0x8000; int e = (h >> 10) & 0x1F; uint16_t m = h & 0x3FF;
    if (e == 0) return 0;
    if (e == 31) return s ? -1e30f : 1e30f;
    union { uint32_t u; float f; } conv;
    conv.u = ((uint32_t)s << 16) | ((uint32_t)(e-15+127) << 23) | ((uint32_t)m << 13);
    return conv.f;
}
#endif

/* Read a single voxel from a frame buffer as float.
 * Declared here, defined after SFA struct (needs SFA* for column info).
 * Usage: sfa_read_voxel_f32(buf, sfa, col, idx) */
/* Full signature declared after SFA struct definition below */

/* ---- Structures ---- */

typedef struct {
    char name[12];
    uint8_t dtype, semantic, component, flags;
    double scale;
} SFA_Column;

typedef struct {
    uint64_t jmpf_offset;
    uint32_t first_frame;
    uint32_t frame_count;
} SFA_L1Entry;

typedef struct {
    double time;
    uint64_t offset;
    uint64_t compressed_size;
    uint32_t checksum;
    uint32_t frame_type;  /* SFA_FRAME_VOXEL=0, SFA_FRAME_VEC_I=1, VEC_P=2, VEC_K=3 */
} SFA_L2Entry;

/* Multi-grid: per-grid definition (stored in GDEF chunk) */
typedef struct {
    uint16_t grid_id;
    uint16_t parent_id;     /* 0xFFFF = root (no parent) */
    uint32_t Nx, Ny, Nz;
    double   Lx, Ly, Lz;   /* half-domain sizes */
    double   cx, cy, cz;   /* center in world coordinates */
    uint16_t ghost;         /* ghost zone width */
    uint16_t kvmd_set;      /* which KVMD set has this grid's params */
} SFA_GridDef;

#define SFA_MAX_GRIDS 16

typedef struct {
    FILE *fp;
    int mode;  /* 0=read, 1=write */

    /* Header */
    uint32_t version, flags;
    uint32_t Nx, Ny, Nz;
    double Lx, Ly, Lz, dt;
    uint32_t n_columns;
    uint32_t total_frames;
    uint64_t first_jtop_offset;
    uint64_t cdef_offset;
    uint32_t jtop_max, jmpf_max;
    uint64_t N_total;
    uint64_t frame_bytes;  /* uncompressed frame size */

    /* Schema */
    SFA_Column *columns;

    /* Write state */
    uint64_t cur_jtop_offset;
    uint32_t cur_jtop_entries;
    uint64_t cur_jmpf_offset;
    uint32_t cur_jmpf_entries;
    uint32_t cur_l1_index;

    /* Split output mode: frames go to individual .sfp files */
    int split_output;           /* 0=single file (default), 1=split .sfp */
    char split_base[512];       /* base path without extension */

    /* Read state / mmap */
    void *mmap_ptr;
    uint64_t mmap_size;

    /* KVMD metadata (optional, written between CDEF and JTOP) */
    #define SFA_MAX_KVMD_SETS 16
    int n_kvmd_sets;
    struct { uint8_t *payload; uint64_t payload_len; } kvmd_sets[SFA_MAX_KVMD_SETS];

    /* GDEF multi-grid (optional, written between CDEF and KVMD) */
    int n_grids;
    SFA_GridDef grids[SFA_MAX_GRIDS];

    /* Temporal model for FRVD P-frame reconstruction */
    float *vec_temporal_mean;   /* n_patches × n_coeffs */
    float *vec_temporal_amp;
    float *vec_temporal_phase;
    float vec_temporal_omega;
    int vec_temporal_valid;     /* nonzero if temporal model loaded */
    float *vec_prev_coeffs;    /* fallback: previous frame coefficients */
    uint32_t vec_n_patches;
    uint16_t vec_n_coeffs;
} SFA;

/* Implementation of sfa_read_voxel_f32 (declared above, needs SFA) */
static inline float sfa_read_voxel_f32(const void *buf, const SFA *sfa, int col, long idx) {
    uint64_t off = 0;
    for (int c = 0; c < col; c++) off += sfa->N_total * sfa_dtype_size[sfa->columns[c].dtype];
    int dt = sfa->columns[col].dtype;
    const uint8_t *src = (const uint8_t*)buf + off;
    if (dt == SFA_F16) return sfa_f16_to_f32(((const uint16_t*)src)[idx]);
    if (dt == SFA_F32) return ((const float*)src)[idx];
    if (dt == SFA_F64) return (float)((const double*)src)[idx];
    return 0;
}

/* KVMD: key-value metadata associated with a range of frames.
 * Multiple KVMD chunks support multi-resolution / multi-parameter runs.
 *
 * Binary layout per KVMD chunk (after standard 12-byte chunk header):
 *   version:     u8   (1)
 *   set_id:      u16  (0-based identifier)
 *   n_pairs:     u16  (number of key-value pairs)
 *   first_frame: u32  (0xFFFFFFFF = applies to all frames)
 *   frame_count: u32  (0xFFFFFFFF = all frames from first_frame onward)
 *   reserved:    u8[3]
 *   --- 16-byte sub-header above ---
 *   pairs:       [key\0 value\0]...  (n_pairs null-terminated string pairs)
 */

typedef struct {
    uint16_t set_id;
    uint32_t first_frame;   /* 0xFFFFFFFF = all */
    uint32_t frame_count;   /* 0xFFFFFFFF = all */
    int n_pairs;
    char keys[128][64];
    char values[128][256];
} SFA_KVMDSet;

/* ---- API ---- */

SFA *sfa_create(const char *path, uint32_t Nx, uint32_t Ny, uint32_t Nz,
                double Lx, double Ly, double Lz, double dt);
void sfa_add_column(SFA *s, const char *name, uint8_t dtype,
                    uint8_t semantic, uint8_t component);
void sfa_add_kvmd(SFA *s, uint16_t set_id, uint32_t first_frame,
                  uint32_t frame_count, const char **keys, const char **values, int n_pairs);
void sfa_add_grid(SFA *s, uint16_t grid_id, uint16_t parent_id,
                  uint32_t Nx, uint32_t Ny, uint32_t Nz,
                  double Lx, double Ly, double Lz,
                  double cx, double cy, double cz,
                  uint16_t ghost, uint16_t kvmd_set);
int  sfa_finalize_header(SFA *s);
int  sfa_write_frame(SFA *s, double time, void **column_data);
int  sfa_write_frame_ex(SFA *s, double time, void **column_data,
                        uint32_t grid_id, uint64_t frame_bytes_override);

/* ---- Cell-native frames (Voronoi foam) ----
 * sfa_write_mesh_frame writes an FMSH chunk carrying cell positions +
 *   volumes (and optionally faces/CSR) for the Voronoi cells.
 * sfa_write_cell_frame writes an FCEL chunk with per-cell column data.
 *   Subsequent FCEL frames are interpreted relative to the most recent
 *   FMSH frame in the file, allowing multi-resolution sequences.
 *
 * flags: bit0=has_faces, bit1=has_csr, bit2=pos_f64 (default f32). */
int sfa_write_mesh_frame(SFA *s, double time,
                         uint32_t N_cells, uint32_t N_faces,
                         double L,
                         const double *cell_pos,
                         const double *cell_vol,
                         const void *face_records,
                         const uint32_t *csr_off,
                         const uint32_t *csr_idx,
                         uint8_t flags);
int sfa_write_cell_frame(SFA *s, double time,
                         uint32_t N_cells, uint32_t n_columns,
                         uint8_t dtype, void * const *column_data);

typedef struct {
    uint32_t N_cells;
    uint32_t N_faces;
    double   L;
    double  *cell_pos;          /* 3 × N_cells */
    double  *cell_vol;          /* N_cells */
    void    *face_records;      /* N_faces × 40 bytes (NULL if absent) */
    uint32_t *csr_off;
    uint32_t *csr_idx;
    uint8_t  flags;
} SFA_Mesh;

int sfa_read_mesh_frame(SFA *s, uint32_t frame_idx, SFA_Mesh *out);
int sfa_read_cell_frame(SFA *s, uint32_t frame_idx,
                        uint32_t *out_N_cells, uint32_t *out_n_columns,
                        uint8_t *out_dtype, void **out_data);
void sfa_mesh_free(SFA_Mesh *m);

/* ---- Cell I-frame with temporal model + sparse P-frame ----
 * Pattern matches the FRVD I/P-frame format. The I-frame embeds a
 * Fourier model (mean + amp·cos(ω·t + phase)) per cell per column. The
 * P-frame stores only cells whose actual value departs from the model
 * prediction by more than the writer's threshold.
 *
 * mean, amp, phase: each n_columns × N_cells f32 arrays, column-major.
 * column_data: same layout as sfa_write_cell_frame (n_columns pointers
 *              to N_cells × dtype-size byte buffers).
 *
 * For the P-frame: cell_ids[n_changed] is a sorted list of cell indices
 * with significant residual; delta_values is n_changed × n_columns f32
 * residuals (column-major). The reader reconstructs the full state via
 *     state[c, k] = mean[c, k] + amp[c, k]·cos(ω·t + phase[c, k]) + delta_or_zero
 * where delta_or_zero is 0 unless c appears in cell_ids. */
int sfa_write_cell_iframe_temporal(SFA *s, double time,
                                    uint32_t N_cells, uint32_t n_columns,
                                    uint8_t dtype, void * const *column_data,
                                    float omega,
                                    const float *mean,
                                    const float *amp,
                                    const float *phase);
int sfa_write_cell_pframe(SFA *s, double time,
                          uint32_t N_cells, uint32_t n_columns,
                          uint32_t n_changed,
                          const uint32_t *cell_ids,
                          const float *delta_values);

/* Reader returns: n_changed, allocated cell_ids (4 × n_changed bytes),
 * allocated delta_values (4 × n_changed × n_columns bytes). Caller frees. */
int sfa_read_cell_pframe(SFA *s, uint32_t frame_idx,
                         uint32_t *out_N_cells, uint32_t *out_n_columns,
                         uint32_t *out_n_changed,
                         uint32_t **out_cell_ids,
                         float **out_delta_values);

/* Reader for the temporal-model fields embedded in an FCEL v2 I-frame.
 * Returns -1 if the frame doesn't carry a temporal model (flags bit2=0).
 * out_omega, out_mean/amp/phase are caller-owned (free after use). */
int sfa_read_cell_iframe_temporal(SFA *s, uint32_t frame_idx,
                                   float *out_omega,
                                   uint32_t *out_N_cells,
                                   uint32_t *out_n_columns,
                                   float **out_mean,
                                   float **out_amp,
                                   float **out_phase);

/* ---- Standalone temporal-model frame (FMTL) ----
 * Stores a per-cell Fourier model as its own frame, separate from any
 * particular FCEL I-frame. Once written, every subsequent FCEP frame
 * uses this model for prediction. A new FMTL supersedes the old.
 *
 * Use this when the model is nearly stationary and you don't want to
 * embed it in every I-frame — typical compression win is 5×–10× for
 * dense temporal captures with many I-frames.
 *
 * Layout: header (28 bytes) + omega (4) + mean (n_columns × N_cells × 4)
 *                          + amp (same)  + phase (same).
 * Total compressed size ≈ 50–100 MB for L=40 / 3.26M cells / 6 columns.
 */
int sfa_write_temporal_model_frame(SFA *s, double time,
                                    uint32_t N_cells, uint32_t n_columns,
                                    float omega,
                                    const float *mean,
                                    const float *amp,
                                    const float *phase);

int sfa_read_temporal_model_frame(SFA *s, uint32_t frame_idx,
                                   float *out_omega,
                                   uint32_t *out_N_cells,
                                   uint32_t *out_n_columns,
                                   float **out_mean,
                                   float **out_amp,
                                   float **out_phase);

void sfa_close(SFA *s);

/* Vector frame types for polynomial patch data.
 * These share the same file, index, and compression as voxel frames.
 * The frame_type field in the JMPF index distinguishes them. */

/* Write a vector I-frame: full set of polynomial patches.
 * Payload: [n_patches(u32)][block_size(u8)][n_coeffs(u16)][pad(u8)]
 *          [origins(i16×3 × n_patches)][coeffs(f32 × n_patches × n_coeffs)] */
int sfa_write_vec_iframe(SFA *s, double time,
                         uint32_t n_patches, uint8_t block_size, uint16_t n_coeffs,
                         const int16_t *origins,  /* n_patches × 3 */
                         const float *coeffs);    /* n_patches × n_coeffs */

/* Write a vector I-frame with temporal model for P-frame prediction.
 * flags byte bit 0 = temporal model present.
 * After coefficients: [omega(f32)][mean(f32×n)][amp(f32×n)][phase(f32×n)]
 * where n = n_patches × n_coeffs.
 * P-frame deltas are then: actual - (mean + amp*cos(omega*t + phase)). */
int sfa_write_vec_iframe_temporal(SFA *s, double time,
                         uint32_t n_patches, uint8_t block_size, uint16_t n_coeffs,
                         const int16_t *origins, const float *coeffs,
                         float omega, const float *temp_mean,
                         const float *temp_amp, const float *temp_phase);

/* Write a vector P-frame: sparse delta patches.
 * Payload: [n_deltas(u32)][n_coeffs(u16)][pad(u16)]
 *          [indices(u32 × n_deltas)][delta_coeffs(f32 × n_deltas × n_coeffs)] */
int sfa_write_vec_pframe(SFA *s, double time,
                         uint32_t n_deltas, uint16_t n_coeffs,
                         const uint32_t *indices,
                         const float *delta_coeffs);

/* Get the frame type for an indexed frame */
uint32_t sfa_frame_type(SFA *s, uint32_t frame_idx);

/* Streaming support:
 *
 * In streaming mode, the JTOP and JMPF are pre-written with zeroed entries.
 * Each sfa_write_frame appends FRMD data and does a seek+overwrite of the
 * 32-byte JMPF slot — no full rewrite needed.
 *
 * Receiver workflow:
 *   1. Receive the header block (SFAH+CDEF+JTOP+JMPF) — write to disk
 *   2. For each incoming FRMD: append to file, then seek to the pre-allocated
 *      JMPF slot and write the 32-byte entry (time, offset, size, checksum)
 *   3. On completion: seek to SFAH offset 68, write total_frames
 *
 * The JMPF slots are at known offsets:
 *   slot_offset = jmpf_chunk_offset + 12 + 8 + (frame_index % jmpf_max) * 32
 *
 * A separate fixup pass can index any .sfa with zeroed entries by scanning
 * FRMD chunks sequentially and filling in the JMPF slots.
 */

/* Compute the byte offset of JMPF entry for a given frame */
static inline uint64_t sfa_jmpf_slot_offset(SFA *s, uint32_t frame_in_table) {
    return s->cur_jmpf_offset + 12 + 8 + (uint64_t)frame_in_table * 32;
}

/* Split output: write frames as individual .sfp files.
 * Call sfa_enable_split() BEFORE sfa_finalize_header().
 * Each sfa_write_frame() then writes to <base>_NNNN.sfp.tmp,
 * renames to <base>_NNNN.sfp on completion (atomic).
 * The .sfa file contains only the header (no frame data).
 *
 * To reconstruct a single .sfa:
 *   sfa_combine("output.sfa", "output_*.sfp", "combined.sfa")
 * or manually:
 *   cat output.sfa output_0001.sfp ... > combined.sfa
 *   sfa_fixup_index combined.sfa
 *
 * SFP file format:
 *   "SFP1" (4 bytes) — magic
 *   uint32 frame_idx  — 0-based frame number
 *   double time       — simulation time
 *   FRMD chunk        — standard "FRMD" + size(8) + compressed data
 */
#define SFA_SFP_MAGIC "SFP1"
void sfa_enable_split(SFA *s);
int  sfa_combine(const char *header_path, const char *output_path,
                 const char **sfp_paths, int n_sfp);

/* Fixup: scan an incomplete .sfa file, fill in zeroed JMPF entries from FRMD chunks */
int sfa_fixup_index(const char *path);

/* Count valid (non-truncated, non-missing) frames given the actual file size.
 * Returns the number of consecutive valid frames from the start. */
int sfa_count_valid_frames(SFA *s, uint64_t file_bytes);

SFA *sfa_open(const char *path);
int  sfa_read_frame(SFA *s, uint32_t frame_idx, void *buf);
double sfa_frame_time(SFA *s, uint32_t frame_idx);
int  sfa_read_kvmd(SFA *s, SFA_KVMDSet *sets, int max_sets);
int  sfa_read_grids(SFA *s, SFA_GridDef *grids, int max_grids);

const char *sfa_dtype_name(uint8_t dtype);

#endif /* SFA_H */


/* ================================================================
   IMPLEMENTATION
   ================================================================ */

#ifdef SFA_IMPLEMENTATION

/* Suppress -Wunused-result for fread/fwrite in the implementation.
 * The SFA format uses fread for sequential binary parsing where partial
 * reads indicate a truncated file, which is handled at a higher level
 * by sfa_count_valid_frames() rather than per-call error checking. */
#if defined(__GNUC__) || defined(__clang__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-result"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <zstd.h>

/* ---- Helpers ---- */

/* ---- Byte Stream Split: transpose bytes for better compression ----
 *
 * Input (N values of S bytes each):
 *   v0[b0 b1 ... bS-1], v1[b0 b1 ... bS-1], ...
 *
 * Output (grouped by byte position):
 *   all b0 bytes, all b1 bytes, ..., all bS-1 bytes
 *
 * This clusters exponent bytes (low entropy), sign bits, and mantissa
 * bytes separately. Exponents compress ~10x, mantissa bytes 1-2x,
 * overall 2-5x improvement over raw zstd on float64 data.
 */
static void sfa_bss_encode(const uint8_t *src, uint8_t *dst,
                           size_t n_values, int elem_size) {
    for (int b = 0; b < elem_size; b++)
        for (size_t i = 0; i < n_values; i++)
            dst[b * n_values + i] = src[i * elem_size + b];
}

static void sfa_bss_decode(const uint8_t *src, uint8_t *dst,
                           size_t n_values, int elem_size) {
    for (int b = 0; b < elem_size; b++)
        for (size_t i = 0; i < n_values; i++)
            dst[i * elem_size + b] = src[b * n_values + i];
}

/* Apply BSS to an entire frame (multiple columns, possibly different dtypes) */
static void sfa_bss_encode_frame(SFA *s, const uint8_t *raw, uint8_t *bss) {
    uint64_t off = 0;
    for (uint32_t c = 0; c < s->n_columns; c++) {
        int es = sfa_dtype_size[s->columns[c].dtype];
        uint64_t col_bytes = s->N_total * es;
        sfa_bss_encode(raw + off, bss + off, s->N_total, es);
        off += col_bytes;
    }
}

static void sfa_bss_decode_frame(SFA *s, const uint8_t *bss, uint8_t *raw) {
    uint64_t off = 0;
    for (uint32_t c = 0; c < s->n_columns; c++) {
        int es = sfa_dtype_size[s->columns[c].dtype];
        uint64_t col_bytes = s->N_total * es;
        sfa_bss_decode(bss + off, raw + off, s->N_total, es);
        off += col_bytes;
    }
}

static uint32_t sfa_crc32(const void *data, size_t len) {
    /* Simple CRC32 (Castagnoli) */
    uint32_t crc = 0xFFFFFFFF;
    const uint8_t *p = (const uint8_t *)data;
    for (size_t i = 0; i < len; i++) {
        crc ^= p[i];
        for (int j = 0; j < 8; j++)
            crc = (crc >> 1) ^ (0x82F63B78 & -(crc & 1));
    }
    return ~crc;
}

static void sfa_write_chunk_header(FILE *fp, const char type[4], uint64_t size) {
    fwrite(type, 1, 4, fp);
    fwrite(&size, 8, 1, fp);
}

static void sfa_compute_frame_bytes(SFA *s) {
    s->N_total = (uint64_t)s->Nx * s->Ny * s->Nz;
    s->frame_bytes = 0;
    for (uint32_t c = 0; c < s->n_columns; c++)
        s->frame_bytes += s->N_total * sfa_dtype_size[s->columns[c].dtype];
}

/* ---- Writer ---- */

SFA *sfa_create(const char *path, uint32_t Nx, uint32_t Ny, uint32_t Nz,
                double Lx, double Ly, double Lz, double dt) {
    SFA *s = (SFA*)calloc(1, sizeof(SFA));
    s->fp = fopen(path, "wb");
    if (!s->fp) { free(s); return NULL; }
    s->mode = 1;
    s->version = SFA_VERSION;
    s->flags = SFA_CODEC_BSS | SFA_FLAG_STREAMING;  /* streaming by default */
    s->Nx = Nx; s->Ny = Ny; s->Nz = Nz;
    s->Lx = Lx; s->Ly = Ly; s->Lz = Lz;
    s->dt = dt;
    s->jtop_max = SFA_JTOP_DEFAULT;
    s->jmpf_max = SFA_JMPF_DEFAULT;
    s->columns = (SFA_Column*)calloc(64, sizeof(SFA_Column));
    s->n_columns = 0;
    return s;
}

void sfa_add_column(SFA *s, const char *name, uint8_t dtype,
                    uint8_t semantic, uint8_t component) {
    SFA_Column *c = &s->columns[s->n_columns++];
    strncpy(c->name, name, 11);
    c->dtype = dtype;
    c->semantic = semantic;
    c->component = component;
    c->flags = 0;
    c->scale = 1.0;
}

void sfa_add_grid(SFA *s, uint16_t grid_id, uint16_t parent_id,
                  uint32_t Nx, uint32_t Ny, uint32_t Nz,
                  double Lx, double Ly, double Lz,
                  double cx, double cy, double cz,
                  uint16_t ghost, uint16_t kvmd_set) {
    if (s->n_grids >= SFA_MAX_GRIDS) return;
    SFA_GridDef *g = &s->grids[s->n_grids++];
    g->grid_id = grid_id; g->parent_id = parent_id;
    g->Nx = Nx; g->Ny = Ny; g->Nz = Nz;
    g->Lx = Lx; g->Ly = Ly; g->Lz = Lz;
    g->cx = cx; g->cy = cy; g->cz = cz;
    g->ghost = ghost; g->kvmd_set = kvmd_set;
}

void sfa_add_kvmd(SFA *s, uint16_t set_id, uint32_t first_frame,
                  uint32_t frame_count, const char **keys, const char **values, int n_pairs) {
    if (s->n_kvmd_sets >= SFA_MAX_KVMD_SETS) return;
    /* Compute payload size: 16-byte sub-header + pairs */
    uint64_t pair_bytes = 0;
    for (int i = 0; i < n_pairs; i++)
        pair_bytes += strlen(keys[i]) + 1 + strlen(values[i]) + 1;
    uint64_t payload_len = 16 + pair_bytes;

    uint8_t *buf = (uint8_t*)malloc(payload_len);
    uint64_t off = 0;
    /* Sub-header: version(1) set_id(2) n_pairs(2) first_frame(4) frame_count(4) reserved(3) */
    buf[off++] = 1;  /* version */
    memcpy(buf + off, &set_id, 2); off += 2;
    uint16_t np = (uint16_t)n_pairs;
    memcpy(buf + off, &np, 2); off += 2;
    memcpy(buf + off, &first_frame, 4); off += 4;
    memcpy(buf + off, &frame_count, 4); off += 4;
    buf[off++] = 0; buf[off++] = 0; buf[off++] = 0;  /* reserved */
    /* Pairs */
    for (int i = 0; i < n_pairs; i++) {
        uint64_t kl = strlen(keys[i]) + 1, vl = strlen(values[i]) + 1;
        memcpy(buf + off, keys[i], kl); off += kl;
        memcpy(buf + off, values[i], vl); off += vl;
    }
    int idx = s->n_kvmd_sets++;
    s->kvmd_sets[idx].payload = buf;
    s->kvmd_sets[idx].payload_len = payload_len;
}

int sfa_finalize_header(SFA *s) {
    sfa_compute_frame_bytes(s);
    FILE *fp = s->fp;

    /* ---- SFAH ---- */
    uint64_t sfah_size = 100;
    sfa_write_chunk_header(fp, "SFAH", sfah_size);
    fwrite(&s->version, 4, 1, fp);
    fwrite(&s->flags, 4, 1, fp);
    fwrite(&s->Nx, 4, 1, fp);
    fwrite(&s->Ny, 4, 1, fp);
    fwrite(&s->Nz, 4, 1, fp);
    fwrite(&s->Lx, 8, 1, fp);
    fwrite(&s->Ly, 8, 1, fp);
    fwrite(&s->Lz, 8, 1, fp);
    fwrite(&s->dt, 8, 1, fp);
    fwrite(&s->n_columns, 4, 1, fp);
    s->total_frames = 0;
    fwrite(&s->total_frames, 4, 1, fp);
    /* first_jtop_offset: will patch */
    uint64_t zero64 = 0;
    long jtop_off_pos = ftell(fp);
    fwrite(&zero64, 8, 1, fp);     /* first_jtop_offset placeholder */
    s->cdef_offset = ftell(fp) + 8 + 4 + 4;  /* after remaining header fields */
    /* cdef_offset placeholder */
    long cdef_off_pos = ftell(fp);
    fwrite(&zero64, 8, 1, fp);
    fwrite(&s->jtop_max, 4, 1, fp);
    fwrite(&s->jmpf_max, 4, 1, fp);
    uint32_t zero32 = 0;
    fwrite(&zero32, 4, 1, fp);     /* reserved */

    /* ---- CDEF ---- */
    s->cdef_offset = (uint64_t)ftell(fp);
    uint64_t cdef_size = 12 + (uint64_t)s->n_columns * 24;
    sfa_write_chunk_header(fp, "CDEF", cdef_size);
    for (uint32_t c = 0; c < s->n_columns; c++) {
        fwrite(s->columns[c].name, 1, 12, fp);
        fwrite(&s->columns[c].dtype, 1, 1, fp);
        fwrite(&s->columns[c].semantic, 1, 1, fp);
        fwrite(&s->columns[c].component, 1, 1, fp);
        fwrite(&s->columns[c].flags, 1, 1, fp);
        fwrite(&s->columns[c].scale, 8, 1, fp);
    }

    /* ---- GDEF (optional, between CDEF and KVMD/JTOP) ---- */
    if (s->n_grids > 0) {
        /* Sub-header: version(1) + n_grids(2) + reserved(5) = 8 bytes
         * Per grid: grid_id(2) + parent_id(2) + Nx/Ny/Nz(12) + Lx/Ly/Lz(24) +
         *           cx/cy/cz(24) + ghost(2) + kvmd_set(2) + reserved(4) = 72 bytes */
        uint64_t gdef_size = 8 + (uint64_t)s->n_grids * 72;
        sfa_write_chunk_header(fp, "GDEF", gdef_size);
        uint8_t gdef_ver = 1;
        uint16_t ng = (uint16_t)s->n_grids;
        uint8_t gdef_reserved[5] = {0};
        fwrite(&gdef_ver, 1, 1, fp);
        fwrite(&ng, 2, 1, fp);
        fwrite(gdef_reserved, 1, 5, fp);
        for (int gi = 0; gi < s->n_grids; gi++) {
            SFA_GridDef *gd = &s->grids[gi];
            uint32_t gd_reserved = 0;
            fwrite(&gd->grid_id, 2, 1, fp);
            fwrite(&gd->parent_id, 2, 1, fp);
            fwrite(&gd->Nx, 4, 1, fp); fwrite(&gd->Ny, 4, 1, fp); fwrite(&gd->Nz, 4, 1, fp);
            fwrite(&gd->Lx, 8, 1, fp); fwrite(&gd->Ly, 8, 1, fp); fwrite(&gd->Lz, 8, 1, fp);
            fwrite(&gd->cx, 8, 1, fp); fwrite(&gd->cy, 8, 1, fp); fwrite(&gd->cz, 8, 1, fp);
            fwrite(&gd->ghost, 2, 1, fp);
            fwrite(&gd->kvmd_set, 2, 1, fp);
            fwrite(&gd_reserved, 4, 1, fp);
        }
    }

    /* ---- KVMD (optional, between GDEF and JTOP) ---- */
    for (int ki = 0; ki < s->n_kvmd_sets; ki++) {
        sfa_write_chunk_header(fp, "KVMD", s->kvmd_sets[ki].payload_len);
        fwrite(s->kvmd_sets[ki].payload, 1, s->kvmd_sets[ki].payload_len, fp);
        free(s->kvmd_sets[ki].payload);
        s->kvmd_sets[ki].payload = NULL;
    }

    /* ---- JTOP ---- */
    s->first_jtop_offset = (uint64_t)ftell(fp);
    s->cur_jtop_offset = s->first_jtop_offset;
    uint64_t jtop_size = 12 + 16 + (uint64_t)s->jtop_max * 16;
    sfa_write_chunk_header(fp, "JTOP", jtop_size);
    /* JTOP header: max, current, next_offset */
    fwrite(&s->jtop_max, 4, 1, fp);
    s->cur_jtop_entries = 0;
    fwrite(&s->cur_jtop_entries, 4, 1, fp);
    fwrite(&zero64, 8, 1, fp);  /* next_jtop_offset = 0 */
    /* JTOP entries: zero-fill */
    for (uint32_t i = 0; i < s->jtop_max; i++) {
        fwrite(&zero64, 8, 1, fp);  /* jmpf_offset */
        fwrite(&zero32, 4, 1, fp);  /* first_frame */
        fwrite(&zero32, 4, 1, fp);  /* frame_count */
    }

    /* ---- First JMPF ---- */
    s->cur_jmpf_offset = (uint64_t)ftell(fp);
    uint64_t jmpf_size = 12 + 8 + (uint64_t)s->jmpf_max * 32;
    sfa_write_chunk_header(fp, "JMPF", jmpf_size);
    fwrite(&s->jmpf_max, 4, 1, fp);
    s->cur_jmpf_entries = 0;
    fwrite(&s->cur_jmpf_entries, 4, 1, fp);
    /* JMPF entries: zero-fill */
    for (uint32_t i = 0; i < s->jmpf_max; i++) {
        double zt = 0; fwrite(&zt, 8, 1, fp);
        fwrite(&zero64, 8, 1, fp);
        fwrite(&zero64, 8, 1, fp);
        fwrite(&zero32, 4, 1, fp);
        fwrite(&zero32, 4, 1, fp);
    }

    /* Register first JMPF in JTOP */
    s->cur_l1_index = 0;
    /* Patch JTOP entry 0 */
    long jtop_entry0 = (long)s->cur_jtop_offset + 12 + 16;  /* after chunk header + jtop header */
    fseek(fp, jtop_entry0, SEEK_SET);
    fwrite(&s->cur_jmpf_offset, 8, 1, fp);
    zero32 = 0; fwrite(&zero32, 4, 1, fp);  /* first_frame = 0 */
    fwrite(&zero32, 4, 1, fp);  /* frame_count = 0 (updated later) */
    /* Update JTOP current_entries = 1 */
    s->cur_jtop_entries = 1;
    fseek(fp, (long)s->cur_jtop_offset + 12 + 4, SEEK_SET);
    fwrite(&s->cur_jtop_entries, 4, 1, fp);

    /* Patch SFAH pointers */
    fseek(fp, jtop_off_pos, SEEK_SET);
    fwrite(&s->first_jtop_offset, 8, 1, fp);
    fseek(fp, cdef_off_pos, SEEK_SET);
    fwrite(&s->cdef_offset, 8, 1, fp);

    /* Seek to end for frame writing */
    fseek(fp, 0, SEEK_END);
    return 0;
}

int sfa_write_frame(SFA *s, double time, void **column_data) {
    return sfa_write_frame_ex(s, time, column_data, 0, s->frame_bytes);
}

int sfa_write_frame_ex(SFA *s, double time, void **column_data,
                       uint32_t grid_id, uint64_t frame_bytes_override) {
    /* Temporarily swap N_total and frame_bytes for variable-size grids */
    uint64_t saved_N_total = s->N_total;
    uint64_t saved_frame_bytes = s->frame_bytes;
    if (frame_bytes_override != s->frame_bytes) {
        /* Compute N_total for this grid from frame_bytes and column dtypes */
        uint64_t bytes_per_point = 0;
        for (uint32_t c = 0; c < s->n_columns; c++)
            bytes_per_point += sfa_dtype_size[s->columns[c].dtype];
        s->N_total = frame_bytes_override / bytes_per_point;
        s->frame_bytes = frame_bytes_override;
    }

    /* Assemble uncompressed frame */
    uint8_t *raw = (uint8_t*)malloc(s->frame_bytes);
    uint64_t off = 0;
    for (uint32_t c = 0; c < s->n_columns; c++) {
        uint64_t col_bytes = s->N_total * sfa_dtype_size[s->columns[c].dtype];
        memcpy(raw + off, column_data[c], col_bytes);
        off += col_bytes;
    }

    /* CRC on original data (before any transform) */
    uint32_t checksum = sfa_crc32(raw, s->frame_bytes);

    int codec = s->flags & 0xF;
    void *comp;
    size_t comp_size;

    if (codec == SFA_CODEC_RAW) {
        /* No compression — write raw data directly */
        comp = raw;
        comp_size = s->frame_bytes;
        raw = NULL;  /* don't free — comp points to it */
    } else if (codec == SFA_CODEC_COLZSTD) {
        /* Per-column BSS+zstd with OpenMP parallel compression.
         * Each column is independently BSS-encoded then zstd-compressed
         * in parallel across available threads. Assembly is sequential.
         * Layout: [n_cols(4)] [comp_size(8)]×n_cols [col0_data...] [col1_data...] ... */
        uint32_t nc = s->n_columns;
        uint64_t *col_comp_sizes = (uint64_t*)calloc(nc, sizeof(uint64_t));
        uint8_t **col_comp_data = (uint8_t**)calloc(nc, sizeof(uint8_t*));
        int col_error = 0;

        /* Compute column offsets and sizes (needed for parallel indexing) */
        uint64_t *col_offsets = (uint64_t*)calloc(nc, sizeof(uint64_t));
        uint64_t *col_bytesv = (uint64_t*)calloc(nc, sizeof(uint64_t));
        {
            uint64_t off2 = 0;
            for (uint32_t c = 0; c < nc; c++) {
                col_offsets[c] = off2;
                col_bytesv[c] = s->N_total * sfa_dtype_size[s->columns[c].dtype];
                off2 += col_bytesv[c];
            }
        }

        /* Parallel BSS-encode + zstd compress each column */
        #ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic)
        #endif
        for (uint32_t c = 0; c < nc; c++) {
            if (col_error) continue;  /* skip if a previous column failed */
            int es = sfa_dtype_size[s->columns[c].dtype];
            uint64_t col_bytes = col_bytesv[c];

            /* BSS-encode: [v0b0,v0b1,...] → [b0v0,b0v1,...,b1v0,...] */
            uint8_t *bss_col = (uint8_t*)malloc(col_bytes);
            uint8_t *src_col = raw + col_offsets[c];
            for (int b = 0; b < es; b++) {
                for (uint64_t v = 0; v < s->N_total; v++) {
                    bss_col[b * s->N_total + v] = src_col[v * es + b];
                }
            }

            /* Compress */
            size_t bound = ZSTD_compressBound(col_bytes);
            col_comp_data[c] = (uint8_t*)malloc(bound);
            col_comp_sizes[c] = ZSTD_compress(col_comp_data[c], bound,
                                               bss_col, col_bytes, 3);
            free(bss_col);
            if (ZSTD_isError(col_comp_sizes[c])) {
                fprintf(stderr, "SFA: colzstd error col %u: %s\n", c,
                        ZSTD_getErrorName(col_comp_sizes[c]));
                col_error = 1;
            }
        }

        free(col_offsets); free(col_bytesv);

        if (col_error) {
            for (uint32_t c = 0; c < nc; c++) free(col_comp_data[c]);
            free(col_comp_sizes); free(col_comp_data); free(raw);
            return -1;
        }

        /* Sequential assembly: [n_cols][sizes...][data...] */
        uint64_t total_comp = 4 + nc * 8;
        for (uint32_t c = 0; c < nc; c++) total_comp += col_comp_sizes[c];

        comp = malloc(total_comp);
        uint8_t *p = (uint8_t*)comp;
        memcpy(p, &nc, 4); p += 4;
        memcpy(p, col_comp_sizes, nc * 8); p += nc * 8;
        for (uint32_t c = 0; c < nc; c++) {
            memcpy(p, col_comp_data[c], col_comp_sizes[c]);
            p += col_comp_sizes[c];
            free(col_comp_data[c]);
        }
        comp_size = total_comp;
        free(col_comp_sizes); free(col_comp_data); free(raw);
    } else {
        /* Apply BSS transform if enabled */
        uint8_t *to_compress = raw;
        uint8_t *bss_buf = NULL;
        if (codec == SFA_CODEC_BSS) {
            bss_buf = (uint8_t*)malloc(s->frame_bytes);
            sfa_bss_encode_frame(s, raw, bss_buf);
            to_compress = bss_buf;
        }

        /* Compress */
        size_t bound = ZSTD_compressBound(s->frame_bytes);
        comp = (void*)malloc(bound);
        comp_size = ZSTD_compress(comp, bound, to_compress, s->frame_bytes, 3);
        free(raw);
        if (bss_buf) free(bss_buf);
    }

    if (comp_size > 0 && ZSTD_isError(comp_size)) {
        fprintf(stderr, "SFA: zstd error: %s\n", ZSTD_getErrorName(comp_size));
        free(comp);
        return -1;
    }

    /* Need new JMPF? */
    if (s->cur_jmpf_entries >= s->jmpf_max) {
        /* Write new JMPF */
        fseek(s->fp, 0, SEEK_END);
        uint64_t new_jmpf = (uint64_t)ftell(s->fp);
        uint64_t jmpf_size = 12 + 8 + (uint64_t)s->jmpf_max * 32;
        sfa_write_chunk_header(s->fp, "JMPF", jmpf_size);
        fwrite(&s->jmpf_max, 4, 1, s->fp);
        uint32_t zero32 = 0;
        fwrite(&zero32, 4, 1, s->fp);  /* current_entries = 0 */
        /* Zero-fill entries */
        uint64_t zero64 = 0; double zt = 0;
        for (uint32_t i = 0; i < s->jmpf_max; i++) {
            fwrite(&zt, 8, 1, s->fp);
            fwrite(&zero64, 8, 1, s->fp);
            fwrite(&zero64, 8, 1, s->fp);
            fwrite(&zero32, 4, 1, s->fp);
            fwrite(&zero32, 4, 1, s->fp);
        }

        /* Register in JTOP */
        s->cur_l1_index++;
        /* TODO: handle JTOP overflow (chain new JTOP) — for >1M frames */
        long jtop_entry = (long)s->cur_jtop_offset + 12 + 16 + s->cur_l1_index * 16;
        fseek(s->fp, jtop_entry, SEEK_SET);
        fwrite(&new_jmpf, 8, 1, s->fp);
        uint32_t ff = s->total_frames;
        fwrite(&ff, 4, 1, s->fp);
        fwrite(&zero32, 4, 1, s->fp);
        /* Update JTOP current_entries */
        s->cur_jtop_entries = s->cur_l1_index + 1;
        fseek(s->fp, (long)s->cur_jtop_offset + 12 + 4, SEEK_SET);
        fwrite(&s->cur_jtop_entries, 4, 1, s->fp);

        s->cur_jmpf_offset = new_jmpf;
        s->cur_jmpf_entries = 0;
        fseek(s->fp, 0, SEEK_END);
    }

    if (s->split_output) {
        /* Split mode: write FRMD to a separate .sfp file with atomic rename */
        char tmp_path[600], final_path[600];
        snprintf(tmp_path,   600, "%s_%04d.sfp.tmp", s->split_base, s->total_frames + 1);
        snprintf(final_path, 600, "%s_%04d.sfp",     s->split_base, s->total_frames + 1);

        FILE *sfp = fopen(tmp_path, "wb");
        if (!sfp) {
            fprintf(stderr, "SFA: cannot open SFP file %s\n", tmp_path);
            free(comp);
            return -1;
        }

        /* SFP header: magic + frame_idx + time */
        fwrite(SFA_SFP_MAGIC, 1, 4, sfp);
        uint32_t fidx = s->total_frames;
        fwrite(&fidx, 4, 1, sfp);
        fwrite(&time, 8, 1, sfp);

        /* Standard FRMD chunk */
        uint64_t frmd_size = 12 + comp_size;
        sfa_write_chunk_header(sfp, "FRMD", frmd_size);
        fwrite(comp, 1, comp_size, sfp);
        fclose(sfp);

        /* Atomic rename: .sfp.tmp -> .sfp */
        rename(tmp_path, final_path);

        free(comp);
        s->total_frames++;
    } else {
        /* Standard mode: append FRMD to the .sfa file and update index */
        fseek(s->fp, 0, SEEK_END);
        uint64_t frmd_offset = (uint64_t)ftell(s->fp);
        uint64_t frmd_size = 12 + comp_size;
        sfa_write_chunk_header(s->fp, "FRMD", frmd_size);
        fwrite(comp, 1, comp_size, s->fp);
        free(comp);

        /* Update JMPF entry */
        long jmpf_entry = (long)s->cur_jmpf_offset + 12 + 8 + s->cur_jmpf_entries * 32;
        fseek(s->fp, jmpf_entry, SEEK_SET);
        fwrite(&time, 8, 1, s->fp);
        fwrite(&frmd_offset, 8, 1, s->fp);
        fwrite(&comp_size, 8, 1, s->fp);
        fwrite(&checksum, 4, 1, s->fp);
        fwrite(&grid_id, 4, 1, s->fp);

        /* Update JMPF current_entries */
        s->cur_jmpf_entries++;
        fseek(s->fp, (long)s->cur_jmpf_offset + 12 + 4, SEEK_SET);
        fwrite(&s->cur_jmpf_entries, 4, 1, s->fp);

        /* Update JTOP frame_count for current L2 */
        long jtop_fc = (long)s->cur_jtop_offset + 12 + 16 + s->cur_l1_index * 16 + 12;
        fseek(s->fp, jtop_fc, SEEK_SET);
        fwrite(&s->cur_jmpf_entries, 4, 1, s->fp);

        s->total_frames++;
        fseek(s->fp, 0, SEEK_END);
    }

    /* Restore original N_total and frame_bytes */
    s->N_total = saved_N_total;
    s->frame_bytes = saved_frame_bytes;
    return 0;
}

/* ---- Vector frame writers ---- */

/* Internal: write a FRVD chunk with the given payload and frame_type */
/* Generalized chunk writer: zstd-compress payload, write a chunk of the
 * given fourcc with [time(8)|compressed_data] body, and append a JMPF entry
 * with the requested frame_type. Used by FMSH and FCEL writers. */
static int sfa_write_typed_chunk(SFA *s, double time,
                                  uint32_t chunk_type, uint32_t frame_type,
                                  const void *payload, uint64_t payload_size) {
    FILE *fp = s->fp;
    if (!fp || s->mode != 1) return -1;

    size_t bound = ZSTD_compressBound(payload_size);
    void *comp = malloc(bound);
    size_t comp_size = ZSTD_compress(comp, bound, payload, payload_size, 3);
    if (ZSTD_isError(comp_size)) { free(comp); return -1; }

    fseek(fp, 0, SEEK_END);
    uint64_t frame_offset = ftell(fp);
    uint64_t chunk_size = 8 + comp_size;
    fwrite(&chunk_type, 4, 1, fp);
    fwrite(&chunk_size, 8, 1, fp);
    fwrite(&time, 8, 1, fp);
    fwrite(comp, 1, comp_size, fp);
    free(comp);

    uint32_t checksum = sfa_crc32(payload, payload_size);
    if (s->cur_jmpf_offset) {
        /* JMPF overflow handling */
        if (s->cur_jmpf_entries >= s->jmpf_max) {
            fseek(fp, 0, SEEK_END);
            uint64_t new_jmpf = (uint64_t)ftell(fp);
            uint64_t jmpf_size = 12 + 8 + (uint64_t)s->jmpf_max * 32;
            sfa_write_chunk_header(fp, "JMPF", jmpf_size);
            fwrite(&s->jmpf_max, 4, 1, fp);
            uint32_t zero32 = 0;
            fwrite(&zero32, 4, 1, fp);
            uint64_t zero64 = 0; double zt = 0;
            for (uint32_t i = 0; i < s->jmpf_max; i++) {
                fwrite(&zt, 8, 1, fp); fwrite(&zero64, 8, 1, fp);
                fwrite(&zero64, 8, 1, fp); fwrite(&zero32, 4, 1, fp);
                fwrite(&zero32, 4, 1, fp);
            }
            s->cur_l1_index++;
            long jtop_entry = (long)s->cur_jtop_offset + 12 + 16 + s->cur_l1_index * 16;
            fseek(fp, jtop_entry, SEEK_SET);
            fwrite(&new_jmpf, 8, 1, fp);
            uint32_t first_frame = s->total_frames;
            fwrite(&first_frame, 4, 1, fp);
            fwrite(&zero32, 4, 1, fp);
            fseek(fp, (long)s->cur_jtop_offset + 12 + 4, SEEK_SET);
            uint32_t new_cur = s->cur_l1_index + 1;
            fwrite(&new_cur, 4, 1, fp);
            s->cur_jmpf_offset = new_jmpf;
            s->cur_jmpf_entries = 0;
            fseek(fp, 0, SEEK_END);
        }
        long jmpf_entry = (long)s->cur_jmpf_offset + 12 + 8 + (long)s->cur_jmpf_entries * 32;
        fseek(fp, jmpf_entry, SEEK_SET);
        fwrite(&time, 8, 1, fp);
        fwrite(&frame_offset, 8, 1, fp);
        fwrite(&comp_size, 8, 1, fp);
        fwrite(&checksum, 4, 1, fp);
        fwrite(&frame_type, 4, 1, fp);
        s->cur_jmpf_entries++;
        fseek(fp, (long)s->cur_jmpf_offset + 12 + 4, SEEK_SET);
        fwrite(&s->cur_jmpf_entries, 4, 1, fp);
        long jtop_fc = (long)s->cur_jtop_offset + 12 + 16 + s->cur_l1_index * 16 + 12;
        fseek(fp, jtop_fc, SEEK_SET);
        fwrite(&s->cur_jmpf_entries, 4, 1, fp);
        s->total_frames++;
        fseek(fp, 0, SEEK_END);
    }
    return 0;
}

/* ---- FMSH (Voronoi mesh frame) ---- */
int sfa_write_mesh_frame(SFA *s, double time,
                         uint32_t N_cells, uint32_t N_faces,
                         double L,
                         const double *cell_pos,
                         const double *cell_vol,
                         const void *face_records,
                         const uint32_t *csr_off,
                         const uint32_t *csr_idx,
                         uint8_t flags) {
    int has_faces = (flags & 0x01) ? 1 : 0;
    int has_csr   = (flags & 0x02) ? 1 : 0;
    int pos_f64   = (flags & 0x04) ? 1 : 0;
    int pos_size  = pos_f64 ? 8 : 4;

    if (has_faces && !face_records) return -1;
    if (has_csr && (!csr_off || !csr_idx)) return -1;
    if (has_csr && !has_faces) return -1;  /* CSR meaningless without faces */

    uint64_t header_size = 28;     /* magic(4)+version(4)+N_cells(4)+N_faces(4)+flags(1)+resv(3)+L(8) */
    uint64_t pos_bytes  = (uint64_t)N_cells * 3 * pos_size;
    uint64_t vol_bytes  = (uint64_t)N_cells * pos_size;
    uint64_t face_bytes = has_faces ? (uint64_t)N_faces * 40 : 0;
    uint64_t csr_off_bytes = has_csr ? ((uint64_t)N_cells + 1) * 4 : 0;
    uint64_t csr_idx_total = has_csr ? csr_off[N_cells] : 0;
    uint64_t csr_idx_bytes = has_csr ? csr_idx_total * 4 : 0;
    uint64_t payload_size = header_size + pos_bytes + vol_bytes
                          + face_bytes + csr_off_bytes + csr_idx_bytes;

    uint8_t *payload = (uint8_t*)malloc(payload_size);
    if (!payload) return -1;
    uint64_t off = 0;
    uint32_t magic   = SFA_CHUNK_FMSH;
    uint32_t version = 1;
    uint8_t  resv[3] = {0,0,0};
    memcpy(payload + off, &magic, 4);     off += 4;
    memcpy(payload + off, &version, 4);   off += 4;
    memcpy(payload + off, &N_cells, 4);   off += 4;
    memcpy(payload + off, &N_faces, 4);   off += 4;
    memcpy(payload + off, &flags, 1);     off += 1;
    memcpy(payload + off, resv, 3);       off += 3;
    memcpy(payload + off, &L, 8);         off += 8;

    if (pos_f64) {
        memcpy(payload + off, cell_pos, pos_bytes); off += pos_bytes;
        memcpy(payload + off, cell_vol, vol_bytes); off += vol_bytes;
    } else {
        float *fpos = (float*)(payload + off);
        for (uint64_t i = 0; i < (uint64_t)N_cells * 3; i++) fpos[i] = (float)cell_pos[i];
        off += pos_bytes;
        float *fvol = (float*)(payload + off);
        for (uint64_t i = 0; i < N_cells; i++) fvol[i] = (float)cell_vol[i];
        off += vol_bytes;
    }
    if (has_faces) {
        memcpy(payload + off, face_records, face_bytes); off += face_bytes;
    }
    if (has_csr) {
        memcpy(payload + off, csr_off, csr_off_bytes); off += csr_off_bytes;
        memcpy(payload + off, csr_idx, csr_idx_bytes); off += csr_idx_bytes;
    }

    int ret = sfa_write_typed_chunk(s, time, SFA_CHUNK_FMSH,
                                     SFA_FRAME_MESH, payload, payload_size);
    free(payload);
    return ret;
}

/* ---- FCEL (cell-data frame) ---- */
int sfa_write_cell_frame(SFA *s, double time,
                         uint32_t N_cells, uint32_t n_columns,
                         uint8_t dtype, void * const *column_data) {
    if (dtype > SFA_F128) return -1;
    int es = sfa_dtype_size[dtype];
    uint64_t header_size = 20;     /* magic(4)+version(4)+N_cells(4)+n_columns(4)+dtype(1)+flags(1)+resv(2) */
    uint64_t col_bytes   = (uint64_t)N_cells * es;
    uint64_t payload_size = header_size + (uint64_t)n_columns * col_bytes;

    uint8_t *payload = (uint8_t*)malloc(payload_size);
    if (!payload) return -1;
    uint64_t off = 0;
    uint32_t magic   = SFA_CHUNK_FCEL;
    uint32_t version = 1;
    uint8_t  flags   = 0;
    uint8_t  resv[2] = {0,0};
    memcpy(payload + off, &magic, 4);       off += 4;
    memcpy(payload + off, &version, 4);     off += 4;
    memcpy(payload + off, &N_cells, 4);     off += 4;
    memcpy(payload + off, &n_columns, 4);   off += 4;
    memcpy(payload + off, &dtype, 1);       off += 1;
    memcpy(payload + off, &flags, 1);       off += 1;
    memcpy(payload + off, resv, 2);         off += 2;

    for (uint32_t c = 0; c < n_columns; c++) {
        memcpy(payload + off, column_data[c], col_bytes);
        off += col_bytes;
    }

    int ret = sfa_write_typed_chunk(s, time, SFA_CHUNK_FCEL,
                                     SFA_FRAME_CELL, payload, payload_size);
    free(payload);
    return ret;
}

/* ---- Reading FMSH / FCEL frames ----
 * Both follow the same flow as sfa_read_frame: locate the JMPF entry via
 * sfa_find_frame, read the compressed payload, decompress, parse. We
 * bypass the voxel frame_bytes machinery since cell/mesh frames carry
 * their own size headers. */
static int sfa_find_frame(SFA *s, uint32_t frame_idx, SFA_L2Entry *out);  /* forward decl */

static int sfa_read_typed_chunk(SFA *s, uint32_t frame_idx,
                                 void **out_payload, uint64_t *out_size,
                                 uint32_t *out_frame_type) {
    if (s->mode != 0) return -1;
    SFA_L2Entry e;
    if (sfa_find_frame(s, frame_idx, &e) < 0) return -1;
    *out_frame_type = e.frame_type;

    /* Chunk layout: type(4) + size(8) + time(8) + compressed_data */
    fseek(s->fp, (long)e.offset + 12 + 8, SEEK_SET);
    uint8_t *comp = (uint8_t*)malloc(e.compressed_size);
    if (!comp) return -1;
    if (fread(comp, 1, e.compressed_size, s->fp) != e.compressed_size) {
        free(comp); return -1;
    }
    unsigned long long content = ZSTD_getFrameContentSize(comp, e.compressed_size);
    if (content == ZSTD_CONTENTSIZE_ERROR || content == ZSTD_CONTENTSIZE_UNKNOWN) {
        content = (unsigned long long)e.compressed_size * 200 + 1024;
    }
    uint8_t *payload = (uint8_t*)malloc(content);
    if (!payload) { free(comp); return -1; }
    size_t got = ZSTD_decompress(payload, content, comp, e.compressed_size);
    free(comp);
    if (ZSTD_isError(got)) { free(payload); return -1; }
    *out_payload = payload;
    *out_size = got;
    return 0;
}

int sfa_read_mesh_frame(SFA *s, uint32_t frame_idx, SFA_Mesh *out) {
    void *payload = NULL; uint64_t psize = 0; uint32_t ftype = 0;
    if (sfa_read_typed_chunk(s, frame_idx, &payload, &psize, &ftype) != 0) return -1;
    if (ftype != SFA_FRAME_MESH) { free(payload); return -1; }
    uint8_t *p = (uint8_t*)payload;
    uint64_t off = 0;
    uint32_t magic, version;
    memcpy(&magic, p + off, 4);   off += 4;
    if (magic != SFA_CHUNK_FMSH) { free(payload); return -1; }
    memcpy(&version, p + off, 4); off += 4;
    memcpy(&out->N_cells, p + off, 4); off += 4;
    memcpy(&out->N_faces, p + off, 4); off += 4;
    memcpy(&out->flags, p + off, 1);   off += 1;
    off += 3;
    memcpy(&out->L, p + off, 8); off += 8;

    int has_faces = (out->flags & 0x01) ? 1 : 0;
    int has_csr   = (out->flags & 0x02) ? 1 : 0;
    int pos_f64   = (out->flags & 0x04) ? 1 : 0;
    int pos_size  = pos_f64 ? 8 : 4;

    out->cell_pos = (double*)malloc(sizeof(double) * 3 * out->N_cells);
    out->cell_vol = (double*)malloc(sizeof(double) * out->N_cells);
    if (pos_f64) {
        memcpy(out->cell_pos, p + off, 3 * out->N_cells * 8);
        off += 3 * out->N_cells * 8;
        memcpy(out->cell_vol, p + off, out->N_cells * 8);
        off += out->N_cells * 8;
    } else {
        const float *fp_pos = (const float*)(p + off);
        for (uint64_t i = 0; i < (uint64_t)3 * out->N_cells; i++)
            out->cell_pos[i] = (double)fp_pos[i];
        off += 3 * out->N_cells * 4;
        const float *fp_vol = (const float*)(p + off);
        for (uint64_t i = 0; i < out->N_cells; i++)
            out->cell_vol[i] = (double)fp_vol[i];
        off += out->N_cells * 4;
    }
    if (has_faces) {
        out->face_records = malloc((uint64_t)out->N_faces * 40);
        memcpy(out->face_records, p + off, (uint64_t)out->N_faces * 40);
        off += (uint64_t)out->N_faces * 40;
    } else { out->face_records = NULL; }
    if (has_csr) {
        out->csr_off = (uint32_t*)malloc(((uint64_t)out->N_cells + 1) * 4);
        memcpy(out->csr_off, p + off, ((uint64_t)out->N_cells + 1) * 4);
        off += ((uint64_t)out->N_cells + 1) * 4;
        uint32_t total = out->csr_off[out->N_cells];
        out->csr_idx = (uint32_t*)malloc((uint64_t)total * 4);
        memcpy(out->csr_idx, p + off, (uint64_t)total * 4);
        off += (uint64_t)total * 4;
    } else { out->csr_off = NULL; out->csr_idx = NULL; }
    free(payload);
    return 0;
}

int sfa_read_cell_frame(SFA *s, uint32_t frame_idx,
                        uint32_t *out_N_cells, uint32_t *out_n_columns,
                        uint8_t *out_dtype, void **out_data) {
    void *payload = NULL; uint64_t psize = 0; uint32_t ftype = 0;
    if (sfa_read_typed_chunk(s, frame_idx, &payload, &psize, &ftype) != 0) return -1;
    if (ftype != SFA_FRAME_CELL) { free(payload); return -1; }
    uint8_t *p = (uint8_t*)payload;
    uint64_t off = 0;
    uint32_t magic;
    memcpy(&magic, p + off, 4); off += 4;
    if (magic != SFA_CHUNK_FCEL) { free(payload); return -1; }
    off += 4;   /* version */
    memcpy(out_N_cells, p + off, 4);    off += 4;
    memcpy(out_n_columns, p + off, 4);  off += 4;
    memcpy(out_dtype, p + off, 1);      off += 1;
    off += 3;   /* flags + reserved */
    int es = sfa_dtype_size[*out_dtype];
    uint64_t data_bytes = (uint64_t)*out_N_cells * *out_n_columns * es;
    *out_data = malloc(data_bytes);
    memcpy(*out_data, p + off, data_bytes);
    free(payload);
    return 0;
}

void sfa_mesh_free(SFA_Mesh *m) {
    if (!m) return;
    free(m->cell_pos); free(m->cell_vol);
    free(m->face_records);
    free(m->csr_off); free(m->csr_idx);
    memset(m, 0, sizeof(*m));
}

/* ---- Cell I-frame with temporal model (FCEL v2) ---- */
int sfa_write_cell_iframe_temporal(SFA *s, double time,
                                    uint32_t N_cells, uint32_t n_columns,
                                    uint8_t dtype, void * const *column_data,
                                    float omega,
                                    const float *mean,
                                    const float *amp,
                                    const float *phase) {
    if (dtype > SFA_F128) return -1;
    int es = sfa_dtype_size[dtype];
    uint64_t header_size = 20;
    uint64_t col_bytes   = (uint64_t)N_cells * es;
    uint64_t cells_bytes = (uint64_t)n_columns * col_bytes;
    uint64_t model_bytes = 4 + 3 * (uint64_t)n_columns * (uint64_t)N_cells * 4;
    uint64_t payload_size = header_size + cells_bytes + model_bytes;

    uint8_t *payload = (uint8_t*)malloc(payload_size);
    if (!payload) return -1;
    uint64_t off = 0;
    uint32_t magic   = SFA_CHUNK_FCEL;
    uint32_t version = 2;     /* version 2: with temporal model */
    uint8_t  flags   = 0x04;  /* bit2 = has temporal model */
    uint8_t  resv[2] = {0,0};
    memcpy(payload + off, &magic, 4);     off += 4;
    memcpy(payload + off, &version, 4);   off += 4;
    memcpy(payload + off, &N_cells, 4);   off += 4;
    memcpy(payload + off, &n_columns, 4); off += 4;
    memcpy(payload + off, &dtype, 1);     off += 1;
    memcpy(payload + off, &flags, 1);     off += 1;
    memcpy(payload + off, resv, 2);       off += 2;

    for (uint32_t c = 0; c < n_columns; c++) {
        memcpy(payload + off, column_data[c], col_bytes);
        off += col_bytes;
    }
    memcpy(payload + off, &omega, 4); off += 4;
    uint64_t mb = (uint64_t)n_columns * N_cells * 4;
    memcpy(payload + off, mean,  mb); off += mb;
    memcpy(payload + off, amp,   mb); off += mb;
    memcpy(payload + off, phase, mb); off += mb;

    int ret = sfa_write_typed_chunk(s, time, SFA_CHUNK_FCEL,
                                     SFA_FRAME_CELL, payload, payload_size);
    free(payload);
    return ret;
}

/* ---- Cell P-frame (FCEP, sparse delta) ---- */
int sfa_write_cell_pframe(SFA *s, double time,
                          uint32_t N_cells, uint32_t n_columns,
                          uint32_t n_changed,
                          const uint32_t *cell_ids,
                          const float *delta_values) {
    /* Header: magic(4)+ver(4)+N_cells(4)+n_cols(4)+n_changed(4)+dtype(1)+flags(1)+resv(2) = 24 */
    uint64_t header_size = 24;
    uint64_t ids_bytes   = (uint64_t)n_changed * 4;
    uint64_t delta_bytes = (uint64_t)n_changed * n_columns * 4;
    uint64_t payload_size = header_size + ids_bytes + delta_bytes;

    uint8_t *payload = (uint8_t*)malloc(payload_size);
    if (!payload) return -1;
    uint64_t off = 0;
    uint32_t magic   = SFA_CHUNK_FCEP;
    uint32_t version = 1;
    uint8_t  dtype   = SFA_F32;  /* deltas always f32 */
    uint8_t  flags   = 0;
    uint8_t  resv[2] = {0,0};
    memcpy(payload + off, &magic, 4);      off += 4;
    memcpy(payload + off, &version, 4);    off += 4;
    memcpy(payload + off, &N_cells, 4);    off += 4;
    memcpy(payload + off, &n_columns, 4);  off += 4;
    memcpy(payload + off, &n_changed, 4);  off += 4;
    memcpy(payload + off, &dtype, 1);      off += 1;
    memcpy(payload + off, &flags, 1);      off += 1;
    memcpy(payload + off, resv, 2);        off += 2;
    if (n_changed > 0) {
        memcpy(payload + off, cell_ids, ids_bytes);   off += ids_bytes;
        memcpy(payload + off, delta_values, delta_bytes); off += delta_bytes;
    }

    int ret = sfa_write_typed_chunk(s, time, SFA_CHUNK_FCEP,
                                     SFA_FRAME_CELL_P, payload, payload_size);
    free(payload);
    return ret;
}

/* Reader for FCEP. Caller frees out_cell_ids and out_delta_values. */
int sfa_read_cell_pframe(SFA *s, uint32_t frame_idx,
                         uint32_t *out_N_cells, uint32_t *out_n_columns,
                         uint32_t *out_n_changed,
                         uint32_t **out_cell_ids,
                         float **out_delta_values) {
    void *payload = NULL; uint64_t psize = 0; uint32_t ftype = 0;
    if (sfa_read_typed_chunk(s, frame_idx, &payload, &psize, &ftype) != 0) return -1;
    if (ftype != SFA_FRAME_CELL_P) { free(payload); return -1; }
    uint8_t *p = (uint8_t*)payload;
    uint64_t off = 0;
    uint32_t magic;
    memcpy(&magic, p + off, 4); off += 4;
    if (magic != SFA_CHUNK_FCEP) { free(payload); return -1; }
    off += 4;   /* version */
    memcpy(out_N_cells,    p + off, 4); off += 4;
    memcpy(out_n_columns,  p + off, 4); off += 4;
    memcpy(out_n_changed,  p + off, 4); off += 4;
    /* dtype + flags + reserved skipped; we always read F32 */
    off += 4;
    uint32_t nc = *out_n_changed;
    *out_cell_ids = NULL;
    *out_delta_values = NULL;
    if (nc > 0) {
        *out_cell_ids = (uint32_t*)malloc((uint64_t)nc * 4);
        memcpy(*out_cell_ids, p + off, (uint64_t)nc * 4);
        off += (uint64_t)nc * 4;
        *out_delta_values = (float*)malloc((uint64_t)nc * (*out_n_columns) * 4);
        memcpy(*out_delta_values, p + off,
               (uint64_t)nc * (*out_n_columns) * 4);
    }
    free(payload);
    return 0;
}

/* ---- FMTL: standalone temporal-model frame ---- */
int sfa_write_temporal_model_frame(SFA *s, double time,
                                    uint32_t N_cells, uint32_t n_columns,
                                    float omega,
                                    const float *mean,
                                    const float *amp,
                                    const float *phase) {
    /* Header: magic(4)+ver(4)+N_cells(4)+n_cols(4)+flags(1)+resv(3)+omega(4)+pad(4) = 28 */
    uint64_t header_size = 28;
    uint64_t mb = (uint64_t)n_columns * N_cells * 4;
    uint64_t payload_size = header_size + 3 * mb;

    uint8_t *payload = (uint8_t*)malloc(payload_size);
    if (!payload) return -1;
    uint64_t off = 0;
    uint32_t magic = SFA_CHUNK_FMTL;
    uint32_t version = 1;
    uint8_t  flags = 0;
    uint8_t  resv[3] = {0,0,0};
    uint32_t pad = 0;
    memcpy(payload + off, &magic, 4);     off += 4;
    memcpy(payload + off, &version, 4);   off += 4;
    memcpy(payload + off, &N_cells, 4);   off += 4;
    memcpy(payload + off, &n_columns, 4); off += 4;
    memcpy(payload + off, &flags, 1);     off += 1;
    memcpy(payload + off, resv, 3);       off += 3;
    memcpy(payload + off, &omega, 4);     off += 4;
    memcpy(payload + off, &pad, 4);       off += 4;
    memcpy(payload + off, mean,  mb);     off += mb;
    memcpy(payload + off, amp,   mb);     off += mb;
    memcpy(payload + off, phase, mb);     off += mb;

    int ret = sfa_write_typed_chunk(s, time, SFA_CHUNK_FMTL,
                                     SFA_FRAME_TEMPORAL_MODEL,
                                     payload, payload_size);
    free(payload);
    return ret;
}

int sfa_read_temporal_model_frame(SFA *s, uint32_t frame_idx,
                                   float *out_omega,
                                   uint32_t *out_N_cells,
                                   uint32_t *out_n_columns,
                                   float **out_mean,
                                   float **out_amp,
                                   float **out_phase) {
    void *payload = NULL; uint64_t psize = 0; uint32_t ftype = 0;
    if (sfa_read_typed_chunk(s, frame_idx, &payload, &psize, &ftype) != 0) return -1;
    if (ftype != SFA_FRAME_TEMPORAL_MODEL) { free(payload); return -1; }
    uint8_t *p = (uint8_t*)payload;
    uint32_t magic;
    memcpy(&magic, p, 4);
    if (magic != SFA_CHUNK_FMTL) { free(payload); return -1; }
    uint32_t Nc, ncols;
    memcpy(&Nc, p + 8, 4);
    memcpy(&ncols, p + 12, 4);
    /* flags p[16], resv p[17..19] */
    memcpy(out_omega, p + 20, 4);
    /* pad p[24..27] */
    *out_N_cells = Nc;
    *out_n_columns = ncols;
    uint64_t mb = (uint64_t)ncols * Nc * 4;
    if (28 + 3 * mb > psize) { free(payload); return -1; }
    *out_mean  = (float*)malloc(mb); memcpy(*out_mean,  p + 28,        mb);
    *out_amp   = (float*)malloc(mb); memcpy(*out_amp,   p + 28 + mb,   mb);
    *out_phase = (float*)malloc(mb); memcpy(*out_phase, p + 28 + 2*mb, mb);
    free(payload);
    return 0;
}

/* Reader for the temporal model embedded in an FCEL v2 I-frame. */
int sfa_read_cell_iframe_temporal(SFA *s, uint32_t frame_idx,
                                   float *out_omega,
                                   uint32_t *out_N_cells,
                                   uint32_t *out_n_columns,
                                   float **out_mean,
                                   float **out_amp,
                                   float **out_phase) {
    void *payload = NULL; uint64_t psize = 0; uint32_t ftype = 0;
    if (sfa_read_typed_chunk(s, frame_idx, &payload, &psize, &ftype) != 0) return -1;
    if (ftype != SFA_FRAME_CELL) { free(payload); return -1; }
    uint8_t *p = (uint8_t*)payload;
    uint64_t off = 0;
    uint32_t magic, version;
    memcpy(&magic, p + off, 4); off += 4;
    if (magic != SFA_CHUNK_FCEL) { free(payload); return -1; }
    memcpy(&version, p + off, 4); off += 4;
    uint32_t Nc, ncols;
    memcpy(&Nc, p + off, 4);    off += 4;
    memcpy(&ncols, p + off, 4); off += 4;
    uint8_t dtype = p[off]; off += 1;
    uint8_t flags = p[off]; off += 1;
    off += 2;   /* reserved */
    if (version < 2 || (flags & 0x04) == 0) { free(payload); return -1; }
    int es = sfa_dtype_size[dtype];
    /* Skip cell values */
    off += (uint64_t)ncols * Nc * es;
    /* Read model */
    if (off + 4 > psize) { free(payload); return -1; }
    memcpy(out_omega, p + off, 4); off += 4;
    *out_N_cells   = Nc;
    *out_n_columns = ncols;
    uint64_t mb = (uint64_t)ncols * Nc * 4;
    if (off + 3*mb > psize) { free(payload); return -1; }
    *out_mean  = (float*)malloc(mb); memcpy(*out_mean,  p + off, mb); off += mb;
    *out_amp   = (float*)malloc(mb); memcpy(*out_amp,   p + off, mb); off += mb;
    *out_phase = (float*)malloc(mb); memcpy(*out_phase, p + off, mb); off += mb;
    free(payload);
    return 0;
}

static int sfa_write_vec_chunk(SFA *s, double time, uint32_t frame_type,
                                const void *payload, uint64_t payload_size) {
    FILE *fp = s->fp;
    if (!fp || s->mode != 1) return -1;

    /* Zstd compress the payload */
    size_t bound = ZSTD_compressBound(payload_size);
    void *comp = malloc(bound);
    size_t comp_size = ZSTD_compress(comp, bound, payload, payload_size, 3);
    if (ZSTD_isError(comp_size)) { free(comp); return -1; }

    /* Write FRVD chunk: type(4) + size(8) + time(8) + compressed_data */
    fseek(fp, 0, SEEK_END);
    uint64_t frame_offset = ftell(fp);

    uint32_t chunk_type = SFA_CHUNK_FRVD;
    uint64_t chunk_size = 8 + comp_size;  /* time(8) + compressed data */
    fwrite(&chunk_type, 4, 1, fp);
    fwrite(&chunk_size, 8, 1, fp);
    fwrite(&time, 8, 1, fp);
    fwrite(comp, 1, comp_size, fp);
    free(comp);

    /* Update JMPF index (same as voxel frames) */
    uint32_t checksum = sfa_crc32(payload, payload_size);
    if (s->cur_jmpf_offset) {
        /* Handle JMPF overflow — allocate new JMPF if full */
        if (s->cur_jmpf_entries >= s->jmpf_max) {
            fseek(fp, 0, SEEK_END);
            uint64_t new_jmpf = (uint64_t)ftell(fp);
            uint64_t jmpf_size = 12 + 8 + (uint64_t)s->jmpf_max * 32;
            sfa_write_chunk_header(fp, "JMPF", jmpf_size);
            fwrite(&s->jmpf_max, 4, 1, fp);
            uint32_t zero32 = 0;
            fwrite(&zero32, 4, 1, fp);
            uint64_t zero64 = 0; double zt = 0;
            for (uint32_t i = 0; i < s->jmpf_max; i++) {
                fwrite(&zt, 8, 1, fp); fwrite(&zero64, 8, 1, fp);
                fwrite(&zero64, 8, 1, fp); fwrite(&zero32, 4, 1, fp);
                fwrite(&zero32, 4, 1, fp);
            }
            s->cur_l1_index++;
            long jtop_entry = (long)s->cur_jtop_offset + 12 + 16 + s->cur_l1_index * 16;
            fseek(fp, jtop_entry, SEEK_SET);
            fwrite(&new_jmpf, 8, 1, fp);
            uint32_t first_frame = s->total_frames;
            fwrite(&first_frame, 4, 1, fp);
            fwrite(&zero32, 4, 1, fp);
            fseek(fp, (long)s->cur_jtop_offset + 12 + 4, SEEK_SET);
            uint32_t new_cur = s->cur_l1_index + 1;
            fwrite(&new_cur, 4, 1, fp);
            s->cur_jmpf_offset = new_jmpf;
            s->cur_jmpf_entries = 0;
            fseek(fp, 0, SEEK_END);
        }

        /* JMPF layout: chunk_header(12) + max_entries(4) + current_entries(4) + entries... */
        long jmpf_entry = (long)s->cur_jmpf_offset + 12 + 8 + (long)s->cur_jmpf_entries * 32;
        fseek(fp, jmpf_entry, SEEK_SET);
        fwrite(&time, 8, 1, fp);
        fwrite(&frame_offset, 8, 1, fp);
        fwrite(&comp_size, 8, 1, fp);
        fwrite(&checksum, 4, 1, fp);
        fwrite(&frame_type, 4, 1, fp);
        s->cur_jmpf_entries++;
        /* Update current_entries in JMPF header */
        fseek(fp, (long)s->cur_jmpf_offset + 12 + 4, SEEK_SET);
        fwrite(&s->cur_jmpf_entries, 4, 1, fp);

        /* Update JTOP frame_count */
        long jtop_fc = (long)s->cur_jtop_offset + 12 + 16 + s->cur_l1_index * 16 + 12;
        fseek(fp, jtop_fc, SEEK_SET);
        fwrite(&s->cur_jmpf_entries, 4, 1, fp);

        s->total_frames++;
        fseek(fp, 0, SEEK_END);
    }
    return 0;
}

int sfa_write_vec_iframe(SFA *s, double time,
                         uint32_t n_patches, uint8_t block_size, uint16_t n_coeffs,
                         const int16_t *origins, const float *coeffs) {
    /* Build payload: header + origins + coefficients */
    uint8_t pad = 0;
    uint64_t origin_bytes = (uint64_t)n_patches * 6;  /* 3 × int16 */
    uint64_t coeff_bytes = (uint64_t)n_patches * n_coeffs * 4;
    uint64_t payload_size = 8 + origin_bytes + coeff_bytes;  /* header(8) + data */

    uint8_t *payload = (uint8_t*)malloc(payload_size);
    uint64_t off = 0;

    memcpy(payload + off, &n_patches, 4); off += 4;
    memcpy(payload + off, &block_size, 1); off += 1;
    memcpy(payload + off, &n_coeffs, 2); off += 2;
    memcpy(payload + off, &pad, 1); off += 1;
    memcpy(payload + off, origins, origin_bytes); off += origin_bytes;
    memcpy(payload + off, coeffs, coeff_bytes); off += coeff_bytes;

    int ret = sfa_write_vec_chunk(s, time, SFA_FRAME_VEC_I, payload, payload_size);
    free(payload);
    return ret;
}

int sfa_write_vec_iframe_temporal(SFA *s, double time,
                         uint32_t n_patches, uint8_t block_size, uint16_t n_coeffs,
                         const int16_t *origins, const float *coeffs,
                         float omega, const float *temp_mean,
                         const float *temp_amp, const float *temp_phase) {
    uint8_t flags = 0x01;  /* bit 0 = has temporal model */
    uint64_t n_total = (uint64_t)n_patches * n_coeffs;
    uint64_t origin_bytes = (uint64_t)n_patches * 6;
    uint64_t coeff_bytes = n_total * 4;
    uint64_t model_bytes = 4 + n_total * 4 * 3;  /* omega + mean + amp + phase */
    uint64_t payload_size = 8 + origin_bytes + coeff_bytes + model_bytes;

    uint8_t *payload = (uint8_t*)malloc(payload_size);
    uint64_t off = 0;

    memcpy(payload + off, &n_patches, 4); off += 4;
    memcpy(payload + off, &block_size, 1); off += 1;
    memcpy(payload + off, &n_coeffs, 2); off += 2;
    memcpy(payload + off, &flags, 1); off += 1;
    memcpy(payload + off, origins, origin_bytes); off += origin_bytes;
    memcpy(payload + off, coeffs, coeff_bytes); off += coeff_bytes;
    memcpy(payload + off, &omega, 4); off += 4;
    memcpy(payload + off, temp_mean, n_total * 4); off += n_total * 4;
    memcpy(payload + off, temp_amp, n_total * 4); off += n_total * 4;
    memcpy(payload + off, temp_phase, n_total * 4); off += n_total * 4;

    int ret = sfa_write_vec_chunk(s, time, SFA_FRAME_VEC_I, payload, payload_size);
    free(payload);
    return ret;
}

int sfa_write_vec_pframe(SFA *s, double time,
                         uint32_t n_deltas, uint16_t n_coeffs,
                         const uint32_t *indices, const float *delta_coeffs) {
    uint16_t pad = 0;
    uint64_t idx_bytes = (uint64_t)n_deltas * 4;
    uint64_t coeff_bytes = (uint64_t)n_deltas * n_coeffs * 4;
    uint64_t payload_size = 8 + idx_bytes + coeff_bytes;

    uint8_t *payload = (uint8_t*)malloc(payload_size);
    uint64_t off = 0;

    memcpy(payload + off, &n_deltas, 4); off += 4;
    memcpy(payload + off, &n_coeffs, 2); off += 2;
    memcpy(payload + off, &pad, 2); off += 2;
    if (n_deltas > 0 && indices && delta_coeffs) {
        memcpy(payload + off, indices, idx_bytes); off += idx_bytes;
        memcpy(payload + off, delta_coeffs, coeff_bytes); off += coeff_bytes;
    }

    int ret = sfa_write_vec_chunk(s, time, SFA_FRAME_VEC_P, payload, payload_size);
    free(payload);
    return ret;
}

static int sfa_find_frame(SFA *s, uint32_t frame_idx, SFA_L2Entry *out);  /* forward decl */

uint32_t sfa_frame_type(SFA *s, uint32_t frame_idx) {
    SFA_L2Entry entry;
    if (sfa_find_frame(s, frame_idx, &entry) != 0) return SFA_FRAME_VOXEL;
    return entry.frame_type;
}

void sfa_enable_split(SFA *s) {
    s->split_output = 1;
    /* The caller must set split_base before calling finalize_header.
     * Typically: strncpy(sfa->split_base, "output", 511) after sfa_create. */
}

void sfa_close(SFA *s) {
    if (s->mode == 1 && s->fp) {
        if (s->split_output) {
            /* Split mode: patch total_frames but leave streaming flag
             * (JMPF entries are empty — fixup needed after combine) */
            fseek(s->fp, 12 + 4 + 4 + 4*3 + 8*3 + 8 + 4, SEEK_SET);
            fwrite(&s->total_frames, 4, 1, s->fp);
        } else {
            /* Patch total_frames in SFAH (offset 68) */
            fseek(s->fp, 12 + 4 + 4 + 4*3 + 8*3 + 8 + 4, SEEK_SET);
            fwrite(&s->total_frames, 4, 1, s->fp);

            /* Clear streaming flag — file is now complete with valid indexes */
            uint32_t flags = s->flags & ~SFA_FLAG_STREAMING;
            fseek(s->fp, 12 + 4, SEEK_SET);  /* offset 16: flags field */
            fwrite(&flags, 4, 1, s->fp);
        }
    }
    if (s->fp) fclose(s->fp);
    if (s->columns) free(s->columns);
    free(s);
}

int sfa_combine(const char *header_path, const char *output_path,
                const char **sfp_paths, int n_sfp) {
    /* Copy header .sfa to output */
    FILE *hdr = fopen(header_path, "rb");
    FILE *out = fopen(output_path, "wb");
    if (!hdr || !out) {
        if (hdr) fclose(hdr);
        if (out) fclose(out);
        return -1;
    }
    /* Copy header file entirely */
    char buf[65536];
    size_t n;
    while ((n = fread(buf, 1, sizeof(buf), hdr)) > 0)
        fwrite(buf, 1, n, out);
    fclose(hdr);

    /* Append each .sfp's FRMD chunk (skip the 16-byte SFP header) */
    for (int i = 0; i < n_sfp; i++) {
        FILE *sfp = fopen(sfp_paths[i], "rb");
        if (!sfp) { fprintf(stderr, "sfa_combine: cannot open %s\n", sfp_paths[i]); continue; }

        /* Verify SFP magic */
        char magic[4];
        fread(magic, 1, 4, sfp);
        if (memcmp(magic, SFA_SFP_MAGIC, 4) != 0) {
            fprintf(stderr, "sfa_combine: bad magic in %s\n", sfp_paths[i]);
            fclose(sfp);
            continue;
        }
        /* Skip frame_idx(4) + time(8) = 12 more bytes */
        fseek(sfp, 12, SEEK_CUR);

        /* Copy FRMD chunk */
        while ((n = fread(buf, 1, sizeof(buf), sfp)) > 0)
            fwrite(buf, 1, n, out);
        fclose(sfp);
    }
    fclose(out);

    /* Run fixup to populate JMPF entries */
    return sfa_fixup_index(output_path);
}

/* ---- Reader ---- */

SFA *sfa_open(const char *path) {
    FILE *fp = fopen(path, "rb");
    if (!fp) return NULL;

    SFA *s = (SFA*)calloc(1, sizeof(SFA));
    s->fp = fp;
    s->mode = 0;

    /* Read SFAH */
    char type[4]; uint64_t size;
    fread(type, 1, 4, fp);
    fread(&size, 8, 1, fp);
    if (memcmp(type, "SFAH", 4) != 0) { fclose(fp); free(s); return NULL; }

    fread(&s->version, 4, 1, fp);
    fread(&s->flags, 4, 1, fp);
    fread(&s->Nx, 4, 1, fp);
    fread(&s->Ny, 4, 1, fp);
    fread(&s->Nz, 4, 1, fp);
    fread(&s->Lx, 8, 1, fp);
    fread(&s->Ly, 8, 1, fp);
    fread(&s->Lz, 8, 1, fp);
    fread(&s->dt, 8, 1, fp);
    fread(&s->n_columns, 4, 1, fp);
    fread(&s->total_frames, 4, 1, fp);
    fread(&s->first_jtop_offset, 8, 1, fp);
    fread(&s->cdef_offset, 8, 1, fp);
    fread(&s->jtop_max, 4, 1, fp);
    fread(&s->jmpf_max, 4, 1, fp);

    /* Read CDEF */
    fseek(fp, (long)s->cdef_offset, SEEK_SET);
    fread(type, 1, 4, fp);
    fread(&size, 8, 1, fp);
    s->columns = (SFA_Column*)calloc(s->n_columns, sizeof(SFA_Column));
    for (uint32_t c = 0; c < s->n_columns; c++) {
        fread(s->columns[c].name, 1, 12, fp);
        fread(&s->columns[c].dtype, 1, 1, fp);
        fread(&s->columns[c].semantic, 1, 1, fp);
        fread(&s->columns[c].component, 1, 1, fp);
        fread(&s->columns[c].flags, 1, 1, fp);
        fread(&s->columns[c].scale, 8, 1, fp);
    }

    sfa_compute_frame_bytes(s);

    /* Streaming files: total_frames in header is 0 until sfa_close().
     * The JTOP/JMPF index IS updated live during writes, so scan it
     * to get the actual frame count. This allows reading streaming files
     * while they're still being written. */
    /* NOTE: Streaming frame count scan disabled due to crashes on large
     * multi-chain JTOP files (>1024 frames). Tools that need frame counts
     * for streaming files should call sfa_count_valid_frames() explicitly.
     * TODO: Fix the streaming scan to handle multi-chain JTOP safely. */

    return s;
}

/* Find L2 entry for a frame */
static int sfa_find_frame(SFA *s, uint32_t frame_idx, SFA_L2Entry *out) {
    uint64_t jtop_off = s->first_jtop_offset;
    FILE *fp = s->fp;

    while (jtop_off != 0) {
        fseek(fp, (long)jtop_off + 12, SEEK_SET);  /* skip chunk header */
        uint32_t max_ent, cur_ent;
        uint64_t next_jtop;
        fread(&max_ent, 4, 1, fp);
        fread(&cur_ent, 4, 1, fp);
        fread(&next_jtop, 8, 1, fp);

        /* Which L2 table? */
        uint32_t l1_idx = frame_idx / s->jmpf_max;
        if (l1_idx < cur_ent) {
            /* Read JTOP entry */
            fseek(fp, (long)jtop_off + 12 + 16 + l1_idx * 16, SEEK_SET);
            SFA_L1Entry l1;
            fread(&l1.jmpf_offset, 8, 1, fp);
            fread(&l1.first_frame, 4, 1, fp);
            fread(&l1.frame_count, 4, 1, fp);

            uint32_t l2_idx = frame_idx % s->jmpf_max;
            if (l2_idx >= l1.frame_count) return -1;

            /* Read JMPF entry */
            fseek(fp, (long)l1.jmpf_offset + 12 + 8 + l2_idx * 32, SEEK_SET);
            fread(&out->time, 8, 1, fp);
            fread(&out->offset, 8, 1, fp);
            fread(&out->compressed_size, 8, 1, fp);
            fread(&out->checksum, 4, 1, fp);
            fread(&out->frame_type, 4, 1, fp);
            return 0;
        }

        /* Move to next JTOP (frames >= jtop_max * jmpf_max) */
        frame_idx -= cur_ent * s->jmpf_max;
        jtop_off = next_jtop;
    }
    return -1;
}

int sfa_read_frame(SFA *s, uint32_t frame_idx, void *buf) {
    SFA_L2Entry entry;
    if (sfa_find_frame(s, frame_idx, &entry) < 0) return -1;

    /* FRVD vector frame: decompress, evaluate patches into buf as if it were voxel data */
    if (entry.frame_type == SFA_FRAME_VEC_I || entry.frame_type == SFA_FRAME_VEC_P) {
        /* Read compressed payload (skip chunk header + time) */
        fseek(s->fp, (long)entry.offset + 20, SEEK_SET);
        void *comp = malloc(entry.compressed_size);
        if (!comp) return -1;
        if (fread(comp, 1, entry.compressed_size, s->fp) != entry.compressed_size) {
            free(comp); return -1;
        }
        /* Decompress — use ZSTD content size hint, fallback to generous estimate */
        unsigned long long content_size = ZSTD_getFrameContentSize(comp, entry.compressed_size);
        size_t bound;
        if (content_size != ZSTD_CONTENTSIZE_UNKNOWN && content_size != ZSTD_CONTENTSIZE_ERROR)
            bound = (size_t)content_size;
        else
            bound = s->frame_bytes > entry.compressed_size * 20 ? s->frame_bytes : entry.compressed_size * 20;
        /* For vec frames with temporal model, payload can be 4x frame_bytes */
        if (bound < s->frame_bytes * 4) bound = s->frame_bytes * 4;
        void *payload = malloc(bound);
        if (!payload) { free(comp); return -1; }
        size_t dec = ZSTD_decompress(payload, bound, comp, entry.compressed_size);
        free(comp);
        if (ZSTD_isError(dec)) { free(payload); return -1; }

        /* Parse header (shared by I and P frames) */
        uint8_t *p = (uint8_t*)payload;

        if (entry.frame_type == SFA_FRAME_VEC_I) {
            /* I-frame: [n_patches(4)][bs(1)][nc(2)][flags(1)][origins(n*6)][coeffs(n*nc*4)]
             * If flags & 1: [omega(4)][mean(n*nc*4)][amp(n*nc*4)][phase(n*nc*4)] */
            uint32_t np; memcpy(&np, p, 4); p += 4;
            int bs = p[0]; p++;
            uint16_t nc; memcpy(&nc, p, 2); p += 2;
            uint8_t flags = p[0]; p++;

            int16_t *origins = (int16_t*)p;
            p += np * 6;
            float *coeffs = (float*)p;
            p += (long)np * nc * 4;

            /* Store coefficients for P-frame reconstruction */
            long n_total = (long)np * nc;
            s->vec_n_patches = np;
            s->vec_n_coeffs = nc;
            free(s->vec_prev_coeffs);
            s->vec_prev_coeffs = (float*)malloc(n_total * sizeof(float));
            memcpy(s->vec_prev_coeffs, coeffs, n_total * sizeof(float));

            /* Parse temporal model if present */
            if (flags & 0x01) {
                float omega; memcpy(&omega, p, 4); p += 4;
                s->vec_temporal_omega = omega;
                free(s->vec_temporal_mean); free(s->vec_temporal_amp); free(s->vec_temporal_phase);
                s->vec_temporal_mean  = (float*)malloc(n_total * sizeof(float));
                s->vec_temporal_amp   = (float*)malloc(n_total * sizeof(float));
                s->vec_temporal_phase = (float*)malloc(n_total * sizeof(float));
                memcpy(s->vec_temporal_mean, p, n_total * 4); p += n_total * 4;
                memcpy(s->vec_temporal_amp, p, n_total * 4); p += n_total * 4;
                memcpy(s->vec_temporal_phase, p, n_total * 4); p += n_total * 4;
                s->vec_temporal_valid = 1;
            } else {
                s->vec_temporal_valid = 0;
            }

            /* Evaluate coefficients to voxel buffer */
            int n_fields = 1, pf_nc = nc;
            if (nc > 64) {
                if (nc % 64 == 0) { n_fields = nc / 64; pf_nc = 64; }
                else if (nc % 27 == 0) { n_fields = nc / 27; pf_nc = 27; }
                else if (nc % 8 == 0)  { n_fields = nc / 8;  pf_nc = 8; }
            }
            int order = 3;
            if (pf_nc <= 27) order = 2;
            if (pf_nc <= 8) order = 1;
            int o1 = order + 1;

            uint64_t col_offsets[64];
            { uint64_t off2 = 0;
              for (uint32_t cc = 0; cc < s->n_columns; cc++) {
                  col_offsets[cc] = off2;
                  off2 += s->N_total * sfa_dtype_size[s->columns[cc].dtype];
              }
            }
            memset(buf, 0, s->frame_bytes);
            int N = s->Nx;
            #ifdef _OPENMP
            #pragma omp parallel for schedule(dynamic)
            #endif
            for (uint32_t pi = 0; pi < np; pi++) {
                int ox = origins[pi*3], oy = origins[pi*3+1], oz = origins[pi*3+2];
                float *pc = coeffs + (long)pi * nc;
                for (int di = 0; di < bs; di++) {
                    int gi = ox + di; if (gi < 0 || gi >= N) continue;
                    double tx = (bs>1) ? (double)di/(bs-1) : 0;
                    double sx = 2*tx-1;
                    double txa[4] = {1, sx, 2*sx*sx-1, 4*sx*sx*sx-3*sx};
                    for (int dj = 0; dj < bs; dj++) {
                        int gj = oy + dj; if (gj < 0 || gj >= N) continue;
                        double ty = (bs>1) ? (double)dj/(bs-1) : 0;
                        double sy = 2*ty-1;
                        double tya[4] = {1, sy, 2*sy*sy-1, 4*sy*sy*sy-3*sy};
                        for (int dk = 0; dk < bs; dk++) {
                            int gk = oz + dk; if (gk < 0 || gk >= N) continue;
                            double tz = (bs>1) ? (double)dk/(bs-1) : 0;
                            double sz = 2*tz-1;
                            double tza[4] = {1, sz, 2*sz*sz-1, 4*sz*sz*sz-3*sz};
                            long idx = (long)gi*N*N + gj*N + gk;
                            for (int f = 0; f < n_fields && f < (int)s->n_columns; f++) {
                                float *fc = pc + f * pf_nc;
                                double val = 0;
                                for (int a = 0; a < o1; a++)
                                for (int b = 0; b < o1; b++) {
                                    double ab = txa[a]*tya[b];
                                    int base2 = a*o1*o1 + b*o1;
                                    for (int c = 0; c < o1; c++)
                                        val += fc[base2+c] * ab * tza[c];
                                }
                                float fv = (float)val;
                                int es = sfa_dtype_size[s->columns[f].dtype];
                                uint8_t *dst = (uint8_t*)buf + col_offsets[f] + idx * es;
                                if (es == 4) memcpy(dst, &fv, 4);
                                else if (es == 2) {
                                    uint32_t bits; memcpy(&bits, &fv, 4);
                                    uint16_t sign = (bits >> 16) & 0x8000;
                                    int exp2 = ((bits>>23)&0xFF) - 127 + 15;
                                    uint16_t mant = (bits >> 13) & 0x3FF;
                                    uint16_t h = (exp2<=0)?sign:(exp2>=31)?(sign|0x7C00):(sign|(exp2<<10)|mant);
                                    memcpy(dst, &h, 2);
                                }
                            }
                        }
                    }
                }
            }
        } else {
            /* P-frame: [n_deltas(4)][n_coeffs(2)][pad(2)][indices(n*4)][deltas(n*nc*4)]
             * Reconstruct: if temporal model, actual = predicted(t) + delta
             *              else actual = prev + delta */
            uint32_t nd; memcpy(&nd, p, 4); p += 4;
            uint16_t nc; memcpy(&nc, p, 2); p += 4; /* nc + pad */

            if (!s->vec_prev_coeffs || s->vec_n_coeffs != nc) {
                /* No baseline — can't reconstruct */
                free(payload);
                memset(buf, 0, s->frame_bytes);
                return 0;
            }

            long n_total = (long)s->vec_n_patches * nc;

            /* Build reconstructed coefficients */
            float *recon = (float*)malloc(n_total * sizeof(float));

            if (s->vec_temporal_valid) {
                /* predicted = mean + amp*cos(omega*t + phase) */
                float omega = s->vec_temporal_omega;
                for (long i = 0; i < n_total; i++)
                    recon[i] = s->vec_temporal_mean[i]
                             + s->vec_temporal_amp[i]
                             * cosf(omega * (float)entry.time + s->vec_temporal_phase[i]);
            } else {
                memcpy(recon, s->vec_prev_coeffs, n_total * sizeof(float));
            }

            /* Apply sparse deltas */
            uint32_t *indices = (uint32_t*)p;
            p += nd * 4;
            float *deltas = (float*)p;
            for (uint32_t d = 0; d < nd; d++) {
                uint32_t pi = indices[d];
                if (pi >= s->vec_n_patches) continue;
                long base2 = (long)pi * nc;
                for (int c = 0; c < nc; c++)
                    recon[base2 + c] += deltas[(long)d * nc + c];
            }

            /* Save as prev for next frame */
            memcpy(s->vec_prev_coeffs, recon, n_total * sizeof(float));

            /* Evaluate to voxels (same logic as I-frame) */
            int n_fields = 1, pf_nc = nc;
            if (nc > 64) {
                if (nc % 64 == 0) { n_fields = nc / 64; pf_nc = 64; }
                else if (nc % 27 == 0) { n_fields = nc / 27; pf_nc = 27; }
                else if (nc % 8 == 0)  { n_fields = nc / 8;  pf_nc = 8; }
            }
            int order = 3;
            if (pf_nc <= 27) order = 2;
            if (pf_nc <= 8) order = 1;
            int o1 = order + 1;
            int bs = 8;  /* block size from I-frame — TODO: store in SFA state */

            uint64_t col_offsets[64];
            { uint64_t off2 = 0;
              for (uint32_t cc = 0; cc < s->n_columns; cc++) {
                  col_offsets[cc] = off2;
                  off2 += s->N_total * sfa_dtype_size[s->columns[cc].dtype];
              }
            }
            memset(buf, 0, s->frame_bytes);
            int N = s->Nx;
            int BN = N / bs;
            #ifdef _OPENMP
            #pragma omp parallel for schedule(dynamic)
            #endif
            for (uint32_t pi = 0; pi < s->vec_n_patches; pi++) {
                int bk = pi % BN, bj = (pi / BN) % BN, bi = pi / (BN * BN);
                int ox = bi * bs, oy = bj * bs, oz = bk * bs;
                float *pc = recon + (long)pi * nc;
                for (int di = 0; di < bs; di++) {
                    int gi = ox + di; if (gi >= N) continue;
                    double tx = (bs>1) ? (double)di/(bs-1) : 0;
                    double sx = 2*tx-1;
                    double txa[4] = {1, sx, 2*sx*sx-1, 4*sx*sx*sx-3*sx};
                    for (int dj = 0; dj < bs; dj++) {
                        int gj = oy + dj; if (gj >= N) continue;
                        double ty = (bs>1) ? (double)dj/(bs-1) : 0;
                        double sy = 2*ty-1;
                        double tya[4] = {1, sy, 2*sy*sy-1, 4*sy*sy*sy-3*sy};
                        for (int dk = 0; dk < bs; dk++) {
                            int gk = oz + dk; if (gk >= N) continue;
                            double tz = (bs>1) ? (double)dk/(bs-1) : 0;
                            double sz = 2*tz-1;
                            double tza[4] = {1, sz, 2*sz*sz-1, 4*sz*sz*sz-3*sz};
                            long idx = (long)gi*N*N + gj*N + gk;
                            for (int f = 0; f < n_fields && f < (int)s->n_columns; f++) {
                                float *fc = pc + f * pf_nc;
                                double val = 0;
                                for (int a = 0; a < o1; a++)
                                for (int b = 0; b < o1; b++) {
                                    double ab = txa[a]*tya[b];
                                    int base2 = a*o1*o1 + b*o1;
                                    for (int c = 0; c < o1; c++)
                                        val += fc[base2+c] * ab * tza[c];
                                }
                                float fv = (float)val;
                                int es = sfa_dtype_size[s->columns[f].dtype];
                                uint8_t *dst = (uint8_t*)buf + col_offsets[f] + idx * es;
                                if (es == 4) memcpy(dst, &fv, 4);
                                else if (es == 2) {
                                    uint32_t bits; memcpy(&bits, &fv, 4);
                                    uint16_t sign = (bits >> 16) & 0x8000;
                                    int exp2 = ((bits>>23)&0xFF) - 127 + 15;
                                    uint16_t mant = (bits >> 13) & 0x3FF;
                                    uint16_t h = (exp2<=0)?sign:(exp2>=31)?(sign|0x7C00):(sign|(exp2<<10)|mant);
                                    memcpy(dst, &h, 2);
                                }
                            }
                        }
                    }
                }
            }
            free(recon);
        }
        free(payload);
        return 0;
    }

    int codec = s->flags & 0xF;

    /* Validate entry before allocating */
    if (entry.compressed_size == 0 || entry.compressed_size > 50000000000ULL) {
        /* Sanity: no frame should be >50 GB compressed */
        return -1;
    }

    /* Read compressed data */
    fseek(s->fp, (long)entry.offset + 12, SEEK_SET);  /* skip FRMD chunk header */
    void *comp = (void*)malloc(entry.compressed_size);
    if (!comp) return -1;
    size_t nread = fread(comp, 1, entry.compressed_size, s->fp);
    if (nread != entry.compressed_size) { free(comp); return -1; }

    if (codec == SFA_CODEC_COLZSTD) {
        /* Per-column zstd: [n_cols(4)] [comp_size(8)]×n_cols [data...] */
        uint8_t *p = (uint8_t*)comp;
        uint32_t nc;
        memcpy(&nc, p, 4); p += 4;
        if (nc != s->n_columns) {
            fprintf(stderr, "SFA: colzstd n_cols mismatch (%u vs %u)\n", nc, s->n_columns);
            free(comp);
            return -1;
        }
        uint64_t *col_sizes = (uint64_t*)p;
        p += nc * 8;

        /* Decompress and BSS-decode each column sequentially
         * (C version — Go version does parallel) */
        uint64_t out_off = 0;
        for (uint32_t c = 0; c < nc; c++) {
            int es = sfa_dtype_size[s->columns[c].dtype];
            uint64_t col_bytes = s->N_total * es;

            /* Decompress into temp buffer (BSS-encoded) */
            uint8_t *bss_col = (uint8_t*)malloc(col_bytes);
            size_t result = ZSTD_decompress(bss_col, col_bytes, p, col_sizes[c]);
            if (ZSTD_isError(result)) {
                fprintf(stderr, "SFA: colzstd decompress col %u: %s\n",
                        c, ZSTD_getErrorName(result));
                free(bss_col); free(comp);
                return -1;
            }

            /* BSS-decode: [b0v0,b0v1,...,b1v0,...] → [v0b0,v0b1,...,v1b0,v1b1,...] */
            uint8_t *dst_col = (uint8_t*)buf + out_off;
            for (int b = 0; b < es; b++) {
                for (uint64_t v = 0; v < s->N_total; v++) {
                    dst_col[v * es + b] = bss_col[b * s->N_total + v];
                }
            }
            free(bss_col);

            p += col_sizes[c];
            out_off += col_bytes;
        }
        free(comp);
        return 0;
    }

    /* Standard codecs (RAW, ZSTD, BSS, etc.) */
    if (codec == SFA_CODEC_RAW) {
        memcpy(buf, comp, s->frame_bytes);
        free(comp);
        return 0;
    }

    /* Decompress (into buf or temp if BSS needs reversal) */
    uint8_t *decomp_buf = (uint8_t *)buf;
    uint8_t *bss_temp = NULL;
    if (codec == SFA_CODEC_BSS) {
        bss_temp = (uint8_t*)malloc(s->frame_bytes);
        decomp_buf = bss_temp;
    }

    size_t result = ZSTD_decompress(decomp_buf, s->frame_bytes, comp, entry.compressed_size);
    free(comp);

    if (ZSTD_isError(result)) {
        fprintf(stderr, "SFA: decompress error: %s\n", ZSTD_getErrorName(result));
        if (bss_temp) free(bss_temp);
        return -1;
    }

    /* Reverse BSS if enabled */
    if (bss_temp) {
        sfa_bss_decode_frame(s, bss_temp, (uint8_t *)buf);
        free(bss_temp);
    }

    /* Verify checksum (on original, un-transformed data) */
    uint32_t check = sfa_crc32(buf, s->frame_bytes);
    if (check != entry.checksum) {
        fprintf(stderr, "SFA: checksum mismatch frame %u (got %08x, expected %08x)\n",
                frame_idx, check, entry.checksum);
        return -2;
    }
    return 0;
}

double sfa_frame_time(SFA *s, uint32_t frame_idx) {
    SFA_L2Entry entry;
    if (sfa_find_frame(s, frame_idx, &entry) < 0) return -1.0;
    return entry.time;
}

int sfa_read_kvmd(SFA *s, SFA_KVMDSet *sets, int max_sets) {
    if (!s || !s->fp) return 0;
    /* KVMD chunks live between CDEF end and JTOP start.
     * CDEF starts at s->cdef_offset. Its total size is 12 (chunk header) + payload.
     * CDEF payload = n_columns * 24. */
    uint64_t cdef_end = s->cdef_offset + 12 + (uint64_t)s->n_columns * 24;
    uint64_t jtop_start = s->first_jtop_offset;
    if (cdef_end >= jtop_start) return 0;  /* no room for KVMD */

    FILE *fp = s->fp;
    long saved_pos = ftell(fp);
    fseek(fp, (long)cdef_end, SEEK_SET);

    int n_found = 0;
    while ((uint64_t)ftell(fp) + 12 <= jtop_start && n_found < max_sets) {
        char type[4]; uint64_t size;
        if (fread(type, 1, 4, fp) != 4) break;
        if (fread(&size, 8, 1, fp) != 1) break;
        if (memcmp(type, "KVMD", 4) != 0) {
            fseek(fp, (long)size, SEEK_CUR);  /* skip unknown chunk */
            continue;
        }
        if (size < 16) { fseek(fp, (long)size, SEEK_CUR); continue; }

        /* Read sub-header */
        uint8_t version; uint16_t set_id, n_pairs;
        uint32_t first_frame, frame_count;
        uint8_t reserved[3];
        fread(&version, 1, 1, fp);
        fread(&set_id, 2, 1, fp);
        fread(&n_pairs, 2, 1, fp);
        fread(&first_frame, 4, 1, fp);
        fread(&frame_count, 4, 1, fp);
        fread(reserved, 1, 3, fp);

        SFA_KVMDSet *out = &sets[n_found];
        memset(out, 0, sizeof(*out));
        out->set_id = set_id;
        out->first_frame = first_frame;
        out->frame_count = frame_count;
        out->n_pairs = 0;

        /* Read pairs from remaining payload */
        uint64_t pair_bytes = size - 16;
        if (pair_bytes > 0 && pair_bytes < 1000000) {
            char *buf = (char*)malloc(pair_bytes + 1);
            fread(buf, 1, pair_bytes, fp);
            buf[pair_bytes] = '\0';

            uint64_t pos = 0;
            for (int i = 0; i < (int)n_pairs && pos < pair_bytes && out->n_pairs < 128; i++) {
                char *key = buf + pos;
                uint64_t kl = strlen(key); pos += kl + 1;
                if (pos >= pair_bytes) break;
                char *val = buf + pos;
                uint64_t vl = strlen(val); pos += vl + 1;
                strncpy(out->keys[out->n_pairs], key, 63);
                strncpy(out->values[out->n_pairs], val, 255);
                out->n_pairs++;
            }
            free(buf);
        }
        n_found++;
    }

    fseek(fp, saved_pos, SEEK_SET);
    return n_found;
}

int sfa_read_grids(SFA *s, SFA_GridDef *grids, int max_grids) {
    if (!s || !s->fp) return 0;
    uint64_t cdef_end = s->cdef_offset + 12 + (uint64_t)s->n_columns * 24;
    uint64_t jtop_start = s->first_jtop_offset;
    if (cdef_end >= jtop_start) return 0;

    FILE *fp = s->fp;
    long saved_pos = ftell(fp);
    fseek(fp, (long)cdef_end, SEEK_SET);

    int n_found = 0;
    while ((uint64_t)ftell(fp) + 12 <= jtop_start) {
        char type[4]; uint64_t size;
        if (fread(type, 1, 4, fp) != 4) break;
        if (fread(&size, 8, 1, fp) != 1) break;
        if (memcmp(type, "GDEF", 4) != 0) {
            fseek(fp, (long)size, SEEK_CUR);
            continue;
        }
        if (size < 8) { fseek(fp, (long)size, SEEK_CUR); continue; }
        uint8_t ver; uint16_t ng; uint8_t res[5];
        fread(&ver, 1, 1, fp); fread(&ng, 2, 1, fp); fread(res, 1, 5, fp);
        for (int i = 0; i < (int)ng && n_found < max_grids; i++) {
            SFA_GridDef *g = &grids[n_found];
            uint32_t pad;
            fread(&g->grid_id, 2, 1, fp); fread(&g->parent_id, 2, 1, fp);
            fread(&g->Nx, 4, 1, fp); fread(&g->Ny, 4, 1, fp); fread(&g->Nz, 4, 1, fp);
            fread(&g->Lx, 8, 1, fp); fread(&g->Ly, 8, 1, fp); fread(&g->Lz, 8, 1, fp);
            fread(&g->cx, 8, 1, fp); fread(&g->cy, 8, 1, fp); fread(&g->cz, 8, 1, fp);
            fread(&g->ghost, 2, 1, fp); fread(&g->kvmd_set, 2, 1, fp);
            fread(&pad, 4, 1, fp);
            n_found++;
        }
        break;  /* only one GDEF chunk expected */
    }

    fseek(fp, saved_pos, SEEK_SET);
    return n_found;
}

const char *sfa_dtype_name(uint8_t dtype) {
    static const char *names[] = {
        "f16","f32","f64","f128","i8","i16","i32","i64","u8","u16","u32","u64"
    };
    return (dtype < 12) ? names[dtype] : "unknown";
}

/* ---- Frame validity check ---- */

int sfa_count_valid_frames(SFA *s, uint64_t file_bytes) {
    if (!s || s->total_frames == 0 || !s->fp) return 0;
    FILE *fp = s->fp;
    long saved = ftell(fp);

    fseek(fp, (long)s->first_jtop_offset + 12 + 4 + 4 + 8, SEEK_SET);
    uint64_t jmpf_off;
    fread(&jmpf_off, 8, 1, fp);
    fseek(fp, (long)jmpf_off + 12 + 4 + 4, SEEK_SET);

    int valid = 0;
    for (uint32_t f = 0; f < s->total_frames; f++) {
        double ftime;
        uint64_t foffset, fcomp_size;
        uint32_t fcrc, freserved;
        fread(&ftime, 8, 1, fp);
        fread(&foffset, 8, 1, fp);
        fread(&fcomp_size, 8, 1, fp);
        fread(&fcrc, 4, 1, fp);
        fread(&freserved, 4, 1, fp);

        if (foffset == 0 && fcomp_size == 0 && f > 0) break;
        if (foffset + 12 + fcomp_size > file_bytes) break;
        valid++;
    }

    fseek(fp, saved, SEEK_SET);
    return valid;
}

/* ---- Streaming fixup: scan FRMD chunks, fill zeroed JMPF entries ---- */

int sfa_fixup_index(const char *path) {
    FILE *fp = fopen(path, "r+b");
    if (!fp) return -1;

    /* Read header basics */
    char type[4]; uint64_t size;
    fread(type, 1, 4, fp); fread(&size, 8, 1, fp);
    if (memcmp(type, "SFAH", 4) != 0) { fclose(fp); return -1; }

    uint32_t version, flags;
    fread(&version, 4, 1, fp);
    fread(&flags, 4, 1, fp);

    uint32_t Nx, Ny, Nz;
    fread(&Nx, 4, 1, fp); fread(&Ny, 4, 1, fp); fread(&Nz, 4, 1, fp);
    double Lx, Ly, Lz, dt_hdr;
    fread(&Lx, 8, 1, fp); fread(&Ly, 8, 1, fp); fread(&Lz, 8, 1, fp);
    fread(&dt_hdr, 8, 1, fp);
    uint32_t n_columns, total_frames_hdr;
    fread(&n_columns, 4, 1, fp);
    fread(&total_frames_hdr, 4, 1, fp);
    uint64_t first_jtop;
    fread(&first_jtop, 8, 1, fp);
    uint64_t cdef_offset;
    fread(&cdef_offset, 8, 1, fp);
    uint32_t jtop_max, jmpf_max;
    fread(&jtop_max, 4, 1, fp);
    fread(&jmpf_max, 4, 1, fp);

    /* Read CDEF to compute frame_bytes */
    fseek(fp, (long)cdef_offset, SEEK_SET);
    fread(type, 1, 4, fp); fread(&size, 8, 1, fp);
    uint64_t N_total = (uint64_t)Nx * Ny * Nz;
    uint64_t frame_bytes = 0;
    for (uint32_t c = 0; c < n_columns; c++) {
        uint8_t col_dtype;
        fseek(fp, (long)cdef_offset + 12 + c * 24 + 12, SEEK_SET);
        fread(&col_dtype, 1, 1, fp);
        frame_bytes += N_total * sfa_dtype_size[col_dtype];
    }

    int codec = flags & 0xF;

    /* Find the first JMPF offset from JTOP */
    fseek(fp, (long)first_jtop + 12 + 16, SEEK_SET);  /* first JTOP entry */
    uint64_t jmpf_offset;
    fread(&jmpf_offset, 8, 1, fp);

    /* Read existing JMPF entries to recover timestamps */
    double *saved_times = (double*)calloc(jmpf_max, sizeof(double));
    for (uint32_t i = 0; i < jmpf_max; i++) {
        fseek(fp, (long)jmpf_offset + 12 + 8 + i * 32, SEEK_SET);
        fread(&saved_times[i], 8, 1, fp);
    }

    /* Scan the file for FRMD chunks */
    fseek(fp, 0, SEEK_END);
    long file_size = ftell(fp);

    /* Find data start: after JTOP + JMPF */
    uint64_t jmpf_size_bytes = 12 + 8 + (uint64_t)jmpf_max * 32;
    uint64_t scan_start = jmpf_offset + jmpf_size_bytes;

    /* Allocate decompression buffer for CRC computation */
    uint8_t *decomp = (uint8_t*)malloc(frame_bytes);
    uint8_t *bss_temp = NULL;
    if (codec == SFA_CODEC_BSS) bss_temp = (uint8_t*)malloc(frame_bytes);

    uint32_t frame_count = 0;
    uint64_t pos = scan_start;

    while (pos + 12 <= (uint64_t)file_size) {
        fseek(fp, (long)pos, SEEK_SET);
        char chunk_type[4];
        uint64_t chunk_size;
        if (fread(chunk_type, 1, 4, fp) != 4) break;
        if (fread(&chunk_size, 8, 1, fp) != 1) break;
        if (chunk_size < 12 || chunk_size > (uint64_t)file_size) break;

        int is_frmd = (memcmp(chunk_type, "FRMD", 4) == 0);
        int is_frvd = (memcmp(chunk_type, "FRVD", 4) == 0);

        if (is_frmd || is_frvd) {
            /* FRMD: chunk_size = 12 + compressed_data_size (includes chunk header)
             * FRVD: chunk_size = time(8) + compressed_data_size */
            double embedded_time = 0;
            uint64_t comp_size;
            if (is_frvd) {
                fread(&embedded_time, 8, 1, fp);
                comp_size = chunk_size - 8;
            } else {
                comp_size = chunk_size - 12;
            }

            /* Read compressed data (skip decompression for FRVD — just index it) */
            void *comp = NULL;
            if (is_frmd) {
                comp = malloc(comp_size);
                fread(comp, 1, comp_size, fp);
            }

            uint32_t checksum = 0;
            if (is_frmd && comp) {
            /* Decompress to compute real CRC32 (FRMD only) */
            uint8_t *target = (codec == SFA_CODEC_BSS) ? bss_temp : decomp;
            size_t dec_size = ZSTD_decompress(target, frame_bytes, comp, comp_size);
            free(comp);

            if (!ZSTD_isError(dec_size)) {
                if (codec == SFA_CODEC_BSS && bss_temp) {
                    /* Reverse BSS to get original data for CRC */
                    /* Manual BSS decode — we don't have an SFA struct here */
                    /* For now, read column layout to do per-column BSS decode */
                    uint64_t off = 0;
                    for (uint32_t c = 0; c < n_columns; c++) {
                        uint8_t col_dtype;
                        fseek(fp, (long)cdef_offset + 12 + c * 24 + 12, SEEK_SET);
                        fread(&col_dtype, 1, 1, fp);
                        int es = sfa_dtype_size[col_dtype];
                        uint64_t col_bytes = N_total * es;
                        sfa_bss_decode(bss_temp + off, decomp + off, N_total, es);
                        off += col_bytes;
                    }
                }
                checksum = sfa_crc32(decomp, frame_bytes);
            }
            } /* end if (is_frmd && comp) */

            /* Determine which JMPF slot this frame goes in */
            uint32_t slot = frame_count % jmpf_max;

            /* Recover time: FRVD has embedded time, FRMD uses saved or computed */
            double time_val;
            if (is_frvd) {
                time_val = embedded_time;
            } else {
                time_val = saved_times[slot];
                if (time_val == 0.0 && frame_count > 0)
                    time_val = dt_hdr * frame_count;
            }

            /* Write the JMPF entry */
            long slot_off = (long)jmpf_offset + 12 + 8 + (long)slot * 32;
            fseek(fp, slot_off, SEEK_SET);
            fwrite(&time_val, 8, 1, fp);
            fwrite(&pos, 8, 1, fp);          /* FRMD chunk offset */
            fwrite(&comp_size, 8, 1, fp);    /* compressed size */
            fwrite(&checksum, 4, 1, fp);
            uint32_t ftype = SFA_FRAME_VOXEL;
            if (is_frvd) {
                /* Read existing frame_type from JMPF if slot was previously written.
                 * The vec writers store the correct I/P type when writing.
                 * Only fall back to heuristic if the slot is empty (0). */
                fseek(fp, slot_off + 28, SEEK_SET);  /* frame_type is at offset 28 in JMPF entry */
                uint32_t existing_ftype = 0;
                fread(&existing_ftype, 4, 1, fp);
                if (existing_ftype == SFA_FRAME_VEC_I || existing_ftype == SFA_FRAME_VEC_P) {
                    ftype = existing_ftype;
                } else {
                    ftype = SFA_FRAME_VEC_I;  /* default for uninitialized entries */
                }
            }
            fwrite(&ftype, 4, 1, fp);

            frame_count++;

            /* Handle JMPF overflow: would need to create new JMPF + update JTOP.
               For simplicity, cap at jmpf_max frames in fixup mode. */
            if (frame_count >= jmpf_max) {
                fprintf(stderr, "sfa_fixup: capped at %u frames (JMPF full)\n", jmpf_max);
                break;
            }
        }

        /* FRMD chunk_size includes the 12-byte header; FRVD chunk_size is payload only */
        if (is_frvd)
            pos += 12 + chunk_size;
        else
            pos += chunk_size;
    }

    free(decomp);
    if (bss_temp) free(bss_temp);
    free(saved_times);

    /* Update JMPF current_entries */
    fseek(fp, (long)jmpf_offset + 12 + 4, SEEK_SET);
    fwrite(&frame_count, 4, 1, fp);

    /* Update JTOP entry 0 frame_count */
    fseek(fp, (long)first_jtop + 12 + 16 + 12, SEEK_SET);  /* first entry + 12 bytes */
    fwrite(&frame_count, 4, 1, fp);

    /* Update SFAH total_frames and clear streaming flag */
    fseek(fp, 12 + 4 + 4 + 4*3 + 8*3 + 8 + 4, SEEK_SET);  /* offset 68 */
    fwrite(&frame_count, 4, 1, fp);

    uint32_t new_flags = flags & ~SFA_FLAG_STREAMING;
    fseek(fp, 12 + 4, SEEK_SET);
    fwrite(&new_flags, 4, 1, fp);

    fclose(fp);
    printf("sfa_fixup: indexed %u frames, cleared streaming flag\n", frame_count);
    return (int)frame_count;
}

#if defined(__GNUC__) || defined(__clang__)
#pragma GCC diagnostic pop
#endif

#endif /* SFA_IMPLEMENTATION */
