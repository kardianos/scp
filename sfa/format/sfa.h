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
    uint32_t reserved;
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
} SFA;

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
void sfa_close(SFA *s);

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
        /* Per-column BSS+zstd: BSS-encode each column, then compress independently.
         * Layout: [n_cols(4)] [comp_size(8)]×n_cols [col0_data...] [col1_data...] ... */
        uint32_t nc = s->n_columns;
        uint64_t *col_comp_sizes = (uint64_t*)calloc(nc, sizeof(uint64_t));
        uint8_t **col_comp_data = (uint8_t**)calloc(nc, sizeof(uint8_t*));
        uint64_t total_comp = 4 + nc * 8;  /* header: n_cols + sizes array */

        uint64_t col_off = 0;
        for (uint32_t c = 0; c < nc; c++) {
            int es = sfa_dtype_size[s->columns[c].dtype];
            uint64_t col_bytes = s->N_total * es;

            /* BSS-encode this column: reorder bytes [v0b0,v0b1,...,v1b0,...] → [v0b0,v1b0,...,v0b1,v1b1,...] */
            uint8_t *bss_col = (uint8_t*)malloc(col_bytes);
            uint8_t *src_col = raw + col_off;
            for (int b = 0; b < es; b++) {
                for (uint64_t v = 0; v < s->N_total; v++) {
                    bss_col[b * s->N_total + v] = src_col[v * es + b];
                }
            }

            size_t bound = ZSTD_compressBound(col_bytes);
            col_comp_data[c] = (uint8_t*)malloc(bound);
            col_comp_sizes[c] = ZSTD_compress(col_comp_data[c], bound,
                                               bss_col, col_bytes, 3);
            free(bss_col);
            if (ZSTD_isError(col_comp_sizes[c])) {
                fprintf(stderr, "SFA: colzstd error col %u: %s\n", c,
                        ZSTD_getErrorName(col_comp_sizes[c]));
                for (uint32_t j = 0; j <= c; j++) free(col_comp_data[j]);
                free(col_comp_sizes); free(col_comp_data); free(raw);
                return -1;
            }
            total_comp += col_comp_sizes[c];
            col_off += col_bytes;
        }

        /* Assemble into a single buffer: [n_cols][sizes...][data...] */
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

void sfa_enable_split(SFA *s) {
    s->split_output = 1;
    /* Derive base path: strip .sfa extension from output path */
    const char *path = NULL;
    /* Get filename from fp — we need the output path stored somewhere.
     * The caller must set split_base before calling finalize_header.
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
    if ((s->flags & SFA_FLAG_STREAMING) && s->total_frames == 0) {
        uint32_t count = 0;
        uint64_t jtop_off = s->first_jtop_offset;
        while (jtop_off != 0) {
            fseek(fp, (long)jtop_off + 12, SEEK_SET);
            uint32_t max_ent, cur_ent;
            uint64_t next_jtop;
            fread(&max_ent, 4, 1, fp);
            fread(&cur_ent, 4, 1, fp);
            fread(&next_jtop, 8, 1, fp);
            for (uint32_t i = 0; i < cur_ent; i++) {
                fseek(fp, (long)jtop_off + 12 + 16 + i * 16 + 12, SEEK_SET);
                uint32_t fc;
                fread(&fc, 4, 1, fp);
                count += fc;
            }
            jtop_off = next_jtop;
        }
        if (count > 0)
            s->total_frames = count;
    }

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
            fread(&out->reserved, 4, 1, fp);
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

    int codec = s->flags & 0xF;

    /* Read compressed data */
    fseek(s->fp, (long)entry.offset + 12, SEEK_SET);  /* skip FRMD chunk header */
    void *comp = (void*)malloc(entry.compressed_size);
    fread(comp, 1, entry.compressed_size, s->fp);

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

        if (memcmp(chunk_type, "FRMD", 4) == 0) {
            uint64_t comp_size = chunk_size - 12;

            /* Read compressed data */
            void *comp = malloc(comp_size);
            fread(comp, 1, comp_size, fp);

            /* Decompress to compute real CRC32 */
            uint8_t *target = (codec == SFA_CODEC_BSS) ? bss_temp : decomp;
            size_t dec_size = ZSTD_decompress(target, frame_bytes, comp, comp_size);
            free(comp);

            uint32_t checksum = 0;
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

            /* Determine which JMPF slot this frame goes in */
            uint32_t slot = frame_count % jmpf_max;

            /* Recover time: use saved JMPF time if nonzero, else use dt*frame */
            double time_val = saved_times[slot];
            if (time_val == 0.0 && frame_count > 0)
                time_val = dt_hdr * frame_count;

            /* Write the JMPF entry */
            long slot_off = (long)jmpf_offset + 12 + 8 + (long)slot * 32;
            fseek(fp, slot_off, SEEK_SET);
            fwrite(&time_val, 8, 1, fp);
            fwrite(&pos, 8, 1, fp);          /* FRMD chunk offset */
            fwrite(&comp_size, 8, 1, fp);    /* compressed size */
            fwrite(&checksum, 4, 1, fp);
            uint32_t zero = 0;
            fwrite(&zero, 4, 1, fp);

            frame_count++;

            /* Handle JMPF overflow: would need to create new JMPF + update JTOP.
               For simplicity, cap at jmpf_max frames in fixup mode. */
            if (frame_count >= jmpf_max) {
                fprintf(stderr, "sfa_fixup: capped at %u frames (JMPF full)\n", jmpf_max);
                break;
            }
        }

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

#endif /* SFA_IMPLEMENTATION */
