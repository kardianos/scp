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
 *    void *buf = malloc(s->frame_bytes);
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

/* Compression codec flags (in SFAH.flags bits 0-3) */
#define SFA_CODEC_RAW     0   /* no compression */
#define SFA_CODEC_ZSTD    1   /* zstd only */
#define SFA_CODEC_BSS     2   /* byte stream split + zstd (default, best for floats) */

/* Chunk type codes (as uint32 for comparison) */
#define SFA_CHUNK_SFAH  0x48414653  /* "SFAH" little-endian */
#define SFA_CHUNK_CDEF  0x46454443  /* "CDEF" */
#define SFA_CHUNK_JTOP  0x504F544A  /* "JTOP" */
#define SFA_CHUNK_JMPF  0x46504D4A  /* "JMPF" */
#define SFA_CHUNK_FRMD  0x444D5246  /* "FRMD" */

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

    /* Read state / mmap */
    void *mmap_ptr;
    uint64_t mmap_size;
} SFA;

/* ---- API ---- */

SFA *sfa_create(const char *path, uint32_t Nx, uint32_t Ny, uint32_t Nz,
                double Lx, double Ly, double Lz, double dt);
void sfa_add_column(SFA *s, const char *name, uint8_t dtype,
                    uint8_t semantic, uint8_t component);
int  sfa_finalize_header(SFA *s);
int  sfa_write_frame(SFA *s, double time, void **column_data);
void sfa_close(SFA *s);

SFA *sfa_open(const char *path);
int  sfa_read_frame(SFA *s, uint32_t frame_idx, void *buf);
double sfa_frame_time(SFA *s, uint32_t frame_idx);

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
    SFA *s = calloc(1, sizeof(SFA));
    s->fp = fopen(path, "wb");
    if (!s->fp) { free(s); return NULL; }
    s->mode = 1;
    s->version = SFA_VERSION;
    s->flags = SFA_CODEC_BSS;  /* default: BSS + zstd */
    s->Nx = Nx; s->Ny = Ny; s->Nz = Nz;
    s->Lx = Lx; s->Ly = Ly; s->Lz = Lz;
    s->dt = dt;
    s->jtop_max = SFA_JTOP_DEFAULT;
    s->jmpf_max = SFA_JMPF_DEFAULT;
    s->columns = calloc(64, sizeof(SFA_Column));
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
    /* Assemble uncompressed frame */
    uint8_t *raw = malloc(s->frame_bytes);
    uint64_t off = 0;
    for (uint32_t c = 0; c < s->n_columns; c++) {
        uint64_t col_bytes = s->N_total * sfa_dtype_size[s->columns[c].dtype];
        memcpy(raw + off, column_data[c], col_bytes);
        off += col_bytes;
    }

    /* CRC on original data (before any transform) */
    uint32_t checksum = sfa_crc32(raw, s->frame_bytes);

    /* Apply BSS transform if enabled */
    uint8_t *to_compress = raw;
    uint8_t *bss_buf = NULL;
    if ((s->flags & 0xF) == SFA_CODEC_BSS) {
        bss_buf = malloc(s->frame_bytes);
        sfa_bss_encode_frame(s, raw, bss_buf);
        to_compress = bss_buf;
    }

    /* Compress */
    size_t bound = ZSTD_compressBound(s->frame_bytes);
    void *comp = malloc(bound);
    size_t comp_size = ZSTD_compress(comp, bound, to_compress, s->frame_bytes, 3);
    free(raw);
    if (bss_buf) free(bss_buf);

    if (ZSTD_isError(comp_size)) {
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

    /* Write FRMD chunk */
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
    { uint32_t z = 0; fwrite(&z, 4, 1, s->fp); }

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
    return 0;
}

void sfa_close(SFA *s) {
    if (s->mode == 1 && s->fp) {
        /* Patch total_frames in SFAH */
        fseek(s->fp, 12 + 4 + 4 + 4*3 + 8*3 + 8 + 4, SEEK_SET);  /* offset 68 */
        fwrite(&s->total_frames, 4, 1, s->fp);
    }
    if (s->fp) fclose(s->fp);
    if (s->columns) free(s->columns);
    free(s);
}

/* ---- Reader ---- */

SFA *sfa_open(const char *path) {
    FILE *fp = fopen(path, "rb");
    if (!fp) return NULL;

    SFA *s = calloc(1, sizeof(SFA));
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
    s->columns = calloc(s->n_columns, sizeof(SFA_Column));
    for (uint32_t c = 0; c < s->n_columns; c++) {
        fread(s->columns[c].name, 1, 12, fp);
        fread(&s->columns[c].dtype, 1, 1, fp);
        fread(&s->columns[c].semantic, 1, 1, fp);
        fread(&s->columns[c].component, 1, 1, fp);
        fread(&s->columns[c].flags, 1, 1, fp);
        fread(&s->columns[c].scale, 8, 1, fp);
    }

    sfa_compute_frame_bytes(s);
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

    /* Read compressed data */
    fseek(s->fp, (long)entry.offset + 12, SEEK_SET);  /* skip FRMD chunk header */
    void *comp = malloc(entry.compressed_size);
    fread(comp, 1, entry.compressed_size, s->fp);

    /* Decompress (into buf or temp if BSS needs reversal) */
    uint8_t *decomp_buf = (uint8_t *)buf;
    uint8_t *bss_temp = NULL;
    if ((s->flags & 0xF) == SFA_CODEC_BSS) {
        bss_temp = malloc(s->frame_bytes);
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

const char *sfa_dtype_name(uint8_t dtype) {
    static const char *names[] = {
        "f16","f32","f64","f128","i8","i16","i32","i64","u8","u16","u32","u64"
    };
    return (dtype < 12) ? names[dtype] : "unknown";
}

#endif /* SFA_IMPLEMENTATION */
