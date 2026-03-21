# SFA — SCP Field Archive Format (v3.1)

A single-file, compressed, seekable, mmap-friendly, columnar container
for time-series volumetric field data. Chunk-based with two-level jump
tables, columnar schema, and multi-resolution AMR patch support.

All multi-byte values are little-endian.

---

## File Structure

The file is a sequence of typed chunks. Every chunk begins with:

```
char[4]   type       four-character chunk type code
uint64    size       total size of this chunk in bytes (including these 12 bytes)
```

Walk the file: read 12 bytes, examine type, skip (size - 12) to next chunk.
Unknown chunk types are safely skippable.

### Chunk Types

| Type | Name | Description |
|------|------|-------------|
| SFAH | File Header | Version, global metadata, pointers |
| GRID | Grid Definition | Named subgrid: dims, extent, center, parent |
| CDEF | Column Definitions | Schema: name, dtype, semantic per column |
| JTOP | Jump Table (top) | L1 index: points to JMPF chunks |
| JMPF | Jump Table (frame) | L2 index: points to FRMD chunks |
| FGRP | Frame Group | Groups FRMDs from one timestep across patches |
| FRMD | Frame Data | Compressed columnar field data for one grid, one timestep |
| META | Metadata | Key-value pairs (optional, user annotations) |

### Typical file layout (single grid, uniform)

```
[SFAH]  File header
[GRID]  Grid 0 (uniform, the whole domain)
[CDEF]  Column schema
[JTOP]  Top-level jump table
[JMPF]  Frame jump table 0
[FRMD]  Frame 0 (grid 0)
[FRMD]  Frame 1 (grid 0)
...
```

### Typical file layout (AMR, multi-resolution)

```
[SFAH]  File header
[GRID]  Grid 0: "global"   (N=256, L=100000, coarse)
[GRID]  Grid 1: "braid"    (N=128, L=5,      fine, parent=0, center=(0,0,0))
[GRID]  Grid 2: "electron" (N=128, L=5,      fine, parent=0, center=(0,0,50000))
[CDEF]  Column schema (shared across all grids)
[JTOP]  Top-level jump table
[JMPF]  Frame jump table 0
[FGRP]  Frame group t=0: {grid0: FRMD@off1, grid1: FRMD@off2, grid2: FRMD@off3}
[FRMD]  t=0, grid 0
[FRMD]  t=0, grid 1
[FRMD]  t=0, grid 2
[FGRP]  Frame group t=dt: {...}
[FRMD]  t=dt, grid 0
...
```

---

## Chunk: SFAH — File Header

Size: 112 bytes total.

```
Offset  Type        Field
0       char[4]     type                "SFAH"
4       uint64      size                112
12      uint32      version             3  (minor=1 implied by GRID chunk presence)
16      uint32      flags               bit 0-3: compression codec
                                          0=raw, 1=zstd, 2=bss+zstd (default)
                                          3=f32+bss+zstd (lossy, viz)
                                          4=f16+bss+zstd (lossy, preview)
                                          5=bq8+zstd (lossy, thumbnail)
                                          6=idelta+bss+zstd (lossless)
                                          7-15=reserved
                                        bit 4: AMR mode (1=multi-grid)
                                        bit 5: streaming (indexes may be incomplete)
                                        bit 6-31: reserved
20      uint32      n_grids             number of GRID definitions (1 for uniform)
24      uint32      n_columns           column count (schema in CDEF)
28      uint32      total_frames        total frame groups written
32      uint64      first_grid_offset   byte offset of first GRID chunk
40      uint64      cdef_offset         byte offset of CDEF chunk
48      uint64      first_jtop_offset   byte offset of first JTOP chunk
56      uint64      last_jtop_offset    byte offset of most recent JTOP
64      uint32      jtop_max_entries    L1 capacity (default 1024)
68      uint32      jmpf_max_entries    L2 capacity (default 1024)
72      float64     dt                  base simulation timestep
80      float64     time_start          simulation time of first frame
88      float64     time_end            simulation time of last frame (updated on close)
96      uint32      reserved[4]         zero
```

---

## Chunk: GRID — Grid Definition

Size: 12 + 112 bytes = 124 bytes.

```
Offset  Type        Field
0       char[4]     type                "GRID"
4       uint64      size                124

--- payload (112 bytes) ---

12      uint32      grid_id             sequential ID (0, 1, 2, ...)
16      char[16]    name                null-terminated (e.g. "braid", "global")
32      uint32      Nx                  grid points in x
36      uint32      Ny                  grid points in y
40      uint32      Nz                  grid points in z
44      uint32      parent_grid_id      ID of parent (coarse) grid, 0xFFFFFFFF if root
48      float64     Lx                  half-width in x
56      float64     Ly                  half-width in y
64      float64     Lz                  half-width in z
72      float64     cx                  center x coordinate (in parent's frame)
80      float64     cy                  center y coordinate
88      float64     cz                  center z coordinate
96      uint32      refine_ratio        refinement vs parent (1=same, 2=2× finer, etc.)
100     uint32      time_refine_ratio   temporal refinement (fine steps per coarse step)
104     float64     dt_local            this grid's timestep (may differ from base dt)
112     uint32      flags               bit 0: static (never changes between frames)
                                        bit 1: boundary_only (only store boundary data)
                                        bit 2-31: reserved
116     uint32      reserved            zero
120     uint32      reserved            zero
```

Grid spacing: dx = 2Lx/(Nx-1), dy = 2Ly/(Ny-1), dz = 2Lz/(Nz-1).
Total points: N_total = Nx × Ny × Nz.

Physical coordinates of grid point (ix, iy, iz):
```
x = cx - Lx + ix × 2Lx/(Nx-1)
y = cy - Ly + iy × 2Ly/(Ny-1)
z = cz - Lz + iz × 2Lz/(Nz-1)
```

For a single uniform grid: one GRID with parent=0xFFFFFFFF, center=(0,0,0).
For AMR: root grid (coarse) + child grids (fine) with parent references.

---

## Chunk: CDEF — Column Definitions

Size: 12 + n_columns × 24 bytes.

```
Offset  Type        Field
0       char[4]     type                "CDEF"
4       uint64      size                12 + n_columns × 24
```

Column entries (24 bytes each):

```
Offset  Type        Field
0       char[12]    name                null-terminated (e.g. "phi_x", "theta_z")
12      uint8       dtype               data type code
13      uint8       semantic            semantic role code
14      uint8       component           axis (0=x, 1=y, 2=z)
15      uint8       flags               bit 0: optional
                                        bit 1: derived
                                        bit 2: grid-specific (not all grids have it)
16      float64     scale               scale factor (default 1.0)
```

### dtype codes

| Code | Name | Bytes |
|------|------|-------|
| 0 | f16 | 2 |
| 1 | f32 | 4 |
| 2 | f64 | 8 |
| 3 | f128 | 16 |
| 4-7 | i8, i16, i32, i64 | 1, 2, 4, 8 |
| 8-11 | u8, u16, u32, u64 | 1, 2, 4, 8 |

### semantic codes

| Code | Role |
|------|------|
| 0 | position (φ) |
| 1 | angle (θ) |
| 2 | velocity |
| 3 | acceleration |
| 4 | energy |
| 5 | binding (|P|) |
| 6 | torsion |
| 7 | metric |
| 8 | mask/labels |
| 255 | custom |

---

## Chunk: FGRP — Frame Group

Groups multiple FRMDs from the same timestep across different grids.
Size: 12 + 16 + n_grids × 24 bytes.

```
Offset  Type        Field
0       char[4]     type                "FGRP"
4       uint64      size                12 + 16 + n_grids × 24

--- header (16 bytes) ---

12      float64     time                simulation time for this group
20      uint32      n_entries           number of grid FRMDs in this group
24      uint32      frame_index         global frame group index

--- entries (n_entries × 24 bytes each) ---

28      uint32      grid_id             which GRID this FRMD belongs to
32      uint64      frmd_offset         byte offset of the FRMD chunk
40      uint64      compressed_size     compressed payload size
48      uint32      checksum            CRC32 of uncompressed data
```

For single-grid files: FGRP has one entry. The JMPF points to the FGRP
(not directly to FRMD). The FGRP then points to its constituent FRMDs.

For AMR files: FGRP has n_grids entries, one FRMD per active grid.
Grids flagged as `static` may be omitted from later FGRPs (read from
the first FGRP that included them).

---

## Chunk: JTOP — Top-Level Jump Table (L1)

Points to JMPF chunks. Size: 12 + 16 + max_l1 × 16.

```
Offset  Type        Field
0       char[4]     type                "JTOP"
4       uint64      size

--- header (16 bytes) ---
12      uint32      max_entries
16      uint32      current_entries
20      uint64      next_jtop_offset    (0 = none)

--- entries (max_entries × 16 bytes) ---
28+     uint64      jmpf_offset
        uint32      first_frame
        uint32      frame_count
```

---

## Chunk: JMPF — Frame-Level Jump Table (L2)

Points to FGRP chunks (not directly to FRMD).
Size: 12 + 8 + max_l2 × 32.

```
Offset  Type        Field
0       char[4]     type                "JMPF"
4       uint64      size

--- header (8 bytes) ---
12      uint32      max_entries
16      uint32      current_entries

--- entries (max_entries × 32 bytes) ---
20+     float64     time
        uint64      fgrp_offset         points to FGRP chunk (not FRMD)
        uint64      reserved
        uint32      checksum
        uint32      reserved
```

---

## Chunk: FRMD — Frame Data

Compressed columnar data for ONE grid, ONE timestep.

```
Offset  Type        Field
0       char[4]     type                "FRMD"
4       uint64      size                12 + 4 + compressed_size
12      uint32      grid_id             which GRID this data is for
16      uint8[...]  payload             BSS+zstd compressed column data
```

### Uncompressed payload

Columns concatenated in CDEF order. Each column is a contiguous array
of N_total values of that column's dtype, where N_total = Nx×Ny×Nz
for the specific grid_id.

```
[column 0: N_total × dtype_size bytes]
[column 1: N_total × dtype_size bytes]
...
```

### Spatial indexing

```
flat_index = ix × (Ny × Nz) + iy × Nz + iz
x = cx - Lx + ix × 2Lx/(Nx-1)
y = cy - Ly + iy × 2Ly/(Ny-1)
z = cz - Lz + iz × 2Lz/(Nz-1)
```

ix varies slowest, iz varies fastest (C row-major).

### Compression

Determined by SFAH.flags bits 0-3. Available codecs:

- **codec 2 (default)**: BSS + zstd level 3. Per-column byte transpose
  groups exponent bytes together for better compression. ~1.5-5x lossless.
- **codec 3**: f32 downcast + BSS + zstd. Lossy (max err ~3e-8, 7 digits).
  ~5-18x. Suitable for visualization data.
- **codec 4**: f16 downcast + BSS + zstd. Lossy (max err ~2.4e-4, 3 digits).
  ~57-89x. Suitable for preview/thumbnail frames.
- **codec 5**: Block-quantize to 8-bit (256-value blocks with per-block
  scale+min as f32) + zstd. Lossy (max err ~3e-3). ~119-163x.
- **codec 6**: Int64 spatial delta + BSS + zstd. Lossless, +5-7% vs codec 2
  at higher zstd levels, but slower decode.

Lossy codecs (3-5) store reduced-precision data; the CDEF dtype remains
unchanged (typically f64). Readers reconstruct to the declared dtype.
See COMPRESSION_STUDY.md for detailed benchmarks.

---

## Chunk: META — Metadata (optional)

Key-value pairs for user annotations.

```
Offset  Type        Field
0       char[4]     type                "META"
4       uint64      size
12      uint32      n_entries
16+     For each entry:
          char[32]  key
          char[64]  value
```

Example: {"m2": "2.25", "mu": "-41.345", "eta": "0.5", "desc": "Cosserat braid T=100"}

---

## Seeking to Frame Group F

```
1. Read SFAH.first_jtop_offset → JTOP
2. l1_index = F / jmpf_max_entries → JTOP.entries[l1_index].jmpf_offset → JMPF
3. l2_index = F % jmpf_max_entries → JMPF.entries[l2_index].fgrp_offset → FGRP
4. Read FGRP → for each grid: FGRP.entries[g].frmd_offset → FRMD
5. Decompress FRMD payload for the desired grid(s)
```

Two reads for the jump tables + one FGRP read + one FRMD decompress per grid.

---

## Backward Compatibility

- v3.0 files (no GRID/FGRP chunks): reader creates an implicit single GRID
  from the SFAH dimensions and treats each FRMD as a single-entry FGRP.
- v3.1 files with n_grids=1: functionally identical to v3.0. The FGRP
  has one entry. Old readers that skip unknown chunks (GRID, FGRP) and
  read FRMD directly still work.

---

## Streaming Mode

When SFAH.flags bit 5 is set, the file is in streaming mode:
- JTOP and JMPF are pre-written with zeroed entry slots
- total_frames may not reflect actual frame count
- Readers should scan for non-zero JMPF entries

### Write protocol (producer)

1. Write SFAH+CDEF+JTOP+JMPF with zeroed entries (header block)
2. For each frame:
   a. Append FRMD chunk at end of file
   b. Seek to the pre-allocated JMPF slot:
      `slot_offset = jmpf_offset + 12 + 8 + (frame_idx % jmpf_max) * 32`
   c. Write 32 bytes: {time, frmd_offset, compressed_size, checksum, 0}
   d. Seek back to end for next FRMD
3. On completion: update total_frames at SFAH offset 68, clear bit 5

The JMPF slots are at KNOWN OFFSETS — no full rewrite needed. Each
frame write is: one append (FRMD) + one seek+overwrite (32-byte slot).

### Receive protocol (consumer, e.g., streaming from remote GPU)

1. Receive header block (SFAH+CDEF+JTOP+JMPF) → write to disk
2. Receive FRMD chunks → append to file
3. After each FRMD: seek to corresponding JMPF slot, write entry
4. The viewer can read frames as they arrive (check JMPF slot != 0)
5. On stream end: update total_frames, clear streaming flag

### Fixup pass (for incomplete files)

If the producer crashes or the stream is interrupted, the .sfa file
has FRMD chunks but zeroed JMPF entries. The fixup function:

```c
int sfa_fixup_index(const char *path);
```

Scans the file for FRMD chunks (identified by "FRMD" type tag),
fills in the zeroed JMPF slots, updates total_frames, and clears
the streaming flag. This is an O(n) scan with O(1) writes per frame.

---

## Summary

| Feature | Value |
|---------|-------|
| Seek cost | O(2) reads for <1M frames |
| AMR patches | unlimited via GRID definitions |
| Frames per timestep | n_grids FRMDs grouped in FGRP |
| Compression | 7 codecs: lossless (raw/zstd/bss+zstd/idelta+bss+zstd), lossy (f32/f16/bq8) |
| Column types | f16-f128, i8-u64 |
| Grid types | uniform, AMR child, static |
| Dependencies | libzstd |
