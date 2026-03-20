# SFA — SCP Field Archive Format (v3)

A single-file, compressed, seekable, mmap-friendly, columnar container
for time-series volumetric field data. Chunk-based with two-level jump
tables for O(2) seek to any frame.

All multi-byte values are little-endian.

---

## File Structure

The file is a sequence of typed chunks. Every chunk begins with:

```
char[4]   type       four-character chunk type code
uint64    size       total size of this chunk in bytes (including these 12 bytes)
```

Walk the file: read 12 bytes, examine type, skip (size - 12) bytes to next chunk.
Unknown chunk types are safely skippable.

### Chunk Types

| Type | Name | Description |
|------|------|-------------|
| SFAH | File Header | Version, grid dims, column count, pointers |
| CDEF | Column Definitions | Schema: name, dtype, semantic per column |
| JTOP | Jump Table (top) | L1 index: points to JMPF chunks |
| JMPF | Jump Table (frame) | L2 index: points to FRMD chunks |
| FRMD | Frame Data | Compressed columnar field data for one timestep |

### Typical file layout

```
[SFAH]  File header (100 bytes)
[CDEF]  Column schema (12 + n_columns × 24 bytes)
[JTOP]  Top-level jump table (12 + 16 + max_l1 × 16 bytes)
[JMPF]  Frame jump table 0 (12 + 8 + max_l2 × 32 bytes)
[FRMD]  Frame 0
[FRMD]  Frame 1
...
[FRMD]  Frame 1023
[JMPF]  Frame jump table 1 (created when frame 1024 is written)
[FRMD]  Frame 1024
...
```

New JMPF chunks are appended lazily — only created when that frame
range is first written. If all JTOP entries fill (>1M frames), a new
JTOP is chained from the previous one.

---

## Chunk: SFAH — File Header

Size: 100 bytes total.

```
Offset  Type        Field
0       char[4]     type                "SFAH"
4       uint64      size                100
12      uint32      version             3
16      uint32      flags               reserved (0)
20      uint32      Nx                  grid points in x
24      uint32      Ny                  grid points in y
28      uint32      Nz                  grid points in z
32      float64     Lx                  domain half-width in x
40      float64     Ly                  domain half-width in y
48      float64     Lz                  domain half-width in z
56      float64     dt                  simulation time step
64      uint32      n_columns           column count (schema in CDEF)
68      uint32      total_frames        total frames written (updated on close)
72      uint64      first_jtop_offset   byte offset of the first JTOP chunk
80      uint64      cdef_offset         byte offset of the CDEF chunk
88      uint32      jtop_max_entries    L1 capacity (default 1024)
92      uint32      jmpf_max_entries    L2 capacity (default 1024)
96      uint32      reserved            0
```

Grid: dx = 2Lx/(Nx-1), N_total = Nx × Ny × Nz.

---

## Chunk: CDEF — Column Definitions

Size: 12 + n_columns × 24 bytes.

```
Offset  Type        Field
0       char[4]     type                "CDEF"
4       uint64      size                12 + n_columns × 24
```

Followed by n_columns entries, each 24 bytes:

```
Offset  Type        Field
0       char[12]    name                null-terminated (e.g. "phi_x")
12      uint8       dtype               data type code
13      uint8       semantic            semantic role code
14      uint8       component           axis (0=x, 1=y, 2=z)
15      uint8       flags               bit 0: optional, bit 1: derived
16      float64     scale               scale factor (default 1.0)
```

### dtype codes

| Code | Name | Bytes | Description |
|------|------|-------|-------------|
| 0 | f16 | 2 | IEEE 754 half |
| 1 | f32 | 4 | IEEE 754 single |
| 2 | f64 | 8 | IEEE 754 double |
| 3 | f128 | 16 | IEEE 754 quad |
| 4 | i8 | 1 | signed 8-bit |
| 5 | i16 | 2 | signed 16-bit |
| 6 | i32 | 4 | signed 32-bit |
| 7 | i64 | 8 | signed 64-bit |
| 8 | u8 | 1 | unsigned 8-bit |
| 9 | u16 | 2 | unsigned 16-bit |
| 10 | u32 | 4 | unsigned 32-bit |
| 11 | u64 | 8 | unsigned 64-bit |

Size lookup: `dtype_size[] = {2,4,8,16, 1,2,4,8, 1,2,4,8}`.

### semantic codes

| Code | Role | Description |
|------|------|-------------|
| 0 | position | displacement field φ_a |
| 1 | angle | rotation field θ_a |
| 2 | velocity | ∂φ/∂t or ∂θ/∂t |
| 3 | acceleration | ∂²/∂t² |
| 4 | energy | energy density |
| 5 | binding | triple product, |P| |
| 6 | torsion | antisymmetric gradient |
| 7 | metric | derived metric |
| 8 | mask | integer labels |
| 255 | custom | interpret by name |

---

## Chunk: JTOP — Top-Level Jump Table (L1)

Points to JMPF chunks. Size: 12 + 16 + max_l1 × 16 bytes.

```
Offset  Type        Field
0       char[4]     type                "JTOP"
4       uint64      size                12 + 16 + max_l1 × 16
```

### JTOP header (16 bytes)

```
12      uint32      max_entries         L1 capacity (matches SFAH.jtop_max_entries)
16      uint32      current_entries     L2 tables created so far
20      uint64      next_jtop_offset    byte offset of next JTOP (0 = none)
```

### JTOP entries (max_entries × 16 bytes each, starting at offset 28)

```
Offset  Type        Field
0       uint64      jmpf_offset         byte offset of a JMPF chunk
8       uint32      first_frame         global index of first frame in that JMPF
12      uint32      frame_count         frames stored in that JMPF
```

Unused entries (index >= current_entries): all zeros.

---

## Chunk: JMPF — Frame-Level Jump Table (L2)

Points to FRMD chunks. Size: 12 + 8 + max_l2 × 32 bytes.

```
Offset  Type        Field
0       char[4]     type                "JMPF"
4       uint64      size                12 + 8 + max_l2 × 32
```

### JMPF header (8 bytes)

```
12      uint32      max_entries         L2 capacity (matches SFAH.jmpf_max_entries)
16      uint32      current_entries     frames stored in this table
```

### JMPF entries (max_entries × 32 bytes each, starting at offset 20)

```
Offset  Type        Field
0       float64     time                simulation time of this frame
8       uint64      offset              byte offset of the FRMD chunk
16      uint64      compressed_size     compressed payload size within FRMD
24      uint32      checksum            CRC32 of uncompressed data
28      uint32      reserved            0
```

Unused entries: all zeros.

---

## Chunk: FRMD — Frame Data

Size: 12 + compressed payload.

```
Offset  Type        Field
0       char[4]     type                "FRMD"
4       uint64      size                12 + compressed_size
12      uint8[...]  payload             zstd-compressed column data
```

### Uncompressed payload layout

Columns concatenated in CDEF order. Each column is a contiguous
array of N_total values of that column's dtype:

```
[column 0: N_total × dtype_size[col0.dtype] bytes]
[column 1: N_total × dtype_size[col1.dtype] bytes]
...
[column n-1: N_total × dtype_size[col_{n-1}.dtype] bytes]
```

Total uncompressed size:
  frame_bytes = Σ (N_total × dtype_size[columns[c].dtype])

### Spatial indexing within a column

A column stores N_total = Nx × Ny × Nz values as a flat array.
The mapping from 3D grid coordinates to flat index:

```
flat_index = ix × (Ny × Nz) + iy × Nz + iz
```

- ix ∈ [0, Nx-1] — x index, varies SLOWEST
- iy ∈ [0, Ny-1] — y index
- iz ∈ [0, Nz-1] — z index, varies FASTEST

C row-major order. Traversal: iz increments first, then iy, then ix.

Physical coordinates:

```
x = -Lx + ix × 2Lx / (Nx - 1)
y = -Ly + iy × 2Ly / (Ny - 1)
z = -Lz + iz × 2Lz / (Nz - 1)
```

### Memory layout example (Nx=4, Ny=3, Nz=2)

```
Index:  0     1     2     3     4     5     6     7     ...
Grid: (0,0,0)(0,0,1)(0,1,0)(0,1,1)(0,2,0)(0,2,1)(1,0,0)(1,0,1) ...
        ---iz---  ---iz---  ---iz---  ---iz---  ---iz---
        --iy=0--  --iy=1--  --iy=2--  --iy=0--  --iy=1--
        ---------ix=0---------         ---------ix=1----------
```

An x-slice (fixed ix): contiguous block of Ny×Nz values.
A y-slice (fixed iy): Nz-sized blocks spaced Ny×Nz apart.
A z-slice (fixed iz): every Nz-th value.

---

## Seeking to Frame F

```
Step 1: Read SFAH.first_jtop_offset → seek to JTOP
Step 2: l1_index = F / jmpf_max_entries
         If l1_index >= JTOP.current_entries:
           follow JTOP.next_jtop_offset (rare, >1M frames)
Step 3: Read JTOP.entries[l1_index].jmpf_offset → seek to JMPF
Step 4: l2_index = F % jmpf_max_entries
         Read JMPF.entries[l2_index].offset → seek to FRMD
Step 5: Decompress FRMD payload
```

Two reads (JTOP + JMPF) for any frame in the first million.
With mmap: both reads are pointer dereferences, zero syscalls.

---

## Writing Sequence

```
1. Write SFAH (total_frames=0, offsets TBD)
2. Write CDEF
3. Write JTOP (empty, max_entries slots)
4. Write JMPF_0 (empty, max_entries slots)
5. Patch SFAH: first_jtop_offset, cdef_offset
6. Patch JTOP.entries[0]: point to JMPF_0

For each frame:
  7. Compress column data → write FRMD
  8. Update JMPF.entries[current]: time, offset, size, checksum
  9. Increment JMPF.current_entries
  10. If JMPF full:
        Write new JMPF_{n+1}
        Add entry in JTOP for new JMPF
        Increment JTOP.current_entries
  11. If JTOP full (>1M frames):
        Write new JTOP, chain from previous

On close:
  12. Update SFAH.total_frames
```

---

## Forward Compatibility

- Unknown chunk types: skip using size field
- Unknown dtype codes: skip column (size computable from code table)
- Unknown semantic codes: treat as custom (255)
- Version > 3: warn but read (chunk structure is self-describing)
- Optional columns (flags bit 0): safe to ignore

---

## Summary

| Property | Value |
|----------|-------|
| Seek cost | O(2) reads for <1M frames, O(3) for <1B |
| Chunk overhead | 12 bytes per chunk |
| L1 entry | 16 bytes (offset + range) |
| L2 entry | 32 bytes (time + offset + size + checksum) |
| Column def | 24 bytes per column |
| Default L1 capacity | 1024 → covers 1024 L2 tables |
| Default L2 capacity | 1024 → covers 1024 frames each |
| Max frames (no chain) | 1024 × 1024 = 1,048,576 |
| Compression | zstd |
| Grid | non-uniform Nx × Ny × Nz |
| Types | f16, f32, f64, f128, i8–i64, u8–u64 |
| Dependencies | libzstd |
