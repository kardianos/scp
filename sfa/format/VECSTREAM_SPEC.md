# VecStream Binary Format Specification

Version 1. Streaming polynomial-patch compression for 3D volumetric field data.

## Overview

VecStream compresses time-series volumetric data by representing each field as a
set of tensor-product polynomial patches over a regular block decomposition. Temporal
coherence is exploited with a video-codec-like I-frame / P-frame / K-frame scheme:

- **I-frame**: Full set of polynomial patches for one field at one time.
- **P-frame**: Sparse delta coefficients against the most recent I-frame or accumulated P-frames.
- **K-frame**: Raw voxel keyframe (zstd-compressed) for random access and error verification.

An optional **Fourier layer** stores per-coefficient temporal fits for ultra-compact
long-range interpolation.

## File Layout

```
[FileHeader]          64 bytes, fixed
[Frame 0]             I-frame or K-frame
[Frame 1]             P-frame, I-frame, or K-frame
...
[Frame N-1]           last frame
[FourierLayer]        optional, per-field
[FrameIndex]          array of index entries
[Footer]              16 bytes, fixed
```

All multi-byte values are little-endian. Offsets are absolute byte positions from
the start of the file.

## FileHeader (64 bytes)

| Offset | Size | Type     | Field       | Description                                    |
|--------|------|----------|-------------|------------------------------------------------|
| 0      | 4    | char[4]  | magic       | `"VECS"` (0x56, 0x45, 0x43, 0x53)             |
| 4      | 4    | uint32   | version     | Format version (currently 1)                   |
| 8      | 4    | uint32   | Nx          | Grid points in x                               |
| 12     | 4    | uint32   | Ny          | Grid points in y                               |
| 16     | 4    | uint32   | Nz          | Grid points in z                               |
| 20     | 8    | float64  | Lx          | Domain half-size in x                          |
| 28     | 8    | float64  | Ly          | Domain half-size in y                          |
| 36     | 8    | float64  | Lz          | Domain half-size in z                          |
| 44     | 8    | float64  | dt          | Time step between frames                       |
| 52     | 2    | uint16   | n_fields    | Number of field columns                        |
| 54     | 2    | uint16   | block_size  | Block edge length (e.g., 8)                    |
| 56     | 4    | uint32   | n_frames    | Total frames (updated at close)                |
| 60     | 4    | uint32   | flags       | Bit 0: has Fourier layer. Bits 1-31: reserved. |

## Frame

Each frame begins with a 32-byte header followed by a zstd-compressed payload.

### Frame Header (32 bytes)

| Offset | Size | Type     | Field             | Description                              |
|--------|------|----------|-------------------|------------------------------------------|
| 0      | 4    | char[4]  | magic             | `"FRVS"` (0x46, 0x52, 0x56, 0x53)       |
| 4      | 1    | uint8    | frame_type        | 0 = I-frame, 1 = P-frame, 2 = K-frame   |
| 5      | 1    | uint8    | field_idx         | Which field column this frame encodes    |
| 6      | 2    | uint16   | padding           | Reserved, must be 0                      |
| 8      | 8    | float64  | time              | Simulation time of this frame            |
| 16     | 8    | uint64   | compressed_size   | Bytes of payload after zstd compression  |
| 24     | 8    | uint64   | uncompressed_size | Bytes of payload before compression      |

### I-frame Payload (before compression)

```
n_patches:  uint32                    -- number of patches in this frame
[PatchHeader + coefficients] × n_patches
```

Each patch:

```
PatchHeader:  16 bytes (see below)
coefficients: n_coeffs × float32     -- tensor-product polynomial coefficients
```

### PatchHeader (16 bytes)

| Offset | Size | Type     | Field    | Description                                     |
|--------|------|----------|----------|-------------------------------------------------|
| 0      | 2    | int16    | origin_x | Block origin in grid coordinates (x)            |
| 2      | 2    | int16    | origin_y | Block origin in grid coordinates (y)            |
| 4      | 2    | int16    | origin_z | Block origin in grid coordinates (z)            |
| 6      | 1    | uint8    | size_x   | Block dimension in x (typically block_size)      |
| 7      | 1    | uint8    | size_y   | Block dimension in y                             |
| 8      | 1    | uint8    | size_z   | Block dimension in z                             |
| 9      | 1    | uint8    | order    | Polynomial order (3 = tricubic, 64 coefficients) |
| 10     | 2    | uint16   | n_coeffs | Number of coefficients: (order+1)^3              |
| 12     | 2    | float16  | max_err  | Max approximation error for this patch           |
| 14     | 2    | float16  | rms_err  | RMS approximation error for this patch           |

Coefficient ordering: tensor-product indices `[a][b][c]` where
`a` varies slowest (x), `c` varies fastest (z), matching
`coeffs[a * (O+1)^2 + b * (O+1) + c]` for polynomial order O.

### P-frame Payload (before compression)

```
n_deltas:      uint32                 -- number of patches with nonzero deltas
coeffs_per_patch: uint16              -- coefficients per patch (must match I-frame)
padding:       uint16                 -- reserved, 0
[patch_index(uint32) + delta_coeffs(coeffs_per_patch × float32)] × n_deltas
```

Only patches whose coefficient delta exceeds the encoder's tolerance are stored.
To reconstruct: accumulate deltas onto the most recent I-frame coefficients.

### K-frame Payload (before compression)

Raw voxel data for the entire field: `Nx × Ny × Nz × sizeof(float32)` bytes
before zstd compression. Stored as float32 regardless of original dtype.

K-frames serve as random-access sync points and verification references. A decoder
can start playback from any K-frame without needing prior frames.

## FourierLayer (optional)

Present when `flags & 1`. Located after the last frame, before the FrameIndex.

### FourierLayer Header (16 bytes)

| Offset | Size | Type     | Field           | Description                        |
|--------|------|----------|-----------------|------------------------------------|
| 0      | 4    | char[4]  | magic           | `"FLVS"` (0x46, 0x4C, 0x56, 0x53) |
| 4      | 1    | uint8    | field_idx       | Which field this layer covers      |
| 5      | 1    | uint8    | n_harmonics     | Fourier harmonics per coefficient  |
| 6      | 2    | uint16   | padding         | Reserved, 0                        |
| 8      | 4    | uint32   | n_patches       | Number of patches                  |
| 12     | 2    | uint16   | coeffs_per_patch| Coefficients per patch             |
| 14     | 2    | uint16   | padding2        | Reserved, 0                        |

### FourierLayer Data

For each patch, for each coefficient, for each harmonic:
```
amplitude:  float32    -- Fourier amplitude
omega:      float32    -- angular frequency
phase:      float32    -- phase offset
offset:     float32    -- DC offset (only for harmonic 0)
```

Total bytes: `n_patches × coeffs_per_patch × n_harmonics × 16`
(for harmonic 0 this includes the DC offset; for harmonics > 0
the offset field is zero).

The Fourier layer enables reconstruction at arbitrary times without
frame-by-frame playback: `c(t) = offset + sum_h amp_h * sin(omega_h * t + phase_h)`.

## FrameIndex

Located after the last frame (or FourierLayer if present).

### FrameIndex Header (16 bytes)

| Offset | Size | Type     | Field       | Description                          |
|--------|------|----------|-------------|--------------------------------------|
| 0      | 4    | char[4]  | magic       | `"IXVS"` (0x49, 0x58, 0x56, 0x53)   |
| 4      | 4    | uint32   | n_entries   | Number of index entries              |
| 8      | 8    | uint64   | reserved    | Must be 0                            |

### FrameIndex Entry (24 bytes)

| Offset | Size | Type    | Field      | Description                            |
|--------|------|---------|------------|----------------------------------------|
| 0      | 8    | float64 | time       | Simulation time                        |
| 8      | 8    | uint64  | offset     | Byte offset of frame header in file    |
| 16     | 1    | uint8   | frame_type | 0=I, 1=P, 2=K                         |
| 17     | 1    | uint8   | field_idx  | Which field                            |
| 18     | 6    | -       | padding    | Reserved, 0                            |

Total index size: `16 + n_entries × 24` bytes.

## Footer (16 bytes)

| Offset | Size | Type     | Field           | Description                        |
|--------|------|----------|-----------------|------------------------------------|
| 0      | 8    | uint64   | index_offset    | Byte offset of FrameIndex header   |
| 8      | 4    | uint32   | checksum        | CRC32 of FileHeader                |
| 12     | 4    | char[4]  | magic           | `"VSND"` (0x56, 0x53, 0x4E, 0x44) |

The footer enables backward seeking: read the last 16 bytes, verify the
`"VSND"` magic, then seek to `index_offset` to load the full frame index.

## Reconstruction Algorithm

To reconstruct field F at frame T:

1. Find the most recent I-frame or K-frame for field F at or before T.
2. If K-frame: decompress raw voxels directly.
3. If I-frame: decompress patches, then apply all P-frame deltas in order
   up to frame T.
4. Evaluate all patches onto the output voxel grid.

For random access, K-frames at regular intervals (e.g., every 100 frames)
bound the worst-case reconstruction cost.

## Compression Notes

- All frame payloads are zstd-compressed (level 3 default).
- I-frames: 64 f32 coefficients per 8^3 block = 256 bytes/block vs 2048 bytes raw.
  With zstd: typically 100-150 bytes/block (14-20x compression before temporal).
- P-frames: sparse, typically 5-30% of patches change per timestep.
  With zstd: 50-200 bytes per changed patch.
- K-frames: raw f32 voxels, zstd-compressed. Typically 2-4x compression.
- Fourier layer: 16 bytes per coefficient per harmonic. With 1-4 harmonics
  and 64 coefficients per patch, ~4 KB per patch for long-duration temporal fit.
