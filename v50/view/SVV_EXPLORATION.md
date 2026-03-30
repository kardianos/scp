# SVV Format Exploration: Compressed Volumetric Viewer Data

Goal: read a 100+ GB SFA file and produce a 100-500 MB file that enables fast
frame scrubbing through thousands of frames in an interactive 3D viewer.

---

## 1. What the Viewer Actually Needs

The volview viewer (`sfa/volview/main.go`) reads 12 f32 columns per frame
(3 phi, 3 theta, 6 velocities) but immediately reduces them to 3-4 scalar
channels for the GPU ray marcher:

**Field view** (default):
- R = |P| = |phi0 * phi1 * phi2| (binding density)
- G = (phi0^2 + phi1^2 + phi2^2) * (1 - sqrt(|P|_norm)) (phi energy, suppressed where P is strong)
- B = theta0^2 + theta1^2 + theta2^2 (theta field)

**Velocity view**:
- R = sqrt(vphi0^2 + vphi1^2 + vphi2^2) (speed)
- G = |vphi0| (x-velocity)
- B = |vphi1| (y-velocity)

**Acceleration view**:
- R = cbrt(|P|) (cube-root binding density for dynamic range)
- G = theta_rms^2
- B = |v| (velocity magnitude)

Each is uploaded as an N^3 RGBA32F 3D texture (4 floats per voxel, alpha=1.0
always). The fragment shader does front-to-back ray marching with hardware
trilinear sampling and a transfer function (brightness, opacity, gamma).

**Key insight**: the viewer never accesses the raw phi/theta/velocity fields
independently after computing the derived channels. We can pre-compute the
derived scalars at conversion time and discard the raw fields entirely.

Reduction from pre-computation alone:
- Raw: 12 columns x 4 bytes = 48 bytes/voxel
- Derived: 3 channels x 4 bytes = 12 bytes/voxel (4x reduction)
- Or: 3 channels x 1 byte = 3 bytes/voxel if using uint8 (16x reduction)

The RGBA32F texture is already overkill -- the data goes through gamma
correction and gets clamped to [0,1] in the shader. The ray marcher uses
200 steps with trilinear sampling. An RGBA8 texture (uint8 per channel)
would be visually indistinguishable from RGBA32F for volume rendering at
this step count.

---

## 2. Spatial Downsampling

For interactive viewing, full simulation resolution is unnecessary. The
physically interesting features (braid cores, theta shells, binding halos)
span 5-20 grid cells. A 2x downsample loses fine detail but preserves all
gross structure.

### Size estimates for N=384 (production grid)

| Viewer resolution | Voxels | RGBA8 per frame | RGBA16 per frame |
|:-:|:-:|:-:|:-:|
| 384 (full) | 56.6M | 226 MB | 452 MB |
| 192 (2x down) | 7.1M | 28 MB | 57 MB |
| 128 (3x down) | 2.1M | 8.4 MB | 17 MB |
| 96 (4x down) | 885K | 3.5 MB | 7.1 MB |
| 64 (6x down) | 262K | 1.0 MB | 2.1 MB |

### Size estimates for N=128 (smaller runs)

| Viewer resolution | Voxels | RGBA8 per frame | RGBA16 per frame |
|:-:|:-:|:-:|:-:|
| 128 (full) | 2.1M | 8.4 MB | 17 MB |
| 64 (2x down) | 262K | 1.0 MB | 2.1 MB |

Downsampling method: box filter (average of 2x2x2 or 3x3x3 blocks). For
derived quantities that are inherently non-negative (|P|, theta^2, |v|),
averaging is fine. For quantities that could cancel (raw phi), we would
need RMS averaging, but since we pre-compute derived scalars first, simple
averaging works.

---

## 3. Temporal Compression

### The problem

The COMPRESSION_STUDY.md found that temporal differencing provides NO benefit
at the typical diagnostic cadence (dt=0.476, writing every 10 steps). But the
use case here is different: burst windows with dt=0.01 per frame and 1000+
frames. At this much finer temporal resolution, adjacent frames are extremely
similar.

### Approach A: Keyframe + delta encoding

Store a full keyframe every K frames. Between keyframes, store the difference
from the previous frame. Since derived quantities are uint8 or uint16, the
deltas are small integers that compress extremely well with zstd.

Expected delta characteristics for derived uint8 channels at dt=0.01:
- Most voxels change by 0 (background) or +/-1 (slow evolution)
- Only voxels near fast-moving features change by more
- Estimated: 90%+ of delta bytes are zero -> zstd ratio 10-50x on deltas

With K=20 keyframes:
- 1 keyframe every 20 frames at full size
- 19 delta frames at ~5-20% of keyframe size
- Effective per-frame cost: ~15-25% of keyframe

### Approach B: Adaptive frame skipping

Compute RMS difference between consecutive derived volumes. Only store a
frame when the RMS change exceeds a threshold. For slowly-evolving periods,
this can skip 10-100x frames. On playback, the viewer holds the last
stored frame until the next keyframe.

Problem: creates temporal jitter during fast events (which are exactly the
frames you want to see). Better suited as a complement to delta encoding,
not a replacement.

### Approach C: Fixed temporal subsampling

Store every Mth frame. For a burst window with 4000 frames at dt=0.01
(total T=40), storing every 4th frame gives 1000 frames at dt=0.04 -- still
visually smooth. Combined with spatial downsampling this may be sufficient.

### Recommendation

Keyframe + delta encoding (Approach A) gives the best compression with full
temporal fidelity. Adaptive skipping can be layered on top as a second pass.
For simplicity, fixed subsampling (every 2nd or 4th frame) is also effective.

---

## 4. Existing Volumetric Formats

### OpenVDB / NanoVDB

**What**: Industry-standard sparse volumetric data structure (VFX, games).
NanoVDB is a GPU-friendly linearized subset.

**Strengths**:
- Sparse representation -- empty regions cost nearly nothing
- GPU-native (NanoVDB runs on CUDA, OpenGL, WebGL, DirectX)
- Hardware-friendly memory layout (breadth-first tree, 32B-aligned nodes)
- Built-in block quantization (2/4/8/16-bit, 4-6x compression)
- Time-varying data support

**Weaknesses for our use case**:
- Our data is NOT sparse. SFA fields have non-zero background everywhere
  (A_bg=0.1 background field). The entire N^3 grid is populated.
- Heavy C++ dependency (OpenVDB requires TBB, Boost, blosc, etc.)
- NanoVDB is read-only; creating NanoVDB requires OpenVDB
- No native temporal delta compression
- Massively overengineered for our uniform-grid, dense data

**Verdict**: Poor fit. The sparse tree structure is wasted on dense data.
The library dependency chain is prohibitive.

### ZibraVDB / NeuralVDB

**What**: AI-based VDB compression. ZibraVDB claims 40-150x compression with
neural decompression on GPU. NeuralVDB (NVIDIA) uses neural networks for up
to 100x memory reduction.

**Strengths**: Extreme compression ratios for VFX data.

**Weaknesses**:
- Proprietary / commercial (ZibraVDB)
- Requires neural network inference for decompression (GPU-heavy)
- Designed for VFX smoke/cloud data, not physics fields
- Dependency on CUDA or specific GPU compute

**Verdict**: Interesting technology but wrong tool. We need simple, fast,
CPU-decodable data, not neural network inference.

### NRRD (Nearly Raw Raster Data)

**What**: Simple N-dimensional raster format with text header + raw data.
Supports gzip/bzip2 compression.

**Strengths**:
- Trivially simple format (text header + compressed blob)
- Supports any dimensionality (3D, 4D for time-varying)
- Wide tool support (Slicer, ParaView, ITK)
- Can use any element type (uint8, float16, float32, etc.)

**Weaknesses**:
- No per-frame random access (single compressed blob for entire volume)
- No temporal compression
- No streaming / progressive loading
- Compression limited to gzip/bzip2 (no zstd, no BSS)
- No index structure for seeking to frame N

**Verdict**: Too simple. Would work for single frames but fails the "fast
scrubbing through thousands of frames" requirement due to lack of random
access.

### KTX2 (Khronos Texture Container)

**What**: GPU texture container format with Basis Universal compression.
Supports 3D textures, mipmaps, cubemaps, texture arrays.

**Strengths**:
- Standard GPU texture format (Vulkan, OpenGL, WebGL)
- Basis Universal transcodes to any GPU format at runtime
- Supports 3D textures natively
- Browser-compatible (WebGL, WebGPU)
- Supercompression (zstd) in the container

**Weaknesses**:
- Basis Universal compression is designed for color textures (YCbCr color
  space, ETC1/UASTC block compression). Our data is 3-channel scientific
  scalar fields, not perceptual color.
- No time dimension -- KTX2 is a single-frame container. Time-varying data
  would need a sequence of KTX2 files or a custom wrapper.
- Block compression formats (ETC1S, UASTC, BC7) are designed for 2D texture
  blocks, not 3D voxel neighborhoods.
- Significant complexity for limited benefit over raw zstd.

**Verdict**: Possible for a WebGL viewer's last mile, but not suitable as
the primary archival/scrubbing format. The 2D texture compression codecs
are a poor fit for 3D scalar data.

### Raw volume sequences

**What**: Directory of raw binary files + JSON metadata.

**Strengths**:
- Dead simple. Each frame is one file.
- Instant random access (open file N, read it).
- Any viewer can load raw volumes.
- Easy to generate, easy to extend.

**Weaknesses**:
- Thousands of small files (filesystem overhead, slow to copy/archive)
- No temporal compression
- No single-file convenience
- Metadata separate from data

**Verdict**: Simple and functional but messy at 4000 frames. A single-file
format with an index is better.

---

## 5. Recommended Approach: Custom SVV Format

None of the existing formats satisfy all requirements:
1. Single file (not thousands of loose files)
2. Per-frame random access (seek to frame N instantly)
3. Pre-computed derived channels (not raw fields)
4. Spatial downsampling
5. Temporal compression (delta encoding)
6. Simple enough to implement in a weekend
7. zstd compression (already a dependency)

A custom format (SVV = SCP Volumetric Viewer) is the right call. It shares
SFA's chunk-based architecture but stores pre-computed visualization data
instead of raw physics fields.

### SVV Format Specification (Draft)

```
SVV file = HEADER + CHANNEL_DEFS + FRAME_INDEX + FRAME_DATA...

HEADER (fixed 64 bytes):
  magic:       "SVV1"           (4 bytes)
  version:     uint32           (1)
  flags:       uint32           (bit 0: has_delta, bit 1: has_mipmap)
  nx, ny, nz:  uint32 x 3      (viewer grid dimensions, after downsampling)
  n_channels:  uint32           (typically 3: |P|, phi2, theta2)
  n_frames:    uint32
  n_views:     uint32           (number of view modes stored, typically 1-3)
  dtype:       uint32           (0=uint8, 1=uint16, 2=float16)
  keyframe_interval: uint32     (K, for delta encoding; 0 = all keyframes)
  t_start:     float64          (simulation time of first frame)
  t_end:       float64          (simulation time of last frame)
  reserved:    8 bytes

CHANNEL_DEF (per channel, 32 bytes each):
  name:        char[16]         (e.g., "abs_P", "theta2", "phi2")
  semantic:    uint8            (0=binding, 1=torsion, 2=energy, 3=velocity, ...)
  norm_mode:   uint8            (0=per-frame auto, 1=global fixed, 2=log scale)
  scale:       float32          (for denormalization: physical = uint8 * scale)
  offset:      float32          (physical = uint8 * scale + offset)
  reserved:    6 bytes

VIEW_DEF (per view mode, 64 bytes each):
  name:        char[16]         (e.g., "field", "velocity", "accel")
  r_channel:   uint8            (index into channel array for red)
  g_channel:   uint8            (index into channel array for green)
  b_channel:   uint8            (index into channel array for blue)
  reserved:    45 bytes

FRAME_INDEX (array of n_frames entries, 32 bytes each):
  sim_time:    float64          (simulation time)
  file_offset: uint64           (byte offset of compressed frame data)
  comp_size:   uint64           (compressed size in bytes)
  flags:       uint32           (bit 0: is_keyframe, bit 1-2: view_id)
  checksum:    uint32           (CRC32 of decompressed data)

FRAME_DATA (per frame, variable size):
  -- If keyframe:
     zstd_compressed(channel0[nx*ny*nz] + channel1[...] + channel2[...] + ...)
  -- If delta frame:
     zstd_compressed(delta_channel0[...] + delta_channel1[...] + ...)
     where delta = (current - previous), stored as signed int8/int16
```

### Design decisions

**uint8 channels**: The viewer normalizes everything to [0,1] anyway, and the
RGBA32F texture is sampled at 200 ray steps with trilinear interpolation. At
8-bit precision, quantization noise is 1/256 = 0.4%, which is invisible in
volume rendering. This is the single biggest size win: 4 bytes -> 1 byte per
channel per voxel.

If more precision is wanted (e.g., for zoom-in on faint features), uint16
doubles the per-voxel cost but gives 1/65536 = 0.0015% quantization. This
is well below what the eye can distinguish even in high-dynamic-range modes.

**Per-frame normalization vs. global**: Per-frame auto-normalization
(current volview behavior) is simple but causes brightness flicker when
scrubbing. Global normalization (find max across ALL frames, normalize once)
gives stable brightness but may waste dynamic range on frames with outlier
values. The SVV format stores both the scale/offset (for global) and a
per-frame flag (for auto). The viewer can choose at runtime.

**Keyframe interval K=20**: At dt=0.01, a group of 20 frames spans dt=0.2,
over which the field changes very little. Delta frames in this regime should
compress 10-50x better than keyframes.

**View modes**: The converter pre-computes field view, velocity view, and
accel view channels. Storing all three as separate channel sets allows the
viewer to switch views without re-reading frames. This costs 3x the
channels but is still small.

**No octree/LOD**: For the target resolution (64-192 voxels per side), an
octree is overengineered. The entire volume fits in a single 3D texture and
can be uploaded in one call. If we ever need N=512+ viewer data, a simple
mipmap chain (full res + half res + quarter res) stored as separate frames
would be simpler than a full octree.

---

## 6. Size Estimates

### Use case: N=384 SFA, 4000 frames, burst window

**Source SFA (f32, 12 columns, colzstd)**:
- Uncompressed frame: 384^3 * 12 * 4 = 2.72 GB
- Compressed frame (~1.4:1): ~1.94 GB
- 4000 frames: ~7.7 TB uncompressed, ~5.5 TB compressed (massive)

This is the motivating problem.

**SVV at 128^3, uint8, 3 channels, field view only**:

| Component | Per frame | 4000 frames |
|:-:|:-:|:-:|
| Uncompressed keyframe | 6.3 MB (128^3 * 3) | -- |
| Compressed keyframe (~3:1 zstd) | 2.1 MB | -- |
| Compressed delta (~15:1 zstd) | 0.42 MB | -- |
| Total with K=20 (200 key + 3800 delta) | -- | **2.0 GB** |
| Total with K=20, every 2nd frame | -- | **1.0 GB** |
| Total with K=20, every 4th frame | -- | **500 MB** |

**SVV at 96^3, uint8, 3 channels, field view only**:

| Component | Per frame | 4000 frames |
|:-:|:-:|:-:|
| Uncompressed keyframe | 2.65 MB (96^3 * 3) | -- |
| Compressed keyframe (~3:1 zstd) | 0.88 MB | -- |
| Compressed delta (~15:1 zstd) | 0.18 MB | -- |
| Total with K=20 (200 key + 3800 delta) | -- | **850 MB** |
| Total with K=20, every 2nd frame | -- | **425 MB** |
| Total with K=20, every 4th frame | -- | **212 MB** |

**SVV at 64^3, uint8, 3 channels, field view only**:

| Component | Per frame | 4000 frames |
|:-:|:-:|:-:|
| Uncompressed keyframe | 786 KB (64^3 * 3) | -- |
| Compressed keyframe (~3:1 zstd) | 262 KB | -- |
| Compressed delta (~15:1 zstd) | 52 KB | -- |
| Total with K=20 (200 key + 3800 delta) | -- | **250 MB** |
| Total with K=20, every 4th frame | -- | **63 MB** |

**SVV at 64^3, uint8, 3 views x 3 channels = 9 channels total**:

| Component | Per frame | 4000 frames |
|:-:|:-:|:-:|
| Total with K=20 | -- | **750 MB** |
| Total with K=20, every 2nd frame | -- | **375 MB** |

### Sweet spot

For the target range of 100-500 MB with 4000 frames from a 384^3 simulation:

- **96^3, uint8, 3ch, K=20, skip 2**: ~425 MB (good quality, single view)
- **64^3, uint8, 9ch, K=20, skip 2**: ~375 MB (lower res, all 3 views)
- **64^3, uint8, 3ch, K=20, all frames**: ~250 MB (lower res, full temporal)
- **96^3, uint8, 3ch, K=20, skip 4**: ~212 MB (minimal, still usable)

Recommendation: **96^3, uint8, 3 channels, K=20** gives the best
quality/size tradeoff. At 96^3 the braid cores are clearly visible (each
core spans ~15-20 voxels at this resolution). Additional views can be stored
as separate channel groups, scaling linearly.

### Use case: N=128 SFA, 4000 frames

**Source**: 128^3 * 12 * 4 * 4000 = ~960 GB compressed at ~1.4:1 ratio.

**SVV at 64^3, uint8, 3ch, K=20**: ~250 MB (done).

---

## 7. WebGL Viewer Possibility

A WebGL/WebGPU viewer is feasible and would pair well with the SVV format.

### Architecture

1. Serve the SVV file over HTTP (static file or range requests)
2. Browser fetches the header + frame index (first ~200 KB)
3. On scrub, fetch the needed keyframe + delta chain (1-20 frame reads)
4. Decompress with wasm-zstd (existing WebAssembly zstd ports)
5. Reconstruct uint8 RGBA volume from channels
6. Upload as WebGL2 3D texture (RGBA8, supported since WebGL2)
7. Ray march with fragment shader (nearly identical to current volview GLSL)

### Size constraints

- WebGL2 3D texture max size: typically 2048^3 (driver-dependent)
- 96^3 RGBA8 texture: 3.5 MB GPU memory (trivial)
- 128^3 RGBA8 texture: 8.4 MB GPU memory (trivial)
- No GPU memory pressure at these sizes

### Bandwidth

At 64^3 with delta encoding, each frame fetch is ~50-250 KB. At 30 fps
scrubbing, that is 1.5-7.5 MB/s -- easily within localhost or LAN speeds.
Over the internet, pre-fetching a window of frames (e.g., 100 frames ahead)
into a ring buffer keeps latency hidden.

### Libraries

- Three.js has `DataTexture3D` for 3D volume textures (since r118)
- Raw WebGL2 `texImage3D` with `RGBA8` is straightforward
- The ray marching shader is ~50 lines of GLSL (already written in volview)
- wasm-zstd: multiple implementations exist (nicolo-ribaudo/pako-zstd, etc.)

### Feasibility

A minimal WebGL2 SVV viewer is roughly:
- 200 lines of JavaScript (fetch, decompress, texture upload, controls)
- 50 lines of GLSL (ray marcher, copied from volview with minor changes)
- 1 dependency (wasm-zstd)

This is a realistic weekend project. The SVV format is designed to make this
easy: the frame index enables seeking, uint8 channels map directly to RGBA8
textures, and zstd handles compression.

---

## 8. Conversion Pipeline

The `sfa2svv` converter would:

1. Open source SFA file
2. Scan all frames to find global normalization ranges (optional 2-pass)
3. For each frame:
   a. Read and decompress the 12 f32 columns
   b. Compute derived channels (|P|, phi^2, theta^2, |v|, etc.)
   c. Downsample to target resolution (box filter)
   d. Normalize to [0, 255] uint8 (or [0, 65535] uint16)
   e. If delta frame: subtract previous keyframe reconstruction
   f. Compress with zstd
   g. Write to SVV file, update frame index
4. Write frame index and header (can be at start with pre-allocated space,
   or at end with a pointer in the header)

Performance estimate: The bottleneck is SFA decompression + derived quantity
computation. At ~500 MB/s decompression and 384^3 frames, each frame takes
~5s to decompress + ~1s to compute + ~0.01s to downsample and compress.
4000 frames at ~6s each = ~7 hours single-threaded. With OpenMP and
pipelining (decompress frame N while computing frame N-1), closer to 2-3
hours.

For a 100 GB source SFA (e.g., N=128, 4000 frames at colzstd), conversion
would be much faster: ~25 MB/frame compressed, ~0.1s decompress each, total
~10 minutes.

---

## 9. Integration with Existing Volview

Two paths:

**Path A: SVV as a separate file, separate viewer**

Write a `sfa2svv` converter (C, ~500 lines) and a `svvview` viewer (modify
volview to accept SVV instead of SFA). The viewer becomes much simpler:
- No SFA parsing
- No derived quantity computation
- Just: read frame from index, decompress, upload RGBA8 texture
- Frame switching becomes near-instant (<10ms instead of 100ms+)

**Path B: SVV as a cache layer inside volview**

Volview generates the SVV cache on first open of a large SFA file (shows
progress bar), then uses the cache for scrubbing. The SVV cache lives next
to the SFA file as `filename.svv`. Subsequent opens skip conversion.

Path B is better UX but more complex. Path A is simpler to implement.

---

## 10. Summary and Recommendation

| Approach | File size (4000 frames, N=384 source) | Random access | Temporal compression | Complexity |
|:-:|:-:|:-:|:-:|:-:|
| Raw SFA (current) | 5.5 TB | Yes (slow) | No | 0 (exists) |
| SFA f16 transcode | ~100 GB | Yes (slow) | No | Low |
| NRRD 4D | ~750 MB | No | No | Low |
| NanoVDB sequence | ~2 GB | Yes | No | High |
| KTX2 sequence | ~1 GB | Per-file only | No | Medium |
| **SVV (uint8, 96^3, K=20)** | **~425 MB** | **Yes (fast)** | **Yes** | **Medium** |
| SVV (uint8, 64^3, K=20) | ~250 MB | Yes (fast) | Yes | Medium |

**Recommendation: Build the SVV format.**

It is the only approach that simultaneously achieves:
- 10,000x size reduction (5.5 TB -> 425 MB)
- Sub-millisecond frame seeking (byte offset from index)
- Sub-10ms frame loading (decompress 0.1-2 MB of zstd)
- Visually faithful rendering (3 derived channels at 96^3 uint8)
- Path to WebGL viewer (uint8 maps directly to RGBA8 textures)
- Simple implementation (no external library dependencies beyond zstd)

The format is intentionally simple: it is a flat array of zstd-compressed
uint8 volumes with a frame index for random access and optional delta
encoding for temporal compression. No trees, no neural networks, no GPU
compute for decompression.

### Implementation order

1. `sfa2svv` converter (C, reads SFA, writes SVV) -- the critical piece
2. Modify volview to load SVV files (or write standalone `svvview`)
3. Optional: WebGL viewer for browser-based scrubbing
4. Optional: 2-pass global normalization, multi-view support
