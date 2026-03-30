# V50/View: SFA Preview Converter and Viewer Improvements

## Goal

Convert large SFA simulation files (100+ GB) into small SFA preview files
(~400 MB) that reuse the existing SFA container format. The preview files
can be loaded in volview for fast scrubbing through thousands of frames.

## Part 1: SFA Preview Converter (`sfa_preview`)

### What it does

Reads a full-resolution SFA file and outputs a smaller SFA file with:
- Fewer columns (3 uint8 derived channels instead of 12 float32 physics)
- Lower spatial resolution (downsampled to 64³ or configurable)
- Same frame count and timestamps
- All original KVMD parameters preserved
- Additional KVMD for reconstruction: quantization ranges, source N, etc.

### Column mapping

    Source (12 × f32):              Preview (3 × uint8):
    phi_x, phi_y, phi_z        →   P_abs    (uint8, |φ₀φ₁φ₂| quantized)
    theta_x, theta_y, theta_z  →   phi_sq   (uint8, |φ|² quantized)
    phi_vx, phi_vy, phi_vz     →   theta_sq (uint8, |θ|² quantized)
    theta_vx, theta_vy, theta_vz   (dropped)

### Quantization

For each channel, scan a sample of frames (first, middle, last) to find
the global min/max. Map [0, max] → [0, 255] linearly. Store the mapping
in KVMD:

    preview = true
    source_file = pp_D12_N100.sfa
    source_N = 100
    preview_N = 64
    quant_P_max = 0.15
    quant_phi_sq_max = 1.0
    quant_theta_sq_max = 0.05

This allows reconstruction: real_value = uint8_value / 255.0 * max.

### Downsampling

Simple box filter: average a 2³ (or NxNxN/preview_N³) block of source
voxels into one preview voxel. The derived quantities (|P|, |φ|², |θ|²)
are computed BEFORE downsampling to preserve nonlinear structure.

Order of operations per frame:
1. Read full-res frame (12 × f32 × N³)
2. Compute derived quantities at full resolution:
   - P_abs[i] = |phi0[i] * phi1[i] * phi2[i]|
   - phi_sq[i] = phi0[i]² + phi1[i]² + phi2[i]²
   - theta_sq[i] = theta0[i]² + theta1[i]² + theta2[i]²
3. Downsample each channel: box-average into preview_N³
4. Quantize to uint8 using global min/max
5. Write 3 × uint8 × preview_N³ to output SFA

### KVMD preservation

The output SFA copies ALL KVMD from the source (physics parameters,
grid info, etc.) as KVMD set 0. The preview-specific metadata goes in
KVMD set 1. This means any tool reading the preview can see the original
simulation parameters.

### Memory

At N=100: one source frame = 48 MB. Three derived channels at N=100 =
12 MB. Preview frame at N=64 = 786 KB. Total working memory ~60 MB.
Constant regardless of frame count — single-pass streaming.

### Estimated output sizes

    Source N=100, 4000 frames:
    Preview N=64, 3 × uint8: 786 KB/frame raw, ~200 KB compressed
    Total: ~800 MB raw, ~400 MB with COLZSTD

    Source N=384, 100 frames:
    Preview N=64, 3 × uint8: same per frame
    Total: ~20 MB

### Usage

    sfa_preview input.sfa output_preview.sfa [--preview_N 64] [--channels P,phi,theta]

### Build

    gcc -O3 -fopenmp -o sfa_preview sfa_preview.c -lzstd -lm

## Part 2: Viewer Improvements

### Fast scrubbing (volview changes)

The current volview viewer supports:
- Home/End: jump to first/last frame
- Left/Right arrow: frame ±1
- Mouse scroll: frame ±1

Add:
- **Shift+Left/Right**: jump ±10 frames
- **Ctrl+Left/Right**: jump ±100 frames
- **Page Up/Down**: jump ±100 frames
- **Number keys 1-9**: jump to 10%, 20%, ..., 90% of total frames
- **Scrub bar**: click anywhere on a horizontal bar at the bottom to
  jump to that frame proportion (mouse position / window width × total)
- **J/K keys**: variable-speed playback (J = faster, K = slower/reverse)

### Preview file detection

When volview opens an SFA with columns named `P_abs`, `phi_sq`, `theta_sq`
of type uint8, it should:
- Skip the derived quantity computation (data is already pre-computed)
- Apply the inverse quantization from KVMD (value = uint8 / 255 × max)
- Map directly to the volume texture

This requires minimal code changes: check column names/types in the
loading path, and branch to a "preview mode" that skips computation.

### Frame caching

For preview files (~200 KB/frame), cache the last N frames in RAM
(say N=100 = 20 MB). This makes backward scrubbing instant — no
re-decompression needed. For full SFA files this is impractical
(48 MB/frame) but for previews it's trivial.

## Part 3: Implementation Order

1. **sfa_preview converter** — standalone C tool. Read SFA, write smaller SFA.
   Test with the C4 162 GB file → target ~400 MB preview.

2. **Viewer scrub keys** — add Shift+arrow and PgUp/PgDn to volview.
   Quick Go code change, immediate usability improvement.

3. **Preview mode in viewer** — detect uint8 columns, skip computation.
   Moderate Go change, enables fast loading of preview files.

4. **Frame cache** — add LRU cache for recently decoded frames.
   Moderate Go change, makes backward scrubbing instant.

## Files

    v50/view/
      PLAN.md            — this file
      SVV_EXPLORATION.md — initial format exploration (from C4)
      sfa_preview.c      — the converter (to be written)
