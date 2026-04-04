# Code Duplication Report

Systematic audit of duplicated code in the SFA toolchain.
Each finding rates severity as HIGH (divergence risk, consolidate), MEDIUM (annoying
but low risk), or LOW (acceptable / intentional).

---

## 1. scp_sim.cu FrameWriter vs sfa.h sfa_write_frame_ex

**Files:** `sfa/sim/scp_sim.cu` (lines 1269-1461) vs `sfa/format/sfa.h` (lines 608-848)

**Finding: HIGH severity -- parallel reimplementation of SFA write path**

The CUDA kernel contains a complete FrameWriter subsystem (`fw_init`, `fw_destroy`,
`fw_write_raw`, `fw_write_voxel`, `fw_write_vec_iframe`, `fw_write_vec_pframe`) that
reimplements the compression and file-writing logic from `sfa_write_frame_ex` and
`sfa_write_vec_chunk` in sfa.h.

Specifically duplicated:
- **BSS+zstd per-column compression** (fw_write_voxel lines 1351-1376 vs
  sfa_write_frame_ex COLZSTD path lines 643-718). Nearly identical: column
  offset computation, BSS byte-shuffle loop, ZSTD_compress per column,
  assembly into [nc][sizes][data] payload.
- **JMPF index update** (fw_write_raw lines 1311-1329 vs sfa_write_frame_ex
  lines 821-841). Same fseek/fwrite pattern for JMPF entry, JMPF entry count,
  and JTOP frame count.
- **Vector frame compression** (fw_write_vec_iframe/pframe lines 1400-1461 vs
  sfa_write_vec_chunk lines 853-900). Same payload assembly + zstd compression.

**Why it exists:** The FrameWriter adds a pthread mutex so the async snapshot
writer thread and the main-thread vector hook can safely interleave writes to
the same SFA file. sfa.h's write functions are not thread-safe.

**What diverges:** fw_write_voxel always uses COLZSTD codec (hardcoded), while
sfa_write_frame_ex supports RAW, BSS, COLZSTD, and split output. fw_write_raw
does not handle JMPF overflow (new JMPF chunk creation), does not support
split output mode (.sfp files), and does not error-check ZSTD_compress return
values for individual columns. The index entry format also differs: fw_write_raw
writes `frame_type` in the last 4 bytes of the JMPF entry where
sfa_write_frame_ex writes `grid_id`.

**Recommendation:** Add a mutex/lock field to the SFA struct itself and make
`sfa_write_frame_ex` thread-safe internally (lock around the fseek/fwrite
section only, compress outside the lock). Then fw_write_voxel can be replaced
by `sfa_write_frame(s, time, cols)` directly. The JMPF overflow and split-mode
bugs in the CUDA path would be fixed for free.

---

## 2. scp_sim.cu vs scp_sim.c -- Host-side code

**Files:** `sfa/sim/scp_sim.cu` (2391 lines) vs `sfa/sim/scp_sim.c` (1317 lines)

**Finding: HIGH severity -- ~700 lines of near-identical host code**

The following are duplicated nearly verbatim (whitespace differences only):

| Section | .c lines | .cu lines | Identical? |
|---------|----------|-----------|------------|
| Config struct | 39-83 | 29-63 | ~95% (cu adds `output_split`) |
| cfg_defaults() | 85-108 | 65-88 | ~95% (cu adds `output_split`) |
| cfg_set() | 110-164 | 91-145 | ~95% (cu adds `output_split`) |
| cfg_load() | 166-191 | 148-170 | ~99% |
| cfg_print() | 193-217 | 172-196 | ~95% (banner text differs) |
| f64_to_f16() | 223-233 | 202-212 | identical |
| f16_to_f64() | 235-245 | 213-223 | identical |
| Grid struct | 251-261 | 311-319 | ~90% (cu omits mismatch/harden) |
| grid_alloc() | 263-291 | 321-338 | ~85% (cu omits mismatch alloc) |
| grid_save_pinned() | 303-314 | 1192-1203 | identical |
| init_oscillon() | 332-345 | 352-363 | identical logic |
| init_braid() | 347-369 | 365-385 | identical logic |
| init_from_sfa() | 371-437 | 387-424 | identical logic |
| init_template() | 461-553 | 429-553 | ~90% similar |
| do_init() | 555-561 | 555-561 | identical |
| compute_energy() | 851-930 | 1736-1802 | ~90% (cu lacks `#pragma omp`) |
| theta_rms() | 932-938 | 1804-1809 | ~95% (cu lacks omp) |
| P_integrated() | 940-945 | 1811-1815 | ~95% (cu lacks omp) |
| sfa_snap() | 953-978 | 1821-1849 | ~85% (cu adds gpu_download) |

Total: approximately 700 lines of host-only C code duplicated between the two files.

**Why it exists:** CUDA compilation requires `.cu` extension and nvcc. Sharing
host-side C code between a .c and .cu file is awkward -- you cannot `#include`
a .c file into a .cu file easily due to compiler differences.

**What diverges:** Several subtle differences have crept in:
- `compute_energy` in .cu does NOT use local `MODE` / `KG` variables; reads
  `c->mode` / `c->kappa_gamma` directly. Also computes curl once unconditionally
  rather than repeating it per-term (minor optimization difference).
- .cu `Grid` omits `mismatch[]` and `harden_Q[]` (GPU has its own device arrays).
- .cu `cfg_set` has `output_split` which .c doesn't have.
- .cu `init_from_exec()` is missing entirely.

**Recommendation:** Extract host-side code (Config, Grid, init functions,
compute_energy, f16 helpers) into `sfa/sim/scp_common.h` as a single-header
library. Both .c and .cu include it. Only the physics loop and GPU kernels
remain in the separate files. This would eliminate ~700 lines of duplication
and prevent divergence.

---

## 3. vecstream.h -- Partially redundant with SFA native FRVD chunks

**Files:** `sfa/format/vecstream.h` (1109 lines) vs SFA FRVD chunks in `sfa/format/sfa.h`

**Finding: MEDIUM severity -- two vector-frame formats coexist**

The CUDA kernel (`scp_sim.cu`) writes vector frames as native SFA FRVD chunks
using `fw_write_vec_iframe` / `fw_write_vec_pframe`. The CPU kernel (`scp_sim.c`)
writes vector frames to a separate `.vecstream` file using `vecstream.h`.

vecstream.h's `vecstream_fit_patch()` function (lines 777-900) contains a tricubic
patch fitting implementation that is also implemented in:
- `sfa/analysis/sfa_vecstream.c` fit_patch() (lines 69-160)
- `sfa/sim/vecstream_gpu_v2.cuh` (GPU kernel version)
- `sfa/format/vecstream.h` itself (lines 777-900)

The .vecstream file format is a separate binary format with its own magic bytes,
header, index, and frame structure. Now that SFA natively supports FRVD vector
chunks, the separate format is redundant for new work.

**Current users of vecstream.h:**
- `sfa/sim/scp_sim.c` -- CPU kernel (writes .vecstream files)
- `sfa/sim/vecstream_hook.h` -- dead code, not included anywhere
- `sfa/analysis/vecstream_verify.c` -- verification tool
- `sfa/tools/vecstream2sfa.c` -- converter
- `sfa/tools/sfa2vecstream.c` -- converter

**Recommendation:** Migrate `scp_sim.c` to write SFA FRVD chunks (like the CUDA
kernel already does), then deprecate vecstream.h. The converter tools and
verification tool become unnecessary. This would remove ~1100 lines of format
code and simplify the build.

---

## 4. vecstream_hook.h -- Dead code

**File:** `sfa/sim/vecstream_hook.h` (348 lines)

**Finding: LOW severity -- entirely dead code**

This header was designed as a portable vecstream hook for both CPU and CUDA
kernels. It is not `#include`d by any file in the codebase. Neither `scp_sim.c`
nor `scp_sim.cu` uses it:
- `scp_sim.c` has its own inline vecstream integration (using vecstream.h directly).
- `scp_sim.cu` uses `vecstream_gpu_v2.cuh` + its own VecHookCtx.

**Recommendation:** Delete. 348 lines of unused code.

---

## 5. vecstream_gpu.cuh (v1) -- Dead code

**File:** `sfa/sim/vecstream_gpu.cuh` (206 lines)

**Finding: LOW severity -- superseded by v2, entirely dead code**

Only referenced within its own comments (line 7: `#include "vecstream_gpu.cuh"`
in a usage example). `scp_sim.cu` includes `vecstream_gpu_v2.cuh` exclusively
(line 1856).

The v1 and v2 both contain the Vandermonde projection matrix computation
(Gauss-Jordan inversion of V^T V), which is the same algorithm in both.
v2 additionally computes a quadratic projection matrix and supports
multi-field patches.

**Recommendation:** Delete. 206 lines of unused code.

---

## 6. Analysis tools -- f16f() and col_f() copy-pasted

**Functions:** `f16f()` and `col_f()` appear identically in 7 analysis files.

**f16f() (f16-to-float conversion) -- 6 identical copies:**
| File | Line |
|------|------|
| sfa/analysis/sfa_vectorize.c | 67 |
| sfa/analysis/sfa_vectorize3d.c | 38 |
| sfa/analysis/sfa_vecstream.c | 38 |
| sfa/analysis/sfa_hybrid_vec.c | 26 |
| sfa/analysis/sfa_temporal_vec.c | 31 |
| sfa/analysis/equation_fit.c | 24 |

**col_f() (extract column value as float) -- 7 identical copies:**
| File | Line |
|------|------|
| sfa/analysis/sfa_vectorize.c | 75 |
| sfa/analysis/sfa_vectorize3d.c | 56 |
| sfa/analysis/sfa_vecstream.c | 45 |
| sfa/analysis/sfa_hybrid_vec.c | 33 |
| sfa/analysis/sfa_temporal_vec.c | 38 |
| sfa/analysis/equation_fit.c | 31 |
| sfa/tools/sfa2vecstream.c | 25 |

These are byte-for-byte identical across all files.

**Finding: HIGH severity for maintenance**

**Recommendation:** Add `sfa_f16_to_f32()`, `sfa_f16_to_f64()`, `sfa_f64_to_f16()`,
and `sfa_col_as_float()` to `sfa.h` in the public API section (before
`#ifdef SFA_IMPLEMENTATION`). This is the natural home since these functions
operate on SFA column data and dtypes. All 19 files with f16 helpers and all 7
with col_f could then use the sfa.h versions directly.

---

## 7. f16/f64 conversion -- 19 files, 3 variants

**Finding: MEDIUM severity -- three slightly different implementations**

There are three variants of f16 conversion scattered across the codebase:

**Variant A: `f16_to_f64` / `f64_to_f16` (double precision)**
Files: scp_sim.c, scp_sim.cu, scp_multi.c, scp_multi.cu, shell_analysis.c,
sfa_extract.c, sfa_render.c, freq_phase.c, accel_analysis.c, breathing.c,
gen_phase_confined.c, gen_proton_analytical.c (12 files)

**Variant B: `f16f` (float precision)**
Files: sfa_vectorize.c, sfa_vectorize3d.c, sfa_vecstream.c, sfa_hybrid_vec.c,
sfa_temporal_vec.c, equation_fit.c (6 files)

**Variant C: `sfa_f16_to_f64` (prefixed version)**
File: sfa/tools/sfa.c (1 file)

All three implement the same bit-manipulation but return different types (double
vs float) and have different edge-case handling (variant B returns +/-1e30f for
infinity, variants A and C return +/-INFINITY).

**Recommendation:** Same as item 6 -- consolidate into sfa.h. Provide both
`sfa_f16_to_f32()` and `sfa_f16_to_f64()` variants for callers that need
different precision.

---

## 8. Tricubic patch fitting -- 3 implementations

**Files:**
- `sfa/format/vecstream.h` vecstream_fit_patch() (lines 777-900)
- `sfa/analysis/sfa_vecstream.c` fit_patch() (lines 69-160)
- `sfa/sim/vecstream_gpu_v2.cuh` GPU kernel (lines 47-130)

**Finding: LOW severity -- different contexts require different implementations**

All three implement the same mathematical operation (separable tricubic
least-squares fit over a BS^3 block) but for different execution contexts:
- vecstream.h: generic C, used by the CPU kernel's vec output path
- sfa_vecstream.c: standalone analysis tool, reads SFA files
- vecstream_gpu_v2.cuh: CUDA kernel, runs on GPU with shared memory

The GPU version is necessarily different (CUDA intrinsics, shared memory,
thread-per-voxel). The two CPU versions are similar but not identical --
sfa_vecstream.c computes its own Vandermonde inverse, while vecstream.h
uses the same algorithm but with different variable names and error handling.

`sfa/analysis/sfa_vectorize3d.c` has a fourth variant (`fit_patch` at line 167)
that uses a general least-squares solver rather than the separable approach,
supporting variable-size blocks and adaptive refinement.

**Recommendation:** Acceptable duplication. The GPU version must be separate.
If vecstream.h is deprecated (see item 3), the duplication reduces to just
sfa_vecstream.c + the GPU kernel, which is unavoidable.

---

## Summary

| # | Item | Severity | Lines | Action |
|---|------|----------|-------|--------|
| 1 | FrameWriter in scp_sim.cu | HIGH | ~200 | Make sfa_write_frame_ex thread-safe |
| 2 | Host code in .cu vs .c | HIGH | ~700 | Extract scp_common.h |
| 3 | vecstream.h format | MEDIUM | ~1100 | Migrate CPU to FRVD, deprecate |
| 4 | vecstream_hook.h | LOW | 348 | Delete (dead code) |
| 5 | vecstream_gpu.cuh v1 | LOW | 206 | Delete (dead code) |
| 6 | f16f/col_f copy-paste | HIGH | ~20 per file x 7 | Add to sfa.h |
| 7 | f16 conversion variants | MEDIUM | ~15 per file x 19 | Consolidate in sfa.h |
| 8 | Tricubic fitting | LOW | n/a | Acceptable (GPU/CPU split) |

**Quick wins (no risk):** Delete items 4 and 5 (554 lines of dead code).

**High-value consolidation:** Items 6+7 (add f16 helpers to sfa.h) eliminates
~25 lines per file across 19 files and prevents the three variant flavors from
diverging further.

**Largest cleanup:** Item 2 (scp_common.h) would eliminate ~700 lines and
prevent the Config/init/diagnostics code from diverging between CPU and CUDA
builds. This is the highest-risk change but also the highest payoff.
