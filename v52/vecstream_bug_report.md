# Vecstream V2 Heap Corruption Bug Report

## Summary

The `malloc(): corrupted top size` crash at t~100 (after ~50 P-frames) is caused
by a buffer overflow in the iframe path of `vec_hook()` in `scp_sim.cu`. The
`VecPatch` struct has a fixed-size `coeffs[64]` array, but the v2 multi-field code
writes 226 coefficients per patch, overflowing by 648 bytes per patch.

## Root Cause

**File**: `/home/d/code/scp/sfa/sim/scp_sim.cu`, line ~1820 (iframe path)

```c
vp->n_coeffs = vh->cpp;  // 226
memcpy(vp->coeffs, vh->h_coeffs + (long)p * vh->cpp, vh->cpp * sizeof(float));
//     ^^^^^^^^^^                                     ^^^^^^^^^^^^^^^^^^^^^^^
//     64 floats (256 bytes)                          226 floats (904 bytes)
```

The `VecPatch` struct is defined in `vecstream.h` (line 54-62):

```c
#define VS_NCOEFFS  64   /* (VS_MAX_ORDER+1)^3 = 4^3 */

typedef struct {
    int16_t  origin_x, origin_y, origin_z;
    uint8_t  size_x, size_y, size_z;
    uint8_t  order;
    uint16_t n_coeffs;
    uint16_t max_err_f16;
    uint16_t rms_err_f16;
    float    coeffs[VS_NCOEFFS];  // <-- ONLY 64 floats
} VecPatch;
```

The v2 code sets `vh->cpp = VS2_MULTI_TOTAL = 226` (64 phi_mag + 3*27 phi_corr +
3*27 theta_corr), but `VecPatch.coeffs` only holds 64 floats. The `memcpy` writes
226 floats (904 bytes) into a 256-byte buffer, overflowing 648 bytes.

## Why It Crashes After ~50 P-frames, Not Immediately

1. Frame 0 is an iframe. The code does `malloc(np * sizeof(VecPatch))` to allocate
   a contiguous array of `np` patches.

2. For each patch `p`, the `memcpy` overflows into the memory occupied by subsequent
   patches in the array. For patches 0 through np-2, this overwrites data that will
   subsequently be overwritten by later iterations anyway, so the corruption is
   "invisible" within the array body.

3. The **last patch** (p = np-1) overflows 648 bytes past the end of the malloc'd
   block, corrupting heap metadata (the "top chunk" size in glibc's malloc).

4. `vecstream_write_iframe` then reads `vp->n_coeffs` (226) from each patch and
   copies `n_coeffs * 4` bytes from `vp->coeffs`. For all patches except the last,
   this reads data that was overwritten by the NEXT patch's `memcpy`, producing
   garbled iframe data. For the last patch, the data is correct (it was the last
   to be written). But this read doesn't cause a crash because the overflowed
   region is still mapped memory.

5. `free(patches)` at line 1823 *might* crash, or the corrupted heap metadata
   might survive if glibc doesn't immediately check the top chunk. Subsequent
   `malloc`/`free` calls (from zstd compression in P-frame writes, etc.) eventually
   detect the corrupted top chunk size, triggering the crash ~50 frames later.

## Secondary Issue: Garbled I-Frame Data

Even if the crash is fixed, there is a data correctness bug in the iframe path.
Because each `VecPatch.coeffs` is only 64 floats, only the first 64 of 226
coefficients are stored correctly per patch. The `vecstream_write_iframe` function
serializes `n_coeffs` (226) floats from the 64-element array, reading 162 floats
of garbage per patch. Any reader that reconstructs fields from iframe data will
get corrupted results.

## Proposed Fix

There are two approaches:

### Option A: Bypass VecPatch for iframes (minimal change to vecstream.h)

In `vec_hook()`, replace the iframe path that uses `VecPatch` structs with a direct
call that builds the binary payload manually, similar to what P-frames already do.
Instead of:

```c
VecPatch *patches = (VecPatch*)malloc(np * sizeof(VecPatch));
for (int p = 0; p < np; p++) {
    // ... fill VecPatch ...
    memcpy(vp->coeffs, vh->h_coeffs + (long)p * vh->cpp, vh->cpp * sizeof(float));
}
vecstream_write_iframe(vh->vs, t, 0, patches, np);
free(patches);
```

Write a raw payload directly:

```c
// Build iframe payload manually with full 226 coefficients
size_t patch_header = 16;
size_t total = 4 + (size_t)np * (patch_header + vh->cpp * 4);
uint8_t *payload = (uint8_t*)malloc(total);
uint8_t *ptr = payload;
memcpy(ptr, &np_u32, 4); ptr += 4;
for (int p = 0; p < np; p++) {
    // write 16-byte patch header (origin, size, order, n_coeffs, errors)
    // write vh->cpp floats from vh->h_coeffs + p * vh->cpp
}
// compress and write as I-frame directly
```

### Option B: Increase VS_NCOEFFS in vecstream.h

Change `VS_NCOEFFS` from 64 to 256 (or make it dynamic). This is simpler but
changes the library ABI and wastes memory for users who only need 64 coefficients.
Given that VS2_MULTI_TOTAL is 226, a value of 256 would suffice.

**Recommendation**: Option A is better because it keeps vecstream.h as-is and
isolates the v2-specific layout to scp_sim.cu. The VecPatch struct was designed
for v1's single-field 4^3 = 64 coefficient patches and shouldn't need to change
for every new encoding scheme.

## Other Checks (No Additional Bugs Found)

### Atomic counter bounds (OK)
`vs2_temporal_residual` uses `atomicAdd(n_nonzero, 1)` with the result as a slot
index into `d_nz_patches[n_patches]`. Since there are exactly `n_patches` threads
and each can add at most once, the maximum slot value is `n_patches - 1`, which
is within bounds.

### Gather buffer bounds (OK)
`h_gather_buf` is allocated as `n_total = n_patches * cpp` floats. The gather loop
indexes `i * cpp` where `i < h_nnz <= n_patches`, so the max write position is
`(n_patches - 1) * cpp + cpp - 1 = n_total - 1`. Within bounds.

### CUDA async synchronization (OK)
The P-frame path has `cudaStreamSynchronize(vh->stream)` at line 1838 before
reading `h_nnz`, and another at line 1844 before the gather loop reads
`h_nz_patches` and `h_residual`. Both barriers are correct.

### Temporal refit division by zero (Low risk)
`vs2_temporal_refit` divides by `count`. The guard at line 1793 ensures
`temporal_count > 2` before calling refit, so `count` is always >= 3. Safe.

### vecstream_write_pframe internals (OK)
The function mallocs `total_size = 8 + n_deltas * (4 + cpp * 4)`, fills it
with a memcpy loop using proper `(size_t)i * coeffs_per_patch` offsets, compresses,
frees, and writes. No overflow here. The `coeffs_per_patch` parameter is `uint16_t`,
and 226 fits in uint16_t. No issues.

## If Static Analysis Is Insufficient

Run with compute-sanitizer (replaces cuda-memcheck in CUDA 11.6+):

```bash
compute-sanitizer --tool memcheck --leak-check full ./scp_sim_cuda config.cfg
```

For host-side heap checking specifically, since the bug is a host `memcpy` overflow:

```bash
# With AddressSanitizer (recompile with -fsanitize=address)
nvcc -O1 -g -Xcompiler -fsanitize=address -arch=sm_70 -o scp_sim_asan scp_sim.cu -lzstd -lm -lpthread
LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libasan.so.6 ./scp_sim_asan config.cfg

# Or with valgrind (slow but no recompile)
valgrind --tool=memcheck --track-origins=yes ./scp_sim_cuda config.cfg
```

ASan would immediately pinpoint the overflow at the exact `memcpy` call with a
stack trace showing the 648-byte out-of-bounds write.
