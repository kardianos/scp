# SFA Vectorization Design: Progressive Polynomial Compression

## Current State

`sfa_vectorize.c` works for 1D slices — compresses field data along
one axis into polynomial vectors (linear → quadratic → cubic) with
merge-on-tolerance. Achieves 5-8× compression at 0.01% error.

## Next Steps

### 1. Full 3D Vectorization

Extend from 1D slices to full 3D:
- Run along all three axes independently for each voxel row/column
- OR: use 3D tensor-product polynomials (patches of 2³, 4³, 8³ voxels)
- The tensor-product approach is more natural: fit a patch of voxels
  with a trilinear/triquadratic/tricubic polynomial

### 2. Reconstruction and Error Scoring

Write a reconstruction pass that:
1. Takes the polynomial vectors
2. Regenerates point data at the original voxel positions
3. Compares against the original data
4. Scores closeness: L∞ (max error), L² (RMS), L¹ (mean abs)
5. Iteratively adjusts tolerance to minimize error at target compression

### 3. SFA Output with Vector Overlays

Write the polynomial vectors as additional SFA columns:
- Each vector becomes a colored line segment in 3D space
- The line's color encodes the field term (R=P, G=φ, B=θ)
- The line's thickness/opacity encodes the polynomial order
- Volview renders these as line overlays on top of the volume

### 4. Iterative Error Minimization Loop

```
target_compression = 16×
tol = initial_guess
for iteration in 1..max_iter:
    vectors = vectorize(data, tol)
    recon = reconstruct(vectors)
    error = score(data, recon)
    compression = n_points / n_vectors
    
    if compression > target and error < target_error:
        tol *= 1.1  # relax tolerance to compress more
    elif compression < target:
        tol *= 0.9  # tighten tolerance
    elif error > target_error:
        tol *= 0.8  # tighten tolerance
    
    if converged: break
```
