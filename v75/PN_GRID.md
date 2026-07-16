# Grid sizing for multi-particle packages (boundary safety)

**Date:** 2026-07-15

## Do we need a larger grid for Z6 / C₁₂?

**Yes for composites larger than F18 Z2.** Not because N (grid points) must grow with particle *count* alone, but because **package radius** (nuclear extent + L shell + Coulomb halo) must stay inside the **live domain** before the absorbing sponge.

## Domain geometry

| Symbol | Meaning |
|--------|---------|
| `L` | Box half-width (domain \(x\in[-L,L]\)) |
| `N` | Grid points per axis; \(dx = 2L/(N-1)\) |
| `damp_width` | Sponge thickness at outer radii |
| \(R_\mathrm{damp}=L-\mathrm{damp\_width}\) | Start of absorbing BC |

**Rule of thumb:**

\[
R_L + 2\,r_\mathrm{core} + \Delta_\mathrm{buf} \;\lesssim\; R_\mathrm{damp}
\]

with \(r_\mathrm{core}\sim 3\)–\(4\), \(\Delta_\mathrm{buf}\gtrsim 6\) (Coulomb/radiation margin).

## Comparison

| Package | N | L | \(dx\) | \(R_\mathrm{damp}\) | \(R_\mathrm{nuc}\) | \(R_L\) | L-margin | Status |
|---------|---|---|--------|---------------------|--------------------|---------|----------|--------|
| F18 Z2 | 128 | 32 | 0.50 | 28 | 5 | 18 | **10** | OK (tight but worked) |
| F16 Z6 (B1) | 192 | 48 | 0.50 | 44 | ~10 | 22 | **22** | PASS_park proven |
| **F19 Z6 (B2)** | **192** | **48** | **0.50** | **43** (dw=5) | **8** (p) | **22** | **~17** | **chosen** |
| F18 box + Z6 L@22 | 128 | 32 | 0.50 | 28 | 8 | 22 | **2** | **unsafe** (sponge eats L) |
| N=256 L=64 | 256 | 64 | 0.50 | 59 | 10 | 28 | 26 | safer, ~16GB tight on V100-16 |

## Boundary effects if undersized

1. **Sponge absorption of L** — massL/Ql secular loss, false “atom death”
2. **Nuclear charge drain** — Q_flux drop not from park but boundary
3. **False Coulomb** — image-like damping of E near sponge
4. **Fission misread** — fragments enter damp zone and vanish

## Memory (V100-16GB)

- Multi-fabric B2: N=192 is **proven** (F16 mz2)
- N=256 multi-fab is **risky** (~full 16GB with staging); use only if R_L must grow past ~28

## F19 freeze

```
N=192  L=48  damp_width=5
p octahedron R=8
n octahedron R=5.5 (inner)
L octahedron R=22  (count = Z)
dx ≈ 0.50  (same resolution as F18/F16)
```

Resolution is held at \(dx\sim 0.5\); **only box size increases** with package scale.
