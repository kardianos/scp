# 05 — Shulga S^7 harmonic kernel replaces Fano-graph compression

**Date**: 2026-05-28
**Artifacts**: `05_7d_shulga_kernel.py` (run completed), this file.
**User directive**: "Yes, replace compression with Shulga harmonic sums."

**What changed from 04_**:
- The "compression" step (information reduction / integration-out of the fast 7D modes) is no longer a projection onto the top-k eigenvectors of the combinatorial Fano-plane graph Laplacian.
- It is now the Shulga kernel G(θ) — the Green function of the Laplacian on S^7, computed via truncated Gegenbauer sums C_l^{(3)} in the v59 7D_Algebra / Shulga work.
- Documented key value used: G(2π/3) ≈ -0.128 (evaluated at the Z3 / triality family-shift angle).
- Self-energy (diagonal) is large positive, regularized by the same UV cutoff used in the v59 derivation.

**Result** (from the run):
- Target: m ≈ [2.879, 0.540, 1.080], φ = 0.2222 rad.
- Shulga kernel output: m ≈ [0.006, 0.497, 0.497], φ ≈ 0.742 rad (Δ ≈ 0.520 from target).

The phase and hierarchy are essentially the same as the Fano-compression result in 04_ (and the flat filter in 03_). The Shulga kernel modulates the overall scale (large positive self-energy + negative cross term) but does not preferentially select the linear combination of the ±i modes that would produce the Brannen offset 2/9 or the observed mass ratios after J-rotation.

**Precise obstruction (updated)**:
The Shulga kernel G(θ) at the documented angles is a global, angle-dependent scalar multiplier on the frequency content. It changes the *amplitude* of the self vs. mixing contributions, but the underlying "which ±i modes from the L_ei are active" is still set by the uniform sum of all seven left-multiplication operators. The kernel does not yet act as a *mode selector* in the representation-theoretic sense (G2-highest-weight projector or character filter on the 7). Until the kernel is applied after (or as part of) a G2-equivariant or Fano-eigenvector-weighted frequency filter, the phase offset remains the same "common" value produced by the current J realizations and base vector.

**What the Shulga kernel *did* contribute**:
- It is the first time the heavyweight, geometrically-derived 7D frequency analysis already present in the repo (the Laplacian Green function on S^7 used to derive λ/μ) was used as the central compression/integration operator in the projection.
- The negative sign at 2π/3 (Z3 shift) is correctly pulling the cross-generation mixing in the opposite direction from the self-energy, which is the correct qualitative behavior for a "family-shift" interaction in the algebra.
- The magnitude of the effect is now tied to a documented internal-geometry calculation rather than an arbitrary graph Laplacian.

**Honest assessment**:
Replacing the toy Fano compression with the Shulga kernel is the right conceptual move (it brings in the actual 7D harmonic analysis the repo already contains). However, at the current level of wiring (uniform L_ei sum + Shulga modulation + J rotation), it does not yet solve the phase-offset or hierarchy problem. The kernel is doing amplitude modulation, not mode selection.

**Next cut (06_, if pursued)**:
1. Apply a G2-highest-weight or Fano-eigenvector filter *first* (select the right combination of ±i modes using the actual automorphism group of the 7), *then* modulate the selected content with the Shulga kernel.
2. Couple the compressed + projected internal 7D layer to the Cl(3,1) spacetime factor (v60/gravity_recast/07) and test whether the surviving long-range tensor modes (the soldered 2-form whose trace recovers the old OBE) emerge naturally.
3. Full regression against all structural integers (Q, 28/3 deviation, 784, gauge 5/2/9, α²¹ magnitude class) rather than just the Brannen kernel.

This cut demonstrates that the Shulga harmonic sums can be wired into the projection operator without breaking the existing test harness. That is incremental progress on "replace compression with Shulga."

*Verification*: The 05_ script ran cleanly; the printed result matches the numbers above. The kernel uses the exact documented G(2π/3) value from the v59 Shulga derivation.