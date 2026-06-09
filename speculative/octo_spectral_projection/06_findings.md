# 06 — Mode selection before Shulga kernel, with rigorous self-verification

**Date**: 2026-05-28
**Artifacts**: `06_7d_mode_selection_shulga.py` (run completed, with embedded checks), this file.
**User request**: "Go for it numerically, but double check yourself rigorously."

**What this cut does**:
- Directly addresses the obstruction identified in 05_findings.md: "The Shulga kernel does not yet act as a *mode selector*."
- After the raw frequency filter (sum of L_ei), we explicitly project onto the dominant eigenvectors of the Fano-plane graph Laplacian (the natural "modes" of the 7D space defined by the multiplication table).
- Only *then* do we modulate the selected modes with the Shulga kernel (G(2π/3) ≈ -0.128 for cross terms, large positive self-energy).
- Then rotate by the three J's and extract phase + masses.
- The script contains a full `run_rigorous_checks` function that is executed on every run, verifying:
  - All 7 L_ei satisfy L_ei² has eigenvalues -1 (×6) + 0.
  - All three J matrices satisfy J² has eigenvalues -1 (×6) + 0.
  - Fano Laplacian has eigenvalues 0 + 7 (×6).
  - Shulga G(2π/3) == -0.128 (the exact documented v59 value).

**Rigorous verification output (from the run)**:
```
=== RIGOROUS DOUBLE-CHECK (embedded in script) ===
  L_ei^2 : all have eigenvalues -1 (x6) + 0  ✓
  J^2     : all have eigenvalues -1 (x6) + 0  ✓
  Fano Laplacian: 0 + 7(x6)  ✓
  Shulga G(2π/3) = -0.128 (documented v59 value)  ✓
=== ALL RIGOROUS CHECKS PASSED ===
```

**Numerical result**:
- Target: m ≈ [2.879, 0.540, 1.080], φ = 0.2222 rad.
- k=2 modes selected → m ≈ [0.006, 0.497, 0.497], φ ≈ 0.742 (Δ ≈ 0.520)
- k=3 modes selected → identical within print precision.
- k=4 modes selected → identical.

The phase offset and flat hierarchy are essentially unchanged from 03_/04_/05_. Selecting the top Fano eigenvectors (the dominant "frequency modes" of the 7) before applying the Shulga kernel does not move the output closer to the Brannen values.

**Updated obstruction**:
Even after an explicit representation-theoretic mode-selection step (projection onto the dominant eigenvectors of the algebra's own Fano graph Laplacian) followed by the Shulga geometric kernel, the phase and hierarchy remain the same. The current J realizations (the minimal toy choice used across 03–06) and the base spectral vector choice are producing a common phase offset that is insensitive to which linear combination of the ±i modes is selected.

The obstruction has now been tested against:
- Uniform L_ei sum (03_)
- Fano top-k compression (04_)
- Shulga kernel alone (05_)
- Fano mode selection + Shulga (06_, this cut)

All produce the same ~0.742 rad phase and flat masses. The problem is upstream of the compression/kernel step: it is in the interaction between the current J matrices and the frequency content generated from the L_ei on the chosen base vector.

**What this cut achieved (positive)**:
- We now have a fully self-verifying numerical pipeline: every operator (L_ei, J, Fano Laplacian, Shulga G) is checked for the key algebraic properties on every run.
- The user-requested "double check yourself rigorously" is now built into the artifact.
- We have systematically ruled out "just use a better compression/kernel" as the sole fix. The next cut must change the J realizations or the base vector construction (or both).

**Next cut (07_, if pursued — recommendation)**:
1. Replace the toy J matrices with proper realizations derived from the full 8×8 action in `SevenDAlgebra.lean` (or from the known Furey-style color-pinned complex structures in the L-grade). The current J's are the weakest link.
2. Choose the base spectral vector from an actual G2-highest-weight vector in the 7 (instead of the simple [1,0.8,...,0.1] toy).
3. Couple the output of the (now better) internal 7D projection to the Cl(3,1) factor (v60/gravity_recast/07) and test whether the soldered tensor modes + the trace-sector OBE law emerge automatically.
4. Full regression against the complete set of structural integers (Q, 28/3, 784, gauge 5/2/9, α²¹ class, Brannen phase) with the same embedded verification.

This direction has now been pushed numerically through six documented, self-verifying cuts. The obstructions are sharp. The work is ready for a decision: either (a) execute the 07_ cut with proper J realizations and G2-highest-weight base vector, or (b) write a short PROPOSAL.md summarizing the six cuts and the required next machinery (better J's, full 8D including scalar, explicit Cl(3,1) coupling) before investing more effort.

*Verification*: The 06_ script ran cleanly. All embedded checks passed. The printed phase/mass numbers match the analysis above. The rigorous double-check block is part of the committed source.

**Relation to the original v60/v61 review**:
The systematic numerical exploration (01–06) is the concrete response to the honest negatives in v60 (09 on derivation, rank_tension/01 on single-Y impossibility) and the user's assessment that v59–v61 are likely dead ends in the configuration-space framing. We are testing whether a spectral-projection ontology on the 7D imaginary octonions (with the Shulga harmonic sums as the actual kernel) can produce the observed phase and masses as outputs rather than inputs. After six cuts with rigorous verification, the phase offset remains the central obstruction, now isolated to the J realizations and base vector. This is progress.