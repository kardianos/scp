# 01 — Toy frequency/phase projection on Z3: negative for this operator

**Date**: 2026-05-27 (speculative experiment after v59–v61 dead-end assessment)
**Artifacts**: `01_toy_frequency_projection.py` (runs clean), this file.
**Hypothesis under test**: The octonionic algebra (octo-space) is the primary information substrate; Brannen masses, Koide Q=2/3, and especially the Brannen phase φ=2/9 are outputs of a frequency/spectral transform + phase projection applied to algebraic information. The "before it comes real" layer contains the actual dynamics; 3+1 observables are the projected result.

**Test design (deliberately minimal to avoid circularity)**:
- Raw "octo-info": real 3-vector v on the sedenion S3 / triality Z3 (the generation 3-cycle). Input v is *not* pre-loaded with φ=2/9 or Brannen eigenvalue ratios.
- Frequency transform: classical DFT on Z3 (DC + two complex conjugate modes at ±1).
- Projection rule (the operator being tested): weighted sum of frequency components using only fixed v59 algebra invariants (dim G2=14, dim Spin(7)=21 → Q=14/21; universal deviation 28/3; selection dims 28/35/63). Extract effective phase φ_eff = arg(weighted sum); effective masses from |weighted frequency components|; compute Q_eff.
- Weight modes: "g2_content", "deviation", "selection", "none" (control).
- Success criterion: a natural (non-tuned) weighting produces φ_eff within ~0.05 rad of 2/9 *and* Q_eff within 1% of 2/3 from at least one non-special input v, without the projection rule containing 2/9 or the target mass ratios as free parameters.

**Result**: **Clean negative for this operator class.**

| weight_mode | φ_eff (typical) | Δ to 2/9 | Q_eff range | closest to 2/3? |
|-------------|------------------|----------|-------------|-----------------|
| g2_content  | 0 (or 2π)       | 0.222    | 0.36–0.85   | no              |
| deviation   | 0 (or 2π)       | 0.222    | 0.34–0.80   | no              |
| selection   | 0 or ~6.28      | 0.222–0.226 | 0.40–0.88 | no              |
| none        | 0 (or 2π)       | 0.222    | 0.34–0.79   | sometimes ~0.69–0.77 |

- No mode ever produces a phase near 2/9 from the test inputs.
- Algebra-derived weights do not help; in several cases they move Q_eff *away* from 2/3 relative to the unweighted control.
- The unweighted ("none") control sometimes gets Q closer to 2/3 than the weighted versions — the invariants are not doing useful spectral work in this setup.
- The DFT phase for any real input v on Z3 is either 0 (pure DC or symmetric) or determined by the input asymmetry; the real scalar weights cannot rotate it into the 2/9 offset.

**Why this formalization failed (precise obstruction)**:
1. **Phase source problem**. A classical DFT on a real 3-vector produces complex coefficients, but the argument of a *real-weighted* sum of those coefficients is fixed by the input v's asymmetry. Real weights (derived from dim ratios) scale amplitudes; they do not introduce an additional phase rotation. To get a preferred offset of 2/9 you must either (a) put it into the input, (b) make the weights complex, or (c) use a different "frequency operator" whose eigenmodes or characters already carry the desired phase.
2. **Scalar invariants are not enough**. dim G2, 28/3, selection dims are *numbers*. They can weight power in frequency bands, but the G7 radian-insert problem is exactly that a *number* (Q) must appear as the *argument* of a cosine (or as a phase in a complex exponential). A purely real weighting cannot transmute a scalar into an angle unless the operator itself has a built-in complex structure or character that supplies the rotation.
3. **Z3 DFT alone is too classical**. The sedenion S3 / triality action is already order-3, so the DFT is in some sense "already there." But the algebra's non-associativity, the L-grade complex structures J (J² = −1, three of them pinned by color), and the derivation algebra of the octonions are not used. A Fourier transform that ignores the multiplication law is probably the wrong transform for this algebra.

**What this does *not* prove**:
- It does not disprove the broader "octo-space is spectral, perception is projection" reframing. It only falsifies one simple concrete realization (classical DFT on the generation Z3 + real scalar weighting from dim ratios).
- It does not touch the richer structures (full L-grade 28 with its three J's, End(L) 784, the actual octonion multiplication table, derivation operators, or representation characters of G2/Spin(7)/Spin(8)).

**Adjacent prior work that already points in a better direction**:
- `brannen_phase_alpha.py` (root): explores φ_d ≈ 14·α(M_Z) and φ_u ≈ −10·α(0), with 14 = dim G2. This is a perturbative (loop) picture, but it already uses the G2 dimension as the structural coefficient in front of α. In a spectral view this could be re-read as "the +1 frequency mode couples with strength set by the G2 content."
- v59 G7 gap (`LeptonPhaseMagnitude.lean` lines 217–221): explicit TODO for a "spectral/character quantity (e.g. a character of the order-3 sedenion ψ weighted by the G₂-content Casimir) that PRODUCES cos(2/3) with the argument = Q." This is almost exactly the object a successful projection operator would need.
- `sfa/analysis/freq_phase.c`: phase extraction via `atan2(phi_1, phi_0)` at cluster centroids + autocorrelation for breathing frequencies. In the new ontology this tool is not post-processing; it is a model of the projection step itself.
- v2 hopfion/field-knot work: phase `arg(Φ)`, winding numbers, Hopf fibration, phase-dependent coupling `cos(φ₁ − φ₂)`. These already treat phase as generative rather than derived.

**What a second, non-toy iteration would require (minimal upgrades)**:
1. **Complex weighting / characters**. Replace real scalar weights with multiplication by a complex character or by one of the L-grade complex structures J (the three J's with J² = −1 that are already central to lepton = L forcing). The G2-content ratio could appear as the magnitude of the character; the J's would supply the phase rotation.
2. **Derivation operator as frequency generator**. Use the actual derivation algebra of the octonions (G2) or the associator to define a "frequency operator" whose eigenvalues or action on the even subalgebra replace the classical DFT. The eigen-"frequencies" would be algebraically natural rather than imposed by a Z3 cycle.
3. **Skewness / higher moments as the phase source** (per the G7 Lean TODO). Instead of arg of a weighted sum, extract the phase from a skewness measure or a character evaluated on the frequency spectrum. The target is then literally `cos(3φ) = cos(Q)` or a similar invariant, produced by the algebra without being inserted.
4. **Full L-grade or End(L) as the space**. The toy used only the generation Z3 factor. A real test must act on the 28-dimensional L-grade (or the 784-dimensional End(L)) and show that the projected spectrum reproduces the Brannen kernel on the grade-selected sector.
5. **Inverse map test**. Given an observed Brannen kernel (or an observed SFA field configuration), can a pre-image in the algebra + a projection operator reproduce it? This is the falsifiable direction for the "what we perceive is the projection" claim.

**Honest assessment**:
This experiment was a necessary first cut. It failed cleanly and for understandable structural reasons (real weights on a classical DFT cannot manufacture the required phase offset from scalar algebra invariants). The failure is useful: it tells us that any viable projection operator must use the *complex* or *non-associative* structure of the algebra (the J's, the derivation algebra, characters with phase) rather than dimension ratios as pure numbers.

The broader reframing (octo-space as spectral substrate, perception as projection) remains live and is arguably the most coherent way out of the G9 derivation obstruction, G1 rank tension, and G7 radian-insert problems that killed the configuration-space field theory program. But it will require new mathematical machinery (an algebraically natural Fourier or character transform on Cl(7)_even or its L-grade) rather than borrowing classical signal-processing tools and dressing them with dim ratios.

**Next step (if pursued)**: a bounded 2–3 week attack on item 1 or 2 above (complex J-weighted characters on the L-grade, or derivation operator as frequency generator), with the explicit success criterion being a machine-checkable statement that a natural algebraic character produces a phase whose cosine is within 1% of the observed Brannen phase or satisfies `cos(3φ) = Q` (or a documented clean negative that rules out this class).

**Relation to v60/v61 closeouts**: This direction is a direct response to the honest negatives in `v60/gravity_recast/09_obe_to_plebanski_findings.md` (Plebański parent is an independent posit) and `v60/gaps/rank_tension/01_findings.md` (single-Y is impossible; two-object reading is required). If the algebra is spectral and only its projection is "real," then asking the OBE to imply the full tensor theory in configuration space was always the wrong question.

*Verification*: `01_toy_frequency_projection.py` runs to completion with no errors; all printed numbers match the analysis above. No Lean module yet (toy level); the experiment is small enough that the Python source is the authoritative artifact.