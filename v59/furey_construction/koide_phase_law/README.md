# koide_phase_law/ — study of  φ = Q/3

A focused study of the charged-lepton Brannen phase law **φ = Q/3** (lepton value
φ = 2/9), where Q = dimG₂/dimSpin7 = 2/3 is the Koide ratio.  Spun off because the
relation is simple, fundamental, and Koide-tight (10⁻⁵), yet has no derived mechanism.

## The claim in one line
The two dimensionless shape parameters of the charged-lepton mass triple are both
functions of the single structural number Q: `t² = (3Q−1)/2` (Koide) and `φ = Q/3`,
equivalently `9·arg ξ = 1 + 2|ξ|²` for the circulant order parameter ξ.

## Contents
- **`CONCEPT.md`** — theory document: the statement, the positive empirical status
  (as tight as Koide), the limits (what the law is NOT), the central question, and the
  "is it alone?" verdict.
- **`physical_mechanisms.md`** — candidate *physical* causes (CP flavon potential,
  geometric/Berry phase of the Z₃ cycle, silent-direction phase-fixing, hierarchy-from-
  phase, RG fixed point), with verdicts.  Recommended path: the J∘Z₃ geometric-phase
  computation (M2/M3).
- **`study.py`** — reproducible investigation (cross-sector, relation forms, precision,
  the electron near-boundary, the family hunt).  `python3 study.py`.
- **`FORWARD_PROPOSAL_v58_dynamics.md`** — brings the v58 unified field equation + its
  constraints into v59: the v58 Higgs-like vacuum manifold `M M̃ = v²` *is* the Koide
  surface `|ξ|²=1/2`; v59's grade-rigidity *derives* v58's hand-imposed scalar(real)/
  bivector(complex) split; and v58's bivector connection `ω ∝ ∇(MM̃)/|M|²` supplies the
  connection the geometric-phase mechanism (M2/M3) needs.  Recasts the next experiment as a
  holonomy calculation of the v58 connection around the generation Z₃ cycle.

## Machine-checked backing (in `../lean/`)
- `LeptonPhaseEmpirical.lean` — **positive**: φ=2/9 at Koide precision; "/3" = N_gen.
- `BrannenPhase.lean` — Q is phase-independent; phase enters only via cos 3φ; T1.2 null.
- `PhaseAmbiguity.lean` — masses fix only cos 3φ; φ=2/9 ≡ 2/9+2π/3; cos(2/3) ≠ cos²θ_W.
- `PhaseExclusions.lean` — φ not π-rational (non-geometric); gauge cBL=2 pinned.
- `TwoNinthsUnification.lean`, `WeinbergPatiSalam.lean` — the gauge-sector Q/3 = 2/9.

## Status
**Positive**: precise (10⁻⁵), lepton, divisor = generations, a genuine "second Koide
relation."  **Negative**: no mechanism; geometry / Weinberg-identity / naive potential
excluded; lepton-only.  **Witt map: DONE** (`witt_map.py`) — 3 disjoint γ-pairs + scalar `i` give 3 fermionic
modes with the exact Fock grading `(0,1,1,1,2,2,2,3)`; 8 states = `Λ•(ℂ³_color)` of one gen.
**Inter-generational connection: BUILT and tested** (`holonomy_result.md`) — the Brannen ring
has Wilson loop `arg(ξ³)=3φ` (a free Peierls/AB flux); **`arg(ξ³)=Q` is NOT forced** — no
structural curvature gives it (v58 connection is exact ⇒ 0; discrete/center fluxes π-rational;
associator integer; color `Z₃` fixes the lepton).  Clincher: `Q=2/3` is not π-rational
(`lean/PhaseExclusions.koide_not_pi_rational`), so it cannot be a geometric holonomy.
**⇒ the geometric-phase mechanism (M2/M3, the former lead) is RULED OUT.**

**Sedenion S₃ analysis** (`sedenion_S3_analysis.md`, `sedenion_s3.py`): per the literature the
generation symmetry lives in `Aut(𝕊)=G₂×S₃` (sedenions), not `Aut(𝕆)=G₂`.  Built `𝕊` and the
explicit `S₃` generators `ε,ψ`; **verified both are genuine automorphisms** (`ε²=ψ³=I`,
`εψ=ψ²ε`).  Result — the `φ=Q/3` mystery splits cleanly:
- **`/3` = `N_generations` = the sedenion `S₃` factor → EMERGENT/structural** (the `Z₃=⟨ψ⟩` is a
  real automorphism, ψ = a 120° rotation, the √3 = sin120°).  *Confirms the hypothesis.*
- **`Q = 2/3` (the phase magnitude `3φ`) = a coupling value → RESIDUAL**: ψ's intrinsic phases
  are `{0,±2π/3}` (π-rational), but `φ=2/9` is not π-rational, so the symmetry cannot fix it.

So `φ=Q/3` reduces to *(generations from S₃) × (a single non-symmetry magnitude Q)*.  The count
is structural; the magnitude is the residual, mechanism-less, Koide-class (`~10⁻⁵`) input whose
only structural target (transcendental `cos 2/3`) matches no candidate.
