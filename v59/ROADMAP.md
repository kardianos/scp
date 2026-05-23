# v59 ROADMAP — Open Questions and Next-Session Targets

**Date**: 2026-05-22
**Status**: Consolidated open-questions list with concrete next-session priorities.
**Parents**: [`SESSION_2026-05-22.md`](SESSION_2026-05-22.md), [`SUMMARY.md`](SUMMARY.md)

This document is the master "where do we go next" list for v59.  It
consolidates every open question across the project as of 2026-05-22.

---

## Big picture: the four open frontiers

```
   1. DYNAMICAL ξ(x) field        ← active investigation, this session
   2. MECHANISM for selection rule (Z₂×Z₂ pattern)
   3. LAGRANGIAN derivation of (5, 21/16) prefactors
   4. SCALE bridges (Brannen a ↔ Higgs VEV ↔ Planck mass)  [PARTIAL CLOSE]
```

Each is a multi-session research effort.  None of them is "do once and
finish".  But each has concrete entry points.

**Status update (later in 2026-05-22 session)**: Frontier 4 / Q4-1
(Brannen-to-Higgs bridge) is largely **closed**. The relation
`v_Higgs = D_lepton² · a_l² = 28² · a_l²` matches at 0.068 %, cascades
to m_W (0.04 %), m_Z (0.02 %) when combined with g_W² = 5√α and the
new `sin²θ_W = 2/9`, `cos²θ_W = 7/9`. The SM tree identity then forces
**α(M_Z) = 25/(324π²)** at 0.03 % match — making α structural at BOTH
ends of the running scale (q=0 from -ln α + 2α = π²/2; q=M_Z from
this session). See `synthesis/FINDINGS_scale_bridge.md`.

The EW sector now has zero empirical α inputs; only a_l is empirical.

Remaining open in Frontier 4: Higgs mass mechanism (m_H ≈ √(7/27)·v
at 0.14 %, no mechanism); Planck-mass bridge; cross-sector Brannen
scales (a_u² ≈ 35·a_d² hint, within quark-mass uncertainty).

---

## Frontier 1: Dynamical ξ(x) field

**Question**: How does the v59 Brannen-kernel parameter ξ ∈ ℍ become a
dynamical field ξ(x) on spacetime?

### Status (as of 2026-05-22 session)

We set up the simplest Lagrangian:
```
   L_ξ  =  (1/2) (∂_μ ξ̄)(∂^μ ξ)  −  (λ/4)(|ξ|² − 1/2)²
   EOM  =  □ ξ_a  =  -λ (|ξ|² − 1/2) · ξ_a
```
Maxima-verified Hessian at lepton vacuum (1/√2, 0, 0, 0):
- 3 massless Goldstone modes (S³ tangent)
- 1 massive radial mode with mass² = λ

**Qualitative SM Higgs-sector match**: 3 Goldstones could be eaten by
gauging the silent SU(2)_L → W±, Z⁰ longitudinals.  Radial mode could
be the Higgs.

### Open questions

**Q1-1** (Sector multi-vacua): How does ONE field accommodate THREE sector
equilibria (1/2, 3/5, 7/9)?

Four candidate options (`synthesis/DYNAMIC_XI.md` § 3):
- **Option A**: Sector-specific `ξ_X(x)` fields (one per fermion type).
  Simplest but multiplies field content.
- **Option B**: Single ξ(x) with multi-minimum potential
  `V = λ · ∏_X (|ξ|² − r_X²)²`.
- **Option C**: Single ξ(x), sector-dependent effective potential
  (matter-content-modulated).
- **Option D**: ξ is the LEPTON field; quarks are composite bound states
  exploiting the additive identity `D_u = D_e + D_d`.

**Q1-2** (Coupling): What sets the scalar coupling λ?
- λ = π²/2 (= v59 S_em)?  Mass would be √(π²/2) ≈ 2.22.
- λ = 1?
- λ = empirical (matched to Higgs mass 125 GeV)?
- No candidate matches if we want Brannen scale a → Higgs VEV.

**Q1-3** (Scale bridge): Brannen scale `a ≈ 17.7 √MeV ≈ 314 MeV²` versus
SM Higgs VEV `v ≈ 246 GeV` — **6 orders of magnitude apart**.  What
structural factor bridges them?
- `(246 GeV)² / (314 MeV²) ≈ 1.92 × 10⁸` — no obvious structural origin.
- v59 has dim Spin(7) = 21, but 21^? doesn't give 10⁸ cleanly.
- Maybe the Higgs is a COMPOSITE of v59-Brannen fields?

**Q1-4** (Gauging): Does gauging the silent SU(2)_L explicitly absorb the
3 Goldstones into W±, Z⁰?
- Need to write `D_μ ξ = ∂_μ ξ + [A_μ, ξ]` with A_μ ∈ su(2).
- Add YM kinetic term `-(1/4 g²) F^a_μν F^{aμν}`.
- Compute gauge boson masses from `⟨ξ⟩`-induced Higgs mechanism.

**Q1-5** (Yukawa): Explicit form of the Yukawa term that gives fermions
mass from `ξ(x)` at vacuum?
```
   L_Yuk  =  Σ_X y_X · ψ̄_X · M_X(ξ) · ψ_X + h.c.
```
- For each sector X, M_X(ξ) is a 3×3 matrix on generations.
- ψ_X is a fermion field with appropriate SM quantum numbers.
- At vacuum, gives masses → m_X^k = y_X · s_X^k (Brannen amplitudes).

**Q1-6** (L⊕F in dynamics): How does the L ⊕ F decomposition of
Cl(7)_even enter the dynamics?
- Pure ξ ∈ ℍ is 4-dim; L ⊕ F is in Cl(7)_even (64-dim).
- For multi-sector dynamics, ξ probably needs to be embedded in
  Cl(7)_even rather than ℍ alone.

### Roadmap for Frontier 1

**Next session (single):**
1. Choose between Option A and Option D (most physical).  Argue from
   field-content economy vs. composite-quark plausibility.
2. Gauge the silent SU(2)_L explicitly.  Compute m_W from `g · ⟨ξ⟩`.
3. Test if v59 scale a × (some specific structural factor) = SM Higgs VEV.

**Multi-session (Frontier 1 completion):**
4. Write the FULL Furey ℂ⊗ℍ⊗𝕆 field Φ(x) Lagrangian.
5. Derive sector-specific equilibria from a single potential.
6. Compute the Higgs mass from structural inputs.
7. Predict W±, Z⁰ masses.

---

## Frontier 2: Mechanism for the selection rule

**Question**: Why does each Furey N-sector couple to its specific
(Bit-L, Bit-F) subset of Cl(7)_even?

### Status

Empirical observation (`cosserat_experiment/16_Z2_decomposition.py`):
```
   N=0 (lepton)  : (Bit-L=1, Bit-F=0)  → L only → 28
   N=1 (d-quark) : (Bit-L=0, Bit-F=1)  → F only → 35
   N=2 (u-quark) : (Bit-L=1, Bit-F=1)  → L⊕F  → 63
   N=3 (neutrino): (?, ?)                 → ? → ?
```
Seven hypotheses tested, none derived.

### Open questions

**Q2-1** (Lagrangian mechanism): Is there a specific v59-natural
Lagrangian whose Yukawa term has sector-dependent Cl-grade projectors?

**Q2-2** (Why d-quark skips L): The most puzzling case — d-quarks
physically transform under gauge fields, so why is Bit-L = 0 for them?

**Q2-3** (Neutrino D_3): What's the predicted D_3 for the neutrino sector?
- Case A: D_3 = 0 (Brannen ansatz inapplicable to neutrinos)
- Case B: D_3 = 35 (same as d-quark)
- Case C: Other v59-structural value

**Q2-4** (Spin(8)/Spin(7)/G_2 ladder): Does the symmetry-breaking pattern
naturally select (Bit-L, Bit-F) for each sector?

### Roadmap for Frontier 2

**Next session:**
1. Try the **Furey Fock-space derivation**: for each fermion state
   |Ω_N⟩, identify which Cl(7)_even operators give non-trivial diagonal
   action (mass term).  Compare to the (Bit-L, Bit-F) pattern.
2. Test the **Spin(8) → Spin(7) → G_2 breaking** pattern.  At each step,
   identify broken generators and see if they map to L vs F.

**Multi-session:**
3. Write the full Brannen Yukawa Lagrangian in ℂ⊗ℍ⊗𝕆.  Check whether
   the L ⊕ F structure emerges automatically from the form.
4. Look at WZW terms on the Spin(8) coset.

---

## Frontier 3: Lagrangian derivation of (5, 21/16) prefactors

**Question**: Why do the v59-tier conjectures have prefactors of 5 (g_W²)
and 21/16 (G_e)?

### Status

- g_W² = 5·√α matches at 0.28 %.
- G_e = (21/16)·α²¹ matches at 0.25 %.
- The 5 has *two* equivalent structural interpretations:
  - Killing-form embedding index of so(3) ⊂ so(7).
  - dim Spin(7) − dim Cl(3,1) = 21 − 16.
- The 21/16 = dim Spin(7) / dim Cl(3,1).
- Both prefactors use ONLY the two integers {21, 16}.

Standard Yang-Mills scenarios (A-E in `10_YM_lagrangian.py`) don't give
either prefactor cleanly.

### Open questions

**Q3-1** (Source of 5): What Lagrangian normalization makes g_W² ∝ √α
rather than g_W² ∝ α or g_W² independent of α?

**Q3-2** (Source of 21/16): What measure/Jacobian factor in the gravity
action produces (21/16) specifically?

**Q3-3** (Unification): Both prefactors involve {21, 16} via different
operations (5 = subtraction, 21/16 = ratio).  Is there a single
mechanism that produces both?

**Q3-4** (Spin(7) Yang-Mills): Does a Spin(7) gauge theory broken to
SU(2)_L give the right embedding ratios with the v59 instanton structure?

### Roadmap for Frontier 3

**Next session:**
1. Explicit YM coupling-flow calculation: Spin(7) gauge with coupling
   `g_Spin(7) = α^(1/4)` → broken to SU(2)_L with `g_W = √5 · g_Spin(7)`.
   Compute the symmetry-breaking pattern that natural produces this.
2. Test the **Pontryagin density** interpretation for G_e: gravity might
   arise from a specific Pontryagin class of a Spin(7) bundle, with the
   (21/16) being a normalization factor.

**Multi-session:**
3. Write the FULL gauge Lagrangian with all sectors (SU(3)_c, SU(2)_L,
   U(1)_em, gravity) and check whether the prefactors are forced.
4. Verify with Lean (formalize the embedding-index calculation in matrix
   Lie algebra terms).

---

## Frontier 4: Scale bridges

**Question**: How do the v59-natural scales relate to the empirical scales?

### Status

We have:
- Brannen scale `a ≈ 17.7 √MeV` (sets lepton masses)
- v59 EM instanton action `S_em = π²/2 ≈ 4.92` (sets α)
- SM Higgs VEV `v ≈ 246 GeV`
- Planck mass `M_P ≈ 1.22 × 10¹⁹ GeV`
- Electroweak/QCD scales

The Brannen scale `a` is set EMPIRICALLY by `(√m_e + √m_μ + √m_τ)/3`.  No
structural derivation in v59.

### Open questions

**Q4-1** (Brannen → Higgs): What structural factor connects `a` to `v_Higgs`?
- Factor needed: `(v/a)² ≈ (246 GeV / 17.7 √MeV)² ≈ 1.92 × 10⁸`.
- Doesn't match any obvious v59-natural quantity.

**Q4-2** (Higgs → Planck): The hierarchy problem.  Standard SM puzzle.

**Q4-3** (v59 lepton mass scale `a`): Can a/M_P be derived from v59
structural inputs?

**Q4-4** (cosmological constant): Where does Λ live in v59?

### Roadmap for Frontier 4

**Next session (any small step):**
1. Test specific candidate bridges:
   - `v/a = √(dim Furey × something)`?
   - `v/a = α^(-some-power)`?
   - `v/a = M_P × structural`?

**Multi-session:**
2. Connect to the full Furey ℂ⊗ℍ⊗𝕆 → SM via a specific scale-setting
   mechanism (dimensional transmutation, anomalous dimension, etc.).
3. Derive Planck mass / dimensionless gravity from v59 structure.

---

## Minor technical items (Lean)

These would tighten the formalisation but don't open new physics.

**Q-L1** (Killing-form general formula): Prove `B_{so(N)}|_{so(n)} = (N-2)/(n-2) · B_{so(n)}` from the actual Killing-form definition, not just stated as a definition.

**Q-L2** (Brannen matrix eigenvalue formula): Lean-formalize the cyclic-shift matrix's eigenvalues `λ_k = a(1 + ξω^k + ξ̄ω^{-k})` directly from the matrix calculation.

**Q-L3** (Quark sector matrix interpretation): Same for quarks — the t² = 3/5, 7/9 fits should be matched to specific Brannen-matrix eigenvalues for each sector.

---

## Speculative / longer-term issues

These need new theory:

**Q-S1**: Full Standard Model gauge group emergence (especially U(1)_Y and Higgs sector).

**Q-S2**: CKM mixing matrix from v59.

**Q-S3**: PMNS matrix (neutrino mixing).

**Q-S4**: Cosmological constant / dark energy in v59 framework.

**Q-S5**: Origin of three generations (currently from Z₃ triality but mechanism not derived).

**Q-S6**: Origin of CP violation.

---

## Priority for the next session (single, focused)

If you have one session to push v59 forward, here's the recommended
sequence (in decreasing order of expected payoff):

### Tier 1 (highest probability of clean result)

1. **Choose Option A vs Option D** for ξ multi-sector formulation.
   - Pro of D: uses additive identity D_u = D_e + D_d structurally;
     elegant; predicts u-quark as bound state of (lepton, d-quark).
   - Pro of A: simplest; allows direct Yukawa per sector.
   - **Decision needed**: which physical interpretation to commit to.

2. **Gauge the silent SU(2)_L explicitly**.  Write D_μ ξ = ∂_μ ξ + [A_μ, ξ];
   compute m_W from ⟨ξ⟩.  This is the explicit verification that the
   v59 silent direction IS the SM SU(2)_L.

### Tier 2 (more speculative but lower-hanging if it works)

3. **Test scale bridge candidates**: try `v/a = (some structural)`.
   Look for ratios with `α^k` or `(dim X)^p` matches.

4. **Test Furey Fock-space derivation of selection rule**.  For each
   |Ω_N⟩ fermion state, compute which Cl(7)_even operators give
   non-trivial mass terms.  Compare to (Bit-L, Bit-F) pattern.

### Tier 3 (deeper, longer)

5. **Write the FULL Furey ℂ⊗ℍ⊗𝕆 field Lagrangian** with all sectors,
   coupling constants, and Yukawa terms.

---

## Down to project: explicit action items for next session

If you want me to start immediately on a specific item next time, here are
the concrete entry points in priority order:

### A. Gauging the silent SU(2)_L explicitly

Files to create:
- `v59/synthesis/gauged_xi_lagrangian.md` — gauged Lagrangian writeup
- `v59/synthesis/gauged_xi.mac` — Maxima for gauge boson mass calculation
- Possibly `lean/GaugedSilentSU2.lean` — Lean-encode the gauge structure

Specific calculation:
- Set up `D_μ ξ = ∂_μ ξ + (g/2) [A^a τ^a, ξ]` where τ^a are Pauli matrices.
- Compute |D_μ ξ|² at vacuum ⟨ξ⟩ = (v_ξ, 0, 0, 0).
- Read off m_A^a from the quadratic gauge-boson terms.
- Expected: m_W ∝ g · v_ξ; compare to v59-predicted g_W = √5 · α^(1/4).

Expected output: confirmation that the silent SU(2)/U(1) is literally
SU(2)_L; gauge-boson mass formula in v59 structural terms.

### B. Scale bridge investigation

Files to create:
- `v59/synthesis/scale_bridge.py` — numerical scan of candidate ratios
- `v59/synthesis/FINDINGS_scale.md` — findings doc

Specific scan:
- Empirical `v_Higgs/Brannen-a` ratio = 246 GeV / (17.7 √MeV) ≈ √(1.92×10⁸)
- Try: `α^p · (dim G_2)^q · (dim Spin(7))^r · ...` for various integer p, q, r.
- Look for ratios within 1 % of the empirical.

### C. Option A vs Option D analysis

Files to create:
- `v59/synthesis/OPTIONS_A_D.md` — pro/con analysis with concrete predictions
- `v59/synthesis/composite_quark.py` — if Option D, test bound-state
  interpretation

Specific test for D:
- u-quark ambient = L ⊕ F.  If u-quark is a bound state of (lepton, d-quark),
  then u-quark mass should come from `m_L + m_F` or `m_L * m_F` or similar.
- Check empirically: m_u / (m_lepton × m_d_quark) — does it give a clean
  structural number?

---

## Cross-references

- **Master session record**: [`SESSION_2026-05-22.md`](SESSION_2026-05-22.md)
- **Prediction summary**: [`SUMMARY.md`](SUMMARY.md)
- **Lean development**: [`furey_construction/lean/README.md`](furey_construction/lean/README.md)
- **Furey program**: [`furey_construction/ALL_FINDINGS.md`](furey_construction/ALL_FINDINGS.md)
- **Cosserat experiments**: [`cosserat_experiment/FINDINGS.md`](cosserat_experiment/FINDINGS.md)
- **Synthesis (v58⊕v59)**: [`synthesis/SYNTHESIS.md`](synthesis/SYNTHESIS.md),
  [`synthesis/DYNAMIC_XI.md`](synthesis/DYNAMIC_XI.md)

## What this roadmap is NOT

- It's NOT a guarantee that each next-session attempt will succeed.
- It's NOT a multi-month research plan — each item could spawn its own
  multi-session investigation.
- The "Tier 1" priority isn't sacred; if you have specific insights or
  new connections, redirect.

## What it IS

- A consolidated list of every open question identified by 2026-05-22.
- Concrete next-session entry points with specific files and
  calculations.
- Honest about uncertainty: many of these may not have clean answers.
