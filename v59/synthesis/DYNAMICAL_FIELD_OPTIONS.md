# Dynamical Field Options — Analysis and Validation Plan

**Date**: 2026-05-22
**Parent**: [`MOTIVATIONS_AND_CONSEQUENCES.md`](MOTIVATIONS_AND_CONSEQUENCES.md),
[`../INTEGRATION.md`](../INTEGRATION.md)
**Related**: `DYNAMIC_XI.md` (initial setup), all Step 1-5 work

This is the analysis document for "what is the v59 dynamical field?" — i.e.,
what is the bosonic structure that the Lagrangian acts on. The answer to
this question gates everything downstream (Lagrangian, simulation, predictions).

---

## The question

In v59, ξ is treated as the Brannen quaternion that gives mass via
$M = a(I + \xi S + \bar\xi S^2)$. But ξ has been a *static parameter* until
now. To write the Lagrangian, we must commit to:

- What's the spacetime-dependent dynamical field?
- What's its algebraic structure?
- What's the relation between sectors (lepton, d-quark, u-quark)?

---

## Five options (originally A-D, now extended to E)

### Option A — Sector-specific quaternion fields

$$
\xi_e(x), \xi_d(x), \xi_u(x), \xi_\nu(x) \in \mathbb{H}
$$

Each fermion sector has its own ξ field, with sector-specific constraint
$|\xi_X|^2 = 1 - 14/D_X$.

**Pros**:
- Simplest field content per sector
- Each ξ_X has its own constraint surface $S^3$ ⊂ ℍ
- Yukawa is straightforward: $y_X \bar\psi_X M_X(\xi_X) \psi_X$
- Each sector's Higgs sector has the right 3 Goldstones + 1 radial spectrum

**Cons**:
- 3 (or 4) separate ℍ-valued fields multiplies field content
- Doesn't unify sectors at the dynamical level
- The Pati-Salam Spin(7) unification (Step 5) doesn't act naturally on independent ξ_X
- The scale bridge $v_H = 28^2 a_\ell^2$ would relate only the lepton ξ

**Validation status**: PARTIAL — works per sector but doesn't unify.

---

### Option B — Single ξ ∈ ℍ with multi-minimum potential

$$
V(|\xi|^2) = \lambda \cdot \prod_X (|\xi|^2 - r_X^2)^2
$$

One field, multiple vacuum branches.

**Pros**:
- Single field — most economical
- All three constraint surfaces are critical loci of one V

**Cons**:
- The "vacuum branch each fermion sees" must be sector-specific — ad hoc
- Goldstone count per branch is the same (3+1 from O(4)/O(3))
- ξ ∈ ℍ is only 4-dim — too small for the 64-dim Furey content
- After Step 1 (16 fermion states per gen in ℂ⊗ℍ⊗𝕆), the Higgs needs more structure

**Validation status**: INSUFFICIENT — too small to carry all SM info.

---

### Option C — Single ξ ∈ ℍ with sector-dependent effective potential

Phenomenological: each sector's local matter content modifies the effective V.

**Validation status**: NOT FUNDAMENTAL — placeholder for refined options.

---

### Option D — Composite quarks from leptons (additive D_u = D_l + D_d)

ξ_lepton is fundamental; quarks emerge as bound states using D_u = D_l + D_d
additive identity.

**Pros**:
- Reduces field content
- Uses an established v59-structural identity

**Cons**:
- Direct test of multiplicative composition: $m_u \approx m_e \cdot m_d$?
  $0.511 \times 4.67 \approx 2.4$ MeV — close to $m_u \approx 2.16$ MeV by coincidence,
  but $m_t \neq m_\tau \cdot m_b$ (1.78 × 4.18 = 7.4 GeV ≠ 173 GeV).
- Additive: $m_t \neq m_\tau + m_b$ (1.78 + 4.18 = 6 GeV ≠ 173 GeV).
- Composite doesn't reproduce quark mass spectrum

**Validation status**: FALSIFIED at the mass-spectrum level.

---

### Option E — Φ ∈ Cl(7)_even with sector-specific Cl-grade projections

The dynamical field $\Phi(x)$ lives in $\mathrm{Cl}(7)_\text{even} \cong \mathbb{C}\otimes\mathbb{O}$ (64 cplx dim).
Each fermion sector couples to a different Cl-grade sub-ambient:
- Lepton: $\Phi_L = P_L \Phi \in L = \Lambda^2 \oplus \Lambda^6$ (28-dim)
- d-quark: $\Phi_F = P_F \Phi \in F = \Lambda^4$ (35-dim)
- u-quark: $\Phi_{L\oplus F} = P_{L\oplus F} \Phi \in L \oplus F$ (63-dim)

**Pros**:
- Unifies all sectors as projections of one Φ
- Pati-Salam Spin(7) acts naturally
- Selection rule (Step 4) is automatic via μ-eigenspace projection
- 64-dim ambient is the natural Furey color algebra
- Compatible with anomaly cancellation (Step 1)

**Cons**:
- 64-dim is HEAVY — many components are "auxiliary"
- Constraint surfaces |ξ_X|² = r_X² must be self-consistently embedded in Φ
- The 4-dim Higgs doublet structure must emerge from a 4-dim slice
- If sectors share an orthogonal-projection decomposition,
  $r_L^2 + r_F^2 = 1/2 + 3/5 = 11/10 \neq 7/9 = r_{L\oplus F}^2$ — **a potential inconsistency**

**Validation status**: REQUIRES TEST — see below.

---

### Option E* — Refined: sector-specific ℍ slices embedded in Cl(7)_even

Each sector has its own ℍ_X ⊂ Cl(7)_even (a 4-dim slice). The slices are
NOT orthogonal projections of a single Φ; instead, they're sector-specific
embeddings related by Pati-Salam Spin(7) rotations.

This is a hybrid of A and E:
- Sector-specific ξ_X ∈ ℍ_X (like Option A)
- All ℍ_X live in Cl(7)_even (Option E unification)
- The Spin(7) acts on the unified Cl(7)_even

**Pros**:
- Each sector's Higgs has correct 3+1 Goldstone+radial spectrum
- Constraint surfaces are sector-local (no inconsistency)
- Spin(7) unifies the structure

**Cons**:
- More fields than Option E (4 ℍ_X fields)
- Still need to specify how ℍ_X are related by Spin(7)

**Validation status**: LIKELY winner after the test below.

---

## Constraints any option must satisfy (from Steps 1-5)

| # | Constraint | Source | Falsifies |
|---|---|---|---|
| 1 | v_H = 28²·a_l² (0.07 %) | Scale bridge | Any option not giving this scale |
| 2 | 3 Goldstones + 1 radial per sector vacuum | Higgs phenomenology | Wrong-dim Higgs sector |
| 3 | sin²θ_W = 2/9 from gauging | Step 5 | Wrong gauge structure |
| 4 | Brannen kernel M = a(I+ξS+ξ̄S²) per sector | Step 1 | Wrong flavor structure |
| 5 | (1-Q_N)·D_N = 28/3 universal | Step 4 | Wrong cross-sector pattern |
| 6 | μ-bisection L⊕F | Step 4 | Wrong grading |
| 7 | Pati-Salam Spin(7) | Step 5 | Wrong gauge unification |
| 8 | Anomaly cancellation | Step 1 | Wrong fermion content |
| 9 | sin²θ_C = 7·α | Step 3 | Wrong inter-sector mixing |
| 10 | y_top ≈ 1 (a_u² = 72·a_l²) | Scale bridge | Wrong sector-scale ratio |
| 11 | m_H²/v² = 7/27 | Higgs mass | Wrong Higgs sector |

A passing option must satisfy ALL 11 constraints.

---

## Validation plan

### Phase 1 (this session): Quick falsification test of Option E

Test the strongest version of Option E: single Φ ∈ Cl(7)_even with the
simplest sector-projection potential. If it fails, Option E* is the
fallback. See `validate_option_E.py` for the calculation.

### Phase 2 (next session): Spectrum verification

Once an option survives Phase 1, compute Hessian at the chosen vacuum;
verify 3 Goldstones + 1 radial; check m_H = √(7/27)·v.

### Phase 3 (later): Lattice simulation

Implement the surviving option in `sfa/sim/scp_sim.cu` and verify field-
theoretic consistency at strong coupling.

---

## Current recommendation (going into Phase 1)

**Test Option E first** (single Φ with sector projections). If it fails the
constraint-simultaneity test, **fall back to Option E***  (sector-specific
ℍ fields embedded in Cl(7)_even).

If Option E passes Phase 1, we have a single unified Higgs field with
clean Cl(7)_even structure. If it fails (likely), Option E* is the
working hypothesis and we proceed with sector-specific fields.

---

## Phase 1 result (DONE — see `validate_option_E.py`)

### Falsification

Simple Option E is **FALSIFIED**:
$$
r_L^2 + r_F^2 = \frac{1}{2} + \frac{3}{5} = \frac{11}{10} \;\neq\; \frac{7}{9} = r_{L\oplus F}^2
$$

For orthogonal projectors $P_L \perp P_F$ in $\mathrm{Cl}(7)_\text{even}$,
the additivity $|P_{L\oplus F}\Phi|^2 = |P_L\Phi|^2 + |P_F\Phi|^2$ forces
the sum of sector constraints — but the v59 cross-sector pattern
$r_X^2 = 1 - 14/D_X$ is NOT additive. The three sector vacua cannot be
simultaneously realized by a single Φ with orthogonal sector projections.

This is a STRUCTURAL fact, not a numerical mismatch — the failure ratio
1.1 vs 0.778 is far from any rounding or normalization adjustment.

### Construction: Option E* is the working hypothesis

Option E* (sector-specific $\xi_X \in \mathbb{H}_X \hookrightarrow \mathrm{Cl}(7)_\text{even}$)
passes 9 of 11 constraints directly:

| ✓/? | Constraint |
|---|---|
| ✓ | v_H = 28²·a_l² (lepton scale bridge) |
| ✓ | 3 Goldstones + 1 radial per sector vacuum (O(4)→O(3) on ℍ_X) |
| ✓ | sin²θ_W = 2/9 from gauging silent SU(2)_L |
| ✓ | cos²θ_W = 7/9 = t²_u (cross-sector) |
| ✓ | Brannen kernel per sector |
| ✓ | Sector constraints \|ξ_X\|² = 1 − 14/D_X |
| ✓ | (1−Q_N)·D_N = 28/3 universal |
| ✓ | L⊕F = μ-bisection embedding |
| ✓ | Pati-Salam G_2 × SU(2)_L × SU(2)_R × U(1)_{B-L} |
| ? | sin²θ_C = 7·α — cross-sector mixing needs explicit derivation |
| ? | m_H²/v² = 7/27 — Hessian computation needs explicit potential |

### Sub-options of Option E*

Two sub-options remain open:

**Option E* (multi-Higgs)**: each fermion sector has its own physical
Higgs-like scalar with VEV $v_X = D_X^2 \cdot a_X^2$. Predictions:
- Lepton Higgs ≈ 246 GeV (= SM Higgs ✓)
- d-quark "Higgs analog" ≈ 800 GeV (HL-LHC reach — likely ruled out if real)
- u-quark "Higgs analog" ≈ 90 TeV (FCC reach)

**Option E** (single-Higgs, structured Yukawas)**: there's only ONE physical
Higgs field (the lepton one at 246 GeV). The "sector ξ_X" are
*Yukawa-matrix-structure parameters*, NOT separate fields. Fermion masses:
$m_f = y_f \cdot v_\text{SM}/\sqrt{2}$ universally; the Brannen-Z₃
structure parametrises the Yukawa matrices $y_X^{ij}$.

Pre-experimental constraints favor Option E**:
- SM Higgs couplings to b-quark/c-quark/τ scale linearly with $m_f$, consistent
  with one universal $v_\text{SM} = 246$ GeV (HL-LHC will tighten this)
- No 800 GeV scalar excess in current LHC searches (constrains Option E*)

### Verdict

**Option E** is the working hypothesis going into the Lagrangian:
- Single physical Higgs $\Phi$ in lepton-sector ambient $L = \Lambda^2 \oplus \Lambda^6$
- $v_\text{SM} = 28^2 \cdot a_\ell^2 = 246$ GeV
- Yukawa matrices $y_X^{ij}$ have Brannen-Z₃ structure with sector-specific
  parameters $(a_X, t_X, \varphi_X)$ tied to $|ξ_X|^2 = 1 - 14/D_X$
- The "extra" structure (Cl(7)_even, Pati-Salam Spin(7), L⊕F μ-bisection)
  provides the algebraic SCAFFOLDING for the gauge group and selection
  rules — but the dynamical fields are LEAN: Φ (Higgs), ψ_X (fermions),
  A^a_μ, B_μ, G^a_μ (gauge), g_μν (gravity).

This is essentially the SM with a STRUCTURED Yukawa sector and an
algebraic origin for the gauge group and Higgs VEV.

---

## Cross-references

- `DYNAMIC_XI.md` — initial framing of the question (pre-Step 5)
- `FINDINGS_scale_bridge.md` — scale bridge constraint
- `MOTIVATIONS_AND_CONSEQUENCES.md` — full conjecture list
- `../INTEGRATION.md` — overall integration picture
- `../furey_construction/FINDINGS_ckm_and_selection.md` — Steps 3+4
- `../furey_construction/FINDINGS_u1y_unification.md` — Step 5
- `validate_option_E.py` (to be created) — Phase 1 test
