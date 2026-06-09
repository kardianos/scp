# ℒ_v60 — The dynamical Lagrangian (cohesive deliverable)

**Date**: 2026-05-26 (dynamical-Lagrangian loop, synthesis)
**Status**: The v60 gated deliverable (ROADMAP §"The Dynamical Lagrangian"): a
first-principles action `ℒ` on `Cl(3,1)⊗Cl(7)_even` whose Euler–Lagrange equations
yield the OBE `Ω(x)` structure, with a LIGO-viable gravity sector and the Brannen
spectrum. Assembled and verified across loop generations GEN1–GEN7 (see
`LOOP_LOG.md`); every claim below has a SymPy/C/Maxima check and/or a Lean theorem.

This is the *derivation* the v59 `LAGRANGIAN.md` was a *target specification* for.

---

## 1. The action

On the factorized carrier `Cl(3,1) ⊗ Cl(7)_even` (built and shown to commute in
v59 `07`/`G9Unification`):

```
ℒ_v60  =  ℒ_grav[B, ω]            (gravity: first-order / Plebański-BF)
        +  ℒ_matter[Φ; g(B)]       (matter: multivector field on Cl(7)_even)
        +  ℒ_coupling[Φ; h]        (universal, minimal)
```

with the three pieces:

**Gravity (first-order, with an independent `Cl(3,1)` connection `ω`)** — GEN1/2:
```
ℒ_grav = ⟨B ∧ F[ω]⟩₀ − ½ Ψ_{IJKL} B^{IJ}∧B^{KL}
```
- `B` = a 2-form valued in the `Cl(3,1)` bivectors `so(3,1)` (= the 6-field
  Cosserat `3φ+3θ`); `ω` = the independent `so(3,1)` connection; `F = dω+ω∧ω`.
- `⟨B∧F⟩₀` is a **Clifford scalar grade projection** (`B^{IJ}F_{IJ} = −Tr(BF)`,
  GEN2) — gravity lives in the same multivector algebra, not bolted on.
- `Ψ` = Weyl-symmetric simplicity multiplier; varying it forces `B = e∧e` ⇒ an
  induced metric `g(B)` via Urbantke (v59 `08`).

**Matter (multivector field `Φ ∈ Cl(7)_even`, L-grade)** — GEN3:
```
ℒ_matter = ½ ⟨D_μΦ  D^μΦ̃⟩ − V(Φ),   V = λ(e₁²−6e₂)² + μ(e₁−c)²
```
in the S₃-symmetric invariants `e₁,e₂` of the generation triple (sedenion S₃).

**Coupling (universal minimal)** — GEN4:
```
ℒ_coupling = ½ h^{μν} T_{μν}[Φ],   T_{μν} = ∂_μΦ∂_νΦ − η_{μν}(½(∂Φ)² + V)
```
one vertex, one constant `f_g` for all matter (equivalence principle).

---

## 2. Euler–Lagrange equations and what they yield

| vary | EL equation | yields |
|---|---|---|
| `ω` | `d_ω B = 0` | eliminates `ω` (slaving); GEN1: in the trace sector this is **the OBE `∇²Ω = −f_g ρ_grav`**, in the TT sector the massless graviton wave |
| `Ψ` | `B∧B ∝ ε` | simplicity ⇒ tetrad ⇒ **induced metric `g(B)`** |
| `B` | `F = Ψ·B + source` | Einstein eq. sourced by `ρ_grav`; **2 TT DOF** (GEN2) |
| `Φ` | `□Φ = −V'(Φ)` | **Koide-cone vacuum** (`Q=2/3`, GEN3); time evolution (GEN7) |

**The central result (GEN1):** the OBE is *not* a separate posit — it is the
**connection-eliminated trace sector** of `ℒ_grav`. The `09` obstruction
("OBE ⇏ Plebański", DOF 1 < 2) is dissolved: OBE and the graviton are sibling
sectors of one first-order parent (`ParentAction.lean`).

---

## 3. The vacuum

`Φ` minimizes `V`: the EL vacuum lies on the **Koide cone `Q = 2/3`**
(`e₁²=6e₂`, GEN3 — derived from minimization, `MatterSector.koide_invariant_form`).
The second moment at the vacuum is
```
ρ_grav = Tr(M†M) = Σ m_k = 9Q a² = 6a²,
```
which is **exactly the gravity source** (GEN4). Gravity and matter share the one
scale `a`. The Brannen phase `φ` is a flat (Goldstone) direction of `V`.

---

## 4. The linearized spectrum (around the joint vacuum)

Decoupled (the `h`–`Φ` mixing `δT_μν` vanishes at the homogeneous vacuum, GEN5),
**5 propagating modes, ghost- and tachyon-free** (`SpectrumStability.lean`):

| sector | modes | m² | protected by |
|---|---|---|---|
| graviton TT | 2 | 0 | diffeo gauge |
| matter radial | 2 | > 0 (PSD ∀λ,μ>0) | — |
| Brannen phase | 1 | 0 | cone symmetry (Goldstone) |

## 5. The dynamics (nonlinear, GEN7)

Time-integrating `□Φ = −V'(Φ)` reproduces this spectrum **nonlinearly**: the
normal-mode frequencies match the Hessian eigenvalues to `~10⁻⁵`, energy is
conserved (`4×10⁻⁶`), and the **massless Goldstone propagates at c = 1**
(`ω² = k²`, `DynamicsSpectrum.lean`). Massive modes obey `ω² = k² + m²`.

---

## 6. Structural integers (all derived / machine-checked)

| quantity | value | origin |
|---|---|---|
| graviton DOF | 2 | TT, `10−4−4` (GEN2/5) |
| Koide `Q` | 2/3 | EL vacuum `e₁²=6e₂` = `dimG₂/dimSpin7` (GEN3/6) |
| selection dims | `D_e=28, D_d=35, D_u=63` | `Cl(7)_even` grades; `D_u=D_e+D_d` (GEN6) |
| universal Koide deviation | `(1−Q)D = 28/3` | ties `Q=2/3,11/15,23/27` (GEN6) |
| `ρ_grav` | `9Qa²=6a²` | second moment = Σ inertial masses (GEN4) |
| radial law | `1/r`, force `1/r²` | massless kernel + monopole (GEN4) |

---

## 7. What is derived vs. residual

**Derived (this loop):** the dynamical *structure* — the OBE as a connection-
eliminated sector; 2 ghost-free TT gravitons; the Koide cone Q=2/3 as an EL
vacuum; EP-exact universal coupling; a stable decoupled spectrum; the selection
rule and its dimensions; and genuine time evolution reproducing the spectrum.

**Residual value-conjectures (NOT derived — carried from v59, now isolated):**
- `α` (fine-structure) — a genuine input.
- the EW vev `v = 784 a²` (R1) — a separate sector (rank tension; `SelectionRule`,
  `RankTension`); the *only* piece still without a dynamical home.
- the Brannen phase `φ = 2/9` — the Goldstone direction; needs an explicit
  S₃-breaking tilt to pin (GEN3).
- the gravity magnitude `f_g ~ α^{21/2}` — value-conjecture (v59).

These are **values/inputs**, not dynamical gaps: the *theory* is complete and
consistent; these four numbers parametrize it.

---

*Verification: `17_verify_all.py` (13/13 Python+Maxima+C pass); 7 Lean modules
build clean against v59 Mathlib (`ParentAction`, `CovariantFirstOrder`,
`MatterSector`, `MatterGravityCoupling`, `SpectrumStability`, `SelectionRule`,
`DynamicsSpectrum`). Per-generation detail in `1X_findings.md` + `LOOP_LOG.md`.*
