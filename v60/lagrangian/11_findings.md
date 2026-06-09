# GEN2 — Full covariant first-order action: `B∧F` grade projection, connection elimination, DOF + ghost audit

**Date**: 2026-05-26 (dynamical-Lagrangian loop, Generation 2)
**Artifacts**:
- `11_covariant_firstorder.py` (runs clean; every assertion passes)
- `11_grade_projection.mac` (independent **Maxima** CAS cross-check of Part 1)
- `../lean/CovariantFirstOrder.lean` (compiles clean against v59 Mathlib;
  headline `gen2_covariant_firstorder` **axiom-free**)
**Builds on**: GEN1 (`10_*`, the sector-by-sector connection elimination).

---

## Verdict

GEN1 eliminated the connection in reduced (helicity-0 and TT) channels. **GEN2
does the full covariant 4D theory at once** and confirms it is linearized Einstein
gravity with exactly **2 ghost-free TT degrees of freedom**, with the gravitational
`B∧F` term identified as a **Clifford scalar-grade projection** — the precise bridge
to the `Cl(3,1)⊗Cl(7)_even` multivector framework of `06`/`07`.

Three new results, each independently verified:

### Part 1 — `B∧F` is a scalar grade projection `⟨B F⟩₀` (SymPy + Maxima)

With the `so(3,1)` bivector generators `(M_{IJ})^a_b = δ^a_I η_{Jb} − δ^a_J η_{Ib}`:

```
Tr(M_{IJ} M_{KL}) = 2(η_{IL}η_{JK} − η_{IK}η_{JL})        [residual 0, both CAS]
```

so for `B = ½B^{IJ}M_{IJ}`, `F = ½F^{KL}M_{KL}`:

```
B^{IJ} F_{IJ} = − Tr(B F) = ⟨B F⟩₀   (scalar grade, up to convention sign)
```

Maxima independently reproduced both (`B^{IJ}F_{IJ} = −Tr(BF) = −8` on a test
pair). **The gravity BF Lagrangian's internal index contraction *is* the grade-0
part of the Clifford geometric product** — the gravity sector lives in the same
multivector algebra as the rest of the theory, not as a bolted-on tensor.

### Part 2 — Full connection elimination (SymPy, 4D, all components)

Linearized Palatini with **independent** connection `Γ^a_{bc}` (40 components,
symmetric lower pair) and `h_{μν}` (10):
```
S₂ = η^{μν}(Γ^λ_{λρ}Γ^ρ_{μν} − Γ^λ_{νρ}Γ^ρ_{μλ})
     − h̄^{μν}(∂_λ Γ^λ_{μν} − ∂_ν Γ^λ_{λμ})
```
The connection EOM `∂S₂/∂Γ = 0` is **linear** in `Γ`: `A·Γ = −b` with `A` the
constant `40×40` Hessian. Result:
- **`rank A = 40/40`, nullity 0** ⟹ the connection is eliminated **uniquely**
  (no projective flat direction survives in *this* action's stationary point).
- Substituting back: `S_eff[h]` is **real** and **second-order** — the
  Fierz-Pauli / linearized-Einstein action. This is the full-tensor version of
  GEN1's `g_i = ∂_iΩ`.

### Part 3 — DOF + ghost audit (SymPy)

On the connection-eliminated `S_eff`, at a **null** momentum `k=(1,0,0,1)`:
- `dim ker O(k) = 6`; all **4** linearized-diffeo gauge directions annihilate `O`
  (verified: `O·ξ = 0`). **Physical DOF `= 6 − 4 = 2`** (the TT graviton).
- **Ghost-free**: the TT polarization (`h_+`) kinetic coefficient is
  `S_eff[TT]/A² = ½(k_z² − w²) = −½k²`, whose ratio to a reference healthy
  massless scalar (computed by the *identical* substitution) is **`+1`** — same
  sign ⟹ not a ghost.
- The conformal/trace mode has the opposite (`−3/2`) sign — the textbook
  "wrong-sign" conformal mode of EH — but it is **non-propagating** (not among
  the 2 DOF), so there is **no physical scalar ghost**. Ghost-freedom is secured
  by the DOF count itself.

---

## Lean (`CovariantFirstOrder.lean`, axiom-free)

| theorem | content |
|---|---|
| `connComps_val` | `40` independent connection components (the eliminated system) |
| `connection_nullity_zero` | `40 − 40 = 0` (unique elimination) |
| `spin_decomposition` | Barnes-Rivers `5+3+1+1 = 10` |
| `physical_dof_fierz_pauli` | `10 − 4 − 4 = 2` |
| `physical_dof_kernel` | `6 − 4 = 2` (the SymPy kernel-minus-gauge count) |
| `no_ghost` | `0 < ttToScalarRatio` (TT kinetic sign healthy) |
| `bf_from_trace` | `bfContractSign·2 = −soTraceCoeff` (grade-projection constant) |
| `gen2_covariant_firstorder` | the bundled headline |

---

## Status table

| claim | status | tool |
|---|---|---|
| `Tr(M_IJ M_KL) = 2(η_IL η_JK − η_IK η_JL)` | **verified** | SymPy + Maxima |
| `B^{IJ}F_{IJ} = −Tr(BF) = ⟨BF⟩₀` | **verified** | SymPy + Maxima |
| connection EOM linear, Hessian rank 40, unique elimination | **verified** | SymPy |
| `S_eff` real, 2nd-order (Fierz-Pauli) | **verified** | SymPy |
| exactly 2 physical TT DOF, diffeo-gauge-invariant | **verified** | SymPy |
| TT mode ghost-free (same sign as healthy scalar) | **verified** | SymPy |
| no physical scalar ghost (conformal mode non-propagating) | **verified** | SymPy |
| integer backbone (counts, signs) | **machine-checked, axiom-free** | Lean |

---

## What GEN2 settled, and what's next

- **Settled**: the gravity sector's dynamics are now the **full covariant**
  first-order action, connection-eliminated to ghost-free linearized GR with the
  correct 2 TT DOF, and `B∧F` is a grade projection in `Cl(3,1)`. Combined with
  GEN1 (`PARENT ⟹ OBE` in the trace sector), the gravity half of the dynamical
  Lagrangian is in hand at linearized level.
- **Next (GEN3, aspect 3)**: the **matter / internal sector**. Write a kinetic +
  potential term for the multivector field `Φ ∈ Cl(7)_even` whose EL vacuum
  enforces the **Brannen / Koide constraint surface** (`t² = 1/2`, `φ = 2/9`) and
  whose second moment is exactly the `ρ_grav = Tr(M†M)` source that GEN1/GEN2
  feed into the gravity trace law. That couples the two halves and is the gateway
  to the full linearized spectrum (GEN5).
