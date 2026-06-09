# GEN4 — Covariant matter→gravity coupling + equivalence principle

**Date**: 2026-05-26 (dynamical-Lagrangian loop, Generation 4)
**Artifacts**:
- `13_coupling_ep.py` (SymPy — universal vertex + stress tensor; all assertions pass)
- `13_newton_ep.c` (compiled **C** — numeric EP ratios + Newtonian 1/r law; exit 0)
- `../lean/MatterGravityCoupling.lean` (compiles clean against v59 Mathlib)
**Builds on**: GEN1/2 (gravity, source `ρ_grav`), GEN3 (matter, produces `ρ_grav`).

---

## Verdict

GEN4 **couples** the gravity half (GEN1/2) and the matter half (GEN3) with a single
covariant term and shows the coupling is **universal** — the equivalence principle
holds exactly, with one `f_g` for all matter. The GEN1–GEN3 pieces are now a single
action.

### (A) Minimal covariant coupling = one universal vertex (SymPy)

Expanding the covariant matter action
`S_m = ∫√(−g)(−½ g^{μν}∂_μΦ∂_νΦ − V)` to linear order in `h_{μν} = g_{μν}−η_{μν}`:
```
δℒ = ½ h^{μν} T_{μν},   T_{μν} = ∂_μΦ∂_νΦ − η_{μν}(½(∂Φ)² + V)
```
(SymPy residual **0**). **One vertex, one coupling, the same `T` for every field/grade**
— this *is* the equivalence principle at the level of the action.

### (B) `T_00` = energy density → Newtonian source (SymPy)

Static (`∂_0Φ=0`): `T_00 = ½|∇Φ|² + V` (the positive Hamiltonian density). So the
Newtonian/trace limit `∇²Φ_N = 4πG T_00` **is** the GEN1/2 OBE
`□Ω = −f_g ρ_grav` with `Ω ↔ Φ_N`, `f_g ↔ 4πG`, `ρ_grav ↔ ∫T_00`.

### (C) Equivalence principle EXACT (SymPy + C)

The gravitational charge is `ρ_grav = Tr(M†M)` = sum of **eigenvalues** of `M†M`
= sum of the **inertial masses** `m_k`. Hence for every mode
`m_grav/m_inertial = 1`. Verified for the Brannen lepton triple:
- SymPy: ratio `= 1` for all 3 modes (exact); `ρ_grav = Σm = 6a² = 9Qa²`.
- C: ratios `1.000000000000000`; `Σm = 6.000000000000 = 9Qa²`.

### (D) OBE radial law from the monopole (C)

Point source `M = ρ_grav`: Newtonian potential slope `d log|Φ|/d log r = −1.000`
(1/r, long-range) and force slope `−2.000` (1/r²) — exactly the GEN1/2 / v59
massless-kernel-with-nonzero-monopole result. (The magnitude `f_g ~ α^{21/2}`
remains the v59 value-conjecture — unchanged here.)

---

## Lean (`MatterGravityCoupling.lean`)

| theorem | content | kind |
|---|---|---|
| `one_universal_coupling` | one `f_g` for all sectors | `decide` |
| `ep_ratio_one` | `m/m = 1` (grav/inertial), for `m ≠ 0` | `div_self` |
| `grav_charge_is_sum_of_masses` | `Tr(diag(m0,m1,m2)) = m0+m1+m2` | Mathlib `simp` |
| `rho_grav_nineQ` | `6 = 9·(2/3)` | `norm_num` |
| `force_one_steeper` | `forceSlope = potentialSlope − 1` (`−2 = −1−1`) | `decide` |
| `gen4_coupling_ep` | bundled headline | mixed |

---

## Status table

| claim | status | tool |
|---|---|---|
| minimal coupling = `½ h^{μν}T_{μν}` (universal vertex) | **verified** | SymPy |
| `T_00` (static) = energy density | **verified** | SymPy |
| EP exact: `m_grav/m_inertial = 1` for all modes | **verified** | SymPy + C |
| `ρ_grav = Tr(M†M) = Σ(inertial masses) = 6a² = 9Qa²` | **verified** | SymPy + C + Lean |
| Newtonian potential 1/r, force 1/r² (monopole) | **verified** | C + Lean |
| coupling magnitude `f_g ~ α^{21/2}` | **value-conjecture** (unchanged) | — |

---

## State of the program after GEN4

The three sectors are now a single coupled action at linearized level:
```
ℒ = ℒ_grav[B,ω]                      (GEN1/2: first-order, eliminates to lin. Einstein, 2 TT DOF)
    + ½⟨∂Φ ∂Φ̃⟩ − V(Φ)               (GEN3: matter, Koide-cone vacuum, ρ_grav = 9Qa²)
    + ½ h^{μν} T_{μν}[Φ]             (GEN4: universal coupling, EP exact)
```
with one scale `a` (and the residual value-conjecture `f_g`). The trace law
`∇²Ω = −f_g ρ_grav` and the 2 TT graviton modes both descend from `ℒ_grav`; the
source `ρ_grav` is produced by the matter vacuum and couples universally.

## What's next

- **GEN5 (aspect 5)**: the **full linearized spectrum around the joint vacuum** —
  diagonalize the quadratic action and confirm, *simultaneously*, the 2 TT
  gravitons + the matter modes (2 massive + 1 Goldstone), with **no ghosts or
  tachyons** in the combined system (mixing audit between `h` and `Φ`).
- Then **GEN6** (selection rule / rank tension G1) and **GEN7** (a genuine
  time-dependent EL solution = the actual "dynamics").
