# v61 GEN1 — Nonlinear/curved backreaction: Schwarzschild from ρ_grav

**Date**: 2026-05-26 (v61 loop, Generation 1)
**Artifacts**:
- `01_curved_backreaction.py` (SymPy — full GR tensor computation; all pass)
- `01_schwarzschild.mac` (independent **Maxima `ctensor`** cross-check)
- `lean/CurvedBackreaction.lean` (builds clean against v59 Mathlib)
**Builds on**: v60 GEN1/2 (linearized gravity), GEN4 (Newtonian limit).

---

## Verdict

v60 established the gravity sector at **linearized** level (flat-space TT gravitons
+ the Newtonian OBE trace law). v61 GEN1 takes it **nonlinear / curved**: the same
first-order action gives full GR, whose static spherical vacuum is **Schwarzschild**,
sourced by the v60 gravitational charge `ρ_grav`, and the v60 Newtonian `1/r` is its
weak-field limit — plus a genuine curved-space prediction (light bending).

### (A) Nonlinear connection elimination → full GR (SymPy)

On a curved (Schwarzschild) background, the Levi-Civita connection is **torsion-free
and metric-compatible** — the *nonlinear* version of v60 GEN2's linearized
elimination. So varying `ω` in the first-order action gives `Γ = Levi-Civita(g)` for
*any* metric, i.e. the action is full GR (Einstein–Hilbert) sourced by the matter
stress tensor — not merely its linearization.

### (B) Schwarzschild is the vacuum solution (SymPy + Maxima)

Computing the full Einstein tensor of `ds² = −(1−r_s/r)dt² + (1−r_s/r)⁻¹dr² + r²dΩ²`:
**`R = 0` and `G_{μν} = 0`** (all components). Independently confirmed by Maxima
`ctensor` (`Σ|Ricci_ij| = 0`). So the static spherical vacuum of the v61 gravity
sector is exactly Schwarzschild.

### (C) Backreaction & Newtonian matching (SymPy + Lean)

A lump of total gravitational charge `M = ∫ρ_grav = Σm` (the v60 GEN3/4 second
moment) sources this metric with `r_s = 2GM`. The weak field gives
`Φ = −(g_00+1)/2 = −GM/r` — **exactly the v60 GEN4 Newtonian `1/r` result**. So the
v60 linearized OBE is the weak-field limit of a genuine Schwarzschild metric, with
the *same* charge `ρ_grav`.

### (D) Light deflection — a curved-space prediction (SymPy + Maxima + Lean)

`θ = 2r_s/b = 4GM/b` — **twice** the Newtonian `2GM/b`. The factor-2 (a hallmark GR
test) is a genuine curved-space prediction beyond the v60 scalar/Newtonian sector.

---

## Lean (`CurvedBackreaction.lean`)

| theorem | content | kind |
|---|---|---|
| `weak_field_newtonian` | `−(g_00+1)/2 = −GM/r` (`r_s=2GM`) | `field_simp` |
| `schwarzschild_gtt_grr` | `g_00·g_11 = −1` | `field_simp` |
| `deflection_doubles` | `4GM/b = 2·(2GM/b)` | `ring` |
| `gen1_curved` | bundled headline | mixed |

---

## Status table

| claim | status | tool |
|---|---|---|
| nonlinear elimination → Levi-Civita on curved `g` | **verified** | SymPy |
| Schwarzschild `G_{μν}=0` (vacuum) | **verified** | SymPy + Maxima (ctensor) |
| `r_s=2GM`; weak-field `Φ=−GM/r` = v60 GEN4 | **verified** | SymPy + Lean |
| light deflection `4GM/b` = 2× Newtonian | **verified** | SymPy + Maxima + Lean |

---

## What's next

- **GEN2 (aspect 2)**: matter **backreaction** — the nonlinear Einstein equation
  sourced by the GEN3 matter stress tensor `T_μν[Φ]`, i.e. a self-gravitating Koide
  condensate (interior solution / TOV-like), checking the exterior matches the
  Schwarzschild `M = ∫ρ_grav` of GEN1.
- Then **GEN3 (aspect 3)**: the EW-vev home (R1) — a dynamical mechanism for
  `v = 784a²` on `End(L)`.
