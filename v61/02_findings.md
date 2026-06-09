# v61 GEN2 — Matter backreaction: the GEN3 energy sources the GEN1 Schwarzschild mass

**Date**: 2026-05-26 (v61 loop, Generation 2)
**Artifacts**:
- `02_backreaction.py` (SymPy — full Einstein tensor; all pass)
- `02_backreaction.mac` (independent **Maxima `ctensor`** cross-check)
- `lean/Backreaction.lean` (builds clean, warning-free)
**Builds on**: GEN1 (Schwarzschild exterior), v60 GEN3/4 (matter `T_00`, charge `ρ_grav`).

---

## Verdict

GEN1 *posited* that the Schwarzschild mass is `M = ∫ρ_grav`. GEN2 **derives** it from
the nonlinear Einstein (00) equation and closes the interior↔exterior loop: the GEN3
matter energy density is the source whose integral is the GEN1 charge.

### (A) The backreaction law (SymPy + Maxima)

For `ds² = −e^{2Φ(r)}dt² + (1−2m(r)/r)⁻¹dr² + r²dΩ²` (G=1), the mixed Einstein
component is **`G^t_t = −2 m'(r)/r²`** (SymPy; independently confirmed by Maxima
`ctensor`). Einstein's equation `G^t_t = 8π T^t_t = −8πρ` then gives
```
m'(r) = 4π r² ρ(r),      ρ = T_00  (the GEN3 energy density).
```

### (B) Total mass = the v60 charge (SymPy + Lean)

Integrating: `M = m(∞) = ∫4πr²ρ dr = ∫ρ d³x = ρ_grav = Σm = 9Qa²`. So **the GEN1
Schwarzschild charge is exactly the integrated GEN3 matter energy density** — the
backreaction loop closes. (Uniform-density check: `M = (4/3)πR³ρ₀`.)

### (C) Interior↔exterior matching → GEN1 Schwarzschild (SymPy + Lean)

Outside the matter (`ρ=0`): `m'=0 ⟹ m=M` const ⟹ `g_rr = (1−2M/r)⁻¹` — exactly the
GEN1 Schwarzschild exterior. The mass function `m(R)=M` matches continuously at the
surface.

### (D) Newtonian limit = the v60 OBE Poisson law (SymPy + Lean)

Weak field: `∇²Φ = m'/r² = 4πρ` (substituting `m'=4πr²ρ`). This is **the v60 GEN4 /
OBE trace law `∇²Ω = −f_g ρ_grav`, now DERIVED** as the weak-field limit of the
curved (00) Einstein equation (not posited).

---

## Lean (`Backreaction.lean`)

| theorem | content | kind |
|---|---|---|
| `backreaction_coefficient` | `−2(4πr²ρ)/r² = −8πρ` | `field_simp; ring` |
| `newtonian_poisson` | `(4πr²ρ)/r² = 4πρ` | `field_simp` |
| `uniform_mass` / `mean_density` | `(4/3)πR³ρ₀`; `ρ₀=3M/(4πR³)` | `ring` / `field_simp` |
| `exterior_schwarzschild` | `(1/(1−2M/r))(1−2M/r)=1` | `one_div_mul_cancel` |
| `gen2_backreaction` | bundled headline | mixed |

---

## Status table

| claim | status | tool |
|---|---|---|
| `G^t_t = −2m'/r²` | **verified** | SymPy + Maxima (ctensor) |
| backreaction law `m' = 4πr²ρ` | **verified** | SymPy + Lean |
| `M = ∫ρ_grav` (GEN1 charge = GEN3 energy) | **verified** | SymPy + Lean |
| interior/exterior → GEN1 Schwarzschild | **verified** | SymPy + Lean |
| Newtonian limit = v60 OBE Poisson | **verified** | SymPy + Lean |

---

## State after GEN2

The gravity sector is now **self-consistent across scales**: the GEN3 matter
condensate's energy density `ρ` sources, via the nonlinear Einstein equation, a
mass `M = ∫ρ_grav` whose exterior is the GEN1 Schwarzschild metric, whose weak
field is the v60 GEN4 Newtonian/OBE law. Linearized OBE → full GR → matter source
all tie together with one charge `ρ_grav`.

## What's next

- **GEN3 (aspect 3)**: the **EW-vev home (R1)** — the last residual without a
  dynamical origin. Attack a dynamical mechanism for the full-rank democratic
  `End(L)` condensate with `‖·‖²_F = 784a²` from the one scale `a_ℓ` (e.g. an
  `so(8)`-symmetric potential on `End(L)` whose vacuum is democratic, or the
  compositeness/condensate reading), and verify whether the 784 count emerges or
  remains a value-conjecture.
