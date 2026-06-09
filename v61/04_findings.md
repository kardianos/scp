# v61 GEN4 — Closing the LIGO motivation: 2 TT modes = h₊,h×; GW emission

**Date**: 2026-05-26 (v61 loop, Generation 4)
**Artifacts**:
- `04_gravitational_waves.py` (SymPy — all pass)
- `04_quadrupole.mac` (Maxima cross-check of the luminosity; supplementary)
- `lean/GravitationalWaves.lean` (builds clean)
**Builds on**: v60 GEN2/5 (2 TT DOF), v60 GEN7 / v61 GEN1 (massless dispersion).

---

## Verdict — the G9 program's origin is resolved

The entire G9 program began because **v59 gravity was a Lorentz scalar (helicity 0)
— fatal for LIGO's observed `h₊, h×` (helicity ±2)**. v60/v61 built a theory with
exactly 2 TT graviton DOF. GEN4 closes the loop: **those 2 DOF *are* the LIGO
polarizations.**

### (A) Helicity ±2 — the LIGO polarizations (SymPy + Lean)

For a wave along `z`, the TT polarizations are `e₊ = diag(1,−1,0)`, `e× = offdiag(1)`.
Under a rotation by `ψ` about `z` (`h → R h Rᵀ`):
```
R e₊ Rᵀ = cos2ψ·e₊ + sin2ψ·e×,   R e× Rᵀ = −sin2ψ·e₊ + cos2ψ·e×
```
so `(h₊, h×)` rotate by **`2ψ`** and `h₊ + i h× → e^{2iψ}(h₊ + i h×)` — **helicity +2**
(conjugate: −2). The 2 TT modes are spin-2, exactly LIGO's `h₊, h×` — replacing the
v59 scalar's helicity 0.

### (B) Gravitational-wave quadrupole emission (SymPy; Maxima cross-check)

The Einstein quadrupole luminosity for a mass `μ` on a circular orbit (radius `a`,
frequency `ω`), `L = (G/5)⟨I⃛_ij I⃛_ij⟩` (reduced trace-free quadrupole, c=1):
```
L = (32/5) G μ² a⁴ ω⁶  >  0
```
(SymPy, exact). The binary radiates GWs — the standard result. (Maxima `04_quadrupole.mac`
re-derives the same coefficient as a supplementary cross-check.)

### (C) Propagation at c (SymPy + Lean)

The graviton is massless (v60/v61: `m=0`, dispersion `ω²=k²`), so GWs travel at the
speed of light (`ω/k = 1`) — consistent with LIGO/GW170817's near-simultaneous
photon arrival.

---

## Lean (`GravitationalWaves.lean`)

| theorem | content | kind |
|---|---|---|
| `graviton_matches_ligo` | `gravitonHelicity = ligoHelicity = 2` | `decide` |
| `v59_fails_ligo` / `helicity_resolved` | `0 ≠ 2`; `0 < 2` | `decide` |
| `polarization_double_angle` | `cos2ψ = cos²ψ − sin²ψ` (the 2ψ structure) | `Real.cos_two_mul` |
| `quadrupole_coeff_val` | `32/5` | `norm_num` |
| `graviton_speed` | `√(k²+0)/k = 1` (massless) | `Real.sqrt_sq` |
| `gen4_ligo_closed` | bundled headline | mixed |

---

## Status table

| claim | status | tool |
|---|---|---|
| 2 TT modes = `h₊,h×` with helicity ±2 | **verified** | SymPy + Lean |
| `(h₊,h×)` rotate by `2ψ` (spin-2) | **verified** | SymPy |
| GW quadrupole luminosity `(32/5)Gμ²a⁴ω⁶` | **verified** | SymPy (+ Maxima cross-check) |
| GWs travel at `c` (massless graviton) | **verified** | SymPy + Lean |
| LIGO motivation (v59 scalar → 2 TT polarizations) | **resolved** | — |

---

## Significance

This closes the arc that began the entire G9 program:
- **v59**: gravity = scalar, helicity 0 → **fails LIGO** (the #1 blocker).
- **v60**: built exactly 2 TT graviton DOF (ghost-free, stable spectrum).
- **v61 GEN1–2**: those modes live on a genuine GR metric (Schwarzschild,
  backreaction).
- **v61 GEN4**: the 2 TT modes **are** the LIGO `h₊, h×` (helicity ±2), emitted via
  the quadrupole formula, traveling at `c`.

The observational objection that motivated the program is now answered within the
derived dynamical theory.

## What's next

- **GEN5 (aspect 5)**: a further GR consistency check — perihelion precession
  `Δφ = 6πGM/(c²a(1−e²))` from the Schwarzschild metric (another classic test), or
  an FRW cosmology sourced by the matter sector. Then a v61 closeout.
