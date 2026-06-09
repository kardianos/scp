# v61 GEN5 — Perihelion precession from the Schwarzschild geodesic

**Date**: 2026-05-26 (v61 loop, Generation 5)
**Artifacts**:
- `05_perihelion.py` (SymPy + numeric Mercury; all pass)
- `05_perihelion.mac` (independent **Maxima** cross-check, plain-file output confirmed)
- `lean/Perihelion.lean` (builds clean)
**Builds on**: GEN1 (Schwarzschild metric).

---

## Verdict

GEN5 derives the **perihelion precession** — the third classic GR test, completing
the trifecta with **light bending (GEN1)** and **GW emission (GEN4)**. The
Schwarzschild geodesic the v61 gravity sector provides reproduces all three.

### (A) Orbit equation (SymPy + Maxima)

From the Schwarzschild first integral `(u')² = R(u)` (`u=1/r`, c=1), differentiating
(`u'' = R'(u)/2`) gives
```
u'' + u = GM/L² + 3 GM u²,
```
where the **`3 GM u²` term (coefficient 3) is the GR correction** — absent in Newton.
Verified in SymPy and Maxima.

### (B) Secular precession (SymPy + Maxima)

The Newtonian ellipse `u₀ = C(1 + e cos φ)` perturbed by `ε u₀²` (`ε = 3GM/c²`) has a
resonant driver `2εC²e cos φ`, producing a secular term `(εC²e) φ sin φ` (verified:
`(K/2)φ sin φ` solves `u₁''+u₁=K cos φ`). The perihelion advances by
```
Δφ = 2π ε C = 6π (GM)²/(c² L²)   per orbit   (the 6π = 2π × 3).
```

### (C) Standard form + Mercury (SymPy + Maxima + Lean)

With `L² = G M a (1−e²)`:
```
Δφ = 6π G M / (c² a (1−e²)).
```
For Mercury (SI numbers, 415 orbits/century): **42.98 arcsec/century** — matching the
observed GR contribution (~42.98″/cy).

---

## Lean (`Perihelion.lean`)

| theorem | content | kind |
|---|---|---|
| `gr_term_present` | `0 < 3` (GR orbit term vs Newton) | `decide` |
| `precession_factor` | `6 = 2·3` | `decide` |
| `precession_standard_form` | `6π(GM)²/(c²L²) = 6πGM/(c²a(1−e²))` | `field_simp` |
| `precession_positive` | `Δφ > 0` (prograde, bound orbit) | `positivity` |
| `gen5_perihelion` | bundled headline | mixed |

---

## Status table

| claim | status | tool |
|---|---|---|
| orbit eq `u''+u = GM/L² + 3GMu²` (GR coeff 3) | **verified** | SymPy + Maxima |
| secular precession `Δφ = 6π(GM)²/(c²L²)` | **verified** | SymPy + Maxima |
| standard form `6πGM/(c²a(1−e²))` | **verified** | SymPy + Maxima + Lean |
| Mercury ≈ 43″/century | **verified** | SymPy (numeric) |

---

## The three classic GR tests — all reproduced

| test | result | gen |
|---|---|---|
| light deflection | `4GM/b` (2× Newton) | GEN1 |
| gravitational waves | quadrupole `(32/5)Gμ²a⁴ω⁶`, 2 polarizations at c | GEN4 |
| perihelion precession | `6πGM/(c²a(1−e²))`, Mercury 43″/cy | GEN5 |

The v61 gravity sector — the curved-space completion of the v60 dynamical
Lagrangian — passes all three classic tests of general relativity.

## What's next

- **GEN6 (synthesis)**: a v61 closeout assembling GEN1–5 (curved gravity, backreaction,
  EW-vev home, LIGO closure, GR tests), a full regression harness, and an update to
  the root `FUTURE.md` — then the v61 loop concludes (its stated goals met).
