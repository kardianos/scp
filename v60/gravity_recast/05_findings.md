# G9 (C) — soldered-tetrad route: exact DOF, weak-field bridge, naturalness

**Date**: 2026-05-25 (advance on scenario C of `04_findings.md`)
**Artifacts**:
- `05_dof_and_weakfield.py` (runs clean; all assertions pass)
- `../lean/G9Soldering.lean` §3b (machine-checked DOF identities; axiom-clean)
**Builds on**: `04_findings.md` (soldering = Minkowski sum; carrier belongs in `Cl(3,1)`)

This pushes the soldered-tetrad route from "±2 is *present*" (04) to "the route gives
*exactly* the 2 LIGO modes **and** contains the v59 scalar law," leaving one
well-posed open item (naturalness).

---

## (C1) Exact degree-of-freedom count = 2 — DONE

The soldered 2-form's metric fluctuation `h_μν` is a symmetric rank-2 (10 comps).
For a massless graviton along `z`, computed entirely by rank/nullspace
(`05_dof_and_weakfield.py`, no number hand-set):

```
symmetric h_μν components ......................... 10
after transverse (4) + traceless (1) constraints .. 5
residual gauge directions preserving TT ........... 3      ← incl. the e₀+e₃ direction
PHYSICAL DOF = 5 − 3 = 2
helicities of the 2 physical modes = [−2, +2]
```

So the count is **exactly 2**, helicity ±2 — the LIGO modes `h_+, h_×`, with **no
extra scalar/vector/ghost** modes. Equivalent covariant count `10 − 4 − 4 = 2`.
Machine-checked integer identities and the helicity multiset (`card = 2`,
`= {±2}`) are in `G9Soldering.lean §3b` (`graviton_dof_covariant`,
`graviton_dof_TT`, `graviton_helicities`), axiom-clean.

*(Subtlety found & fixed: a naive count gave 3 because, for null `k`, the
transverse-preserving gauge directions are the 3-dim space `{k·ξ=0}` — which
includes the `e₀+e₃` combination, not just the two coordinate `ξ`'s. Counting all
three recovers the correct 2.)*

---

## (C2) Weak-field bridge: the v59 scalar law is the trace sector — DONE

Linearized Einstein in harmonic gauge, `□ h̄_{μν} = −16πG T_{μν}`. The static 00
component with `h_{00} = −2Φ`, `T_{00} = ρ` gives `∇²Φ = 4πG ρ`, i.e.
`□Φ = −4πG ρ`. Identifying

| v59 object | linearized-GR object |
|---|---|
| `Ω_grav` | trace/Newtonian potential `Φ` |
| `f_g` | `−4πG` (gravitational coupling) |
| `ρ_grav = Tr(M†M) = Σ m_k` | `T_{00}` (second-moment source) |

the established v59 law `□Ω = f_g ρ_grav` is **exactly the trace (helicity-0,
Newtonian) sector** of the one tensor equation, whose **transverse-traceless sector
is the 2 graviton modes of (C1)**. So moving the carrier to the `Cl(3,1)` soldered
2-form **contains the v59 scalar result** (as its trace) and **adds** the LIGO ±2
modes — nothing is lost, and the radial/magnitude results sit in the trace sector
unchanged at this level.

---

## (C3) Naturalness — UPGRADED from "plausible" to "derived (mod one v59 gap)"

The resolution requires the soldering `so(3,1)` to **commute with** the full
internal `Spin(7)` (which must stay unbroken to keep lepton = L, Koide, gauge).
This is now a *derivation*, not a hope (`06_lorentz_commutant.py`, runs clean):

**DERIVED (Schur, verified with v59's octonion structure).** `Spin(7)` is
irreducible on the octonion 8, so its commutant inside `End(ℝ⁸) = Cl(7)_even` is
**1-dimensional (scalars)** — computed exactly (commutant dim = 1, via the 21
explicit `so(7)` generators `L_iL_j`). Since `dim so(3,1) = 6 > 1`, **no Lorentz
`so(3,1)` can commute with the internal symmetry inside the algebra v59 used for
gravity** (`Λ²(V⁸)=so(8) ⊂ End(V⁸)`). So:
- the v59 carrier had **no room** for a Lorentz structure → forced scalar
  (independent confirmation of the no-go), and
- the spin-2 carrier is **forced out** of the internal arena into a tensor
  extension — exactly the spacetime `ℂ⊗ℍ / Cl(3,1)` factor of 04/05. This is *why*
  the internal-soldering routes 01–03 had to fail, now rigorous.

Machine-checked consequence: `G9Soldering.no_internal_lorentz`
(`commutantDimInternal = 1 < 6 = dimSO31`), axiom-free.

**VERIFIED (constructive).** In the tensor extension (spacetime⊗octonion, 32-dim),
`so(3,1)` on the spacetime factor commutes **exactly** with the full internal
`Spin(7)` on the octonion factor (`max|[·,·]| = 0`), and acts trivially on the
octonion sector (`max|[so(3,1)⊗I, I⊗(any octonion op)]| = 0`) — so `G₂`/triality/
color literally cannot see it. The `so(3,1)` is genuinely Lorentzian (non-compact
boosts: real eigenvalues ±1; rotations imaginary ±i; preserves `η`). Compatibility
is therefore **automatic**, and the A2/B1 worry does not arise for the `Cl(3,1)`
carrier.

**RESIDUAL (honest — a v59 gap, not a flaw in this argument).** That this tensor
extension *is* v59's single intended algebra `ℂ⊗ℍ⊗𝕆`, with a unique/forced
embedding, is exactly the spacetime–internal integration that v59's own
`SevenDAlgebra.lean` header lists as "a future refinement." Pinning that map (an
explicit `Cl(3,1) ↪` the full algebra commuting with the realized `Spin(7)`)
closes C3 completely. Until then, C3 is **derived modulo v59's unfinished merge** —
a substantial upgrade from "plausible."

---

## Status of (C)

| item | status |
|---|---|
| ±2 present iff soldered | **proved** (04; `G9Soldering`) |
| exactly 2 TT DOF, no extras | **proved** (C1; rank computation + Lean) |
| recovers v59 `□Ω=f_g ρ` as trace sector | **shown** (C2; symbolic) |
| no internal Lorentz (Schur) ⇒ carrier forced to spacetime factor | **derived** (C3 Q1; commutant=1, Lean `no_internal_lorentz`) |
| `so(3,1)`↔internal commutation in tensor extension | **verified** (C3 Q2; `max|[·,·]|=0`, octonion-inert) |
| the extension IS v59's `Cl(spacetime)⊗Cl(internal)` | **resolved** (06–07; v59 SYNTHESIS §8 — see `07_findings.md`) |
| full Plebański action + EOM = OBE | open (next) |
| magnitude `α²¹` in emergent metric | open (sits in trace sector; likely unchanged) |

---

## Other parts we can move forward (menu)

1. **(C3) embedding + commutation** — the natural next domino for gravity; needs the
   `ℂ⊗ℍ⊗𝕆` embedding choice, then a decidable/numeric `[Lorentz, color]=0` check.
2. **Full Plebański action** — write `S[B]` for the `Cl(3,1)` 2-form + simplicity
   constraint, derive EOM, show weak-field limit = (C1)+(C2). Builds on C3.
3. **G1 (rank tension)** — entirely unstarted. Self-contained first step: classify
   the maximal-subalgebra chains of `so(8)` and test whether any leaves a rank-3
   light sector while preserving the 784 Frobenius scale (numeric + rep theory).
4. **Magnitude** — verify `G_e=(21/16)α²¹` survives as the trace-sector coupling of
   the emergent metric (should be inherited from the v59 scalar law via C2).

Items 1–2 deepen gravity; item 3 opens the second blocker; item 4 is a consistency
check. 1 and 3 are the most independent and highest-value.
