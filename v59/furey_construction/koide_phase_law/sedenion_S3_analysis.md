# Analysis in the sedenion S₃ factor: the "/3" is emergent, the "Q" is residual

*Proceeding past the octonion `G₂` layer into the sedenion `S₃` factor of
`Aut(𝕊) = Aut(𝕆) × S₃ = G₂ × S₃` (Gresnigt et al. 2024; Eakin–Sathaye).  Computation in
`sedenion_s3.py`; the π-rationality clincher is machine-checked
(`lean/PhaseExclusions.koide_not_pi_rational`).*

## Why move to the sedenions

Our holonomy test failed *inside* the octonions because the generation `Z₃` is not an
octonion automorphism — `Aut(𝕆) = G₂` has no `S₃` factor.  The literature locates the
generation symmetry one Cayley–Dickson step up: in the **sedenions** `𝕊`,
`Aut(𝕊) = G₂ × S₃`, and the `S₃` is "the algebraic source for exactly three generations."
So the right place to test the `/3` is the `S₃` factor.

## What was built and verified (all computed in `sedenion_s3.py`)

The sedenions `𝕊` (16-dim, Cayley–Dickson `(a,b)(c,d) = (ac − d̄b, da + bc̄)`) and the explicit
`S₃` generators:

- `ε : A + Bs₈ ↦ A − Bs₈`  (order 2),
- `ψ : A + Bs₈ ↦ ¼[A + 3A* + √3(B − B*)] + ¼[B + 3B* − √3(A − A*)]s₈`  (order 3).

**Verified:** `ε² = I`, `ψ³ = I`, `εψ = ψ²ε` (so `⟨ε,ψ⟩ ≅ S₃`), and — checked on all `16×16`
basis products — **both `ε` and `ψ` are genuine automorphisms** `f(xy) = f(x)f(y)`.

`ψ` acts as a **120° rotation** `[[−½, √3/2],[−√3/2, −½]]` of the two 7-dim imaginary-octonion
blocks `(Im A, Im B)` (the `√3 = sin 120°`), fixing the real parts `e₀, e₈`.  So `⟨ψ⟩` realises
the `Z₃` Fourier split `1, ω, ω²` — exactly the Brannen generation cycle `S`.

## The two-part verdict

**(1) The "/3" is emergent — confirmed and now rigorous.**  The generation `Z₃ = ⟨ψ⟩` is a
*genuine sedenion automorphism*, not an external add-on.  So the divisor in `φ = Q/3` — the
number of generations — **is structural, emergent from (extended) octo-space**, precisely the
hypothesis.  This is real progress: the `/3` is no longer assumed, it is the `S₃` of `Aut(𝕊)`.

**(2) The "Q" (the phase magnitude) is NOT fixed by the symmetry — residual.**
- `ψ`'s *intrinsic* phases are the cube roots `{0, ±2π/3}` — **π-rational**.
- The Brannen phase is the *coupling* phase `φ = arg(ξ)` in `M = a(I + ξψ + ξ̄ψ²)`, **not** an
  eigenphase of `ψ`.  `ψ`-covariance only forces `M` to be circulant — **any `ξ` is allowed**.
- `φ = 2/9` is **not π-rational** (`koide_not_pi_rational`), so it is not *any* symmetry or
  holonomy phase.  Hence the `S₃` does **not** fix `φ = 2/9`.

## The clean decomposition this buys

`φ = Q / 3` now splits into an explained part and a residual:

| piece | meaning | status |
|---|---|---|
| `/3` | `N_generations` = the sedenion `S₃` automorphism factor | **emergent / structural ✓** |
| `Q = 2/3` | the phase *magnitude* `3φ` = the coupling value | **residual** — not symmetry-fixed; not π-rational; orthogonal to the generation symmetry |

So moving to the `S₃` factor **confirms the `/3` is emergent** and **isolates the entire
residual mystery to one sharp question**: *why is the coupling magnitude exactly
`3φ = Q = dimG₂/dimSpin7 = 2/3`?*  That is a mass-sector/dynamical input, decoupled from the
(now-explained) generation count.  The `S₃` symmetry gives the *structure* of three
generations; it does not — and provably cannot, by π-rationality — give the *value* of the
generation phase.

## Net

- **Answered (the user's hypothesis):** yes, the `/3` is emergent from octo-space — precisely,
  from the sedenion `S₃ ⊂ Aut(𝕊)`, verified as a genuine automorphism here.
- **Residual:** `Q = 2/3` (the phase magnitude `3φ`) is a coupling value, not a symmetry phase,
  and provably not a holonomy (π-rational obstruction).  The `φ=Q/3` law is thereby reduced to:
  *(generations from `S₃`) × (a single non-symmetry magnitude `Q`)*.
- This is the cleanest available statement of where the derivation stands: the count is
  structural, the magnitude is not — and the magnitude's only structural target, the
  transcendental `cos(2/3)`, matches no candidate (`PhaseAmbiguity`).
