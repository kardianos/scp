# Integration plan: put v58 dynamics under the v59 dependency relations

*2026-05-24.  Goal: replace the conjectural/fitted dependency relations (v_Higgs=28²a_l²,
g_W²=5√α, α, …) with **physical** relations derived from the v58 field equation + its vacuum
constraint, using the v58↔v59 dimension map.  Status: a research program; the first reframing
(v_Higgs) is concrete and tested, the rest are targets.*

## The linchpin identification

The v58 density `ρ_M = ½(M M̃ − v²)` defines a **Higgs-like vacuum manifold** `M M̃ = v²` with
"vacuum expectation v."  **Identify `v = v_Higgs`** (the electroweak VEV).  Then the entire v59
dependency stack is a statement about *excitations on the EW vacuum manifold*, and the dim-map
(even/ℍ sector; `L=Λ²⊕Λ⁶` = the complex/chiral grades = the mass-bearing directions, proven in
`lean/BladeSquareSign`+`LeptonRealityForcing`) tells us *which* directions carry what.

## (1) v_Higgs — reframed into a physical mode-count (the main win here)

`v_Higgs = 28²·a_l²` is dimensionally `v_Higgs = (dim L · a_l)²`, i.e.

> **√v_Higgs = dim(L) · a_l**,   equivalently   **a_l = √v_Higgs / dim(L)**   (verified, 0.03%).

Physical reading: the lepton Brannen amplitude is the **EW-vacuum scale per L-direction** —
`√v` shared (equipartitioned) over the `dim(L)=28` mass-bearing (chiral/complex) directions.
This is no longer a fitted "784"; it is a mode-count over exactly the grade we *proved* carries
the complex structure / mass (`L=Λ²⊕Λ⁶`).

**It is lepton-specific** (tested): the naive `a_q = √v/D_q` fails badly (`a_d` off 1.8×, `a_u`
off 19×).  This *matches* the picture — the lepton is the **color singlet / the ℍ-sector "density
achiever"** of v58 (`PARTICLES_AS_DENSITY_ACHIEVERS.md`); the EW vacuum couples to it, not to the
colored quarks the same way.

**v58 derivation target.**  Show the v58 dynamics on `M M̃ = v²` give a lepton (color-singlet)
excitation amplitude `√v / dim(L)`:
- (a) the vacuum sits at `M M̃ = v² = v_Higgs²` — set by the equation's `μ⟨Ω,Ω⟩` + `v` (the
  Higgs-like minimum, the analog of v59's proven `XiVacuum` |ξ|²=1/2);
- (b) mass excitations are confined to the `L`-grade (28 directions) — **already proven** (v59:
  the complex structure / mass term lives in `L`, not `F`);
- (c) equipartition of `√v` over those 28 directions gives the per-mode amplitude `a_l = √v/28`.
(a) and (c) are the v58-dynamical pieces to derive; (b) is in hand.  If (a)+(c) follow from the
field equation, **v_Higgs=28²a_l² is derived**, and the "1 input" headline becomes earned:
`a_l` and `v_Higgs` are the *same* scale (one dimensionful input) tied by the mode-count `dim L`.

> **CORRECTION (2026-05-24, `integration_v58/03_higgs_bridge_result.md`).**  Step (c) as
> written is **wrong**: the Cl(7) computation shows the 28 `L`-blades are *Frobenius-orthogonal*,
> so equipartition of `√v` over them (vacuum-as-a-vector-in-`L`) gives `√v = √28·a_l ≈ 5.3 a_l`
> under the physical energy norm — **not** `28 a_l`.  The factor `28` is correct only under the
> **bilinear** reading: the lepton Yukawa is a *matrix* on the 28-dim `L`, with `dim(L)²=784`
> components, and `v_Higgs = ‖Y‖²_Frobenius = dim(L)²·a_l²` with `so(8)`-democratic entries
> (natural `L²` norm).  See `03_higgs_bridge_result.md` for the computed result and the
> remaining residual.

## (2) g_W² = 5√α — reframe via the v58 bivector connection

In v58 the gauge field is the **bivector (grade-2 = chiral/EM) connection** `Ω`; couplings are its
normalizations.  Reframe:
- `5 = h∨(Spin7)` (dual Coxeter, **proven**) = the connection's β-function / Casimir
  normalization for the `SU(2)_L` sub-bundle (`so(3) ⊂ so(7)`);
- `√α`: the EM (bivector) coupling; the **square-root** is the signature of the v58 **quadratic**
  self-interaction `λΩ² + μ⟨Ω,Ω⟩` (a coupling that enters EM linearly enters the squared/√ form
  in the cross-channel).

**v58 target.**  Derive `g_W²` as (dual-Coxeter normalization `5`) × (EM bivector coupling, whose
`√α` form comes from the `Ω²` term).  This would turn `g_W²=5√α` from "ansatz with no Lagrangian"
into a projection of the v58 equation; and since `α(M_Z)=25/(324π²)` is *derived from* `g_W²=5√α`
+ `sin²θ_W=2/9` (proven) + the SM relation, **closing `g_W²=5√α` closes `α(M_Z)` too.**

## (3) α — the v58 EM/bivector self-coupling

`α` is the strength of the bivector (chiral/EM) sector.  In v58 that is the coefficient of the
self-interaction in the EM channel.  The α(0)/α(M_Z) conjectures would become the **fixed-point /
normalization condition of the bivector self-coupling** on the vacuum manifold.  Lowest-priority,
most speculative — but it is the only dimensionless input, so deriving it is the endgame.

## The dimension map (the tool)

| v59 integer | what it is | v58-physical role |
|---|---|---|
| `dim L = 28` | `Λ²⊕Λ⁶` = complex/chiral grades (mass-bearing) | # directions the EW vacuum equipartitions over → `√v=28 a_l` |
| `5 = h∨(Spin7)` | dual Coxeter | normalization of the `SU(2)_L` bivector connection → `g_W²=5√α` |
| `dimG₂=14, dimSpin7=21` | `Q = 14/21 = 2/3` | Koide ratio (the lepton amplitude shape) |
| `D_d=35, D_u=63` | `F`, `L⊕F` ambient dims | quark sectors (the looser, scale-dependent Koide) |
| `√v = v_Higgs` | v58 vacuum norm | the EW vacuum manifold `M M̃ = v²` |

## Honest status & order of work

This is a **program**, not a result.  Each item above is a *target* — a way the v58 equation
*could* derive a relation currently held as conjecture.  What is genuinely new and solid here:
- the **`v = v_Higgs` identification** (gives v58's `v` a physical meaning), and
- the **`√v = dim(L)·a_l` reframing** — the worst offender (`28²`) is now a physical, tested,
  lepton-specific mode-count over exactly the proven mass-bearing grade.

Order (load-bearing first, mirrors `RIGOR_AUDIT.md`):
1. **v_Higgs / the equipartition** — derive (a)+(c) above from the v58 equation (the vacuum
   minimum and the `√v/dim L` mode-amplitude).  Highest value: earns the "1 input" claim.
2. **g_W²=5√α** — derive the bivector-connection normalization; this also closes `α(M_Z)`.
3. **α** — the bivector self-coupling fixed point.
4. quark Koide / selection rule — lowest (soft, RG-dependent); treat as pattern.

The caveat from `RIGOR_AUDIT.md` stands until these are done: the "≈1 parameter" claim is
conditional.  But the v58 integration is the right vehicle to make it unconditional, because it
supplies the one thing the v59 numerology lacks — **a dynamical law on a vacuum manifold** — and
the dim-map tells that law which grades do what.
