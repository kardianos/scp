# v63 — Trying Both Integration Routes (a clean convergence)

**Date**: 2026-05-28
**Status**: both readings of "integrate according to local density" carried out.
Neither delivers a parameter-free `2/3` — but they converge, with the v59
`PhaseExclusions` theorems, on a single sharp dichotomy that *explains why no such
derivation exists* via integration, and points to the one route that survives.

---

## Reading A — weight by the dynamical local density (`integrate_density.py`)

**A1 (generation-space density).** The Brannen masses `mₖ` are a density over the
`Z₃` generation cycle. The density-weighted **circular mean** of the generation
angle is closed-form:

```
    Σₖ mₖ e^{iθₖ} = a²(6t·e^{−iφ} + 3t²·e^{2iφ}),   θₖ = 2πk/3
    Φ = arg(Σ) = −φ + atan2( 3t² sin3φ ,  6t + 3t² cos3φ )
```

Verified numerically (closed form matches to `1e-12`). At the lepton point
(`t²=1/2`, `φ=2/9`) this gives `Φ = −0.0528 rad`. **Observations:**
- The correction carries `cos(3φ)=cos(2/3)` — the transcendental invariant — so
  density-weighted integration naturally produces a **transcendental** phase
  observable (the right number-type).
- But `φ` entered *via the amplitudes* `mₖ(φ)`. This is a valid phase-**extraction**
  ("integrate to get the form/phase"), **not** a derivation of the value `2/9`.

**A2 (spatial soliton density).** Weighting by the actual soliton profile `ρ(r)`
(v3 data): the dimensionless, scale-free observable `⟨r²⟩/⟨r⟩²` is **1.2514
(massless)** vs **1.2184 (massive)** — a 2.6 % difference. Any dimensionless
spatial density-weighted quantity inherits the soliton's parameter-dependence
(same as `v_min/c = 0.670 vs 0.625`), so it **cannot be a universal constant**.
Same failure mode as the `c_eff` lead.

> Reading A verdict: density-weighting is a legitimate phase-extraction that lands
> in the transcendental number-type, but does **not** derive a parameter-free phase.

---

## Reading B — weight by the canonical measure: DH/character localization (`dh_localization.py`)

The G₂ 7-rep character `χ₇(X) = Σ_weights e^{i⟨μ,X⟩}` (a Duistermaat–Heckman /
Weyl localization: a finite sum over the 6 short roots + 0).

- **Canonical (triality) point** `X=(2π/3,0)`: every weight angle is a rational
  multiple of `2π` (commensurate), so each `e^{i⟨μ,X⟩}` is a root of unity and
  `χ₇ = 2` — an **algebraic** number. No transcendental can appear.
- **The phase `2/3`** needs `⟨μ,X⟩ = 2/3`, i.e. `2/3 ∈ 2π·ℚ`. But `2/3` is **not**
  a rational multiple of `π` (`min|2/3 − 2π·p/r|` over `r≤50` is `~6e-3`, bounded
  away from 0; proven exactly in v59). So `2/3` appears only at a **generic,
  off-lattice** `X` — a free continuous parameter chosen by hand, not a canonical
  element.

> Reading B verdict: this is exactly v59 `PhaseExclusions` E1/E1′ — geometric
> holonomies/characters give π-rational angles; `2/3` is not π-rational, so it is
> "a free Aharonov–Bohm phase, unforced by the algebra." **Pre-closed in v59.**

---

## The unifying dichotomy (both readings + the whole v62→v63 arc)

```
    TRANSCENDENTAL phase 2/3  ⟺  OFF-LATTICE (non-π-rational) angle
                              ⟺  a free CONTINUOUS parameter.
```

Geometry, characters, holonomies, and the canonical/uniform measure **all** hand
you commensurate (π-rational) angles → algebraic cosines. The phase `2/3` is a
**bare radian ratio** (not π-times-rational). That is precisely v59's E1 reading:
*"the phase must be a RATIO (like Q/3) inserted as a radian value, not an angle
read off a rotation."* A bare rational radian is the signature of a **dynamical
loop-ratio** (the only number-type that produces non-π transcendentals), not a
geometric/topological angle (which would be a `π`-multiple).

**Machine-checked** (`lean/PhaseOffLattice.lean`, no `sorry`): the triality angle
`2π/3` is π-rational with `cos = −1/2` (algebraic); the phase argument `2/3` is
**not** π-rational (`phase_arg_not_pi_rational`). Headline `phase_off_lattice`
bundles the dichotomy.

---

## Double-check (before / after)

**Before**: confirmed the relevant objects — Brannen masses as a `Z₃` density;
the G₂ 7-rep weights = 6 short roots + 0; and that v59 already proved
`phase_not_pi_rational`/`koide_not_pi_rational` (read `PhaseExclusions.lean`).

**After**:
- Reading A: closed-form circular mean matches numeric to `1e-12`; spatial
  observable parameter-dependent (`1.2514` vs `1.2184`). 4/4 checks.
- Reading B: `χ₇` at the triality element `= 2` (algebraic), all weight phases
  roots of unity; `2/3` off-lattice (`gap ~6e-3`); character a continuous dial.
  5/5 checks.
- Lean `PhaseOffLattice.lean` builds clean.

---

## Net result

"Integrate according to local density" does not rescue a parameter-free `2/3`:
- the **dynamical** density (Reading A) is parameter-dependent (and only extracts,
  doesn't derive);
- the **canonical** measure (Reading B) gives π-rational/algebraic values, never
  the bare-rational `2/3`.

Both confirm the same thing from opposite ends: **`2/3` is not a geometric angle —
it is a dynamical loop-ratio.** That is not a dead end for the *theory*; it is a
definite answer about the *category* of the phase. The one surviving home is the
continuous, dynamically-pinned loop-ratio `cos(3φ) = −c₃/(4c₆)`
(`v62/residual_audit/phase_dynamical_home.py`). Everything else — c_eff, holonomy,
canonical localization — has now been tried and closed.

**Artifacts**: `integrate_density.py` (4/4), `dh_localization.py` (5/5),
`lean/PhaseOffLattice.lean` (no `sorry`).
