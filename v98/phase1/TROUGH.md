# Phase 1B — the 3+2 fiber and the trough

**Status:** executed 2026-08-13. Paper + CAS. No kernel.
Scripts: `cas/scale.mac` (Maxima), `cas/wind.mac` (Maxima).
Logs: `cas/out/scale.out`, `cas/out/wind.out`.

The charter forbids writing `V(ρ)` by hand. This note asks what
conservation plus “partial cycles convert nothing” actually
force, and records where they do **not** force a unique radius.

---

## 0. Verdict in four lines

1. A 3-D link has a unique transverse 2-plane. That is the
   first fiber (3+2). Winding `W ∈ ℤ` lives there. `W = ½`
   does not close. That is the no-half-proton statement,
   topological, not dynamical.
2. Pure kinetics on that fiber have **no trough**. Scaled
   Dirichlet energy is `R`-invariant (S1). Area cost drives
   `R → 0` (S2). Log-IR + area still has no real critical
   point (S3; Maxima’s formal roots are imaginary). 3-D
   unlocked lumps shrink as `E ~ R` (S4). This is Face C,
   proved.
3. A **discrete** cycle `C_n` plus a site cost `β n` *does*
   have a finite minimizer, but `n*` slides with `ρ = β/α`
   (W3–W5). At `ρ = 1`, even `n` ties `n* = 4` and `6`.
4. Unique size is therefore **not forced** until `ρ` is
   forced. Forcing `ρ` from one energy, without importing
   a well, is the remaining 1B gap. It is not closed here.

---

## 1. The fiber, with no extra stage

A 3-D direction `û` (a link) has a 2-plane `û⊥`. Those are
the chiral planes of `v89/PRINCIPLE.md` §3. They are not a
global `R³ × K`. They exist as the internal geometry of a
conversion and die with it.

| occupation | meaning |
|---|---|
| `W = 0` | φ-independent on the fiber; unlocked travel in 3-space (the slit) |
| `W = ±1` | first lock; ψ single-valued only if the cycle closes |
| `W = ±2, …` | extra wraps; candidates for decaying excitations (`UNFOLD.md` §5) |
| non-integer | not a lock; partial cycle converts nothing |

Proton / neutron, if they live here, are two **occupancies** of
the same `W = ±1` lock, not two windings. The elementary gap
(Higgs sector) is whatever keeps a `W = 0` traveling mode from
falling into `W = ±1`. Those are identity claims. They are not
yet dynamics.

3+4 is still only the named widening: two independent integers
on a 4-plane, if 3+2 cannot host both handedness and the n/p
slot. Nothing in this note forced that widening.

---

## 2. Continuum: no trough from kinetics (Maxima)

Setup: `ψ = f(ρ/R) e^{i W φ}` on the transverse plane.

**S1.** Substitute `ρ = R u`. Both the radial and angular
Dirichlet integrals lose `R`. Fiber `|∇|²` is scale-invariant.

**S2.** Add a space-area cost `E_s = k R²` (space is a mode;
area is space-energy). Then `dE/dR = 2kR > 0`. No interior
critical point. The pinch wants `R = 0`.

**S3.** Grain-fixed core of width `a`, angular IR
`2π W² log(R/a)`, plus area:

```
dE/dR = 2 k R + 2π W² / R
```

Every term is positive for `R > 0`. `solve(dE=0)` returns a
pure-imaginary pair. **No real trough.** The minimum is
`R = a⁺` — the grain.

**S4.** An unlocked 3-D lump `φ_R(x) = φ(x/R)` has
`E ~ R` and shrinks. That is the Face C hoard, as a theorem
of kinetics.

So: a continuum 3+2 field with only gradient energy and a
space-area cost cannot select a composite size. It either
shrinks to nothing or sits on the grain. Sitting on the grain
is unique, and it is the wrong uniqueness (every lock would
be one grain; nuclei would have to be joint copies, which is
Phase 4, but the *nucleon* would have no interior).

A log/area trough was a tempting mistake. The CAS refused it.
That is why this arm used Maxima.

---

## 3. Discrete cycle: a trough exists, and it slides (Maxima)

The fiber, discretised, is `C_n`. Amplitude 1 per site
(cap-shaped, not a new law):

```
E_w(n, W) = n · λ_W = 4 n sin²(π W / n)
```

Checked: `(n,W) = (3,1) → 9`, `(4,1) → 8`, `(6,1) → 6`.
As `n → ∞`, `E_w ~ 4π² / n → 0`. Winding wants to spread.
That is the combinatorial hoard.

Add a site cost — one space-unit per vertex, the cheapest
conservation-geometry term that is not an imported well:

```
E(n) = α E_w(n,1) + β n ,   ρ = β/α
```

Integer argmin on `n = 3..24` (W4):

| ρ | n* (all n) | n* (even only) |
|---|---|---|
| 0 | 24 (the cutoff) | 24 |
| 0.1 | 20 | 20 |
| 0.3 | 11 | 12 |
| 0.5 | 8 | 8 |
| 1 | **5** | **4 and 6, tied** |
| ≥ 2 | 3 | 4 |

At `ρ = 1` the values are almost flat: `E(3)=E(4)=E(6)=12`,
`E(5)≈11.91`. A one-percent tilt picks any of `{3,4,5,6}`.

**W3:** `ρ = 0` has no finite minimizer. That is Face C
again (no space-leg, winding spreads).

**W4:** `ρ > 0` produces a trough. The trough’s *location*
is not forced.

**W6:** `W` is an integer or the cycle does not close. No
stable half-winding. That part is clean.

---

## 4. What would force ρ

`α` is pattern-energy (gradients / winding). `β` is
space-energy per occupied contact. They are one energy, so
in principle `ρ` is a conversion identity, not a knob.

This note does **not** derive that identity. Candidates,
named so they are not invented later as knobs:

- the grain: `β/α` is a pure number built from `2π` and the
  two-atom credit. Then `n*` is an integer near that number.
  Suspicious until the algebra is written.
- occupancy, not geometry: qgraph M3 (`Δ/T` on a *given*
  chart) selects `N`, not `n`. It does not replace ρ.
- 1A’s strict filter: if F-odd+F-3 is adopted, the smallest
  allowed `n` is 6 (`K_{3,3}`, not `C₆`). Combined with
  `ρ = 1` (even, tie 4 and 6) this *picks 6* because 4 is
  not 3-connected. That is a **joint** 1A+1B observation,
  not a derivation of ρ. It is the least-imported way a
  unique smallest robust even lock currently sits at 6.

Writing `ρ` by hand is the same crime as writing `V(ρ)`.

---

## 5. What this does to Phase 1C

- Do not open a 5-D or 7-D cell lattice. The fiber is the
  2-plane of a 3-D link. Continuum kinetics on it do not
  help.
- Do not expect a trough from `|∇|²` plus area. Proved absent.
- A discrete lock plus a space-per-contact cost is the only
  trough this arm found. It needs `ρ` forced, or a filter
  from 1A that cuts the slide down to one `n`.
- p/n as two occupancies of `W = ±1`, Higgs as the `W = 0`
  gap, stay the identity claims. They do not depend on ρ.

Reproduce:

```
maxima --very-quiet -b v98/phase1/cas/scale.mac < /dev/null
maxima --very-quiet -b v98/phase1/cas/wind.mac  < /dev/null
```
