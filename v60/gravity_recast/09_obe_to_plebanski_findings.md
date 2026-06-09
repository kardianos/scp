# Can the Plebański action be DERIVED from the OBE? — the precise obstruction

**Date**: 2026-05-25 (the open item left by `08_plebanski_action.md` §5: "action
**derived from the OBE** — open")
**Artifacts**:
- `09_obe_to_plebanski.py` (runs clean; all assertions pass)
- `../lean/G9ObeToPlebanski.lean` (compiles clean against v59 Mathlib; the
  obstruction theorems `dof_strict_increase`, `fields_strict_increase`,
  `no_content_preserving_derivation`, `strict_structure_increase` depend on **no
  axioms at all** — pure `decide`)
**Builds on**: `04`–`08` (soldering → 2 DOF → carrier forced into `Cl(3,1)` →
unified algebra → posited Plebański action + verified induced metric).

---

## Verdict (up front, honest)

**NOT DERIVABLE.** The Plebański 2-form action of `08` is an **independent posit**,
not a consequence of the OBE. The OBE supplies its trace/Newtonian *sub*sector (a
consistency anchor that the posited action was engineered to match), but it
provides **neither** of the two ingredients that make Plebański a *gravity* theory:
the **BF kinetic term** `B∧F` (no independent connection in the OBE) and the
**simplicity constraint** (no quadratic spacetime constraint, no rank‑4 Weyl
multiplier in the OBE).

The implication that actually holds runs the **opposite** way:

> **Plebański action** ──(linearize, take the 00/trace component)──▶ **OBE scalar law**
> `□Ω_grav = f_g ρ_grav`.

So the OBE is a *consequence* (the trace subsector) of the posited action, **not its
origin**. This is exactly what `08`'s own status table already flagged as open; `09`
turns "open" into a *characterized* negative result with a precise obstruction.

---

## What the OBE actually is, vs. what Plebański is

The single decisive fact is a **mismatch in dynamical content**, made into integers
and machine-checked.

| | OBE gravity sector | Plebański action |
|---|---|---|
| independent field(s) | **1** — `Ω_grav`, a scalar potential | `B`(36) + `ω`(24) + `Ψ` (multiplier) |
| equation order | **2nd** (`□Ω = f_g ρ`) | **1st** (`B∧F`, `F=dω+ω∧ω`) |
| propagating DOF | **1** (helicity 0, Newtonian/trace) | **2** (helicity ±2, TT graviton) |
| connection | **none** (`Ω` is *slaved*: `Ω=f_g∇(K∗ρ)`) | independent `ω`, varied freely |
| pairing | convolution kernel `K(x,x′)` — **metric-dependent** | wedge `B∧F` — **metric-INDEPENDENT** |
| metric | **fixed background** (built into `K`, into `□`) | **induced** `g(B)` via Urbantke |

The OBE is the **already-integrated, connection-eliminated, second-order scalar**
form. Plebański is a **first-order two-field tensor** action that *produces* the
metric. The target has strictly more dynamical content — you cannot manufacture two
propagating tensor DOF and an emergent metric out of one slaved scalar mode.

`G9ObeToPlebanski.lean` machine-checks the integer backbone, axiom-free:
`obeDOF = 1 < 2 = plebDOF` (`dof_strict_increase`), `obeFields = 1 < 60 =
plebFields` (`fields_strict_increase`), and the resulting obstruction
`¬(plebDOF ≤ obeDOF)` (`no_content_preserving_derivation`).

---

## Ingredient-by-ingredient (the body of `09_*.py`)

### (1) The BF kinetic term `B∧F` — **ABSENT from the OBE**

`B∧F` has three irreducible parts, none present in the OBE as independent
structure:

- **(i) an independent connection `ω`.** The OBE connection `Ω_grav` is *not*
  independent: by `NEW_OBE_FORMULATION §3a` it is the solved convolution
  `Ω = f_g ∇(K∗ρ)`, algebraically slaved to the source through a **fixed** kernel.
  There is no free `ω` to vary; `F[ω]` has no OBE analogue.
- **(ii) a metric-independent wedge.** The OBE's pairing is the Green's-function
  kernel `K`, which is built from `□` on a **fixed** background metric — the
  *opposite* of the topological, metric-free wedge `B∧F`.
- **(iii) first-order EOM** `d_ω B = 0` and `F = Ψ·B + source`. The OBE is
  second-order (`□Ω = ρ`). A first-order BF system *reduces to* a second-order wave
  equation only **after** eliminating `ω`; the OBE starts already-eliminated.

The 1D toy in `09_*.py` makes the irreversibility explicit: the first-order pair
`(φ′=p, p′=s)` (two fields) eliminates to `φ″=s` (one field, residual 4.7e−2,
consistent), but the inverse map `φ″=s ↦ (φ,p)` is **not unique** — the connection
`p` is an extra datum the scalar equation does not carry. Re-introducing `ω`, the
wedge, and the multiplier is **added** structure, not derived.

### (2) The simplicity constraint `Ψ` — **ABSENT from the OBE**

Simplicity (`B^{IJ}∧B^{KL} ∝ ε^{IJKL}`, enforced by the Weyl-symmetric multiplier
`Ψ_{IJKL}`) is the part that makes a metric *emerge* from `B` (verified in `08`).
The only OBE object that could conceivably play `Ψ`'s role is the v59 selection-rule
projector `Π_N` (`§4`). It fails on all three counts:

- **(a) order.** `Π_N` is a **linear** idempotent projector (`Π² = Π`, verified);
  simplicity is **quadratic** in `B`. Different algebraic order.
- **(b) sector.** `Π_N` acts on the **internal** `Cl(7)_even` (L/F grades);
  simplicity acts on the **spacetime** `Cl(3,1)` 2‑form `B`. By `06`/`07` these are
  *commuting* tensor factors, so an internal projector is **blind** to the spacetime
  2-form — it cannot impose `B∧B ∝ ε`.
- **(c) symmetry.** The OBE source `ρ_grav = Tr(M†M)` is a single Lorentz
  **scalar** (verified: a real number). `Ψ_{IJKL}` is a rank‑4 Weyl-symmetric
  tensor. A scalar cannot generate / be a rank‑4 multiplier; **no OBE object has
  `Ψ`'s symmetry**.

So the simplicity constraint — the keystone of the induced-metric mechanism — has
**no OBE origin**.

### (3) The source `f_g ρ_grav` — **PRESENT, but only as a consistency check, and in the wrong logical direction**

This is the one leg that closes. The linearized Einstein 00/trace component
(harmonic gauge `□h̄_{μν} = −16πG T_{μν}`, static `h_00 = −2Φ`, `T_00 = ρ`) gives
`∇²Φ = 4πG ρ`, i.e. `□Φ = −4πG ρ`. With `Ω_grav ↔ Φ`, `f_g ↔ −4πG`, `ρ_grav ↔
T_00`, the OBE scalar law **is** the trace sector — exactly `05` (C2), re-verified
on numbers in `§4`.

**But this is `Plebański ⟹ OBE`, not `OBE ⟹ Plebański`.** The action's *trace
sector* reproduces the scalar law; matching a necessary trace is a consistency
check the posit passes, **not** a derivation of the posit. (And the `α²¹` magnitude
sits in this same trace coupling `f_g`, where it remains a value-conjecture, not
derived — unchanged from the v59 status.)

---

## The precise obstruction (one sentence)

> The OBE is an **integrated, connection-slaved, second-order scalar** equation on a
> **fixed** background metric, carrying **1 propagating (helicity-0) DOF**; the
> Plebański action is a **first-order, two-field (`B`,`ω`) tensor** theory that
> **produces** the metric and carries **2 propagating (helicity ±2) DOF** plus a
> quadratic **simplicity constraint** with a **rank-4 Weyl multiplier** — none of
> which (independent connection, metric-free wedge, quadratic constraint, Weyl
> multiplier) exists in the OBE. The target has strictly more dynamical content, so
> no content-preserving map `OBE ⟹ Plebański` can exist; only the trace-subsector
> map `Plebański ⟹ OBE` holds.

This is the same shape of honest result the project has reached elsewhere (e.g.
"`g_top` is FREE at all levels tested"; "`α` is a genuine value-conjecture, not
derivable"): a structure that *matches* on the sector where it can be checked, but
is **posited** where it adds content beyond what the fundamental equation supplies.

---

## Status table

| claim | status |
|---|---|
| trace source `□Ω = f_g ρ_grav` = Plebański's 00/trace sector | **verified** (`§4`; = `05` C2) — direction is Plebański⇒OBE |
| OBE contains a `B∧F` kinetic term | **NO** (`§2`) — OBE is the `ω`-eliminated 2nd-order scalar; obstruction is genuine |
| OBE forces / contains the simplicity constraint `Ψ` | **NO** (`§3`) — `Π_N` is linear, internal, scalar-sourced; no rank-4 Weyl object |
| OBE source carries `Ψ`'s Weyl symmetry | **NO** (`§3c`) — `ρ_grav` is a Lorentz scalar |
| content-preserving derivation `OBE ⟹ Plebański` | **impossible** — `dof_strict_increase` / `no_content_preserving_derivation` (Lean, axiom-free) |
| Plebański action **derived from the OBE** | **NO — independent posit** (the headline answer to the `08` open item) |
| field/DOF/component integer backbone of the obstruction | **machine-checked** (`G9ObeToPlebanski.lean`, axiom-clean) |

---

## What this does and does not change for the G9 arc

- It does **not** invalidate `04`–`08`. Those establish that *if* you adopt a
  spin-2 carrier, it is forced into the `Cl(3,1)` factor, gives exactly 2 TT DOF,
  commutes with the internal sector, and (via the posited Plebański action) induces
  a metric whose trace reproduces the OBE scalar law. All of that stands.
- It **does** correct any reading that the Plebański action is *implied by* the OBE.
  It is not. The OBE under-determines the gravity sector: it fixes the scalar/trace
  law and the source, but the **tensor (BF) structure and the simplicity constraint
  are additional postulates**, chosen (legitimately, as the minimal Lorentz-covariant
  completion) but not forced by the OBE.
- The honest closing statement for the gravity sector: **the OBE fixes the
  Newtonian/trace law and the gravitational charge `ρ_grav`; promoting that scalar
  law to a full tensor (LIGO-compatible, metric-producing) gravity theory requires
  the Plebański BF + simplicity structure as an independent input.** Whether a
  *more fundamental* OBE (e.g. a genuinely first-order multivector action whose
  Euler–Lagrange equations are the OBE force law, with an independent `Cl(3,1)`
  connection) could contain `B∧F` is a well-posed next question — but it would
  require writing that action first, which is precisely the "construct the
  Lagrangian whose EL equations yield `Ω(x)`" step that `NEW_OBE_FORMULATION §5`
  itself lists as not-yet-done.
