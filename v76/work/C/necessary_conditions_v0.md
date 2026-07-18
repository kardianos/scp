# Necessary Conditions for a Monist Medium (v0)

**Approach:** C — reverse theoretical / relational  
**Date:** 2026-07-18  
**Status:** Round 1 draft — constraints on media, not a dynamics  
**Inputs (phenomenological targets, not GR ontology):**  
local \(c=\mathrm{const}\) in lab frames; \(E\leftrightarrow m\) via \(c^2\); weak-field light deflection \(\propto M/b\); Shapiro-like excess delay  
**Non-circularity rule:** Do **not** assume Einstein’s equation as a premise. Use only operational targets and monist eligibility (`PROBLEM.md` §7). Where GR formulas appear, they are **fit targets** that any medium must approximately reproduce in the weak-field, slow regime — not axioms.

---

## 0. What “necessary” means here

A condition \(N\) is necessary if:

> every monist-eligible medium that realizes local \(c\), shared \(E\leftrightarrow m\), and weak light-bending/delay **must** satisfy \(N\), or else one of those targets fails or dualism re-enters.

These are filters for Approaches A/B/D. They are not yet a unique PDE.

Notation (working):

| Symbol | Meaning |
|--------|---------|
| \(\mathcal{M}\) | medium state (the only continuum) |
| free / bound | partitions of \(\mathcal{M}\) (states, not substances) |
| \(c\) | free-field locality bound in a local frame |
| lock | stable (or long-lived) mass-form region |
| \(E_{\mathrm{bound}}\) | unlock ledger of a lock (rest frame) |
| \(m\) | inertial mass of the lock |
| path cost | excess free-signal travel relative to uniform free medium |

---

## 1. From local \(c = \mathrm{const}\) (C1)

### NC-1.1 — Locality is a free-medium property

\(c\) must be defined as the **supremum of free-update causal speed** in a local frame built from the medium itself (rods/clocks of free or weakly bound oscillators), **not** as a constant of a pre-given metric.

**Fails if:** \(c\) is read off a fixed background chart while locks live “on top.”

### NC-1.2 — Frame calibration identity

There must exist an operational procedure that constructs local orthonormal frames such that free signals are isotropic at speed \(c\), **including near locks** (to leading order in lock strength, outside horizons/traps).

**Implication:** Global charts that hold free density uniform will **not** keep free null cones trivial once locks exist. Warp is the mismatch between “local \(c\) enforcement” and “global uniform-medium charts.”

### NC-1.3 — Rods and clocks are medium-made

Any claim about constant local \(c\) requires that rods/clocks are states of \(\mathcal{M}\). Therefore:

1. Lock formation that changes free budget **can** change nearby clock rates and rod lengths (operationally).
2. “Local \(c=\mathrm{const}\)” is a **normalization** of free signalling against those rods/clocks, not an independent geometric postulate.

**Fails if:** rods/clocks are external classical instruments immune to free-budget change (dualist measurement).

### NC-1.4 — Free radiation is free medium in flight

Packets that propagate at \(c\) are free (or free-carrying) medium configurations, not a second species. Their speed bound is the **same** \(c\) that enters \(E=mc^2\).

### NC-1.5 — Bound regions reduce free participation

Wherever mass-form exists, free degrees of freedom available for causal update are reduced (amplitude, channels, unlocked budget — model-dependent).  
**Not required yet:** a unique formula for how they are reduced.  
**Required:** reduction is dynamical identity of lock-in, not a painted “matter density” on full free continuum.

---

## 2. From \(E \leftrightarrow m\) via \(c^2\) (C2)

Detail note: `reverse_Emc2_bound_v0.md`.

### NC-2.1 — Single ledger for unlock and inertia

For an isolated lock in its rest frame, one number \(E_\star\) must satisfy:

\[
E_{\mathrm{unlock}} = E_\star,
\qquad
m_{\mathrm{inertial}} = \frac{E_\star}{c^2},
\qquad
E_{\mathrm{rest}} = m_{\mathrm{inertial}}\, c^2 = E_\star.
\]

“Bound” means: medium capacity counted in \(E_\star\) is **unavailable** for free propagation until unlock work \(\ge E_\star\) is supplied.

### NC-2.2 — \(c^2\) is forced by locality as the only free velocity scale

Dimensional and causal necessity: conversion between energy ledger and inertial resistance, if mediated only by free-medium response limited by \(c\), must scale as \(E/c^2\). No independent mass unit.

### NC-2.3 — Inertia is relational free-medium resistance

Accelerating a lock must cost free-medium rearrangement work whose low-velocity expansion is \(\frac12 m v^2\) with the **same** \(m = E_\star/c^2\).  
If unlock energy and kinematic inertia use different ledgers, monist \(E=mc^2\) fails.

### NC-2.4 — Formation depletes free budget by \(E_\star\) (integral sense)

Creating a lock of unlock energy \(E_\star\) must remove at least \(E_\star\) of free reconfiguration capacity from the continuum (globally, up to radiation losses accounted in the same ledger).  
Local form \(\rho_{\mathrm{free}}+\rho_{\mathrm{bound}}=\rho_0\) is **sufficient** but **not necessary** (see NC-3.4).

### NC-2.5 — No dualist mass primitive

Mass is not an attribute of a point particle on a stage. It is \(E_\star/c^2\) of a lock. Worldlines, if used, are effective descriptions of lock centers.

---

## 3. From weak-field light deflection and delay (C3)

Detail note: `path_cost_profile_v0.md`.

Phenomenological targets (weak, slow, isolated lock mass \(M = E_\star/c^2\), impact parameter \(b\gg\) lock size):

| Target | Rough form (fit target) |
|--------|-------------------------|
| Light deflection | \(\Delta\theta \sim \mathcal{O}(1)\times \dfrac{G_{\mathrm{eff}} M}{c^2 b}\) with \(\mathcal{O}(1)\sim 4\) preferred (Einstein) |
| Shapiro-like delay | excess coordinate time \(\Delta t \sim \dfrac{G_{\mathrm{eff}} M}{c^3}\ln(\cdots)\) |
| Newtonian limit (orbits) | \(a \sim G_{\mathrm{eff}} M / r^2\) for slow locks (optional Round-1 soft target) |

\(G_{\mathrm{eff}}\) is **emergent** from medium coupling, not a second ontology.

### NC-3.1 — Long-range free-signal path cost \(\propto M/r\) class

Holding **local** free speed \(= c\), global free-signal travel past an isolated lock must acquire an excess path cost whose weak-field leading term is in the **\(M/(c^2 r)\) multipole class** (isotropic monopole of lock strength), so that integrated deflection scales as \(M/b\) and delay as \((M/c^3)\ln\).

### NC-3.2 — Path cost is free-medium structure, not a second field

The excess cost must be a functional of free/bound arrangement of \(\mathcal{M}\) alone. Forbidden as ontology: independent metric field sourced by \(T_{\mu\nu}[\text{matter}]\).

### NC-3.3 — Compact support free-density alone is insufficient (local optics)

If free-signal “index” is a **local** function \(n=n(\rho_{\mathrm{free}})\) only, and free deficit is **compact** (co-located with the lock under local budget \(\rho_f+\rho_b=\rho_0\)), then \(n-1\) has compact support and **cannot** produce long-range \(1/r\) path-cost tails.

**Therefore any monist medium that keeps local algebraic optics and long-range weak lensing must either:**

1. allow **extended** free-medium rearrangement (non-compact free contrast), or  
2. define path cost **non-locally / relationally** (connectivity, constraints, global free-path calculus) so compact locks still induce \(1/r\)-class costs, or  
3. weaken local budget identity (integral / topological / constraint budget).

This is a hard reverse filter on A/B designs.

### NC-3.4 — Integrated free deficit vs \(1/r\) density profile

A free-density contrast \(\delta\rho_{\mathrm{free}}(r)\propto 1/r\) for all large \(r\) has **divergent** spatial integral.  
If the medium insists on

\[
\int (\rho_0 - \rho_{\mathrm{free}})\,dV = E_\star < \infty
\]

with \(\rho_0\) finite, then exterior free-density **cannot** be pure \(1/r\) at infinity.

**Split forced by reverse theory:**

| Object | May be long-range \(1/r\) class? | Must have finite integral \(=E_\star\)? |
|--------|----------------------------------|----------------------------------------|
| Bound lock energy density | no (compact OK) | yes (total \(E_\star\)) |
| Free **energy-density** deficit | only if integrable (e.g. faster than \(1/r^{3+\varepsilon}\) if isotropic, or screened) | yes if local budget |
| Free **path-cost / effective index** contrast | **yes — required for weak lensing** | no — path cost is not the same integral as energy budget |

**NC-3.4 (statement):** Path-cost profile and free-energy-density profile are **not** freely interchangeable. Media that set \(n-1 \propto \delta\rho_{\mathrm{energy}}\) with local budget face the compact-support vs long-range tension (NC-3.3–3.4).

### NC-3.5 — Coefficient \(G_{\mathrm{eff}}\) fixed by lock–free coupling

The map \(E_\star \mapsto\) path-cost amplitude must be universal for isolated locks (equivalence-principle-like): same \(E_\star\) ⇒ same weak field, independent of lock micro-structure to leading order.

### NC-3.6 — Local \(c\) preserved under weak warp

Whatever produces path cost must leave **local** free-signal speed \(= c\) in local frames (NC-1.2). Coordinate speed may vary; local proper free speed may not (outside extreme lock-up).

### NC-3.7 — No second gravity pass

Ray deflection and delay must be readable from free-medium evolution / free-path calculus alone. Scoring a medium by feeding \(E_\star\) into Poisson/Einstein is dualist diagnostics, not monist proof (`PROBLEM.md` §3).

---

## 4. Cross-cutting monist eligibility (from PROBLEM + reverse)

### NC-0.1 — One continuum

No independent stage metric as distance-carrier; null/timelike structure is medium-defined (or operationally equivalent).

### NC-0.2 — One evolution

Formation of locks **is** rearrangement of \(\mathcal{M}\). Not: evolve matter on fixed \(g\), then curve \(g\).

### NC-0.3 — Discriminating observables

At least one observable ties lensing/inertia to free-medium deficit or free-path structure, not only to \(\int T_{00}\) on flat space.

### NC-0.4 — Particle stability is lock stability

Stability criteria refer to medium lock-in, not only VK of a field on \(\mathbb{R}^3\).

---

## 5. Minimal necessary-condition checklist (short list)

Any candidate medium must:

1. Define \(c\) from free updates; calibrate local frames so free signals see that \(c\) (NC-1.1–1.3).  
2. Identify unlock energy and inertial mass via one ledger \(E_\star\) with \(m=E_\star/c^2\) (NC-2.1–2.3).  
3. Ensure lock formation removes free capacity of order \(E_\star\) (NC-2.4).  
4. Induce free-signal path costs in the weak-field \(M/(c^2 r)\) class without a second gravity sector (NC-3.1–3.2, 3.7).  
5. Resolve the **local-optics vs long-range** tension explicitly (NC-3.3–3.4) — pick extended free rearrangement, nonlocal path cost, or non-local budget.  
6. Keep local free speed \(c\) under weak warp (NC-3.6).  
7. Make \(G_{\mathrm{eff}}\) emergent and universal in \(E_\star\) (NC-3.5).  
8. Strip dualist stage/actor structures (see `dualist_strip_list_v0.md`).

---

## 6. Explicitly not claimed (Round 1)

- Uniqueness of medium PDE.  
- Derivation of full Einstein equation.  
- Galactic DM residuals (C4 deferred).  
- Horizon microphysics.  
- That \(\rho_f+\rho_b=\rho_0\) is mandatory locally.

---

## 7. Pointers

| Artifact | Role |
|----------|------|
| `path_cost_profile_v0.md` | Forced path-cost / free-budget sketches (C3) |
| `reverse_Emc2_bound_v0.md` | What “bound” must mean (C2) |
| `dualist_strip_list_v0.md` | Structures to strip (C5) |
| Log `C_reverse_theory.log` | C-001+ provenance |

---

## 8. Revision rule

Strengthen only by: (i) sharper operational targets, (ii) no-go lemmas (as NC-3.3–3.4), (iii) kill results from B/D. Never by assuming \(G_{\mu\nu}=8\pi G T_{\mu\nu}\) as monist premise.
