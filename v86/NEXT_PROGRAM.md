# v86 — what to run next, after Part 0
**Date:** 2026-07-26 · Triage of the unrun census rungs, the new work Part 0
forces, and the proposed proof architecture. Companion to `PART0_RESULTS.md`.

---

## Part 1 — Census triage: keep, redesign, kill

### KILL outright

**HC-4 (line-width spectroscopy) — do not run, in any form resembling its
design.** Three independent measurements each kill it:

1. **There is no line to measure where the program works.** HC-1's
   box-converged solve finds zero bound internal modes for ω ≥ 1.36. The
   entire production era lives there. A width measurement needs a narrow
   line; there isn't one.
2. **The pre-registered observable tests the wrong power.** "Width ∝
   amplitude² = golden-rule confirmation" assumes direct (n = 1) emission.
   Every monopole candidate in the catalog has its fundamental below Ω_c, so
   the leading process is n ≥ 2 and the exponent differs per mode.
3. **The target lines were probably never linear modes.** The QRK-1 lines
   (0.018–0.126) do not appear anywhere in the linear catalog at the ω where
   they were seen, and they sit below the box IR cutoff — the council's own
   box-IR hazard (GROUNDING §2 fix #4). They are candidates for nonlinear
   frequency-pulling or box artifacts.

The halving list's demotion of HC-4 was economical; it is now *earned*. What
survives is a single cheap confirmation, not a campaign — see HC-4-lite below.

**HC-5 (overload test) — already merged into HC-6; keep it dead.** Seeding
2-frequency content in one component tests decay-to-sector-minimum, which HC-6
tests better with GSS-unstable partitions. No separate rung.

**EX-3 (full boost-invariance audit) — stays cut.** The 1–5% lattice
group-velocity anomaly already bounds "is movement simulated correctly", and
N7's dx pair now measures the same systematic on a *dynamical* observable for
free. A dedicated lab-vs-boosted dx set is expensive certification, not
discovery.

### KEEP, unchanged

**HC-6 (converse decay) — but gated on a cheap prerequisite.** HC-3 found
n(D) = 1 at every scanned partition, so there are no decay targets on that ray.
Worse, HC-1 showed analytically that L_x^flav = L₀ − A with A = 2P₀ < 0, so the
flavour channels can *never* contribute a negative direction in this potential.
It is therefore a live possibility that **the flavoured branch has no
GSS-unstable region at all**, in which case HC-6 is unrunnable — and that is
itself a result worth stating. Decide it first with HC-3-volume (below); run
HC-6 only if a target exists.

### REDESIGN

**HC-3 → HC-3-volume (free shooter work, do this first).** The present scan is
one detuning ray at one total charge: it maps a stable *tube*, not the
partition space. Extend D_ab = ∂Q_a/∂ω_b over the full accessible
(Q₀,Q₁,Q₂) volume at several total Q and look for a signature change. Two
outcomes, both useful: it finds HC-6's targets, or it establishes that n(D) = 1
everywhere the flavoured branch exists, which turns the GSS matching from a
local check into a global statement about this sector.

**HC-1 → HC-1-gauged. This is the single biggest gap in the census.** Everything
measured so far is ungauged, which GROUNDING §1 names as the *anchor* — but the
production regime is g = 0.05, Coulomb phase, where GSS is explicitly heuristic.
A gauged BdG needs the A₀ perturbation solved self-consistently with the matter
perturbation (the Gauss constraint couples the channels), which is a real build
but a bounded one: the same tridiagonal machinery plus one Poisson solve per
mode. Until it exists, **every census statement about production objects is
heuristic**, and that is the weakest brick in the Part-0 wall.

**HC-4 → HC-4-lite (one run, not a campaign).** Exactly one bound internal mode
exists in the whole scanned branch (ω ≈ 1.33–1.34, Ω ≈ 0.094 and 0.135). HC-2
predicts its leading radiative order (n = 2–3) from the harmonic table. A single
N=64 run at ω = 1.34, kicked in the monopole channel, either shows the predicted
mode at the predicted frequency with a width consistent with the predicted
order, or it doesn't. That is a genuine falsification test of the HC-1/HC-2
chain and costs one CPU run.

---

## Part 2 — New work that Part 0 forces

### D8 — the discrete momentum-balance theorem  *(highest value, and the reason
### N7's force metrology is only good to 4.4%)*

**The finding.** N7's momentum balance exposed that the measured stress flux
integrates to only 95.6% of the momentum the ball actually gained, and a
plane-to-plane scan shows the force varying **1.4–8.4% between surfaces one
unit apart, in a vacuum gap where a divergence-free momentum flux must be
flat**. The error is in the discrete surface integral, not the physics.

**Why this matters beyond N7.** The kernel has an *exact* discrete Gauss law —
residual pinned at 10⁻¹³, by construction, and it is tripwire #1 for
implementation bugs. There is no equivalent for momentum. Yet a
translation-invariant discretisation should possess an exactly conserved lattice
momentum and an exactly local discrete balance
d/dt P_i(V) = −Σ_∂V (discrete flux), where "discrete flux" is whatever falls out
of the kernel's *own* difference stencil — not the continuum T^{ji} sampled at
voxel centres, which is what I integrated. Getting this right would move N7's
force metrology from 4.4% to the integrator floor and promote F/a from
"excluded, biased" to a first-class independent estimator.

**The rung, in three parts:**

- **D8a — RUN, 2026-07-26. ANSWERED.** One boosted ball (v = 0.05, ω = 1.430),
  *periodic* box, ungauged, no sponge, N = 64, L = 16, T = 60.

  | quantity | drift over the run |
  |---|---|
  | Q | 1.3×10⁻⁹ (exact by construction) |
  | E | 1.1×10⁻⁴ |
  | **P_x** | **3.9×10⁻⁴ (max excursion 4.9×10⁻⁴)** |

  **Total lattice momentum is conserved to 5×10⁻⁴ — ninety times better than the
  4.4% shortfall in the surface flux.** So the defect is *not* lattice
  non-conservation: it is my surface object, and it is fixable. The 5×10⁻⁴
  figure is the ceiling a corrected discrete flux could reach, which is 70×
  below ε — comfortably enough to promote F/a to a first-class estimator and to
  tighten D5′ well below its present 0.5%.

  Note the hierarchy this exposes: Q and the Gauss law are exact *by
  construction* (10⁻⁹, 10⁻¹³); E and P are merely well-conserved by the
  integrator (10⁻⁴). Momentum has no structural guarantee in this kernel the way
  charge does — that asymmetry is itself worth stating in the theory document,
  and closing it is what D8b/D8c are for.
- **D8b (symbolic derivation).** Derive the exact discrete balance from the
  kernel's stencils — covariant link differences, the compact plaquette, the
  leapfrog time-stepping — rather than by discretising the continuum answer.
  This is mechanical symbolic algebra over a finite stencil, and it produces the
  correct discrete flux as its output.
- **D8c (Lean).** State and machine-check the identity. See Part 3.

### D9 — the identity suite, formalised

Part 0 produced several exact statements that are currently backed only by
numerics. They are finite algebraic facts and are genuinely provable:

| identity | status now | provable? |
|---|---|---|
| ωQ = 2(E_kin + E_g) (Gauss) | numeric 5×10⁻⁹ | yes — continuum + discrete |
| Σ = E_∇+E_m+E_V−E_kin−E_g | derived, numeric 5×10⁻⁸ | yes |
| Derrick: (E_∇−E_g)+3(E_m+E_V−E_kin)=0 | derived, numeric 1.4×10⁻⁷ | yes (continuum) |
| Σ = (2/3)(E_∇−E_g) | derived | yes, follows |
| linear mode E = ω\|Q\| + Surf(r_cut) | derived, used in N10 | yes |
| dev = Var(ω)/⟨ω⟩² for a superposition | derived, used in N10 | yes, trivially |
| Ω_c = m − ω (co-rotating continuum edge) | derived | yes |
| n(H_ω) = n(L_x) + n(L₀), and L₀ ≥ 0, L_x^flav ≥ L₀ for μ<0 | derived, numeric | yes — the μ<0 sign argument is a one-liner |
| discrete momentum balance | **open — D8** | the target |

### N8 — ring mass ladder, un-deferred

N8 was parked because D5′ was open: "does M(n) move at fixed Q" is meaningless
without knowing what M is. **D5′ is now closed (M = E/c²)**, so N8 becomes a
clean test of whether winding buys mass at fixed charge — the first
foundation-level crack in the Q-degeneracy. Threshold already pre-registered:
claim a crack only if ΔE/E > 2ε over the accessible n. Reuse the v73 ring
machinery; cheap.

### EX-1 — now unblocked, with one honest caveat about the protocol

D5′ closed means a boosted composite's mass is no longer ambiguous, so EX-1 can
be interpreted. Two arms, and they differ in feasibility:

- **Ramped (the sanctioned protocol).** Amendment 1 requires ≥4 steps with
  spacing ≫ 1/ε₁ = 473 t.u. for X10c. That is ~2000 t.u. of ramp plus a
  measurement window — **~30 h of V100, not the ~20 h the proposal assumed.**
  It is implementable without touching the kernel by restarting from the final
  frame and applying each Δv externally, which needs a new tool (`sfa_boost`,
  in build).
- **Sudden kick at v = 0.02 (also sanctioned).** One seed, one run, and it is
  explicitly retained in Amendment 1 "as a stripping-fraction measurement"
  (pre-registered ~13% strip). Roughly a third of the cost.

Pre-registered from EX-2: the radiative losses should be dominated by above-gap
matter waves, not gauge radiation.

**Seed status: BUILT AND VALIDATED** (`/space/scp/v86/ex1/ex1_seed_v002.sfa`).
A new tool was required — every seed generator in the repo builds a boosted
object from a radial *profile*, and EX-1 needs to boost a relaxed *end-state*.
`n_battery/sfa_boost.c` does it, and three things had to be got right:

1. **The rotation must be exact, not expanded.** The rest-time offset is γvx, so
   the induced phase is γωvx — at v = 0.02 and |x| = 40 that is 1.1 radians. A
   first-order expansion of it lost exactly half the momentum in testing.
2. **The clock must be local.** X10c carries two clocks of opposite sign (ball
   and opposite-charge cloud); a single global ω boosts one correctly and
   corrupts the other. ω_loc is measured per voxel per component, with an
   amplitude floor (the near-vacuum halo's ω_loc is noise, and rotating by a
   noisy ω_loc(x) injects real spurious momentum through x·∂_xω_loc).
3. **The gauge sector is not transformed**, by choice: a correct boost gives the
   links an A₀ piece and breaks the kernel's temporal gauge. E is zeroed and the
   kernel's init Gauss projection rebuilds it.

Validation and calibration:

| check | result |
|---|---|
| static ball vs the purpose-built `gen_qball_boost` | P matches to **0.1%** |
| momentum gain on the X10c state | **1.100 at v=0.01, 1.104 at v=0.02** — linear and stable to 0.4% |
| where the excess sits | the core, not the halo (x>20 carries 0.1% of P) |
| origin | the untransformed gauge sector: in temporal gauge the local phase rate is w_eff(r) = ω − χ(r), not the bare ω |
| calibrated request | v = 0.018764 → **v_eff = P/E = 0.01998**, 0.09% off target |
| charge / energy preserved | −0.13% / −0.05% beyond the deliberate E_em removal |

The 10% gain is a *fixed, linear* property of the configuration, so calibrating
against it is legitimate; it is documented here rather than hidden, and the
run's actual velocity will be measured from the frames regardless. v_eff sits at
**38% of the adiabatic threshold** α_f = 0.053.

---

## Part 3 — Proof architecture: Lean, and what Lean cannot do

The instinct that the momentum-integral error "deserves a dedicated theorem" is
right, and so is the judgement that **Lean alone is insufficient**. Being precise
about why, and about what the missing piece actually is:

### What Lean can carry here, today

Every identity in the D9 table is a *finite* statement: an algebraic relation
among sums over a finite lattice, or a continuum identity provable by parts.
These are squarely in reach of the existing `lean/` setup (Lean 4.29). Proving
them converts a row of **[D]** tags into theorem-grade, and — more valuable —
forces the *hypotheses* into the open. The Derrick identity, for instance,
silently assumes the scaled pair need not satisfy Gauss; the linear-mode
identity assumes a boundary term vanishes. Formalising surfaces those.

### What Lean cannot carry

Lean can prove a theorem about an operator. It cannot certify that
`scp_sim.c` *implements* that operator. That gap is unbridgeable by proof alone
and needs either extraction/codegen (rewriting the kernel — out of scope and
policy-blocked) or a **differential-testing bridge**: the Lean statement and the
numerical residual gate must be emitted from *one* symbolic source, so a drift
between them is a build failure rather than a silent divergence.

### The missing piece, stated concretely

What is wanted — a symbolic engine carrying high-level physical relations
(Maxwell, SR, Newton-with-bounds) as rewrite rules with **validity domains
attached**, so that an obligation discharges as "Newtonian dynamics applies here
because v/c < 0.03 and the bound is tracked through the derivation" — does not
exist off the shelf. Proof assistants have the logic but no physics; CAS systems
have the algebra but no bounds tracking; neither carries the applicability
side-conditions as first-class objects.

The tractable version, which is worth building and is *not* a research project:

```
   kernel stencils  ──►  symbolic derivation (SymPy)  ──┬──►  Lean statement
   (one source of                                       │     (machine-checked)
    truth, extracted                                    └──►  numeric residual
    from scp_sim.c)                                           gate (CI tripwire)
```

with each derivation step annotated by the domain on which it is valid
(v ≪ c, dx → 0, |Φ| below saturation, ω in the branch window). The annotations
are not decoration: they are what turns "the virial identity holds" into "the
virial identity holds on the VK-stable branch of a potential with μ < 0", which
is the actual content.

**Recommended first target: D8.** It is the smallest problem with all three
ingredients — a real numerical defect (4.4%), a derivable exact answer, and a
statement worth machine-checking — and unlike the ε identities it is not already
confirmed to 10⁻⁸, so the pipeline would be *earning* its keep rather than
re-certifying something known.

---

## Recommended order

1. **D8a** — periodic-box momentum conservation test (minutes, CPU). Decides
   whether the defect is my surface object or the lattice.
2. **EX-1** — GPU, authorised. Runs while the rest proceeds.
3. **HC-3-volume** — free shooter; decides whether HC-6 exists.
4. **D8b/D8c** — derive the discrete balance, formalise it, re-run N7's flux.
5. **HC-1-gauged** — the census's biggest gap.
6. **N8** — ring ladder, now interpretable.
7. **HC-4-lite** — one run, falsifies or confirms the HC-1/HC-2 chain.
8. **D9** — the rest of the identity suite in Lean, once the pipeline exists.
