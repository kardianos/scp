# v65 FINDINGS — Mechanism (A): intrinsic kernel self-tuning (implemented; substrate-blocked)

**Date**: 2026-06-05
**Goal**: realize mechanism (A) of [`SELF_TUNING.md`](SELF_TUNING.md) *inside* the kernel —
κ a dynamical variable that self-organizes to the stability edge at fixed conserved charge,
as a term in the time evolution (not an outer loop). **Kernel change authorized by user.**

---

## What was built (and works)

Added an opt-in `self_tune` mode to `scp_sim.cu` (+ config fields in `scp_config.h`),
gated by `self_tune=1` so the default physics is byte-for-byte unchanged when off:
- κ promoted to a dynamical variable, updated every `st_dt` time units in the main loop;
- **SOC law**: `κ *= (1−st_eps)` when the core is stable, `κ *= (1+st_gamma)` on collapse
  (collapse sensed via `P_max > st_pcrit` or NaN) — the action pull down + stability push up;
- **charge projection** (`st_project=1`): rescale φ each tune step to hold `Q=∫|φ|²=2·E_mass/m²`
  fixed (the mechanism-A Lagrange constraint), reusing the existing diagnostic reduction;
- a `rescale_phi_kernel`; the `[self-tune]` trajectory streams to stdout.

**Verified functional**: compiles (nvcc), κ is dynamical, the SOC update and charge
projection run, the log streams. The *implementation* is correct and ready.

## What it revealed (the blocker) — no persistent soliton to self-tune against

Across **six** runs (oscillon & neutron seeds; bc absorbing & periodic; A_bg 0 & 0.1; κ₀ 4–50;
slow & fast descent; project 0 & 1), the same wall:

- **The field always decays/disperses.** The oscillon disperses; the neutron template *melts*
  even at its own equilibration κ=50 with energy-conserving (periodic) BC (`P_max 0.10→10⁻⁹`,
  `Q 7.8→0.003`). A_bg=0.1 slows but does not stop it.
- Because there is **no persistent dense core**, the collapse the feedback needs never
  develops: `P_max` decreases monotonically as the field melts, the sensor never fires, and
  κ drifts to the floor (≈0) with no edge.
- **Charge projection (`st_project=1`) over-stabilizes**: holding `∫|φ|²` fixed caps the
  amplitude and *prevents* collapse entirely — the wrong charge (it suppresses the very
  instability that defines the edge). Q bounced 1.5–50 re-inflating a dispersing field.
- **Fast descent got closest**: P_max rose to **1.89 at κ≈1.9** (right at the X1-measured
  κ_crit≈1.5) before the field dispersed and P_max fell back — the κ-response is real but
  the collapse cannot complete on a non-persistent soliton.

## Root cause (a known project fact)

These solitons are **metastable lattice saddles, not stable attractors** (MEMORY: *"Topology
loss: gradient flow loses Q at all tested N; Skyrmion is a lattice saddle point"*). So there
is no persistent dense object for a *continuous* self-tuning to push against. Two coupled
problems:
1. **No persistent soliton** — every seed melts/disperses on the descent timescale.
2. **Collapse is irreversible within a run** — once a core spikes it cannot un-spike, so even
   if triggered, the SOC feedback cannot *oscillate* around the edge the way X1 does.

**Why X1 (mechanism B) worked and this does not**: X1 re-seeds a *fresh compact oscillon
every step*, so each evaluation has a dense core that genuinely collapses at low κ
(P_max→14.7), and each step is independent/reversible. The continuous intrinsic version loses
both properties.

## Status and path forward

- **Mechanism validated**: by X1 (`FINDINGS_X1.md`) — SOC self-tuning to κ_crit≈1.46 is real.
- **Intrinsic kernel mechanism**: *implemented and ready*, but blocked on the soliton-
  stability substrate. It is NOT a refutation of (A) — it is a prerequisite that is unmet.
- **To make the intrinsic version work, one of**:
  1. **A genuinely stable soliton** (solve the lattice-saddle problem first) — then adiabatic
     κ-descent on a persistent object would reach κ_crit cleanly. This is a real open problem.
  2. **A reversible stability-margin sensor** instead of irreversible collapse — e.g. tune κ
     on the soliton's *lowest fluctuation eigenvalue* (the normal-mode / NTK operator the
     project already computes): self-tune κ so the marginal mode sits at ω²=0 (the edge),
     which is reversible and needs no collapse. This is the cleaner intrinsic formulation and
     reuses existing machinery (`normal_modes`, BLV).
  3. **A re-compacting term** (periodic re-projection of the field onto a compact ansatz) —
     essentially X1's re-seed made intrinsic; least principled.

**Recommendation**: option (2) — drive κ by the **fluctuation-spectrum margin** (ω²_min→0),
not by collapse. It is reversible, action-grounded (the Hessian *is* δ²E/δΦ²), and sidesteps
the persistence problem. Requires adding a (cheap) lowest-eigenvalue estimate to the tune step.

**Artifacts**: kernel diff in `scp_sim.cu`/`scp_config.h` (opt-in `self_tune`); six run logs.

---

## Addendum — bidirectional density homeostasis (user's idea): implemented, same wall

The unidirectional flaw above (dissolution misread as "stable") motivated a **bidirectional
density feedback** (added: `st_ptarget`, `st_gain`): regulate the core density `P_max` to a
setpoint by `κ *= (1 + gain·(P_max−target)/target)` — higher κ lowers equilibrium density, so
this is stable negative feedback that should correct *both* dissolution (P_max<target → κ
down → re-concentrate) and collapse (P_max>target → κ up → spread).

**It works in the collapse direction.** From a compact oscillon (P_max=1.1), the feedback
correctly raised κ and drove P_max to **exactly the target (0.33 at t=6)**.

**It cannot prevent dissolution.** Once P_max fell below target, lowering κ (all the way to
10⁻³) did **not** re-concentrate the field — P_max kept crashing to ~0. Two BCs bracket the
failure:
- **periodic** (energy-conserved): the dispersed field sloshes forever and is never
  re-gathered (`Q~2` stays in the box but spread out, `P_max→2×10⁻⁵`);
- **broad damping** (global cooling, `damp_width=9`): the field is *absorbed away* entirely
  (`Q 17→10⁻⁹`, `P_max→10⁻¹⁹`).

**Root cause (two compounding facts):**
1. **κ has asymmetric authority** — it spreads a dense core (raise κ) but cannot re-gather a
   dispersed one (lower κ): the restoring force `∝V'(P)=μP/(1+κP²)²` is negligible once P is
   small. Control is effectively one-directional in *authority* even when bidirectional in
   *intent*.
2. **Re-gathering requires removing the dispersal (kinetic) energy, but removing it removes
   the field** — there is no operating point where dissipation re-concentrates instead of
   deleting.

Underneath both: **no stable soliton attractor exists** for the feedback to hold (the
lattice-saddle fact). A controller can hold a system *at* a stable fixed point; it cannot
*manufacture* one. The user's bidirectional concept is correct; it is defeated by the
actuator (κ) authority and the missing attractor, not by the control logic.

**Conclusion (updated)**: self-tuning κ — by collapse (SOC) *or* by density homeostasis —
cannot work on a substrate with no persistent soliton. The remaining viable routes are
unchanged: **(2) the reversible eigenvalue-margin sensor** (tune κ so ω²_min→0; reads the
stability margin of whatever quasi-static soliton is present, needs neither collapse nor
re-gathering), or first **solve the soliton-stability problem** (a genuine stable attractor).
The density-feedback fields (`st_ptarget`, `st_gain`) remain in the opt-in `self_tune` mode.
