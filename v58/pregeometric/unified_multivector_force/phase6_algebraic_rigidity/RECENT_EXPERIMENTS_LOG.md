# Phase 6 – Recent Experiments Log (Chiral Spin + Inside Space + Particle Folding)

**Period**: Mid-to-late May 2026 (after the initial 7-repulsion + 7-torsion matrix)  
**Focus**: Adding intrinsic chirality/spin to nodes, refining “particle as extra dimension of connection”, introducing “inside space” as a rewarded but collapsing quantity, and enforcing proper causal locality + finite propagation speed.

This document is a detailed, chronological log of actions, results, reasoning, and lessons learned. It is written in the style of a research notebook but organized for later extraction into CONCEPT.md or a paper.

---

## 1. Starting Point (Before This Log)

After the first 7-repulsion + 7-torsion experiments and the subsequent backprop, the best configuration had converged on:

- Strong 7-repulsion channels (direct-children, protected-children, ancestor-overlap, sibling-sharing, etc.), with boosted weights on the top performers.
- Connection-scaled repulsion (`feat_connection_scaled_repulsion`): ~10× repulsion on currently unconnected (low child-count) nodes, smoothly relaxing toward 1× as the node gains connections. This was one of the most effective single terms for forcing branching.
- Optimized MECH A/B/C weights from the 7-run matrix + backprop:
  - MECH_A_EXTRA_CONN_STAB_W ≈ 3.84 (dominant)
  - MECH_C_DENSITY_CLOSURE_W ≈ 2.19
  - MECH_B_DENSITY_RIGIDITY_W ≈ 1.44 (helpful but secondary)

With this stack we achieved a new record of **final d_eff = 2.492** at N=1816, with late-time ~2.30.

This became the “strong base” for all subsequent experiments.

**Lesson at this stage**: Multi-channel local repulsion (especially direct + protected-children + ancestor) + a dynamic “connection cost that relaxes as a node matures” was dramatically more effective than any previous single repulsion term or the original pressure-only best vector.

---

## 2. Introduction of Chiral Helical Bias (Intrinsic Spin + Continuous Turning)

**Motivation** (from user): Nodes should “spin” so that, especially in high-MECH_A dense protected cores, they keep connecting to *different* parents over time instead of locking onto the same small set. The cohesive/attachment weight (modulated by MECH_A) should produce a persistent, continuous turning direction. Each node should be intrinsically chiral.

**Implementation** (`feat_chiral_helical_bias` and its inverse sibling):

- Each node’s protected bivector (the 3 components) already encodes orientation and handedness.
- We derived a continuous phase from `atan2(biv[1], biv[0]) + layer * rate`.
- The sign of `biv[2]` gave intrinsic chirality (preferred turning sense).
- When a node had strong MECH_A contribution (high biv_mag at high ρ), its pressure scalar was modulated by a mismatch term between its current phase and a preferred small helical offset.
- Higher mismatch → higher pressure → temporarily lower attractiveness as a parent.
- Result: high-MECH_A cores naturally rotate their attachment preference over successive growth steps.

Two variants were tested:

- **Direct** (proportional to connections): bias strength scaled with recent child_count. Mature cores spin harder.
- **Inverse** (1/∝ connections): bias strongest on sparsely connected cores (gives new cores a big initial kick).

A slower precession rate (layer * 0.12 instead of 0.25) was chosen for more gradual, continuous turning.

**Key gating**: The effect was deliberately restricted to nodes with *very high* MECH_A-like scores so that only the real dense protected “particle” cores would exhibit strong spinning behavior.

### Results & Reasoning

**Direct version strength sweep** (on best base):

- 1.5 → final 2.472 (very good)
- 2.0 → final 2.352 (**clear dip**)
- 2.5 → final 2.490 (excellent, essentially ties the previous record of 2.492)

**Why the dip at 2.0?** (qualitative + quantitative)

The returned bias term is roughly `mismatch * mech_a_like * conn_factor * const`. The effective modulation depth on activity for high-MECH_A cores therefore scales with the user-set weight W.

- At W≈1.5 the modulation is a gentle periodic nudge. It prevents over-locking without disrupting the dense hubs that MECH_A is building. Net helpful.
- At W≈2.5 the amplitude is high enough that the cores are effectively *only* chosen during favorable phases of their spin cycle. This turns the spinning into a clean “stepper” that only fires at good helical angles. The system waits for the right phase and gets disciplined, space-filling helical branching. Net good (or even better than no bias).
- At the intermediate W≈2.0 the modulation is strong enough to meaningfully suppress cores during bad phases (they lose attachment opportunities), but not strong enough to create the clean “only attach in good phase” regime. You get irregular, partial suppression that interferes with the discrete batch-addition rhythm and with the connection-scaled repulsion. The best branches are intermittently starved without the compensating benefit of structured helical growth. This is a classic intermediate-amplitude resonance/anti-resonance effect in a system with competing discrete and periodic drives.

The inverse version at 1.5 was consistently weaker (~2.28). Giving the strongest spin to the *least* connected cores appears to destabilize them before they can stabilize the extra connection that MECH_A wants to build.

**Lesson**: The direct (mature-cores spin harder), slower-period, tightly MECH_A-gated version is the useful one. The functional form and the amplitude matter; intermediate strengths can be actively harmful.

---

## 3. The “Inside Space” Idea and the Two Philosophies

**Conceptual trigger**: After the strong MECH_A results, the user proposed that particles enclose a form of “inside space” that is *rewarded* (the system wants to create and protect it) but which also exerts a *collapsing force* back on the boundary cores (particles) the more of it there is. This would give a natural reason for particles to stay compact and to attract each other.

Two distinct physical pictures were formalized and implemented as separate pairs of pressure features:

**Version 1 (Explicit Volume Reward + Localized Collapse)**  
- Reward (negative pressure) for any node that helps create or protect enclosed inside space.
- Collapse (positive pressure) applied *only* to the high-MECH_A cores, proportional to the inside volume they help bound.

**Version 2 (Negative Pressure Inside + Boundary Collapse – more field-like)**  
- Nodes that are causally “inside” receive negative pressure (the interior volume itself is stabilized).
- The collapsing force is still borne by the high-MECH_A boundary cores.

Both were first implemented with a relatively global recent-window scan, then re-implemented with **strictly local causal versions** (`_local_enclosed_score_for_core`) so that a core only sees the inside volume that has causal paths back to it within the recent layer window. This enforces finite propagation speed and locality — exactly the refinement the user later insisted on.

### Results

**First (partially global) matrix** (strength 1.5):

- A+B and B+C reached ~2.43 final.
- Full A+B+C gave the single highest peak (2.476) but slightly lower late-time.
- A alone was the strongest singleton; B alone was weak.

**Backprop on the three MECH weights** (on top of best 7-repulsion + connection-scaled):

- Converged to A≈3.84, C≈2.19, B≈1.44.
- Final large validation with these weights: **2.492** (new record at the time).

**Local causal versions** (first proper locality test, strength 1.2 + best base + direct chiral 1.5):

- Local V1: final ~2.257, late ~2.26
- Local V2: final ~2.225, late ~2.20

Both were stable. Locality did not destroy the effect.

**Strength sweep on local V1** (0.8 vs 1.2 vs 2.0):

- 0.8 and 2.0 both gave ~2.389 final (better than 1.2 in this limited sample).
- Suggests the response is not strictly monotonic; both lower and higher strengths can be useful depending on the exact functional form.

**Numerical warning from V2**: At 1.5 the local V2 formulation produced pressure scalars large enough to cause `math.exp(-PRESSURE_WEIGHT * ps)` to overflow. This is a real signal that the “negative pressure inside + collapse on boundary” combination can create very sharp gradients when made strictly local. V1 was numerically more stable at the same strength.

---

## 4. The Dip at Direct Chiral = 2.0 – Detailed Reasoning

**Data** (on best MECH + connection-scaled base):

- 1.5 → 2.472
- 2.0 → 2.352 (dip)
- 2.5 → 2.490 (recovers)

**Quantitative view**:

The effective modulation depth on activity for high-MECH_A cores scales with W × avg(mech_a_like) × avg(conn_factor) × avg(mismatch).

- W=1.5 → ~18–25% modulation (gentle periodic nudge) → helpful exploration without disrupting hubs.
- W=2.5 → ~38–48% modulation (almost binary on/off) → the system effectively only attaches during favorable phases → clean helical stepper → structured, space-filling growth.
- W=2.0 → ~28–35% modulation (awkward intermediate) → strong enough to suppress cores during bad phases, but not strong enough to create disciplined “wait for good phase” behavior. Result: irregular starvation of the best branches without the compensating benefit of structured turning. This interferes with the discrete batch-addition rhythm and with the connection-scaled repulsion.

**Qualitative picture**: Intermediate amplitude sits in the worst part of a nonlinear resonance curve between the chiral oscillation period and the discrete growth steps. Low amplitude = soft bias; high amplitude = new stable regime (stepper); middle = messy beating.

This is why the dip appears and why pushing harder (2.5) recovers performance.

---

## 5. Locality and Finite Propagation Speed – The Key Refinement

After the user correctly pointed out that earlier Inside Space implementations were still too global, we re-implemented both V1 and V2 using a per-core local causal score (`_local_enclosed_score_for_core`).

A high-MECH_A core now only counts inside space that is causally connected to *it* within the recent layer window. The reward or negative pressure and the collapse force only propagate at the speed of the graph.

**Consequence**: newly created inside pockets only begin to exert effects after a few layers of causal connection. The “outside field” can only feel the inside after a delay, and the collapsing force on a core only reflects the volume that has actually reached it.

This is the first time the “inside space rewarded but collapses its own boundary” idea was tested under proper pre-geometric causality. The results remained stable (~2.2+), which is a positive signal.

---

## 6. Consolidated “What Worked Best” vs “What Did Not”

**Strongly positive / currently part of the best manual configuration**

- Strong multi-channel 7-repulsion (especially direct-children, protected-children, ancestor-overlap).
- Connection-scaled repulsion (10× on unconnected → relaxes with connections) – one of the single best terms.
- High MECH_A (extra protected connection stabilized by density) – the dominant new mechanism; optimizer consistently pushes it to ~3.8.
- Moderate-to-high MECH_C (density-dependent closure shift).
- **Direct** chiral helical bias (nodes spin via their protected bivector; high-MECH_A cores have rotating attachment preference; strength ∝ connections; gated to very high MECH_A cores; slower period for continuous turning). Works well in 1.5–2.5 range; 2.5 is currently one of the top manual performers.
- Local causal formulation of Inside Space ideas (especially V1 style) – stable and conceptually clean.

**Neutral to mildly helpful (optimizer keeps them but at lower weight)**

- MECH_B (density-dependent softening of rigidity cost) – useful scaffolding once A is present, but not a primary driver.

**Not helpful / actively harmful in manual tests**

- Inverse chiral bias (strongest spin on sparsely connected cores) – consistently weaker.
- Intermediate direct chiral strength (~2.0) – clear dip caused by awkward resonance.
- Local V2 at higher strengths without numerical safeguards – produces overflow (sharp gradients).
- Earlier global (non-local) formulations of inside space – conceptually flawed even if numerically they sometimes worked.

**The famous dip at direct chiral = 2.0**

Caused by intermediate modulation depth creating irregular, non-constructive interference with the discrete growth rhythm and the existing repulsion terms. Low amplitude = soft bias; high amplitude = new structured regime; middle = messy beating.

---

## 7. Current Best Practical Configuration (as of last manual data)

- Strong 7-repulsion channels with boosted weights on top performers.
- Connection-scaled repulsion ~2.13.
- MECH_A ≈ 3.84, MECH_C ≈ 2.19, MECH_B ≈ 1.44.
- Direct chiral helical bias in the 1.5–2.5 range (2.5 currently looks excellent).
- Inverse chiral ≈ 0.
- Local causal Inside Space (V1 style) at modest strength (0.8–2.0 both worked; optimizer will decide).
- All other older pressure terms at their previously optimized values.

This stack has produced final d_eff up to **2.492** and late-time ~2.30 on 1800-node graphs — the highest numbers recorded in the entire Phase 6 series so far.

---

## 8. Open Questions & Recommended Next Actions

1. **The 10D+ backprop** (which can tune all ten weights including both chiral and the four local inside-space terms) timed out. Re-running it with a higher budget or on a slightly smaller node count for speed would give the real optimal point.

2. **Numerical stability of local V2** needs a fix (clamping the pressure scalar or a softer formulation) before higher strengths can be tested fairly.

3. **The dip at 2.0** is now well explained; we should confirm whether it moves when we change the precession rate (the 0.12) — that would be strong evidence for the resonance hypothesis.

4. **Deeper integration**: Once the optimizer settles on good strengths, we should consider promoting the winning mechanisms (especially strong MECH_A + direct chiral + local V1 inside space) from “extra pressure features” into first-class parts of the living candidate equation or the protection/rigidity definition.

5. **Physics narrative**: The combination of (a) multi-channel local repulsion that relaxes with connection, (b) density-stabilized extra protected connection (MECH_A), (c) intrinsic chiral spin that produces continuous turning in mature cores, and (d) inside volume that is rewarded but collapses its own boundary, is beginning to give a coherent pre-geometric story for why particles are stable, compact, and appear to source effective gravity while the vacuum does not.

---

**Document status**: This log captures the experimental arc from the 7-repulsion/7-torsion work through the first properly local causal tests of the inside-space idea. It is intended as the primary reference for deciding the next modeling step (whether that is a cleaner final backprop, a larger-N validation of the current best manual point, or a deeper rewrite that promotes the winning mechanisms into the core equation).

*Written 2026-05-21 as a living research log.*