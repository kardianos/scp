# V48 Proposal — Structure-First Binding: Time-Averaged Topology-Dependent Coupling

## The Problem (Restated)

V45 showed no nuclear binding with constant η=0.5 — the energy is minimized
at infinite separation. V47 showed that amplifying η with η(P) causes
structural COLLAPSE (merger) rather than binding. Neither approach produces
what real nuclei are: two distinct particles bound at finite separation
without merging.

The force-balance / parameter-tuning approach has failed. The problem is
structural, not parametric.

## Insight 1: Structure Creates Binding, Not Forces

Previous approaches treated binding as a force-balance problem: find the
right coupling strength where attraction overcomes Derrick's gradient
penalty. This is the wrong framing.

Binding is not about balancing forces. It is about CONVERGENT STRUCTURE —
a field configuration that maintains itself as a dynamical attractor.

Evidence:
- A single braid is a structural attractor. It maintains itself at T=500+.
- A 3-braid baryon with phase offsets {0, 2π/3, 4π/3} is a structural
  attractor. The phase cancellation P→0 at overlap CREATES the confinement.
  No force holds it together — the geometry prevents it from dissolving.
- A two-baryon bound state must be a structural attractor of the SAME kind.
  The field naturally converges to it, maintains it, and resists perturbation.

If such a structural attractor doesn't exist for the current equation, no
amount of parameter tuning will create binding. If it does exist, we need
to find it — not engineer it.

## Insight 2: It's Not 3D Geometry — It's 4D Structure

The braid is not a static 3D shape. It's a 4D pattern: a helix in space
TIMES an oscillation in time. The particle IS the oscillation. A snapshot
showing a "breathing sphere" is a cross-section of a process, not the
thing itself.

Two baryons binding is not two 3D objects touching. It is two 4D
oscillation patterns interlocking. Breathing synchronization isn't a
side effect of binding — it IS the binding. The spatial structure follows
from the temporal coherence.

This means:
- T=500 (≈3 breathing cycles) is far too short to see binding form
- The binding structure involves PHASE relationships in time, not just
  proximity in space
- The carrier phase, breathing phase, and θ radiation phase are all
  structural coordinates — not just the spatial position

## Insight 3: The V47 Blowup Was a Rate Problem

V47 used η(P) = η₀ + η₁P²/(1+κP²) with instantaneous P. The coupling
changed at every timestep, tracking the carrier oscillation (~4t period).
This created a feedback loop running at the CARRIER frequency — orders of
magnitude faster than the STRUCTURAL timescale (~150t breathing period).

The physics of the blowup:
- t=0: P oscillates normally at carrier frequency
- η tracks P → oscillates at carrier frequency
- Oscillating η amplifies the carrier oscillation
- Amplified carrier → larger P peaks → larger η peaks
- Runaway in ~10 carrier cycles (40t)

The structural coupling should respond to the TIME-AVERAGED P, not the
instantaneous P. The time-averaged ⟨P²⟩_τ represents the persistent
binding density — the STRUCTURAL content at that point. The instantaneous
P includes the carrier oscillation, which is not structure — it's the
medium vibrating.

The V46 analytical value η₁≈115 may be CORRECT for the equilibrium. The
blowup happened because we applied the equilibrium coupling at the
wrong timescale.

## Insight 4: Parameters May Be Field-Dependent (Recursive)

If η = η(⟨P²⟩_τ), why not also:
- m² = m²(⟨ρ⟩_τ) — effective mass depends on local energy density
- μ = μ(⟨P²⟩_τ) — potential strength depends on binding density
- κ = κ(⟨P²⟩_τ) — saturation scale depends on local structure

Each parameter responds to the local field state. The field state depends
on the parameters. The physical solution is the FIXED POINT where everything
is self-consistent.

This is not parameter tuning — it's finding the attractor of a
self-referential system. The physics IS the fixed point. Constant
parameters are the zeroth-order approximation.

Implement ONE parameter at a time. η first (directly addresses binding).
Others later, after η stabilizes.

## Insight 5: What Structure Would Work?

The current braid is 1D (helical along one axis). The 3-braid composite
is assembled from 1D pieces. These might be metastable structures — they
survive for simulation timescales but may not be the TRUE attractors.

A proton-like structure might need:
- Genuinely 3D topology (not just composited 1D pieces)
- Chirality (handedness — distinguishes particle from antiparticle)
- Use of ALL 6 DOF (φ₀,φ₁,φ₂,θ₀,θ₁,θ₂), not just the φ sector
- The θ sector participating in the STRUCTURE, not just radiating from it
- A topology that supports LINKING (two structures interlocked without
  merging) rather than just PROXIMITY (two structures placed nearby)

Candidates from topology:
- Hedgehog (Skyrmion-like): φ wraps S² with winding number
- Hopfion: linked field lines (Hopf fibration)
- Vortex ring: braid closed into a circle
- Linked rings: two vortex rings passing through each other (one φ, one θ)

The linking number of two rings is a topological invariant — it provides
binding that CANNOT be undone by continuous deformation. This is binding
through geometry, not through forces.

BUT: these are still 3D geometric objects. The full picture requires
4D structure — the time evolution is part of the topology. A "particle"
might be a vortex ring that PRECESSES, with the precession frequency
being an additional structural parameter.

This is a hard open problem. The topological analysis can be done with
algebra (homotopy groups, linking numbers), but constructing explicit
solutions requires either analytical insight or numerical search.

## Insight 6: The Binding Limit

There may be a fundamental limit to binding strength per timestep. Just
as c limits information propagation, there may be a maximum rate at which
the coupling can change. The V47 blowup occurred because we violated this
limit — we allowed infinite-rate coupling change.

The time-averaged coupling η(⟨P²⟩_τ) naturally enforces this limit.
The parameter τ IS the binding rate limit. Physics cannot bind faster
than one breathing cycle.

## Proposal: Time-Averaged Topology-Dependent Coupling

### The equation modification

Add one new field: P_avg(x,t) — the exponentially smoothed time-average
of P²(x,t) at each point.

Update rule (each timestep):
```
P_avg[idx] = (1 - α) × P_avg[idx] + α × P(x,t)²
```
where α = dt/τ, τ ≈ breathing period ≈ 150 code time units.

Effective coupling:
```
η_eff(x,t) = η₀ + η₁ × P_avg[idx] / (1 + κ × P_avg[idx])
```

This uses the SAME saturating form as V47 but with the TIME-AVERAGED P
instead of instantaneous P. The time averaging:
- Rate-limits the coupling change (maximum rate = 1/τ)
- Filters out the carrier oscillation (period ~4t << τ=150)
- Responds to persistent binding density (structural signal)
- Prevents the instantaneous feedback that caused the V47 blowup

### Implementation

1. Add `P_avg[N³]` array to the Grid struct (one extra array, ~1 GB at N=512)
2. Initialize P_avg = 0 (builds up from zero as structure forms)
3. Each timestep: update P_avg with exponential smoothing
4. Compute η_eff from P_avg instead of instantaneous P
5. Config parameters: `eta1` (coupling amplification), `tau` (smoothing time)

Cost: one multiply-add per voxel per timestep + one array. Negligible
compared to the force computation.

### Test plan

1. **Stability test**: η₁=375, τ=150 at N=64 T=100 (5 min CPU)
   - If stable: the rate limit works
   - If unstable: try τ=500 (slower response)

2. **Single baryon structure**: η₁=375, τ=150 at N=192 T=500
   - Does the baryon maintain its structure?
   - Does P_avg build up to the expected equilibrium value?
   - What is η_eff at the core after equilibration?

3. **Binding test**: two baryons at D=10,15,25,80 with η₁=375, τ=150
   - Does binding appear at finite D?
   - Does the structure maintain two distinct objects (not merge)?
   - Does P_avg form a "bridge" between the two baryons?

4. **Phase exploration**: same baryons but with different carrier phase
   offsets (Δ = 0, π/6, π/3, π/2, π)
   - Which relative phase creates the most stable bound pair?
   - Does the system naturally evolve toward a preferred Δ?

### What to look for

The signature of structural binding would be:
- P_avg forms a SADDLE between two baryons (high at each core, with a
  nonzero bridge between them)
- η_eff is elevated in the bridge region (coupling stronger BETWEEN the
  baryons, not just AT them)
- The θ field forms a standing wave in the bridge region
- The breathing phases of the two baryons converge over time (synchronization)
- E_total(D_bound) < E_total(D=80) for some finite D_bound

The signature of structural failure would be:
- P_avg collapses to a single peak (merger)
- P_avg stays as two isolated peaks with no bridge (no interaction)
- The system blows up (still unstable even with time averaging)

### Future: recursive parameters

After η(⟨P²⟩_τ) is validated, the same mechanism can be applied to
other parameters:
- κ_eff(⟨P²⟩_τ) — saturation responds to local binding density
- μ_eff(⟨ρ⟩_τ) — potential depth responds to local energy density

Each adds one smoothed field and one parameter. The system of smoothed
fields evolves toward a self-consistent fixed point — the parameters
are not tuned but DERIVED from the field's own structure.

### Connection to the topological question

The time-averaged P_avg is a DETECTOR of topological content. Regions
with persistent nonzero ⟨P²⟩ have braids or binding structures.
Regions with ⟨P²⟩≈0 are background. The coupling η(⟨P²⟩) naturally
partitions space into "structured" (strong coupling) and "unstructured"
(weak coupling) regions.

The topology question (what 3D/4D structures exist?) is separate from
the coupling question (how strongly do they interact?). But the time-
averaged coupling ENABLES the topology by providing the energy to
maintain complex structures that constant coupling cannot sustain.

The topological exploration (hedgehogs, hopfions, linked rings) should
proceed in parallel as an analytical/seed-generation effort, while the
time-averaged coupling provides the numerical framework to test whether
these structures are dynamically stable.

## Summary

| Idea | Mechanism | Status |
|------|-----------|--------|
| Constant η (V45) | No amplification | No binding (Derrick) |
| Instantaneous η(P) (V47) | Fast feedback | Unstable (blowup/merger) |
| **Time-averaged η(⟨P²⟩_τ)** | **Rate-limited structural coupling** | **Proposed (V48)** |
| Recursive parameters | Self-consistent fixed point | Future |
| Topological seeds | Geometric interlocking | Future (analytical) |

The core hypothesis: binding fails not because η₁ is wrong, but because
the coupling responded at the wrong timescale. The structural coupling
must evolve at the STRUCTURAL timescale (breathing period), not at the
carrier oscillation timescale (4t). Time-averaging enforces this
naturally and may be all that's needed.
