## Timeline reading (5–10 bullets)

- **Stages 0–2A are genuinely done and reusable.** Z-carbon parks, the fusion branch, ledger discipline (gauss_max ~1e-13 as tripwire culture) are assets. Nothing in the Stage-3 failure invalidates the nuclear stack.
- **The Stage-3 failure is intrinsic to the representation, not a tuning miss.** In the Q-ball framework, "light" means small Q, which means sitting near the bottom of the stability window — a light multiplet hump is *definitionally* the most fragile object in the theory. Multi-fab L inherited the same Cosserat potential and therefore the same fragility floor. An electron analog built this way will always shred (hard v_t) or evaporate (long T), because marginal solitons do that.
- **v79's E_em nulling is diagnostic, not incidental:** a wide light multiplet acts as an *absorbing bath* for the shared gauge field. Any hump-electron screens the nucleus rather than orbiting it. This is the multiplet-shadow of "L never opened a sequestration channel" (DUAL_CHARTS.md).
- **Force physics is real; product physics is not.** F PASS (monotone a_rel(D)), E-lite charge bookkeeping, and S_pair E_em survival (v80 G=0.62) show the *free gauge medium does the Coulomb job correctly*. What fails is only the light *matter carrier*. The redesign should keep the medium and replace the carrier.
- **v80's S4 "hold" being a COM/window artifact is a second lesson: diagnostics on Eulerian windows (Q_core, Q_flux) cannot certify bound products.** Any new design must make binding a property of tracked objects, not of a fixed box.
- **The two-stacks honesty (v78) is the crack the redesign walks through:** the monist free/bound language could *score* everything but never *was* the state. v80's kill gates are correctly aimed at preventing a third round of that.
- **v76's GRIN kill still binds:** whatever we build, c stays fixed in the physical (C-chart) law; sequestration thickens the measure, never the local speed.
- **Timeline pressure favors a standalone sandbox.** Kernel modification is authorized, but the last three versions spent GPU-nights on runs whose *state type* was the problem. Two CPU weeks on the right state type is cheaper than one more V100 night on the wrong one.

## Three redesign options (ranked)

### Option A (primary) — Lock-carrier engine: discrete typed locks on a real free medium

The direct implementation of the 4+3+6+10 blend: keep the gauge medium (which works), replace the light matter carrier with a first-class lock object.

**State.**
- *Free medium (fields, as view-compatible substrate — allowed per CONTEXT.md §"not claiming"):* Yee-lattice gauge sector (E on links, B on plaquettes), plus a scalar **capacity field** n(x) with per-site split n_free + n_seq = n_tot (shape-3 content: path measure ℓ is a readout of n_seq thickening).
- *Locks (array of structs, not fields):* `{id, type ∈ {N,e}, q = ±1, E_i (energy ledger), p_i, X_i (continuous), a_i (footprint radius), S_i (sequestered capacity), clock φ_i, ω_i}`. Rest energy is **derived**: E_i(rest) = μ·S_i. Mass is never stored — M_i ≡ E_i/c² is an M-chart readout.

**Step (C-chart law, c = 1 lattice unit/tick, fixed).**
1. FDTD Maxwell update; lock currents deposited with a charge-conserving (Esirkepov-class) smeared scheme → discrete Gauss law at machine floor by construction.
2. Lock push: relativistic, v = pc²/E_i (so |v| < c is a ledger consequence, not a clamp — shape 6 partially literal), force = q(E + v×B) interpolated at the footprint.
3. Sequestration force: each lock maintains footprint w(x−X_i) claiming n_seq; overlap beyond n_tot costs energy → medium-mediated contact/saturation core. **No pairwise potential anywhere in the code.**
4. Ledger: dE_i = F·v dt exchanged exactly against field energy via the discrete Poynting theorem; boundary flux accounted.

**How it maps v80 shapes.** 4 = architecture verbatim (locks are objects; particle ≠ hump definitionally). 3 = capacity field is the medium content; ℓ readout from n_seq. 6 = field update is exactly light-cone limited and lock speed is budget-limited by its own ledger. 10 = charge lives on the gauge links, sourced by locks, Gauss-exact. Charts: C = the step law itself; M-chart = post-processing view holding M_i(0) fixed and reporting residual as per-lock energy-table/c_eff drift; E-chart = fixed region budgets with Poynting honesty. All three are re-partitions of one state → passes representation kill gate 3.

**Kernel delta.** *None for two weeks.* Standalone `v81/sandbox/lock_medium.c` (CPU, OpenMP). SFA only as `export_grid(state)` for volview; lock tracks as TSV. Kernel promotion (a `locks` block in scp_sim, possible sfa.h particle-track chunk) is a *later, separate* authorization spend, only after the sandbox passes.

**First experiment.** 96³–128³, CPU. (i) Heavy lock (S=100) + light lock (S=1, opposite q) at rest, D ∈ {8,12,16,20,24} → a_rel(D) must be Coulomb-class *with no potential term in code*. (ii) Tangential kick → radiative inspiral → settle at the saturation-core separation; (iii) positronium proper: two S=1 locks, ±q, same protocol. T ~ 3000–5000 ticks. Observables: X_i(t), E_i(t), field energy, boundary Poynting, Gauss residual, E_bind = E_final − 2E_free.

**Success metrics.** Coulomb-class force from medium only; durable bound final state (≥20 orbital periods, or a static minimum with two-sided restoring force); total ledger closure ≤1e-6 relative including boundary flux; Gauss at floor; M- and E-chart readouts implemented and mutually consistent.

**Kill gates.** Dead if: attraction requires adding an explicit 1/r² term; the bound state exists only with radiation switched off; ledger closure needs a fudge sink; the capacity field turns out dynamically inert (pure bookkeeping, no force) — then the free medium is fake and this is PIC.

**Honest risk.** Medium-on-grid is acceptable (grid as view), but the lock could degenerate into a PIC macro-particle with a dressed-up contact potential. The three commitments that keep it honest: rest energy *from* sequestered capacity, mass always derived, core force always through the shared capacity field. If any of the three has to be faked, say so and kill.

### Option B (charge-honesty upgrade) — Integer-flux substrate: matter as flux-string endpoints

Shape 10 taken as primary. **State:** integer (or quantized) electric flux on links + a derived (not stored) endpoint list; charge = divergence of flux, an integer per site — topological, impossible to hand-place as a hump. Optional per-site budget b(x) for shape 6. **Step:** div-preserving plaquette moves plus endpoint hops, deterministic steepest-descent with an inertia ledger; c = one hop per tick, *exactly* literal. **Binding for free:** opposite endpoints are connected by flux; strong coupling → linear confinement (always bound), weak coupling → spreading flux → Coulomb-class; the coupling dial interpolates. Positronium = a rotating flux dipole. **Kernel delta:** none — pure sandbox, and it should stay one for a while. **First experiment:** 64³, dipole at D ∈ {8,16}, force vs D across coupling; search for a regime with both attraction and non-instant annihilation of adjacent endpoints. **Kill gates:** no coupling regime gives attraction *and* persistence; dynamics rules need per-experiment tuning; inertia/mass cannot be defined from the state. **Honest risk:** real-time dynamics for discrete flux has no canonical form — this is invented dynamics and may behave like a cellular-automaton toy. That is why it is ranked second: highest structural honesty, highest research risk, worst week-scale fit. Its role: if Option A's Maxwell sector feels like "imported EM," B is the replacement charge sector, slotting under A's lock layer.

### Option C (fallback / Stage-5 seed) — Measure-primary with free↔bound exchange

Shape 3 taken as primary. **State per site:** (n_free, n_bound, gauge links); locks are *emergent* — connected components of n_bound, with a local condensation/evaporation exchange rule (direct descendant of the v73 uptake/layment ledger). **Step:** hyperbolic transport of n_free at fixed c; exchange driven by local intensity; bound blobs advect. **Unique win:** locks can form and dissolve dynamically — this is the only option with a native path to Stage 5 (spontaneous production). **Kernel delta:** none; sandbox. **First experiment:** seed two opposite bound blobs, watch whether sequestration persists and binds. **Kill gate:** this is the option most exposed to v80 kill gate 1 — n_bound(x) *is* a field on N³, and the design degenerates into two-fluid plasma. Use it only if Option A succeeds and its capacity field asks to be promoted to carry the locks itself. Do not start here.

## Recommended primary path + first 2 weeks of work

**Primary: Option A**, sandbox-first, with B held as the later charge-sector upgrade and C as the Stage-5 follow-on. Rationale: the Coulomb medium is the one thing v75–v80 proved works; the light carrier is the one thing that provably cannot be a multiplet hump; Option A changes exactly the broken part, is implementable on CPU in days, and satisfies all five v80 kill gates by construction rather than by argument.

- **Days 1–2:** `v81/` scaffolding; FDTD medium + charge-conserving deposition. Unit gate: one static lock → Gauss at floor, correct exterior Coulomb profile.
- **Days 3–4:** two static opposite locks; measure force via ledger momentum transfer; a_rel(D) for D=8..24. This is the *F-repro gate*: reproduce multi-fab's F PASS with zero multiplet fields. If this fails, stop.
- **Days 5–6:** capacity field + sequestration: rest energy from S_i, saturation core; verify a static equilibrium separation exists.
- **Days 7–8:** relativistic push + radiation; inspiral run with full boundary Poynting accounting; ledger-closure gate.
- **Days 9–10:** positronium proper (two light locks, ±q): small v_t capture scan; persistence ≥20 periods or two-sided static minimum; measure E_bind.
- **Days 11–12:** charts as code: `view_M(state)`, `view_E(state)` post-processors; `export_grid(state)` → SFA for volview morphology.
- **Days 13–14:** `v81/FINDINGS.md`; run every kill gate explicitly; only then decide whether to spend the kernel authorization (locks block in scp_sim) or pivot to B's charge sector.

Deliberately deferred: hydrogenoid with the heavy lock replaced by an actual C-nucleus export, GPU port, any sfa.h change.

## What NOT to do

- **No multi-fab reruns** — no Z6+L6, no soft-v_t orbit grids, no third L-species tuning pass (per brief; the failure is representational).
- **Don't spend the kernel authorization first.** Authorization ≠ obligation. Modifying scp_sim before a sandbox demonstrates the state type works would weld the new ontology onto the old one under schedule pressure — the exact "extra arrays" failure mode.
- **Don't let c vary as ontology.** The v76 GRIN kill stands: c_eff belongs only in M-chart readouts.
- **Don't define the electron as a small Q-ball** — including through the back door of "initialize the lock from a radial_qball profile." Lock rest energy comes from its sequestration ledger, full stop.
- **Don't add a pairwise Coulomb formula** to make early force tests pass. The moment F(D) is computed from X_i−X_j instead of from the medium, the design is dead and should be reported dead.
- **Don't trust Eulerian-window diagnostics for binding** (the S4 lesson). Bound/free status is read off lock tracks and ledgers, never off Q_core in a fixed box.
- **Don't target carbon.** Positronium, then hydrogenoid, then nothing else until both are green.
- **Don't build the shape catalog.** One blend, one sandbox, kill gates run honestly, two weeks.
