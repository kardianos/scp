## Timeline reading (5–10 bullets)

- v74→v78 progressively engineered *recipes* (Q-ball profiles, ladders, R1–R10 freeze) on an unchanged state type: continuum multiplets on fixed N³. The R-matrix captured relations; the simulator representation never moved. v79/v81-B1 is where that debt came due.
- The B1 F-pass is real and worth keeping as a *target spec*: opposite-sector attraction with monotone \(a_{\mathrm{rel}}(D)\) exists in this kernel family. Force was never the blocker.
- All three O-failures are one failure at different speeds: the light charge is a field hump whose charge \(Q=\int|\phi|^2\) and \(E_{\mathrm{em}}\) are *soft integrals* — nothing structural stops them leaking. Soft \(T\sim4000\) evaporation and \(v_t=0.5\) shredding in \(<0.5\) rev are the same defect (no sequestration) at two energy inputs.
- v79 is the sharpest signal: hand-placed L6 didn't merely fail to bind — it *nulled* \(E_{\mathrm{em}}\) while C/Q tracked the L=0 control exactly. The third fabric wasn't inert emulation; it actively damaged the EM sector. That rules out "multi-fab is a passable shadow" as a fallback position.
- S4's real short/medium hold (once the COM artifact was removed) shows the attractive core dynamics is salvageable physics — but it lived in a window where the light hump hadn't yet radiated. Durability, not binding strength, is the quantity the next state type must buy.
- The v80 thesis is correct about the bottleneck being state type, but a wholesale "new substrate" is not week-scale implementable. The blend (4+3+6+10) has one piece that is both literal and cheap: **charge as Gauss defect on the free gauge medium (shape 10)**. That is also precisely the piece that fixes durability — an integer defect charge cannot evaporate.
- Chart machinery (C/M/E) is not a substrate; it is ledger discipline over one state. It can and should be implemented as readout from day one at near-zero cost — it is what keeps any candidate honest about the free/bound split and prevents quietly reintroducing parallel arrays.
- Sequencing consequence: every recorded failure mode is a *durability-of-light-charge* failure. The first experiment of any redesign must therefore test **lock durability and ledger honesty**, with binding as the second metric — not the reverse.
- Two revolutions at once (new charge representation *and* new geometry/path-measure dynamics) is how this program previously produced untestable layers. One literal change per stage; path-measure \(\ell\) enters only after locks are proven.

## Three redesign options (or ranked paths)

### Option 1 (rank 1) — Locks on a link-gauge free medium: defects + Maxwell substrate

  Shapes 10 + 4, with 6 as the literal step coin; 3 deferred as an extension slot.

- **State vector.** Two *kinds* of state, one world:
    - Free medium: gauge links — \(A_i(x)\) (or \(E_i,B_i\) on a Yee lattice), one shared U(1) sector, evolved hyperbolically at fixed \(c\).
    - Locks: a short list \(L_k=\{\mathbf{x}_k,\mathbf{p}_k,q_k\in\{\pm1\},E^\star_k,\text{type}\}\). The light charge *is* a lock carrying integer Gauss charge with a conserving string/current deposit — not a hump, not a Q-ball, not a third fabric.
- **Step law.** Per tick, \(\mathrm{d}t\) set by the hop budget \(c\) (CFL promoted to engine law, shape 6): E/B updated from curl; locks pushed by locally interpolated fields (Boris-type); current deposited with an exactly charge-conserving scheme (Esirkepov-type) so Gauss' law holds to machine precision every tick; lock ledger exchange satisfies \(\Delta E^\star+\Delta KE = -\int\!\mathbf{J}\cdot\mathbf{E}\). No pairwise force is ever inserted.
- **Conserved quantities.** Total charge (exact, integer per lock, Gauss-law protected); total energy = field energy + lock ledgers (to discretization error); approximately momentum. The free/bound split is a re-partition of one ledger — charts C/M/E fall out as readouts (C: fixed-\(c\) medium view; M: fixed lock labels; E: fixed region budgets).
- **v80 mapping.** 10 literally (charge = defect flux, no first-class charged \(\Phi\)); 4 literally (two kinds of state, isolation by *type*); 6 literally (per-tick causal spend); 3 as an empty slot (locks can later source an \(\ell\) field). Multi-fab intent survives as *lock species* (bag/charge/light) on one medium — exactly the REPRESENTATION.md architecture sketch. Passes all five v80 kill gates structurally: state is not \(\phi\) copies; particle is sequestration of budget with exact ledger; charts are re-partitions; no Q-ball precursors; \(c\) is the step coin.
- **Kernel delta.** New module `scp_locks.h/.c`: lock array (SoA), field→lock interpolation, lock→link current deposit, pusher, ledger diagnostics appended to `diag.tsv`; hooks into the existing Maxwell/U(1) tick after the E-update; config keys `n_locks` + per-lock `{type,q,Estar,x,v}`. Existing Cosserat \(\phi/\theta\) machinery untouched; SFA export becomes a chart projection of one state.
- **Sandbox first.** Yes — build the lock/pusher/deposit numerics standalone (plain C, ~1 week) where time-stepping, interpolation order, and subcycling can be tuned without kernel weight; port into `scp_sim` once the kill-gate scorecard is green. The medium update itself can reuse known-good kernel code at port time.
- **First experiment (minimum).** 3D, \(N=64^3\), \(L=32\), \(\mathrm{d}t=0.25\) CFL, \(T=2000\). Two locks \(q=\pm1\), \(E^\star\) equal, placed at rest, separation \(D=16\). Observables: \(D(t)\), \(E_{\mathrm{em}}(t)\), lock \(KE+E^\star\) ledgers, \(E_{\mathrm{total}}\) drift, flux through a detector shell. Success metrics:
    1. Rest-pair attraction reproduces the F-pass spec: monotone \(a_{\mathrm{rel}}(D)\) over \(D\in\{8,12,16,20,24\}\), with sign and scaling consistent with the emergent field force.
    2. Soft long-\(T\) run: inspiral proceeds with locks structurally incapable of evaporation; \(E_{\mathrm{em}}\) floors at the Gauss field of the surviving charges (never nulls); \(E_{\mathrm{total}}\) drift \(<2\%\).
    3. Hard \(v_t=0.5\) orbit: pair survives \(\ge 5\) revs with measurable radiation-reaction inspiral; no shredding is *possible* by construction, so the observable becomes the honest inspiral rate.
- **Kill gate.** Dies if binding/capture requires an inserted pair potential or hand damping (medium is fake — shape 4's own named risk); dies if grid-Cherenkov/discretization noise dominates the two-body dynamics at the scales of interest *after* standard remedies (higher-order shapes, subcycling) — that would mean the budget law isn't controlling the physics.
- **Honesty clause.** This is the oldest split in computational physics — fields plus particles — and the philosophy-paint risk is real: one could call any PIC code "free medium + locks." What makes it a v80 move rather than paint: (a) the *failed* object (Noether-charge hump) is removed as the definition of the light sector; (b) one shared medium, species as lock types, no parallel arrays; (c) C/M/E ledgers computed from one state as first-class output; (d) Stage-4 path replaces bare charge-locks with bag-locks sourced by the existing Cosserat core rather than adding a fourth fabric. The residual imported element is the lock's \(E^\star\) inertia label — that is the sequestered budget, and it is explicitly ledgered rather than hidden in a mass parameter of a hump.

### Option 2 (rank 2) — Dynamic path-measure medium: \(\ell(x,t)\) as stored structure

  Shape 3 primary, with 10 for charge and 6 for the step law.

- **State vector.** Free medium: scalar path-cost/capacity field \(\ell(x,t)\) plus the same gauge links as Option 1. Locks as Option 1, but additionally *sourcing* \(\ell\) (thickening). Signals and gauge dynamics propagate on the \(\ell\)-weighted metric at fixed coordinate \(c\) — the C-chart is native; the M-chart optical redraw (\(c_{\mathrm{eff}}=c/n\)) exists only as readout.
- **Step law.** \(\ell\) obeys a *hyperbolic* sourced law, e.g. \(\ddot\ell - c^2\nabla^2\ell = -S(\rho_E)\) with \(S\) from lock+field energy density (hyperbolic is mandatory: an elliptic per-tick Poisson solve is the named shape-3 collapse). Gauge links evolve with \(\ell\)-weighted curls; locks feel \(\ell\)-gradients in addition to Lorentz-type force.
- **Conserved quantities.** As Option 1, plus an \(\ell\)-sector energy flux that must be included in the E-chart ledger or the budget book will not close.
- **v80 mapping.** The only option where the *stored object* is monist path-cost (R4 primary, R1 native). Directly implements "warp at fixed local \(c\)" as data.
- **Kernel delta.** Everything in Option 1 plus one scalar hyperbolic field and weighted discrete curls — moderate in code, large in validation.
- **First experiment (minimum).** Static spine first: pin one heavy lock at box center, \(N=64^3\), \(T=500\); let \(\ell\) relax hyperbolically; check the radial \(\ell(r)\) profile is stationary and non-runaway; then release a light test lock and measure trajectory deflection versus the fixed-\(c\) budget prediction. Only if the spine passes: two-lock positronium repeat of the Option-1 experiment with \(\ell\) active.
- **Kill gate.** Dies if \(\ell\) requires an elliptic solve per tick to remain sane (it is then a side-Poisson after \(\rho_b\) — the exact collapse SHAPES.md warns of); dies if removing \(\ell\) leaves dynamics statistically unchanged (decorative field = \(\phi(x)\)-with-extra-arrays in a new costume); dies if \(\ell\) back-reaction runaways (attractive collapse) without a mechanical stabilization that isn't a fudge.
- **Honesty clause.** Highest conceptual fidelity to v80, but it doubles the number of unproven changes versus Option 1, and its failure modes (passive scalar, or GRIN rerun) look exactly like places this program has already been. v76's kill is respected only if \(\ell\) makes *no* gravity claim — it is medium state, and the optical redraw stays in the M-chart.

### Option 3 (rank 3) — Token/currency cellular substrate: \(c\) as literal budget, particles as medium defects

  Shape 6 literal + 4 architecture; the only option where no field is imported at all.

- **State vector.** One conserved currency: token density \(u(x)\) and link fluxes \(f_i(x)\) with a hard per-tick hop cap \(c\). Locks are not objects in the state — they are *self-sustaining circulation patterns* (vortex defects) of the token flow; charge = circulation handedness \(\pm\); mass/inertia = sequestered update budget (tokens stuck in the pattern). Positronium = vortex–antivortex pair of one substrate.
- **Step law.** Conservative local transfer rule (lattice-gas style): tokens move along links, total tokens exactly conserved, per-link per-tick spend \(\le c\); vortices advect under token-pressure gradients.
- **Conserved quantities.** Token total (the budget itself); vortex winding (charge, topological); approximately energy.
- **v80 mapping.** The purest realization of "free measure as world": free/bound is literally a re-count of the same tokens; kill gates 1, 2, 5 are passed by construction.
- **Kernel delta.** None — this must be a pure 2D sandbox spike; it shares no numerics with `scp_sim`.
- **First experiment (minimum).** 2D, \(N=512^2\), \(T=10^5\) ticks. Nucleate a vortex–antivortex pair; measure interaction-force sign and radial scaling versus separation; check long-lived bound gyration with token-conservative radiation.
- **Kill gate.** Dies at two weeks regardless of promise if the defect pair does not show clean attraction with a power-law-ish profile, or if the circulation patterns evaporate under their own radiation — which would reproduce, in the new substrate, the exact hump-durability failure it was built to escape.
- **Honesty clause.** This is a research direction, not a Stage-3 deliverable: getting Coulomb-type forces and radiation at \(c\) out of a discrete currency is unproven, calibration to the R-matrix is far, and the probability it becomes the atom path on a week scale is low. Its value is as a bounded probe of whether *any* single-currency substrate supports durable opposite defects — a question no amount of multiplet tuning can answer.

## Recommended primary path + first 2 weeks of work

  **Primary: Option 1** (locks on the link-gauge medium), sandbox first, kernel port second. Option 2's \(\ell\) field is pre-allocated as an extension slot in the state design and is the natural Stage-3.5 once locks are green; Option 3 runs only as a time-boxed 2D spike, in parallel at most, never in the kernel.

  Week 1 — sandbox core:
- Days 1–3: standalone 3D Yee E/B solver, CFL budget as step law, lock list, Boris push, charge-conserving deposit. Sanity: single accelerated lock radiates; radiation exits the box; E-chart ledger (field + lock + shell flux) closes to \(<1\%\). If the ledger doesn't close here, nothing downstream means anything.
- Days 4–7: rest-pair scan \(D\in\{8,12,16,20,24\}\): monotone \(a_{\mathrm{rel}}(D)\) from emergent field force only; then the soft long-\(T\) run (\(T=2000\), the configuration that evaporated at \(t\sim1500\)): locks persist by construction, \(E_{\mathrm{em}}\) floors at the surviving Gauss field, inspiral tracked, drift budget written down.

  Week 2 — the shredding test and the port decision:
- Days 8–10: \(v_t=0.5\) circular-orbit initialization, \(\ge 5\) revs; measure radiation-reaction inspiral rate against classical expectation; tune interpolation order/subcycling if grid noise is non-negligible (remedy, not death).
- Days 11–12: implement the C/M/E chart readout module as first-class output (C: medium view at fixed \(c\); M: lock ledger; E: region budgets with shell-flux honesty — the v79 \(E_{\mathrm{total}}\) drift lesson).
- Days 13–14: kill-gate scorecard written against the three success metrics; if green, port as `scp_locks.h/.c` into `scp_sim`, reusing the kernel's existing Maxwell tick and SFA export as the chart projection; verify parity against sandbox on the \(D=16\) case.

  Day-14 deliverable: a positronium-analog two-lock bound/inspiraling system with intact locks and an honest energy ledger, or a written kill naming which gate fired. No carbon, no Stage-4 hooks until this is green.

## What NOT to do

- Do not re-run multi-fab Z6+L6, soft \(v_t\) orbit grids, or any "park one more L variant" campaign — that state type has been falsified for this purpose three times (evaporation, shredding, \(E_{\mathrm{em}}\) nulling).
- Do not add a fourth fabric, per-species gauge copies, or per-species \(\ell\) — parallel arrays are the emulation habit in new clothes; isolation must come from lock *type* on one medium.
- Do not hand-place Q-balls, Gaussians, or any field hump as the *definition* of the light sector — that is v80 kill gate 4 and the direct cause of all three O-failures.
- Do not inject pairwise Coulomb forces, springs, or damping "to get positronium" — the moment the medium is bypassed, the experiment tests nothing (shape 4's named risk, and Option 1's kill gate).
- Do not put \(\ell\) (Option 2) into the kernel before locks are proven — two unproven state changes at once is how failures become undiagnosable.
- Do not let \(\ell\) become an elliptic per-tick solve or a gravity claim — side-Poisson collapse and the v76 GRIN kill both still hold; the optical \(c_{\mathrm{eff}}\) redraw stays M-chart readout only.
- Do not build the CA/currency substrate inside `scp_sim` — it shares no numerics with the kernel; sandbox only, two-week box, hard gate.
- Do not spend the redesign budget on chart formalism alone — C/M/E is a readout module riding on a real state change, not a deliverable by itself; ledgers without a new state type are v78's gap restated.
- Do not aim the first experiments at binding energy spectroscopy or carbon compatibility — the first metric is durability + ledger honesty, because that is what actually failed.

