# KIMI_FEEDBACK.md — External review of v80 (kimi CLI, model k3)

## A. Steelman

The thesis, restated: **a simulator can only say what its state type can store.** v78 froze a ledger of relations (locality \(c\), \(E \leftrightarrow Mc^2\), free/bound budgets, dual sources, Coulomb) that presupposes a two-kind ontology — a free medium, and bound *sequestrations* of it. Every implementation so far stores exactly one kind of thing: field amplitude at grid sites, times N fabrics. In that state type, "bound" can only be superposition, "isolation" can only be weak coupling, "particle" can only be a hump. The central monist claim — matter is free structure sequestered — has no data structure anywhere in the program. v79's null atom is then the predictable symptom, not a surprise: a hand-placed hump in a third copy of the same type cannot be an "opposite-sector lock"; it can only cancel \(E_{\mathrm{em}}\) through the shared \(A\), which is what it did.

Is it coherent? **As a critique, yes.** It is the standard representational argument (Eulerian density vs. particles in fluids; collective coordinates vs. raw field for solitons), applied to SCP's own gap. The strongest evidence for it is not v79 but the admission in `REPRESENTATION.md` itself: "v78 can *score* free/bound; the code always *is* floats at sites." That is checkable and true independent of any run.

**As a program, not yet.** The load-bearing nouns are undefined: *free substrate* (what object? graph? measure? links?), *sequestration* (what dynamical event creates a lock?), *chart* (view, constraint, or frame?). v80 is a coherent thesis with placeholder vocabulary.

## B. Agree / disagree

**Strong:**

- Multi-fab is N copies of one representation. Correct and worth saying hard: adding replicas of a type is not a new type. v79's failure mode (opposite sector can only couple through the shared gauge field it was supposed to be isolated from) is structural, not parametric.
- Relations ≠ substrate. R1–R10 are constraints; they do not define a step operator. Correct.
- Charts must be re-partitions of one state, and \(c_{\mathrm{eff}}\) must stay a computed readout, never an evolved field. Correct, and consistent with the v76 kill.

**Wrong or overstated:**

1. **"Fixed grid" is a red herring; state-kind is the real axis — and the docs keep conflating them.** Every candidate shape will be implemented on a discrete lattice or graph: shape 3's \(\ell\) will be an array, shape 6's budget will be an array, shape 10's links *are* a lattice. Any simulator is discrete; "not a fixed Eulerian grid" as written filters *implementations*, not *ontologies*. The real filter is: does the state contain typed, budget-bearing objects with identity, or only densities? And: does the coupling law permit sectors that interact with the medium but not directly with each other? You could fail every kill gate off-grid and pass all of them on one (shape 4 on a fixed mesh + discrete lock objects). Rewrite the filter accordingly or shape 3/6 will "pass" by renaming.
2. **The v79 → representation inference is weak.** One hand-placed recipe failed. That is also consistent with wrong initial conditions, wrong coupling strengths, no capture dynamics, or insufficient T. v80 admits this ("fair test of the recipe, weak test of monist field nature") and then pivots the whole program on it anyway. The honest claim is "v79 is *consistent with* representation inadequacy," not evidence for it.
3. **R1–R10 underdetermine all ten shapes.** "Relations→field fails" is right, but relations→*substrate* is equally unconstrained: ten shapes fit the same ledger. Without a quantitative anchor (a new substrate must reproduce v78 R-matrix numerics on export), shape selection is taste. Say this plainly in the docs; it is the sharpest cargo-cult guard available.
4. "Inexpressible" is too strong. Locks *can* be represented on top of fields (moduli-space/collective-coordinate reductions do exactly this). The true claim is "unnatural and structureless," not impossible.

## C. Shapes 3, 4, 6, 10

**Shape 3 — path-measure primary.** Promise: makes "warp at fixed \(c\)" the stored object; eikonal/travel-time is the right language for a free medium. Nearest existing ideas: Regge calculus (edge lengths as the degrees of freedom — the strongest match), Finsler/optical (Gordon) metrics, analog-gravity acoustic metrics, fast-marching eikonal solvers. Fatal risk: \(\ell(x)\) as a scalar field on fixed N³ with hyperbolic evolution *is* \(\phi(x)\) with the units renamed — highest collapse-to-\(\phi\) risk of the four. It escapes only if the state is genuinely relational (lengths on evolving connectivity, i.e. shape 1), not a measure-field on a fixed mesh. **Collapses back: yes, by default, unless implemented as edge data on dynamical connectivity.**

**Shape 4 — dual substrate (free medium + typed locks).** Promise: gives particles identity, type, and ledger; recasts multi-fab intent correctly (species of lock, one medium); matches how simulators are actually architected. Nearest existing ideas: PIC / hybrid particle-mesh methods (this *is* a dual substrate — mesh fields plus discrete particles with cross-forcing); soliton collective-coordinate reductions; moving punctures in numerical relativity (punctures = locks in a metric field). Fatal risks, plural: (a) becomes N-body with hand potentials if the medium is decorative; (b) two-way coupling consistency (self-force, lattice Cherenkov artifacts, energy drift) is exactly where hybrid schemes die — the PIC literature is the cautionary file; (c) **locks have no origin**: v80 asserts "a lock exists if free substrate supports stable sequestration" and provides zero sequestration dynamics. Locks are assumed, not produced. If locks are inserted by fiat, you have re-imported chemistry through the type system. Also note a real tension: a structureless lock coupled to a continuum medium radiates and decays unless topology protects it — which pushes 4 into needing 10. **Collapses back: no — it escapes \(\phi(x)\) by construction; its failure mode is dualism-by-API instead.**

**Shape 6 — update-budget monism.** Promise: takes R1 literally as a computational coin; "inertia as stuck updates" is a genuinely different intuition. Nearest existing ideas: cellular automata / quantum CA / digital physics — and, uncomfortably, *your own explicit integrator*: the CFL condition already is a per-tick causal spend limit with \(c\) as hop radius. Fatal risk: either trivial (rebranding explicit timestepping) or false (if "budget" is meant as more). The burden is to exhibit one phenomenon the budget language explains that the scheme didn't already have — e.g. it must say something predictive about lock acceleration. If it can't, kill it first; it is the most metaphorical of the four and the most likely to produce numerology. **Collapses back: medium-high — as a per-site budget array it is \(\phi(x)\) again; it survives only as edge/event bookkeeping, at which point it has merged into shapes 1/5.**

**Shape 10 — gauge/defect first.** Promise: charge as topology (winding, frustration, holonomy) is the one place where "not a hump in \(\phi\)" is already achieved in known physics — defects are classified by homotopy, not amplitude; Gauss becomes counting. Nearest existing ideas: Wilson-link lattice gauge theory (connections on links *are* the primary state — and note it lives on a fixed lattice, further proof the grid filter misfires), topological defects/Kibble, BKT vortices, compact QED monopoles, discrete exterior calculus. Fatal risk: alone it supplies charge without mass, inertia, or lock dynamics; getting stable, mobile, inertial defect particles from pure gauge is unproven and hard (the known success case is confinement → hadrons, i.e. reinventing QCD). Honest caveat for the toy scale: in 2D, gauged vortices need an order parameter — a complex scalar on the lattice, i.e. \(\phi(x)\) returns through the Higgs door — while pure-gauge charges are hand-placed external sources, which violates kill gate 4. This fork is not a detail; it is the experiment (see F). **Collapses back: low for the charge sector itself, but watch the order-parameter backdoor.**

**On the blend (4+3+6+10):** it is a wish-list composite, not a design. Each shape's fatal risk is "covered" by another shape's promise, which is how incoherent monoliths get born. The defensible minimal pair is **4+10** — topology supplies lock stability, locks supply the typed ledger — with 3 demoted to a readout (travel-time analysis of the medium) and 6 to a naming convention. Pick one load-bearing idea; the union of liked features is not an architecture.

## D. Dual charts (C/M/E)

"Reality as real + POV" is **useful as a software doctrine, trite as physics unless formalized**. As software, it is just state + views — and it is exactly the discipline that would have prevented the v79 error (reading an export as the world). Keep it for that.

Concrete legitimacy test: a chart system is coherent iff chart transitions are **invertible maps on one state**, adding and losing nothing. M-chart passes *if* \(c_{\mathrm{eff}}\) is always derived from the C-chart state and never evolved — the optical/Gordon-metric rewrite of a fixed-\(c\) theory is standard and fine as numerics. The moment any chart residual (\(M_{\mathrm{eff}}\), \(c_{\mathrm{eff}}\)) acquires its own update law, it is a new field and the doctrine has failed. State this rule in the docs; it is one sentence and it kills the "variable \(c\) as gravity" licensing risk outright.

**E-chart is the soft spot.** "Hold a region's energy budget constant" is not a coordinate choice; it is a constraint on an open subsystem — thermodynamics, not a frame. Treating it as a peer of C and M invites bookkeeping gauge to masquerade as physics. Demote E to an analysis view with explicit coarse-graining labels.

**The lock triple \(E_\star = M_\star c_\star^2\) is numerology as written.** If \(c_\star\) is the universal free \(c\), then \(M_\star\) and \(E_\star\) are one dial, not two, and the "only two free" framing is vacuous. If \(c_\star\) is per-lock, you have reintroduced per-particle light speed through the back door — the exact thing v76 killed. Fix the table: one global \(c\), locks carry one ledger number plus type.

And a blunt correction to `DUAL_CHARTS.md`'s v79 diagnosis: "L nulling \(E_{\mathrm{em}}\)" is consistent with "extra multiplet on shared A," yes — and equally consistent with "right ontology, wrong coupling." Charts do not diagnose which layer failed. Only a state-level experiment can.

## E. Multi-fab

**Demote to export/projection and benchmark harness. Do not redesign it; do not run new atom parks on it.**

Concretely:

1. **Freeze all monism claims from multi-fab outputs.** v79 showed the state type has no slot for "opposite lock," so no atom produced there can constitute evidence about monism. It remains a real product API for stable nuclei — as engineering, not as ontology.
2. **Keep the v79 T=800 core/control data as the regression set.** Define the export contract now: any future substrate must emit \((Q_\phi, Q_{\mathrm{flux}}, E_{\mathrm{em}}(t), r_{\mathrm{core}})\) comparable to the L=0 control within a stated tolerance. This is the quantitative anchor B(3) says you currently lack.
3. **One allowed product use:** acceptance target for a new substrate's `export_grid()`. If a product atom is still wanted in parallel, it belongs to the product line (recipe work: dynamical capture instead of hand-placement), explicitly quarantined from the representation line.

Redesigning multi-fab to "fix" v79 would be the worst option: more effort into the representation you just argued is an emulation.

## F. What to build next

**One toy: 2D U(1) link medium + two typed locks, standalone, ~a few hundred lines, no `scp_sim` dependency.**

- **Free medium:** U(1) link variables on a 2D lattice, free Hamiltonian/wave update at fixed \(c\) (shape 10's free part).
- **Locks:** exactly two typed objects — id, charge ±, position, phase clock, budget ledger (shape 4). They source link flux (Gauss constraint) and are forced by the medium. Charge from winding/orientation, never from the sign of an amplitude.
- **Charts:** implement `view_C` (travel times, flux at fixed \(c\)), `view_M` (derived \(c_{\mathrm{eff}}\) redraw), `view_E` (budget ledger), and **numerically test chart-transition invertibility** — this directly tests the D doctrine's legitimacy, not just the physics.
- **The experiment is the fork named in C:** attempt stable mobile ±locks *without* (a) hand-placed external sources and *without* (b) an order-parameter hump field. Pre-register the kill gates: locks need hand-tuned mutual potentials → medium is decorative, dead; attraction requires reintroducing a site-field soliton → collapsed to \(\phi(x)\); charts not losslessly invertible → chart doctrine incoherent. Any of these outcomes is a *result*, which is what makes this a test rather than a demo.

This is the smallest artifact that touches the actual gap (typed non-field state + isolation + charts), needs no carbon, no multi-fab park, and no kernel work. **Kernel rewrite is not required now** — the open question is "what is the state type," which a toy answers for days of effort; rewriting `scp_sim` before the toy falsifies or survives the 4+10 pair would be the premature move. The toy's output is precisely the API contract (`state / step / view / export`) that a later kernel rewrite would implement.

## G. Kill list — v80 should explicitly refuse to

1. Run any carbon/atom matrix, on any fabric count, as evidence for the representation thesis.
2. Touch `scp_sim`/`sfa.h` before the toy passes or fails its pre-registered gates.
3. Permit per-lock \(c_\star\) anywhere in any ledger. One \(c\), in C-chart; \(c_{\mathrm{eff}}\) derived only.
4. Permit chart residuals (\(M_{\mathrm{eff}}\), \(c_{\mathrm{eff}}\)) to have their own dynamics.
5. Accept any candidate whose state is an amplitude array on fixed N³ renamed "measure" or "budget." Renaming is not re-representation; shapes 3 and 6 die here by default unless proven otherwise.
6. Accept hand-placed locks as the *definition* of a light sector in any future atom claim (make kill gate 4 absolute).
7. Revive GRIN/optical-metric gravity, including rebranded as "M-chart gravity" (v76 kill stands).
8. Spend effort on shape 7 (spectral/label-primary): geometry-secondary toys maximize cargo-cult return for this program.
9. Interpret chart information loss as physics. A non-invertible view is a bug or a coarse-graining; label which.
10. Treat the 4+3+6+10 blend as "the design," and refuse to admit new shapes beyond the ten until one focus shape has died or passed.

## H. Score

| Item | Score 1–5 | Note |
|------|-----------|------|
| Thesis clarity | 4 | The diagnosis (relations≠substrate; multi-fab = N copies; free/bound has no data structure) is sharp and quotable. Docked one for the grid-vs-kind conflation and undefined core nouns (substrate, sequestration, chart). |
| Novelty-with-teeth | 2 | Every piece exists elsewhere (Regge calculus, lattice gauge, PIC, CA, optical metric); novelty is local to this program — which v80 allows, but teeth require one chosen shape plus a falsifiable test, and v80 ships a menu instead. |
| Risk of cargo-cult | 4 | High: renaming \(\phi(x)\) as \(\ell\) or "budget," per-lock \(c_\star\), E-chart metaphysics, and "reality as POV" language can all cosplay depth. The kill gates are genuinely good; the lock triple and E-chart are where it would sneak in. |
| Path to a runnable sandbox | 3 | The MVE is the right instinct but under-specified (no update rule, no pass/fail numbers, locks' origin unaddressed). Path is open precisely because kernel work is forbidden. |
| Fit to carbon atom goal | 3 | Indirect by design: a detour that could correct the question "what would an atom even be," but nothing here shortens the distance to C₁₂; best case is a better-posed Stage 3, not an atom. |

**Bottom line:** the thesis is right about *what* is missing (a typed, two-kind state) and wrong about *where* the enemy lives (the grid is innocent; density-only state is not). Build the 2D gauge-medium + two-locks toy with pre-registered kill gates, anchor everything to the v79 control data as export contract, and let exactly one shape die or survive before writing any atlas as if the blend were settled.
