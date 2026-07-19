# SCP Next-Step Recommendations (v80 era)

## A. Timeline reading

- **Stage 2A is genuinely closed.** Z-carbon parks; the nucleus template (gauged Q-ball) is solid. This is the strongest result in the campaign and nothing recent challenges it.
- **v78's atom package was ADOPTED but the C₁₂ product was already BLOCKED** — the campaign entered v79 knowing the naive assembly path was broken. v79 confirmed *why*: a hand-placed L6 shell acts as a neutralizing/radiation bath, collapsing E_em 6.13→0.36 while C/Q bookkeeping is untouched. The gauge ledger (Gauss at 1e-13 floor) is reliable throughout, so bookkeeping integrity is not the bottleneck — field dynamics is.
- **v79 vs v80 S4 is the single most informative contrast in the brief:** fat Z6+L6 kills E_em, minimal C+L hydrogenoid-class holds it (0.71→0.75 over T=800). Shell *loading*, not the multi-fab architecture per se, is what kills the EM sector. This is a mechanism lead, not yet a mechanism.
- **v80's G = 0.62 clears the continue threshold but is soft everywhere it matters physically.** S_ledger = 1.00 is a consistency check, not physics. The load-bearing components — force (0.55), orbit (0.50) — rest on indirect proxies (E_em comparisons at fixed D, E drift) because the SFA pair tracks were pruned and multi-rev was never measured.
- **E drift is a red flag in both "positive" v80 results:** +14% at rest, −27% at vt=0.08. Held E_em with drifting total E means the "bound" configurations are exchanging energy with something (radiation, grid, each other) — the campaign currently cannot say which.
- **Orbit remains unestablished.** v74–75: partial classical-ish capture, head-on merges. v80: mild low-vt behavior, vt=0.08 bleeding energy, no multi-rev measurement. Nothing since v75 has moved this needle.
- **The v80 representation thesis is untested design work.** The external review narrowed it sensibly (grid innocent; density-only state the problem; 2D toy first), but zero lines of substrate code exist. It is a paper tiger until the 2D toy runs.
- **S_morph = 0.25 and the pruned SFAs are self-inflicted evidence gaps.** The campaign scored itself on morphology it never captured and force curves it deleted. Instrumentation discipline, not physics, cost the most information per dollar in v80.

## B. Long-term goal status

- **Stage 3 (light opposite-charge sector) is the active blocker and is partially open.** The C+L two-body problem shows: long-lived coexistence (S4), energy leaks of unknown sink, no measured Coulomb acceleration curve, no orbit. Positronium-first remains the right gate and is not passed.
- **Stage 4 (C + 6 bound light charges) is not approachable yet, and v79 shows the shortcut is fatal.** Skipping Stage 3's two-body mechanics to hand-place shells produces neutralization, not atoms. The Stage 4 gate should be explicit: measured attractive a_rel(D) over a D-range + at least one multi-rev bound orbit with bounded E drift.
- **Stage 5 (spontaneous production) is untouched and should stay untouched.**
- **What would count as real progress toward carbon:** in priority order — (1) a measured, gauge-clean attractive force law a_rel(D) for C–L with error bars; (2) one closed multi-revolution orbit with |ΔE/E| bounded and the sink identified; (3) a principled explanation of why L6 collapses E_em while L1 holds it (this determines whether "6 electrons" is even reachable by loading); (4) only then, staged shell loading 1→2→…→6 with E_em survival as the per-step gate. Note that none of these require the representation substrate — they are product-line mechanics.

## C. Recommended next step (primary)

**Instrument, then measure the C–L two-body force and orbit properly — one overnight GPU run, no kernel changes.**

Step 1 (short GPU, ~1–2 h sanity): Re-run a minimal C–L pair at D ∈ {8, 12, 16, 24} from rest, N=192, T≈200, **with `mf_pair_track` SFA output retained and volview morphology capture enabled** for at least one D. Confirm SFA writing works end-to-end and storage fits (`/space/scp/`).

Step 2 (overnight GPU, the main event): Two grids on the V100-class instance:
- **Force grid:** C–L from rest, D ∈ {6, 8, 10, 12, 16, 20, 24, 32}, T long enough for measurable displacement; extract a_rel(D) from pair tracks and test against Coulomb (1/D²) and alternatives (Yukawa-like, constant). Include same-sign controls at the same D values (repulsion baseline) and an L=0 control (fabric-only background drift subtraction).
- **Orbit grid:** vt ∈ {0.02, 0.04, 0.08, 0.12} at the two D values where step-1 shows the mildest behavior, T ≥ one estimated classical period at the measured a_rel (use step-1 force law to set this — do not guess vt from v80's failed 0.08). Track separation vs time for multi-rev classification.

Step 3 (CPU, next morning): Fit force law; classify each orbit run as bound / merging / escaping; decompose E drift into EM, kinetic, and radiation-proxy channels; produce the volview morphology for the best orbit candidate.

**Why this before alternatives:** Every fork in the campaign's decision tree — product ladder vs representation pivot, Stage 3 pass/fail, shell-loading feasibility — depends on whether C–L actually attracts with a clean inverse-square law and can hold an orbit. v80 tried to score exactly this and couldn't, because the tracks were pruned and morphology never captured. Re-running the same physics *with the instrumentation fixed* converts v80's soft 0.5–0.55 scores into a hard yes/no. The representation substrate (line B) cannot be validated without a product-side ground truth to emulate anyway, so it is strictly downstream of this measurement.

**Success criteria:** a_rel(D) fit with R² > 0.9 to a monotone attractive law for opposite sign (and repulsive for same sign); ≥1 vt/D combo with ≥2 full revolutions and |ΔE/E| < 10% with identified dominant sink; E_em held (no v79-style collapse) throughout.
**Kill criteria for multi-fab two-body atoms:** a_rel flat/noise across D, or every orbit run merging/escaping with E_em drain — that *would* be "clear force death" and triggers the soft-kill clause and a full pivot to line B.
**Effort:** ~1 h config/seed work + short GPU sanity + overnight GPU (~6–8 h) + half-day CPU analysis. **Kernel auth: no** — everything runs through existing generators (`gen_mf_pair_boost`, `gen_mf_shell_orbit`, `mf_pair_track`) and configs. The only change is operational: stop pruning SFAs and budget disk on `/space/scp/`.

## D. Explicitly defer or kill

- **Do not re-run Z6+L6 (or any multi-L shell park) as the next experiment.** v79 already answered it: fat shells neutralize. The only legitimate version is *staged* loading gated on per-step E_em survival, and that belongs after the C–L force law is measured. The brief's instruction stands.
- **Do not start the kernel rewrite for the representation substrate.** The external review said 2D toy first; writing kernel code before the two-body ground truth exists is building an emulator of an unknown.
- **Do not run more G-score composite campaigns.** G = 0.62 demonstrates the composite's failure mode: it averaged away the missing a_rel and missing morphology behind a perfect ledger score. Replace with hard per-component gates; retire G as a decision device.
- **Do not run morphology-free or SFA-pruned GPU jobs again.** Any run that cannot produce pair tracks and at least sampled volview frames is wasted wall-clock at this stage.
- **Defer Stage 2B (multi-center carbon) and anything Stage 5.** Neither unblocks Stage 3.
- **Defer the RC1 co-field / dual-channel extensions** unless the force grid shows a specific pathology they address.

## E. Product vs representation split

Over the next 1–2 calendar weeks: **~80% product line A, ~20% representation line B**, structured as:

- **Line A (product, GPU-critical-path):** the force/orbit grid of section C, plus one targeted diagnostic the v79/v80 contrast demands — a staged-loading ladder C+L₁, C+L₂, C+L₄ at fixed N and T=400–800 recording E_em per step. This isolates where E_em death begins as a function of shell count and is cheap (three short runs). Line A produces the ground-truth dataset.
- **Line B (representation, CPU-only, no kernel):** implement the 2D free-medium + typed locks toy exactly as the external review scoped it — path-measure state, dual free+defect decomposition, c-as-update-budget — as a standalone Python/NumPy toy, not in `scp_sim`. Its explicit validation target is to *reproduce the a_rel(D) curve measured by line A*. If line A hasn't delivered the curve yet, line B builds the toy skeleton and unit-tests conservation bookkeeping, but claims nothing physical.
- Hard rule: line B gets zero GPU and zero kernel time this fortnight; line A gets zero refactoring. Re-evaluate the split only when the force law is in hand — a clean Coulomb result strengthens the emulation thesis (representation becomes the scaling path to shells); a dead force law makes line B the main line.

## F. Confidence

**3/5.**

The recommendation direction (measure the two-body force law properly before anything else) is robust — it follows from the campaign's own evidence gaps and needs no new physics. What caps the score:

- The pruned SFAs mean I cannot verify that "elite opp Q≈0 / same Q~230" actually implies the force asymmetry attributed to it; the entire
