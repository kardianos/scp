# SCP v80-era Next-Step Recommendations — Independent Advisor Report

## A. Timeline reading

- **The nuclear side is genuinely closed.** Stage 2A (Z-carbon liquid drop) is done, free A=12 super-criticality at g=0.05 is a measured constraint, and nothing since v74 has contradicted it. Nuclei are not the bottleneck.
- **v79 established a real negative result, not a failed run:** a hand-placed Z6+L6 shell is a neutralizing radiation bath. It nulls E_em (6.13 → 0.36) and Gauss flux (−68%) while leaving the C/Q core inventory identical to the L=0 control. "Place 6 light charges around a nucleus" is dead as a construction method at this coupling.
- **v80's single most important datum is the S4 contrast:** minimal C+L hydrogenoid-class holds E_em (0.71 → 0.75) over T=800 where the fat atom collapsed. Whatever binds atoms here, it exists — if at all — in the one-electron sector first. This validates the Stage-3 ladder ordering (positronium/hydrogenoid before multi-electron).
- **But the force law is unmeasured.** S_force=0.55 rests on completion + qualitative E_em ordering; the a_rel pair tracks were pruned with the SFAs. There is currently **no quantitative evidence that the C–L attraction is long-range Coulomb (∝1/D²)** rather than a short-range screened interaction. This is the load-bearing unknown for Stage 4.
- **Orbits are not demonstrated.** vt=0.08 shows −27% E drift and no multi-rev was measured. "Orbit/capture" remains a claim, not a result.
- **The S4 bookkeeping anomaly is unresolved:** Q_flux 114→4.7 and Q_core 114→0.46 "likely COM/window" is a hypothesis, not a check. If that drop is real charge loss, the S4 hold is much weaker than G=0.62 implies.
- **Ledger discipline is solid everywhere** (gauss_max at 1e-13 floor across v79 and all 16 v80 jobs). The gauge implementation is not in question; the physics of binding is.
- **The representation thesis (v80) is a design document with one external review, zero simulation states.** The kimi k3 recommendation (2D free-medium + typed-locks toy before any kernel rewrite) is the correct containment of that line.
- **Morphology is a blind spot** (S_morph=0.25, no volview on the entire campaign). Cheap to fix, and it matters: "E_em held" without knowing whether the L is orbiting, smeared, or sitting inside the core is ambiguous.

## B. Long-term goal status

**Stage 3 is partially open, better than pre-v80 but not passed.** What exists: opposite-attract/same-repel (E-lite, F8–F9); a single light L that coexists with a C core for T=800 without killing E_em. What does not exist: a measured long-range force law, a bound orbit surviving multiple revolutions, any multi-L configuration that is not a radiation bath.

**Stage 4 is blocked on exactly two measurements**, in order:
1. **F(D) ∝ 1/D² with the right sign** between C-core and L over at least a factor of ~4–8 in D. Without this, "Coulomb-bound" is a slogan.
2. **A hydrogenoid bound state with ≥2 revolutions** and bounded energy drift. That is the fabric hydrogen atom. Carbon-atom work before hydrogen works is skipping the stack — the same discipline that made Stages 0–2A honest.

**Stage 5 is not meaningfully closer** and should not be discussed until Stage 4 has a hydrogen.

**Real progress toward the carbon atom now means:** a measured force-law exponent, a multi-rev hydrogenoid orbit with morphology, and a resolved S4 bookkeeping story. Nothing else counts.

## C. Recommended next step (primary)

A tightly ordered 3-step sequence — inseparable because each gates the next.

**Step 1 (do first, cheap): Resolve the S4 Q_flux/Q_core windowing question.**
- Re-analyze existing S4 diagnostics with a COM-tracked integration window (centroid-following sphere, not a fixed box region). If the diag TSVs are insufficient, re-run S4 short (T≈200, N=192) with tracked-window diagnostics and volview snapshots at start/mid/end.
- **Why first:** if the drop is real charge loss, S_Lhold collapses and the whole G=0.62 continue-signal weakens; every downstream run inherits this measurement method.
- **Success:** COM-tracked Q_core stays ~constant → windowing confirmed. **Kill:** charge genuinely leaving → treat S4 as slow evaporation, tighten timescales before orbits.
- Effort: CPU analysis, at most short GPU. **Kernel auth: no.**

**Step 2: Re-run the pair force matrix with tracks retained — measure F(D).**
- `gen_mf_pair_boost` + `mf_pair_track`: C-core + single L at rest, D ∈ {8, 12, 16, 24, 32} (as large as N=192/L=48 allows; go to N=256 if edges bite), plus a same-charge control at 2–3 D values. **Keep the a_rel tracks and diag TSVs; do not prune before extraction.** Fit a_rel(D) to a power law and to Yukawa; add volview snapshots per run.
- **Why before orbits:** the pruned tracks are v80's biggest self-inflicted evidence gap, and the exponent decides everything. If the interaction is screened at D ~ core radius, no orbit parameter scan will ever find a hydrogen — you'd learn that fact expensively via Step 3 otherwise.
- **Success:** exponent within ~±0.3 of −2 over the range, opposite-sign attraction, magnitude scaling sanely with Q_L. **Kill:** clean exponential screening or exponent ≪ −3 ⇒ multi-fab L does not produce a Coulomb sector at g=0.05 ⇒ Stage 4 blocked in this architecture; escalate to the user (gauge-parameter discussion or weight shift to line B). Do not soft-pedal this outcome.
- Effort: short GPU (~5 runs × modest T; half a day wall). **Kernel auth: no.**

**Step 3 (only if Step 2 passes): Hydrogenoid multi-rev orbit.**
- `gen_mf_shell_orbit`, single L, initial D and vt chosen from the *measured* F(D) for a circular orbit (not scanned blind — that's why Step 2 comes first). N=192–256, T long enough for ≥3 predicted periods, COM-tracked diagnostics from Step 1, volview renders mandatory.
- **Success = the fabric hydrogen atom:** ≥2 revolutions, E_em drift bounded (say <10%/rev), D(t) bounded, ledger at floor. This is the single biggest possible advance toward Stage 4.
- **Kill:** inspiral/ejection within 1 rev at F(D)-derived parameters ⇒ radiation reaction dominates ⇒ next question is drag scaling, not more parameter scans.
- Effort: overnight GPU. **Kernel auth: no.**

## D. Explicitly defer or kill

- **Kill: any Z6+L6 or multi-L re-park.** v79 was decisive and no new mechanism exists. Multi-L returns only after a one-L orbit works, and then via sequential capture into a measured potential, not hand-placement.
- **Kill: composite goal-function campaigns as the primary instrument.** G=0.62 blends strong evidence (ledger) with absent evidence (pruned tracks, no morphology) into one reassuring number. The next runs should answer single questions with kept data.
- **Defer: monist kernel rewrite / primary free-substrate state.** Needs kernel auth, and the thesis has zero simulation backing. The 2D toy earns it or kills it first.
- **Defer: Stage 5, A=12 substructure (2B), abundance work.** Nothing upstream is ready.
- **Don't** silently drop SFAs again on any force/orbit run; pruning the pair tracks converted a completed campaign into qualitative evidence. Track files are tiny — keep them locally as policy.

## E. Product vs representation

**~75% product / ~25% representation over the next 1–2 weeks.**

- **Line A (product):** Steps 1–3 above are the critical path and are GPU-shaped. If Step 3 passes, follow with orbit robustness (vary vt, Q_L) and a first two-L sequential capture attempt. If Step 2 kills, Line A pauses at a clean, documented negative and the split inverts.
- **Line B (representation):** exactly the kimi-k3-endorsed **2D free-medium + typed-locks CPU toy**, standalone code (no `scp_sim` changes, no kernel auth), timeboxed to ~2–3 days of effort. Deliverable: does a "free substrate + gauge defects" state reproduce even one v78 relation (e.g., a piece of R1–R10) that the density-only state cannot? If not, the thesis stays a design note and stops consuming attention.
- The lines don't contend for resources (GPU vs CPU), so the split is about analysis/attention, not hardware. Do not let Line B's conceptual appeal displace the concrete Stage-3 measurements — the carbon goal runs through hydrogen, not through a rewrite.

## F. Confidence

**3.5 / 5.**

High confidence in the ordering (windowing check → force law → orbit): it directly targets the two measurements Stage 4 provably needs, requires no kernel auth, and every branch has a decisive outcome. The main deductions: (1) S_force and S_Lhold rest on pruned/ambiguous data, so my premise that "minimal C+L is promising" could dissolve at Step 1; (2) I cannot verify from the brief how far D can be pushed before box-image/absorbing-BC systematics (a known v70 failure mode) contaminate the force fit — grid sizing in Step 2 may need one iteration; (3) the representation-line timebox is a judgment call on a design I've only seen summarized. None of these change the recommendation's structure, only its first-attempt parameters.
