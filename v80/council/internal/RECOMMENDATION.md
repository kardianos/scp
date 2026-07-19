# Strategic recommendation — SCP after v80 overnight campaign

**Advisor role:** independent; did not author campaigns.  
**Date:** 2026-07-19  
**Sources:** shared brief, v79 RESULTS, v80 campaign RESULTS/GOAL, optional `v80/` design notes.  
**Rule:** no invented simulation numbers; no Z6+L6 re-park as primary without a new mechanism.

---

## A. Timeline reading

What the last few versions actually established vs still open:

1. **Stage 2A is banked.** Liquid-drop Z-carbon parks; free A=12 stays super-critical at g=0.05. Nuclear baseline is not the bottleneck for atoms.

2. **Stage 3 architecture exists on paper and in multi-fab seeds (C/Q/L), but is not “solved.”** Opposite-charge light sector is the atom blocker; the program has a working multi-fab product path *and* a representation-critique line, not a cold electron analog.

3. **Same-fabric E-lite force (v75 F8–F9) and multi-fab charge bookkeeping (v80 S1) both work at the ledger level:** opposite ≈ neutral net Q, same ≈ 2Q; Gauss stays on the 1e-13 floor. Gauge health is not in doubt for these runs.

4. **v79 closed the hand-placed multi-electron atom recipe.** Z6N6+L6 at T=800: core inventory tracks the L=0 control (identical Q_phi / Q_core / r_core end state); L **nulls E_em** and crushes flux. That is a neutralizing/radiation bath, not a parked multipole. Re-running the same park is not science.

5. **v80 correctly shifted product work from “carbon park” to “force / pair / orbit graph.”** 16/16 jobs OK, G = 0.62 (≥ 0.55 continue). Soft-kill of multi-fab atoms is **not** justified.

6. **Hydrogenoid-class seed is the real product win of the night.** S4 T=800 held E_em (0.71→0.75) and Q_phi — the opposite of v79’s fat-shell death. Minimal C+L is the only multi-fab atom-adjacent object that looks alive on long T.

7. **Still unproven (explicit non-claims):** measured Coulomb \(a_{\mathrm{rel}}\); multi-rev orbit/capture; cold multi-L atom; free-substrate representation as primary state. Force score is soft because SFAs were pruned and pair tracks were never computed.

8. **Orbit is partial and energy-fragile.** Low \(v_t\) mild; \(v_t=0.08\) hard E drift (−27%). Without tracks, “orbit” is a vibes score, not a capture class.

9. **Diagnostics trap:** S4 \(Q_{\mathrm{flux}}\) / \(Q_{\mathrm{core}}\) collapse while \(Q_\phi\) holds is flagged as COM/window, not charge death. Product decisions must not treat flux alone as annihilation.

10. **Representation thesis (v80 design) is coherent critique, not GPU-proven.** Multi-fab may be N copies of one Eulerian density type; external feedback says density-only state is the problem, grid is innocent — but that is design pressure, not a kill of the product G score.

---

## B. Long-term goal status

| Stage | Status after v79+v80 |
|-------|----------------------|
| 0–2A | Done / banked |
| 2B multi-center nucleus | Optional; do not touch now |
| **3 electron / opposite light** | **Open and blocking.** Product line has force/sign and minimal C+L E_em hold; **no** demonstrated cold bound positronium analog with multi-rev orbit + measured attraction law |
| **4 carbon atom** | **Not started.** v79 multi-L park failed; stacking L on Z6 without capture dynamics is dead |
| 5 spontaneous production | Far future |

**What would count as real progress toward a carbon atom (structural):**

1. **Measured opposite-charge attraction law** on multi-fab C–L (or Q–L): \(a_{\mathrm{rel}}(D)\) from kept SFAs + `mf_pair_track` / equivalent, same-sign control, Gauss floor held.
2. **At least one multi-rev (or clear capture) orbit** with track TSV — not energy-drift vibes alone — for a **hydrogenoid-class** seed (1 light opposite, not 6).
3. **Stationary-ish bound ledger:** E_em held over long T (already S4-class), plus COM-aware flux/core that does not fake charge death.
4. Only **after** 1–3: dynamical multi-L assembly (capture / sequential binding), **not** hand-placed L shells on a Z=6 octa. Carbon atom = Stage 4; skipping the hydrogenoid force/orbit graph is how you re-earn v79.

**Honest product today:** fabric nuclei + multi-fab mechanics graph with G=0.62. **Not** an atom. **Not** Stage 4.

---

## C. Recommended next step (primary)

### Primary sequence (tightly ordered; do not reorder)

**C1 → C2 → C3 gate.** One product arc: close the missing force/orbit measurements on hydrogenoid-class multi-fab, then decide.

#### C1. Force law with kept SFAs + pair tracks (do this first)

- **What:** Re-run a **minimal** force matrix only:
  - Multi-fab C–L (or campaign B1) at **D = 12, 16, 20** (same distances as S1).
  - Controls: opposite vs same-sign E-lite if still needed for sign check.
  - **Keep SFAs** on remote/local `/space/scp/v80_tracks/` (do not prune mid-queue).
  - Run `mf_pair_track` / `sfa_qball_track` → \(D(t)\), \(a_{\mathrm{rel}}(D)\), optional COM-corrected cores.
- **Why before alternatives:** S_force = 0.55 is the weakest high-weight score **because tracks were never produced**. Every atom claim (positronium, capture, carbon) presupposes opposite attraction of known sign and rough scale. Without \(a_{\mathrm{rel}}\), orbit scans and L-hold extensions are decoration.
- **Success:** Opposite \(a_{\mathrm{rel}}\) toward each other; same-sign repel or weaker; monotonic-ish falloff over the three D’s; Gauss floor; no hard-fail disk/prune.
- **Kill / retune:** Opposite does not attract under track analysis; or attraction is merge-only at all D with no soft force regime; or tracks impossible due to tool/seed failure (fix tooling before theory).
- **Effort:** Short–medium GPU (order of hours, not overnight carbon matrix). Analysis mostly CPU on diags+SFA.
- **Kernel auth:** **No.**

#### C2. Orbit with tracks on the gentlest hydrogenoid seed

- **What:** Only if C1 shows soft opposite attraction. Re-run low-\(v_t\) orbits (0.03 / 0.05 first; **skip or demote** 0.08 until low-v works) with **kept SFAs** + track TSVs; publish multi-rev count or explicit non-capture.
- **Why:** v80 S3 proved “jobs complete”; it did **not** prove orbit. High-\(v_t\) already radiates. Capture class is the Stage-3 success proxy (positronium analog), not Z=6 parking.
- **Success:** ≥1 clear multi-rev or stable capture window with E_em not collapsing and track radius bounded; Gauss floor.
- **Kill:** All soft \(v_t\) either merge head-on or fly apart with no intermediate; E_em death on minimal C+L when dynamics turn on (contradicts S4 rest-ish hold — investigate seed, not stack L).
- **Effort:** Short–overnight GPU depending on T and N.
- **Kernel auth:** **No.**

#### C3. Gate only after C1–C2

- If tracks show Coulomb-class opposite force **and** at least one orbit/capture class: **then** T-extend hydrogenoid / mild multi-L (2, not 6) **with dynamical assembly**, still no Z6 hand shell.
- If force fails hard: soft-kill product multi-fab atoms (G path dies on physics, not on pruned SFAs) and **shift primary** to representation toy (Section E).
- If force OK but orbit always radiates/merges: product path becomes “force-only phenomenology”; do **not** claim atoms; representation work becomes the only route to Stage 4 language.

**Why this over re-running Z6+L6, kernel rewrite, or pure design week:**  
G≥0.55 says product has legs; v79 says multi-L hand-park is dead; S4 says minimal C+L is the viable object; the **only** high-value missing measurement is tracks. Representation is legitimate but does not replace closing the force graph you already paid GPU for.

---

## D. Explicitly defer or kill

| Action | Verdict | Why |
|--------|---------|-----|
| **Re-run Z6+L6 hand-placed park as main experiment** | **Kill as primary** | v79 already null: L kills E_em; core = L=0 control. No new mechanism proposed → pure waste. |
| **Stage 4 carbon atom matrix** | **Defer** | Stage 3 hydrogenoid force+orbit not closed. |
| **Stage 2B multi-center nucleus** | **Defer** | Orthogonal to atom blocker; dilutes focus. |
| **Kernel rewrite / free-substrate into `scp_sim`** | **Defer** | Needs explicit user auth; kimi/v80 design still undefined nouns (substrate, lock, sequestration). Toy first. |
| **Overnight “more multi-fab parks” without tracks** | **Kill** | Repeats v80’s blind spot: scores without \(a_{\mathrm{rel}}\). |
| **Soft-kill multi-fab atoms now** | **Do not** | G=0.62; hard_fail false; S4 E_em hold. Soft-kill only if C1 force dies or G≪0.3. |
| **Treating S4 flux/core drop as charge annihilation** | **Kill as interpretation** | Brief: COM/window; \(Q_\phi\) held. Fix diagnostics or tracking, don’t rewrite theory on window bug. |
| **Volview/morph-only package as main next** | **Defer** | S_morph=0.25 is low weight; pretty pictures after tracks. Optional add-on on kept SFAs. |
| **η>0 / process-form reopening for atom** | **Defer** | Atom path is multi-fab charge + gauge; don’t reopen η-drain without a specific atom mechanism. |

---

## E. Product vs representation

**Split for the next 1–2 calendar weeks** (assuming one person or one primary agent):

| Line | Share | What | When |
|------|------:|------|------|
| **(A) Product multi-fab** | **~70%** week 1, **~50%** week 2 | C1 tracks+force, then C2 orbit tracks; write force/orbit RESULTS with numbers from tracks only | Immediately; uses existing seeds/generators |
| **(B) Representation substrate** | **~30%** week 1, **~50%** week 2 if C1 OK; **primary** if C1 force dies | **CPU-only** 2D free-medium + typed locks toy (or smallest dual-state export) per v80/Kimi: define state type, one sequestration demo, one `export_grid()` that could later accept product metrics | Daytime design/code; **no** scp_sim change; **no** carbon park as “evidence” for the thesis |

**Rules of quarantine (do not mix):**

- Product runs **must not** be scored as proof of free-substrate ontology.
- Representation toys **must not** be sold as carbon progress until they export a measurable force/orbit proxy.
- If product C1 succeeds: keep product as Stage-3 ladder; representation remains strategic insurance and theory honesty.
- If product C1 fails: flip to **~80% representation / ~20% autopsy** (one controlled re-seed, not Z6 park).

**Do not** redesign multi-fab inside the kernel “to fix v79.” That spends effort on the representation you already called an emulation. Fix product recipes (tracks, dynamical capture) or change primary state in a **toy**, not a silent multiplet tweak.

---

## F. Confidence

**Score: 4 / 5**

| Strength | Weakness |
|----------|----------|
| v79 null is clean and controlled (atom vs L=0) | No actual \(a_{\mathrm{rel}}\) numbers — force recommendation is structural, not data-driven beyond sign/E_em patterns |
| v80 G=0.62 + S4 E_em hold cleanly separates “minimal C+L alive” from “fat L dead” | Orbit multi-rev unknown; S4 flux/core ambiguity |
| Next step is the campaign’s own stated gap (kept SFA + tracks) | Representation effort split is judgment; 30/70 could flip if user prioritizes ontology over Stage-3 ladder |
| Constraints respected (no kernel, no Z6 re-park primary) | External Kimi critique may underweight dynamical-capture as alternative to representation failure |

**Why not 5:** The highest-weight product claim (Coulomb-class force) is still **unmeasured**. The recommendation is the correct experiment to close that gap, not a claim that the force already exists at the needed quality.

**Bottom line:** Ship **tracked force on hydrogenoid-class multi-fab** before anything that looks like an atom. Park carbon shells only after opposite light charges show attraction + orbit/capture. Representation work in parallel as CPU design — not as an excuse to ignore the missing tracks, and not as a premature `scp_sim` rewrite.
