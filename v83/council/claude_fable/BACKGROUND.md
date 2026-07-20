# Advisor background: SCP Stage 2B / anti-lock substructure / Stage 3 unblock

**To:** Claude Fable (independent design critic)  
**From:** Program orchestrator  
**Date:** 2026-07-20  
**Task:** Critique the hypothesis that Stage 2B (multi-center substructure with high lock : single anti-lock core) is required to unblock Stage 3. See also `MY_HYPOTHESIS.md` written *before* your response.  

**You are not asked to implement code.** Give structured opinion: agree/disagree, risks, first kill experiments, alternative framings.

---

## 0. Standing goal

Produce a **carbon atom structural analog** from fabric alone (gauged complex Cosserat + free gauge medium + typed locks as evolved) — not imported chemistry particles.

Carbon is a **scale stack**. Do not skip stages.

| Stage | Target | Status |
|-------|--------|--------|
| 0 | Carbon mapping (\(Z,A\), observables via \(Q,Q_a,\omega_a\), Gauss) | Done — `v74/CARBON_MAP.md` |
| 1 | Nucleon template (gauged Q-ball; η-soliton when η>0) | Largely done v66–v73 |
| 2A | Liquid-drop Z-carbon nucleus | Done primary — c6_light parks Q→650 |
| 2B | True multi-center carbon (retain \(A\) substructure) | **Named; not delivered** — was labeled “optional research” |
| 3 | Light opposite-charge stable objects; positronium analog first | **Blocking for atoms** |
| 4 | C-nucleus + 6 light opposite charges, Coulomb-bound via \(A\) | After 2+3 |
| 5 | Spontaneous production / abundance near carbon | After engineered seeds |

Success is **structural** (ledgers, branch, stability), not quantitative chemistry.

---

## 1. Ontology (current best)

### Heavy sector (nuclear)
- Complex Cosserat multiplet \(\Phi_a\), product potential \(V(s)\), \(s=\prod|\Phi_a|^2\).
- Optional gauge \(A\) (photon analog); charge conservation + discrete Gauss exact by construction when implemented correctly.
- Stable objects: **Q-balls / balls / baryons** — fused charge droplets on a measured branch.
- “Quarks” = component fields \(\Phi_a\): fractional \(Q/3\), no isolated existence, only triple-complete binds (v71 QUARK).
- Flavor = \((\omega_a,Q_a)\) partition; flavored baryons dynamically stable (fb1).
- Nuclei in practice: **liquid-drop fusion** (co-phase multi-ball → one droplet). He/Li/Z-carbon measured. Super-critical evaporates toward \(Q_{\max}\), not fission by default.

### Light sector (Stage 3 redesign — monist locks)
- Multi-fabric product path (C/Q/L multiplet L as electron) **failed orbit / product** (v80 O FAIL). Do not resume product multi-fab.
- **P1 locks:** discrete typed carriers on free U(1)/Yee medium (PIC: Boris, charge-conserving J, Gauss). Kernel: `scp_locks.h`, `n_locks=0` safe.
- Charge locks (`type=0`): \(q=\pm1\), \(E^\star\) sequestered budget → inertia.
- **Anti-locks (`type=1`):** \(q=0\), no ρ/J; bag attract charge locks (pairwise or grid co-field \(B(x)\)). Called “lock-Higgs” — **not** Cosserat multiplet Higgs.
- Cosserat multiplet Higgs for light free sector: **wrong tool** (vacuum Φ; decorative or harmful).

### Process-form / capacity (theory thread)
- Free substrate + sequestration; short-range capacity well as second scale (v82 A) for free pairs.
- Pair-overlap Gaussian is a **proxy** for depletion, not final monist co-field.

---

## 2. Stage 2A facts (v74)

- **c6_light:** 6× light nucleons → stable on-branch **single droplet** (Z-carbon map).  
- **c12_light:** 12× → merge → evaporate toward \(Q_{\max}=921\) (A=12 hot control).  
- Multi-center static bound state with \(A\) retained: **not available**.  
- Explicit non-goals of 2A: electrons, chemistry, kernel mods.

Implication: current “carbon nucleus” is a **featureless (or soft multipole) droplet** for Stage 4 purposes unless 2B exists.

---

## 3. Substructure and interlock (v71) — the reason parts exist

### 3.1 QUARK campaign
- Single-component lump: \(s\equiv0\), disperses (no binding).  
- Two components: same.  
- Three components within ~2 core radii: baryon self-assembles.  
- Flavored (unequal \(\omega_a\)): stationary + 3D-stable clocks.  
- No meson sector (structural difference from QCD).  
- All three quarks same sign of diagonal U(1) charge — no \(+2/3\) vs \(-1/3\) inside one baryon via gauge alone.

**Theory reason substructure exists:** force on \(\Phi_a \propto \prod_{b\neq a}|\Phi_b|^2\). Incomplete color → free KG dispersal. Parts are not optional decoration.

### 3.2 Flavored collision taxonomy (COLLIDE)
Per-component relative phase is the nuclear-force knob (g=0):

| phases | outcome |
|--------|---------|
| all aligned | merge |
| all anti | exterior bounce (cores never touch; node in tail) |
| mixed (1 anti, 2 aligned) | collide-and-re-emerge; flavor-differential charge transfer |

Nuclear force ≈ per-flavor \(\cos\Delta\phi_a\) channel sum weighted by tail range \(\mu_a=\sqrt{m^2-\omega_a^2}\).

### 3.3 Phase-interlock molecule campaign (im1, im2, mc1, mc2)

Design intent: static two-ball standoff — repulsion short, attraction long.

- New baryon: two-low / one-high \(\omega=1.38/1.38/1.42\).  
- **im1** (2 short anti, 1 long aligned): nonlinearly **repulsive** in overlap; no lock.  
- **im2** (**1 anti, 2 aligned**): collapses to **one composite**: aligned flavors fuse to central peak; frustrated flavor = **two lobes around a node** (renders cyan). Metastable: frustrated flavor burns off ~−33%/500 t.u. Lifetime hundreds of t.u.  
- **mc1/mc2** crash two im2 molecules: **total fusion**; interlock survives as **momentum** (thrust along frustration axis), not retained multi-center geometry.

**Gauged interlock** (Coulomb + channels): listed FUTURE.md #11; **never run**.

Moral written in lab notes: “phase frustration is never load-bearing at equilibrium — ejects or burns off.” User now proposes anti-lock may supply the missing **load-bearing pin**.

---

## 4. Stage 3 path history (compressed)

| Era | Attempt | Result |
|-----|---------|--------|
| v75–v80 multi-fab | Opposite light multiplet L | F can pass; **O FAIL**; \(E_{\mathrm{em}}\) issues with light humps |
| v81 P1 locks | Free medium + charge locks | Durable carriers; Coulomb from medium; Gauss fix critical |
| N=256 + anti-lock bag | Second scale for free ± pair | Expansion slowed; **no multi-rev park** |
| S1 grid bag from anti-locks | Co-field bag | Collapse pass-through; revs ≪ 1 |
| Cosserat Higgs for light | Multiplet cavity | Wrong for vacuum light |
| v82 A capacity well | \(F_{\mathrm{core}}\) repel short + EM attract | Well PASS; circular seed bugs fixed; **analytic multi-rev PASS**; live PIC ≥1 band PASS (fragile) |
| v82 B magnetic | Uniform \(B_z\) | FAIL as primary scale |
| v82 C rigid composite | Fixed internal \(D\), free COM | PASS\* structural durability (defines away free \(r\)) |

Program had labeled 2B “optional.” User correction: **2B is required to unblock Stage 3.**

---

## 5. Anti-lock implementation facts (kernel)

From `scp_locks.h` / config:

- `type=1`: anti-lock; skipped in ρ/J deposit.  
- Bag modes: (1) pairwise attract charge locks toward anti-locks for \(r<r_{\mathrm{bag}}\); (2) grid \(B(x)\) CIC deposit from anti-locks, \(F=-\kappa\nabla B\).  
- Soft core: short-range repel among charge locks.  
- `locks_medium_only=1`: skip multiplet force (needed large N CPU).

Never systematically tested: **many charge locks + exactly one anti-lock as a single composite unit** intended to park as substructure (nucleon-analog or carbon-center).

---

## 6. User directive (verbatim spirit)

- Stop treating 2B as optional.  
- Substructure has a reason it exists (see §3).  
- Anti-lock Higgs useful **in the substructure**.  
- Ratio: **high lock count with only a single anti-lock** in the substructure.  
- Look for such a construction; free-pair A is not the atom path center of mass.

---

## 7. Orchestrator hypothesis (summary — full text in MY_HYPOTHESIS.md)

**H1:** 2B retained multi-center / composite anatomy is on the critical path to Stage 3/4; droplet-only 2A is insufficient as the atom’s nuclear target.  

**H2:** Anti-lock is the interior bag pin with stoichiometry \(N_{\mathrm{lock}}\gg 1\), \(N_{\mathrm{anti}}=1\) per unit; Stage 3 light partner is exterior opposite charge, not a second anti-lock.  

**H3:** Order of work: unit (locks+1 anti) → multi-center carbon → exterior light capture → six lights.  

**H4:** Falsifiers listed in MY_HYPOTHESIS.md.

Confidence: H1 medium–high; H2 medium (motivated by im2 anatomy + failed free-pair bag, not yet built).

---

## 8. Open technical gaps

1. No parked multi-center carbon.  
2. No high-lock : 1-anti composite park in monist locks.  
3. Gauged Cosserat interlock untested.  
4. Live free-pair capacity orbit fragile (A3 self-force / SOR).  
5. How Cosserat nuclear core couples to lock bag (OP foreshadowed; not designed).  
6. 2D log potential vs 3D Coulomb for exterior Stage 3.  
7. Ledger: anti-lock \(E^\star\) meaning if \(q=0\).

---

## 9. Constraints for advice

- Do not recommend resuming multi-fab product as Stage 3 delivery.  
- Do not recommend Cosserat multiplet Higgs as free light electron.  
- Kernel edits for Stage 3 were authorized earlier; prefer clear minimal deltas.  
- Prefer kill experiments that are cheap and decisive.  
- Distinguish: (a) theory necessity of 2B, (b) necessity of anti-lock stoichiometry, (c) monist locks vs Cosserat fields as the 2B substrate.

---

## 10. Questions for you (answer all)

1. **H1:** Is 2B required to unblock Stage 3, or can droplet 2A + free-pair Stage 3 still deliver a carbon atom structural analog?  
2. **H2:** Is “high lock : single anti-lock core” the right substructure stoichiometry? Better ratios? Anti-lock vs capacity co-field vs gauged interlock as the pin?  
3. **im2:** Is the cyan metastable composite the right geometric target to stabilize with anti-lock, or a red herring?  
4. **Substrate:** Should 2B be built first in (A) monist locks only, (B) Cosserat flavored multi-ball + gauge, or (C) hybrid bag-lock on Cosserat core?  
5. **First two weeks:** ordered experiment list with PASS/FAIL kill criteria.  
6. **What to stop doing** immediately.  
7. **Biggest risk** that this whole reframing is cargo-cult (importing “nucleus has parts” without fabric necessity).  
8. **Stage 3 definition:** after 2B, what is the minimal positronium-analog acceptance test?

Write as a design recommendation memo. Be skeptical of the orchestrator. Prefer null results and kill gates.
