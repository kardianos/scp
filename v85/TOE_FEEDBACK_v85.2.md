# TOE v85.2 FEEDBACK — adversarial + constructive review

**Date:** 2026-07-21
**Agent:** Claude Fable 5 (same session as `v85/ANALYSIS.md`; that document's data and scratchpad reused)
**Target:** `v85/TOE_v85.2.md` (Steps 0–11, pillars, locks, experiments)
**Constraints honored:** kernel sources read-only, no lattice runs; new scratch script `toe_scales.py` in the session scratchpad; every number below traces to a shown command/output or to ANALYSIS §2.

---

## 0. Verdict summary

| Item | Verdict |
|---|---|
| Fidelity to ANALYSIS corrections | **High.** All eight major corrections are incorporated (ε>0, H5-weak, integral D8, redirect-not-block, 3-regime detuning, no interior Q\*, kinetic park, −Q_N neutrality). Residuals: one misattributed pillar number (99.75%), one wrong scale number (166 → **157.4**), two over-tags, dropped X10 error budget. |
| Step 9.1 (E_rad = ΔQ·Δω) | **Stronger than it looks — with new support and new limits.** I verify E = ωQ *exactly* for bound KG-Coulomb eigenmodes (canonical-Hamiltonian proof, §2.2), and classical radiation reaction *does* pin endpoints asymptotically. But at finite carrier charge the levels self-shift by **~29% of transition energies** (computed §3 step E): the Planck *ratio* survives; universal sharp lines do not. Interrupted transitions legally trap fractions — an unregistered anti-QM discriminator. |
| Step 9.2 (carriers required) | Logic valid as a *correspondence constraint* (nature-input, not fabric measurement) — mis-tagged [D]. Carrier *existence* has no kernel mechanism (T2 absent, admitted) and no designed rung. |
| Step 9.3 (ω³ bath) | **Fails on the lattice as stated.** Preferred frame + cutoff void the uniqueness argument; at any candidate amplitude the bath outweighs matter **30–8000×** (§3 step F); absorbing BC cannot sustain it; classical evolution thermalizes it (the known SED disease). |
| Step 9.4 (bath-equilibrium pins q₀) | **Incoherent as stated:** the only scale-free (massless) sector is chargeless — absorption/emission balance can pin E at fixed Q but can never move Q; VK stability supplies no restoring force toward any special Q. Constructive replacement offered (§2.5): one-way **bath erosion** with survivor floor at Q_thr(E/Q = m) — coherent, and config-testable as a *transient*-agitation assay. |
| Scale arithmetic | Doc's "≈166" corrects to **157.43** (g=0.05); thresholds now exact per branch (144.80 / 146.65 / 157.43 / 232.42 at g=0/0.02/0.05/0.10). Pairwise "derive g" fixed points exist: **g\*=0.0674 (q₀=93.2)** for the fold pairing, **g\*=0.0412 (q₀=152.5)** for the threshold pairing; exact identity **α = g/2** if q₀=2π/g; **no triple coincidence** (Q_thr/Q_min = 1.67–1.76 always); **no g reproduces α = 1/137** (would need q₀=431 vs Q_thr≈146 there). New bound found: the absolutely-stable window closes at **g ≈ 0.10–0.11** (Q_thr ≥ Q_max). |
| In-box testability of Step 9 | **Severe, unregistered:** shell-transition "photons" have λ_rad = 5,000–11,300 — **45–60× larger than any feasible box** (§3 step G). BATH-1 and any shell-scale emission-quantization test are not runnable as designed. Feasible salvage: core-clock-difference lines (Δω ≈ 0.05–0.09, λ ≈ 70–125). |
| New structural gaps | Carrier count for Z-carbon neutrality is **not 6** for any candidate q₀ (5.17–7.25; exactly 6 needs q₀=108.3 — closest scale is the *light nucleon* 114.13); carrier mass ≈ nucleon mass (**no m_e/m_p hierarchy**: 1/6.8 vs nature's 1/22,000 for the carbon system). The v75 light-sector problem re-enters through the carrier door. |
| Composite universality (Step 8) | Correctly stated — but the roadmap consequence is not drawn: **Stage 2A's single-clock droplets are disqualified as nuclei-for-atoms; 2B becomes mandatory** (currently "optional" in CLAUDE.md; v83 multi-center attempts were FAIL/NULL). |

---

## 1. Fidelity check (ANALYSIS → TOE)

**Correctly incorporated (verified against ANALYSIS §§1–2):** ε = +0.9–4.3% and E=Qω as thin-wall asymptote (Step 3); H5-weak with the continuum-branch refutation quoted (Step 2); M ≡ ∫T₀₀/c² with D5′ (Step 3 + locks); integral-form D8 with the 33% product-form failure (Step 3); three-U(1) verification and redirect-not-block (Step 4); the 3-regime detuning law with rectified term and lock bound (Step 5); F(Q) frozen form, F(650)=0.980, no interior iron peak, fullness-can't-stop-annihilation (Step 5); kinetic park + episodic shed + F6 AMBIGUOUS (Step 6); shell spectra quoted digit-exact vs ANALYSIS §2D, α ceiling 0.17–0.21 (Step 7); −Q_N neutrality (Step 7); λ = 2πQ(1+ε)/p (Step 8); honest D6 (Step 8); pillar precisions match their sources.

**Residual misstatements to fix (quoted):**

1. Step 2 / pillars: *"[M — v73 ledger: 99.75% cycle closure measured]"* and *"Process-form ledger closure | 99.75%"*. **Misattributed.** 99.75% is the **pass-through layment completeness of a moving ball at fixed probes** (v73 §3 E2). The stationary cycle-closure measurement is E3's flux profile (|P(r)| ≲ 2×10⁻³, alternating sign). And ANALYSIS §1.4 flags the 0.25% residue as a predicted **ARTIFACT(dx)** pending X2 — it should not sit unqualified in the pillar table.
2. Step 9.4: *"evaporation threshold Q(E/Q = m) ≈ 166 [D, interpolated]"*. **Wrong number** (+5.4%): exact values computed in §3 step A — **157.43** at g=0.05 (144.80 ungauged, 146.65 at g=0.02, 232.42 at g=0.10). The doc's self-flag "exact value is free analysis" is honored: here is the analysis.
3. Step 2: *"existence window … (1.40624, 1.4975) at g=0.05"*. 1.4975 is a **scan endpoint**, not a physical edge: the existence window top is m = 1.5 (CLAUDE.md concurs); the VK-**stable** segment ends at ω ≈ 1.4848 (Q_min fold). Two different physical points, neither at 1.4975.
4. Step 7: the X10 predictions are quoted without ANALYSIS §2D's bound caveats (core-overlap shift O(0.4–0.8%), l-splitting unresolvable, box-mode discretization). A prediction stripped of its error budget will false-fail or false-pass.
5. Step 8: *"detection is a whole-carrier transaction (all-or-nothing by H5)"*. H5 gives closure-or-death for *processes*; all-or-nothing *transfer* is T2 — which the doc's own lock table declares absent. The sentence silently substitutes the [C] carrier hypothesis for a measured axiom. Reword as "by carrier integrity (Step 7, [C])".

---

## 2. New-content critique (Step 7 carriers, Step 9 chain, BATH-1)

### 2.1 Steelman first — what genuinely works

The backward chain is a real tightening of the program. Its skeleton: (light is quantized) ⇒ (matter transfers charge in universal units) ⇒ (shell occupants are discrete carriers) ⇒ (something pins the carrier scale on a continuous branch) ⇒ (H1's moving vacuum is the only available pinning agent). Each implication is honest about its direction (working *backward* from nature), the carrier picture **repairs the ħ-universality defect** I raised in ANALYSIS §1.2 (per-object ℏ_eff=Q becomes universal iff all matter is q₀-closures and composites retain constituent clocks — one coherent architecture with Step 8), and the un-carried alternative (continuous cloud ⇒ classical light) is correctly identified and rejected. Registering LIGHT-STAT and Q-PIN as open locks is the right epistemic move.

### 2.2 Step 9.1 — E_rad = ΔQ·ω_rad: verified further than claimed, then bounded

**New support (constructive).** The doc tags "linear modes obey E = ωQ exactly" as [D] without proof for *bound* modes. I close that gap: for a charged KG eigenmode φ(x)e^{−iωt} in a static external Coulomb potential, the **canonical Hamiltonian** (which books the −qA₀ρ interaction term correctly) gives, using the eigenvalue relation ∫|∇φ|² = ∫[(ω+α/r)²−m²]|φ|²:

  H = ∫2(ω+α/r)[(ω+α/r) − α/r]|φ|² = ω·∫2(ω+α/r)|φ|² = **ω·Q exactly**.

(The naive energy ∫|D_tψ|²+|∇ψ|²+m²|ψ|² over-counts by ∫2(α/r)(ω+α/r)|φ|² and would *fail* E=ωQ — the identity holds for the conserved H, not the field-energy integral. Worth stating in the doc; it is the actual theorem Step 9.1 leans on.)

**Endpoint pinning — partly classical, no carriers needed.** In a fixed linear problem, the only non-radiating configurations are single-eigenmode (or degenerate-mixture) occupancies; any cross-term beats and radiates until exhausted. So radiation reaction *classically* drives complete transfer, and E_rad → ΔQ·Δω asymptotically by conservation. The doc's "endpoints pinned by carrier integrity, whatever the intermediate chirp" is thus *stronger than needed* for the isolated ideal case — carriers are needed for three other things: (a) **universality of ΔQ across events** (without carriers each event has its own continuous ΔQ); (b) **interruption protection** — a collision or boundary event mid-decay legally traps a *fraction* of ΔQ·Δω on this kernel (nothing forbids it): the fabric predicts fractional-energy events under interruption, real QM forbids them — an unregistered, in-principle-measurable discriminator; (c) absorber-side thresholds.

**New limit (computed, §3 step E).** At finite carrier charge the linearization that makes E=ωQ exact breaks: the carrier's self-Coulomb energy difference between modes is ~(gq₀)²/8π·(1/⟨r₁ₛ⟩−1/⟨r₂ₛ⟩) = 0.0248 vs transition energy q₀Δω = 0.0849 — **29%**. Consequences: level positions dress at O(30%) (the §2D/X10 tables are the q₀→0 limit); the emitted line chirps between initial- and final-dressed beat frequencies. The Planck **ratio** E_rad/ω̄_rad = ΔQ survives (the source always radiates at its actual instantaneous beat), but **universal sharp spectral lines do not** — atomic-spectroscopy correspondence degrades to ~30% at q₀ ≈ 90–157 unless screening cancellations rescue it (not computed anywhere; should be lock-listed).

### 2.3 Step 9.2 — carriers: valid constraint, phantom mechanism

The implication "continuous cloud ⇒ continuous ΔQ ⇒ classical light" is sound. But: (i) its premise ("light must be quantized") is a **correspondence requirement imported from nature**, not a fabric measurement or derivation — [D] is the wrong tag (§4); (ii) the kernel provides **no mechanism** that holds a bound-mode occupancy at fixed q₀: a carrier *in* a shell mode is not a branch soliton (E/Q = ω_n ≈ 1.4994 < 1.5283 = branch value at Q=90 — it is *off-branch*), so VK stability does not protect it, and the phase-blind linear sector lets any fraction of its charge move to another mode continuously. Carrier integrity in a *bound mode* is exactly T2 — absent. The doc's own lock table says so ("carrier-level integrity is the working substitute") — but a substitute needs a mechanism or at least a rung, and none is designed: X8 tests flavor conservation, not integrity; F6-cluster tests *free* lump quantization. **Missing rung:** seed a shell mode with 1.5 q₀ of charge; if carrier physics is real, relaxation must expel 0.5 q₀ rather than sit at 1.5 — currently nothing in the theory forbids 1.5.

### 2.4 Step 9.3 — the ω³ bath on a lattice: four independent failures

1. **The uniqueness argument does not survive the lattice.** ρ(ω) ∝ ω³ is the unique *boost-invariant* spectrum of a massless field in continuum 3+1D. This lattice has a preferred frame and measured Lorentz violation at object scales (v70 group-velocity lag 1–5%), and a cutoff ω_c = π/dx ≈ 8–11 where the ω³ weight concentrates its energy. "Boost invariance forces" is false here; at best "boost invariance of the continuum limit motivates" — and the bath's dominant content sits precisely where the continuum limit is worst.
2. **Energy budget is catastrophic** (§3 step F): ρ_bath = ħ_bath·ω_c⁴/8π². At ħ_bath = q₀ = 90: **8,300×** the ball's mean energy density on the pathfinder grid (2,800× on the v74 grid). Even ħ_bath = 1 gives 31–92×. Any bath at a candidate-scale amplitude does not perturb the matter sector — it erases it.
3. **Absorbing BC cannot sustain it.** The sponge damps all velocities in the boundary shell (scp_sim.c:1440–1470); an injected bath decays on the box-crossing time (~L/c ≈ 50–100 t.u.). "Sustained bath" requires continuous boundary injection — machinery the kernel does not have; adding it is a kernel change requiring explicit authorization. BATH-1 as listed ("park balls … in a sustained ω³-shaped bath") is **not config-expressible**.
4. **Classical evolution destroys the spectrum.** A classical nonlinear medium in a structured bath thermalizes toward equipartition (Rayleigh–Jeans) — the known instability of stochastic electrodynamics, which the doc's LIGHT-STAT lock acknowledges for *statistics* but not for the bath's own **stationarity**, which Q-PIN needs. On this kernel the ball itself is the nonlinearity that scrambles the spectrum.

### 2.5 Step 9.4 — bath-equilibrium selection: incoherent as stated; a coherent replacement exists

**Why the stated mechanism cannot work.** "Absorption balances emission" must balance *charge*, not energy, to select Q. The only sector that could plausibly carry a scale-free bath is the massless gauge field — which is **chargeless**: absorbing/emitting it changes E at fixed Q. The charged sectors (Φ, and Θ in the gauged kernel — both link-transported, scp_sim.c:1151–1154) are **gapped** (m=1.5, m_θ=1.6): their bath is thermal-like, not ω³, and dies on the sponge. Meanwhile VK stability holds at *every* branch Q — there is no restoring force toward any special charge. Net: a neutral bath pins nothing; a charged bath isn't the postulated object. Also note the sign structure: bath heating *enables* shedding (charged, above-gap emission — the v74 evaporation channel), while nothing in a neutral bath ever *adds* charge — the flow in Q is **one-way down**.

**Constructive replacement (offered, [C]):** take the one-way flow seriously. Under sustained agitation, complete evaporation is energetically open only below the threshold **Q_thr(E/Q = m)** and closed above it — Q_thr is the unique erosion-respecting boundary on the branch. "Bath erosion" then predicts a **survivor floor**: agitated populations converge onto Q ≳ Q_thr(g) (= 157.4 at g=0.05), with sub-threshold balls slowly disintegrating. This mechanism (i) uses only measured energetics, (ii) needs no detailed balance, no stable ω³ spectrum, and no sustained bath — a **transient**-bath survival assay (inject bath as initial condition; repeat; track the Q census) is config-expressible today with `gen_qball_bath`-class seeding, absorbing BC and all; (iii) hands Q-PIN a definite candidate: **q₀ = Q_thr(g)**, and with it the fixed point of §3 step C: g\* = 0.0412, q₀ = 152.5. It remains [C] — but it is a *coherent* [C]. Recast BATH-1 accordingly (§7, R6).

### 2.6 Step 7 carriers — the two unregistered structural gaps

1. **Carrier count ≠ 6.** Neutralizing the Z-carbon core (Q_N = 650) takes 650/q₀ carriers: **7.25 / 5.17 / 4.13** for q₀ = Q_min / 2π/g / Q_thr (§3 step D). Exactly 6 requires q₀ = **108.3** — none of the candidate scales; the nearest object in the theory is the **light nucleon** Q_N = 114.13 (giving 5.70). The "Z=6" chemistry mapping does not currently survive carrier quantization with any proposed q₀.
2. **No mass hierarchy.** A free q₀=90 carrier weighs E = 138 vs the carbon core's 934 — ratio **1/6.8** (nature: m_e/m(¹²C) ≈ 1/22,000); per nucleon, carriers weigh ≈ 0.8 nucleon masses. The "electron" of this theory is nucleon-heavy. This resurrects, inside the carrier picture, exactly the light-sector problem that motivated v75 multi-fabric — and nothing in v85.2 addresses it. It belongs in the locks table.

---

## 3. Scale arithmetic (commands + outputs)

All from `toe_scales.py` (scratchpad; inputs: `branch_g{000,002,005}.tsv` regenerated in ANALYSIS §2A from the verbatim v69 shooter, plus `v69/theory/gscan.tsv` for g=0.10).

**Step A — exact evaporation thresholds** (E/Q crossing m = 1.5 on the VK-stable side, linear interp between adjacent solutions):
```
g=0.00: Q_min= 86.59  Q_thr=144.80 (ω=1.43795)  Q_max= 352.3*  2π/g=  ∞
g=0.02: Q_min= 87.06  Q_thr=146.65 (ω=1.43875)  Q_max=5199.2   2π/g=314.16
g=0.05: Q_min= 89.69  Q_thr=157.43 (ω=1.44317)  Q_max= 923.9   2π/g=125.66
g=0.10: Q_min=102.82  Q_thr=232.42 (17-row gscan block; sparse)  Q_max=260.3  2π/g=62.83
```
(*g=0 table regenerated only to Q=352; the ungauged branch continues to Q=1.44×10⁵ per gscan. The g=0 threshold 144.80 confirms v72's "continuum η=0 branch puts Q\* ≈ 145".)

**Step B — the doc's "≈166":** computed Q_thr(g=0.05) = 157.43 → the doc's number is **+5.4% off** and matches no branch (g=0.02 gives 146.65). Correct it.

**Step C — "derive g" fixed points** (quadratic fits through the four g values; scipy brentq):
```
Q_min(g) fit [1955.25, −34.75, 86.71];  2π/g = Q_min(g) at g* = 0.0674 → q₀ = 93.2
Q_thr(g) fit [11994.5, −338.8, 146.04]; 2π/g = Q_thr(g) at g* = 0.0412 → q₀ = 152.5
Identity: q₀ = 2π/g ⟹ α = g²q₀/4π = g/2 exactly.  At g*=0.0412: α = 0.0206.
Nature's α = 1/137.04 would need g = 2/137 = 0.0146 → q₀ = 2π/g = 431 — but
Q_thr(0.0146) ≈ 146 ≠ 431: no g reproduces the real fine-structure constant.
Triple coincidence impossible: Q_thr/Q_min = 1.67 (g=0) … 1.76 (g=0.05) — never 1.
```
Verdict on the speculation: **realizable but under-determined** — demanding flux-quantum/charge-quantum consistency does pin g, but the answer depends on *which* charge scale is nominated (0.0674 vs 0.0412), the two pairings disagree, all three scales never merge, and neither fixed point lands near nature's α. The clean by-product worth keeping: **α = g/2** under any flux=charge identification. Caveats: the g=0.10 row rests on a 17-row gscan block (fold under-resolved); and comparing Noether charge to a flux quantum needs a Dirac-type commensurability argument — on this kernel matter is fixed charge-1 under the compact group, so "charge quantum = flux quantum" is a category leap to flag, not a lattice identity.

**New structural bound (from the same fits):** Q_thr(g) rises steeply while Q_max(g) ≈ (2.1–2.6)/g² falls; at g = 0.10 the absolutely-stable segment [Q_thr, Q_max] = [232, 260] is an 11% sliver, and the fits close it entirely at **g ≈ 0.10–0.11**. Above that coupling the theory has *no* evaporation-proof matter. (Fit-extrapolated; one dense shooter block at g=0.10–0.12 would confirm — analysis-only, cheap.)

**Step D — carrier counts** (Z-carbon, Q_N=650): 650/q₀ = 7.25 (Q_min), 5.17 (2π/g), 4.13 (Q_thr), 5.70 (light nucleon); exactly 6 ⟺ q₀ = 108.3. **Step E — self-energy chirp:** Δ(self)/(q₀Δω) = 0.0248/0.0849 = **29%** (pathfinder 1s→2s, q₀=90; point-Coulomb estimate, ±factor 2). **Step F — bath density:** ρ_bath = ħω_c⁴/8π² = 16,595 (ħ=90) / 184 (ħ=1) vs ball mean ≈ 2 → **8,300× / 92×** (pathfinder grid; 2,807×/31× on the v74 grid). **Step G — box cavity:** λ_rad = 2π/ε₁ = **4,998** (pathfinder) and **11,313** (carbon) vs boxes of 110 and 190 → **45–60× too small**; feasible radiation window is ω_rad ≳ 2π/2L ≈ 0.03–0.06, i.e. core-clock-difference lines (Δω ≈ 0.05–0.09, λ ≈ 70–125), not shell lines.

---

## 4. Epistemic-tag audit

| Location | Tag | Should be | Why |
|---|---|---|---|
| Step 2 "99.75% cycle closure" | [M] | [M, re-described] | Right number, wrong measurement (§1-1); residue is a flagged X2 target |
| Step 3 "linear eigenmodes E=ωQ exactly" | [D] | [D] ✓ (keep, add proof + limits) | Now actually proven for bound modes (§2.2); holds in the linearization/external-field limit; breaks at O(30%) for finite q₀ |
| Step 5 hard-core row | [M] | [D/G] | κ-knee suppression ×12 is potential-form arithmetic at measured s; capacity fold is [M]; but contact **repulsion between saturated objects has never been run** (X6b unrun) |
| Step 7 shell geometry | [D, G] | ✓ | Appropriate; restore the error budget |
| Step 7 carriers | [C→G] | **[C]** | "G = designed, numbers exist, not yet run" — no carrier-integrity rung exists; the only numbers are modes (which are the *continuous* piece). Design the 1.5 q₀ rung (§2.3) to earn G |
| Step 9.1 | [D] | [D conditional on carriers + completion] | Per-mode identity is [D]; "Planck's relation with ħ=ΔQ" inherits carrier-[C]; state the conditionality |
| Step 9.2 | [D] | **[C-corr]** (correspondence constraint) | Valid inference *from an empirical target*, not from fabric measurements; the tag system needs this category or [C] |
| Step 9.3, 9.4 | [C] | [C] ✓ | Correctly tagged; but 9.3's inner sentence "boost invariance **forces**" is stated as [D]-grade — weaken on-lattice (§2.4-1) |
| Step 9.4 "≈166 [D, interpolated]" | [D] | [M-corrected: 157.43] | Number was wrong; now measured from the branch |
| Pillars "ledger closure 99.75%" | pillar | re-describe | Same as row 1 |
| Step 8 "all-or-nothing by H5" | — | "by carrier integrity [C]" | H5 does not supply transfer quantization; T2 does, and it is absent |

---

## 5. QFT-backward

**Maps well.**
- *Emission bookkeeping:* E_rad = ΔQ·Δω ↔ E = ħω with ħ = q₀, including the ratio's robustness to level dressing (§2.2). Absorption/stimulated processes = classical driven response — natural.
- *Spontaneous emission:* classical radiation reaction reproduces the A-coefficient structure (ω³|d|²); the bath's added role (no-delay photoelectric, ground-state stabilization) is exactly the SED correspondence — the one regime where SED historically works.
- *AB effect, flux phases:* native to the compact-U(1) sector; already in the canon table.
- *ħ-universality:* the carrier picture is the first proposal in the program that could make ℏ_eff universal across matter, light, and angular momentum — **iff** all matter is q₀-composites retaining constituent clocks (Step 8's 2B architecture). Internally coherent.

**Conflicts / gaps.**
1. **Field-level mode quantization is absent.** In QED the radiation field itself is quantized; here light quanta are parasitic on matter transitions. Free-field consequences follow: classical blackbody equilibrium is Rayleigh–Jeans (UV catastrophe) unless the ω³ bath is *imposed and stable* — and §2.4-4 argues it is not stable on this kernel. Planck's spectrum — the doc's own Step-9 motivation — is not yet derivable.
2. **Photon statistics:** antibunching/sub-Poissonian light is impossible for a positive-P classical field + bath. LIGHT-STAT registers this; it should also register that *this* is where the "carriers make light quantized" claim finds its sharpest empirical ceiling: quantized *exchange* ≠ quantized *field*.
3. **Vacuum energy:** QED normal-orders the zero-point energy away; on a lattice the bath is real energy, 30–8,000× matter density (§3 F). Under Step 10's monism this should gravitate catastrophically — the cosmological-constant problem at full strength. **Constructive out (worth writing into Step 10):** monist gravity attributes curvature to *locked structure* (free-capacity deficit), not to free motion — if only closures gravitate, the bath is weightless by construction. That is a falsifiable structural stance unique to this ontology, and it is currently unstated.
4. **Universality of angular momentum:** QM quantizes L in the same ħ for everything; here L_z = nQ per object — universal only under the same 2B/carrier architecture (conflict resolved *only* if point 4 of §6 is solved).
5. **Mass hierarchy:** QED's electron is 1/1836 of a nucleon; the fabric's carrier is ≈0.8 nucleons (§2.6-2). No mechanism in sight on one fabric; this is the strongest single quantitative disagreement with nature in the whole edifice.
6. **In-sim spontaneous emission** is Purcell-modified: the box is far below cavity cutoff at shell-line frequencies (§3 G) — an instrument conflict, not a theory conflict, but it silently invalidates naive in-box tests of the QED mapping.

---

## 6. Top-5 weakest load-bearing points (ranked)

1. **Q-PIN via bath equilibrium (Step 9.4).** Chargeless massless bath cannot select charge; VK gives no restoring force; ω³ is neither unique nor stable on-lattice; absorbing BC forbids sustained baths; candidate amplitudes swamp matter 10²–10⁴×. Everything called "one constant" hangs on this. (Salvage path exists: erosion floor at Q_thr — §2.5.)
2. **Carrier existence without a mechanism (Steps 7/9.2).** T2 is admitted absent; bound-mode occupancies are off-branch and unprotected; no rung tests integrity. Quantized light, atom occupancy, and the ħ-unification all route through this single unsupported node.
3. **In-box untestability of the quantum-light sector.** Shell photons 45–60× larger than any feasible box; BATH-1 unrunnable as designed; the program's decisive Step-9 experiments cannot currently discriminate quantum from classical light at the scales where the theory makes its claims.
4. **2B-mandatory + missing hierarchy + carrier count.** Composite universality disqualifies Stage-2A droplets as nuclei for atoms (2B has only FAIL/NULL history); carriers weigh ≈ nucleons (no m_e/m_p analog); no candidate q₀ gives 6 carriers for Z-carbon (need 108.3; candidates 89.7/125.7/157.4). Three independent misses converging on the same architectural node.
5. **Finite-q₀ line dressing (29%) + interruption fractions.** The Planck ratio survives, but universal sharp spectra don't, and interrupted transitions legally trap fractional quanta — both are concrete anti-QM phenomenology the doc neither registers nor turns into discriminating predictions.

(Runner-up: the scale numerology's under-determination — two incompatible fixed points, no triple coincidence, no route to α = 1/137 — §3 C.)

---

## 7. Recommended revisions (concrete edits)

- **R1 (Step 2 + pillars):** re-describe 99.75% as "pass-through layment completeness (moving ball, v73 E2)"; add "0.25% residue = predicted ARTIFACT(dx), X2 target."
- **R2 (Step 9.4):** replace "≈166" with the measured table: Q_thr = 144.80 / 146.65 / 157.43 / 232.42 at g = 0 / 0.02 / 0.05 / 0.10 [M].
- **R3 (Step 9.1):** add the canonical-H proof of E=ωQ for bound modes; state Planck-relation conditionality (carriers + completion); register the two new discriminators (fractional-energy interruption events; 29% dressed-line chirp at finite q₀) as predictions-or-problems with a lock.
- **R4 (Step 9.3):** demote "boost invariance forces ω³" to a continuum-limit motivation; add the lattice energy-budget numbers (30–8,300×) and the BC-sustainment impossibility as stated constraints on any bath claim.
- **R5 (Step 9.4):** either adopt the erosion-floor mechanism (q₀ = Q_thr(g), one-way flow, survivor census) as the Q-PIN default [C], or mark Q-PIN "no coherent mechanism on this kernel" — the current detailed-balance wording is not repairable as written.
- **R6 (BATH-1):** redesign as a **transient-agitation survival assay**: seed ball + band-limited bath (initial condition only; amplitude ≤ O(ball); avoid the ω³ label), repeat over Q ∈ {90…250} spanning Q_thr = 157, score survival/disintegration census vs Q. Config-expressible now; tests the erosion floor, which is the testable core of Q-PIN.
- **R7 (Step 7):** restore the X10 error budget (core-overlap %, l-splitting unresolvable, box-mode spacing); retag carriers [C]; add the carrier-integrity rung (seed 1.5 q₀ in a mode; integer expulsion = pass) — that rung, once designed with numbers, earns the G.
- **R8 (Step 8):** draw the roadmap consequence explicitly: 2B is now load-bearing for atoms (Stage-2A droplets disqualified as atomic nuclei); flag for the user's CLAUDE.md stage-table decision.
- **R9 (locks):** add HIER (carrier/nucleon mass ratio ≈ 0.8 vs nature's 1/1836 — no mechanism) and COUNT (Z=6 needs q₀=108.3; candidates miss it; nearest object is the light nucleon 114.13). Add the α = g/2 identity and the g ≈ 0.10–0.11 no-stable-matter bound to the measured-structure section (fit-extrapolated; confirmable by one cheap shooter block).
- **R10 (Step 10):** add the constructive monist answer to the vacuum-energy conflict — only locked structure gravitates, free bath motion does not — as an explicit [C] with the note that it is this ontology's native resolution of the cosmological-constant problem, and a falsifier if Step-10 work ever shows bath energy curving the fabric.

---

*Method note: new numerics in `toe_scales.py` (session scratchpad), consuming the ANALYSIS §2A branch tables and `v69/theory/gscan.tsv`; the KG numbers reuse `kg_coulomb2.py` outputs quoted in ANALYSIS §2D. No repo file other than this one was created; kernel sources untouched. The g=0.10 rows inherit gscan's 17-row sparseness — flagged wherever used.*
