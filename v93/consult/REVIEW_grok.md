I'll read the review request and every file it points to, then write a rigorous critical review with concrete file:line citations.Next I'll pull the design docs, findings, and the pass-U / arg(ψ) kernel regions the brief names.Now reading the kernel pass-U, Givens hop, arg(ψ) door, and pass-F regions in detail.Checking QUENCH logs and whether L1_FINDINGS is missing the arg-door results claimed in the brief.Re-running winding analysis properly and checking a few more code details that affect the conclusions.# v93 Critical Review — Theory, Code, Results, Conclusions

Reviewer stance: agree that unitary dense transport is the right primitive and that L1-B / L1-C are structurally solid. Disagree that L1-A and the arg(ψ) door are as confirmatory as the write-up claims. Several code/design mismatches and analysis flaws weaken the stronger conclusions.

---

## 0. Executive verdict

| Claim | Verdict | One line |
|-------|---------|----------|
| Pass U is norm-preserving Givens hop (cousin of F) | **Sound** | Algebra matches pass F; pairwise conservation is a theorem |
| In-pass precession ordering | **Sound for dense history; not identical to F’s sign** | Dense matches old `th2 += w2e·dt`; field uses opposite chirality |
| Conservation theorem (L1-B) | **Confirmed** | drift_rel = 0 at amp_nat 0/2.5/5 |
| Anti-ignition (L1-C) | **Confirmed and well-diagnosed** | Bath byte-identical; gate self-starves on random phase |
| L1-A “mechanism confirmed / thesis confirmed” | **Overclaimed** | Directional bias is real; L1-A bar is not met; free-packet centroid motion is not yet “translation of a dense object” |
| Narrow resonance as dispersion | **Plausible, under-tested** | Consistent with tight-binding + on-site ω; measurement and binding confounds not ruled out |
| arg(ψ) door energy-conserving | **Ledger-conserving, not unitary-merge** | Phase-only composition; \|ψ\| forced by Em ledger |
| arg(ψ) implements §II.7 / IV.6 | **Does not** | Still a **cell-clock** write; not slot-borne |
| Early-time R2d “3× retention” | **Weak / overstated** | Partial support; nlive small; m-peak not locked to +2; pipeline bug in suite script |
| “Matter churn is the wall” | **Incomplete** | Field RA2 dies; imprint quality fails; reverse door loses phase |
| Packet needs binding | **Agree** | Free unitary packet must diffract; L1-A should be reframed or binding added |

---

## 1. Theory / code soundness

### 1.1 Pass U hop algebra — correct

The Givens update is identical in structure to pass F:

```1561:1565:v93/kernel/freecell.c
                double m1i = dm1_[i], m2i = dm2_[i], m1j = dm1_[j], m2j = dm2_[j];
                dm1_[i] = cc * m1i + ss * m2j;
                dm2_[i] = cc * m2i - ss * m1j;
                dm1_[j] = cc * m1j + ss * m2i;
                dm2_[j] = cc * m2j - ss * m1i;
```

vs field hops at `freecell.c:2086–2089`. Pairwise norm preservation is immediate (same cross-term cancellation as F). C and Go mirrors match (`step.go:488–491`). No sign error in the hop itself relative to the field template.

For small τ, energy flux is
\[
\Delta E_i \approx 2\tau\,\mathrm{Im}(\psi_i^*\psi_j) = 2\tau(m_{1i}m_{2j}-m_{2i}m_{1j}),
\]
so **J ∝ Im, not Re**. Code comments at `freecell.c:1522–1523` and L1_FINDINGS get this right. The design doc is inconsistent:

- `v93/README.md:176–179` — `2·Re(ψ_i*·w·ψ_j)`
- `v93/README.md:211` — `J_s = Im(...)` (correct)
- `v92/consult/SUBQUARK_synthesis.md:25` — Re again

**Fix:** normalize the design language to `J_s = 2 sin(τ) Im(ψ_i* ψ_j)` (exact first-order flux for this Givens form) and drop Re.

### 1.2 τ construction — intentional Face-2 choice, but incomplete vs design

```1477:1480:v93/kernel/freecell.c
            double gsym = sqrt(g_ij * g_ji);
            double tau = P.amp_nat * base * gsym;
            if (tau > 0.5) tau = 0.5;   /* integrator guard: deep-overlap cap */
            shau_[s] = tau;
```

with `base = kd * dt * geo * gpl * res` (`freecell.c:1391`).

What is good:
- `dt` is inside `base` → τ scales as a rate×dt, like pass F’s `field_J * w/√(si sj) * dt` (`freecell.c:2082`).
- Dropping head after Face 1 is correct: head≈0 at cap freezes the core (`L1_FINDINGS.md:35–43`).
- Mobility is implicit in ψ via √Em (energy that can move), not double-counted in τ — that is the right pass-F analogy.

What is fragile:
1. **τ clamp at 0.5 rad/step** is a hard nonlinearity. At amp_nat ≳ few on deep-overlap links, the map saturates → amp_nat is not a pure linear hopping strength. That alone can create “resonance windows” and direction flips without true continuum dispersion.
2. **gsym = √(g_ij g_ji)** is even in the closure phase. Directional transport must come entirely from Im(ψ_i*ψ_j) inside the hop, not from an asymmetric gate. That is fine in principle, but it means any residual directed bias in g_ij − g_ji (the old want) is discarded. You may be throwing away a useful odd channel (κ_reac already derived this for the additive path at `freecell.c:1395–1410` and is zero under V3a).
3. When `amp_nat>0`, **additive want + inflight (passes 3–5) are fully bypassed** (`freecell.c:1514–1515`). That is a bigger physics change than “replace debit with rotation”: flight load, rough conversion on links, and parcel-borne phase windows (`amp_tau` shadow deposits) all go dark for dense. The design said “transport primitive swap”; the implementation also **kills the flight sector** for matter.

### 1.3 In-pass precession — right ordering for dense; not “exactly” pass F

```1544:1550:v93/kernel/freecell.c
        for (int i = 0; i < NC; i++) {
            double ang = w2e[i] * dt;
            double cc = lut_cos(ang), ss = lut_sin(ang);
            double a = dm1_[i], b = dm2_[i];
            dm1_[i] = cc * a - ss * b;
            dm2_[i] = ss * a + cc * b;
        }
```

This is multiplication by **e^{+i w2e dt}**, matching the historical out-of-pass advance `th2 += w2e·dt` (`freecell.c:2143–2145` skipped when amp_nat>0).

Pass F does the opposite chirality:
```2057:2062:v93/kernel/freecell.c
        double ang = w1e[i] * dt;
        ...
        fa1[i] = cc * a1 + ss * a2;
        fa2[i] = -ss * a1 + cc * a2;
```
i.e. **e^{−i w1e dt}**.

So Face 3 correctly removes the dense Trotter split relative to the *dense* clock, but the claim “exactly mirroring pass F” (`L1_FINDINGS.md:107–109`, `freecell.c:1538`) is false on sign. Not necessarily a bug — field and dense may be conjugate sectors — but it should not be sold as “byte-identical ordering to F.” If you ever want a single complex-structure law across modes, this sign needs a design decision.

**Skip of pass-6 th2 when amp_nat>0 is correct** for avoiding double precession.

### 1.4 Sequential Givens product — shared Trotter with F

Hops are applied in roster order from the i-side only (`freecell.c:1551–1555`). Non-commuting rotations ⇒ order-dependent effective Hamiltonian. Same as F, so A/B-safe, but for a free packet at long T this is another wobble source beyond on-site/hop Trotter. Worth a sweep with randomized slot order (or checkerboard) to bound the artifact.

### 1.5 Empty-cell phase writeback

```1576:1581:v93/kernel/freecell.c
        for (int i = 0; i < NC; i++) {
            double a = dm1_[i], b = dm2_[i];
            Em[i] = a * a + b * b;
            double ph = atan2(b, a);
            if (ph < 0) ph += TWO_PI;
            th2[i] = ph;
        }
```

`atan2(0,0) → 0`: every empty cell’s clock is forced to 0 each step under amp_nat>0. That is a **global phase sink / artificial phase reference** the additive path does not have. On a bath that is “inert” this is invisible (τ≈0), but any cell that briefly holds then loses Em gets phase reset to 0 — a systematic phase-slip machine relevant to QUENCH. Prefer: if Em < ε, leave th2 unchanged (or leave undefined).

### 1.6 arg(ψ) door — energy-safe, not the designed retention fix

```2178:2193:v93/kernel/freecell.c
                    if (P.amp_door > 0) {
                        double aph = atan2(fa2[i], fa1[i]);
                        double Em_old = Em[i] - (d1 + dsp);
                        if (Em_old < 0) Em_old = 0;
                        double ro = sqrt(Em_old), rc = sqrt(d1);
                        double s1 = ro * lut_cos(th2[i]) + rc * lut_cos(aph);
                        double s2 = ro * lut_sin(th2[i]) + rc * lut_sin(aph);
                        th2[i] = fmod(atan2(s2, s1) + 8.0 * TWO_PI, TWO_PI);
                    }
```

**What is legitimate:**
- Em/Ee/Es ledger is unchanged: magnitude of ψ is still forced by Em; only th2 is rewritten. Total mode energy conservation is intact (L1-B still holds).
- aph is read after real scaling of fa (`fac = √((Ee−d1)/Ee)`), which preserves arg — OK.
- Coherent arg merge is strictly better than qp_phase’s partial pull (`mix * wrap(aph−th2)`).

**What is not legitimate as §II.7 / IV.6:**

1. **Still a cell-clock write.** Design: *“Phase is carried slot-borne (protected from delivery churn)”* (`README.md:238–239`, `IV.6` at `README.md:390–394`). Code: writes `th2[i]`. The forbidden m·th2 pattern is replaced by a better cell write, not by a slot-borne atom. IV.6 is not implemented.

2. **Not a unitary amplitude merge.** True coherent sum would set
   `Em_new = |√Em_old e^{iθ} + √d1 e^{i aph}|² = Em_old + d1 + 2√(Em_old d1) cos Δ`,
   which **breaks** the Em+=d1 ledger. You correctly refuse that (IV.4 / conservation). But then the cross-term is **discarded from energy and only used for phase**. Calling this “carries the interference cross-term = the current” (`freecell.c:253–255`) overstates it: the current is not deposited; only a phase estimate is.

3. **Space footprint `dsp` has no phase partner.** Merge uses `√d1` not `√(d1+dsp)` while Em gains d1+dsp. Space is magnitude-only — consistent with “space is a mode,” but then arg is defined on a slice of the new Em, not the full matter parcel.

4. **`amp_door` is a boolean.** Any `amp_door>0` fully enables the path; the value is never a mix weight (`amp_door` only appears in `> 0` tests). Naming implies a continuous knob; it is a switch.

5. **Evaporation is not phase-faithful in reverse.** `field_inject` (`freecell.c:922–934`): if Ee already >0, multiplies fa by a real factor (field phase unchanged); only vacuum inject uses th2. Matter→field does **not** compose arg the way field→matter does under amp_door. Bidirectional “atoms carry phase” is one-way.

### 1.7 Current / P1 meter

P1 uses ΔEm_j × link displacement (`freecell.c:1566–1572`). That is the right pairwise-conserved energy moment. It is **not** equal to `2τ Im(ψ_i*ψ_j)` after a finite Givens (higher orders in τ), but sign and leading order match. No bug found.

### 1.8 Stale docs

- `v93/README.md:1–7` still says “docs-only, NOT opened… awaiting sign-off” while code and CLAUDE.md say v93 is active with full authorization. Status banner is false.
- `L1_FINDINGS.md` ends with “arg(ψ) door is the registered next step” (`L1_FINDINGS.md:296–299`) while `REVIEW_REQUEST.md` already claims door results. Findings file is behind the brief.

---

## 2. L1-A narrow resonance — real physics vs artifact

### 2.1 What is solid

At amp_nat=2.0, 5/5 seeds give **+x** cos ∈ [0.75, 0.995] (`L1_FINDINGS.md:161–166`). Face 2 already had seed-robust **sign** within each amp_nat. Face 3 cleaned the Trotter wobble enough that one seed meets speed≥2.6e-3 and cos≈0.99. That is a real coherent directional bias from a phase-tilted dense packet under unitary hops — a first for the programme. I do not dismiss that.

### 2.2 What is not L1-A

L1-A as staged (`README.md:306–309`): speed ≥ 2.6e-3 **and** cos → 1, **seed-robust**.  
Measured (`L1_FINDINGS.md:168–171`): cos≥0.95 on 2/5; speed bar on 1/5. Calling the thesis “CONFIRMED” (`L1_FINDINGS.md:172–173`, `REVIEW_REQUEST.md:81–83`) confuses **mechanism engaged** with **acceptance met**. Honest label: **NEAR / mechanism-positive, bar open**.

### 2.3 Is the amp_nat≈2 window dispersion?

**Plausible yes**, for a load-dependent tight-binding caricature:
- on-site: ω_i = w2 / (1 + q_detune · x_i), x=(Em+flload)/cap (`freecell.c:1311–1314`)
- hop rate: t ~ amp_nat · kd · geo · gpl · res · gsym  
Group velocity of a wave packet with spatial frequency set by seed kx can reverse when t/ω crosses Brillouin-like regimes, or when load shear warps the local ω across the blob.

**But several artifacts can fake the same table:**

| Artifact | Why it matters | How to tell |
|----------|----------------|-------------|
| **Centroid of a diffracting free packet** | Asymmetric melt / bath coupling moves COM without bulk translation | Track Em-weighted COM **and** Em_tag, rms, peak density, second-moment velocity; require COM motion with rms not exploding |
| **Second-half linear fit only** | `blob_drift` fits only nfsamp/2…end (`freecell.c:4524–4541`); early transient or late melt dominate | Fit full trajectory; report R² of linear fit; split early/late |
| **cos = vx/\|v\| assumes kdir=+x** | `freecell.c:4541`: `vx/sp` — correct only if seed tilt is +x | OK for e3b if kx>0 always; document assumption |
| **τ saturation (0.5)** | Nonlinear map vs amp_nat | Histogram τ; if many links at clamp near “resonance edges,” not continuum dispersion |
| **gsym gate + seed phase noise** | Seeds with better core phase alignment get higher effective t | Correlate per-seed mean gsym / ρ_coh with cos |
| **bath=0 → zero motion** | Medium-required (`L1_FINDINGS.md:69–70`) | Translation may be **drag against bath**, not free ballistic packet — different physics than “unitary dense momentum” |

**Recommended discriminants (cheap):**
1. **Uniform clock probe** (already attempted as h3 with q_detune=0 on QUENCH; do it on e3b): if resonance **widens** and cos seed-variance falls, load-shear is causal.
2. **kx sweep at fixed amp_nat=2:** v_g(kx) should look like a sin-like dispersion branch if tight-binding; random scatter if centroid artifact.
3. **dt half / double at fixed amp_nat·dt product:** true rate physics is invariant; Trotter/clamp artifacts move.
4. **Report Em_tag(t), rms(t), peak Em, and COM(t)** for every e3b row — not only cos/speed.

### 2.4 Widening / pinning magnitude

Ranked:
1. **q_detune=0 on e3b** (uniform dense clock) — direct test of shear hypothesis in `L1_FINDINGS.md:185–189`.
2. **Softer gate or explicit odd channel:** τ ∝ gsym is even; optional `τ_signed ~ κ · base · (g_ij − g_ji)` or reactivate κ_reac-shaped odd piece as a hop phase (careful: must stay real rotation angle or complex SU(2)-style — don’t break norm).
3. **Do not reintroduce head in τ**; if anything, use a **soft floor** on gsym in the core (or gate on phase match only, not headroom).
4. **Separate amp_nat (overall scale) from a chart factor** so resonance is not the only knob that sets speed.
5. **Binding companion** (next section) — free diffraction guarantees seed-variant \|v\|.

---

## 3. QUENCH-3 retention wall

### 3.1 Steady-state result is real: no retained m=+2

Reanalysis of the FCS files (not the broken awk in `quench_argdoor.sh:13`):

| arm | R2d (t=20…300) | RA2 field |
|-----|----------------|-----------|
| h0 no door | ~0.01–0.04 | 1.0 → ~0.15 |
| h1 amp_door | ~0.01–0.04 | similar |
| h2 door+uni | ~0.02–0.05 | similar |
| h3 +q_detune=0 | ~0.02–0.06 | better field RA2 at mid times, still no matter lock |

Steady-state R2d≈0.02–0.03 for all arms is confirmed. Unitary transport alone does not retain winding — agreed.

### 3.2 Early-time “3× retention” — overstated

Early FCS (h0e / h1e / h2e):

| t | h0 R2d (no door) | h1 R2d (arg door) | h2 R2d (door+uni) | nlive |
|---|------------------|-------------------|-------------------|-------|
| 6 | **0.285** pk=m+2 | 0.196 pk=**m−1** | 0.195 pk=m−1 | 21 |
| 8 | 0.081 | 0.125 | 0.110 | 71 |
| 10 | 0.020 | 0.087 | 0.064 | ~138 |
| 12 | 0.023 | **0.143** | **0.133** | ~220 |
| 14 | 0.035 | **0.105** | 0.098 | ~316 |
| 20 | 0.041 | 0.010 | 0.018 | ~530 |

Comments:
1. At t=6, **no-door R2d is higher** than door. With nlive=21, 1/√N ≈ 0.22 — no-door 0.285 is **noise-floor consistent**, not a real winding imprint (no phase write without door).
2. Door arms do keep R2d ≳ 0.1 longer (t=12–14), then collapse by t=20 like everyone else. “Sustains ~3× longer” is a fair soft description of the **decay tail**, not a clean controlled effect size.
3. Door m-peak is often **not m=+2** (m−1, m+0, m+3…). You are not retaining the seeded texture; you are retaining **some** phase coherence that is not topologically the m=+2 vortex.
4. Field RA2 is still high (~0.8) through t=14 while matter already fails — so early loss is **not** only “field died in transit.” Imprint + matter dynamics fail while the field still carries spin.

### 3.3 Is “matter churn” the right diagnosis?

Partially, but incomplete.

**For churn:** births~deaths ~9e3 over T=300 on h0 (`h0_nodoor.log`) means many door cycles; if each condensation is a phase-slip site (design §II.8), r < 1 ⇒ Γ_p large. Agreed as a contributor.

**Against churn-as-sole-wall:**
1. **Field decoherence** is large (RA2 1→0.15) — design already notes this.
2. **Imprint is not m=+2-faithful** even early (peak mode wanders).
3. **Door is cell-th2, not slot-borne** — delivery churn of th2 (precession shear, empty-cell reset to 0, multi-atom overwrites) remains.
4. **h3 q_detune=0 exploded births** (33581 vs ~9e3) with **no** R2d gain — uniform clock without fixing imprint/topology made metabolism worse. That undercuts “just reduce shear.”
5. **Unitary arm halves births** (h2: 4580) but R2d unchanged — less churn did **not** buy retention. That is strong evidence against pure churn diagnosis.

**Better diagnosis (ranked):**
1. **No topological protection:** winding is protected only if \|ψ\|>0 on a cycle (`README.md:250–252`). Condensation into a sparse cloud with Em_LIVE threshold holes is full of phase-slip sites by construction.
2. **Wrong carrier:** cell th2, not slot-borne parcel phase.
3. **Field transit decoherence** + slit/bath geometry.
4. **Churn** as multiplier once 1–3 fail.

### 3.4 Experiments that distinguish

| Experiment | Churn-limited | Imprint-limited | Field-limited |
|------------|---------------|-----------------|---------------|
| Freeze geometry + disable birth/death (or par_tau identity only on long-lived cells) | R2d rises | R2d flat | R2d flat |
| Inject matter phase by hand at t=0 (seed th2 = 2φ, no field) + amp_nat only | retention tests transport | if dies, matter dynamics fail | N/A |
| Strong continuous spin drive (keep RA2 high with source) | if R2d still 0, imprint/hold fail | | if R2d follows RA2, field was limit |
| Slot-borne phase atom (real §II.7): deposit phase on slem/shadow, not th2 | should help if design right | | |
| Measure phase-slip rate per condensation event (pre-register §II.8) | | primary meter | |

### 3.5 Code/method path to real steady retention

1. **Implement slot-borne phase for the fired atom** (actual IV.6): on condensation, write phase into a slot buffer (existing `slph` / shadow `sre_/sim_` infrastructure) that decays with `amp_tau`, and let pass U / deposits read that — do **not** slam th2 every atom.
2. **Fix empty-cell th2 writeback** (above).
3. **Symmetric reverse door** on evaporation: compose field amplitude with √d2 e^{i th2} (with ledger renorm of Ee), not `field_inject` real scale when Ee>0.
4. **Accept torus language:** require a closed dense ring (nv parcel) with \|ψ\|>ε everywhere before claiming winding retention; QUENCH cloud is the wrong object.
5. **Stop claiming early R2d until analyze_winding m-peak is +2 and nlive ≳ 100.**

---

## 4. Packet binding — L1-A framing

**Agree with the free-packet diagnosis** (`L1_FINDINGS.md:74–79`, `REVIEW_REQUEST` item 4). A pure unitary hop with no confining potential is a discrete Schrödinger evolution: wave packets disperse; COM can still drift at v_g. Cap/head/inflight in the additive model were an accidental binder; Face 2–3 removed them from the hop.

**Therefore either:**

**A. Reframe L1-A (recommended short-term)**  
- Meter: **group velocity of a coherent dense packet** (phase-gradient current integrated, P1 dense moment, Fourier peak of Em e^{i th2}).  
- Drop “bound blob cos→1 seed-robust” as the first bar.  
- Keep “direction seed-robust in a window” as L1-A0; ballistic bound object becomes L1-A1 after binding.

**B. Add a binding companion (medium-term)**  
Candidates that respect safeguards:
- **Cap + radiance fixed point** already binds energetically (x̂*=0.62) but does not bind **phase texture**.
- **Do not** put head back into τ (Face 1 lesson).
- Possible: weak **attractive bond on phase-matched dense pairs** (existing κ_bond is geometric; a dense-phase-matched cohesion on positions or on Es footprint).
- Or: **inflight not fully killed** — keep a small additive flight channel for localization while unitary hop carries the coherent current (hybrid). Risky but battery-testable behind a knob.
- Soliton from **nonlinear on-site ω(Em)** (q_detune) + hopping is the classical DNLS mechanism; you already have load-detuning. The narrow resonance may be a hint that **self-trapping is adjacent** — map the existence region carefully rather than treating diffraction only as failure.

**bath=0 → zero motion** is important: current “translation” may require a bath to push against. A true mover should show dipole of J with balanced monopole **in vacuum or thin bath**. That experiment is missing and should gate any “translation IS the current” slogan (`README.md:216–219`).

---

## 5. Missed aspects / alternative explanations

1. **L1-C “byte-identical” is stronger than “self-gating τ≈0”.** With p_gate=8, cos^8 is extremely sharp — random bath truly has τ~0. Coherence-runaway risk is **not** fully tested: the armed case is **ρ_coh≈0.77 coherent bath**, not random. Design §II.9 named that risk; L1-C measured the unarmed bath. **Gap:** run a high-coherence bath (or Phase-M ρ_coh-seeded) at amp_nat=2–5 and watch births/glow.

2. **Driven QUENCH births drop at amp_nat>0** (`L1_FINDINGS.md:216–218`) — channel reshapes condensation. Not ignition, but not “inert on structure” either. Metabolism change without retention is a new regime to document.

3. **L2 ggm 0→0.13 at amp_nat=2 then 0.015 at 3** (`L1_FINDINGS.md:267–272`) has the **same resonance shape as L1-A**. That is either beautiful (one coherent mechanism) or a red flag (one amp_nat accident). Joint amp_nat sweep of e3b cos and tri ggm on the same law table is mandatory.

4. **C≡Go at strength:** design III.3 warned tolerance abx; battery green is at amp_nat=0. Production sweeps at amp_nat=2 need explicit C vs Go on e3b cos/speed, not only inert abx.

5. **Analysis bug:** `quench_argdoor.sh:13` awk `$(NF-3)` does not extract R2d (see `quench_argdoor.out` printing `R2d=M:` / `R2d=pk=m+4`). **Do not trust suite stdout; re-derive from analyze_winding.**

6. **README Re/Im inconsistency** and **stale STATUS banner** erode programme hygiene.

7. **Alternative to “dispersion” for direction reverse between amp_nat 2 and 2.5:** second-half COM fit on a **wobbling/melting** packet can flip fitted vx when the late-time blob fragments. Without trajectory plots this is not excluded.

8. **Synthesis insight still stands:** additive ledger rejecting the interference term is the right post-mortem of v92 feedback. Unitary hops are the correct response. That does **not** automatically deliver L1-A or spin retention.

---

## 6. Concrete next work (ranked)

### P0 — correctness / honesty (do before more physics claims)
1. Fix empty-cell `th2` writeback under pass U (leave phase if Em≈0).
2. Fix `quench_argdoor.sh` R2d extraction; re-publish early/late tables with nlive and m-peak.
3. Update `L1_FINDINGS.md` with arg-door results **and** downgrade “thesis CONFIRMED” to “mechanism engaged; L1-A open.”
4. Fix design doc: J = Im not Re; STATUS = active programme; document that amp_door is cell-level, not slot-borne.

### P1 — decide what L1-A is
5. e3b **instrumentation pack**: COM(t), Em_tag(t), rms, peak Em, linear-fit R², mean τ, clamp fraction.
6. **q_detune=0** and **kx sweep** at amp_nat=2 (dispersion vs artifact).
7. **bath=0 / thin bath** dipole-of-J test.
8. Reframe bars: L1-A0 group velocity; L1-A1 bound object later.

### P2 — retention (real §II.7)
9. Slot-borne fired-atom phase (use shadow / slph), th2 only as bulk clock.
10. Symmetric coherent merge on evaporation.
11. Hand-seeded m=+2 matter ring (no field) + unitary transport — isolates hold vs imprint.
12. Phase-slip rate meter at doors (pre-registered acceptance language).

### P3 — binding / fifth
13. Search DNLS-like self-trapping (q_detune × amp_nat map) rather than only “add a binder.”
14. Joint L1-A / L2 amp_nat fine sweep; seed a live fifth for L2 maintenance (IV.8).
15. Coherent-bath anti-ignition (armed L1-C).

### P4 — cleanup
16. Optional odd component in τ (κ_reac-shaped) as a controlled experiment.
17. Tolerance C≡Go at amp_nat=2 on physics columns.
18. DEC δJ refactor only after L1-A0 is clean (design III.4 step 5 is rightly last).

---

## 7. Answers to the six brief questions (compressed)

1. **Pass U sound?** Yes for Givens + pairwise norm. Precession ordering right for dense history; not same sign as F. arg(ψ) door is Em-ledger-safe but **not** unitary merge and **not** slot-borne. No hop sign bug; design Re/Im docs are wrong. Empty-cell th2→0 is a real defect.

2. **Narrow resonance?** Likely mix of load-dependent tight-binding dispersion **and** free-packet + fit artifacts. Discriminate with q_detune=0, kx sweep, dt invariance, full COM trajectories. Widen via uniform clock and/or binding; don’t chase amp_nat as a magic constant.

3. **QUENCH wall?** Churn is **not** sufficient (unitary arm reduces births, R2d unchanged). Imprint quality + cell-clock carrier + field transit + lack of topological support dominate. Slot-borne phase + ring substrate + reverse door are the code path; churn reduction alone will not yield steady R2d.

4. **Binding?** Necessary for a **bound translating blob**. For a first success bar, reframe L1-A as wave-packet group velocity / current dipole. Parallel track: exploit nonlinear ω(Em) for self-trapping.

5. **Missed:** armed coherence L1-C; bath-required motion; L2 resonance co-located with L1-A; analysis script bug; amp_door boolean; evaporation phase asymmetry; flight sector disabled under amp_nat.

6. **Next:** P0 honesty fixes → P1 L1-A diagnostics → P2 real slot-borne door → P3 binding/fifth.

---

## 8. Bottom line

v93’s core engineering move — **replace additive Em want with pairwise unitary hops** — is sound, battery-safe at amp_nat=0, conservation-theoretic (L1-B), and self-gating on incoherent matter (L1-C). Face 3’s in-pass precession is a real improvement and produces the first seed-robust **directional** dense bias.

It does **not** yet confirm “coherent dense translation” as an accepted L1-A object, does **not** implement the designed arg(ψ)/slot-borne retention fix, and the QUENCH early-time story is softer than the brief states. The programme should tighten meters and close the design–code gap on the door before claiming structural victory beyond conservation and anti-ignition.

Disagreement with the author’s conclusion (1): **partial** — structural wins on B/C yes; “translation emerges / thesis confirmed” is one notch too strong.  
Agreement with conclusion (2): door helps only weakly; wall is not primarily churn.  
Agreement with conclusion (3): resonance, seed variance, diffraction are the open L1-A gaps — add measurement confounds and bath-dependence to that list.
