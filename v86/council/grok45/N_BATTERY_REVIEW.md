# N-battery review — adopted Part 0 vs intent
**Seat:** grok45 · **Date:** 2026-07-26  
**Inputs:** `v86/GROUNDING.md` v2 (§0, §5), `v86/PROPOSAL.md` Amendment 1, prior `FOUNDATIONS.md` §3 + `FEEDBACK.md`, `v85/STATE_OF_THE_THEORY.md`  
**Charge:** (1) confirm/correct N1–N10 as adopted; (2) order-of-work check; (3) open invitation on strategy.

Ground rules: quote attacked text; number findings; CANNOT VERIFY where applicable.

---

# 1. Confirmation / correction of N1–N10

## 1.0 Global adoption issues

**Attacked block (GROUNDING v2 §5 header):**
> "All shooter/analysis-level (zero GPU spend) except N7:"

followed by ten one-line bullets.

**Finding 1.** Adoption preserves **names and scientific intent** of every rung. It **strips** the load-bearing operational content: exact observables, formulas, discriminator tables, pass/fail (kill) criteria, data paths, and costs. As written, §5 is a **todo list**, not an executable protocol. Implementers will invent incompatible operationalizations; results will not be comparable to the FOUNDATIONS design or to each other.

**Finding 2.** PROPOSAL Amendment 1 cites "N1–N10 as specified in GROUNDING v2 §5" — so the thin §5 text is now the **sole specification**. Restoring detail is not optional documentation polish; it is the difference between a battery and a slogan.

**Finding 3.** Amendment 1 order **omits N9 entirely** from the revised schedule (lists N1–N4, N10, N5/N6/N8, N7 — never N9). Either N9 was silently dropped or the order is incomplete. Flag as a planning bug.

Below: for each rung — **faithful?** / **sufficient to discriminate?** / **restored protocol** / **ambiguity flags**.

---

## N1 — Component decomposition of E / Σ

| | |
|---|---|
| **Adopted text** | "component decomposition of E on existing gscan + branch tables → Σ split into E_∇, E_m, E_V, E_g vs Q (is Σ the virial imbalance as claimed?)" |
| **Faithful to intent?** | **Yes** (core question). |
| **Sufficient as written?** | **No** — no formulas, no first-cut vs full-cut, no discriminators, no kill. |

**Restored protocol**

- **Cost:** analysis, hours–1 day; zero new runs for first cut. Full radial split may need a short shooter instrument pass (CPU, free).
- **Data:** `/home/d/code/scp/v69/theory/gscan.tsv` columns `omega, Q, E_matter, E_field, E_total, r_half` (header verified present).
- **Definitions (single-clock):**
  \[
  \Sigma \equiv E - \omega Q,\qquad \varepsilon \equiv \frac{\Sigma}{\omega Q} = \frac{E}{\omega Q}-1.
  \]
  Target partition (from radial energy density, v66/THEORY):
  \[
  \begin{aligned}
  E_{\mathrm{kin}} &= \tfrac12 \omega Q &&\text{(throughput; check vs integrated }\tfrac32\omega^2 f^2\text{)},\\
  E_{\nabla} &= \tfrac32\int (f')^2\,\mathrm{d}V,\\
  E_m &= \tfrac32\int m^2 f^2\,\mathrm{d}V,\\
  E_V &= \int V_t(f^6)\,\mathrm{d}V,\\
  E_g &= \text{gauge energy} = E_{\mathrm{field}}\text{ in gscan}.
  \end{aligned}
  \]
  Consistency identity to report: \(E \stackrel{?}{=} E_{\mathrm{kin}}+E_{\nabla}+E_m+E_V+E_g\) (convention-matched to kernel factors).
- **Phase A (today, gscan only):** \(\Sigma\) vs \(E_{\mathrm{field}}\); \(\varepsilon(g{=}0)\) vs \(\varepsilon(g{>}0)\) gap attributed to \(E_g\).
- **Phase B (shooter emit):** full \(\{E_{\nabla},E_m,E_V,E_g\}\) fractions of \(\Sigma\).

| Outcome | Interpretation |
|---|---|
| \(\Sigma\) tracks \(4\pi R^2\) (+ Coulomb \(\propto g^2Q^2/R\)) | Surface/Coulomb origin |
| \(\Sigma\) dominated by bulk \(E_V\) | Potential-shape / thick-wall |
| \(\varepsilon\to 0\) as \(R/\xi\to\infty\) on g=0 thin wall | Thin-wall asymptote confirmed |
| \(E_g\) explains most of g>0 vs g=0 ε gap | Coulomb is the g-piece |
| Phase-B sum of parts fails to rebuild \(E\) at ≫ shooter tol | **Kill / alarm:** energy bookkeeping broken before ontology talk |

**Ambiguity if left compressed:** Phase A alone (matter vs field) **cannot** answer "is Σ the virial imbalance?" — that requires \(E_{\nabla}\) vs \(E_V\) vs kinetics. Adopted text implies full split; gscan alone does not deliver it. **Require Phase B** for a non-ambiguous N1 verdict.

---

## N2 — Integrated virial identity

| | |
|---|---|
| **Adopted** | "integrated virial identity of the radial ODE (parameter-free) → ε(Q) as a derived function; TOE-grade if it lands." |
| **Faithful?** | **Yes.** |
| **Sufficient?** | **Almost** — missing residual tolerance and fail modes. |

**Restored protocol**

- **Cost:** theory derivation half-day + evaluation on existing shooter profiles (free).
- **Work:** (i) multiply radial profile ODE by \(f\) and/or \(r f'\); integrate; obtain parameter-free relation among \(E_{\nabla},E_m,E_V,\omega,Q\) (and \(E_g\) if gauged). (ii) Evaluate residual \(R_{\mathrm{vir}}\) on every gscan/shooter point.
- **Pass:** \(|R_{\mathrm{vir}}|/E \ll \varepsilon\) (target: at or below shooter ODE residual, typically ≪ 10⁻³ relative).
- **Fail / kill:**
  - Residual ~ ε but systematic → identity wrong or incomplete (missing gauge/surface term).
  - Residual ≫ ODE residual → energy definition inconsistent with Lagrangian used to define \(E\) — **foundational alarm for D5′** (cannot trust \(M=E/c^2\) if \(E\) is not the Hamiltonian that generates the ODE).
- **TOE-grade landing:** \(\varepsilon(Q)\) expressed as a function of measurable profile moments with no free fit parameters.

**Ambiguity:** "TOE-grade if it lands" without residual threshold invites declaring victory on a 1% residual that is the whole ε signal.

---

## N3 — Scaling collapse

| | |
|---|---|
| **Adopted** | "scaling collapse: ε vs ξ/R across g after Coulomb subtraction." |
| **Faithful?** | **Yes**, slightly incomplete (second axis). |
| **Sufficient?** | **No** without axes, g-set, and subtraction recipe. |

**Restored protocol**

- **Cost:** shooter re-drive or re-use gscan rows for g∈{0, 0.02, 0.05, 0.10}; analysis free.
- **Definitions:** \(\xi = 1/\sqrt{m^2-\omega^2}\); \(R = r_{\mathrm{half}}\) (or charge RMS — **freeze one** and report sensitivity).
- **Coulomb subtraction (gauged):** form \(\varepsilon' = (E - E_g^{\mathrm{(Coulomb\ est)}})/(\omega Q) - 1\) with \(E_g^{\mathrm{est}} = g^2 Q^2/(8\pi R)\) (or measured \(E_{\mathrm{field}}\) from N1); also plot raw \(\varepsilon\) vs \(g^2 Q/R\).
- **Plots:** (1) \(\varepsilon\) vs \(\xi/R\); (2) \(\varepsilon\) vs \(g^2 Q/R\); (3) \(\varepsilon'\) vs \(\xi/R\) after subtraction.

| Outcome | Interpretation |
|---|---|
| Collapse on \(\xi/R\) alone (all g) | Pure geometry / surface |
| Needs \(g^2 Q/R\) or only collapses after Coulomb subtraction | Coulomb self-energy essential |
| No collapse | Bulk / thick-wall / potential-shape dominates (still OK; changes closed-form target) |

**Ambiguity:** "after Coulomb subtraction" without specifying measured \(E_{\mathrm{field}}\) vs point-Coulomb estimate can absorb arbitrary residuals. Prefer measured \(E_g\) from N1.

---

## N4 — Multi-flavor Legendre excess

| | |
|---|---|
| **Adopted** | "multi-flavor Legendre excess on the flavored branch." |
| **Faithful?** | **Partially** — name only. |
| **Sufficient?** | **No.** |

**Restored protocol**

- **Cost:** flavored shooter scan (CPU free); partitions as in v71/QRK (e.g. symmetric vs 76.7/103.7/103.7-class).
- **Observable:** \(\Sigma(\vec\omega) = E - \sum_a \omega_a Q_a\) at fixed total \(Q=\sum Q_a\), vary partition \(\vec Q\).
- **Discriminator:**

| Outcome | Interpretation |
|---|---|
| \(\Sigma\) ≈ constant vs partition at fixed \(Q\) | Excess is geometric / total-charge |
| \(\Sigma\) strongly partition-dependent | \(V(s)\) interference is a leading ε channel → feeds HC-3 / census |
| \(\Sigma < 0\) somewhere | Sign convention or non-minimizing branch — investigate before GSS |

- **Kill:** none for foundation; **gate for HC-3** if partition dependence is large — then single-charge ε formulas must not be naively reused.

**Ambiguity if compressed:** "Legendre excess on the flavored branch" can be satisfied by printing one number at one partition — zero discrimination.

---

## N5 — Throughput vs charge (𝒜 ≠ 2πQ)

| | |
|---|---|
| **Adopted** | "throughput-vs-charge (𝒜 vs 2πQ) on flavored/breathing objects — the action≠charge split made empirical." |
| **Faithful?** | **Yes.** |
| **Sufficient?** | **No** — missing observables and data paths. |

**Restored protocol**

- **Cost:** analysis of existing diags; free if kinetic energy is in diag or recoverable from SFA.
- **Data:** `/space/scp/v85/topo1/out/qrk*`; `/space/scp/v85/topo1/gpu/x10c*` (breathing).
- **Observables (report all three):**
  1. \(T \equiv 2 E_{\mathrm{kin}}\) (cycle-averaged if nonstationary),
  2. \(S_\omega \equiv \sum_a \omega_a Q_a\) (or \(\omega Q\) monochromatic),
  3. \(2\pi Q\) (and \(\hbar_{\mathrm{eff}}\) candidates \(T/\omega\), \(E/\omega\)).
- **Process action definition to freeze:** \(\mathcal{A} \equiv \oint 2 E_{\mathrm{kin}}\,\mathrm{d}t\) over one clock period (or \(T/\omega\) for stationary). Compare \(\mathcal{A}/(2\pi)\) to \(Q\).

| Outcome | Interpretation |
|---|---|
| \(T = S_\omega\) within diag noise on flavored | Throughput generalizes (R4 soft) |
| Systematic \(T \neq S_\omega\) on flavored/breathing | Action splits from charge (R4 hard) |
| Cloud sector: field energy \(\approx \omega_{\mathrm{cl}}|Q_{\mathrm{cl}}|\) | Supports linear-mode pillar (ties N10) |

**Ambiguity flags (high):**
1. **Does current diag expose \(E_{\mathrm{kin}}\)?** If not, N5 requires SFA post-processing or a diag extension — cost not "free." **CANNOT VERIFY** without reading current diag columns in topo1 files.
2. Breathing X10c is nonstationary; cycle average must be defined (window ≫ breath period ~500–600 t.u.) or the residual is phase noise, not ontology.

---

## N6 — ħ_eff triple comparison

| | |
|---|---|
| **Adopted** | "ħ_eff triple comparison (reopen v70 calibration against ε)." |
| **Faithful?** | **Yes.** |
| **Sufficient?** | **Borderline** — needs the three quantities named. |

**Restored protocol**

- **Cost:** re-analysis of v70 boost series if archived; else one small boost set (CPU/GPU-light).
- **Triple (same object, same time window):**
  \[
  \hbar_E \equiv \frac{E}{\omega},\qquad
  \hbar_{pk} \equiv \frac{p}{k},\qquad
  \hbar_Q \equiv Q,
  \]
  with \(E=\int T_{00}\) **only** (D5′; do not use \(Q\omega\) as \(E\)).
- **Report:** relative residuals \(\hbar_E/Q-1\), \(\hbar_{pk}/Q-1\), \(\hbar_E/\hbar_{pk}-1\); overlay vs ε from rest-frame branch.
- **Pass for "identity":** all residuals below claim precision (program was quoting 3–5%).  
- **Expected (FOUNDATIONS):** residuals same family as ε (1–4%); **kill the identity language**, keep measured ratios.

**Ambiguity:** "against ε" without requiring \(E=\int T_{00}\) can re-introduce circular Qω-based ħ_eff.

---

## N7 — N-INV weak acceleration (empirical D5′)

| | |
|---|---|
| **Adopted** | "N-INV weak-acceleration inertia measurement — decides D5′ empirically (M_inertial vs E/c² vs Qω/c²). GPU-light." |
| **Faithful?** | **Yes** — this is the load-bearing foundation experiment. |
| **Sufficient?** | **No** — highest ambiguity risk of the battery if protocol stays one line. |

**Restored protocol**

- **Cost:** CPU N=64 first; GPU-light only if force control needs resolution. Design-only until force protocol chosen.
- **Setup:** stationary gauged ball on branch; apply **known** weak force:
  - Option A: very slow boost ramp (global velocity field increase), measure COM acceleration during ramp; or
  - Option B: uniform external E-field analog if kernel supports / slow gauge bias; or
  - Option C: two-ball distant Coulomb pull with measured \(g^2 Q_1 Q_2/4\pi R^2\) as \(F\).
- **Observable:** \(M_{\mathrm{inertial}} = F / a_{\mathrm{COM}}\) in the linear-response window (small \(a\), no radiation-dominated drag).
- **Compare to:** \(E/c^2\), \(Q\omega/c^2\), \(Q\omega(1+\varepsilon)/c^2\) with \(c=1\) lattice units.

| Closest match | Consequence |
|---|---|
| \(E\) | D5′ confirmed; H4 dead as mass definition |
| \(Q\omega\) | Demotion narrative partially revived — force SR rewrite |
| Neither (O(ε) off both) | Medium-added mass / lattice form factor — new physics, report \(M_{\mathrm{inertial}}/E\) |

- **Kill / invalid run:** sponge-dominated drag (force not equal to applied \(F\)); force large enough to excite monopole breath at O(1); COM defined on density that includes radiation wake.

**Ambiguity flags (critical):**
1. Lattice preferred-frame and 1–5% group-velocity anomaly (STATE/v70) mean \(M_{\mathrm{inertial}}\) may disagree with \(E\) at the **same scale as ε**. Pre-register: agreement within anomaly budget = pass for D5′; do not require sub-percent.
2. Without a force calibration independent of assuming \(M=E\), option C (Coulomb) is partially circular if \(g\) was fit assuming energy units. Prefer kinematic boost-ramp (option A) where \(F\) is not assumed.
3. **Adopted form alone would produce an ambiguous result** — restore protocol before scheduling GPU.

---

## N8 — Ring spin ladder M(n) at fixed Q

| | |
|---|---|
| **Adopted** | "ring spin-ladder M(n) at fixed Q (v73 ring machinery) — mass without charge; the first foundation-level degeneracy crack test." |
| **Faithful?** | **Yes.** |
| **Sufficient?** | **Mostly** — add scan and kill numbers. |

**Restored protocol**

- **Cost:** reuse v73 ring seed/tools; shooter or lattice stationary states at fixed \(Q\), scan integer winding \(n\).
- **Measure:** \(E(n)\), \(\omega(n)\), \(\Sigma(n)=E-\omega Q\), \(L_z=nQ\).
- **Discriminator:** fractional mass ladder \(\Delta E / E\) per \(\Delta n=1\) at fixed \(Q\).

| Outcome | Interpretation |
|---|---|
| \(\partial E/\partial n\) large (mass moves ≫ ε at fixed Q) | Degeneracy cracked (R5.2 win) |
| \(E\approx\omega Q\) still, \(\partial E/\partial n\sim 0\) beyond spin kinetic | Degeneracy intact; rings do not buy hierarchy |
| No stable ring branch at needed \(n\) | Test blocked — report, do not claim crack |

**Ambiguity:** "mass without charge" rhetoric if \(E\) barely moves — require a quantitative threshold, e.g. **claim crack only if** \(\Delta E/E > 2\varepsilon\) over accessible \(n\).

---

## N9 — Soft-window (μ,κ) scan

| | |
|---|---|
| **Adopted (GROUNDING §5)** | "soft-window shooter scan (μ,κ) — how much can M/Q move in one sector?" |
| **In Amendment 1 order?** | **MISSING** (Finding 3). |
| **Faithful?** | **Yes** if restored to schedule. |
| **Sufficient?** | **No** without numeric HIER bound. |

**Restored protocol**

- **Cost:** 2D shooter grid, free CPU.
- **Scan:** \((\mu,\kappa)\) about standards \((-41.345, 50)\); map existence window width, \(\varepsilon\) range, and \(E/Q\) dynamic range **at fixed Q** (or fixed mid-branch ω).
- **Discriminator (FOUNDATIONS):**

| Outcome | Interpretation |
|---|---|
| \(E\) at fixed \(Q\) varies by ≳ **2×** inside VK-stable region | Single-sector HIER softens; second sector less mandatory for hierarchy alone |
| Window always ≲ **10%** in \(E/Q\) | Single-sector hierarchy blocked; §D.2 second sector mandatory for HIER |

**Ambiguity:** "how much can M/Q move" without thresholds cannot decide second-sector necessity — the whole point of N9.

---

## N10 — Shell-mode E=ωQ exactness

| | |
|---|---|
| **Adopted** | "shell-mode E = ωQ exactness on X10 profiles (the linear-mode anchor, atom-facing)." |
| **Faithful?** | **Yes.** |
| **Sufficient?** | **No** — core subtraction is load-bearing. |

**Restored protocol**

- **Cost:** analysis of X10 pathfinder/b/c diags + optional KG eigenmode from ANALYSIS §2D; free if energy partitions exist.
- **Data:** `/space/scp/v85/topo1/gpu/x10*`.
- **Must measure cloud-only energy:** \(E_{\mathrm{cloud}} = E_{\mathrm{total}} - E_{\mathrm{core\ alone}}\) (matched Q_N ball without cloud seed), **or** integrate energy density outside a core radius \(r > r_{\mathrm{contact}}\), **or** use linear KG mode normalization.  
  **Do not** compare full ball+cloud \(E\) to \(\omega_{\mathrm{cl}}|Q_{\mathrm{cl}}|\).
- **Pass:** \(|E_{\mathrm{cloud}} - \omega_{\mathrm{cl}}|Q_{\mathrm{cl}}|| / |E_{\mathrm{cloud}}| \ll \varepsilon_{\mathrm{soliton}}\) (target: order of numerical noise / linearization error, not 1–4%).
- **Fail at soliton-ε level:** cloud nonlinear, or subtraction wrong, or identity misapplied.

**Ambiguity flags (high):** Adopted one-liner almost guarantees people will ratio total energy to cloud charge and get garbage. **Without core subtraction, N10 is invalid.**

---

## 1.1 Summary table

| Rung | Intent preserved? | Discriminating if §5-only? | Ambiguity risk | Priority restore |
|---|---|---|---|---|
| N1 | Yes | No (need Phase B) | Medium | formulas + Phase A/B |
| N2 | Yes | Almost | Medium | residual ≪ ε |
| N3 | Yes | No | Medium | axes + Coulomb recipe |
| N4 | Partial | No | High if one-shot | partition scan table |
| N5 | Yes | No | **High** (diag/kin) | observables + windows |
| N6 | Yes | Borderline | Medium | triple + E=∫T₀₀ |
| N7 | Yes | **No** | **Critical** | full force protocol |
| N8 | Yes | Borderline | Medium | ΔE/E > 2ε threshold |
| N9 | Yes (text) / **dropped from order** | No | Medium | thresholds + schedule slot |
| N10 | Yes | **No** | **Critical** | core subtraction |

**Finding 4.** Rungs most likely to produce **ambiguous results in adopted form:** **N7, N10, N5** (protocol-sensitive), then **N4, N9** (under-specified discriminators). N1–N3 are safe if Phase B and residual cuts are restored.

**Finding 5.** GROUNDING §0 foundation freeze is **faithful** to FOUNDATIONS R1/R2/R4 (primary E, M=E/c², Σ excess, three-way split, degeneracy as regime). No correction required to §0 content — only to §5 executability.

---

# 2. Future direction / order of work

## 2.1 Attack the still-unamended PROPOSAL body

**Attacked (PROPOSAL I.1, still in file above Amendment 1):**
> "Prediction: the number of independently stable harmonics = the dimension of the conserved-charge lattice; everything else decays."

**Attacked (PROPOSAL I.1 resonance):**
> "n·ω_i ± m·ω_clock stays below the radiation continuum (gap m=1.5)."

**Finding 6.** Amendment 1 **does not rewrite Part I** — it appends. The body still carries the slogan and continuum formula that FEEDBACK and GROUNDING v2 §2 already killed (multipole-first; BdG edge, not m=1.5 on QRK Ω). **Correction required:** either strike/replace I.1–I.4 in place or mark body **SUPERSEDED by Amendment 1 + GROUNDING v2**. Leaving both invites implementers to follow the wrong text.

## 2.2 Confirm or correct revised order

**Attacked (Amendment 1 order):**
> "Part 0 (N1–N4, N10) → HC-1/HC-2 → N5/N6/N8 → HC-3/HC-6 → HC-4 (hardened) → EX-2 → N7 → EX-1/EX-4 (ramped) → EX-3."

**Finding 7 — Directionally correct, wrongly sequenced in three places.**

| Issue | Why |
|---|---|
| **N7 too late** | Empirical D5′ is foundational; EX-1/4 and any mass-defect claim consume it. Parking N7 after HC-4 + EX-2 delays the freeze §0 promised. N7 can run **in parallel** after N1 once E is trusted. |
| **N9 omitted** | Soft-window HIER bound is the single-sector vs second-sector fork; belongs after N1–N3, before second-sector design talk. |
| **N8 before HC-3** | Acceptable, but N8 needs stationary rings; do not block census if ring tools are sticky — **non-blocking parallel**. |
| **HC-4 after HC-3/6** | Correct (need multipole catalog + unstable control). |
| **EX-2 before EX-1** | Correct (mine existing flux first). |
| **EX-3 last** | Correct as audit; see cuts below. |

**Corrected order (recommended):**
1. **N1 Phase A → N2 → N1 Phase B → N3 → N10** (ε/linear-mode freeze; days, $0)
2. **N4 + N9** (flavored Σ + window HIER bound; shooter free) — *N9 restored*
3. **HC-1** (BdG + multipoles + continuum edges + n(H_ω))
4. **HC-2** (multipole-first audit)
5. **N5, N6** (existing data; action/ħ honesty) — *parallel with HC-1 if staffed*
6. **N7** (D5′ lock) — *start design at step 1; run as soon as force protocol ready, before EX-1*
7. **HC-3 + HC-6** (GSS region + converse)
8. **HC-4 hardened** (only after HC-1 edges + width floor)
9. **EX-2**
10. **EX-1 / EX-4 ramped** (consumes N7 + EX-2)
11. **N8** non-blocking whenever ring tools free
12. **EX-3** only if EX-1 claims "movement correct" at precision better than known 1–5% anomaly

**Finding 8 — Highest information per cost path (critical path):**
\[
\textbf{N1A}\to\textbf{N2}\to\textbf{N1B/N3}\to\textbf{N10}\;\Big\|\;\textbf{HC-1}\to\textbf{HC-2}\to\textbf{HC-3/HC-6}
\]
then **N7** (single decisive dynamics rung) → **EX-2** → **minimal EX-1** (one ramped v≤0.03).  
Everything else is second-order for the standing carbon goal.

**Finding 9 — If forced to halve the program, cut entirely:**

| Cut | Rationale |
|---|---|
| **EX-3** | Known 1–5% lattice anomaly already bounds "simulated correctly"; full lab-vs-boost dx set is expensive certification, not discovery |
| **N9** (conditional) | Cut only if N1–N3 + existing gscan already show ω-window ≲10% forever under standard potential — else keep (cheap!) |
| **HC-5** | Merge into HC-6 (unstable partitions already test overload-to-minimum); dedicated 2-frequency seed is redundant |
| **EX-4 full cloud-mode transport** | Already deferred in GROUNDING §6; keep ball-clock only |
| **EX-1 v=0.05 arm** | Keep single ramped v=0.02–0.03; upper arm is stripping theater |
| **N8** (defer, not kill) | Degeneracy crack is high value but not on carbon critical path; park after Part 0 freeze |
| **N6 full new boost series** | Re-analyze v70 only; no new campaign |

**Do not cut:** N1–N3, N7, N10, HC-1, HC-3/6, EX-2, one EX-1.

**Finding 10.** Part 0 **before** census is correct and should stay. Census **before** heavy exchange is correct. Foundation battery **before** Stage-4 carbon is correct. The only serious ordering error is **late N7** and **ghost N9**.

---

# 3. Open invitation — strategy from the measured state

## 3.1 Strongest alternative strategies

**Alt-A — Carbon endgame first (STATE §3A).**  
Skip deep C1 census and much of Part 0. Run 2B multi-center + D12-clean cloud at neutrality (2–3 GPU). Success = structural Z-mapped atom.  
*Strength:* matches standing goal; physics for shells already PASSED (X10).  
*Fatal weakness:* inherits Q-degeneracy, continuous cloud electrons, no transport audit, no ε ontology freeze — every "atom" claim will be re-litigated. EX rule ("parked atom that cannot move is not an atom") is right; pure Alt-A violates it.

**Alt-B — Second-sector kernel first (STATE §3B / TOE §D.2).**  
Design light-gap sector sharing gauge U(1); authorize kernel surgery; aim at true HIER + orbital carriers.  
*Strength:* only structural fix for nature-like mass hierarchy and universal ħ chain.  
*Fatal weakness:* unauthorized until design doc; does not remove need for Part 0/census in *each* sector; high risk of multi-fab redux (v75/v83 scars) without tighter success criteria.

**Alt-C — Numerics hygiene / D7 first (STATE §4).**  
X2 undone; async SFA corruption; seed-shape mismatch; unaudited dx. Certify one ball CLOSED vs (dt,dx,L) before any width or boost science.  
*Strength:* every HC-4/EX result is otherwise confounded by the same systematics the archive already flags.  
*Weakness:* pure hygiene does not move carbon or foundation theory; can become infinite regression.

**Alt-D — Spectroscopy-only theory export (C1 without carbon).**  
Finish multipole BdG catalog + GSS + widths; publish "fabric particle physics" without Stage-4.  
*Strength:* cleanest intellectual export; matches QRK opening.  
*Weakness:* abandons standing carbon goal; program identity is atom-from-fabric, not Q-ball spectroscopy.

## 3.2 Would I pursue a fundamentally different direction?

**Finding 11 — Not fundamentally different; materially reweighted.**

The **current plan's skeleton is right**:
1. Freeze energy–mass foundation from measurement (Part 0 / §0) — mandatory after the E=Qω overclaim.
2. Mode census with multipole-first + corrected GSS — the right theory target for "why stable objects look like lowest states of sectors."
3. Exchange audit before calling anything an atom — STATE and EX framing are correct.

I would **not** replace this with Alt-B first (kernel surgery without Part 0) or pure Alt-D (abandon carbon).

I **would** reweight toward a hybrid that Alt-A and Alt-C force:

### Recommended strategy (seat proposal)

**Name:** *Foundation freeze → certified ball → minimal carbon → census depth as theory export*

| Phase | Content | Success criterion |
|---|---|---|
| **P0** | N1–N3, N10, N7 (and N9 if cheap) | TOE Step 3 rewritten to §0; ε story closed or bounded; D5′ empirical |
| **P0b** | **D7-lite on one object** (not currently in Amendment 1): ball drift vs dx and L only (skip full 4-axis theater) | Drift decomposition: residual vs ARTIFACT(dx), BC-limited(L) |
| **P1** | HC-1, HC-2, HC-3, HC-6 | C1 at "organizing claim" strength, not theorem theater |
| **P2** | EX-2 + one ramped EX-1 | Moving ball+cloud co-motion number; strip fraction at v=0.02 |
| **P3** | **Carbon structural** (STATE A): 2B droplet or multi-center + D12 cloud, g=0.023/Q_N=650 class | Neutral cloud retained at shell radii on Z-mapped core; process ledger; **no claim of 6 electrons or orbits** |
| **P4** | Deep HC-4 / EX-3 / N8 / second-sector design doc | Only if P3 succeeded and program still wants spectroscopy export or HIER fix |

**Finding 12 — What I discard or demote vs current plan:**

1. **Demote C1 "theorem-shaped" (GROUNDING §4)** to **C1-weak:** measured organization of long-lived spectrum by GSS + multipole + arithmetic, with listed exceptions (monopole breath). Do not hold Stage-4 hostage to full HC-4 ε² theater.
2. **Demote full HC-4** until P0b D7-lite exists — otherwise width science is box fiction (GROUNDING §2 already admits IR cutoff hazard).
3. **Promote D7-lite** into the program explicitly (Amendment gap).
4. **Promote one carbon GPU package earlier** (after P2, not after endless EX-3) — standing goal is carbon, not perfect exchange metrology.
5. **Keep second sector as design-only** until N9 + carbon attempt show single-sector ceiling; do not let HIER rhetoric drive unauthorized kernel work mid-v86.
6. **PROPOSAL body I.1 rewrite** — delete stable-harmonics slogan and m=1.5 arithmetic; Amendment-only patching is unsafe (Finding 6).

**Finding 13 — Defense against strongest alternative (Alt-A pure carbon now):**  
Pure carbon-now maximizes short-term "atom" screenshots and minimizes intellectual control. The program already burned v84–v85 on ontology thrash (ENT/LQ → Qω → ε retreat). Running carbon on an unfrozen mass definition and uncertified boxes will produce another RESULTS.md that cannot say what was measured. Part 0 is **days and zero GPU**; delaying carbon by that amount is mandatory. Delaying carbon by a full hardened HC-4 + EX-3 campaign is **not** mandatory — that is where I diverge from maximal Amendment 1.

**Finding 14 — Defense against strongest theory alternative (Alt-B second sector first):**  
Q-degeneracy is real (STATE §4, TOE §D). Jumping to kernel surgery without (i) ε closed form, (ii) single-sector window bound (N9), (iii) proof that cloud-carbon is the wrong atom kind, repeats multi-fab product failure modes. Second sector should be a **designed response to a measured ceiling**, not a reset button.

## 3.3 Bottom line on strategy

**The current plan is right in architecture (Part 0 → census → exchange → carbon), wrong in emphasis (late N7, missing N9/D7-lite, C1 theorem inflation, PROPOSAL body not superseded, carbon delayed behind EX-3).**  
I would not abandon it for a different fundamental direction. I would **compress** it to the critical path in Finding 8, **insert D7-lite**, **pull N7 forward**, **restore N9**, and **allow structural carbon after one honest EX-1**, with spectroscopy depth as the parallel theory export rather than the gate.

---

# 4. Required document edits (concrete)

1. Expand GROUNDING §5 to include restored protocols (this file §1) or point to a `v86/N_BATTERY_SPEC.md` extracted from them.  
2. Amendment 1 order: insert N9; move N7 to parallel-after-N1 / before EX-1.  
3. Mark PROPOSAL §§I.1–I.4 and old "Order of work" **superseded** by Amendment 1 + GROUNDING v2.  
4. Add D7-lite as Part 0b or EX-0.  
5. Soften GROUNDING §4 C1 "iff / theorem-shaped" to match FEEDBACK Finding 15 / this Finding 12.

---

# 5. Findings index

| # | Content |
|---|---|
| 1 | §5 adoption strips operational content |
| 2 | Thin §5 is sole spec — must expand |
| 3 | N9 missing from Amendment 1 order |
| 4 | Highest ambiguity: N7, N10, N5 |
| 5 | §0 freeze faithful to FOUNDATIONS |
| 6 | PROPOSAL body still carries killed claims |
| 7 | Order: N7 late; N9 ghost; else OK |
| 8 | Critical path N1–N3–N10 ‖ HC-1–3/6 → N7 → EX-2 → mini EX-1 |
| 9 | Halve-program cuts: EX-3, HC-5, EX-1 high-v, EX-4 cloud, defer N8 |
| 10 | Part 0 before census before exchange: confirmed |
| 11 | Strategy: not different, reweighted |
| 12 | Demote C1 theorem / full HC-4; promote D7-lite + earlier carbon |
| 13 | Reject pure carbon-now |
| 14 | Reject second-sector-first |

---

*End of N-battery review — grok45.*
