# Bell rung B0 — working file: GROK seat
(Owned by the grok-4.5 agent. Fable may READ this but must not write it.)

**Start:** 2026-07-26 12:22 CDT  
**Budget:** 4 h wall; VERDICT by ~15:22 even if unfinished.  
**Tools:** Maxima 5.46.0, SymPy 1.12. Worksheets under `v87/work/grok/`.

Status: **PRIMARY STACK COMPLETE** (12:45 CDT). All 7 brief deliverables landed.
Door A fully quantified (I_*=0.046274 bits agreed with FABLE); Tsirelson half
failed for pure MD; crux = postulated not derived; Doors B/C closed for kernel.
Continuing light engagement / polish within budget.

---

## 0. Orientation

Bell's theorem is correct: **(R) ∧ (L) ∧ (MI) ⇒ |S| ≤ 2**. The question is which
assumption O1∧O2 force us to relax, what explicit model lives in that relaxation,
what |S| it yields (and why ≤ 2√2 not 4), and whether the structure is *derived*
from fabric dynamics or only *postulated*.

| Tag | Meaning |
|-----|---------|
| **[D]** | Derived (algebra or theorem from stated assumptions) |
| **[M]** | Measured / numerically verified to stated precision |
| **[P]** | Postulated (structure put in by hand) |
| **[C]** | Conjecture |
| **[UNVERIFIED]** | Literature claim not checked against primary source this session |

**Ontology pre-selection (brief + THEORY A1/A2/A9):**
- O1 → keeps **(R)**; forbids fundamental randomness / branching.
- O2 → **(MI)** is the suspicious axiom (choosers share causal past with source).
- A2 preferred frame 1–5% group-velocity **[M, v86]** → necessary for relativistic Door B, not sufficient (kernel remains local).
- A9: Bell closed by theorem for structural-analog programme unless a relaxation is admitted — this rung makes the relaxation precise.

---

## Worksheet index

| file | content |
|------|---------|
| `work/grok/g1_chsh_bound.py` | Exhaustive s(λ)∈{±2}; CHSH≤2 under R∧L∧MI |
| `work/grok/g2_pr_via_md.py` | Unrestricted MD ⇒ PR box S=4 |
| `work/grok/g3_hemisphere_classical.py` | Classical E=−1+2θ/π; \|S\|=2 at Tsirelson angles |
| `work/grok/g4_doorA_chsh_md.py` | Discrete two-point MD model; S=2√2; MC; I≈0.556 bits (suboptimal) |
| `work/grok/g5_min_MI_chsh.py` | GOOD-tilt + numerical min I → **I=0.046274 bits** at S=2√2 |
| `work/grok/g6_maxima_chsh.mac` | Maxima: s∈{±2}, S_qm=−2√2, S_cl=−2 |
| `work/grok/g7_hall_continuous.py` | Continuous circle MD tilt; E=−cos θ; MC; I |
| `work/grok/g8_doorB_nonlocal.py` | Nonlocal toy model; S=2√2 and S=4; cost = foliation |
| `work/grok/g9_doorC_and_tsirelson.py` | Door C ⊂ Door A operationally; Tsirelson not fabric-derived |
| `work/grok/g10_I_closedform_and_family.py` | Closed-form I; family t: S=2→2√2→4 |
| `work/grok/g11_singlet_MI_lower.py` | Continuum planar MI estimates (achievable) |
| `work/grok/g12_fabric_MD_obstruction.py` | O2⇒MD capacity easy; structure hard |
| `work/grok/g13_fable_pfamily_verify.py` | Independent verify FABLE p-family + MI |
| `work/grok/g14_I_star_maxima.mac` | Maxima: I_*=0.04627384685…, D=√2−1 |
| `work/grok/g15_Imin_formula.py` | I_min(S)=D_KL((1+S/4)/2‖3/4) specials |

---

## Findings

### G1. CHSH bound under (R)∧(L)∧(MI)  **[D]**

For fixed λ, A,A′,B,B′∈{±1}:

    s(λ) = A(B+B′) + A′(B−B′),

and exactly one of |B+B′|, |B−B′| equals 2. Hence **|s(λ)|=2**, so
|S|=|∫ρ s|≤2 for any probability ρ independent of settings.

- Worksheet: `g1_chsh_bound.py` (16/16 assignments); Maxima `g6` unique s ∈ {−2,+2}.
- **Consequence:** any |S|>2 must drop at least one of (R), (L), (MI).

### G2. Unrestricted measurement dependence reaches |S|=4  **[D]**

With (R)∧(L) but free ρ(λ|a,b), the PR box is realisable: E00=E01=E10=+1,
E11=−1 ⇒ S=4. Construction: condition λ on the setting pair; responses stay
local functions of (own setting, λ).

- Worksheet: `g2_pr_via_md.py`.
- **Consequence:** Door A without an extra principle **overshoots to 4**.
  Stopping at 2√2 is a separate half of the problem.

### G3. Classical local hemisphere model (full MI)  **[D,M]**

    λ ~ Unif[0,2π),   A=sgn cos(λ−a),   B=−sgn cos(λ−b).

Closed form: **E(θ) = −1 + 2θ/π** for θ=∠(a,b)∈[0,π].

At Tsirelson angles (a,a′,b,b′)=(0, π/2, π/4, −π/4):

| pair | θ | E_cl | E_qm=−cos θ |
|------|---|------|-------------|
| ab   | π/4 | −1/2 | −1/√2 |
| ab′  | π/4 | −1/2 | −1/√2 |
| a′b  | π/4 | −1/2 | −1/√2 |
| a′b′ | 3π/4 | +1/2 | +1/√2 |

**S_cl = −2**, **S_qm = −2√2**. Numeric integral of sign product matches closed form to ≲10⁻⁵ (grid).

- Worksheets: `g3_hemisphere_classical.py`, Maxima `g6`.

### G4. Door A — discrete local det. model with MD, exact Tsirelson  **[D,M]**

**State space:** λ=(A₀,A₁,B₀,B₁)∈{±1}⁴ (16 deterministic local strategies).  
**Responses (R)+(L):** A(a,λ)=A_a, B(b,λ)=B_b.  
**Event source:** for each run, settings S=(a,b) drawn (here uniform on 4 pairs);
then λ ~ ρ(·|S); outcomes read from the table. Between runs, λ is **not** held
fixed independently of settings — that is the MI relaxation. Epistemic p is
over uncontrolled fabric degrees of freedom that correlate with chooser state.

**Two-point construction (pedagogical, suboptimal MI):** for target correlator t,
mass (1+t)/2 on one λ with product +1 and (1−t)/2 on one with product −1.

With t=±1/√2 in CHSH sign pattern: **S = 2√2 exactly**.  
I(Λ:S) = **0.556483 bits** (H(Λ)=1.157, H(Λ|S)=0.601).

**Monte Carlo** (N=200000, seed 42):

| E | analytic | MC ± SE |
|---|----------|---------|
| E00 | +0.70711 | +0.70602 ± 0.00316 |
| E01 | +0.70711 | +0.70522 ± 0.00317 |
| E10 | +0.70711 | +0.70622 ± 0.00317 |
| E11 | −0.70711 | −0.71048 ± 0.00315 |
| **S** | **2.82843** | **2.82794** |

- Worksheet: `g4_doorA_chsh_md.py`.

### G5. Minimal MI for four-setting CHSH at S=2√2 in the 16-strategy model  **[D,M]**

**GOOD set:** the 8 strategies with s(λ)=+2. Uniform on GOOD gives E=(½,½,½,−½),
**S=2** (classical saturation).

**Tilt construction:** for target T, reweight  
w_i ∝ 1_GOOD(i)·(1 + α·prod_i) with α=(T−m)/(1−T m), m=½.  
At |T|=1/√2: **S=2√2**, **I(Λ:S)=0.04627384685 bits**.

Numerical SLSQP minimisation over all conditionals on 16 atoms (12 random
seeds): **best I = 0.046274 bits** in every successful run — matches the tilt
construction. **FABLE F3 (Blahut–Arimoto) upgrades this from “apparent” to
global:** I_min(S)=D_KL((1+S/4)/2 ‖ 3/4) for all S∈[2,4], so at Tsirelson
**I_min = 0.046274 bits is exact** for 4-setting CHSH over the 16-strategy
model. GROK constructed the achieving family; FABLE certified optimality.
Cross-seat agreement **[D]**.

**Closed form [D]:** with η=(1+1/√2)/2 and μ₊=3/4 (fraction of GOOD with
product +1 on a CHSH correlator),

#if
I_* = D_{KL}\!\left(\eta \,\middle\|\, \tfrac{3}{4}\right)
    = \eta\log_2\frac{\eta}{3/4} + (1-\eta)\log_2\frac{1-\eta}{1/4}
    = 0.04627384685\ldots \text{ bits.}
#endif

(SymPy exact expression in `g10`; equals the numeric I to 1e−16.)

**Compare to brief / Hall:** brief cites ≈0.0663 bits (~1/15) for the *full
singlet* (continuum of settings). Our **0.0463 bits is CHSH-only** (4 setting
pairs) — smaller, as expected. Barrett–Gisin PRL 106, 100406 (2011) and Hall
PRL 105, 250404 (2010) / PRA 84, 022102 (2011) are the literature anchors; the
0.0663 figure is **[UNVERIFIED]** as a primary-source recomputation here, but
is consistent in magnitude with our CHSH result.

- Worksheets: `g5_min_MI_chsh.py`, `g10_I_closedform_and_family.py`.

### G6. Continuous Door-A model (Hall-style circle tilt)  **[D,M]**

Same responses as G3; ρ(φ|θ) is the **level-set exponential tilt** enforcing
⟨σ⟩=−cos θ (σ=AB). Closed form: masses p±=(1∓cos θ)/2 on {σ=±1}, uniform
inside each level set (μ₊=θ/π under classical geometry).

At Tsirelson angles: **|S|=2√2 exactly**.  
Per-setting D_KL(ρ‖uniform)=**0.046274 bits** (same number!).  
Discretised I(Λ:S) to mixture marginal ≈ **0.0246 bits** (four conditionals
are not independent of each other — shared circle, related angles).

MC (N≈3e5): S_MC=−2.825 vs analytic −2.828.

- Worksheet: `g7_hall_continuous.py`.

### G7. One-parameter family: 2 → 2√2 → 4  **[D]**  ★ Tsirelson half

GOOD-tilt with correlator magnitude t ∈ [½, 1]:

| t | S=4t | I(Λ:S) bits |
|---|------|-------------|
| 0.50 | 2.000 | 0.000 |
| 0.60 | 2.400 | 0.010 |
| 1/√2 | **2.828** | **0.046** |
| 0.80 | 3.200 | 0.105 |
| 0.90 | 3.600 | 0.208 |
| 0.99 | 3.960 | 0.378 |
| 1.00 | **4.000** | **0.415** |

**There is no critical point, fixed point, or kink at t=1/√2.** Door A with a
smooth MD budget produces a smooth S(I) curve through Tsirelson toward PR.

**Therefore:** measurement dependence alone **does not** explain why nature
stops at 2√2. An independent principle is required (information causality,
macroscopic locality, exclusivity, or emergent Hilbert-space structure).
None of these is derived from Cosserat fabric dynamics in this rung — they
remain **[P/C]** imports. (Anchors: Pawłowski et al. Nature 461, 1101 (2009);
Navascués & Wunderlich Proc. R. Soc. A 466, 881 (2010); Fritz et al.
Nat. Commun. 4, 2263 (2013) — identifiers from brief.)

- Worksheet: `g10_I_closedform_and_family.py`.

### G8. Door B — nonlocal element  **[D,M]**

**Toy model NB1:** λ=(α,u), α∈{±1} fair, u~Unif[0,1];  
A=α (local); **B=α if u<(1+E(a,b))/2 else −α** (depends on distant setting a).  
(MI) holds; (R) holds; **(L) fails**. MC recovers S≈2√2 for quantum targets
and **S=4** for PR targets.

**Cost:** preferred foliation for the nonlocal influence (who depends on whose
setting). No-signalling can still hold (marginals fair in NB1).

**SCP kernel contact:**
- Kernel is a **local** lattice field theory (finite stencil, link-covariant
  derivatives, local potential). **No Door-B channel exists in the present
  dynamics.** **[M/context from kernel structure]**
- Preferred frame at 1–5% group velocity **[M, THEORY A2]** is a *dispersion /
  anisotropy* fact, not a spacelike influence channel. **Necessary for a
  relativistic Door-B story, nowhere near sufficient.**

**Agreement with FABLE F0** (quote): FABLE writes *“Door B requires a
genuinely nonlocal guidance structure … a preferred frame is necessary for
that but nowhere near sufficient. The kernel has the frame and lacks the
channel.”* — **I agree entirely.** See § Cross-seat.

- Worksheet: `g8_doorB_nonlocal.py`.

### G9. Door C — retrocausality / all-at-once  **[D]**

Operational content of retrocausal correlation of λ with future settings **is
measurement dependence** (Door A statistics). Wharton & Argaman RMP 92, 021002
(2020) — identifier from brief.

**O1 readings:**
1. **Dynamical / Cauchy (kernel-native):** leapfrog IVP; future BC do not
   constrain source λ. **Door C closed for the present realization.**
2. **Block universe (one definite spacetime history):** all-at-once constraints
   allowed in principle; still requires a *principle* selecting quantum not PR.

FABLE F0: *“Door C is closed for this realization — not for field monism in
general.”* — **Agree on kernel-IVP closure.** I add: block reading of O1 keeps
Door C philosophically open but does not help the Tsirelson half.

- Worksheet: `g9_doorC_and_tsirelson.py`.

### G10. The crux — derive vs postulate  **[C/D]**

| Structure required | Status under O1∧O2 + present kernel |
|--------------------|--------------------------------------|
| Local det. responses A(a,λ), B(b,λ) | Compatible with O1 **[D]** |
| ρ(λ\|a,b) ≠ ρ(λ) with I≈0.046 bits (CHSH) or ~0.066 bits (singlet, Hall) | **Motivated** by O2 (shared causal past) but **not derived** from fabric PDE **[P]** |
| Bound forbidding S>2√2 | **Not derived**; pure MD allows S=4 **[D]** |
| Nonlocal guidance (Door B) | **Not present** in kernel; would be new physics **[P if added]** |
| Two-time BC (Door C) | **Not present** in kernel IVP |

**Verdict on crux:** With present tools we can **postulate** a Hall-type MD
model that reproduces CHSH = 2√2 at ~0.046 bits of mutual information, and we
can **prove** that unrestricted MD gives 4. We **cannot** derive the required
ρ(λ|a,b) from fabric dynamics, nor derive Tsirelson.  
A postulated correlation between λ and settings is a **parametrisation of the
answer, not an explanation.**

What would count as derivation:
1. From fabric causal structure, compute the joint distribution of (source
   degrees of freedom, chooser degrees of freedom) and show I(Λ:S) lands at
   the Hall value without free functions; **and**
2. Show the same dynamics forbids PR (S=4), e.g. by an emergent information-
   causality inequality or an emergent operator algebra.

Neither is available today. Topology is the common ingredient of the other
v86 walls (NEXT_PROGRAM); it is **not obviously** the Bell ingredient — Bell
is about **setting–source correlation structure and composition of instruments**,
not about winding numbers.

### G11. Falsifiable consequences  **[C, quantitative]**

At least one fingerprint per open door:

1. **Door A fingerprint — residual measurement dependence.**  
   If I(Λ:S) ≥ I_* ≈ 0.046 bits for CHSH-optimal experiments, then in a
   sufficiently complete model of the chooser+source fabric there exists a
   classical predictor of settings from source-side variables (or vice versa)
   with that mutual information.  
   **Operational test (in-fabric):** build two chooser subsystems whose
   internal state is scrambled by independent high-entropy baths for time
   τ_scram; measure CHSH vs τ_scram. Prediction: **S(τ) → 2** as the baths
   destroy setting–source correlation, with a characteristic scale set by
   fabric correlation length / bath coupling.  
   Quantitative target: to keep |S|≥2.7 requires I≳0.04 bits on our family
   curve (`g10`); full Tsirelson needs ≥0.046 bits in the 16-strategy model.

2. **Door B fingerprint — preferred-frame causal order.**  
   Nonlocal influence respects a foliation. Signature: CHSH value or
   signalling residuals depending on the absolute orientation of the
   Alice–Bob axis relative to the lattice rest frame, at the same 1–5%
   level as the group-velocity anomaly **[M, A2]** — *if* Door B is realised
   by promoting that frame to a Bohm-like slice. (Current kernel: no such
   nonlocal term, so this tests a *modified* kernel, not the present one.)

3. **Negative fingerprint for pure local + MI fabric:**  
   If the kernel is run as-is with spacelike-separated “chooser” patches and
   a shared entangled-analog seed **without** engineered MD, CHSH must
   respect |S|≤2 (G1). Any observed |S|>2 in pure local evolution would
   falsify the locality of the implementation (bug or hidden channel).

---

### G12. FABLE p-family (independent confirmation)  **[D,M]**

Re-derived and checked in `g13_fable_pfamily_verify.py`:

    ρ_p(λ|b) = ((p+1)/(4π)) |λ·b|^p ,   A=sgn(a·λ), B=−sgn(b·λ)

    E_p(θ) = −1 + ((p+1)/π) B((p+2)/2, 1/2) ∫_0^θ sin^p φ dφ

- **p=1:** B(3/2,1/2)=π/2 ⇒ E=−cos θ **exact** **[D]**
- MC (N=2e5) agrees to ≲0.002
- I(Λ:b)|_{p=1,φ=π/2} ≈ 0.202 bits (suboptimal vs I_*=0.046)

Citation Degorre–Laplante–Roland PRA 72, 062314 (2005): **[UNVERIFIED]** as
primary source this session; the *mathematics* is re-derived, not citation-dependent.

### G13. TV cost saturates at Tsirelson  **[D]**

FABLE's bound S ≤ 2+2D (D=Σ_xy TV(ρ_xy,m)) is **tight**. GOOD-tilt achieves
D=√2−1 and S=2√2 simultaneously — equality in the bound. So the *minimal*
total-variation MD budget for Tsirelson CHSH is exactly √2−1.

### G14. Fabric derivation obstruction  **[D/C]**

O2 ⇒ common causal past of source and choosers ⇒ generic I(Λ:S)>0
(capacity: ~0.046 bits needs only eps≈0.25 soft bias on a common-past bit —
`g12`). Locality ⇒ only common-past-mediated MD. The **structure** that hits
quantum correlators and not PR is **not** supplied by the local Cosserat PDE
as currently formulated. Capacity is easy; architecture is the wall.

---

## Deliverable checklist (brief §5)

| # | Deliverable | Status |
|---|-------------|--------|
| 1 | Symbolic model + closed-form E + symbolic/exact S | **Done** G3–G6, Maxima+SymPy |
| 2 | Bound max\|S\| | **Done:** ≤2 under MI; =2√2 when quantum targets imposed; =4 unrestricted MD |
| 3 | Cost quantified | **Done:** I_*=0.046274 bits (CHSH); family to 0.415 bits at PR; Hall ~0.066 singlet noted |
| 4 | Monte Carlo | **Done** g4, g7, g8 |
| 5 | Derive vs postulate | **Done:** structure postulated; not derived from fabric |
| 6 | Falsifiable consequence | **Done** G11 (≥1, quantitative) |
| 7 | Honest negatives | **Done:** Door B/C closed for kernel; MD fails Tsirelson half alone |

---

## Cross-seat notes (FABLE)

### Agreement zone

**FABLE F0 (quoted):**  
> “O1 ∧ O2 admit Door A naturally, admit Door B only by adding structure the
> kernel does not have, and do not admit Door C at all for the present dynamics.”

**GROK:** **Agree** on Door A primary; Door B (frame ≠ channel); Door C closed
for leapfrog IVP. Block-universe O1 keeps Door C philosophically open but does
not help Tsirelson.

**FABLE F2 p-family** ρ_p(λ|b)∝|λ·b|^p on S² with hemisphere responses:
- p=0 → E=−1+2θ/π, S=2
- **p=1 → E=−cos θ exact singlet**, S=2√2
- p→∞ → S→4

**GROK G13 independently confirms:** algebraic identity B(3/2,1/2)=π/2 gives
E_1=−cos θ exactly **[D]**; MC on the sphere matches closed form to ~10⁻³
(N=2e5); S(p=0)=−2, S(p=1)=−2√2 at Tsirelson angles; I(Λ:b) at p=1, φ=90°
≈ **0.202 bits** (GROK grid) vs FABLE **0.2018 bits** — **agreement**.

**FABLE F4:** S_max(p) smooth through 2√2 — Tsirelson imported. **Agree**
(GROK G7/G10 family t is the same conclusion in a different parametrisation).

**FABLE F3 TV bound:** S ≤ 2 + 2D with D=Σ TV(ρ_xy, m), tight at
D=(S−2)/2. **Agree and strengthen:** GROK GOOD-tilt **saturates** it:
D=√2−1≈0.414214 at S=2√2 exactly (see computation 12:43).

### I_min at Tsirelson — dispute **resolved** (FABLE update)

FABLE initially reported I_min≈0.0638 (local convex optimum). GROK claimed
≤0.046274 via GOOD-tilt feasibility. **FABLE F3 update confirms GROK:**
Blahut–Arimoto gives

    I_min(S) = D_KL((1+S/4)/2 ‖ 3/4)

hence **I_min(2√2)=0.046274 bits** exactly, and I_min(4)=log₂(4/3)≈0.415 bits.
FABLE further states the tilt family is the **global optimum at every S**.

| source | quantity | value | status |
|--------|----------|-------|--------|
| **GROK G5/G10 + FABLE F3 BA** | I_min(S=2√2) | **0.046274 bits** | **agreed [D]** |
| Closed form | D_KL(η‖3/4), η=(1+1/√2)/2 | 0.04627384685… | Maxima g14 |
| TV budget | D=√2−1 | 0.414214 | tight, both seats |
| FABLE p=1 geometric | I(Λ:b) | 0.2018 bits | suboptimal (form price) |
| Hall full singlet | I | ~0.0663 bits | **[UNVERIFIED]** primary; ≥ CHSH-only OK |

**Why Hall ~0.066 > CHSH-only 0.046:** continuum of settings vs four pairs.

### Labour / convergence (updated 12:46)

FABLE owns the continuous p-family (best geometric singlet model) and the
BA global-optimality proof for I_min(S). GROK owns the GOOD-tilt construction
that first hit 0.046274 and the TV saturation check. FABLE F10 credits the
correction — collaboration as designed.

**FABLE F6–F8:** Door B closed (frame≠channel); Door C closed for IVP;
structure postulated; mixing dynamics work *against* MD derivation. **Agree
in full.** FABLE's point that the postulate must correlate λ with the *entire
algebra of setting generators* (not one fixed chooser) is sharper than
GROK G12 — adopt it: superdeterministic Door A needs a conspiracy under
arbitrary fresh-entropy injection, which the measured fabric (Adler averaging,
radiation of unclosed flow) does not provide **[C, FABLE-anchored]**.

**FABLE F9 falsifiables:** compatible with GROK G11; their S_max(I) inversion
of the price list is a clean quantitative form. Adopt as shared prediction:

| available I (bits) | S ≤ |
|--------------------|-----|
| 0.04627 | 2√2 |
| 0.02 | ~2.55 |
| 0.01 | ~2.40 |
| 0.001 | ~2.13 |
| 0 | 2 |

### Final joint structural conclusion

A deterministic sub-quantum field theory obeying O1∧O2 **can** produce
CHSH∈(2,2√2] by relaxing (MI) with ≥0.046 bits of source–setting mutual
information (or by adding nonlocal structure, Door B). The **present SCP
kernel does neither**: it is local+IVP, so |S|≤2 for cone-separated
fresh-entropy choosers. Reaching 2√2 requires **postulated** structure not
derived from fabric dynamics, and pure MD **does not** select Tsirelson over
PR. Bell remains a **wall for explanation**.

---

## FINAL VERDICT (locked 12:46 CDT)

| Question | Answer | Tag |
|----------|--------|-----|
| Which door does O1∧O2 motivate? | **Door A (MI)** primary; Door B only with new nonlocal structure; Door C closed for kernel IVP | [D] |
| Can \|S\|>2 under O1∧O2? | **Yes**, via MD (Door A) or nonlocal responses (Door B) | [D] |
| Explicit model with S=2√2? | **Yes** — GOOD-tilt (I-optimal); FABLE p=1 sphere (exact singlet, I=0.20); circle tilt | [D,M] |
| Cost? | **I_min=0.046274 bits** (global, CHSH); TV D=√2−1 tight; full singlet ~0.066 **[UNVERIFIED]** Hall | [D,M] |
| Why not 4? | **Not answered by Door A alone.** I_min(4)=log₂(4/3)≈0.415 bits. Tsirelson imported | [D] |
| Derived from fabric? | **No.** Structure **[P]**; O2 motivates existence of *some* MD, not the quantum tilt; kernel mixing works against it | [C/D] |
| Falsifiable fingerprint? | In-kernel: \|S|≤2 for fresh-entropy cone-separated choosers. Door-A completion: S_max(I) price list | [C/D] |
| **Bottom line** | O1∧O2 **admit** CHSH above 2 by relaxing (MI), with an explicit local deterministic model at **0.046 bits/run**. They do **not** derive that model from fabric dynamics, nor the bound 2√2. The present kernel realises none of the required structure: **Bell is a wall for explanation.** | |

---

## Work log

| time (CDT) | action |
|------------|--------|
| 12:22 | Read BELL_BRIEF, COLLAB, THEORY A1/A2/A9; start file. |
| 12:23–12:30 | G1–G4 worksheets; CHSH bound, PR-via-MD, hemisphere, discrete MD+MC. |
| 12:27–12:32 | g5 min-I (I_*=0.046274); Maxima g6; g7 continuous; g8 Door B; g9 Door C. |
| 12:30–12:33 | g10 family t: 2→2√2→4; closed-form I=D_KL(η‖3/4). |
| 12:33–12:40 | g11 continuum MI; g12 fabric obstruction; read FABLE worksheets. |
| 12:40–12:45 | g13 verify FABLE p-family (p=1 exact singlet); TV saturation; I_min dispute. |
| 12:45 | FABLE confirms I_*=0.046274 via Blahut–Arimoto; dispute resolved. Maxima g14. |
| 12:46 | g15 I_min(S) formula; FABLE F6–F10 engagement; FINAL VERDICT locked. |
| 12:46 | **Primary stack complete.** All 7 deliverables landed. |

---

## Appendix A — explicit GOOD-tilt formulae (Door A, primary model)

**λ-space:** 16 points, of which GOOD = {λ : A₀B₀+A₀B₁+A₁B₀−A₁B₁ = +2} (8 pts).

**Prior base (no MD):** ρ₀(λ)=1/8 on GOOD, 0 else. Then  
E₀₀=E₀₁=E₁₀=+½, E₁₁=−½, **S=2**.

**With MD, target T_ab ∈ {+c,+c,+c,−c}, c=1/√2:**  
α_ab = (T_ab − m_ab)/(1 − T_ab m_ab), m_ab=⟨A_a B_b⟩_{ρ₀}=±½,  
ρ(λ|a,b) ∝ ρ₀(λ)(1 + α_ab A_a(λ) B_b(λ)).

**Result:** E_ab=T_ab exactly, **S=2√2**, **I=D_KL(η‖3/4)** with η=(1+c)/2.

**Event source:** each experimental run draws fabric configuration from
ρ(·|settings); settings may be chosen by fabric subsystems correlated with
the source through the common past (O2). Uncontrolled DOF ⇒ epistemic
frequencies matching the integrals above.

**What is postulated:** the map (a,b)↦α_ab (or equivalently the targets T_ab).  
**What is derived:** once that map is given, E and S and I are fixed algebraically.

---

## Appendix B — what would change the verdict

1. A derivation of information causality (or exclusivity) from fabric instrument
   composition → would close the Tsirelson half under Door A.
2. A measured or derived nonzero I(source:chooser) in a concrete dual-chooser
   simulation → would promote MD from [P] to [M].
3. A kernel modification adding nonlocal guidance on a preferred slice → would
   reopen Door B as programme physics (requires explicit user authorisation for
   kernel change — not done here).
