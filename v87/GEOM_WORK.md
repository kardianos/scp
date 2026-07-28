# Bell rung B0-G — working file: GEOM seat
(Owned by the geometric grok-4.5 agent. FABLE and GROK may READ this, not write it.)

**Start:** 2026-07-26 ~12:49 CDT  
**Primary stack complete:** ~12:55 CDT  
**Budget:** 4 h wall; VERDICT locked below.  
**Tools:** Maxima 5.46.0, SymPy 1.12. Worksheets under `v87/work/geom/`.  
**Extra brief:** `GEOM_BRIEF.md` — geometric hypothesis + **two hard constraints**.

Status: **PRIMARY STACK COMPLETE** (all seven deliverables; geometric hypothesis
scored as Partial / Strong-with-gap — bound conditional on inner-product form
is full; bound from c² alone is an honest **negative**).

---

## 0. Orientation and labour split

**Hypothesis (GEOM_BRIEF §1), held separate at every step:**

| half | assigned to | GEOM task |
|------|-------------|-----------|
| **Reachability** — how get *above* 2? | O2 field monism: "recursive non-locality by geometric construction" | Identify the *mechanism* that breaks a Bell assumption (not mere holism) |
| **Bound** — why stop at 2√2 not 4? | metric quadratic form / same algebra as c² | Derive \|S\|≤2√2 from a 2-norm / inner-product structure |

**Hard constraints (GEOM_BRIEF §3) — ignoring either invalidates the work:**

1. **Light cone alone ≠ 2√2.** PR boxes are no-signalling and reach S=4. So
   "locality produces a hard stop" cannot mean "no-signalling produces 2√2."
   The survivor form is: the *quadratic form / 2-norm / parallelogram law*,
   not causality *per se*.
2. **Field monism alone ≠ Bell escape.** λ = entire Cauchy-slice field still
   gives local responses A(a,λ), B(b,λ) and CHSH≤2. "Holism" is not enough.
   Recursive non-locality must break *Cauchy-slice freeness* of λ.

**What FABLE and GROK already own (read at start; re-read at ~12:54):**
- Door A explicit models, I_min(2√2)=0.046274 bits, p-family, Tsirelson
  *imported* under pure MD, Doors B/C closed for kernel IVP. Both seats
  **COMPLETE** with mutually agreed VERDICT.
- GEOM owns the **BOUND half** and Constraint-2 sharpening of reachability.
  Ideal join tested in Geo7: their Door-A p=1 model lands on E=−a·b; geometry
  then caps S — but geometry does **not** select p=1.

### Worksheet index

| file | content | status |
|------|---------|--------|
| `work/geom/gm1_tsirelson_geometry.py` | Parallelogram + CS ⇒ \|S\|≤2√2 when E=−a·b | **pass** |
| `work/geom/gm1b_tsirelson.mac` | Maxima: bound at φ=π/2 is 2√2; \|S_can\|=2√2 | **pass** |
| `work/geom/gm2_constraint1_PR.py` | NS permits S=4; MD realises PR | **pass** |
| `work/geom/gm3_constraint2_cauchy.py` | Cauchy-slice monism ⇒ CHSH≤2; \|s(λ)\|=2 | **pass** |
| `work/geom/gm4_operator_S2.py` | S²=4I−[A,A′]⊗[B,B′]; ‖S‖=2√2 on Pauli CHSH | **pass** |
| `work/geom/gm5_why_inner_product.py` | Classical sawtooth ≠ cosine; IP form = selection principle | **pass** |
| `work/geom/gm6_lovasz_sketch.py` | ϑ as vector geometry; almost-quantum caution | **pass** |
| `work/geom/gm7_join_p1_geometry.py` | FABLE p=1 → E=−cosθ; S(p) overshoots 2√2 | **pass** |
| `work/geom/gm8_montecarlo.py` | MC N=2e6: p=1 → \|S\|≈2.8288; p=8 → 3.94 | **pass** |
| `work/geom/gm9_fabric_obstruction.py` | Signature mismatch; derive checklist; scorecard | **pass** |

---

## Findings

(numbered Geo0…; tags [D]/[M]/[P]/[C]; [UNVERIFIED] for unchecked literature)

### Geo0. Pre-selection under O1∧O2, with Constraint 2 applied  **[D]**

O1 keeps **(R)**: definite λ, epistemic ρ.  
O2 makes **(MI)** non-axiomatic (choosers = fabric).  
Kernel locality (finite stencil) keeps **(L)** unless new structure is added.

**Constraint 2 applied to the "one field structure" slogan:**

> Quote GEOM_BRIEF §3: "Given local dynamics, one may always take λ to be the
> entire field configuration on a Cauchy slice in the common past. The responses
> A(a, λ), B(b, λ) are then local functions of it, (R) and (L) both hold, and
> CHSH ≤ 2 follows regardless of how holistic the field is."

So O2 *as holism* does **not** open Door B or relax CHSH. What O2 *does*
open is Door A: ρ(λ|a,b) may fail equality with ρ(λ) because the joint
(source, chooser) configuration is one field with common past. That is the
only reachability route consistent with present kernel locality + O1.

**Agreement with FABLE F0 / GROK G8–G9:** Door A primary; Door B needs a
channel the kernel lacks; Door C closed for IVP. GEOM adds: the geometric
hypothesis's "recursive non-locality" is **not yet a mechanism** — it is a
word. Candidate mechanisms (brief §3): (i) global self-consistency / Door C,
(ii) topological/holonomic over-determination of slice data, (iii) Door-A MD
*derived* from geometry. FABLE/GROK show (iii) is motivated but not derived;
(i) is closed for kernel; (ii) remains **[C]** unexplored — no construction
this rung.

---

### Geo1. Conditional geometric bound: E=−a·b ⇒ |S|≤2√2  **[D, M]**

**Theorem (vector route).** Let a,a′,b,b′ be unit vectors in a real
inner-product space, and E(x,y)=−x·y. Then

    S = E(a,b)+E(a,b′)+E(a′,b)−E(a′,b′)
      = −a·(b+b′) − a′·(b−b′),

hence |S| ≤ ‖b+b′‖ + ‖b−b′‖. Parallelogram: ‖b+b′‖²+‖b−b′‖² = 4.
Cauchy–Schwarz on (1,1)·(‖b+b′‖, ‖b−b′‖):

    ‖b+b′‖ + ‖b−b′‖ ≤ √2 · √4 = 2√2,

with equality when b⊥b′ and a ∥ (b+b′), a′ ∥ (b−b′).

**Machine checks** (`gm1`, `gm1b`):
- Symbolic max of 2(|cos(t/2)|+|sin(t/2)|) = 2√2 at t=π/2 **[D]**.
- Maxima: bound at φ=π/2 equals 2√2; |S_canonical|=2√2 **[D]**.
- 20k random unit pairs in ℝ³: max(‖a+a′‖+‖a−a′‖) ≤ 2√2 (gap ~5×10⁻¹⁰) **[M]**.
- 80-start numeric maximisation of S over planar angles: |S|_max = 2√2 to 1e−15 **[M]**.
- Canonical angles: E=−1/√2 three times, +1/√2 once; S=−2√2 exactly **[D]**.

**Status:** this is the *easy* half of the geometric hypothesis. The brief is
explicit: "the real question is not why 2√2 given inner products — that is
elementary — but WHY ARE THE CORRELATIONS INNER PRODUCTS?"

---

### Geo2. Constraint 1 formalised: light cone / NS ⇏ 2√2  **[D]**

Worksheet `gm2`.

1. **PR box:** E₀₀=E₀₁=E₁₀=+1, E₁₁=−1 ⇒ **S=4**. Standard bit-form
   P(A,B|a,b)=½ if A⊕B=a∧b else 0 is **no-signalling** (marginals = ½).
2. **Door A unrestricted** realises the same correlators with local
   deterministic responses and ρ(λ|a,b) free (re-derives GROK G2 / FABLE F4
   endpoint).
3. **Gap:** 4 − 2√2 ≈ 1.172. Relativistic causality-as-NS stops at 4, not
   Tsirelson.

**Consequence for the hypothesis:** any claim that "the c² light cone produces
the 2√2 ceiling" is **false as stated**. The survivor is weaker: a *definite
quadratic form of the same algebraic class* (parallelogram + CS / C*-norm)
produces 2√2 *once correlations are unit-vector bilinears*. That form lives on
**outcome/setting space**, not on spacetime.

**Lorentzian vs Euclidean (Geo9):** spacetime metric is indefinite; the
Tsirelson 2-norm is positive-definite. Same *class*, different *instance*.
No bridge fabric→Euclidean correlator geometry is constructed this rung.

---

### Geo3. Constraint 2 formalised: monism / holism ⇏ CHSH>2  **[D, M]**

Worksheet `gm3`.

1. **Elementary CHSH fact:** for every ±1 assignment of A,A′,B,B′,
   s = A(B+B′)+A′(B−B′) has **|s|=2** exactly (16/16). Under (R)∧(L)∧(MI),
   |S|=|∫ρ s|≤2.
2. **Cauchy-slice model:** λ = slice data; A depends on (a, λ), B on (b, λ)
   only — holism of the *shared past region* of the slice is allowed and still
   yields |S|≤2 under MI. 200 random local models: max |S| ≤ 2 **[M]**.
3. **Howard 1985** separability vs locality (identifier from brief;
   **[UNVERIFIED]** primary this session) names the distinction: O2 attacks
   *separability of substances*, not the *factorisability of the HV measure*.

**What "recursive non-locality by geometric construction" must mean:**
it must break one of: free Cauchy data, (L), or (MI). Mere "one field" does
not. GEOM therefore **rejects** the naive monism slogan as a Bell solution.

**Reachability mechanisms that *work* (list, not constructions):**
| mechanism | door | status under present kernel |
|-----------|------|------------------------------|
| ρ(λ\|a,b)≠ρ(λ) | A | Motived by O2; **postulated** not derived (FABLE/GROK) |
| Nonlocal response | B | Kernel lacks channel; frame ≠ channel |
| Two-time BC | C | Kernel is IVP; closed for realization |
| Holonomic over-determination of slice | C-ish / topology | **[C]** no construction |

---

### Geo4. Operator route — same bound, same class of object  **[D, M]**

Worksheet `gm4`. For dichotomic tensor-product observables and the CHSH sign
pattern S = AB+AB′+A′B−A′B′:

    S² = 4 I − [A,A′] ⊗ [B,B′]     (sign fixed by the minus on A′B′; **[D]**
                                     machine-checked residual = 0 on Pauli CHSH)

Hence ‖S‖ ≤ √(4 + ‖[A,A′]‖‖[B,B′]‖) ≤ √8 = 2√2. Pauli CHSH operators achieve
op-norm exactly 2√2 **[M]**.

This is not a new number: it is the same 2-norm estimate as Geo1 in operator
language. **Neither route uses the light cone.**

Landau, Phys. Lett. A 120, 54 (1987) — identifier from brief; elementary
route consistent with what is re-derived here; primary not re-read → partial
**[UNVERIFIED]** as citation, mathematics independent.

---

### Geo5. The crux of the bound half: why E = −a·b?  **[D, M, P]**

Worksheet `gm5`.

**Classical local hemisphere** (full MI, uniform λ on the circle):

    E_cl(θ) = −1 + 2θ/π ,   θ=∠(a,b) ∈ [0,π]

Numeric check max |err| ≲ 5×10⁻⁶ on 2×10⁵-point grid **[M]**. At θ=π/4:
E_cl=−½ vs E_qm=−1/√2. At Tsirelson angles |S_cl|=2.

**Gram test:** for E=−cos(θ_i−θ_j), the matrix G_ij=−E is PSD of rank 2
(embedding into ℝ²) **[M]**. Classical sawtooth gives a PSD Gram of *high*
numerical rank (~35 on 36 angles) — it is *not* the low-rank Euclidean
embedding E=−a·b of setting vectors.

**Selection principle IP:**
> Correlations are (minus) inner products of unit vectors assigned to
> measurement settings in a real inner-product space.

- IP ⇒ |S|≤2√2 **[D]** (Geo1).
- IP is **not** implied by the light cone **[D]** (Geo2).
- IP is **not** implied by O2 holism **[D]** (Geo3).
- IP is **not** implied by the Cosserat kernel as measured **[C/D]** (Geo9).
- IP *is* output of QM Born rule, and of FABLE's **postulated** p=1 tilt.

**Verdict on bound half of the hypothesis:** the geometric derivation of
2√2 is real and elementary **conditional on IP**. Deriving IP from fabric
geometry (or from c²) is the missing step — and the place the hypothesis
**fails as a full derivation**.

---

### Geo6. Exclusivity / Lovász as "geometry produces a hard stop"  **[C, partial D]**

Worksheet `gm6`. Cabello–Severini–Winter PRL 112, 040401 (2014) — from brief:
quantum bound = Lovász ϑ of the exclusivity graph, defined via orthonormal
representations in a **real vector space**. That is the sharpest existing
example of geometry-as-hard-stop; it is *event* geometry, not spacetime
geometry.

- C5: ϑ=√5 sits strictly between independence number α=2 and clique cover —
  the classic quantum-like gap from vector geometry **[D standard]**.
- Full CSW↔CHSH identification: not re-derived from scratch this session;
  CHSH already has an elementary complete proof (Geo1), so the ϑ story is
  structural context not load-bearing. **[UNVERIFIED]** as full theorem check.
- **Almost-quantum caution** (Navascués et al. Nat. Commun. 6, 6288 (2015),
  from brief): sets strictly larger than quantum satisfy many principles.
  Landing "near" 2√2 is not unique success. **[UNVERIFIED]** primary details;
  caution adopted.

**Fabric contact:** Cosserat fabric supplies (approx) Lorentzian cone +
preferred-frame defect (THEORY A2). It does **not** supply Hilbert-space
event geometry or Born projectors. Topology (A3/A9) discretizes mode indices,
not correlators.

---

### Geo7. Join with FABLE p-family: geometry caps, does not select  **[D, M]**

Worksheet `gm7`. Re-derived FABLE closed form at p=1:

    B(3/2, 1/2) = π/2 ,
    E_1(θ) = −1 + (2/π)·(π/2)·∫₀^θ sin φ dφ = −cos θ

max |E_1 + cos| < 1e−15 on 25 angles **[M/D]**.

| p | \|S\| at Tsirelson angles | vs 2√2 |
|---|---------------------------|--------|
| 0 | 2.0000 | −0.828 |
| 0.5 | 2.4849 | −0.344 |
| **1** | **2.8284 = 2√2** | **0** |
| 1.5 | 3.0817 | +0.253 |
| 2 | 3.2732 | +0.445 |
| 3 | 3.5355 | +0.707 |
| 10 | 3.9725 | +1.144 |

**Join statement:**
- At p=1, E lands on the unit-vector bilinear form; Geo1 then **derives**
  |S|≤2√2 with equality attained. Bound half succeeds *conditionally*.
- S(p) is smooth and strictly increasing through p=1 (agrees FABLE F4).
  **Nothing geometric about spacetime selects p=1.** The glue between
  reachability (MD) and the bound (IP form) is **postulated**.
- For p>1 the model **overshoots** Tsirelson toward PR — failure mode
  reported, not hidden (brief §5.2).

**Engagement — FABLE F4 quote:**
> "S_max(p) is smooth and strictly increasing through p = 1. Nothing in the
> relaxation mechanism distinguishes the quantum point."

**GEOM:** Agree entirely. Geometry of the *correlator* distinguishes the
quantum point *once E=−a·b is given*; geometry of *spacetime / MD mechanism*
does not put you there.

**Engagement — GROK G7/G10:** same conclusion in the t-family (GOOD-tilt):
no kink at t=1/√2. **Agree.**

---

### Geo8. Monte Carlo of the geometric singlet model  **[M]**

Worksheet `gm8`. Explicit PRNG (numpy PCG64, seed 20260726). Density
ρ∝|λ·b|^p on S² via inverse-transform on u=λ·b; deterministic responses
A=sgn(a·λ), B=−sgn(b·λ). N=2×10⁶ per correlator.

| run | S (MC) | se | target |
|-----|--------|-----|--------|
| p=1 quantum angles | −2.828761 | 0.00100 | −2√2 ≈ −2.828427 (0.33 σ) |
| p=0 classical | −2.000748 | 0.00122 | −2 |
| p=8 overshoot | −3.940214 | 0.00024 | >2√2 confirmed |

Symbolic / closed-form / MC agree within stated SE. Matches FABLE F5 and
GROK g13 within method differences (N, seed).

---

### Geo9. Fabric obstruction to *deriving* the geometric form  **[D, C]**

Worksheet `gm9`. Obstruction chain:

1. **[D]** Local lattice PDE + free Cauchy data + fresh chooser entropy
   ⇒ (R)∧(L)∧(MI) ⇒ |S|≤2 (Geo3).
2. **[D]** NS / light cone ⇒ |S|≤4 only (Geo2).
3. **[D]** Preferred frame 1–5% (THEORY A2) is cone anisotropy, not a
   spacelike channel (agree FABLE F6 / GROK G8).
4. **[D]** Lorentzian spacetime form ≠ Euclidean correlator 2-norm
   (signature mismatch). Shared *class*, not shared *instance*.
5. **[C]** Topology quantizes mode indices (A3, A9), not Born correlators.
6. **[D]** Door A can *postulate* a tilt landing on E=−a·b (p=1); geometry
   then caps. Cap derived; landing not.

**What a derivation would require** (checklist for a future rung):
1. Identify fabric DOF that play the role of Bloch/setting vectors.
2. Derive E(a,b)=−u(a)·u(b) for a unit map u from dynamics.
3. Show the same dynamics forbids PR without importing IC/ML/exclusivity.
4. If Door A: derive ρ(λ|a,b) *and* p=1 selection.
5. If topology: exhibit holonomy over-determining Cauchy data (Constraint 2).
6. Survive fresh-entropy injection at choosers (FABLE F8 algebra-of-generators
   point — **adopted**).

**Cost table (adopted from FABLE/GROK joint result, not re-solved):**

| quantity | value | tag |
|----------|-------|-----|
| I_min(S=2√2) | 0.046274 bits = D_KL((1+1/√2)/2 ‖ 3/4) | [D] FABLE+GROK |
| I_min(S=4) | log₂(4/3) ≈ 0.415 bits | [D] |
| TV budget D at Tsirelson | √2−1 ≈ 0.4142 (tight) | [D] |
| p=1 geometric model I | ≈ 0.202 bits (suboptimal form price) | [M] FABLE |

GEOM does not dispute these numbers; independent re-solution was not the
geometric seat's primary task. Spot-check: MC recovers |S|=2√2 for p=1 (Geo8).

---

### Geo10. Falsifiable consequences (geometric + adopted)  **[C, quantitative]**

1. **In-kernel (sharpest, shared with FABLE F9.1 / GROK G11.3).** Present
   kernel, cone-separated fresh-entropy choosers, no engineered MD:
   **|S| ≤ 2**. Observing |S|>2 falsifies locality of the implementation
   (Gauss-law-class tripwire).

2. **Geometric fingerprint if IP form is ever derived.** If fabric forces
   E=−u(a)·u(b), then:
   - Full singlet correlator (not just CHSH) must track −cosθ, not the
     classical sawtooth −1+2θ/π, across a continuum of angles.
   - Residual: any measured E(θ)+cosθ must vanish within experiment SE;
     classical alternative predicts systematic +0.207 at θ=π/4.

3. **Door-A price list (adopt FABLE/GROK).** Available mutual information I
   after last common-past contact bounds |S|:

   | I (bits) | S ≤ |
   |----------|-----|
   | 0 | 2 |
   | 0.001 | ~2.13 |
   | 0.01 | ~2.40 |
   | 0.02 | ~2.55 |
   | 0.04627 | 2√2 |
   | 0.415 | 4 |

4. **Signature-mismatch fingerprint.** If someone claims "c² alone gives
   Tsirelson," demand a derivation that (a) does not smuggle Euclidean
   correlator structure and (b) explains why PR is excluded. Absence of
   either is a failed claim — this rung's Constraint 1 is the formal
   obstruction.

---

### Geo11. Cross-seat engagement (by name)

#### With FABLE

**Agree:**
- Door A primary under O1∧O2; Doors B/C closed for kernel (F0, F6, F7).
- p=1 ⇒ E=−cosθ exact; family overshoots (F2, F4). **Re-derived** in gm7.
- I_min=0.046274 bits global; TV tight at D=√2−1 (F3). **Adopted.**
- Structure postulated not derived; mixing dynamics oppose MD (F8).
- F9 falsifiables — adopted as Geo10.1 and Geo10.3.

**Quote and respond — FABLE F4:**
> "The measurement-dependence door, by itself, has no Tsirelson wall."

**GEOM:** Correct. The geometric hypothesis hoped the *bound* half would
supply that wall from c². **It does not**, unless IP form is independently
forced. Conditional on IP, the wall is 2√2; fabric does not force IP.

**Quote and respond — FABLE verdict on Tsirelson:**
> "Half the problem is solved; the other half is honestly open, for us and
> for the field."

**GEOM:** The "other half" is exactly GEOM's brief. Result: **conditional
success + honest negative on c²-alone.** Not a full geometric derivation of
Tsirelson from fabric; a clean theorem *once* E=−a·b, and a proof that
spacetime light cone cannot be the source.

**Disagreement?** None on numbers or doors. **Emphasis difference:** FABLE
prices Door A; GEOM proves the geometric bound is real but *orthogonal* to
spacetime causality, and that monism alone is not reachability.

#### With GROK

**Agree:**
- G1–G3 classical facts; G5/G10 I_*=0.046274; G7/G10 family through PR.
- G8 Door B frame≠channel; G9 Door C operationally ⊂ Door A for statistics.
- G12 capacity easy, architecture hard.

**Quote — GROK final joint structural conclusion:**
> "Reaching 2√2 requires **postulated** structure not derived from fabric
> dynamics, and pure MD **does not** select Tsirelson over PR. Bell remains
> a **wall for explanation**."

**GEOM:** Agree. GEOM sharpens "wall for explanation" into two walls:
1. **Reachability wall:** deriving ρ(λ|a,b) / breaking Cauchy freeness.
2. **Bound wall:** deriving IP form (or equivalent event geometry) rather
   than spacetime c².

**On GROK G6 "discretised I≈0.0246":** FABLE F10 says this is impossible by
DPI and that the circle tilt actually saturates 0.046274. **GEOM does not
re-compute GROK's g7 mixture**; adopts FABLE's correction as the consistent
picture. Not GEOM's dispute to re-litigate beyond noting the bound
I≥0.046274 for any model achieving S=2√2 (DPI from 16-strategy optimum).

#### Labour join (the design intent of three seats)

| half | owner | result |
|------|-------|--------|
| Reachability above 2 | FABLE/GROK Door A | Yes, at ≥0.046 bits MI, **postulated** |
| Bound at 2√2 | GEOM | Yes **if** E=−a·b; **not** from c² alone |
| Glue (why E=−a·b from fabric) | **nobody this rung** | Missing; postulated (p=1 or QM import) |

The two halves **do not fully join**. That is a result.

---

## Deliverable checklist (brief §5)

| # | Deliverable | Status |
|---|-------------|--------|
| 1 | Symbolic model + closed-form E + symbolic S-max | **Done** — geometric model E=−a·b (conditional); FABLE p=1 as realizing MD model; gm1/gm1b/gm7 |
| 2 | Bound max\|S\| | **Done** — ≤2√2 under IP **[D]**; =4 under NS alone **[D]**; overshoot of p-family reported |
| 3 | Cost quantified | **Done** — adopt I_min=0.046274 bits; TV D=√2−1; p=1 pays ~0.20 bits |
| 4 | Monte Carlo | **Done** — gm8, seed 20260726, N=2e6 |
| 5 | Derived vs postulated | **Done** — bound derived *conditional on IP*; IP and MD tilt **postulated**; c²-alone **fails** |
| 6 | Falsifiable consequence | **Done** — Geo10 (≥1 quantitative) |
| 7 | Honest negatives | **Done** — Constraints 1–2 as theorems; recursive non-locality not a mechanism; no fabric→IP bridge |

---

## Hypothesis scorecard (GEOM_BRIEF §5 success levels)

| level | description | GEOM outcome |
|-------|-------------|--------------|
| 1 Full | fabric forces IP form + metric gives 2√2 + reach derived | **NOT achieved** |
| 2 Strong | bound from c²-form + reach identified but postulated | **NOT quite** — bound is from *Euclidean correlator* form, not c²/light cone; reach = Door A postulated |
| 3 Partial | quadratic form gives *a* bound at wrong value | N/A — conditional bound is *right* value |
| 4 Negative valuable | proof c² alone cannot produce 2√2 | **ACHIEVED [D]** (Constraint 1 + signature mismatch) |

**Reported level: Partial / Strong-with-gap**, with a clean negative on the
strongest reading of "c² produces Tsirelson."

---

## FINAL VERDICT (GEOM seat, locked ~12:55 CDT)

| Question | Answer | Tag |
|----------|--------|-----|
| Which Bell assumption does O1∧O2 relax? | **(MI)** primarily. O1 keeps (R). Kernel keeps (L). Door C closed for IVP. | [D] |
| Does field monism alone give CHSH>2? | **No.** Cauchy-slice λ ⇒ CHSH≤2 (Constraint 2). "Recursive non-locality" is not yet a mechanism. | [D] |
| Does the light cone alone give 2√2? | **No.** NS max is 4 (Constraint 1). Lorentzian ≠ Euclidean correlator 2-norm. | [D] |
| Can geometry give 2√2? | **Yes, conditional on E=−a·b (IP form).** Parallelogram + CS / operator S² identity. Elementary and machine-checked. | [D,M] |
| Can the model exceed 2 and hit ≤2√2? | **Yes**, via Door A (FABLE p=1 / GROK GOOD-tilt). At p=1, geometry *explains the cap*. Cost ≥0.046 bits MI. | [D,M] |
| Is Tsirelson *derived from fabric*? | **No.** IP form and MD tilt are **postulated**. Spacetime c² does not select them. | [D/C] |
| Do the two halves join? | **Not fully.** Reach = Door A (postulated). Bound = IP geometry (conditional). Glue missing. | [D] |
| Falsifiable fingerprint? | In-kernel \|S\|≤2 without engineered MD; E(θ)+cosθ residual test if IP claimed; S_max(I) price list. | [C/D] |
| **Bottom line** | The geometric hypothesis **survives only in weakened form**: 2√2 is a theorem about **Euclidean / operator 2-norms on correlators**, not about the spacetime light cone. Field monism supplies **motivation** for Door A, not an automatic Bell escape. A deterministic sub-quantum field theory **can** sit at 2√2 (explicit MD models), but fabric dynamics as measured **derive neither** the required measurement dependence nor the inner-product correlator form. **Bell remains a wall for explanation; geometry explains the *shape* of the ceiling once you are under it, not how the fabric gets there or why it may not climb to 4.** | |

---

## Work log

| time (CDT) | action |
|------------|--------|
| 12:49 | Read BELL_BRIEF, COLLAB, GEOM_BRIEF; THEORY A1/A2/A9; FABLE/GROK complete stacks. |
| 12:49–12:53 | Write gm1–gm9 worksheets; formalise Constraints 1–2; geometric Tsirelson; operator S²; join p=1. |
| 12:53–12:55 | Fix abs(S) checks; Maxima gm1b; MC gm8 pass; full findings + VERDICT. |
| 12:55 | **Primary stack complete.** All 7 deliverables landed. |

---

## Appendix A — explicit geometric model (bound half, conditional)

**Assumption IP [P if from fabric; D if from QM/p=1 landing]:**  
settings map to unit vectors a,b∈ℝⁿ; E(a,b)=−a·b.

**CHSH:** S=−a·(b+b′)−a′·(b−b′) ⇒ |S|≤‖b+b′‖+‖b−b′‖≤2√2 **[D]**.

**Realizing Door-A model (FABLE; re-checked Geo7/Geo8):**
- λ∈S²; A=sgn(a·λ), B=−sgn(b·λ) **(R)+(L)**;
- ρ(λ|b)=|λ·b|/(2π)  (p=1) **(MI) relaxed**;
- E=−a·b exactly; |S|_max=2√2; I≈0.20 bits (not MI-optimal).

**Event source:** each run draws uncontrolled fabric DOF λ from a
setting-dependent distribution; frequencies are epistemic (O1). What
"produces" the tilt is **not derived** from Cosserat dynamics.

**What geometry contributes:** the cap at 2√2 once IP holds.  
**What geometry does not contribute:** selection of IP, exclusion of PR
from light-cone structure alone, or Cauchy-slice violation.
