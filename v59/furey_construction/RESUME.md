# RESUME — v59 Furey ℂ⊗ℍ⊗𝕆 program

*Snapshot 2026-05-24.  Pick-up point for the parameter-reduction / v58-integration effort.
**START HERE for the 2026-05-24 dynamics/couplings session: `integration_v58/00_SUMMARY.md`** —
the capstone indexing docs `01`–`20` + the new Lean, with the honest status of every claim.
Companion docs: `PARAMETER_BUDGET.md`, `RIGOR_AUDIT.md`, `INTEGRATION_PLAN_v58.md`,
`koide_phase_law/`, `integration_v58/`.  All Lean builds clean; axioms = standard trio.*

**Bottom line of the 2026-05-24 session (`00_SUMMARY.md`):** the octonion framework + the v58/OFE
dynamics fix the DISCRETE data and the *meaning* of the couplings, but NOT their continuous VALUES
(`t²=½`, `φ=2/9`, `α`) — every derivation route closed with a computed reason (kinematics/dynamics
divide).  The values are inputs/residuals (the flavor problem).  New Lean (axiom-clean): `OctoHalf`
(the `½` = L-grade complex root-product), `ChiralPhaseWindow` (M4), `MaximalMixingKoide`,
`HiggsVevReframe`.

## 1. What this is

A program identifying Standard-Model parameters with structural quantities of the Furey color
algebra `Cl(7)_even ≅ ℂ⊗𝕆` (dims `14,21,16,28,35,63`; Koide `Q=dimG₂/dimSpin7`; grade split
`L=Λ²⊕Λ⁶` / `F=Λ⁴`), now being put on a **dynamical** footing via the **v58 multivector field
equation** `⟨DΩ+λΩ²+μ⟨Ω,Ω⟩⟩_{0,2} = f_g(ρ)J_ρ + f_em(ρ)J_χ`, `ρ_M=½(MM̃−v²)`.

## 2. Current status — the one distinction that matters

- **Theorem-grade (Lean-checked):** the *algebraic skeleton* — Koide `Q=14/21`, grade structure
  (lepton=L, `J_c∈Λ²`, `BladeSquareSign`, `LeptonRealityForcing`), color `su(3)` + classification
  (`ColorSU3`), `5=h∨(Spin7)`, `sin²θ_W=2/9` given `(5,2)`, the sedenion `S₃` automorphism
  (the "/3"), and the v_Higgs reframing (`HiggsVevReframe`).
- **Empirical / conjecture-dependent:** every *physical prediction* (masses, `α`, `v_Higgs`,
  `g_W²=5√α`, …) is a match of a structural number to data, several hanging on underived
  ansätze.  **The "≈1 free parameter" headline is conditional on those conjectures.**

## 3. Free variables (genuine inputs)

| input | type | note |
|---|---|---|
| **`a_lepton`** | dimensionful scale | **irreducible** — sets the unit; cannot be "derived," only its dimensionless ratios |
| **`α`** | dimensionless | conjecturally structural (`α(M_Z)=25/(324π²)`); if rejected, free.  *The* legitimate derivation target |
| quark scales `a_u,a_d`; quark phases `φ_u,φ_d` | — | free (phases shown ≠ `Q_q/3`) |
| CKM (4), neutrinos (≥7), `θ_QCD`, `α_s` | — | not engaged |

**Count:** lepton+EW+Higgs block ≈ **1 dimensionful input (`a_lepton`) + `α`** (optimistically 1 if
the `α` conjecture holds).  Full SM: still **~10–12 free** (quark flavour / CKM / neutrino / strong
sectors untouched).

## 4. Dependent variables (turned structural), with status

`[thm]`=Lean identity · `[emp≈X]`=empirical match · `[conj]`=undrived ansatz

| quantity | relation | status |
|---|---|---|
| lepton amplitude `t` | `t²=1/2 ⟺ Q=14/21` (Koide) | `[thm]`+`[emp 10⁻⁵]` |
| lepton phase `φ` | `φ=Q/3`; the `/3`=sedenion `S₃` | `[emp 10⁻⁵]`; `/3` structural, magnitude residual |
| **`v_Higgs`** | `=dim(L)²a_l²` ⟺ `Σ√m/√v=3/28` ⟺ Yukawa sum rule `Σ√y=(3/28)2^{1/4}` | `[thm-equiv]`+`[emp 0.03%]`; mechanism = geometric dilution over `L≅so(8)` (uniformity proven); **one overlap integral open** |
| `sin²θ_W` | `=2/9` (Pati-Salam, `cBL=2` pinned) | `[thm given (5,2)]` |
| `m_W,m_Z` | from `v_Higgs,α,sin²θ_W` | `[emp 0.02–0.04%]` |
| gauge `5` | `=h∨(Spin7)` | `[thm]` |
| quark Koide `Q_d=11/15,Q_u=23/27` | `t²=1−dimG₂/D_N` | `[emp 0.3%, soft/RG-dependent; selection rule UNDRIVEN]` |
| `α(M_Z)` | `=25/(324π²)` | `[conj]` — derived *given* `g_W²=5√α` |

**Load-bearing conjectures still undrived:** `g_W²=5√α` (the `√α` form), `v_Higgs` overlap
integral, `α(0)/α(M_Z)`, `G_e=(21/16)α²¹`, the quark selection rule.

## 5. Active work: v58 integration (replace conjectures with physical derivations)

Linchpin: **identify v58's vacuum norm `v` (in `ρ_M=½(MM̃−v²)`) with `v_Higgs`.**  Then the
dim-map says which Clifford grade carries which physics.

- **#1 v_Higgs — BRIDGE DONE (mechanism pinned), two dimensionful residuals remain.**  Clean form
  locked (`3/28`); reframed as a lepton **Yukawa sum rule** `Σ√y=(3/28)2^{1/4}` (machine-checked).
  **Cl(7) bridge computed** (`integration_v58/bridge_vhiggs_cl7.py`, `03_higgs_bridge_result.md`,
  Lean `HiggsVevReframe.{frobeniusSq_democratic,bilinear_reading,vector_reading,readings_differ}`):
  the 28 `L`-blades are Frobenius-**orthogonal**, so the prior **"equipartition over 28 directions"
  mechanism is WRONG** (it gives `√v=√28·a≈5.3a`, not 28).  The factor 28 is correct only under
  the **bilinear** reading — the lepton Yukawa is a *matrix* on `L`, `v_Higgs=‖Y‖²_Frobenius =
  dim(L)²·a_l²`, `28²` = **#components** of the `L`-bilinear, natural `L²` norm, `so(8)` democracy.
  **Lepton-specificity derived** (quarks `D=35,63` miss `√v=D·a` by 1.8×, 19× — only the singlet's
  `Y` is wholly `L`-grade).  *Irreducible residuals (no dynamical home in v58 or v59-XiVacuum):*
  **(R1)** why `v=‖Y‖²_F` (= 246 GeV); **(R2)** why each component `=a_l`; + the `28`-dim-`L`-vs-`3`-
  generations consistency.  **Recommend: bank #1 here, move to `α`.**  Docs: `01_higgs_vev.md`,
  `02_higgs_iii_yukawa.md`, `03_higgs_bridge_result.md`.
- **α — SCOPED (2026-05-24, `ALPHA_SCOPING.md`).**  Three *separate* conjectures, not one:
  - **C1 `g_W²=5√α` → α(M_Z)=25/(324π²) (0.032%) — ATTEMPTED, NEGATIVE (`04_gW_sqrt_alpha_result.md`).**
    `5=h∨(Spin7)` **proven**; `sin²θ_W=2/9` **proven**.  But the **`√α` form is NOT a derivable law**:
    given the definitional `e=g_W sinθ_W` it is `α(M_Z)=25/(324π²)` *in disguise* (the line `g_W²=18πα`
    meets `5√α` at a **single point** = the EW scale; they disagree at all other scales — verified).
    The `Ω²` mechanism is **excluded 3 ways**: (1) v58 derives EM as *primary* (`DΩ=J_χ`), not sourced
    by a weak `Ω²`; (2) v58 has **no weak sector** and the steelman needs `λ≈1/h∨=0.2`, 40× outside
    v58's safe band `|λ|≤0.005`; (3) standard `U(1)⊂Spin(7)` embedding gives `g_W²∝α` (**linear**),
    never `√`.  **Honest restatement:** `√α(M_Z)=h∨·sin²θ_W/(4π)` is a *value-match* for `α(M_Z)`
    (0.032%), not a gauge law — so **`α` stays a genuine conjectural input** (the "1-input" headline
    does NOT get `α(M_Z)` for free).  **α derivation target: closed negative for v58/v59 as built.**
  - **C2 `−ln α + 2α = π²/2` → α(0) — mechanism-poor.**  Pure BPST-instanton form `−ln α=π²/2`
    (`g²=16=dimCl31`) is only **1.4%**; the `3.6×10⁻⁵` ("0.004%") match needs the **`+2α` term,
    which `04_alpha_prediction.py` Part 3 REVERSE-ENGINEERED** (fitted, not derived).  Instanton
    ungrounded (`π₃(S⁷)=0`).  *(Fixed the wrong "2.4×10⁻³" gap comment in `AlphaZero.lean`.)*
  - **C3 `G_e=(21/16)α²¹` — numerology** (0.25%, wild exponent, no home).  Flag and drop.
- **#4 quark Koide / selection rule — LOW priority** (soft, RG-dependent).

## 6. Short-term goals  *(rewritten 2026-05-24 after the Higgs/α/Koide session)*

**Status:** #1 banked (Higgs mechanism pinned, 2 residuals); α closed negative (C1/C2/C3 all
mechanism-poor — `ALPHA_SCOPING.md`, `04_gW_sqrt_alpha_result.md`).  **The session's real find:
Koide=2/3 and `v_Higgs=28²a²` are ONE mystery — equipartition / maximal mixing of the order
parameter (`EQUIPARTITION_PRINCIPLE.md`).**

> **DYNAMICS LINE — CLOSED NEGATIVE (2026-05-24, `18_`→`20_…md`).**  The hope: lift v58 to Cl(7)
> {0,2,4,6}, use the new **grade-4 (F) channel** (non-flat on a generic 28-dim L-config, `⟨M²⟩₄≠0`) to
> select the alignment.  **CORRECTED/closed:** the lepton is the **color singlet**, confined to the
> 4-dim **`ℍ`-slice** (`{1,e₀₁,e₀₂,e₁₂}`, closed `≅ℍ`, grades `{0,2}`), where the **F-channel is
> IDENTICALLY ZERO** (verified `0.000000`).  So Cl(7) adds no working channel for the lepton; on `ℍ`
> the OFE = the `XiVacuum` `{0,2}` dynamics = exactly v58's Cl(3,0) even part → `|ξ|²=½` only as the
> Mexican-hat parameter `μ²/λ` (INPUT), phase a flat Goldstone (FREE).  **⟹ the constrained OFE does
> NOT fix `t²=½` or `φ`.**  Crack closed by color-singlet-ness itself (`ℍ`-confinement ⇒ F=0).
> Confirms the kinematics/dynamics divide (`17_…md`): dynamics gives vacuum STRUCTURE, not MAGNITUDE.
> The constraint-stack reduction (`19_…md`) stands — sector = 1 number `t²`, `φ` follows — but **no
> in-hand dynamics fixes that number**; it remains a potential-parameter input.

The earlier central target (still relevant as the kinematic backbone):

1. **THE PRINCIPLE — derive maximal mixing.**  Show the fermion order parameter is the
   maximally-mixed (G-symmetric, fixed-norm) config on the vacuum manifold — "democratic mass
   matrix."  Forces `t²=(D−dimG₂)/D` (⇒ Koide 2/3,11/15,23/27) AND the democratic v_Higgs
   bilinear, together.
   - **v58-ENERGY ROUTE = CLOSED (2026-05-24, `integration_v58/05_v58_vacuum_alignment.md`).**
     PROVED: `M∈L=skew ⟹ M²=symmetric ⟹ ⟨M²⟩₂≡0`, so the v58 grade-`{0,2}` quadratic fixes the
     vacuum NORM but leaves ALIGNMENT a **flat modulus** — energy-min neither fights nor selects
     maximal symmetry.  ⟹ **maximal symmetry must be a PRIMITIVE selection principle** (the hunch:
     flatness ⟹ energy is silent ⟹ symmetry is more fundamental).  v58 (quadratic-only) is
     *structurally insufficient* to select the vacuum; a quartic alignment term would be needed.
   - **SYMMETRY ROUTE = BUILT (2026-05-24, `integration_v58/06_koide_from_G2_maximal_mixing.md`,
     `g2_koide_derivation.py`, Lean `MaximalMixingKoide.lean` — axiom-clean, in `AxiomCheck`).**
     Constructed `G₂=Aut(𝕆)` from scratch as the derivation algebra (**dim 14 verified** — grounds
     the `14`, was only `dimGL7−dim3Forms7`); decomposed `L=so(8)=G₂(14)⊕nonG₂(14)` (non-G₂ fraction
     **exactly 1/2**); every so(7) bivector uniformly 2/3 in G₂.  RIGOROUS PARTS STAND: `dimG₂=14`
     built, `L=so(8)=14⊕14`, the arithmetic `(D−dimG₂)/D ⇒` Koide ratios (Lean `MaximalMixingKoide`,
     axiom-clean).
   - **GENERATION MAP = BUILT, BRIDGE FAILS (2026-05-24, `08_generation_map_result.md`).**  Built
     `Aut(𝕊)=G₂×S₃`, verified `ψ`=Z₃ generation cycle.  **DECISIVE NEGATIVE:** S₃-covariance forces
     `M` circulant but leaves the amplitude `ξ`/`t²` **FREE** (any `ξ` allowed); max-mixing over the
     generation structure gives `t²∈{0,1}`, **never 1/2**.  So `(D−14)/D=1/2` is a **dimensional
     coincidence in the so(8) Clifford grade, NOT a generation-level mechanism** — P2's bridge has
     **no realization**.  ⟹ **`06`'s "derivation modulo P1+P2" is RETRACTED.**  Koide `2/3` (`t²=1/2`)
     is a **residual coupling magnitude, exactly parallel to the phase `φ=2/9`**: both are couplings
     the `G₂×S₃` symmetry provably leaves free (form & count fixed, amplitude & phase not).  Back to
     "Kepler not Newton," now with a *proven* reason.  (P2 rep-theory scrutiny: see `07_…md`.)
2. **THE LEADING DIRECTION — the missing MULTIVECTOR representation (`09_structure_magnitude_audit.md`).**
   Structure/magnitude audit (all physical core quantities vs all octo-space structures, mapped):
   **every mapped-with-mechanism quantity is DISCRETE** (charges←SU(3), gen-count←S₃, gauge-grp←Spin7,
   mass-grade←J∈Λ², sin²θ_W←Pati-Salam); **every FAILED one is a continuous MAGNITUDE** (`t²,φ,α,v/a²`).
   Reason: we used only octo-space's **symmetry** content (`G₂,S₃` = magnitude-blind automorphisms) +
   grade-counts, **never its multivector representation** (the geometric PRODUCT / associator `[x,y,z]`
   / invariant forms `φ₃,∗φ` / norm — the magnitude-carrying tensors B9–B11).  Symmetry can't fix a
   magnitude; that's why every attempt failed (incl. v58, whose flatness came from projecting to
   symmetry-invariant scalars `⟨M²⟩₀,₂`).  **The two residues type-match: missing constraint is
   magnitude-carrying; octo's only unused magnitude-content is the multivector product.**  NEXT:
   build the lepton order parameter as a full Cl(7) **multivector** (all grades + geometric product),
   impose a geometric-product constraint (idempotent `M²=M` / norm `MM̃` / associator), test whether it
   fixes `t²=½`, `φ=2/9`.  Either finds the constraint, or tests the last octo-structure class (⇒ magnitudes
   genuinely external).  Re-legitimises v58 as the multivector dynamical home.
   - **FIRST RESULT (2026-05-24, `10_multivector_grade_balance.md`): GRADE-BALANCE gives `t²=½`.**
     `M/a=I+ζ`, `ζ=ξS+ξ̄S²`; scalar(grade-0)=`I` norm²=1, non-scalar `ζ` norm²=`2t²`.
     **Grade-equipartition `1=2t² ⟺ t²=½ ⟺ Q=2/3`** (verified; balanced ONLY at ½).  FIRST
     type-correct condition (magnitude from geometric-product grades, not symmetry) to give the
     RIGHT value — vindicates the audit (symmetry gave 0/1).  Lepton-specific (quarks not balanced).
     **CAVEAT: algebraically `≡ t²=½`, so a type-correct REFRAMING, not yet forced.**
   - **FORCING TESTS — NONE force ½ (2026-05-24, `11_grade_balance_forcing_tests.md`).**  Tested
     (1) idempotent, (2) criticality, (3) color singlet, in isolation + all combos.  **All NO:**
     C1 `ζ²=I` unachievable (closest=⅓); C2 no functional critical at ½ (extrema 0.54/0.87/0.95);
     C3 **vacuous** (color⊥generation; reinforces color-3≠generation-3, neg. for Task #1).  Combos:
     C3 vacuous ⇒ 1&3=C1, 2&3=C2; natural joint **1&2 (min idempotent-defect) ⇒ t²=⅓, NOT ½.**
     Recurring value when a condition bites = **⅓** (=`1−Q_lep`, the `−1/3` thread), not ½.
     **NET: nothing in octo-space (symmetry `09` OR multivector grades `11`) forces `t²=½`; weight of
     evidence ⇒ `t²=½` is likely a genuinely EXTERNAL magnitude.**  THE ONE SURVIVING THREAD = C1
     conditional on the **generation↔internal map** carrying the internal idempotent ½ to the
     generation `t²` (same map unbuilt since `08`).  Make-or-break: build it (⇒ derive Koide) or show
     it can't carry the ½ (⇒ `t²` external, question closed).  Phase `φ=2/9` still separate/untouched.
   - **OCTOMATH — WHAT the ½ IS (2026-05-24, `13_octomath_half.md`, Lean `OctoHalf.lean` axiom-clean).**
     Half-element law (any ring): `(1+u)²=2(1+u)+(u²−1)`, so `P=(1+u)/2` is controlled by `u²` (fixed by
     GRADE via `BladeSquareSign`): **`u²=−1` (L, COMPLEX) ⇒ `P²=P−½`, root-product ½; `u²=+1` (F, REAL)
     ⇒ idempotent, root-product 0.**  Since `mass∈L` (`u²=−1`) is PROVEN, **the ½ is the L-grade
     complex-structure root-product — NOT arbitrary; an F-grade (real) mass would give 0.**  Explains the
     recurring `⅓`/`0` (the `u²=+1` real-idempotent branch).  The internal ½ is now *grounded & Lean-locked*;
     the open part is unchanged (transport to generation `t²` = the map), but its job is now sharp: carry
     the L-grade complex-structure root-product ½ to `t²`.
   - **THE MAP — PARTIAL BUILD (2026-05-24, `14_generation_internal_map.md`).**  Identify the coupling
     `ξ` with the L-grade complex half-element `(1+u)/2`: reversion-norm `ξξ̃=(1−u²)/4=½` = a SCALAR =
     `|ξ|²=t²`.  **So `t²=½` IS carried to the amplitude — the Koide ½ = the L-grade complex
     reversion-norm** (L: scalar ½; F: `(1+w)/2`, NOT a scalar — no norm; `mass∈L` proven ⇒ ½).  Strongest
     amplitude anchor yet.  **BUT: (i) the PHASE is FREE** — the half-element's `45°` is the AMPLITUDE
     balance (`|Re ξ|=|Im ξ| ⇒ |ξ|²=½`), NOT a phase; `|ξ|²=½` holds for any `arg(ξ)`, so `φ` is a free
     rotation on the `|ξ|²=½` circle (residual, `φ=2/9` near the chiral limit, not π-rational, far from
     `45°`); **(ii) not a clean forcing** — `|ξ|²=½` holds for any
     phase (= the grade-balance, now *carried by* the half-element, value ½ structural, but "why
     half-norm" unchanged); **(iii) HARD obstruction:** generation `ζ=ξS+ξ̄S²` is NEVER a complex
     structure (`⟨ζ²⟩₀=+2t²≥0`), so the carrier is `ξ`, not `ζ`.  **VERDICT: amplitude ½ grounded as the
     L-grade complex reversion-norm (proven-anchored); phase `2/9` separate open residual; clean forcing
     not achieved.**  `t²=½` is now "the L-grade complex half-norm," not a bare number.
   - **PHASE `φ=2/9` — M4 chiral window, Lean-locked (2026-05-24, `15_chiral_phase_window.md`,
     `ChiralPhaseWindow.lean` axiom-clean).**  In the Koide amplitude, the electron-massless point is
     `φ=π/12` (`chiral_massless_point`); physical `φ=2/9 < π/12` (`physical_phase_below_chiral`), electron
     light-but-massive (`electron_light_but_massive`).  So `φ=2/9` sits just inside the light-electron
     window (`√m_e/a≈π/12−φ≈0.04`).  **CONSTRAINT, not derivation:** chiral protection ⇒ `φ` near `π/12`
     (window), not the point; `2/9=Q/3` is the separate phase law; offset `π/12−2/9` set by `m_e` (input);
     `π/12` π-rational, `2/9` not.  **Phase `φ=2/9` = free residual, near the chiral edge; NOT derived.**
3. Task #1: color-3 (SU(3) charge ⅓) vs generation-3 (S₃) — same origin?  Plus: the "why half-norm"
   amplitude forcing and the phase `φ=2/9` value (beyond the chiral window) remain open.
4. (Lower) Higgs residuals R1/R2; quark Koide selection rule (undriven `D_N∈{28,35,63}`).

## 7. Long-term goals

1. Make the "≈1 input" claim **unconditional**: every lepton/EW/Higgs relation *derived* from the
   v58 equation, not fitted — leaving only the irreducible scale `a_lepton`.
2. **Derive the residual `Q=2/3`** (the magnitude of `φ=Q/3` and the Koide ratio) — currently
   mechanism-less; provably not a holonomy/symmetry phase; only structural target is `cos(2/3)`.
3. **Engage the untouched SM sectors:** quark scales/phases (the selection rule), CKM, neutrinos,
   strong sector.
4. **The deep open question:** is the family-`3` (sedenion `S₃`) the *same* `3` as the spatial-`3`
   (`Im ℍ`)?  (Ties to v58's emergent-3D program.)
5. Full v58↔v59 unification: the field equation as the dynamical law whose vacuum/constraint
   structure *is* the v59 algebra (even/ℍ sector + grade↔complex law, both established).

## 8. Standing caveats (do not lose)

- Only the **algebraic skeleton is proven**; physics predictions are empirical/conjectural matches.
- `a_lepton` is **dimensionful** → never "derivable" as a number; derive dimensionless ratios + `α`.
- The **octonion `Aut=G₂`** does NOT contain the generation `Z₃`; the "/3" lives in the **sedenion
  `S₃`** (`Aut(𝕊)=G₂×S₃`) — confirmed, and it's why the intra-octonion holonomy test was trivial.
- v58 (`Cl(3,0)`) ≠ v59 (`Cl(7)`) as full algebras; the bridge is the **even/ℍ sector + the
  signature-independent grade↔complex law**, NOT a full identification.  (Cf. the retracted
  octonionic "living-candidate crossover.")

## 9. Key files
- Budget/audit/plan: `PARAMETER_BUDGET.md`, `RIGOR_AUDIT.md`, `INTEGRATION_PLAN_v58.md`.
- φ=Q/3 study: `koide_phase_law/` (`CONCEPT.md`, `holonomy_result.md`, `sedenion_s3.py`,
  `LITERATURE_three_from_octonions.md`).
- v58 integration: `integration_v58/01_higgs_vev.md`, `02_higgs_iii_yukawa.md`.
- Lean: `lean/` (build all + axiom audit via `lean/AxiomCheck.lean`); new this phase —
  `BrannenPhase`, `WeinbergPatiSalam`, `TwoNinthsUnification`, `PhaseExclusions`, `PhaseAmbiguity`,
  `LeptonPhaseEmpirical`, `ColorSU3`, `HiggsVevReframe`.
- Memory: `~/.claude/projects/-home-d-code-scp/memory/project_v59_octonionic_kernel.md`.
