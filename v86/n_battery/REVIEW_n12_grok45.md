I'll read the v86 program docs, the N-battery review protocols, the gauged shooter conventions, and the N1/N2 instrument in order, then stress-test the derivation and code.Next I'll check the shooter's energy conventions and the theory docs that fix the signs and scaling.# v86 N1/N2 implementation review (adversarial)

**Seat:** v86 implementation-review · **Scope:** derivation I1–I4, `v86/n_battery/n12_decomp.py`, protocols in `N_BATTERY_REVIEW.md` §1 · **Files not modified.**

Smoke-check (production `gauged_shooter`, H=0.004): at ω=1.45, g=0 and g=0.05, I1 residual ≲ 10⁻⁹, |R_vir|/E ~ 1.4×10⁻⁷, |I4|/E ~ 5×10⁻⁸, ε ~ 4.0–4.3%. FAST grid (default) is similar at ~10⁻⁶–10⁻⁷. So the algebra is not vapor — but protocol and foundation-prose issues remain.

---

## 1. Derivation I1–I4 vs shooter / theory

### Finding 1 — I1 (Gauss) is correct; tail treatment is consistent

**Attacked text (docstring):**
> `rho_Q = 3 wt f^2 = -lap(chi)/g^2`, so  
> `int chi rho_Q dV = 2 E_g`, and `w Q = 2(E_kin + E_g)`.

Matches `gauged_shooter.py` conventions: `wt = w - chi`, `Q = ∫ 3 wt f² dV`, `E_kin = ∫ (3/2) wt² f² dV`, `E_g = ∫ χ'²/(2g²) dV`, Gauss `∇²χ = −g² ρ_Q`.

With \(w = wt + χ\):
\[
wQ = \int 3\,wt^2 f^2\,\mathrm{d}V + \int χ\,\rho_Q\,\mathrm{d}V = 2E_{\mathrm{kin}} + 2E_g
\]
on solutions (and under integration by parts with correct surface terms).

**Tail:** Docstring says exterior ρ contributes nothing, so use tail-corrected \(E_g\). Slightly loose wording, but the **implementation is right**: truncated-domain IBP leaves a surface term at \(R_{\max}\) equal to \(2E_g^{\mathrm{tail}}\) for \(χ\sim C/r\). Adding
`4π C²/(2 g² rN)` matches both shooter’s `observables()` and I1. Smoke: I1/E ~ 10⁻¹⁶ (g=0) and ~10⁻⁹ (g=0.05).

### Finding 2 — Minus sign on \(E_g\) in \(F\) is correct (saddle)

**Attacked text:**
> `F[f,chi] = E_grad + E_m + E_V - E_kin - E_g`  
> (gauge energy with a MINUS: electrostatics is a saddle in \(a_0\)).

Free variations of this \(F\) recover **both** shooter ODEs:
- \(\delta F/\delta χ = 0\) → Gauss  
- \(\delta F/\delta f = 0\) → matter radial ODE with \(W = 2 V_t' f^5\)

Hamiltonian energy keeps **\(+E_g\)**; the **reduced free functional** that is stationary in \((f,χ)\) has **\(−E_g\)**. That is standard for electrostatics (constraint / saddle in \(A_0\)). **The minus sign is right.**

### Finding 3 — I2 is exact but not independent of I1

On the Gauss constraint, \(Σ = E − wQ = F\) is algebraic:
\[
Σ - F = 2(E_{\mathrm{kin}}+E_g) - wQ = -(\text{I1 residual}).
\]
So **I2_rel ≡ −I1_rel**. Reporting both as two independent identities overstates the test suite. I2 is a rewrite check, not a second physical constraint.

### Finding 4 — Derrick scalings and I3 are correct for free \(F\)

**Attacked text:**
> \(E_{\mathrm{grad}}, E_g \sim L\); \(E_m, E_V, E_{\mathrm{kin}} \sim L^3\);  
> \(R_{\mathrm{vir}} = (E_{\mathrm{grad}}-E_g) + 3(E_m+E_V-E_{\mathrm{kin}}) = 0\).

Under the free-functional family \(f(r)\mapsto f(r/L)\), \(χ(r)\mapsto χ(r/L)\) at fixed \(ω\), those scalings hold in 3D, and \(dF/dL|_{L=1}=0\) yields exactly I3. Scaled \((f,χ)\) need **not** obey Gauss; Derrick only needs a variation of \(F\).

**Caveat (not a kill):** constrained scaling (eliminate \(χ\) via Gauss, then dilate \(f\)) has different naive powers for \(E_g\). The free-functional route is the right one for critical points of \(F\). Numerics already show \(|R_{\mathrm{vir}}|/E \sim 10^{-7}\).

### Finding 5 — I4 follows; ungauged positivity follows

From \(Σ = (E_{\mathrm{grad}}-E_g) + (E_m+E_V-E_{\mathrm{kin}})\) and I3:
\[
Σ = \tfrac23(E_{\mathrm{grad}}-E_g),\qquad
\varepsilon = \frac{E_{\mathrm{grad}}-E_g}{3(E_{\mathrm{kin}}+E_g)}.
\]
Ungauged: \(Σ = \tfrac23 E_{\mathrm{grad}} > 0\). Correct.

### Finding 6 — FOUNDATIONS / GROUNDING prose has the wrong gauge sign for \(Σ\)

**Attacked (FOUNDATIONS R2):**
> \(Σ = E_∇ + E_m + E_V + E_g - E_{\mathrm{kin}}\)  
> with \(E_{\mathrm{kin}}=ωQ/2\).

**Attacked (GROUNDING / THEORY_v86):**
> “gradient + potential + **gauge** − kinetic imbalance”

With the shooter’s **true** kinetic \(E_{\mathrm{kin}}=\int\tfrac32 wt^2 f^2\):
\[
Σ = E_∇ + E_m + E_V - E_{\mathrm{kin}} - E_g.
\]
With throughput \(E_{\mathrm{kin}}^{\mathrm{tp}} ≡ ωQ/2 = E_{\mathrm{kin}}+E_g\):
\[
Σ = E_∇ + E_m + E_V - E_{\mathrm{kin}}^{\mathrm{tp}}
\]
(**no** \(+E_g\)). FOUNDATIONS’ \(+E_g\) is wrong once gauged (error \(2E_g\) if mixed with true kinetic).  

**Instrument derivation supersedes FOUNDATIONS R2 on the gauge sign.** Any N1 writeup must not quote FOUNDATIONS’ formula uncorrected.

### Finding 7 — N1 protocol’s \(E_{\mathrm{kin}}=\tfrac12ωQ\) vs instrument

**Attacked (N_BATTERY_REVIEW N1):**
> \(E_{\mathrm{kin}} = \tfrac12 ωQ\) (throughput; check vs integrated \(\tfrac32 ω^2 f^2\))

Instrument correctly uses integrated \(\tfrac32 wt^2 f^2\). At \(g>0\), \(\tfrac12 ωQ = E_{\mathrm{kin}}+E_g\), not \(E_{\mathrm{kin}}\). Using protocol \(E_{\mathrm{kin}}\) inside the F-decomposition would **corrupt** Σ fractions and I4. Instrument choice is right; **label columns** so they are not compared naively to FOUNDATIONS/N1 “throughput” \(E_{\mathrm{kin}}\).

---

## 2. Code review (`n12_decomp.py`)

### Finding 8 — Energy pieces mirror `G.observables()` (good)

Same extend-to-boundary, `np.gradient(..., edge_order=2)`, `trapz(..., dx=H)`, Coulomb tail, \(Q\)/matter factors. Smoke: piece sum equals `Em+Ef` to machine epsilon. Docstring claim is empirically true.

### Finding 9 — No hard assert against `G.observables()` (gap, not math bug)

Docstring says pieces are “guaranteed” to sum to `observables()`’s \(E\), but `decompose()` never calls it. One line
`assert abs(E - G.observables(...)['Et']) < tol` would catch silent drift if either path changes.

### Finding 10 — Default `N12_FAST=1` is coarser than production / gscan

| | Production shooter | FAST (default) |
|--|--|--|
| H | 0.004 | 0.01 |
| RMAX | 150 | 100 |

Identities still look good on FAST (~10⁻⁶–10⁻⁷), and ω=1.45 matches gscan at displayed precision — fine for a scout. **Archival / “TOE-grade” N2 should set `N12_FAST=0`** (or pin that residuals are grid-limited).

### Finding 11 — `np.gradient` / r=0 (acceptable, same as shooter)

Neither path enforces \(f'(0)=χ'(0)=0\) in the energy gradient. \(r^2\) weight kills the origin; residual level is consistent with discretization, not a showstopper.

### Finding 12 — g-continuation is sound

Same \(χ \propto g^2\) guess and dg-halving as the shooter. Rebuilding each \(g\) from the g=0 base is slower but correct. No bug found.

### Finding 13 — I1/I2 double-counting in the summary table

Max |I1_rel| and |I2_rel| will always match. Prefer one Gauss residual + R_vir + I4.

### Finding 14 — Composition table will look “wild” (expected, not a bug)

\(E_m\) and \(E_{\mathrm{kin}}\) are large and cancel inside \(Σ\). Fractions of \(Σ\) of order tens are normal. Primary N1/N2 diagnostics should be **I4** and \(\tfrac23(E_{\mathrm{grad}}-E_g)\), not raw \(E_m/Σ\).

### Finding 15 — No ODE residual emitted

N2 pass is “at or below shooter ODE residual.” Newton residual exists inside `solve` but is discarded. Without it, “≪ 10⁻³” is a soft slogan, not a per-row gate.

---

## 3. Protocol compliance (N1 restored / N2)

### Finding 16 — N1 Phase B: substantially met; Phase A: missing

| Protocol element | Status |
|--|--|
| Full split \(E_∇, E_m, E_V, E_g, E_{\mathrm{kin}}\) | **Yes** |
| \(Σ\), ε vs Q/g | **Yes** |
| Consistency \(E=\sum\) pieces | **Implicit**, not asserted |
| Phase A on existing `gscan.tsv` | **No** |
| Outcome table (surface / bulk \(E_V\) / Coulomb gap / thin-wall ε→0) | **Not automated** — only raw fractions |
| Kill: parts fail to rebuild \(E\) ≫ shooter tol | **Not coded** |

Phase B is what actually answers “is \(Σ\) the virial imbalance?” Phase A is still a free cross-check against the frozen gscan artifact and should not be skipped if the goal is comparability to the FOUNDATIONS design.

### Finding 17 — N2: identity form good; pass/fail packaging incomplete

| Protocol element | Status |
|--|--|
| Parameter-free relation among energy pieces | **Yes** (Derrick of \(F\); equiv. class to \(f/rf'\) multiply) |
| Evaluate \(R_{\mathrm{vir}}\) on branch | **Yes** |
| Pass: \(\|R_{\mathrm{vir}}\|/E \ll ε\), ≲ ODE residual | **Partially** — residual printed; no ODE floor; no pass/fail bit |
| Fail modes (systematic ~ε → wrong identity; ≫ODE → bookkeeping alarm) | **Not pre-registered in script output** |
| TOE-grade: ε(Q) as derived function of moments | **I4 is exactly that** if residuals land |

### Finding 18 — What is MISSING for non-ambiguous FOUNDATIONS comparison

1. Explicit convention block: true \(E_{\mathrm{kin}}\) vs \(\tfrac12ωQ\); \(F\) has \(−E_g\); correction of FOUNDATIONS R2.  
2. Phase A: reload `v69/theory/gscan.tsv`, compare \(E,Q,E_{\mathrm{field}}\) at shared \((g,ω)\).  
3. Per-row Newton residual + automated PASS/FAIL vs \(\max(|R_{\mathrm{vir}}|/E, |I4|/E) \ll \min(ε)\).  
4. N1 outcome columns: e.g. \(Σ\) vs \(4π R_{\mathrm{half}}^2\), \(E_g\) vs \(g^2 Q^2/(8π R)\), \(E_V\) fraction of \(Σ\) after I4 reduction, ε(g=0) vs ε(g>0) gap vs \(E_g\).  
5. Production grid for the published table (`N12_FAST=0`).  
6. Do not treat I1 and I2 as two independent pillars.

Without (1)–(4), two readers can both “pass N1/N2” and still disagree with FOUNDATIONS’ written \(Σ\) formula and with each other’s surface/bulk verdict.

---

## 4. Pre-registration

### Finding 19 — What falsifies “\(Σ\) is the virial imbalance”

Declare **falsified** if any of the following hold on the VK-stable branch (away from edges where the solver is soft):

1. \(\max |R_{\mathrm{vir}}|/E \gtrsim \varepsilon\) **and** systematic (same sign, smooth in \(ω\)) — I3 wrong/incomplete.  
2. \(\max |Σ - \tfrac23(E_{\mathrm{grad}}-E_g)|/E \gtrsim \varepsilon\) while pieces still sum to \(E\) — I4/algebra broken.  
3. \(\max |wQ - 2(E_{\mathrm{kin}}+E_g)|/E \gg\) ODE/integration floor — Gauss or energy convention broken (**N1 kill / D5′ alarm**).  
4. \(Σ\) composition after I4 is **not** carried by \((E_{\mathrm{grad}}-E_g)\), but by an independent bulk that refuses the virial reduction — “imbalance” narrative fails even if a residual is small by accident.

**Does not falsify:** large cancelling \(E_m\) and \(E_{\mathrm{kin}}\) in the raw partition (expected).

### Finding 20 — If I4 holds to ~10⁻⁶ but ε is 1–4%

That is the **success case**, not a tension.

| Quantity | Meaning |
|--|--|
| \(\varepsilon \sim 1\)–\(4\%\) | Physical Legendre excess (signal) |
| \(\|I4\|/E \sim 10^{-6}\)–\(10^{-7}\) | Discretization / integration floor (noise on the identity) |

**Report:**
- \(ε(Q,g)\) and \(Σ=\tfrac23(E_{\mathrm{grad}}-E_g)\) as the **derived** closed form.  
- Identity residuals **separately**, with pass iff residual \(\ll \min ε\) on the sample (target: ≲ ODE residual, typically ≪ 10⁻³).  
- **Do not** require residual ~ ε, and **do not** claim “ε is numerical error.”  
- Narrative fix for the program: \(Σ = E_∇+E_m+E_V-E_{\mathrm{kin}}-E_g = \tfrac23(E_∇-E_g)\); FOUNDATIONS’ \(+E_g\) wording is superseded.

Near the capacity fold, pre-register watching \(E_g\) vs \(E_{\mathrm{grad}}\): I4 allows \(Σ<0\) if \(E_g>E_{\mathrm{grad}}\) — interesting if seen; not automatically a kill until VK/existence edges are checked.

---

## VERDICT

### **RUN WITH FIXES**

Derivation I1–I4 matches the shooter’s conventions and the variational structure of the gauged radial system. The **minus on \(E_g\) in \(F\) is correct**; Derrick powers are correct; I4 follows; tail-corrected \(E_g\) is consistent with I1. Code is largely a faithful split of `observables()` and is already numerically healthy.

**Required before treating the run as battery-grade (not before a scout run):**

1. **Convention banner in output/results:** true \(E_{\mathrm{kin}}=\int\tfrac32 wt^2f^2\); \(Σ=F=E_∇+E_m+E_V-E_{\mathrm{kin}}-E_g\); note FOUNDATIONS R2 \(+E_g\) is wrong when gauged.  
2. **Assert** `E` vs `G.observables(...]['Et']` (and optionally \(Q\)) per row.  
3. **Phase A:** diff key rows against `v69/theory/gscan.tsv` (\(E_{\mathrm{total}}\), \(E_{\mathrm{field}}\), \(Q\)).  
4. **Record Newton residual**; gate PASS on \(\|R_{\mathrm{vir}}\|/E\) and \(\|I4\|/E \ll ε\) and ≲ that floor.  
5. **Drop or mark I2** as duplicate of I1.  
6. For the frozen Part-0 table: **`N12_FAST=0`** (or document FAST as scout-only).  
7. **Outcome-table hooks** (even crude): \(Σ\) vs \(R_{\mathrm{half}}^2\), \(E_g\) vs point-Coulomb estimate, g=0 vs g>0 ε gap vs \(E_g\).

**Optional scout:** run as-is with `N12_FAST=1` is acceptable for a first look; do not promote those residuals as the TOE-grade N2 landing without (2)–(6).

**Not DO NOT RUN:** no derivation-breaking bug found; no integration convention that would invent a false I4 at the ε scale.