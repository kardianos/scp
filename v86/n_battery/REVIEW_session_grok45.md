I'll adversarially review all four results: read sources, re-derive the math by hand, re-run cheap checks, and hunt for mistakes rather than agreement.Next I'll pull the logs/TSVs, the GSS statement, and the kernel/stencils that these claims rest on.The log may predate the window filter; I'll analyze the TSV, re-derive the key expansions, and re-run the cheap checks.Checking a few more details while the independent Q(ω) scans run.# Adversarial review — four session results

Re-derived expansions by hand, re-ran cheap numerics, compared code to claims. Effort is on what is wrong or overstated.

---

## 1. HC-1-gauged

### Finding 1.1 — Linearised c-coupling  
**Severity: OK**  
**Claim:** “expanding \((D-ic)^2(f+\eta)\) … gives the c-coupling \(-i\dot c\,f + 2\,\mathrm{wt}\,c\,f\)”  
**Check:** With \(D=\partial_t+i\,\mathrm{wt}\) and \(\mathrm{wt}\to\mathrm{wt}-c\),
\[
(D-ic)^2(f+\eta)=D^2(f+\eta)-i\dot c\,f+2\,\mathrm{wt}\,c\,f+O(2).
\]
EOMs then give \(-2\,\mathrm{wt}\,f\,c\) in \(x\) and \(+f\dot c\) in \(y\). Sign and structure match the file and the implemented operators.

### Finding 1.2 — Linearised Gauss  
**Severity: OK**  
**Claim:** \((-\Delta+3g^2 f^2)c=2g^2\,\mathrm{wt}\,f\sum_a x_a\)  
**Check:** Background \(\Delta\chi=-3g^2\,\mathrm{wt}\,f^2\). Static linearisation of \(\rho=\sum_a(\mathrm{wt}-c)(f+x_a)^2\) yields screening \(3g^2 f^2\) (three components) and source \(2g^2\,\mathrm{wt}\,f\sum x_a\). Factors 3, 2, and sign all check out.

### Finding 1.3 — Coefficient 12 and PSD  
**Severity: OK**  
**Claim:** symmetric elimination gives \(L_x^\mathrm{sym}+12g^2(\mathrm{wt}f)K^{-1}(\mathrm{wt}f)\)  
**Check:** \(\sum x_a=3\xi\), \(Kc=6g^2(\mathrm{wt}f)\xi\), static EOM \(L_x\xi+2\,\mathrm{wt}f\,c=0\) ⇒ coefficient **12**. \(K>0\) ⇒ the Coulomb piece is PSD and **cannot add** a negative direction (it can only raise eigenvalues / possibly remove negatives). Numerics: all Coulomb shifts \(+0.029\)–\(+0.068>0\).

### Finding 1.4 — Headline: window / \(Q_\max\) “not GSS”  
**Severity: MAJOR**  
**Claim:** “\(Q_\max=921\) and the narrowed window \((1.406,1.5)\) are the BACKGROUND ceasing to exist, NOT a GSS index change.”  
**What is wrong:** This conflates **two different edges**.

| Edge | What happens | GSS content |
|------|----------------|-------------|
| Low-\(\omega\) / \(Q_\max\) | No solution for \(\omega\le 1.405\); \(Q\to 921\) at \(\omega\approx 1.4062\) (gscan) | Existence only; \(n(H_\omega)=1\) on surviving points |
| High-\(\omega\) | Branch **exists** past a **VK turning point** | GSS via \(n(D)\), not \(n(H)\) |

Independent gauged continuation (same shooter):
- \(Q_\min\approx 89.70\) at \(\omega\approx 1.485\)
- \(\mathrm{d}Q/\mathrm{d}\omega\): negative through 1.48 (\(Q=90.95\)), positive at 1.49 (\(Q=91.86\), \(\mathrm{d}Q/\mathrm{d}\omega\sim +91\))

HC-1’s own table has \(n(H_\omega)=1\) at **both** 1.48 and 1.49. Classic VK: \(n(H)\) stays 1 while \(n(D)\) flips when \(\mathrm{d}Q/\mathrm{d}\omega\) crosses zero. So half the window story **is** a GSS index mismatch on the \(D\) side; the file only measured \(n(H)\) and then over-claimed.

**Corrected statement:** Low-\(\omega\) termination and \(Q_\max=921\) are existence limits with \(n(H_\omega)=1\). High-\(\omega\) has a real VK turning (\(Q_\min\sim 89.7\) near \(\omega\sim 1.485\)); \(n(H_\omega)=1\) on both sides means the GSS change is carried by \(n(\partial Q/\partial\omega)\), which this rung did not compute.

### Finding 1.5 — Goldstone floor  
**Severity: MINOR**  
**Claim:** count \(L_0\) negatives against \(10\times|\lambda_\mathrm{Goldstone}|\)  
**What is wrong:** Not circular for the main claim (\(n(H)\) is carried by \(L_x^\mathrm{sym}\) at \(O(10^{-1})\), far above the \(O(10^{-6})\) floor). It **can** hide a genuine tiny \(L_0\) negative comparable to discretisation. At \(\omega=1.49\), gold residual rises to \(3.2\times 10^{-5}\) (thin-wall / box stress)—floor is looser there.

**Corrected:** Acceptable for “no second \(L_0\) direction at resolved scale”; not a proof that \(L_0\) is free of all negatives below discretisation.

### Finding 1.6 — \(l=0\) only  
**Severity: MAJOR**  
**Claim:** “\(n(H_\omega)=1\) everywhere the gauged branch exists”  
**What is wrong:** Only the spherical sector is built (`base = M2 - wt^2 + P0  # l=0`). Index is **not** complete. \(l=1\) should be translational (near-zero); \(l\ge 2\) can host fission/multipole negatives, especially near large \(Q\) (\(\omega=1.41\), \(Q\sim 529\)–\(921\)).  

**Corrected:** \(n(H_\omega)^{(\ell=0)}=1\) on the scanned gauged branch; full \(n(H_\omega)\) needs \(\ell\ge 1\).

### Finding 1.7 — “BdG spectrum” vs Hessian  
**Severity: MINOR**  
Docstring claims a gauged BdG spectrum; code only diagonalises static \(L_0,L_x^\mathrm{flav},L_x^\mathrm{sym}+\mathrm{Coulomb}\). That is the right object for the GSS Morse index, but it is **not** the dynamical BdG frequency spectrum.

---

## 2. HC-3-volume

### Finding 2.1 — Parametrisation / fundamental domain  
**Severity: OK** (one pedantic note)  
\(\delta_a=\rho\cos(\theta-120^\circ a)\): \(\sum\delta_a=0\) to machine precision.  
\(S_3\cong D_3\) on the traceless plane ⇒ fundamental domain length \(60^\circ\): **[0,60] OK**.  
Old ray \(\propto(-1,+1/2,+1/2)=\theta=180^\circ\); cycle \((d_1,d_2,d_0)\) maps it to \(\theta=60^\circ\): **correct**.  
Edges \(\theta=0\) (one up, two down) vs \(\theta=60\) (two up, one down) are physically distinct: **correct**.

### Finding 2.2 — VK turning / \(n(D)=0\) headline  
**Severity: MAJOR** (turning real; interpretation overstated)  
**Claim:** 20 window-valid points with \(n(D)=0\) near \(\bar\omega=1.48\), due to VK turning at \(\omega\sim 1.481\) with \(Q_\mathrm{min}\sim 86.7\).

**Independent re-solve** (`n4_hc3_flavored` symmetric seed, \(H=0.01\), \(R_\max=60\)):

| \(\omega\) | \(Q_\mathrm{tot}\) | \(\Delta Q/\Delta\omega\) | \(n(D)\) (sym) |
|----------|-------------------|--------------------------|----------------|
| 1.480 | 86.727 | −250 | 1 (ev₀=−54) |
| 1.482 | **86.590** | −68 | **0** (ev₀=+13) |
| 1.484 | 86.902 | +156 | 0 |
| 1.490 | 92.336 | +1425 | 0 |

Turning is **real**, not a continuation artifact. Minimum is \(\omega\approx 1.482\), \(Q\approx 86.59\) (not exactly 1.481 / 86.7, but the same feature).

**What is overstated:**
1. At \(\bar\omega=1.48\), **symmetric** \(n(D)\) is still 1; the 19 points at 1.48 are **detuned** partitions that cross zero because the symmetric eigenvalue is already \(O(10)\)–\(O(50)\) from the turn.
2. Past the turn, \(n(D)=0\) with provisional \(n(H)=1\) is the **ordinary monochromatic VK-unstable branch**, not a novel multi-flavor HC-6 discovery. Framing all 20 as “partition-volume index boundary / HC-6 targets” oversells flavoured novelty.
3. **9/20** of the \(n(D)=0\) points have some \(w_a>1.495\); several show soft pathology (e.g. \(\mathrm{ev}\sim 10^5\)–\(10^6\), large \(Q\) jumps when \(\rho\) increases). Window cut at 1.5 is necessary but not sufficient.

**Corrected:** Symmetric branch has a real VK turn at \(\omega\approx 1.482\), \(Q_\min\approx 86.59\). Detuned points near \(\bar\omega=1.48\) inherit that crossing. Clean interior \(n(D)=0\) points are VK-side unstable states; near-edge points need a tighter quality cut before HC-6 seeding.

### Finding 2.3 — Window filter  
**Severity: MINOR**  
Hard cut \((1.3087,1.5)\) correctly kills the pre-filter \(n(D)=2\) / \(|\lambda|\sim 10^6\) junk (log had those; filtered TSV has **zero** \(n(D)=2\)). Survivors at \(w\sim 1.498\)–\(1.4997\) are still soft (large FD eigenvalues, \(Q\) sensitivity).  

**Corrected:** Filter is right as a first cut; claim “20 genuine branch points” needs residual/edge diagnostics, not only \(w_a\in(w_\mathrm{lo},w_\mathrm{hi})\).

### Finding 2.4 — GSS inference  
**Severity: OK** with scope flag  
**Claim:** \(n(D)=0\) while \(n(H_\omega)=1\) ⇒ GSS mismatch ⇒ HC-6 target.  

`GROUNDING.md` §1: stability needs \(n(H_\omega)=n(\partial Q_a/\partial\omega_b)\). Single-charge check: \(\mathrm{d}Q/\mathrm{d}\omega<0\Rightarrow n(D)=1=n(H)\).  
So \(n(D)=0\neq 1=n(H)\) ⇒ **unstable** — not inverted.  

Caveat the file understates: \(n(H)=1\) is still the **ungauged** one-hump assumption (HC-1-gauged only did \(\ell=0\), \(g=0.05\), monochromatic). HC-3-volume is ungauged flavoured; using gauged HC-1 as if it closed \(n(H)\) for this scan is a stretch.

### Finding 2.5 — Asymmetry \(2.03\times 10^{-3}\)  
**Severity: OK** for signature, **MINOR** near zero crossing  
Filtered max \(|D-D^T|/\max|D|=2.03\times 10^{-3}\) (matches TSV; pre-filter log’s \(1.28\times 10^{-1}\) is the unfiltered run).  
For clear eigenvalues (\(|\lambda|\gtrsim 10\)) this will not flip signs. For the softest positive \(\mathrm{ev}_0\sim 4.3\) (θ=60, ρ=0.015), absolute error \(\sim 10^{-4}\times\max|D|\sim 0.2\) still leaves the sign stable. Trustworthy enough for \(n(D)\) counts on the reported targets; not for precise eigenvalue magnitudes near the turn.

---

## 3. D8b

### Finding 3.1 — Telescoping identities  
**Severity: OK**  
Independent random-field check of the \(d\neq i\) identity: max relative residual \(2.6\times 10^{-16}\). \(d=i\) similarly. Symbolic re-run: 7/7 PASS.

### Finding 3.2 — Mass term  
**Severity: OK**  
\(\sum_n m^2 f_n(f_{n+i}-f_{n-i})=0\) by periodic relabelling (any dx). Confirmed.

### Finding 3.3 — “No exact lattice momentum”  
**Severity: OK** (standard subtleties acknowledged)  
**Claim:** no exactly conserved lattice momentum in general; translation broken to a discrete subgroup ⇒ no Noether current (unlike continuous U(1) for charge).  

**Check:** Continuous Noether requires a continuous symmetry. Lattice shift \(\mathbb{Z}^d\) is discrete ⇒ no continuous Noether current. The continuum-like \(\sum \pi\,D_c\phi\) is **not** conserved once \(U'\) is nonlinear; the residual is the PN force. Special integrable discretisations (Ablowitz–Ladik, etc.) can conserve extra quantities—this stencil is not one of them.  

I did **not** find a missed exact continuum-style \(P\) for this kernel stencil. Claim stands.

### Finding 3.4 — Peierls–Nabarro factor 6 and “not pinning”  
**Severity: OK** for factor; **MINOR** for production-dx wording  
\(V=V_t(s)\), \(s=\prod|\Phi_a|^2=f^6\), equal components:  
\(\partial V/\partial u_a=2V'(s)f^5\) per real component × 3 ⇒ **\(6\,V_p f^5\)**. Correct.  

Measured R/E: \(4.27\times 10^{-5}\) (dx=0.60) → \(2.03\times 10^{-9}\) (dx=0.20). Far below 4.4%.  

**Minor:** “at production dx” is vague; N7 defect was on a specific mesh. Ordering (PN ≪ 4.4%) is robust across the scanned dx.

### Finding 3.5 — Match to `sfa_momentum.c`  
**Severity: MAJOR**  
**Claim (implicit):** derived flux is the correct object; 4.4% defect is the wrong surface object.  

`sfa_momentum.c` still samples continuum  
\(T^{ji}=\sum \partial_j\phi\,\partial_i\phi+\delta_{ji}L_\mathrm{matter}\)  
at cell centres with central differences—not the link flux of D8b.  

So:
- PN ≪ 4.4% **does** rule out pinning as the defect source.
- D8b did **not** re-integrate the link flux on the N7 run and show the residual collapse.
- “Defect is the surface object and Part 2a gives the corrected one” is a derivation claim, **not** a closed re-measurement.

**Corrected:** Pinning is not the 4.4% defect. Whether the link flux restores balance on N7 data is untested.

---

## 4. v87 B1 CHSH

### Finding 4.1 — Algebraic \(|S|\le 2\) identity  
**Severity: OK** (wording nearly airtight)  
**Claim:** offline per-object dichotomic readout ⇒ per-frame CHSH bracket \(\in\{-2,+2\}\) ⇒ sample mean \(|S|\le 2\) algebraically.  

Enumeration over all \(A,A',B,B'\in\{\pm 1\}\): bracket only \(\pm 2\) (8+8). If \(A=A(a,\lambda)\), \(B=B(b,\lambda)\) only, this is identity, independent of the PDE. Gauge-invariant Wilson line is setting-independent ⇒ bound intact.  

**Corrected nuance:** \(S>2\) would still be evidence of an **analysis bug**, not of fabric nonlocality. The writeup already scopes the run as a tripwire; “no evidence about locality in either direction” is right for **locality**.  

CRANK2’s recommended test is not falsifiable in the \(S>2\) direction under this protocol: **correct**.  
“Settings as fabric DOF is a prerequisite, not a parallel gap”: **correct** for a genuine locality test.

### Finding 4.2 — `bell_grid` arm 1 can exceed 2  
**Severity: OK**  
Code builds one histogram \(h[j]\), then  
`E[k]` by cyclic shift of the response against that single histogram (`bell_grid.c` ~103–110). That assumes translation-invariant \(\rho(\lambda)\). Finite-\(N\) histogram noise ⇒ four CHSH terms are not from consistent joint samples ⇒ per-bin bracket not confined to \(\pm 2\) ⇒ search can go above 2.  

Reported \(+3.95/\sqrt{N}\), 12/12 positive on triangle, ~0 on cosine: consistent with this mechanism (and with plateau vs isolated-max statistics).

### Finding 4.3 — Search bias “property of the waveform”  
**Severity: OK**  
Triangle: max on a plateau of many near-optimal quadruples ⇒ \(\mathbb{E}[\max]>\mathrm{truth}\).  
Cosine: isolated Tsirelson point ⇒ searched ≈ fixed.  
Factor ~30 in replica \(c\) and ~200 at production \(N\): sound and load-bearing for Layer 2.

### Finding 4.4 — Arm 2 weighted case  
**Severity: OK** (already caveated in the doc)  
Pointwise \(\{\pm 2\}\) argument fails with \(b\)-dependent weights; null band is required. Measured excess 0.008 vs null band ~0.13: correctly read as consistent with \(2\sqrt{2}\), detection loophole open at \(\eta\approx 0.637\).

---

## Overall verdicts

| Item | Verdict | One-line reason |
|------|---------|-----------------|
| **1. HC-1-gauged** | **SOUND WITH CORRECTIONS** | Linearisation, Gauss, coeff 12, PSD: solid. Headline wrong on high-\(\omega\) VK; \(n(H)\) is \(\ell=0\) only. |
| **2. HC-3-volume** | **SOUND WITH CORRECTIONS** | Geometry OK; VK turn real (indep. \(Q_\min\approx 86.59\) at \(\omega\approx 1.482\)); GSS reading not inverted. Overstates 20 points as clean HC-6 partition targets; near-edge junk remains after filter. |
| **3. D8b** | **SOUND WITH CORRECTIONS** | Algebra and PN factor solid; “not pinning” sound. Link flux **not** shown to fix N7’s 4.4% on actual data; `sfa_momentum` still continuum-centred. |
| **4. v87 B1 CHSH** | **SOUND** | Algebraic bound and prerequisite claim hold; arm-1 excess mechanism verified in code; waveform-dependent search bias correct. |

---

### Highest-priority corrections (if these enter census language)

1. **HC-1-gauged §3/§5:** Split \(Q_\max\)/low-\(\omega\) existence from high-\(\omega\) VK turning; do not call the whole window “not GSS.”  
2. **HC-1-gauged index:** Label \(n(H_\omega)^{(\ell=0)}\); do not claim full GSS index.  
3. **HC-3-volume targets:** Separate monochromatic VK-unstable neighborhood from genuine off-symmetric partition mismatches; tighten quality cuts beyond \(w_a<1.5\).  
4. **D8b → N7:** Recompute momentum balance with the link flux on the same SFA before closing D8.
