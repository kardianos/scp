I'll re-derive the multipole Coulomb test and stress-check the HC-6 target/control setup without modifying anything.Re-deriving the multipole Gauss/translation identities and stress-checking the HC-6 geometry next.# Adversarial review — A (hc1_gauged_l) and B (HC-6 setup)

No files modified. Analytic re-derivation + independent numerics.

---

## A. Finding 1.6 closure (`hc1_gauged_l`)

### Finding A1 — Coulomb coefficient / residual as a self-test  
**Severity: MINOR**  
**Claim:** recovering the l=1 translation zero mode “pins the Coulomb coefficient, sign and Green’s function independently of any l=0 statement.”

**What holds:** The construction is the right one. With \(u=r\phi\), multipole Gauss is
\[
K_l u_c = 2g^2\,(\mathrm{wt}\,f)\sum_a u_{x_a},\quad
K_l=-\partial_r^2+\frac{l(l+1)}{r^2}+3g^2 f^2,
\]
and symmetric elimination still produces **coefficient 12** for every \(l\). Centrifugal on \(K_l\) is required and is present in the code.

Independent continuum check on the shooter:  
\((-\Delta_1+3g^2 f^2)\chi' \stackrel{?}{=} 6g^2\,\mathrm{wt}\,f\,f'\)  
relative residual \(\sim 4\times 10^{-5}\)–\(1\times 10^{-3}\) on \(\{f>0.01 f_0\}\) — identity holds.

**What is weak:** Residual improvement alone does **not** uniquely pin 12. Any \(\alpha\in(0,12]\) that partially cancels \(L_x u\) improves the residual. Distinguisher (do this, not only residual ratios):

1. **Direct:** solve \(K_1 u_c=6g^2(\mathrm{wt}f)(r f')\) and compare \(u_c\) to \(r\chi'\) (shape + scale).  
2. **h-refinement:** only the correct \(\alpha\) drives residual/eigenvalue → discretisation floor.  
3. **Wrong-α scan:** residual vs \(\alpha\) should minimize at 12.

I measured (1): discrete eliminated \(u_c\) vs \(r\chi'\) still has **~8–10% relative error** (grid / interp / outer BC), so the self-test is supportive but not a machine-precision pin of 12. Eigenvalues are sharper: l=1 \(L_x^\mathrm{sym}\) goes \(\sim-10^{-2}\to\sim-10^{-5}\) with Coulomb.

**Corrected claim:** residual+eigenvalue tests support sign, \(K_l\), and \(\alpha\approx 12\); they do not by themselves prove uniqueness of 12 without an \(\alpha\)-scan or \(u_c\leftrightarrow r\chi'\) match under refinement.

---

### Finding A2 — Translation ⇒ \(\delta A_0\) / constrained \(c=\chi'Y_{1m}\)  
**Severity: OK**  
**Claim:** constrained elimination must reproduce \(\chi'(r)Y_{1m}\).

**Derivation (sketch):** Background \(\Delta\chi=-3g^2\,\mathrm{wt}\,f^2\). Differentiate radially and rearrange the l=1 radial operator:
\[
\bigl(-\partial_r^2-\tfrac{2}{r}\partial_r+\tfrac{2}{r^2}+3g^2 f^2\bigr)\chi'
=6g^2\,\mathrm{wt}\,f\,f',
\]
which is exactly Gauss on the symmetric translation mode \(x_a=f'\cos\theta\) (all three components). So yes: constrained \(c\) for \(x=f'\) is \(\chi'\) (times the same angular factor). That is the real content of the test, and it is analytically clean.

---

### Finding A3 — \(\omega=1.49\), ratio \(=1.0\)  
**Severity: MINOR** (not a Coulomb bug; residual metric fails)  
**Claim (worry):** ratio 1.0 at \(\omega=1.49\) might mean the Coulomb term is wrong there.

**What the numbers say (re-run):**

| ω | \|Mc u\|/\|L u\| | res with C | res no C | Lsym+Coul min |
|---|-----------------|------------|----------|----------------|
| 1.41 | ~1.0 | 7.9e-7 | 8.7e-6 | −4.4e-5 |
| 1.48 | ~0.85 | 5.5e-7 | 3.6e-6 | −2.5e-5 |
| **1.49** | **~0.24** | **1.13e-5** | **1.13e-5** | **−2.3e-5** |

At 1.49 the **eigenvalue test still works** (Coulomb lifts −6×10⁻³ → −2×10⁻⁵). The residual ratio fails because:

- Profile is thick-wall: \(r(f>0.01f_0)\approx 13.7\) (vs ~10–11 below).  
- Past the VK turn (\(\omega_\mathrm{turn}\approx 1.485\)).  
- \(u=r f'\) from `np.gradient` on an interpolated soft profile is a poor discrete null-vector; residual floor rises above the Coulomb correction.  
- \|Mc u\|/\|L u\| drops to 0.24, so cancellation is weak in max-norm even when the **spectral** projection onto the true Goldstone is fine.

**Corrected reading:** ratio 1.0 is a **failed residual diagnostic on a soft profile**, not evidence that Coulomb is off. Prefer l=1 eigenvalue under h-refinement; do not treat residual ratio as pass/fail at \(\omega\gtrsim 1.485\).

---

### Finding A4 — \(l_\max=6\) enough?  
**Severity: OK** (with a soft caveat)  
**Claim:** no negative at \(l\ge 2\) ⇒ full \(n(H_\omega)=1\).

**Argument:** For this massive theory the effective potential is compact (supported where \(f=O(1)\)). Centrifugal \(l(l+1)/r^2\) is strictly positive and, once \(\min_r l(l+1)/r^2\) exceeds the depth of the attractive well, the ground eigenvalue of each channel is forced upward in \(l\). Numerics already show, for \(l\ge 3\), \(L_0\approx L_x^\mathrm{flav}\approx L_x^\mathrm{sym}\) (potential corrections irrelevant) and **monotone increase** with \(l\). Coulomb is PSD in every multipole, so it cannot create new negatives.

Not logically airtight for every conceivable potential, but airtight **for this family** given the measured spectrum up to 6 and the compact support of \(f\). Optional cheap insurance: one point at \(l=8,10\) on the large-Q end (\(\omega=1.41\)).

**Caveat (scope):** this is still the **static Hessian** Morse index in radial multipoles, not the full dynamical BdG frequency count (Krein collisions, etc.). Same scope as HC-1-gauged.

---

### Finding A5 — \((2l+1)\) degeneracy for Morse index  
**Severity: OK**  
**Claim:** \(n(H_\omega)=\sum_l(2l+1)\,n_\mathrm{neg}(l)\).

Correct for the dimension of the negative eigenspace when each radial multipole is \((2l+1)\)-fold degenerate under SO(3). Zero modes (l=1 translation) must **not** be counted as negatives — code uses `tol=1e-7*max|ev|` so l=1 eigenvalues \(\sim-10^{-5}\) are zeros, not negatives. Only l=0 contributes one negative (VK dilational). **OK.**

---

### Finding A6 — Headline \(n(H_\omega)=1\) over l  
**Severity: OK** with prior caveats  
Log: negatives by l = `{0:[1], 1:[0], 2…6:[0]}`, \(n(H)=1\) on all six \(\omega\). Combined with A1–A5: **Finding 1.6 is closed at the level claimed**, modulo residual-test weakness at \(\omega=1.49\) and static-vs-dynamical scope.

**Verdict A: SOUND WITH CORRECTIONS**  
(close the residual-ratio rhetoric; keep eigenvalue / analytic translation test as the real pin.)

---

## B. HC-6 setup (before treating the run as a result)

### Finding B1 — Fairness: VK vs sponge vs ω  
**Severity: MAJOR**  
**Claim:** target \(n(D)=0\) vs control \(n(D)=1\) at comparable \(Q_\mathrm{tot}\) is a fair converse-GSS test.

**Confounds (stacked):**

| | Target | Control |
|--|--------|---------|
| \(w\) | (1.495, 1.495, 1.45) | (1.465, 1.465, 1.42) |
| \(\bar\omega_\mathrm{eff}\) | **1.480** | **1.450** |
| vs VK turn (~1.482) | **on/past** | **below** |
| \(r(f_\mathrm{max}>0.01)\) | **15.0** | **9.5** |
| undamped half-box | 19 (L−damp=22−3) | same |
| margin to sponge | **~4** | **~9.5** |
| \(\kappa=\sqrt{m^2-w_\mathrm{hi}^2}\) | 0.122 (w=1.495) | 0.322 (w=1.465) |

So target differs from control in **at least three** ways that all bias toward “target looks worse”:

1. **Monochromatic VK:** past the turn, even a symmetric ball has \(n(D)=0\) and is classically VK-unstable. Detuning at \(\bar\omega=1.48\) largely *pushes effective frequency over that turn* — not a pure flavored GSS novelty.  
2. **Sponge exposure:** extended thick-wall tail sits much closer to the absorber.  
3. **Mean ω / binding:** softer object, more radiation, more absorption.

**Early diag (setup red flag, not a result claim):** by \(t\sim 25\), target has already lost \(\Delta Q\approx -1.0\) (~0.9%) with **partition fractions almost frozen** (0.438→0.437); control \(\Delta Q\sim 10^{-3}\). That pattern is exactly what sponge+radiation predicts, **not** “redistribution to the sector minimum.”

**Corrected design (minimum to call the test fair):**

- Matched **mean ω** or matched **extent** (same \(r(f>0.01)\), same margin to sponge), **or**  
- Extra arm: **monochromatic** seed past the turn (\(w,w,w\) with \(n(D)=0\)) as a second control — if it dies like the target without partition drift, the failure mode is VK/box, not flavored GSS.  
- Prefer **periodic or much larger L** (e.g. L ≳ 30–35 half-domain) so both tails are deep in vacuum.

---

### Finding B2 — \(dx=0.3465\) resolution  
**Severity: MINOR** (opposite of the fear)  
**Claim / worry:** target more extended ⇒ under-resolved first.

Skin depth \(1/\kappa\): target high-\(w\) component has **more** points per skin (~24) than control (~9). Local grid resolution is **finer relative to curvature** for the target. What fails first is **box/sponge**, not dx.

**Corrected:** dx is adequate for both at this ω; do not refine dx first — enlarge L / kill sponge asymmetry first.

---

### Finding B3 — Observables  
**Severity: MINOR**  
**Claim:** \(Q_{p0,1,2}\) is the right “decay to sector minimum” readout.

**OK as primary** for flavored redistribution (kernel appends true per-component Noether charges).  
**Not sufficient alone.** Pre-register jointly:

| Observable | Role |
|------------|------|
| \(Q_a/Q_\mathrm{tot}\) | partition drift (GSS flavored signature) |
| \(Q_\mathrm{tot}(t)\) | sponge / radiation loss (must stay flat for a clean test) |
| \(E_\mathrm{total}\) | binding / dispersal |
| \(s_\mathrm{max}\), \(r_\mathrm{core}\) | collapse vs inflate |
| cluster count / fragmentation | fission channel |

If \(Q_\mathrm{tot}\) drops while fractions stay fixed → **do not** score as GSS converse success.

---

### Finding B4 — False positives to pre-register  
**Severity: MAJOR** (exclusion list)

Treat as **void / not GSS confirmation** if any of:

1. **Sponge false decay:** \(Q_\mathrm{tot}\) falls by ≳ few × control’s \(\Delta Q\), while \(\max_a|\Delta(Q_a/Q)|\) stays ≲ 10⁻³.  
2. **VK / dilational death:** radius grows, \(s_\mathrm{max}\) drops, energy radiates, **fractions frozen** — monochromatic thick-wall dispersal.  
3. **Seed projection junk:** t=0 already far from shooter \((Q,E)\) (target already \(Q=116.03\) vs profile 116.44; watch early transient).  
4. **Absorbing BC as one-way valve:** outgoing charge never returns; any radiation looks like decay.  
5. **Control also sick:** control \(\Delta Q\) or fraction drift comparable to target → void (already in your script).  
6. **Wrong n(H) assumption:** \(n(H_\omega)=1\) is monochromatic **gauged** HC-1; HC-6 is **flavoured ungauged**. If flavoured \(n(H)\neq 1\), index logic changes. Label verdicts “provisional under \(n(H)=1\)”.  
7. **Equal-component coincidence:** target has \(Q_0=Q_1\) by construction (θ=60°); “redistribution” must be specified as e.g. \(Q_2/Q_\mathrm{tot}\to 1/3\) or approach to a measured \(n(D)=1\) neighbor — not vague “components equalize.”

**Success criterion (sharpened):**  
target shows **significant partition motion** toward a lower-index region **while** \(Q_\mathrm{tot}\) conserved to the control’s floor **and** control fractions fixed → then GSS converse is supported.  
Anything else → bound or void, not confirmation.

---

### Finding B5 — Pre-registered reading in `hc6_profiles.py`  
**Severity: MINOR**  
Three-way reading (target dies / both live / both die) is right in spirit. It under-specifies **how** the target dies. Add the partition-vs-\(Q_\mathrm{tot}\) split above before locking the result.

---

## Verdicts

| Item | Verdict |
|------|---------|
| **A. hc1_gauged_l** | **SOUND WITH CORRECTIONS** — multipole ops, \(K_l\), coeff 12, translation analytic identity, \((2l+1)\), and \(n(H_\omega)^{(l\le 6)}=1\) hold. Soften residual-ratio claim; ω=1.49 ratio=1 is soft-profile floor, not a Coulomb failure; eigenvalue test still clean. |
| **B. HC-6 setup** | **GO-WITH-CHANGES** — runnable as a **tripwire / bound**, **not** as a clean converse-GSS confirmation with the present target/control pair. |

### Minimum changes before trusting a positive HC-6 reading
1. Pre-register exclusions B4 (especially sponge + frozen fractions).  
2. Add monochromatic past-turn control **or** enlarge L so both objects sit with comparable vacuum margin.  
3. Score success only on **partition drift at conserved \(Q_\mathrm{tot}\)**, not on \(Q\) loss or \(E\) drop alone.  
4. Label \(n(H)=1\) as imported from monochromatic gauged HC-1, not measured on these seeds.

With those, a careful negative or “bound only” outcome is scientifically usable; a triumphant “GSS converse confirmed” on the present pair would be **overclaimed**.
