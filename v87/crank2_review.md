I'll act as an independent adversarial reviewer: read the brief and results in order, re-run both worksheets, then attack the closed form, dimension-invariance, uniqueness, detection loophole, and fabric fit.Running both worksheets and reading the seat files plus the A6/A2/A9 theory sections.Crank 1 is still running (heavy MC). I'll derive the closed form and check the singular cases independently while it finishes.Checking a few remaining edge cases: the derivative formula cleanly, hybrid source+detection, and whether monochromatic + uniform λ really yields the triangle wave.# Adversarial review: `v87/CRANK2_RESULTS.md`

Independent seat. Worksheets re-run end-to-end; closed form, N-invariance, Fourier table, efficiencies, and hybrids recomputed here. Printed numbers that match my runs are accepted; inferences that do not follow are attacked.

---

## 1. The closed form

**Quote attacked:**  
> `dE_p/dθ = (4/Z_p)|sin θ|^{p−1} sin θ`  
> `E_p(θ) = −1 + (4/Z_p) ∫₀^θ sin^p u du`  
> `Z_p = 2√π Γ((p+1)/2)/Γ(p/2+1)`  
> and “g is odd about φ=π/2 so the full-circle integral vanishes.”

### Finding 1.1 — Closed form is correct **[D]**

Independent derivation:

1. \(Z_p = \int_0^{2\pi}|\cos\phi|^p\,d\phi = 4\int_0^{\pi/2}\cos^p\phi\,d\phi = 2\,B\bigl(\tfrac{p+1}{2},\tfrac12\bigr) = 2\sqrt{\pi}\,\Gamma\bigl(\tfrac{p+1}{2}\bigr)/\Gamma\bigl(\tfrac p2+1\bigr)\).  
   Numerics: error \(\lesssim 10^{-8}\) for \(p\in\{0,0.5,1,1.5,2\}\).

2. Write \(g(\phi)=|\cos\phi|^p\,\mathrm{sgn}(\cos\phi)\). Then  
   \(E(\theta)=-\frac1{Z}\int_0^{2\pi} g(\phi)\,\mathrm{sgn}(\cos(\phi-\theta))\,d\phi\).

3. On the half-circle where \(\cos(\phi-\theta)>0\), the sign is \(+1\). If \(\int g=0\),  
   \(E(\theta)=-\frac2Z\int_{\theta-\pi/2}^{\theta+\pi/2} g\).

4. Leibniz:  
   \(g(\theta+\pi/2)-g(\theta-\pi/2)=-2|\sin\theta|^{p-1}\sin\theta\),  
   hence \(\frac{dE}{d\theta}=\frac4Z|\sin\theta|^{p-1}\sin\theta\).

5. Integrate from \(E(0)=-1\): the incomplete-sine form follows. Landmarks:  
   \(p=0\to(2\theta-\pi)/\pi\), \(p=1\to-\cos\theta\), \(p=2\to(2\theta-\sin 2\theta-\pi)/\pi\). **Match.**

FD of the **closed-form** \(E\) recovers the analytic derivative to \(\sim10^{-8}\). (FD of the **sgn-integrand** trapezoid is biased by grid staircasing; that is a numerics artifact, not a formula bug.)

### Finding 1.2 — Oddness and integrability for \(0<p<1\) **[D]**

- **Oddness:** \(g(\pi/2+t)=-g(\pi/2-t)\) holds for every \(p\) where \(g\) is defined (algebra of absolute value and sign). Full-circle \(\int g\) is numerically \(\sim10^{-6}\) (grid residual), consistent with exact zero.
- **Singularity:** near \(\cos\phi=0\), \(|g|\sim|\cos|^p\). Integrable iff \(p>-1\). So for all \(p>0\), including \(0<p<1\), \(g\in L^1(\mathbb{S}^1)\).
- **Endpoint product:** \(|\sin|^{p-1}\sin\to0\) as \(\theta\to0^+\) for \(p>0\), so \(\frac{dE}{d\theta}\) extends continuously by \(0\) at the endpoints even though the bare factor \(|\sin|^{p-1}\) diverges for \(p<1\).
- Leibniz for an \(L^1\) integrand with jump endpoints is standard; no extra caveat is required beyond \(p>-1\).

**Verdict on §1:** The closed form stands. The singularity worry does not break the argument.

---

## 2. Dimension-invariance

**Quote attacked:**  
> “E(θ) is EXACTLY independent of the dimension N … the radial integral cancels.”

### Finding 2.1 — True for the stated model, including \(N=2\) **[D,M]**

For \(N\ge3\), with \(A=\mathrm{sgn}(a\cdot\lambda)\), \(B=-\mathrm{sgn}(b\cdot\lambda)\), \(\rho\propto|\lambda\cdot b|^p\), and \(a,b\) coplanar: the integrand factors as (radial)×(angular). Radial pieces cancel. Worksheet MC: \(N=2\ldots25\), \(p=1\), \(\theta=60^\circ\), \(|E+\cos\theta|\le4\times10^{-4}\). My MC agrees.

For **\(N=2\)** the writeup’s marginal formula \((1-r^2)^{(N-4)/2}\) is for \(N\ge3\), but the claim survives by degeneration: \(\lambda\) lives on \(S^1\), \(r\equiv1\), and \(E\) is exactly the angular integral that every higher \(N\) reduces to. MC: \(N=2,3,4\) agree for \(p=0,1,2\).

### Finding 2.2 — Invariance is an artifact of 2-plane responses, not “physics of any hidden-variable body” **[D,M]**

Replace \(\mathrm{sgn}\) by soft responses \(A=\tanh(k\,a\cdot\lambda)\), \(B=-\tanh(k\,b\cdot\lambda)\), \(k=3\), \(p=1\), \(\theta=60^\circ\):

| \(N\) | \(E\) (soft) |
|------|----------------|
| 2 | −0.470 |
| 3 | −0.440 |
| 5 | −0.392 |
| 10 | −0.313 |

Strong \(N\)-dependence. Same for \(p=0\).

**Scope (exact):**  
\(E(\theta)\) is \(N\)-independent **if and only if** both the weight and the product of responses depend on \(\lambda\) only through the **angular coordinate in the plane of \((a,b)\)** (up to a common radial factor that cancels). Hemisphere/sgn models satisfy that. General dichotomic or continuous responses that retain radial magnitude do **not**.

Crank 2 inherits crank 1’s claim without stating this restriction. That is an overclaim of generality, not a false computation inside the model.

---

## 3. Uniqueness / monochromaticity

**Quote attacked:**  
> “p = 1 is the unique monochromatic rung.”  
> “p = 1 ⟺ single harmonic ⟺ E=−cos θ ⟺ S_max=2√2.”  
> “A single rotating phase, projected on a setting, yields a single harmonic… p=1 is what a monochromatic object gives.”

### Finding 3.1 — Fourier table verified **[M]**

Re-run of `why_p1.py` part C matches the printed coefficients to the digits shown; at \(p=1\), \(\|\mathrm{higher}\|/c_1\sim2\times10^{-10}\).

Within the **p-family only**, uniqueness is sharp **[D]**:  
on \((0,\pi)\), \(\frac{dE}{d\theta}=\frac4{Z_p}\sin^p\theta\). Equality with \(\sin\theta\) (the derivative of \(-\cos\theta\)) for all \(\theta\) forces \(p=1\) and \(Z_1=4\).

### Finding 3.2 — The fabric inference is a pun on “monochromatic” **[D]** (pun) / **[C]** (false bridge)

Two different meanings:

| Sense | Object | What is single? |
|--------|--------|------------------|
| **Temporal** | \(\Phi_a=f(r)e^{i\omega t}\); HC-1 modes | one **clock frequency** |
| **Angular** | \(E_p(\theta)=\sum c_k\cos(k\theta)\) | one **harmonic in setting angle** |

They are not the same.

**Counterexample (decisive):** Take a monochromatic rotator \(\lambda\) (one phase), responses \(A=\mathrm{sgn}\cos(\lambda-a)\), \(B=-\mathrm{sgn}\cos(\lambda-b)\), **uniform** sampling \(p=0\). Object is temporally monochromatic. Correlation is the **triangle wave** \(E=-1+2\theta/\pi\), with \(c_1\approx0.81\), large odd harmonics, \(S_{\max}=2\) — **not** the cosine rung.

So:

- “Monochromatic object ⇒ \(p=1\)” is **false**.  
- “Monochromatic object + uniform λ ⇒ classical triangle” is what the model actually gives.  
- What selects the single **angular** harmonic is the **sampling weight** \(|\lambda\cdot b|\) (or an equivalent MD tilt), not the single temporal clock.

The rectifier identity \(|x|\mathrm{sgn}(x)=x\) at \(p=1\) **[D]** is a real mechanism for why that weight reconstitutes a linear projection. It does **not** license “HC-1 monochromaticity ⇒ \(p=1\)”.

Section 5.2.3 half-admits this (“monochromaticity is a property of the object, not a derivation of the sampling weight”) while §3 and the status table still sell “What singles out \(p=1\)? **Monochromaticity** **[D]**”. That is inconsistent self-marketing.

---

## 4. Detection-loophole fork

**Quote attacked:**  
> Reading 2 is the detection loophole and is refuted because \(\langle|\cos|\rangle=2/\pi\), \(\langle|\lambda\cdot b|\rangle=1/2\), both \(\ll 2(\sqrt2-1)\approx0.828\).

### Finding 4.1 — Efficiencies and threshold numbers **[D,M]**

| Quantity | Document | Independent |
|----------|----------|-------------|
| \(\langle|\cos\phi|\rangle\) on circle | 0.63662 | \(2/\pi\) exact; MC matches |
| \(\langle|\lambda\cdot b|\rangle\) on \(S^2\) | 0.50010 | \(1/2\) exact; MC 0.4999 |
| \(2(\sqrt2-1)\) | 0.82843 | exact; \(=2/(1+\sqrt2)\) |

The CHSH critical efficiency \(\eta>2(\sqrt2-1)\) for symmetric detectors / maximally entangled target with no-detection assigned in the CHSH bookkeeping is the **standard CHSH detection-loophole number** (Garg–Mermin / related CHSH analyses). **[UNVERIFIED against primary paper this session; identity and algebra of the threshold are standard.]**  
Clausers–Horne / Eberhard thresholds are lower (\(\sim2/3\)); pure Reading 2 still sits below those too (\(2/\pi\approx0.637<2/3\), \(1/2<2/3\)).

### Finding 4.2 — Pure Reading 2 is correctly killed **[D]**

If detection probability is \(|\lambda\cdot b|\) (max-normalized to 1), mean efficiency is ½ (sphere) or \(2/\pi\) (phase). Loophole-free runs report high η **and** \(S>2\). A model whose intrinsic fire rate is ~50–64% cannot be the data-generating process for η≳90% experiments. Identification of pure Reading 2 with a detection-loophole / fair-sampling model is correct in substance.

### Finding 4.3 — Hybrids weaken the cleanliness of the fork **[D,M]**

Document presents a binary: source MD (alive) vs detector sampling (dead). Missing case: **base efficiency + alignment-dependent boost**,

\[
P(\mathrm{detect}\mid\lambda,b)=c+(1-c)\,|\lambda\cdot b|.
\]

My MC (B-side only; Tsirelson angles; post-select):

| \(c\) | \(\eta_B\) | \(S\) |
|------|------------|-------|
| 0.00 | 0.50 | ≈2.83 |
| 0.70 | 0.85 | ≈2.15 |
| 0.85 | 0.93 | ≈2.07 |
| 1.00 | 1.00 | 2.00 |

Both detectors biased: at coincidence η≈0.90 still \(S\approx2.09>2\).

So:

- Pure \(|\lambda\cdot b|\) detection is dead as an explanation of **high-η, near-Tsirelson** data.  
- Mild CHSH violation at high η via **partial** detection bias is still a local loophole model; the 0.828 threshold is about reaching the **quantum** point with ideal visibility, not about forbidding every \(S>2\).  
- A hybrid (partial source MD + partial detection bias) is **not** refuted by §4; source MD is already the surviving reading, so hybrids do not reopen a closed door so much as show the fork is less binary than written.

**§4 is directionally right and load-bearing; it overstates cleanliness.**

---

## 5. Fabric fit (A6)

**Quote attacked:**  
> A6 force \(\propto\cos(\Delta\phi)=|\cos\Delta\phi|\cdot\mathrm{sgn}(\cos\Delta\phi)\) “is precisely weight × response at \(p=1\). That is not a fitted analogy; it is A6 read literally.”  
> Table: “weight = magnitude… **[M]**”, “Does the fabric have the right algebraic form? **Yes, measured**.”

### Finding 5.1 — Force ≠ dichotomic correlation **[D]**

Document §5.2.1 already states this. It is fatal to the “literally” rhetoric: the measured object is a **near-field force** between solitons, not \(E(a,b)=\langle A B\rangle\) for spacelike settings. Same algebraic skeleton, different observable class. Tagging the Bell-side identification **[M]** is false; A6 is **[M]**, the map to CHSH weight×response is **[C]** or **[P]**.

### Finding 5.2 — Factorization is a tautology, not evidence **[D]**

For any real number \(x\), \(x=|x|\,\mathrm{sgn}(x)\). Writing \(\cos=|\cos|\mathrm{sgn}(\cos)\) does not show that nature samples with weight \(|\cos|\) and reports \(\mathrm{sgn}(\cos)\). Any force \(\propto f(\Delta\phi)\) admits the same split. This is reverse-engineered bookkeeping, not a measured sampling architecture.

### Finding 5.3 — \(\Delta\phi\) is not \(\lambda\cdot b\) **[D]** (category error)

| | \(\Delta\phi\) (A6) | \(\lambda\cdot b\) (Bell model) |
|--|---------------------|----------------------------------|
| What | Relative **dynamical** clock phase of two field lumps | Hidden variable dotted with **setting** |
| Control | Evolves under PDE; not a free knob | \(b\) is (ideally) free / chooser output |
| Role | Argument of a force law | Argument of a conditional measure \(\rho(\lambda\mid b)\) |

Substituting a dynamical relative phase for a freely chosen analyzer direction is the central illicit step. O2 makes choosers fabric, so settings are *also* dynamical — but then one must **construct** the chooser functional and show the joint \((\lambda,b)\) law; one may not silently rename \(\Delta\phi\leftrightarrow\lambda\cdot b\).

### Finding 5.4 — No detector setting exists in the kernel yet **[C/gap]**

For a “setting” to exist one needs at least:

1. a fabric degree of freedom that plays the role of analyzer orientation / phase reference;  
2. quasi-independent choice (or controlled preparation) of that DOF on Alice and Bob;  
3. a dichotomic readout map \(A(a,\lambda)\in\{\pm1\}\);  
4. spacelike separation with local response.

None of this is in A6. A6 is contact physics at small \(D\). Exponential death of contact (\(e^{-\kappa D}\sim2\times10^{-3}\) at \(D=12\), document’s own numbers, consistent with κ≈0.5117) is correctly stressed; it blocks using A6 *as* the Bell channel and forces common-past architecture — which returns to GROK G14 / seats’ “capacity easy, structure hard,” not to a one-line A6 reading.

### Finding 5.5 — Exponential gap and G14 correctly flagged **[M/D]**

§5.2.2 is sound: direct locking is numerically dead at Bell scales; correlation must be prepared in the past and survive transport; kernel mixing works against frozen MD. This is one of the document’s honest strengths.

---

## 6. Overclaims and omissions

### Finding 6.1 — Overclaims

1. **Status table:** “What singles out \(p=1\)? Monochromaticity **[D]**” — Fourier uniqueness inside the p-family is **[D]**; fabric monochromaticity link is **[C]** and, as argued, wrong as an implication.  
2. **“Fabric has the right algebraic form? Yes, measured **[M]**”** — measures a force, not Bell \(E(a,b)\).  
3. **“Not a fitted analogy; A6 read literally”** — tautological factorization + category error (Finding 5.2–5.3).  
4. **“The gap has narrowed… to a single, well-posed question rather than a research programme”** — rhetorical compression (Finding 8).  
5. **N-independence** stated without the sgn/2-plane scope restriction.  
6. **Detection fork** presented as fully decisive without hybrids / CH vs CHSH threshold nuance.

### Finding 6.2 — What it should say and does not

1. Explicit scope: N-invariance holds only for responses/weights that factor through the \(a\)–\(b\) plane angle.  
2. Explicit **pun warning** on monochromatic (temporal vs angular).  
3. That \(p=0\) monochromatic already gives triangle / \(S=2\); MD weight is the quantum-selecting structure.  
4. That other MD models (not power weights) also achieve \(E=-\cos\theta\) (Hall, Degorre et al. style); \(p=1\) is not unique outside the family.  
5. Hybrid detection models; distinguish “cannot explain loophole-free Tsirelson” from “cannot produce any \(S>2\)”.  
6. That the Tsirelson half remains open: pure MD dials to 4 (seats already proved this); crank 2 does not add a bound principle — it only prefers the cosine rung *inside an already postulated family*.  
7. Operational definition of a fabric **setting** — absent.  
8. Confusion risk between A6 **Coulomb/gauge** (long-range, internal-state-blind) and A6 **phase-coherent contact** (short-range, \(\cos\Delta\phi\)) — see §7.

What it **does** well and should keep: closed form; rectifier at \(p=1\); pure detection reading killed; exponential suppression numbers; source-side design constraint; honest §5.2 gaps list (if the marketing around it were deleted).

---

## 7. Closing chain in §5.3

**Quote attacked:**  
> gauge couples to \(\rho_Q=\mathrm{Im}(\bar\Phi\dot\Phi)\), bilinear ⇒ rate ∝ projection ⇒ \(p=1\).

### Finding 7.1 — Non-sequitur **[D]** (logic) / **[C]** (as research hope)

Link-by-link:

1. **\(\rho_Q\) bilinear** — true as a local expression **[M/kernel]**. For a monochromatic Q-ball, \(\Phi=fe^{i\omega t}\) gives \(\rho_Q=\omega|\Phi|^2\): a **scalar charge density**, not a projection onto an external setting direction.

2. **Gauge sector vs contact** — THEORY A6 splits:  
   - Coulomb / gauge: long-range, **internal-state-blind**;  
   - phase-coherent contact: \(\propto\cos(\Delta\phi)e^{-\kappa D}\).  
   The \(\cos\Delta\phi\) structure crank 2 wants lives in **contact**, not in the gauge bilinear. Citing \(\rho_Q\) as the source of \(p=1\) **points at the wrong A6 force**.

3. **“Bilinear ⇒ sampling rate ∝ projection”** — does not follow. Bilinearity of a charge density does not determine the measure \(\rho(\lambda\mid b)\) over hidden variables / phases conditioned on chooser state.

4. **Source-side requirement** — correctly restated from §4, but not advanced.

**What would actually establish \(p=1\):**

- Define settings \(a,b\) as fabric functionals (chooser subsystems).  
- From the joint IVP (source + two choosers + common past), compute or prove  
  \(\rho(\lambda\mid b)\propto|\lambda\cdot b|\) (or an equivalent tilt that yields \(E=-\cos\theta\)), **at production**, not at detection.  
- Show the tilt **survives** free evolution / scrambling to spacelike measurement (addresses G14).  
- Map continuous fabric readouts to dichotomic \(A,B\) without reintroducing a detection loophole.  
- Show the dynamics **forbid** nearby rungs that give \(S>2\sqrt2\) (or accept that Tsirelson is still imported).

That is a multi-step derivation programme, not a one-liner from bilinearity.

---

## 8. “The gap is now one exponent wide”

**Quote attacked:**  
> gap narrowed to “derive the \(|\lambda\cdot b|\) sampling weight, source-side, from the gauge bilinear… one exponent wide.”

### Finding 8.1 — Rhetorical compression of a much larger gap **[D]** as assessment

Remaining walls, none of which is “just fix \(p\)”:

| # | Gap | Status |
|---|-----|--------|
| 1 | What is a **setting** in the kernel? | Undefined |
| 2 | Source-side \(\rho(\lambda\mid a,b)\) of quantum form | Postulated |
| 3 | Survival under transport / mixing (G14) | Open obstruction |
| 4 | Continuous force/phase → dichotomic CHSH outcomes | Unbuilt |
| 5 | Why **exactly** \(p=1\) (not 0.9 or 1.1) from dynamics | Open |
| 6 | Why not \(S\to4\) (Tsirelson principle) | Still imported |
| 7 | Spacelike separation with local responses | Only by construction of MD |
| 8 | Gauge vs contact mis-ID in the proposed closing chain | Wrong target coupling |

Calling this “one exponent wide” is **not honest**. The seats already had capacity \(\ge0.046\) bits and a working \(p=1\) model. Crank 2 usefully (i) closed the pure detection reading, (ii) gave a clean ladder/closed form, (iii) isolated the rectifier identity. It did **not** reduce the derivation wall to a single exponent; it **localized one algebraic preference inside a postulated family** and then oversold that localization as the residual gap.

---

## Cross-checks run (this review)

- `nd_shape_limit.py`: exit 0; N-independence MC, waveform ladder, CHSH ceilings (\(p=1\to2.82845\)), chained Bell ratios — match document.  
- `why_p1.py`: exit 0; landmarks, Fourier table, efficiencies — match.  
- Independent SymPy/NumPy: \(Z_p\), \(E_p\), \(dE_p\), \(N=2\), soft-response N-break, hybrid detection, monochromatic+\(p=0\) triangle.

---

## VERDICT

### **SOUND WITH CORRECTIONS**

**Keep (sound):**

- Closed form of the p-ladder and landmarks **[D]**.  
- Rectifier mechanism \(|x|\mathrm{sgn}(x)=x\) at \(p=1\) **[D]**.  
- Fourier uniqueness of \(p=1\) **within the p-family** **[D,M]**.  
- Pure detector-side \(|\lambda\cdot b|\) reading incompatible with loophole-free efficiencies **[D,M]** (as a hard design constraint against pure Reading 2).  
- Exponential death of A6 contact at Bell scales **[M]**; common-past requirement.  
- Explicit admission in §5.2 that force ≠ correlation and monochromaticity ≠ derived weight (if elevated over the marketing).

**Corrections required (without these, several central claims are not sound):**

1. **Scope N-invariance** to sgn / 2-plane factorizing responses; state failure for radial-sensitive responses.  
2. **Retract** “monochromatic object ⇒ \(p=1\)”. Call the temporal/angular clash a **pun**. Tag fabric monochromaticity link **[C]** or delete it as a derivation.  
3. **Demote** “A6 is literally weight×response **[M]**” to “algebraic analogy **[C]**”; factorization is tautological.  
4. **Do not identify** \(\Delta\phi\) with \(\lambda\cdot b\) without a constructed setting.  
5. **Rewrite §5.3**: gauge/\(\rho_Q\) chain is a **non-sequitur** and targets the wrong A6 sector; replace with a real derivation checklist.  
6. **Qualify §4**: pure Reading 2 dead; hybrids and “any \(S>2\)” vs “Tsirelson at high η” distinguished.  
7. **Retract** “gap is one exponent wide” / “single well-posed question rather than a research programme.” Residual gap remains multi-dimensional (settings, source MD derivation, transport, dichotomic readout, Tsirelson bound).  
8. Status table: “fabric has right form? **Yes [M]**” → **No as Bell correlator; yes as force algebra [M/C]**.

**What would make the document NOT SOUND overall:** treating §3 monochromaticity + §5.1 “literally A6” + §5.3 gauge chain + “one exponent” as established results rather than conjectural framing. As geometry of a postulated MD family plus a useful experimental fork, the core math holds; as a claim that the fabric gap has collapsed to deriving one exponent from \(\rho_Q\), it does not.
