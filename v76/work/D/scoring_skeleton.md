# Approach D — Scoring / Inversion Skeleton (Round 1)

**Status:** design + minimal demo (aligned with C/A FOR_D)  
**Code:** `invert_demo.py`, `gen_results_data.py`  
**Not claimed:** recovery of Einstein, real lensing, or galactic DM.

---

## 1. Data (what we invert from)

Synthetic observables from a **hidden monist path-cost kernel** (C Class C —
compact lock + nonlocal free response). Primary target is **ℓ / rays**, not
\(\delta\rho_{\mathrm{free}}\) (A-008 / C-001 FOR_D).

| Datum | Symbol | Role |
|-------|--------|------|
| Impact-parameter ray delays | \(\delta t(b)\) | Born path-time excess from \(\delta\ell\) |
| Impact-parameter deflections | \(\alpha(b)\) | Born bend from \(\partial_y\delta\ell\) |
| Bound mass ledger | \(M = E_\star/c^2\) | Unlock/bound ledger (known in synthetic) |
| Path-cost profile samples | \(\delta\ell(r)\) | Direct C3 target \(\sim M/r\) class |

Rays: straight-line Born (impact \(b\)). Curved rays Round 2+.

---

## 2. Medium parameters (what we fit)

### 2.1 Monist kernel (preferred Class C — **one sector**)

Compact lock ledger \(M\); free-response path-cost excess (medium law, not T→G):

\[
\delta\ell(r) = \frac{\beta\, M}{\sqrt{r^2+\varepsilon^2}}.
\]

Fit: \(\{M,\beta,\varepsilon\}\). Local free \(c\) fixed; warp = chart expression
of path cost. **No separate gravity field.**

### 2.2 Local-optics monist (wrong class for long-range — C no-go)

\[
\rho_{\mathrm{bound}}=A e^{-r^2/(2\sigma^2)},\quad
\rho_{\mathrm{free}}=\rho_0-\rho_{\mathrm{bound}},\quad
n-1=\kappa\,\rho_{\mathrm{bound}}/\rho_0.
\]

Fit: \(\{A,\sigma,\kappa\}\). One sector + budget link, but **compact** ray
signature — expected to lose \(L_{\mathrm{fit}}\) on kernel truth.

### 2.3 Dualist adversarial baseline (two sectors — C strip S2/S6/S9)

1. Matter \(\rho_m\) on fixed stage (mass \(= \int\rho_m\)).
2. Independent Plummer \(\Phi = -G_{\mathrm{eff}} M_{\mathrm{tot}}/\sqrt{r^2+a^2}\).
3. Observables from \(\Phi\), not free–bound identity.

**Important Round-1 fact:** Plummer Born rays are **algebra-isomorphic** to the
monist kernel when \(G_{\mathrm{eff}}=\beta\), \(a=\varepsilon\),
\(M_{\mathrm{tot}}=M\). Pure \(L_{\mathrm{fit}}\) can **tie**. Discrimination is
**ontology Occam** (D3), not fit alone. That is intentional pedagogy — dualism
can fake the map; the second sector still costs.

Parameters: \(\{A_m,\sigma,G_{\mathrm{eff}}\}\).

---

## 3. Loss / score

### 3.1 Fit term

\[
L_{\mathrm{fit}}
=
w_t\,\mathrm{MSE}(\delta t)
+
w_\alpha\,\mathrm{MSE}(\alpha)
+
w_M\big((M^{\mathrm{pred}}-M^{\mathrm{data}})/M^{\mathrm{data}}\big)^2.
\]

### 3.2 Ontology Occam (penalize second sectors)

\[
L_{\mathrm{occ}}
=
\lambda_{\mathrm{sec}}\,(N_{\mathrm{sectors}}-1)
+
\lambda_{\mathrm{link}}\,\mathbf{1}\{\text{no free–bound link}\}.
\]

| Model | \(N_{\mathrm{sectors}}\) | free–bound link |
|-------|--------------------------|-----------------|
| Monist kernel | 1 | yes (ledger \(M=E_\star/c^2\)) |
| Local optics | 1 | yes (pointwise budget) |
| Dualist Plummer | 2 | no |

### 3.3 Extra penalties

\[
L_{\mathrm{xtra}}
=
\lambda_{\mathrm{softE}}\,\mathbf{1}\{\text{Phi/Einstein labeled monist}\}
+
\lambda_{\mathrm{nogo}}\,\mathbf{1}\{\text{local budget claims long-range }\rho\sim 1/r\}.
\]

### 3.4 Combined score

\[
S = L_{\mathrm{fit}} + L_{\mathrm{occ}} + L_{\mathrm{xtra}}.
\]

**Round-1 outcome (zero noise):**  
\(S_{\mathrm{MK}}=0\), \(S_{\mathrm{LO}}\approx 0.19\), \(S_{\mathrm{dual}}=1.5\),
\(S_{\mathrm{softE}}=100\). **score_winner = monist_kernel.**  
Pure-fit MK ties dualist; Occam decides.

---

## 4. Inversion algorithm

1. Generate kernel truth → rays + \(M\).
2. Fit monist kernel, local optics, dualist (pattern search multi-start).
3. Score all; ablate \(\lambda_{\mathrm{sec}}\).
4. Write TSV/JSON/SVG under `results/`.

---

## 5. Explicit non-goals (Round 1)

- Soft-coding Einstein into monist map as “success.”
- Real GR / galactic D5.
- `scp_sim` / Q-balls as monist proof.

---

## 6. Deliverables

| Artifact | Path |
|----------|------|
| Skeleton | `work/D/scoring_skeleton.md` |
| Demo | `work/D/invert_demo.py` |
| Offline tables | `work/D/gen_results_data.py` |
| Results | `work/D/results/*` |
| Log | `logs/D_reverse_numeric.log` |
