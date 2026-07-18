# Forced Free-Budget / Path-Cost Profile Sketch (C3) — v0

**Approach:** C — reverse theoretical  
**Date:** 2026-07-18  
**Status:** approximate reverse sketch; targets not dynamics  
**Rule:** Use weak-field light delay/deflection as **phenomenological targets**. Do not take Einstein’s equation as a monist axiom. Derive what free-medium bookkeeping must look like if  
**warp = constant local \(c\) around locks.**

---

## 1. Operational setup

### 1.1 Local law (seed)

In every local frame constructed from free-medium rods/clocks:

\[
\text{free-signal speed} = c = \mathrm{const},\quad \text{isotropic}.
\]

### 1.2 Global chart (lab / asymptotic free chart)

Choose coordinates that look Minkowski in the asymptotic free region far from locks. In that chart, free signals need **not** travel in coordinate-straight lines at coordinate speed \(c\). The mismatch is the warp bookkeeping.

### 1.3 Isolated lock

A single compact lock with

\[
M \;=\; \frac{E_\star}{c^2},
\qquad
E_\star \;=\; E_{\mathrm{unlock}} \;=\; E_{\mathrm{rest}}.
\]

Spherical symmetry, static, weak field: dimensionless potential strength \(\varepsilon \ll 1\).

---

## 2. What must be reproduced (fit targets)

Denote an emergent coupling \(G_{\mathrm{eff}}\) with dimensions of Newton’s \(G\). Prefer coefficients that match classical weak-field light tests when \(G_{\mathrm{eff}}\to G\):

**Deflection** (grazing ray, impact parameter \(b\)):

\[
\Delta\theta \;=\; \frac{4 G_{\mathrm{eff}} M}{c^2 b}
\quad\text{(Einstein target; Newtonian light would give half)}.
\]

**Shapiro-like one-way excess** (schematic logarithmic form):

\[
\Delta t \;=\; \frac{2 G_{\mathrm{eff}} M}{c^3}
\ln\!\left(\frac{4 r_E r_R}{b^2}\right)
\;+\; \cdots
\]

**Interpretation under constant local \(c\):** these are not “gravity acting on photons as particles on a stage.” They are excess free-path costs and direction changes of free medium signals when global charts are used.

---

## 3. Path-cost language (primary monist object)

Define the **excess optical path** of a free ray \(\gamma\) relative to the uniform free background:

\[
\Delta\mathcal{L}[\gamma]
\;=\;
\int_\gamma \bigl(\ell[\mathcal{M}] - \ell_0\bigr)\,ds_{\mathrm{chart}},
\]

where \(\ell[\mathcal{M}]\) is local path-cost density (inverse free capacity, index, relational cost — model class below), and \(\ell_0\) is the asymptotic free value.

### 3.1 Forced weak-field monopole

For an isolated lock, the leading isotropic excess cost density in the exterior chart must behave as

\[
\boxed{
\ell(r) - \ell_0
\;=\;
\frac{2 G_{\mathrm{eff}} M}{c^2\, r}
\;+\;
\mathcal{O}\!\left(\frac{1}{r^2}\right)
\;+\;
\text{(near-zone lock structure)}
}
\]

so that:

- \(\int (\ell-\ell_0)\,dl\) along a straight-line approximation yields Shapiro-class logarithms;
- transverse gradient \(\nabla_\perp \ell\) yields deflection \(\sim 4 G_{\mathrm{eff}} M/(c^2 b)\) when both “time” and “space” pieces of path cost are present at this order (see §5).

**Units:** \(\ell\) dimensionless if \(ds\) is coordinate length and \(\Delta\mathcal{L}\) is excess coordinate length; excess time \(\Delta t = \Delta\mathcal{L}/c\).

### 3.2 Slogan form

> **Forced path-cost profile (weak field): \(1/r\) monopole with amplitude \(2 G_{\mathrm{eff}} M/c^2\).**  
> Any forward medium must reproduce this **path-cost** tail (or an observationally equivalent nonlocal functional), not necessarily a \(1/r\) free-**energy-density** tail.

---

## 4. Model class A — Local refractive index (optical monism)

### 4.1 Ansatz

In the asymptotic chart, free waves obey an eikonal with coordinate speed \(c/n(r)\), while local frames renormalize so proper free speed is \(c\):

\[
n(r) \;=\; 1 + \delta n(r),
\qquad
|\delta n|\ll 1.
\]

Excess path density \(\ell - \ell_0 = \delta n\) (to first order for isotropic index).

### 4.2 Forced index for Shapiro-class delay

Matching \(\delta n\) to §3.1:

\[
\boxed{
\delta n(r)
\;\approx\;
\frac{2 G_{\mathrm{eff}} M}{c^2\, r}
\quad (r \gg \text{lock radius})
}
\]

(Full Einstein light deflection typically needs an effective spatial contribution of the same order — often written as isotropic \(n \approx 1 - 2\Phi/c^2\) with \(\Phi = -G_{\mathrm{eff}}M/r\), i.e. \(\delta n = 2 G_{\mathrm{eff}}M/(c^2 r)\).)

### 4.3 Link to free density (if local)

Suppose a linear response of index to free-budget contrast:

\[
\delta n
\;=\;
\beta\,
\frac{\rho_0 - \rho_{\mathrm{free}}}{\rho_0}
\;=\;
\beta\,
\frac{\delta\rho_{\mathrm{free}}}{\rho_0},
\qquad
\beta = \mathcal{O}(1)\ \text{medium constant}.
\]

Then the **forced free-budget contrast under local optics** is

\[
\boxed{
\frac{\delta\rho_{\mathrm{free}}(r)}{\rho_0}
\;\approx\;
\frac{2 G_{\mathrm{eff}} M}{\beta\, c^2\, r}
}
\]

i.e. **exterior free-density depletion \(\propto 1/r\)** if one insists that path cost is local in \(\rho_{\mathrm{free}}\).

### 4.4 No-go / tension (must be faced)

Local budget identity

\[
\rho_{\mathrm{free}} + \rho_{\mathrm{bound}} = \rho_0
\]

with compact \(\rho_{\mathrm{bound}}\) implies compact \(\delta\rho_{\mathrm{free}}\) if the identity is **pointwise** and \(\rho_{\mathrm{bound}}=0\) outside the lock. That **contradicts** §4.3’s long-range \(1/r\) free-density tail.

Even if one **allows** long-range \(\delta\rho_{\mathrm{free}}\propto 1/r\),

\[
\int_{R}^{\infty} 4\pi r^2 \delta\rho_{\mathrm{free}}(r)\,dr
\;\propto\;
\int_R^\infty r\,dr
\;=\; \infty
\]

so the free-energy deficit is **not** finite and cannot equal \(E_\star\).

**Conclusion for Class A:**  
Local \(n(\rho_{\mathrm{free}})\) + pointwise budget + finite \(E_\star\) **cannot** all hold if weak lensing requires long-range \(\delta n\propto 1/r\).

Forward theories must break at least one:

| Escape hatch | Meaning |
|--------------|---------|
| **A-escape-1** | Path cost nonlocal: compact free deficit, long-range \(\ell\) via constraints/Green functions |
| **A-escape-2** | Extended free rearrangement not equal to bound support (history halo; budget integral still finite only if \(\delta\rho\) decays fast enough — then \(\delta n\) cannot stay purely local-\(\rho\)) |
| **A-escape-3** | Budget only integral/topological; local \(\rho_f+\rho_b=\rho_0\) dropped |
| **A-escape-4** | Index depends on gradients / strain / free-path anisotropy, not only \(\rho_{\mathrm{free}}\) |

**Reverse recommendation:** treat **path cost \(\ell(r)\)** as the primary forced profile; treat \(\rho_{\mathrm{free}}(r)\) as secondary and model-dependent.

---

## 5. Model class B — Split “time” and “space” path costs

Weak-field isotropic metrics that fit light tests often look like:

\[
ds^2
\;\approx\;
-\,c^2\bigl(1+2\Phi/c^2\bigr)\,dt^2
\;+\;
\bigl(1-2\Phi/c^2\bigr)\,d\mathbf{x}^2,
\qquad
\Phi = -\frac{G_{\mathrm{eff}} M}{r}.
\]

**Caution:** writing \(g_{\mu\nu}\) is allowed as **diagnostic bookkeeping** for free-signal structure, not as an independent Einstein-sourced field.

For null rays, both \(g_{00}\) and \(g_{ij}\) contribute ~ equally to deflection (factor 2 each → total 4).

Monist reading:

- \(g_{00}\) piece ↔ free-clock / free-update rate contrast in the chart;
- \(g_{ij}\) piece ↔ free-rod / free-path spatial capacity contrast;

both must arise from **one** medium state around the lock.

**Forced sketch:**

\[
\boxed{
\frac{\delta\tau_{\mathrm{clock}}}{\tau}
\sim
\frac{|\Phi|}{c^2}
=
\frac{G_{\mathrm{eff}} M}{c^2 r},
\qquad
\frac{\delta L_{\mathrm{rod}}}{L}
\sim
\frac{G_{\mathrm{eff}} M}{c^2 r}
}
\]

in the weak field — same \(1/r\) class, same \(M\), operationally from free medium.

---

## 6. Model class C — Relational path cost (preferred monist shape)

Avoid local \(n(\rho)\) as fundamental. Define free influence relation: event \(p\) can affect \(q\) if there is a free-medium chain with update cost \(\le T\).

**Sketch functional (illustrative, not unique dynamics):**

\[
\Delta\mathcal{L}[\gamma]
\;=\;
\alpha
\int_\gamma
\!\!
\int
\!
\frac{K(\mathbf{x},\mathbf{y})\,\rho_{\mathrm{bound}}(\mathbf{y})}{c^2}
\,dV_y\, ds_x
\]

with \(K\sim 1/|\mathbf{x}-\mathbf{y}|\) (or the medium’s actual free Green function). Then a compact lock automatically produces

\[
\ell(\mathbf{x})-\ell_0
\;=\;
\frac{\alpha E_\star}{c^2 |\mathbf{x}|}
\;+\;\cdots
\;=\;
\frac{\alpha M}{|\mathbf{x}|}
\;+\;\cdots
\]

and \(G_{\mathrm{eff}}\) is fixed by \(\alpha\) and medium normalization:

\[
\frac{2 G_{\mathrm{eff}}}{c^2}
\;\longleftrightarrow\;
\alpha
\quad\text{(match §3.1)}.
\]

**Why this class is reverse-favored:**

1. Bound energy can stay compact with finite \(E_\star\).  
2. Long-range path cost is forced by the free medium’s **linear response kernel**, not by non-integrable free-energy density.  
3. Local \(c\) can be preserved while global chart costs gain \(1/r\) tails.  
4. Matches NC-3.3 escape hatch A-escape-1.

**B/D target:** implement smallest free-response kernel that yields \(\ell-\ell_0 \propto M/r\) when locks form and free budget drops inside the lock.

---

## 7. Free-budget profile — what *is* forced vs optional

| Quantity | Forced (weak field)? | Profile |
|----------|----------------------|---------|
| Path-cost contrast \(\ell(r)-\ell_0\) | **Yes** (for light targets) | \(\propto M/(c^2 r)\) monopole |
| Clock/rod fractional shifts | **Yes** (if full Einstein light factor desired) | same class \(\propto M/(c^2 r)\) |
| Free **energy** density \(\rho_{\mathrm{free}}(r)\) | **Not uniquely** | compact deficit preferred under finite budget; long-range \(1/r\) energy density **disfavored** |
| Bound density \(\rho_{\mathrm{bound}}\) | **Yes** (integral) | \(\int\rho_{\mathrm{bound}}=E_\star\); support ~ lock size |
| Relation \(n(\rho_{\mathrm{free}})\) local | **Optional** | conflicts with long-range if budget local |

### 7.1 Minimal free-budget statement that survives

\[
\boxed{
\begin{aligned}
&\int \mathcal{E}_{\mathrm{bound}}\,dV_{\mathrm{op}} = E_\star = M c^2, \\
&\rho_{\mathrm{free}}\ \text{reduced inside lock (free channels locked)}, \\
&\ell[\mathcal{M}]-\ell_0\ \text{exterior monopole}\ =\ \frac{2G_{\mathrm{eff}}M}{c^2 r}+\cdots
\end{aligned}
}
\]

---

## 8. Matching \(G_{\mathrm{eff}}\) (reverse)

From Class C kernel amplitude \(\alpha\):

\[
G_{\mathrm{eff}}
\;=\;
\frac{\alpha\, c^2}{2}
\quad\text{(schematic; normalize to exact eikonal convention)}.
\]

Equivalence-principle-like demand: \(\alpha\) depends on free-medium constants only, **not** on lock species — otherwise \(E_\star\) and lensing mass disagree (kills monist single-ledger story).

---

## 9. What to hand other approaches

### FOR_A

- Axioms must imply exterior path-cost monopole \(\propto E_\star/(c^2 r)\) from locks + local \(c\), **or** state a different weak-field light law and accept killing Einstein-factor targets.  
- Prefer nonlocal free-response over local \(n(\rho)\) as fundamental.

### FOR_B

- Measure \(\Delta\mathcal{L}\) and \(\Delta\theta\) vs \(E_\star\), not only \(\int T_{00}\).  
- Kill candidate if lock forms with free deficit but rays stay straight in the free-path metric.  
- Try: compact lock + free Green kernel → \(1/r\) cost; do **not** soft-code Poisson for rays.

### FOR_D

- Loss: match \(\delta n_{\mathrm{eff}}(r)\) or \(\ell(r)\) to \(2G M/(c^2 r)\) tail given known \(M=E_\star/c^2\).  
- Adversarial dualist: fixed grid + separate Poisson; monist score must penalize second sector.  
- Penalize recovered media with non-integrable free-energy \(1/r\) if they claim local budget identity.

---

## 10. Falsifiers (medium-side)

1. Stable locks with \(E_\star>0\) but \(\ell-\ell_0=o(1/r)\) (no monopole) → cannot fit deflection \(\propto M/b\).  
2. Monopole amplitude not proportional to \(E_\star\) → inertia/lensing split (dualist mass).  
3. Local \(c\) violated in free regions → fails seed.  
4. Only works when Einstein/Poisson is inserted by hand → ineligible.

---

## 11. Deferred

- Galactic residual free-medium histories (C4).  
- Exact coefficient group (PPN \(\gamma\)) beyond “factor ~2 on space + time.”  
- Strong field / horizons.  
- Vector/tensor free-medium anisotropies (frame dragging analogs).

---

## 12. Bottom line

**Forced (primary):** weak-field free **path-cost** profile

\[
\ell(r)-\ell_0 \;\approx\; \frac{2 G_{\mathrm{eff}} M}{c^2 r}.
\]

**Not forced:** free energy density \(\propto 1/r\).  
**In tension:** local \(n(\rho_{\mathrm{free}})\) + pointwise budget + finite \(E_\star\) + long-range lensing.  
**Preferred reverse shape:** compact bound ledger + free-medium response kernel producing \(1/r\) path cost while local free speed remains \(c\).
