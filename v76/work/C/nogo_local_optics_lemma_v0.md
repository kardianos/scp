# No-Go Lemma: Local Optics + Local Budget ↛ Long-Range Path Cost (v0)

**Approach:** C — reverse theoretical  
**Date:** 2026-07-18 (Round 2)  
**Status:** formal writeup of Round-1 NC-3.3 / NC-3.4  
**Purpose:** Kill local GRIN monism as a channel for Einstein-class weak-field exterior path cost, without assuming GR as ontology.

---

## 0. Setup and notation

Work in a fixed asymptotic chart \(\mathbb{R}^d\), \(d\in\{2,3\}\), used only as a storage scaffold (dualist residue allowed for *coordinates*; not for *physics ontology*).

| Symbol | Meaning |
|--------|---------|
| \(\rho_b,\rho_f\ge 0\) | bound / free budget densities |
| \(\rho_0>0\) | constant total density (local budget) |
| \(E_\star\) | unlock ledger of an isolated lock; \(E_\star<\infty\) |
| \(n(x)\) | chart refractive index for free signals |
| \(\delta n = n-1\) | excess path-cost density (isotropic optics) |
| \(\mathrm{supp}\) | essential support |

**Local budget (pointwise):**
\[
\rho_f(x)+\rho_b(x)=\rho_0
\quad\text{for a.e. }x.
\tag{LB}
\]

**Compact lock:**
\[
\mathrm{supp}(\rho_b)\subset B_R(0)
\quad\text{for some }R<\infty,
\qquad
\int \rho_b = E_\star < \infty.
\tag{CL}
\]

**Local optics (path cost algebraic in free density):**
\[
n(x)=F\bigl(\rho_f(x)\bigr)
\quad\text{with }F:\,(0,\rho_0]\to[1,\infty)\ \text{Borel measurable},
\tag{LO}
\]
and \(F(\rho_0)=1\) (asymptotic free vacuum index). Equivalently
\(\delta n(x)=f(\delta\rho_f(x))\) with \(\delta\rho_f=\rho_0-\rho_f=\rho_b\) under (LB).

**Weak-field long-range path-cost target (C3 monopole class):** there exist
\(C>0\), \(r_*>R\), and \(\eta\in(0,1]\) such that for all \(r\ge r_*\),
\[
\delta n(r\,\hat e)
\;\ge\;
\frac{C}{r^\eta}
\quad\text{for all unit }\hat e
\quad\text{(isotropic lower bound)},
\tag{LR}
\]
or, weaker but sufficient for “not compact GRIN”:
\[
\liminf_{r\to\infty}\, r^\eta \inf_{|x|=r}\delta n(x)
\;>\; 0.
\tag{LR′}
\]
(The Einstein/Shapiro target uses \(\eta=1\), \(C\sim 2G_{\mathrm{eff}}M/c^2\).)

---

## 1. Lemma A — Compact free contrast under (LB)+(CL)

**Lemma A (support).**  
Assume (LB) and (CL). Then
\[
\rho_f(x)=\rho_0
\quad\text{and}\quad
\delta\rho_f(x)=0
\quad\text{for a.e. }|x|>R.
\]
In particular \(\mathrm{supp}(\delta\rho_f)\subset B_R(0)\).

**Proof.**  
From (LB), \(\delta\rho_f=\rho_b\). From (CL), \(\rho_b=0\) a.e. outside \(B_R\). ∎

---

## 2. Lemma B — Local optics cannot long-range if free contrast is compact

**Lemma B (local map of compact support).**  
Assume (LO) and \(\mathrm{supp}(\delta\rho_f)\subset B_R(0)\). Then
\[
\mathrm{supp}(\delta n)\subset B_R(0).
\]
In particular (LR) and (LR′) fail for every \(\eta>0\), \(C>0\), \(r_*>R\).

**Proof.**  
Outside \(B_R\), \(\rho_f=\rho_0\), so \(n=F(\rho_0)=1\) and \(\delta n=0\).  
For \(r>R\), \(\inf_{|x|=r}\delta n(x)=0\), so the liminf in (LR′) is \(0\). ∎

---

## 3. Theorem 1 — Local GRIN no-go (NC-3.3)

**Theorem 1 (local optics + local budget + compact lock ↛ long-range path cost).**  

Assume (LB), (CL), and (LO). Then the exterior path-cost contrast has compact support:
\[
\delta n(x)=0\quad\text{for a.e. }|x|>R.
\]
Hence no isotropic long-range monopole of class (LR)/(LR′) exists. In particular, free-signal excess delay for impact parameters \(b\gg R\) vanishes at Born order in \(\delta n\) (straight-line integral of zero), and there is no \(\Delta\theta\sim M/b\) tail from this channel alone.

**Proof.**  
Lemma A ⇒ compact free contrast. Lemma B ⇒ compact \(\delta n\). ∎

**Corollary 1.1 (numerical form).**  
Any B2-lite medium with \(n=\rho_0/\rho_f\) (or any \(n=n(\rho_f)\) only) and pointwise \(\rho_f+\rho_b=\rho_0\) with compact \(\rho_b\) is at most a **compact GRIN**. It may bend rays near the lock; it cannot supply Einstein-class exterior \(\ell\sim M/r\).

**Corollary 1.2 (escape necessity).**  
To keep compact locks, finite \(E_\star\), and long-range path cost, a medium must violate at least one of (LB), (CL), (LO). Preferred escape (C Class C): replace (LO) by a **nonlocal** free-response path-cost functional
\[
\delta n(x)=\bigl(K*\rho_b\bigr)(x)
\quad\text{or more generally}\quad
\delta\ell = \mathcal{F}[\rho_f,\rho_b]
\text{ nonlocal in space}.
\]

---

## 4. Theorem 2 — Non-integrable free-energy tail (NC-3.4)

Now drop (CL) and ask whether long-range **free-energy** density can carry the optics.

**Assume** a linear (or bi-Lipschitz near vacuum) local optics map:
\[
\delta n(x)
\;=\;
\beta\,\frac{\delta\rho_f(x)}{\rho_0}
\;+\;
o\!\bigl(\delta\rho_f/\rho_0\bigr),
\qquad
\beta\neq 0,
\tag{LO-lin}
\]
and an exterior isotropic tail
\[
\delta n(r)
\;\sim\;
\frac{C}{r^\eta},
\qquad
C\neq 0,\ \eta\le d,
\quad r\to\infty.
\tag{Tail}
\]

**Theorem 2 (divergent free-energy deficit).**  
Under (LO-lin) and (Tail) with \(\eta\le d\),
\[
\int_{|x|>r_*}\delta\rho_f(x)\,d^dx
\;=\;\infty
\]
for every \(r_*<\infty\). Hence no finite unlock ledger can equal the integrated free-energy deficit if that deficit is identified with \(\int\delta\rho_f\).

**Proof.**  
From (LO-lin), \(\delta\rho_f\sim (\rho_0/\beta)\,C\, r^{-\eta}\) for large \(r\). In polar coordinates,
\[
\int_{r_*}^\infty r^{d-1}\, r^{-\eta}\,dr
\]
diverges at infinity whenever \(\eta\le d\). For physical targets \(\eta=1\), \(d\in\{2,3\}\), divergence is immediate. ∎

**Corollary 2.1.**  
One cannot simultaneously hold:

1. local linear optics \(\delta n\propto\delta\rho_f\),  
2. exterior \(\delta n\sim C/r\) (or any \(\eta\le d\)),  
3. \(\int\delta\rho_f = E_\star <\infty\) as the monist free-energy identity.

**Corollary 2.2 (profile split).**  
Path-cost profile \(\delta\ell\sim M/(c^2 r)\) and free-energy density profile are **not** interchangeable. Finite \(E_\star\) forces free-energy deficit to be integrable (if it is the budget); long-range path cost must then be **non-energy** structure (response kernel, strain, connectivity) or non-integrable free energy must be abandoned as the budget identity.

---

## 5. What is *not* claimed

| Not claimed | Why |
|-------------|-----|
| Full GR is false | Targets are phenomenological; lemma is about media classes |
| Nonlocal kernels are monist | Separate stress test (`kernel_dualism_stress_v0.md`) |
| No medium can lens | Only this *class* fails |
| (LB) is always wrong | (LB) may hold if path cost is nonlocal |
| 2D vs 3D changes the no-go | Only changes \(\eta\le d\) threshold in Theorem 2 |

---

## 6. Empirical alignment (Round 1)

| Source | Observation |
|--------|-------------|
| B-003 | Local channel: free deficit co-located; exterior \(n-1\sim 10^{-5}\) at \(r=6\); compact GRIN only |
| D-002 | `local_optics` loses \(L_{\mathrm{fit}}\) on kernel truth (\(\sim 0.19\)) |
| A-008 | Ax4 demotes (B-local)+local \(n\) for long-range; adopts no-go |

---

## 7. Formal status

- **Logical strength:** Theorems 1–2 are elementary measure-support + integral tests. Suitable for Lean later (indicator of support; integral divergence of \(r^{d-1-\eta}\)).
- **Physics strength:** Conclusive against local-GRIN-as-gravity; does not settle free-response monism.

---

## 8. Bottom line

\[
\boxed{
\text{(LB)}+\text{(CL)}+\text{(LO)}
\;\Longrightarrow\;
\delta n\text{ compact}
\;\Longrightarrow\;
\text{no long-range path-cost monopole}.
}
\]
\[
\boxed{
\text{(LO-lin)}+\delta n\sim r^{-\eta}\ (\eta\le d)
\;\Longrightarrow\;
\int\delta\rho_f=\infty.
}
\]

**Escape required** for monist weak-field light targets: nonlocal free path-cost (or abandon local budget / compact lock / Einstein-class tail).
