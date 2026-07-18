# Free-Response Kernel: Monist vs Dualist Stress Test (v0)

**Approach:** C — reverse theoretical  
**Date:** 2026-07-18 (Round 2)  
**Trigger:** O-001 (kernel postulated, not derived); A-008 FOR_C  
  (“does (B-response) smuggle dualism S2?”); B-003 honesty on \(\Phi=\alpha\int\rho_b/R\).  
**Status:** stress-test verdict + necessary conditions NC-K.*

---

## 0. What is under test

**Object K:** a free-response map from lock / bound budget to exterior free path cost,
schematically
\[
\delta\ell(x)
\;=\;
\int K(x,y)\,\rho_b(y)\,dV_y
\quad\text{or}\quad
\delta\ell = \mathcal{K}[\rho_b,\rho_f].
\]

**B Round-1 instance:**
\[
\Phi(x)=\alpha\int\frac{\rho_b(y)}{|x-y|+\varepsilon}\,dA_y,
\qquad
n=1+\Phi.
\]

**Algebraic twin (dualist Plummer / soft Poisson):**
\[
\Phi_{\mathrm{dual}}=-G\int\frac{\rho_m}{R_a},
\quad\text{observables from }\Phi_{\mathrm{dual}}.
\]
D-003: Born rays **isomorphic** when \(G=\alpha\) (or \(\beta\)), same soft core, same \(M\).

---

## 1. Two notions of “monist kernel”

| Sense | Meaning | Ray fit? | PROBLEM.md §7? |
|-------|---------|----------|----------------|
| **Phenomenological monism** | One declared continuum; free/bound labels; no *named* Einstein solve | Yes possible | **Not sufficient** |
| **Ontological monism** | Path cost is *identical* to free-medium dynamics / free connectivity; no independent potential DOF; \(K\) derived or defined as free law of that continuum | Yes if dynamics right | **Required** |

O-001 correctly: Round-1 package is a workable *phenomenology*, not yet a complete monist *theory*.

---

## 2. Dualism smuggling criteria (when is K dualist?)

\(K\) (or \(\Phi=K*\rho_b\)) is **dualist (S2 / soft \(T\to G\))** if any of the following hold:

| ID | Smuggle test | Fail means |
|----|--------------|------------|
| **D-K1** | \(\Phi\) is an **independent field** evolved or solved separately from free medium state | Second sector (metric/potential actor) |
| **D-K2** | Source of \(\Phi\) is **matter density** with no free-budget identity linking source to free capacity | Stage+actor |
| **D-K3** | \(K\) is chosen **only** to match gravity phenomenology, with **no** free-update derivation path | Soft Poisson / softE |
| **D-K4** | Removing the \(\Phi\) pass leaves free medium with **no** exterior path structure, and \(\Phi\) cannot be rewritten as free-path cost of that medium | Bolt-on gravity |
| **D-K5** | Ledger \(M=\int\rho_b\) and lensing amplitude use **independent** couplings (fit \(G\) free of \(E_\star\)) | Split mass (fails C2–C3 joint) |
| **D-K6** | Code/theory has a **“gravity_solver”** or Poisson residual equation as ontology | Explicit S2 |

\(K\) is **monist-eligible** only if **all** of D-K1…D-K6 are avoided **and** NC-K.* below hold.

### 2.1 Algebraic form is not the crime

The Green function \(\sim 1/r\) (3D) or \(\sim\log r\) (2D Laplacian) appears in:

- dualist Poisson gravity, **and**
- monist free-medium linear response (diffusion steady state, screened Coulomb of free charge, constraint projection, graph resistance distance, …).

**Verdict rule:**  
> **Same formula, different ontology.**  
> Dualism = independent potential sourced by matter.  
> Monism = path cost / free capacity response of the **same** continuum that carries free signals and locks.

Isomorphism of rays (D-003) ⇒ **fit cannot decide**. Ontology + derivation + free–bound link must.

---

## 3. Answer to A-008: Does (B-response) smuggle dualism?

### Short verdict

| Claim | Verdict |
|-------|---------|
| Ax4 **(B-response)** as *axiom schema* (“free medium has a linear response producing exterior \(\ell\propto E_\star/r\)”) | **Not dualist by itself** if \(K\) is stipulated as a property of free continuum, not a second field |
| B Round-1 **implementation** \(\Phi=\alpha\int\rho_b/(R+\varepsilon)\) with hand \(\alpha,\varepsilon\) | **Dualist-suspect / soft Poisson** under D-K3 (and borderline D-K1 in code: \(\Phi\) computed as a separate array from \(\rho_b\) only) |
| Same formula **if derived** as free Green of free-medium law | **Monist-eligible** |

### Detailed answer for A

1. **Schema vs instance.**  
   (B-response) as “free continuum responds nonlocally so free-path cost feels compact \(E_\star\)” is monist-compatible. It is *not* automatically S2.

2. **Smuggle risk.**  
   The moment one writes \(\nabla^2\Phi = 4\pi G\rho_b\) (or posts \(\Phi=\alpha\int\rho_b/R\) **as gravity**) with \(\rho_b\) as matter and free medium not generating \(\Phi\), one has S2. B’s honesty note (“postulated, not derived”) is exactly D-K3.

3. **What would clear A’s claim.**  
   A must eventually **derive** \(K\) from free-update axioms (Ax3+Ax7+dynamics), or define path cost **purely** as free influence structure (relational), so that \(\ell-\ell_0\propto E_\star/r\) is a **theorem**, not a second postulate isomorphic to Newton.

4. **Current status (Round 2 start):**  
   - (B-response) is a **legitimate monist research target**.  
   - Round-1 kernel numerics are **not yet a monist proof**.  
   - Calling them “one sector” in D’s score is **provisionally OK for Occam pedagogy** only if tagged `postulated_K=true` and not counted as PROBLEM.md §7 pass.

### One-line for A

> **(B-response) does not smuggle dualism if \(K\) is free-medium law; it does if \(K\) is an underived Poisson twin.** A’s claim is correct *conditionally*; Round-1 B does not yet meet the condition.

---

## 4. Necessary conditions for monist \(K\) (NC-K.*)

Any free-response kernel counted as monist for goal condition (2) must satisfy:

### NC-K.1 — Single dynamical sector
State variables are free/bound (or equivalent) of **one** continuum. \(\delta\ell\) is a **functional of that state**, not an independent field with its own kinetic term unless that field *is* free medium.

### NC-K.2 — Free-origin of \(K\)
Either:
- **(Derive):** \(K\) is the Green function / influence kernel of free-medium evolution (diffusion, wave, relaxation, graph update, constraint solve on free DOF only), **or**
- **(Define):** path cost is defined as free-signal travel cost on free connectivity, and the exterior monopole is **proved** from lock boundary conditions on free medium,

not: “insert \(1/R\) because gravity.”

### NC-K.3 — Bound is lock of the same medium
\(\rho_b\) is bound capacity of the continuum (C2 ledger \(E_\star\)), not a foreign dust species.

### NC-K.4 — Budget / depletion link
Forming \(\rho_b\) reduces free capacity in the integral (or dynamical) sense (NC-2.4). Free dynamics that produce \(K\) must **feel** that reduction (locks as boundary / source of free stress), not ignore free deficit.

### NC-K.5 — Local \(c\) preserved
Free signals in local frames still measure speed \(c\) (NC-1.2, NC-3.6). \(K\) affects chart path cost / global null structure, not a second local speed law for free medium.

### NC-K.6 — Universality in \(E_\star\)
Weak-field amplitude \(\propto E_\star\) only (to leading order), independent of lock micro-structure species → \(m_{\mathrm{lens}}=E_\star/c^2\) (joint C2–C3).

### NC-K.7 — No independent \(G\) fit
\(G_{\mathrm{eff}}\) (or \(\alpha\)) is a **medium constant** fixed by free law, not a free parameter fitted per object against orbits while \(E_\star\) floats separately (D-K5).

### NC-K.8 — Kill dualist rewrite
If the theory is rewritten as “matter \(\rho_m\) + Poisson \(\Phi\)” with free medium **idle**, physics must **change** (missing free-bound dynamics, missing unlock, missing free deficit). If rewrite is faithful and free medium never does work, monism failed.

### NC-K.9 — Tagging for numerics
Exports must tag:
```text
sector_count, free_bound_link, K_origin ∈ {derived, defined_relational, postulated},
gravity_solver ∈ {none, soft_poisson, poisson, einstein}
```
Postulated \(K\) ⇒ not monist-complete even if `gravity_solver=none`.

---

## 5. Kill conditions for B’s \(\alpha\int\rho_b/R\)

The Round-1 kernel is **killed as monist proof** (may remain as diagnostic toy) if:

| Kill ID | Condition |
|---------|-----------|
| **KB1** | \(\alpha,\varepsilon\) remain free phenomenological fits with no free-dynamics derivation by the time theory claims completeness |
| **KB2** | Free medium (\(\rho_f\)) can be set uniform with **no** change to \(\Phi\) (Φ depends only on \(\rho_b\)) **and** no relational free-path story — free DOF idle ⇒ dualist matter+potential |
| **KB3** | Inertia \(m_{\mathrm{push}}\) fails to match \(E_\star/c^2\) while rays still see \(\alpha E_\star\) (split ledger) |
| **KB4** | A dualist code path (Poisson on \(\rho_b\)) is byte-equivalent and scored as monist without Occam penalty (softE) |
| **KB5** | Local \(c\) violated in free regions when \(\Phi\) large |
| **KB6** | Only works when \(\rho_b\) is hand-placed forever; spontaneous free→bound formation cannot source \(\Phi\) continuously from free law |

**Partial pass (diagnostic only):**  
KB not triggered for *ray demos* if: compact lock, free deficit measured, rays from \(n=1+\Phi\), tagged `K_origin=postulated`, not claimed as §7 monism.

**Full pass (monist K):**  
Need free evolution such that steady or causal free response **reproduces** \(\sim\alpha\int\rho_b/R\) (or equivalent monopole) with \(\alpha\) predicted; free DOF essential (KB2 fails to hold).

---

## 6. Monist-friendly derivation shapes (targets for A/B)

These are **eligible attempts**, not theorems:

1. **Free diffusion / relaxation**  
   Free capacity \(u\) relaxes with sinks at locks; steady \(\nabla\cdot(D\nabla u)=-\Gamma\rho_b\) ⇒ Coulomb/log Green. Path cost from \(u\) or \(\nabla u\).  
   *Risk:* equation looks Poisson — must identify \(u\) as free medium, not metric potential.

2. **Graph resistance / commute time**  
   Free edges carry conductance; locks remove free conductance. Effective resistance distance ⇒ path cost. Exterior multipoles from network theory.

3. **Constraint projection**  
   Free updates project onto constraints; locks = frozen constraints; influence kernel of projector.

4. **Wave / condensate eikonal**  
   Free field with local \(c\); bound = condensate; effective index from order parameter **nonlocal** (integral kernel in hydro closure).

5. **Relational only**  
   No chart Green: define influence count through free medium; prove asymptotic cost \(\propto E_\star/r\) in embedding chart.

---

## 7. Scoring implication for D (Occam is not enough)

| Model | L_fit vs monist kernel | Dualism? |
|-------|------------------------|----------|
| Postulated monist \(K\sim 1/R\) | Can be 0 | Incomplete monism (D-K3) |
| Dualist Plummer | 0 (isomorphic) | Dualist (D-K1–2) |
| Derived free Green \(\sim 1/R\) | ~0 if law right | Monist-eligible |
| SoftE (Poisson labeled monist) | 0 | Hard fail (NC strip) |

**New reverse necessary condition (from D-003 + this stress):**

> **NC-obs:** Observational degeneracy between free-response kernel and dualist \(\Phi\) ⇒ monism **requires** ledger/link/derivation criteria beyond light bending alone.

(Already in D-003 FOR_C; adopted as NC-obs.)

**Recommended D tags:**  
`K_origin`, `N_sec`, `link`, plus **derivation score** \(L_{\mathrm{deriv}}\):  
\(0\) if derived/defined relational, \(>0\) if postulated (even if N_sec=1).

---

## 8. Congruence with O-001

| O-001 claim | C Round-2 stance |
|-------------|------------------|
| Local GRIN dead | **Confirmed** — formal Theorems 1–2 |
| Preferred shape compact lock + nonlocal response | **Confirmed as target** |
| Kernel postulated ≠ monist theory | **Confirmed** — D-K3 / KB1 |
| Goal (2) not met | **Confirmed** until NC-K.2 holds in theory+numerics |

---

## 9. Bottom line

**Verdict on Round-1 kernel:**  
**Phenomenologically useful, ontologically dualist-suspect.** Same algebra as soft gravity; monist only if/when free-medium dynamics own \(K\).

**Verdict on (B-response) axiom schema:**  
**Conditionally monist** — allowed as research target; not a dualist smuggle *if* NC-K.1–K.9 are acceptance tests.

**Kill line for B:**  
Hand \(\alpha\int\rho_b/R\) without free-origin derivation **cannot** satisfy goal condition (2), even with beautiful rays and Occam N_sec=1.

**Pass line for B:**  
Evolve free medium; measure emergent path-cost monopole \(\propto E_\star\); \(\alpha\) predicted or fixed by free law; free deficit essential; tag `K_origin=derived`.
