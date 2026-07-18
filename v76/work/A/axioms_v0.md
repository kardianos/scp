# A1 — Locality-First Axiom Set (v0)

**Approach:** A (forward theoretical / relational)  
**Date:** 2026-07-18  
**Status:** draft Round 1 (v0.1) — killable claims; not a finished theory  
**Depends on:** `v76/PROBLEM.md` §7, `v76/APPROACHES.md` §1–2, A1 sketch  
**Cross-check:** C reverse notes (`work/C/path_cost_profile_v0.md`, necessary conditions)  
**Eligibility stance:** monist — no fixed background metric as independent distance carrier; no \(T\to G\) bolt-on as ontology.

---

## 0. Design intent

Postulate the **least** structure that could still yield:

1. free vs bound continuum states;
2. \(c\) as free-field locality (not an external light postulate);
3. a budget identity (free + bound);
4. lock stability;
5. rest mass \(m = E_{\mathrm{bound}}/c^2\) as a **derived** ledger relation;
6. warp as what **constant local \(c\)** looks like around locks in a global chart.

We do **not** start with \(\mathbb{R}^{1,3}\) and paint fields on it. Coordinates and metrics are **operational summaries** of free-field signalling after locks exist.

**Slogan (seed hinge):**  
> Mass is field locked; \(c\) is free-field locality, constant in our frame; energy is the continuum ledger; warp is what constant local \(c\) looks like around locks.

---

## 1. Primitive vocabulary (pre-axiomatic)

| Term | Working meaning |
|------|-----------------|
| **Event** | An elementary occurrence of the continuum (point in the event set). |
| **State** | Local configuration at an event: free / bound / mixed (see Ax2). |
| **Influence** | Possibility that the state at \(e_1\) can depend on the state at \(e_0\). |
| **Free update** | Reconfiguration of unlocked medium that can propagate influence. |
| **Lock** | A self-sustaining bound region (particle-form / mass-form). |
| **Observer frame** | A local choice of rods/clocks built from free medium (Ax7). |
| **Budget density** | Scalar ledger of continuum reconfiguration capacity (not a second fluid). |

These words are fixed by the axioms below; they are not dualist “matter” vs “space.”

---

## 2. Axioms

### Ax1 — One continuum (event set)

There exists a set \(\mathcal{E}\) of **events**. Every physical content of the theory is a structure **on** \(\mathcal{E}\) (states, influence, budget). There is no second carrier of distance independent of that structure.

**Kill condition:** if any later construction treats a fixed metric \(g_{\mu\nu}\) as primitive input (not defined from free-field locality + locks), Ax1 is violated.

### Ax2 — Local states: free / bound / mixed

At each event \(e \in \mathcal{E}\) the continuum has a **state** \(s(e)\) in a state space \(\mathcal{S}\) admitting a decomposition into **free** and **bound** ledger components (densities or measures):

\[
\rho_{\mathrm{free}}(e),\quad \rho_{\mathrm{bound}}(e),\quad
\rho_{\mathrm{free}}(e)\ge 0,\quad \rho_{\mathrm{bound}}(e)\ge 0.
\]

- **Free:** medium available to propagate free updates (Ax3).
- **Bound:** medium locked in mass-form (Ax5); not free to leave without unlocking.
- **Mixed:** both components positive at the same event (allowed).

“Empty space” is not a second substance: it is free medium in a nearly uniform unlocked state (\(\rho_{\mathrm{bound}}\approx 0\), \(\rho_{\mathrm{free}}\) nearly constant).

**Optional refinement (not required in v0):** radiation packets = free field *in flight* (still free component).

### Ax3 — Free-field locality (definition of \(c\))

There is a **local causal bound** on free updates. For any local observer frame \(F\) (Ax7):

\[
c_F \;\equiv\;
\sup\{\, v : \text{a free disturbance can influence at coordinate speed }v\text{ in }F\,\}.
\]

**Axiom content:** \(c_F\) is finite and positive wherever free medium exists and free updates are possible. We write \(c\) for this bound when the frame is fixed by convention (Ax7).

**Reading:** \(c\) is **not** “speed of light as an external constant of empty spacetime.” Light (or the fastest free signal) **is** free-field locality. Same fact, two names.

**Kill condition:** defining \(c\) by embedding into a priori Minkowski space and reading off null cones of that metric (circular dualism).

### Ax4 — Budget identity (shared medium ledger)

There is a **total budget density** \(\rho_{\mathrm{tot}}\) such that, at every event (or as an integral identity over every region \(R\subset\mathcal{E}\)):

\[
\boxed{
\rho_{\mathrm{free}}(e) + \rho_{\mathrm{bound}}(e)
\;=\;
\rho_{\mathrm{tot}}(e)
}
\tag{B}
\]

v0 **strong form (local):** \(\rho_{\mathrm{tot}}\) is a conserved scalar density of the *same* medium — not a fixed *metric*. Simplest closed model:

\[
\rho_{\mathrm{tot}}(e) = \rho_0 = \mathrm{const}
\quad\Rightarrow\quad
\rho_{\mathrm{free}} = \rho_0 - \rho_{\mathrm{bound}}.
\]

**Weaker / preferred forms (v0.1 after C reverse):**

| Form | Statement | Status |
|------|-----------|--------|
| **(B-local)** | pointwise \(\rho_f+\rho_b=\rho_0\) | pedagogically clear; **insufficient alone** for long-range path-cost monopole if path cost is local in \(\rho_f\) (C no-go) |
| **(B-integral)** | \(\int(\rho_f+\rho_b)\) conserved; compact lock removes free capacity \(\sim E_\star\) | **required** (C NC-3) |
| **(B-response)** | free medium has a **linear response / Green kernel** so exterior free-path cost feels compact \(E_\star\) as \(\ell-\ell_0\propto E_\star/(c^2 r)\) | **preferred forward target** (C Class C); not a second gravity sector — free continuum response |

**Monist content of (B):** forming mass-form **uses** free budget (at least integrally). Free fabric capacity participating in free updates is reduced by locks. Geometry of free signals is free-medium arrangement, not \(T\to G\).

**C no-go (adopted as design constraint):**  
local refractive law \(n=n(\rho_{\mathrm{free}})\) **plus** compact support of \(\rho_{\mathrm{bound}}\) under **(B-local)** **cannot** produce exterior \(\delta n\propto 1/r\) with finite \(E_\star\) (divergent free-energy if \(\delta\rho_f\propto 1/r\)). Escape: nonlocal path cost / free response kernel; extended free rearrangement; weaken to (B-integral); or index from strain/gradients not only \(\rho_f\).

**Kill condition:** mass defined as \(\int T_{00}[\phi]\,dV\) on fixed flat measure with no free-budget term; **or** claiming monist long-range lensing from (B-local)+local \(n(\rho_f)\) alone without an escape hatch.

### Ax5 — Lock operation and stability

A **lock** is a compact (or effectively localized) region \(L\subset\mathcal{E}\) together with a bound configuration such that:

1. \(\int_L \rho_{\mathrm{bound}}\,d\mu > 0\) (positive bound budget; measure \(d\mu\) from free rods — Ax7);
2. a **stability predicate** \(\mathrm{Stable}(L)\) holds: free updates in a neighborhood cannot unlock \(L\) on timescales short compared to free transit across \(L\);
3. unlocking requires a finite **unlock energy** \(E_{\mathrm{unlock}}(L)\) drawn from free budget (or conversion of bound → free).

**Formation:** free medium may enter bound form when a stability condition is met (self-sustaining defect / phase lock / topological trap — mechanism deferred; only existence of some lock rule is axiomatic here).

**Stability sketch (relational, not PDE):** \(\mathrm{Stable}(L)\) iff every free perturbation of amplitude below a threshold \(\varepsilon_L\) leaves the bound integral and the lock’s causal exterior unchanged up to gauge of free oscillations. This is the monist stand-in for “particle stability” (medium lock, not only VK stability of \(\phi\) on \(\mathbb{R}^3\)).

### Ax6 — Energy as ledger of reconfiguration capacity

**Energy** is a frame-dependent integral of the continuum ledger — not a fluid poured into space.

In a rest frame \(F\) of a lock \(L\) (Ax7), define:

\[
E_{\mathrm{bound}}(L;F)
\;=\;
\int_L \mathcal{E}\bigl[\rho_{\mathrm{bound}}\bigr]\,d\mu_F,
\]

where \(\mathcal{E}\) is a positive functional of bound density (v0 simplest: \(\mathcal{E}[\rho]=\rho\), so energy density = budget density in chosen units).

Total energy in a region is free + bound ledgers (plus free packets in flight). Conservation, when it holds, is conservation of **medium budget under the dynamics**, not conservation of a substance on a fixed grid.

### Ax7 — Observer frames enforce local \(c = \mathrm{const}\)

An **observer frame** \(F\) is a local system of rods and clocks **made of free medium**:

- clock ticks = periods of a free local oscillator (or free-signal round trips);
- rod lengths = free-signal path costs calibrated so that free signals have isotropic speed \(c_F\) **in that frame**.

**Axiom:** For every event where free medium exists, there exist local frames in which:

\[
c_F = c = \mathrm{const}
\quad\text{(isotropic, independent of free-signal source motion within the usual relativity statements)}.
\]

This is **operational**: we *define* rods/clocks so free-field locality looks constant and isotropic **here**. It is not “flat space pre-exists and light is weird.”

### Ax8 — Influence structure from free medium (emergent geometry)

Between events, **influence** for free signals is constrained by free medium arrangement and locks:

- free signals cannot outrun local \(c\) in any free frame (Ax3+Ax7);
- free-signal **path cost** \(\ell\) is a functional of free-medium state and locks — **not** required to be a strictly local function of \(\rho_{\mathrm{free}}(e)\) alone (v0.1; C reverse).

**Preferred structure (forward target):** a free linear response around a compact lock with ledger \(E_\star\),

\[
\ell(r)-\ell_0 \;\sim\; \frac{\alpha E_\star}{c^4\, r}
\quad\text{(weak field exterior)},
\]

with \(G_{\mathrm{eff}}\) identified from kernel amplitude \(\alpha\) (emergent, not a second sector). Local frames still measure free speed \(c\) (Ax7).

A **global chart** that pretends free medium is still uniform will assign nontrivial path geometry to free signals. That assignment is **warp** (Ax9). The metric \(g_{\mu\nu}\) is **defined** as whatever makes free null structure \(g_{\mu\nu}k^\mu k^\nu=0\) with local speed \(c\) true — not an independent Einstein-sourced field.

### Ax9 — Warp theorem schema (to be proved, not bolted)

**Claim schema (target theorem, not extra axiom):**

\[
\underbrace{c\ \mathrm{constant\ locally\ (Ax7)}}_{\text{free-field law}}
\;+\;
\underbrace{\text{locks with }E_{\mathrm{bound}}>0\ (Ax5,Ax6)}}_{\text{rearrangement}}
\;+\;
\underbrace{\text{budget identity (Ax4)}}_{\text{free less near mass}}
\;\Longrightarrow\;
\underbrace{\text{nontrivial free-signal geometry around locks}}_{\text{warp}}
\]

in any global chart that covers free and locked regions.

**Corollaries intended:**

1. Mass without warp is inconsistent under monism (locks rearrange free paths; constant local \(c\) forces curved global description).
2. Lensing / Shapiro-like delay are free-path cost effects, not \(T\to G\) as ontology.
3. Inertia of a lock is resistance to reconfiguring bound budget relative to surrounding free medium (relational; see §4).

---

## 3. Derived: mass from locality (sketch → see `mass_from_locality.md`)

### Definition (rest mass)

In the rest frame \(F_0\) of a stable lock \(L\):

\[
\boxed{
m(L) \;\equiv\; \frac{E_{\mathrm{bound}}(L;F_0)}{c^2}
}
\tag{M}
\]

where \(c\) is free-field locality in \(F_0\) (Ax3+Ax7).

### Why this is not dualist \(E=mc^2\)

| Dualist habit | Here |
|---------------|------|
| \(m\) = amount of matter-stuff | \(m\) = how much continuum is locked |
| \(E\) = energy *of* that stuff | \(E_{\mathrm{bound}}\) = ledger of bound medium |
| \(c\) from Minkowski light cones | \(c\) = free update bound of the same continuum |
| Geometry separate; \(T\) sources \(G\) | Geometry = free-signal bookkeeping under (B)+(locks) |

### Why \(c^2\) (dimensional / operational sketch)

1. Free-signal locality sets the conversion between **time** and **length** for free rods/clocks: \([c]=\mathrm{L}/\mathrm{T}\).
2. Inertial resistance of a lock under free-medium “push” involves accelerating a bound structure whose causal response is limited by \(c\) (signals that reconfigure the lock’s free envelope travel at most at \(c\)).
3. Unlock energy \(E_{\mathrm{unlock}}\) and the energy cost of boosting the lock both scale with bound ledger; matching them forces \(E_{\mathrm{bound}} = m c^2\) once \(m\) is defined as the inertial coefficient of the lock against free-medium drag (full derivation in `mass_from_locality.md`).

v0 status: **(M) is a definition + consistency requirement**; full inertia derivation is open but the dualist alternative (mass as independent substance) is excluded by Ax1–Ax4.

---

## 4. Relational package (optional strengthening)

Not required for eligibility, but natural with Ax1–Ax9:

- **Mach-like inertia:** \(m(L)\) is only meaningful relative to free medium elsewhere (resistance of lock vs free web).
- **Path cost:** “distance” between free events = infimum free-signal time × \(c\), with locks as obstacles / index gradients via \(\rho_{\mathrm{free}}\).
- **No empty container:** relation counts and free budgets *are* geometry.

---

## 5. Explicit non-axioms (rejected smuggling)

| Rejected | Why |
|----------|-----|
| Fixed \(g_{\mu\nu}\) as stage | Fails PROBLEM §7.1; dualist |
| \(G_{\mu\nu}=8\pi G T_{\mu\nu}/c^4\) as ontology | Secondary sourcing |
| Q-ball on fixed lattice as monist proof | Core is high \(\|\Phi\|\), not free-budget depletion as geometry |
| Defining \(c\) from a priori Minkowski | Circular (Ax3 kill condition) |
| Locks that do not reduce free budget | No monist gravity channel |

---

## 6. Open choices (flagged for later rounds)

1. **Budget form:** (B-local) pedagogical only if path cost nonlocal; (B-integral)+(B-response) preferred after C. Kernel \(K\sim 1/r\) vs strain-based index still open.
2. **\(\mathcal{E}[\rho_{\mathrm{bound}}]\):** identity vs nonlinear functional (saturation).
3. **Signature:** does local \(c=\mathrm{const}\) force Lorentzian signature, or only approximately?
4. **Rod/clock circularity:** clocks made of free medium that itself feels locks — Approach C1 critical path.
5. **Dynamics:** Ax5 states stability existence; free response kernel + lock rule are Approach B’s job.
6. **Einstein recovery:** exterior \(\ell-\ell_0\propto M/(c^2 r)\) is a *theorem target* for free response (C forced profile), not an axiom and not \(T\to G\).
7. **Inertia expansion:** derive \(\tfrac12(E_\star/c^2)v^2\) from free locality + path cost (C FOR_A; still open — `mass_from_locality.md`).

---

## 7. Success criteria for this axiom set (A — theoretical)

| Criterion | v0 status |
|-----------|-----------|
| \(c\) from field locality alone | **Stated** (Ax3); not circular if Ax7 operational |
| \(m=E_{\mathrm{bound}}/c^2\) derived / defined without dualist matter | **Defined (M)**; inertia match **sketched** |
| Warp from local \(c\) + locks | **Schema (Ax9)**; proof needs dynamics + free-path cost |
| No fixed background metric as actor | **Stated** (Ax1, Ax8) |
| Depletion identity | **Stated** (Ax4) |

---

## 8. Cross-approach handoffs

| Tag | Content |
|-----|---------|
| **FOR_B** | Implement budget (B) + lock flag + free updates only within locality radius; measure free drop where bound rises; rays = fastest free paths under local \(c\). Do not fix Euclidean distance as truth. |
| **FOR_C** | Reverse-check: is strong (B) forced by weak lensing, or only integral free deficit \(\propto M\)? Rod/clock from free medium — necessary conditions. |
| **FOR_D** | Score media by: free drop co-located with bound; ray deflection without Poisson solve; penalize second gravity sector. Use (M) as consistency \(m_{\mathrm{ray}}\leftrightarrow E_{\mathrm{bound}}/c^2\). |

---

## 9. Change log

| Ver | Note |
|-----|------|
| v0 | Round 1 draft: Ax1–Ax9, (B), (M), warp schema, kill conditions |
| v0.1 | Adopt C reverse: (B-local) no-go for local optics+long-range 1/r; prefer free response kernel; Ax8 path-cost monopole target; strip-list alignment |
| v0.2 | Round 2: free-response \(K\) defined in `free_response_kernel_v0.md` (F1 preferred dynamical principle; F5 hand kernel = dualist residue / target only). Congruence package for B/C/D. |
| v0.3 | Round 3: B-M2 free Laplace endorsed monist; 3D F1 in `free_response_3d_v0.md`; \(G_{\mathrm{eff}}=\gamma s c^4/(8\pi\sigma_0)\); congruence T2 dimension-aware; hand 1/r still dead. |
