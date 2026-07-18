# Approach B — Minimal Medium Design v0

**Date**: 2026-07-18  
**Track**: Forward numerical (B2 primary; B1 graph deferred)  
**Status**: sandbox-eligible design (not a full theory)

---

## 0. Goal of this medium

Exhibit, in one evolution / state, that:

1. **Free budget drops where bound (mass-form) rises** — shared continuum budget.
2. **Ray deflection / delay appears without a second gravity solver** — geometry of free medium is the warp.
3. **Local free-signal speed \(c\)** is fixed by the free-field update law (locality seed).

This is a **killable sandbox**, not a claim that GR is derived.

---

## 1. Why B2 first (not B1)

| Option | Pros | Cons for Round 1 |
|--------|------|------------------|
| **B1** causal graph | Distance = path cost naturally | Graph geodesic machinery heavier; formation rules underspecified |
| **B2** free-density continuum | Budget identity trivial; ray eikonal cheap | Chart is fixed 2D plane (dualist residue — flagged below) |
| **B3** metric-from-locks | Closest to Einstein diagnostic | Needs B2 refractive structure first |
| **B4** frame enforcement | Full monist hinge | Needs rod/clock model from A/C |

**Choice:** **B2-lite optical monism** on a fixed lab chart, with explicit honesty that the *chart* is scaffolding; the *degrees of freedom* for mass and warp are the same \((\rho_{\mathrm{free}},\rho_{\mathrm{bound}})\).

---

## 2. Continuum content (one medium)

State on a 2D lab chart \((x,y)\):

\[
\rho_{\mathrm{bound}}(x,y) \ge 0,\qquad
\rho_{\mathrm{free}}(x,y) \ge 0,
\]

with **strict budget identity** (v0 choice — killable if too strong):

\[
\rho_{\mathrm{free}} + \rho_{\mathrm{bound}} = \rho_0 = \mathrm{const}.
\]

| Symbol | Meaning |
|--------|---------|
| \(\rho_0\) | Uniform total medium budget density |
| \(\rho_{\mathrm{bound}}\) | Medium locked into mass-form (particle-form) |
| \(\rho_{\mathrm{free}}\) | Medium free to carry free signals |
| \(c\) | Local free-signal locality bound (set = 1 in code units) |

**Mass (ledger, rest frame of lock):**

\[
m \;=\; \frac{1}{c^2}\int \rho_{\mathrm{bound}}\,dA
\quad\text{(v0: \(\mathcal{E}[\rho_b]=\rho_b\); units where \(c=1\))}.
\]

No separate matter field. No \(T_{\mu\nu}\to G_{\mu\nu}\) pass.

---

## 3. Locality-\(c\) and optical structure

### 3.1 Local free speed

Wherever \(\rho_{\mathrm{free}} > 0\), free signals have **local** speed \(c\) when measured with free-medium rods/clocks. In the **lab chart** used by the sandbox, free signals obey an eikonal with **chart speed**

\[
v_{\mathrm{chart}}(x) \;=\; \frac{c}{n(x)},
\]

where the refractive index is **defined from free density only**:

\[
n(x) \;=\; \frac{\rho_0}{\rho_{\mathrm{free}}(x)}
\;=\; \frac{\rho_0}{\rho_0 - \rho_{\mathrm{bound}}(x)}
\quad\bigl(\rho_{\mathrm{free}} \ge \varepsilon_{\min}\bigr).
\]

Interpretation (v0 slogan, not theorem):

- Depleted free medium \(\Rightarrow\) higher chart index \(\Rightarrow\) longer optical path.
- Local \(c\) remains the free-medium locality; the lab chart sees Shapiro-like delay and deflection **as identity of the medium arrangement**, not as a second field.

Floor \(\varepsilon_{\min}\) prevents \(n\to\infty\) (horizon analog deferred).

### 3.2 What is *not* done

- No Poisson solver \(\nabla^2\Phi = 4\pi G\rho\).
- No Einstein integrator.
- No Q-ball / scp_sim kernel.
- Rays do **not** feel \(\rho_{\mathrm{bound}}\) except through the budget \(\rho_{\mathrm{free}}=\rho_0-\rho_{\mathrm{bound}}\).

---

## 4. Lock placement and optional formation

### 4.1 Hand-placed lock (Round 1 primary)

Gaussian bound clump:

\[
\rho_{\mathrm{bound}}(r) \;=\; A\,\exp\bigl(-r^2/(2\sigma^2)\bigr),\quad
A < \rho_0 - \varepsilon_{\min},
\]

then \(\rho_{\mathrm{free}}=\rho_0-\rho_{\mathrm{bound}}\).

### 4.2 Formation rule (optional, B2 exchange)

Discrete relaxation (for later rounds):

\[
\partial_t \rho_{\mathrm{bound}} \;=\; \gamma\,\Theta(\rho_{\mathrm{free}}-\rho_*)\,\rho_{\mathrm{free}}
\;-\; \mu\,\rho_{\mathrm{bound}},
\]

with \(\partial_t\rho_{\mathrm{free}}=-\partial_t\rho_{\mathrm{bound}}\) (budget). Round 1 may only **hand-place** to isolate deficit + ray metrics.

---

## 5. Ray / path metrics (no gravity sector)

### 5.1 Eikonal rays (Hamiltonian form)

In 2D, ray Hamiltonian with \(H = |p| - n(\mathbf{x}) = 0\) (or \(|p|/n = 1\) in \(c=1\)):

\[
\dot{\mathbf{x}} = \frac{\mathbf{p}}{|\mathbf{p}|}\,v_{\mathrm{chart}}
\;=\; \frac{\mathbf{p}}{|\mathbf{p}|}\,\frac{c}{n},\qquad
\dot{\mathbf{p}} = c\,\nabla\log n
\quad\text{(equivalent form used in code: optical length \(L=\int n\,ds\))}.
\]

Implementation v0: integrate rays with optical metric \(ds_{\mathrm{op}} = n\,ds_{\mathrm{eucl}}\) via gradient of \(n\) (Snell / eikonal ODE).

### 5.2 Observables

| Symbol | Definition |
|--------|------------|
| \(\Delta\rho_{\mathrm{free}}\) | \(\rho_0 - \langle\rho_{\mathrm{free}}\rangle_{\mathrm{core}}\) vs exterior (deficit) |
| \(\delta_{\mathrm{defl}}(b)\) | asymptotic deflection angle vs impact parameter \(b\) |
| \(\Delta t_{\mathrm{Sh}}(b)\) | excess travel time vs vacuum ray at same endpoints |
| \(m\) | \(\int\rho_{\mathrm{bound}}\,dA / c^2\) |
| \(\delta_{\mathrm{defl}}\) vs \(m\) | scaling check (weak-field-like \(\propto m/b\) *aspirational*, not enforced) |

Success for Round 1:

1. \(\rho_{\mathrm{free}}\) lower where \(\rho_{\mathrm{bound}}\) higher (identity).
2. \(\delta_{\mathrm{defl}} \ne 0\) and/or \(\Delta t_{\mathrm{Sh}} > 0\) for nonzero lock — **without** Poisson/Einstein.
3. Vacuum control (\(A=0\)): zero deflection, zero excess delay (to machine noise).

---

## 6. Dualist residue (honest)

| Residue | Mitigation plan |
|---------|-----------------|
| Fixed Euclidean chart for storage | Later: path-cost only (B1) or reconstruct \(g\) from local-\(c\) (B3/B4) |
| \(n=\rho_0/\rho_{\mathrm{free}}\) is a *choice* | C3 target profile can kill or force this map |
| Hand-placed lock | Formation dynamics Round 2 |
| No inertia test of locks | Measure \(m\) vs acceleration response later |

Monist-eligible **core**: free and bound share one budget; rays see only free structure; no second gravity pass.

---

## 7. Code map

| Path | Role |
|------|------|
| `work/B/medium_design_v0.md` | This design |
| `work/B/sandbox_b2_lite.py` | Single-file Python sandbox |
| `work/B/outputs/` | Figures + `results.json` |

---

## 8. Round 1 acceptance

PASS if sandbox runs and reports:

- budget residual max \(\lvert\rho_f+\rho_b-\rho_0\rvert < 10^{-12}\),
- core free deficit \(> 0\) for \(A>0\),
- max \(|\delta_{\mathrm{defl}}| > 0\) for \(A>0\),
- vacuum control \(\max|\delta_{\mathrm{defl}}| \approx 0\).

FAIL / kill if rays only bend after adding an explicit Poisson step, or if free density does not drop under the lock.

---

## 9. Amendment after C reverse no-go (Round 1)

**C no-go** (`path_cost_profile_v0.md`): local \(n(\rho_{\mathrm{free}})\) + pointwise
budget + compact \(E_\star\) **cannot** produce long-range \(\ell-\ell_0\propto 1/r\).

| Channel | File | Long-range path cost? | Free deficit? |
|---------|------|----------------------|---------------|
| **B2-lite local optics** | `sandbox_b2_*.py` | No (compact GRIN only) | Yes (identity) |
| **B2-response kernel** | `sandbox_b2_kernel.py` | Yes: \(\Phi=\alpha\int\rho_b/(R+\varepsilon)\) | Yes (ledger still local) |

Round 1 **local** sandbox still proves monist bookkeeping + ray without second
solver. It does **not** claim Einstein-class exterior light bending.

Round 2 primary: evolve/derive free-response kernel from medium dynamics
(not postulate \(\Phi\)), keep single sector, match C monopole target
\(\ell-\ell_0 \approx 2G_{\mathrm{eff}}M/(c^2 r)\).

Code map addendum:
| `sandbox_b2_pure.py` | stdlib RK2 + Born |
| `sandbox_b2_kernel.py` | nonlocal free-response path cost |
| `_compute_born_results.py` | Born series writer |
