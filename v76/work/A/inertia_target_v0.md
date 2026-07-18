# Inertia Target v0 — Non-Tautological \(m_{\mathrm{inertial}}=E_\star/c^2\)

**Approach:** A  
**Round:** 4  
**Audience:** B (measure), C (EQ check), D (triad score)  
**Status:** protocol specification — not yet a closed theorem  
**Blocks:** goal2_FULL / J5 (O-005)

---

## 0. Why this is hard

Identity \(m = E_\star/c^2\) is **ledger definition** of rest mass.  
**Inertia** is dynamical resistance of a lock to acceleration relative to free medium.

If B sets \(a = F/m_{\mathrm{ledger}}\) by construction, the test is **tautological** (B R2 honesty).  
J5 requires an independent operational \(m_{\mathrm{inertial}}\) from free-medium response.

---

## 1. Target equality (triad)

\[
\boxed{
m_{\mathrm{inertial}}
\;\stackrel{?}{=}\;
\frac{E_\star}{c^2}
\;\stackrel{?}{=}\;
M_{\mathrm{ray}}
}
\tag{T5}
\]

| Symbol | Operational meaning |
|--------|---------------------|
| \(E_\star\) | Unlock / bound ledger \(\int\rho_b\,dV\) (or measured unlock work) |
| \(M_{\mathrm{ray}}\) | Mass inferred from exterior path-cost amplitude via \(G_{\mathrm{eff}}\) / Born slope, with \(G_{\mathrm{eff}}\) **fixed by free constants** (not free-fit \(M\)) |
| \(m_{\mathrm{inertial}}\) | Low-acceleration dynamical mass from free-drag protocol below |

**Pass band:** declare \(\varepsilon_{\mathrm{triad}}\) (e.g. 20–30% on coarse grid; tighten with \(N\)).

---

## 2. Forbidden protocols (tautology / dualism)

| Forbidden | Why |
|-----------|-----|
| Set \(a=F_{\mathrm{ext}}/m_{\mathrm{ledger}}\) and report \(m=m_{\mathrm{ledger}}\) | Tautology |
| Fit \(m\) from dualist \(\Phi\) force \(F=-m\nabla\Phi\) with foreign \(G\) | Dualist sector |
| Prescribe lock trajectory and back out \(m\) from prescribed Newtonian law | Circular |
| Use only ray slope as “inertia” | That’s \(M_{\mathrm{ray}}\), not independent \(m_{\mathrm{inertial}}\) |

---

## 3. Allowed protocol A — Free-\(\psi\) force, free motion of lock

### 3.1 Setup

1. Form compact lock; measure \(E_\star=\int\rho_b\,dV\).  
2. Free state: solve F1-3D for \(\psi[\rho_b]\) (same monist solver as R3).  
3. Define **force density on bound medium from free capacity**:

\[
\mathbf{f}
= -\rho_b\,\nabla\psi
\quad\text{(or }-\rho_b\nabla(\gamma\psi)\text{ if path-cost scaled)},
\]

motivated by free energy \(U=\frac{\sigma_0}{2}\int|\nabla\psi|^2\) and virtual work on \(\rho_b\) placement (same structure as electrostatic force on charge density in a potential — here \(\psi\) is free fabric capacity, not dualist gravity).

4. Total free-medium force on lock:

\[
\mathbf{F}_\psi = \int \mathbf{f}\,dV.
\]

For a centered static lock, \(\mathbf{F}_\psi=\mathbf{0}\) by symmetry (equilibrium).

### 3.2 Acceleration test (two variants)

**Variant A1 — displaced lock (static spring):**  
Shift lock center by small \(\delta\mathbf{x}\) relative to free medium / box; re-solve \(\psi\); measure restoring \(F_\psi(\delta x)\). Effective spring \(k=-dF/d(\delta x)\).  
This tests **confinement / self-force**, not inertial mass directly — useful control, **not** \(m_{\mathrm{inertial}}\).

**Variant A2 — slow boost energy (preferred non-tautological):**  
Give lock a small velocity \(\mathbf{v}\) (\(v\ll c\)) by translating \(\rho_b(x-vt)\) (or advect bound density under free dynamics).  
At each time, re-solve/relax \(\psi\) (or use quasi-static co-moving solve).  
Measure free energy

\[
U(v)=\frac{\sigma_0}{2}\int|\nabla\psi_v|^2\,dV.
\]

Fit

\[
U(v)=U(0)+\frac12 m_{\mathrm{inertial}} v^2 + O(v^4).
\]

Compare \(m_{\mathrm{inertial}}\) to \(E_\star/c^2\).

**Why non-tautological:** \(m_{\mathrm{inertial}}\) comes from free-field energy expansion, not from dividing \(F\) by assumed \(m\).

**Caveat:** self-energy \(U(0)\) is cutoff-sensitive; use \(\Delta U=U(v)-U(0)\) only; renormalization may be needed (kinetic_inertia_v0).

### 3.3 Variant A3 — free-drag under external free drive

1. Apply a **known external free-medium drive** that is **not** defined as \(m_{\mathrm{ledger}}a\): e.g. impose a large-scale free potential ramp \(\psi_{\mathrm{ext}}=\mathbf{g}\cdot\mathbf{x}\) with small \(|g|\) as outer BC / body term in free law  
   \(-\sigma\nabla^2\psi = s\rho_b + \rho_f\,\eta_{\mathrm{drive}}\) carefully documented.  
2. Allow lock center of mass \(\mathbf{X}(t)\) to move under \(\mathbf{F}_\psi=\int(-\rho_b\nabla\psi)\,dV\) **and** free relaxation of \(\psi\).  
3. Integrate dynamics \(\dot{\mathbf{P}}=\mathbf{F}_\psi\) only if \(\mathbf{P}\) is free-field momentum of the lock–wake system — **or** measure early-time \(\ddot X\) before walls matter.  
4. \(m_{\mathrm{inertial}}=F_{\mathrm{net}}/\ddot X\) **only if** \(F_{\mathrm{net}}\) is measured from \(\psi\) integrals and \(\ddot X\) from \(\rho_b\) COM motion under the **same** free law — no insertion of \(m_{\mathrm{ledger}}\).

A3 is easier to code wrong; prefer A2 for Round 4 if short on time.

---

## 4. Protocol B — Unlock energy cross-check

Independently:

1. Measure work to dissolve lock into free medium: \(E_{\mathrm{unlock}}\).  
2. Require \(E_{\mathrm{unlock}}\approx E_\star=\int\rho_b\).  
3. If unlock dynamics absent, report \(E_\star\) only and tag `unlock=ledger_proxy`.

Triad uses \(E_{\mathrm{unlock}}\) when available; else \(E_\star\).

---

## 5. Protocol C — Ray mass with fixed \(G_{\mathrm{eff}}\)

1. Fix free constants \((\sigma_0,s,\gamma,c)\) from constitutive choice (same as multipole run).  
2. Predict \(G_{\mathrm{eff}}=\gamma s c^4/(8\pi\sigma_0)\).  
3. From exterior \(\psi\) or delay slope, infer \(M_{\mathrm{ray}}\) with **that** \(G_{\mathrm{eff}}\), not free \(M\).  
4. Compare to \(E_\star/c^2\).

**Pass:** \(M_{\mathrm{ray}}\approx E_\star/c^2\) tests multipole–ledger link (partial triad).  
**Not sufficient alone for J5** — still need \(m_{\mathrm{inertial}}\).

---

## 6. Export package for D

```text
E_star
E_unlock          # or null + unlock=ledger_proxy
m_inertial
M_ray
c
protocol          # A2_boost_energy | A3_free_drag | ...
epsilon_triad
tautology_flag    # must be false
G_eff_predicted
G_eff_used_for_M_ray
```

Score:

\[
L_{\mathrm{triad}}
= \left|\frac{m_{\mathrm{inertial}}}{E_\star/c^2}-1\right|
+ \left|\frac{M_{\mathrm{ray}}}{E_\star/c^2}-1\right|.
\]

---

## 7. Theory expectation (A)

- Scaling \(m\sim E/c^2\) expected from free-response energy / \(c^2\) (EM-mass analogy).  
- **Coefficient exactly 1** is **not** proven (kinetic_inertia_v0).  
- J5 **PASS** if numerics hit band; **FAIL** if systematic split; **DEFER** if only multipole–ledger holds and inertia marked optional for `goal2_PC3D_workable`.

O-005 allows: if J5 rigorously deferred as independent optional, orchestrator may still declare `goal2_PC3D_workable` with residuals.

---

## 8. Minimal Round-4 success for B

| Priority | Deliverable |
|----------|-------------|
| P0 | Attempt A2: \(\Delta U(v)\) vs \(\tfrac12(E_\star/c^2)v^2\) for 2–3 small \(v\) |
| P1 | \(M_{\mathrm{ray}}\) with fixed \(G_{\mathrm{eff}}\) vs \(E_\star/c^2\) |
| P2 | Export triad JSON; `tautology_flag=false` |
| Optional | A3 free-drag if A2 unstable on coarse grid |

---

## 9. Bottom line

**Non-tautological inertia** = free continuum energy (or free force integrals) respond to lock motion/acceleration **without** inserting \(m=E_\star/c^2\) into the integrator by hand.

Until that passes or is explicitly deferred, **goal2_FULL** stays open; **goal2_PC3D_workable** theory package remains MET (THEORETICAL_PACKAGE_v1).
