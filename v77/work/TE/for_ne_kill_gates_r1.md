# FOR_NE — concrete numeric kill-gates (Round 1–2)

**Agent:** TE  
**Date:** 2026-07-18  
**Round:** 1  
**Theory ref:** [`maxwell_monist_v0.md`](maxwell_monist_v0.md)  
**Constitutive:** [`constitutive_table_v0.md`](constitutive_table_v0.md)  
**Audience:** NE (primary); NM/ND for dual-channel later

---

## 0. Implementation posture

- Sandboxes **only** under `v77/work/NE/` (Python local default).  
- May cite/adapt v76 free-capacity code patterns; do **not** edit v76 logs.  
- Units: recommend \(\varepsilon=\mu=c=1\).  
- Tag exports: `sector=1`, `A_origin=free_gauge` or `E_origin=free_maxwell_lite`, `gravity_solver=none`, `em_solver=free_maxwell_lite`.  
- Dualist adversary must be **labeled** dualist (`sector=2` or `softE_em`), never soft-monist.

---

## 1. Round 1 kill-gates (must attempt)

### KG1 — Vacuum control (`D-EM-vacuum-control`)

| | |
|--|--|
| **Setup** | \(\rho_Q\equiv 0\), \(\mathbf{J}_Q\equiv 0\); quiescent IC \(\mathbf{E}=\mathbf{B}=\mathbf{0}\) (or pure gauge) |
| **Evolve / solve** | Static solve and/or short dynamical idle |
| **PASS** | \(\max|\mathbf{E}|\), \(\max|\mathbf{B}|\) ≤ floor (e.g. \(10^{-12}\) relative to a unit-charge reference run, or absolute \(<10^{-10}\) in unitless grid units) |
| **FAIL / KILL signal** | Spontaneous Coulomb field without source (implementation or ontology leak) |

### KG2 — Discrete Gauss law (`D-EM-gauss-coulomb`)

| | |
|--|--|
| **Setup** | Compact charge blob total \(Q=\int\rho_Q=1\) (smooth kernel preferred) |
| **Solve** | Electrostatic: \(-\varepsilon\nabla^2\Phi=\rho_Q\), \(\mathbf{E}=-\nabla\Phi\) **or** free relaxation of Gauss-compatible projection |
| **PASS** | Residual \(r_G=\|\nabla\cdot(\varepsilon\mathbf{E})-\rho_Q\|_\infty\) ≤ \(10^{-8}\) (SOR/FFT) or ≤ \(10^{-6}\) (low-order FD) of \(\|\rho_Q\|_\infty\); **and** flux \(\oint\varepsilon\mathbf{E}\cdot dA = Q\) within 1% on a sphere outside support |
| **FAIL** | Gauss residual \(O(1)\); flux not \(\approx Q\) |

### KG3 — Coulomb multipole \(E\sim 1/r^2\) (`D-EM-gauss-coulomb`)

| | |
|--|--|
| **Setup** | Same as KG2; sample exterior shell \(r\in[r_{\min},r_{\max}]\) outside charge support, away from BC |
| **Metric** | Fit \(\log|\mathbf{E}|\) vs \(\log r\); or \(\mathcal{R}_E(r)=4\pi\varepsilon r^2|\mathbf{E}|/|Q|\) |
| **PASS** | Prefer \(1/r^2\) over \(1/r\) and over \(1/r^3\): e.g. \(R^2_{1/r^2}>0.95\) on clean shell **or** mean \(|\mathcal{R}_E-1|<0.15\) on mid exterior shell; vacuum control still holds |
| **FAIL** | Exterior flat; wrong power; only hand two-body \(F=Q_1Q_2/r^2\) without field solve |

### KG4 — Wave speed \(c=1/\sqrt{\varepsilon\mu}\) (`D-EM-wave-c`)

| | |
|--|--|
| **Setup** | Source-free or compact pulse; set \((\varepsilon,\mu)\) so \(c_{\mathrm{th}}=1/\sqrt{\varepsilon\mu}\) known (try \(c_{\mathrm{th}}=1\) and one off-unit case e.g. \(\varepsilon=4,\mu=1\Rightarrow c=1/2\)) |
| **Evolve** | Hyperbolic Maxwell or wave equation for transverse packet / plane-wave train |
| **PASS** | Measured group or front speed \(c_{\mathrm{meas}}\) within **5%** of \(c_{\mathrm{th}}\) on both unit and off-unit constitutive points |
| **FAIL** | \(c_{\mathrm{meas}}\) independent of \(\varepsilon,\mu\); or locked to grid Courant number irrespective of constitutive \(c\) |

---

## 2. Round 1 optional / Round 2 preferred

### KG5 — Continuity / charge conservation (`D-EM-continuity`)

| | |
|--|--|
| **Setup** | Time-dep \(\rho_Q,\mathbf{J}_Q\) with prescribed continuity-satisfying current **or** lock motion proxy |
| **PASS** | \(|dQ_{\mathrm{tot}}/dt|\) ≤ tol when BC non-absorbing for charge; local \(\partial_t\rho_Q+\nabla\cdot\mathbf{J}_Q\) residual at discrete floor |
| **FAIL** | Net charge drifts without boundary flux |

### KG6 — Dualist Occam adversary (`D-EM-dualist-occam`)

| | |
|--|--|
| **Setup** | Twin: hand Coulomb \(\Phi=\sum q_i/(4\pi\varepsilon r_i)\) on fixed stage with free DOF idle (`sector=2`) vs monist free solve (`sector=1`) |
| **PASS** | Ray/force **fit may be isomorphic** (expected); monist wins only with free-origin tags + idle-DOF penalty on dualist (v76 softE lesson) |
| **FAIL as monist proof** | If NE claims monism **solely** from \(1/r^2\) fit without tags — TE will reject as softE_em |
| **Numeric deliverable** | Export both; document isomorphism residual R-EM2 |

### KG7 — Sibling ≠ ψ (`D-EM-sibling-psi`) — **V77-2 critical**

| | |
|--|--|
| **Setup** | Joint or sequential: (i) \(\rho_b>0\), \(\rho_Q=0\) neutral mass blob; (ii) \(\rho_b>0\), two opposite \(\rho_Q\) with same \(\rho_b\) pattern if possible |
| **PASS** | (i) \(\psi\) exterior nonzero (F1-3D) **and** \(\mathbf{E}\approx\mathbf{0}\); (ii) \(\mathbf{E}\) shows dipole/opposite structure while \(\psi\) remains same-sign monopole from total \(E_\star\) |
| **FAIL / K-EM5** | Neutral mass sources Coulomb \(\mathbf{E}\); or opposite charges force opposite \(\psi\) monopoles under joint laws |

### KG8 — Shared \(c\) consistency (`D-EM-sibling-psi` / Dynamics)

| | |
|--|--|
| **Setup** | Same free medium constants: EM wave \(c_{\mathrm{EM}}\) vs Dynamics free-signal \(c\) (if ND has free wave) |
| **PASS** | \(|c_{\mathrm{EM}}-c_{\mathrm{dyn}}|/c < 5\%\) when constitutive shared |
| **FAIL** | Two unrelated speeds with no frame story |

---

## 3. Demo ID ↔ gate map

| Demo ID | Gates |
|---------|-------|
| D-EM-vacuum-control | KG1 |
| D-EM-gauss-coulomb | KG2, KG3 |
| D-EM-wave-c | KG4 |
| D-EM-continuity | KG5 |
| D-EM-dualist-occam | KG6 |
| D-EM-sibling-psi | KG7, KG8 |
| D-EM-maxwell-lite | theory; numeric congruence when KG1–4 PASS |

---

## 4. Suggested NE Round-1 deliverables

1. `work/NE/sandbox_em_static.py` — KG1–KG3  
2. `work/NE/sandbox_em_wave.py` — KG4 (± off-unit \(\varepsilon\))  
3. `work/NE/outputs/round1_*.json` or `.tsv` with pass flags  
4. Log `NE-*` findings; tag `FOR_TE` if equations need revision  

**Provisional OK:** if TE late, NE may implement from this gate sheet + Maxwell textbook form with monist tags — mark `provisional` until TE congruence check.

---

## 5. Explicit non-goals for NE Round 1

- Real \(\alpha_{\mathrm{EM}}\), SI units.  
- Full particle radiation reaction.  
- Modifying `scp_sim` / `sfa.h`.  
- Replacing v76 \(\psi\) channel with \(\varepsilon(\rho_b)\) GRIN.

---

## 6. Pass bar for EM focus toward V77-1

**EM theory+numeric draft PASS** if:

- TE package exists (done), **and**  
- NE achieves **KG1 ∧ KG2 ∧ (KG3 ∨ KG4)** with honest FAIL logs if not.

**Toward V77-2:** add KG7 (sibling) with ψ co-present.
