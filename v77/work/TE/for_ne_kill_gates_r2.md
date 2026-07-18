# FOR_NE — full Maxwell kill-gates (Round 2)

**Agent:** TE  
**Date:** 2026-07-18  
**Round:** 2  
**Theory:** [`full_maxwell_monist_v0.md`](full_maxwell_monist_v0.md)  
**Joint:** [`dual_channel_joint_v0.md`](dual_channel_joint_v0.md)  
**R1 gates (still valid):** [`for_ne_kill_gates_r1.md`](for_ne_kill_gates_r1.md) — KG1–4 **PASS** stamped congruent

---

## 0. R1 stamp (read first)

TE confirms NE R1 equations match Maxwell-lite subset:

- FG-Q / M3: \(-\varepsilon\nabla^2\Phi=\rho_Q\), \(\mathbf{E}=-\nabla\Phi\)  
- Wave: \((\partial_{tt}-c^2\nabla^2)A=0\), \(c=1/\sqrt{\varepsilon\mu}\)  
- Gauss continuum identity  

→ set `te_equation_match=yes` for KG1–4 scope.  
→ keep `full_maxwell_claim=false` until gates below PASS.

---

## 1. Implementation posture R2

| Item | Recommendation |
|------|----------------|
| Solver | Yee leapfrog 2D TE/TM **or** 3D; tags `em_solver=free_maxwell_yee` |
| Units | \(\varepsilon=\mu=c=1\) default; **one off-unit** wave-with-B run |
| BC | Periodic for divB/Faraday loops; PEC/absorbing for radiation if needed |
| Charge | Cont-satisfying \(\mathbf{J}_Q\) or static \(\rho_Q\) with \(\mathbf{J}_Q=0\) |
| Dual-channel | Prefer coordinate with NM for KG7; NE may do sequential ψ import |
| Dualist twin | sector=2 hand Biot–Savart / Coulomb; document iso residual |
| No | `scp_sim` edits; real \(\alpha\) |

---

## 2. New Round-2 kill-gates (full Maxwell)

### KG-F1 — Divergence-free \(\mathbf{B}\) (`D-EM-divB`)

| | |
|--|--|
| **Setup** | Yee IC with \(\mathrm{div}_h\mathbf{B}=0\) (e.g. pure transverse pulse or \(\mathbf{B}=\mathrm{curl}_h\mathbf{A}\)) |
| **Evolve** | ≥ few hundred steps source-free or with Cont-safe \(\mathbf{J}_Q\) |
| **PASS** | \(\max|\mathrm{div}_h\mathbf{B}| / (\|\mathbf{B}\|_\infty/h+\epsilon) \le 10^{-8}\) (double) or \(\le 10^{-5}\) (float / low order) **and** does not secularly grow |
| **FAIL** | divB grows as free dynamical error (not BC artifact alone) |

### KG-F2 — Faraday (`D-EM-faraday`)

| | |
|--|--|
| **Setup** | 2D or 3D loop \(C=\partial S\); measure \(\mathcal{E}=\oint_C\mathbf{E}\cdot d\boldsymbol{\ell}\) and \(\Phi_B=\int_S\mathbf{B}\cdot d\mathbf{A}\) |
| **PASS** | \(\bigl|\mathcal{E}+\partial_t\Phi_B\bigr| / (|\mathcal{E}|+|\partial_t\Phi_B|+\epsilon) \le 5\%\) time-averaged on smooth pulse (or residual at discrete Stokes identity floor for pure Yee) |
| **FAIL** | Faraday violated at \(O(1)\) while claiming free Maxwell |

### KG-F3 — Ampère–Maxwell / wave with \(\mathbf{B}\) (`D-EM-wave-EB`, `D-EM-ampere`)

| | |
|--|--|
| **Setup** | Source-free transverse EM packet (coupled \(\mathbf{E},\mathbf{B}\)), not scalar \(A\) proxy alone |
| **PASS** | (i) \(c_{\mathrm{meas}}\) within **5%** of \(1/\sqrt{\varepsilon\mu}\); (ii) off-unit e.g. \(\varepsilon=4\Rightarrow c=1/2\) also within 5%; (iii) \(|\mathbf{B}|\sim |\mathbf{E}|/c\) (SI-like) or \(|\mathbf{B}|\sim|\mathbf{E}|\) when \(c=1\) within 15% peak for plane-wave polarization |
| **FAIL** | Speed ignores \(\varepsilon,\mu\); or only E without B dynamics |

### KG-F4 — Continuity + Gauss preservation (`D-EM-continuity`)

| | |
|--|--|
| **Setup** | Compact \(\rho_Q\) with \(\mathbf{J}_Q\) satisfying discrete Cont (moving blob or prescribed divergence-free flow of charge) |
| **PASS** | \(|Q_{\mathrm{tot}}(t)-Q_{\mathrm{tot}}(0)|/|Q_0| \le 10^{-4}\) interior (no BC leak); Gauss residual \(\|\mathrm{div}(\varepsilon\mathbf{E})-\rho_Q\|\) stays ≤ R1 thresholds or after cleanse |
| **FAIL** | Charge drifts without flux; Gauss residual \(O(1)\) secular growth |

### KG-F5 — Static Coulomb still holds under full solver (`D-EM-gauss-coulomb` regression)

| | |
|--|--|
| **Setup** | Static \(\rho_Q\), \(\mathbf{J}_Q=0\); evolve Yee to relaxation **or** show static projection satisfies M3 |
| **PASS** | Exterior \(E\sim 1/r^2\) (3D) or consistent 2D log \(E\sim 1/r\) **labeled as 2D**; vacuum control preserved |
| **FAIL** | Full Maxwell dynamic breaks quasistatic Coulomb that lite already passed |

### KG-F6 — Displacement current necessity (anti-magnetostatics-only)

| | |
|--|--|
| **Setup** | Compare free Maxwell vs **killed** \(\partial_t\mathbf{E}\) term (magnetostatics / incomplete Ampère) on radiating packet |
| **PASS** | Full M4 radiates / propagates at \(c\); incomplete Ampère fails wave or causality diagnostic |
| **FAIL as monism kill** | Only if **incomplete** law is required for monist tags to work — then K-EM9 |

---

## 3. Dual-channel gates (V77-2; with NM)

### KG7 — Sibling ≠ ψ (repeat R1, now mandatory for V77-2)

| Case | PASS |
|------|------|
| \(\rho_b>0,\rho_Q=0\) | \(\psi\) exterior nonzero (F1); \(\max|\mathbf{E}|\) ≤ floor vs unit-charge ref |
| Two opposite \(Q\), both \(\rho_b>0\) | \(\mathbf{E}\) dipole-like; \(\psi\) same-sign monopole from total \(E_\star\) |

### KG8 — Shared \(c\)

| | |
|--|--|
| **PASS** | \(c_{\mathrm{EM}}\) from KG-F3 equals Dynamics/path \(c\) within 5% when constitutive shared |
| **Cross** | ND free-wave already \(v\approx c\); report joint number |

### KG9 — Joint dual multipole (DUAL-0)

| | |
|--|--|
| **Setup** | One grid: solve F1-3D for \(\psi\) and M3 for \(\mathbf{E}\) with Supp |
| **PASS** | Both exterior multipoles PASS their individual bars; tags `dual_channel=1` |
| **Owner** | NE and/or NM; one export JSON is enough for TU |

---

## 4. Pass bars for claims

| Claim | Required gates |
|-------|----------------|
| **full_maxwell_claim=true** | KG-F1 ∧ KG-F2 ∧ KG-F3 ∧ (KG-F4 ∨ static Cont) ∧ KG-F5 |
| **V77-2 dual-channel numeric** | KG7 ∧ KG9 (+ KG8 preferred) |
| **D-EM-maxwell-full LIVE** | full_maxwell_claim + TE congruence stamp R2 |
| **Regression** | R1 KG1–4 must not regress |

---

## 5. Suggested NE R2 deliverables

1. `work/NE/sandbox_ne_r2_yee.py` — 2D TE or TM Yee minimum  
2. Optional: Faraday loop diagnostic; divB monitor  
3. `work/NE/outputs/r2_result.json` with gate flags  
4. Dual-channel: either NE joint or consume NM dual0 fields with Maxwell Φ tags  
5. Log NE-*; `FOR_TE` if Yee identities force equation revision  

**2D honesty:** if Faraday/wave-B are 2D, say so; keep 3D Coulomb from R1 as complementary, not pretend 3D full Yee until done.

---

## 6. Demo ID ↔ gate map

| Demo ID | Gates |
|---------|-------|
| D-EM-divB | KG-F1 |
| D-EM-faraday | KG-F2 |
| D-EM-wave-EB / D-EM-ampere | KG-F3, KG-F6 |
| D-EM-continuity | KG-F4 |
| D-EM-gauss-coulomb | KG-F5 + R1 KG2–3 |
| D-EM-maxwell-full | full bar §4 |
| D-EM-sibling-psi | KG7, KG8 |
| D-DUAL-channel | KG7, KG9 |

---

## 7. Non-goals R2

- Real \(\alpha_{\mathrm{EM}}\), SI  
- Radiation reaction on accelerating locks (R-EM3)  
- Non-abelian / Cosserat kernel  
- Replacing C-ψ with \(\varepsilon(\rho_b)\) GRIN
