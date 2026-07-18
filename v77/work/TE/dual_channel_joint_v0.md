# Dual-channel constitutive joint note (C-ψ + C-A)

**Agent:** TE  
**Date:** 2026-07-18  
**Round:** 2 (written) · **3 (composition stamped)**  
**For:** TU (WORLD/G2), NM (joint sandbox), NE (sibling gates), TM (lock dual-source)  
**Refs:** `full_maxwell_monist_v0.md` (**FROZEN**), `constitutive_table_v0.md`, TM `lock_ontology_v0.md` DUAL-0, TU WORLD  
**Milestone:** **V77-2** dual-channel free medium — **MET** (parent O-004; NM joint + NE full Maxwell)  
**Composition:** see [`dual_channel_composition_r3.md`](dual_channel_composition_r3.md) — full Yee + DUAL-0 lite compose as one C-A channel

---

## 0. Purpose

State **one joint constitutive package** so:

1. TU can score multi-focus unity without three private ontologies.  
2. NM/NE can implement a **single medium** with \(\rho_b\to\psi\) and \(\rho_Q\to\) gauge.  
3. TE-IA1 sibling tests are concrete kill-gates, not slogans.

---

## 1. Joint state (minimal monist dual-channel)

On one continuum region (box):

\[
\bigl\{\,
\rho_f,\;\rho_b,\;\rho_Q,\;\mathbf{J}_Q,\;
\psi,\;
\mathbf{E},\;\mathbf{B},\;
\sigma,\;\varepsilon,\;\mu
\,\bigr\}
\]

| Symbol | Channel | Meaning |
|--------|---------|---------|
| \(\rho_f,\rho_b\) | budget | free / bound occupancy |
| \(\rho_Q,\mathbf{J}_Q\) | C-A source | gauge ledger + current on locks |
| \(\psi\) | **C-ψ** | free-capacity potential |
| \(\mathbf{E},\mathbf{B}\) | **C-A** | free-gauge stress |
| \(\sigma\) | C-ψ constitutive | free conductivity |
| \(\varepsilon,\mu\) | C-A constitutive | free permittivity / permeability |

**Not joint state:** independent metric \(g_{\mu\nu}\) as ontology; foreign \(G\) and \(\alpha\) as second-sector constants.

---

## 2. Joint laws (static DUAL-0 — V77-2 minimum)

\[
\boxed{
\begin{aligned}
\rho_f+\rho_b&=\rho_{\mathrm{tot}},\\
\mathrm{supp}(|\rho_Q|)&\subseteq\mathrm{supp}(\rho_b),\\
-\nabla\cdot\bigl(\sigma(\rho_f)\nabla\psi\bigr)&=s\,\rho_b,\\
\nabla\cdot(\varepsilon\mathbf{E})&=\rho_Q,\\
\mathbf{E}&=-\nabla\Phi\quad\text{(static)},\\
c&=\frac{1}{\sqrt{\varepsilon\mu}}=c_{\mathrm{path}}.
\end{aligned}
}
\tag{DUAL-0}
\]

Vacuum linear: \(\sigma=\sigma_0\), \(\varepsilon=\varepsilon_0\).

Exterior (compact isolated sources):

\[
\psi\sim\frac{s E_\star}{4\pi\sigma_0 r},\qquad
\Phi\sim\frac{Q}{4\pi\varepsilon_0 r},\qquad
E\sim\frac{|Q|}{4\pi\varepsilon_0 r^2}.
\]

**TM force sketch (quasistatic):**

\[
\mathbf{F}_i=\mathbf{F}_i^{\psi}[M_j,\nabla\psi]+\mathbf{F}_i^{C}[Q_j,\mathbf{E}].
\]

Path-cost: same-sign attraction. Coulomb: \(Q_i Q_j\) sign structure.

---

## 3. Joint laws (dynamic DUAL-1 — full Maxwell branch)

When NE full Maxwell is live:

\[
\boxed{
\begin{aligned}
&\text{C-A: M1–M4 + Cont (full\_maxwell\_monist\_v0)},\\
&\text{C-ψ: F1-3D quasistatic (or TD time-dep free response)},\\
&\mathbf{J}_Q=\rho_Q\mathbf{v}_{\mathrm{lock}}+\mathbf{J}_{\mathrm{int}},\\
&c=1/\sqrt{\varepsilon\mu}=c_{\mathrm{path}}=c_{\mathrm{free-wave}}.
\end{aligned}
}
\tag{DUAL-1}
\]

**R2 priority:** DUAL-0 numeric first (static dual Poisson/F1 on one grid). DUAL-1 after Yee green.

---

## 4. Constitutive joint constraints (JC — freeze for TU)

| ID | Constraint | FAIL means |
|----|------------|------------|
| **JC1** | Single free locality \(c\): path-cost \(c\) \(=\) \(1/\sqrt{\varepsilon\mu}\) \(=\) free-wave \(c\) | Split ontology speeds → K2 path |
| **JC2** | Budget \(\rho_f+\rho_b=\rho_{\mathrm{tot}}\) (or integral form); \(u_{\mathrm{EM}}\) in free ledger | Double-count / ghost energy |
| **JC3** | Sources independent: \(\rho_b\not\equiv\rho_Q\); Supp default | Sibling collapse or pure dualist charge fluid |
| **JC4** | \(G_{\mathrm{eff}}\) from \((\gamma,s,\sigma_0,c)\); Coulomb from \(\varepsilon\) — both free constitutive | Foreign second-sector constants as ontology |
| **JC5** | Sector tags: `phi_origin=free_*`, `E_origin=free_maxwell_*`, `gravity_solver=none` | softE / soft dualist |
| **JC6** | TE-IA1: \(\psi\not\equiv\Phi_{\mathrm{EM}}\) as fields | Forced identification kills dual-channel world |
| **JC7** | Neutral mass sources \(\psi\) not \(\mathbf{E}\); opposite \(Q\) do not reverse \(\psi\) monopole | K-EM5 / DS4 |

---

## 5. Recommended joint sandbox units (NE+NM)

| Choice | Value | Note |
|--------|-------|------|
| \(\varepsilon,\mu,c\) | \(1,1,1\) | Match Dynamics |
| \(\sigma_0,s\) | as v76 sandbox | Do not retune to fake Coulomb |
| \(\gamma\) | as path-cost demos | Independent of \(\varepsilon\) |
| Grid | same Cartesian mesh for \(\psi\) and \(\Phi\)/\(\mathbf{E}\) | One medium |
| Locks | ≥2 compact blobs; cases: ++, +−, neutral+\(Q\) | Sibling tests |
| Tags | `dual_channel=1`, `sector=1` | Export JSON |

**Hierarchy scan (TM Ξ):** vary \(k_C/\kappa_\psi\) or \(|Q|/M\) — do **not** delete a channel.

---

## 6. Numeric kill matrix (who owns)

| Gate | Owner | PASS criterion (short) |
|------|-------|------------------------|
| KG7 sibling | NE/NM | Neutral: \(E\approx 0\), \(\psi\neq 0\); opposite \(Q\): \(E\) dipole, \(\psi\) monopole |
| KG8 shared \(c\) | NE+ND | \(c_{\mathrm{EM}}\approx c_{\mathrm{dyn}}\) |
| G1 dual0 | NM | both multipoles on one grid |
| G2 forces | NM | same-sign \(\psi\) attract; Coulomb sign flip |
| Faraday/divB | NE | full Maxwell R2 gates |
| Occam | NE/TU | monist tags vs sector=2 twin |

Detail: TE `for_ne_kill_gates_r2.md`; TM `FOR_NM_gates_v0.md`.

---

## 7. Dualist strip (joint-only)

| ID | Forbidden joint pattern |
|----|-------------------------|
| DJ1 | Solve \(\psi\) with gravity_solver, then EM with foreign stage, stitch forces |
| DJ2 | Identify \(\psi=\Phi\) and call it “unification” |
| DJ3 | Two grids with different \(c\) as ontology |
| DJ4 | Charge without Supp and without free-charge free mode theory |
| DJ5 | Path-cost from local \(n(\rho_f)\) while EM is free Maxwell (mixed dead branch) |

---

## 8. V77-2 scorecard (for TU) — R3 update

| Item | Theory | Numeric |
|------|--------|---------|
| C-ψ LIVE | v76 | v76 |
| C-A lite (static M3) | full_maxwell §3.5 **FROZEN** | NE R1 LIVE_PASS |
| C-A full M1–M4 | full_maxwell **FROZEN** | NE R2 LIVE_PASS; `te_equation_match_full=yes` |
| DUAL-0 joint | TE+TM **FROZEN** | NM R2 KG7/8/9 **PASS** |
| JC1–JC7 written | **YES** | NM/NE stress **PASS** (JC1–5) |
| V77-2 | — | **MET** (O-004 parent) |
| DUAL-1 single-box Yee+ψ | sketch | optional residual RC1 |

---

## 9. Bottom line

**Dual-channel monism = one continuum, two free constitutive laws, two ledgers, one \(c\), Supp on charge, TE-IA1.**  
Static DUAL-0 is the V77-2 minimum (**MET**). Full Maxwell dynamics LIVE as same C-A channel (NE R2).  
Composition with \(\psi\): [`dual_channel_composition_r3.md`](dual_channel_composition_r3.md). No ontology fork.
