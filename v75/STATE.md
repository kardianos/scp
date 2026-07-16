# v75 Current State — Multi-fabric atoms, P/N, C₁₂ path

**As of:** 2026-07-16  
**Era:** U(1) complex Cosserat + multi-fabric Option B (B1 locked → **B2 unlock for P/N**)  
**Goal doc:** `C12_ATOM_GOAL.md` · **Lab notebook:** `FINDINGS.md` (F11–F19) · **Status:** `STATUS.md`

This file is the **theory + ops snapshot** for this version of the work: equations, charge assignment, what was measured, and what is frozen. Someone reading only this should understand *what the model is now* and *what experiments established it*, without replaying the whole trial-and-error trail.

---

## 1. Field content (frozen package)

### 1.1 Matter

Three complex Cosserat sectors (fabrics), each with the standard 12 real fields per fabric (φ + θ, Re/Im, and velocities), when `n_fabrics=3`:

| Fabric | Symbol | Role |
|--------|--------|------|
| **C** | \(\Phi_C^a\), \(\Theta_C^a\) | Nuclear **bag** / binding medium |
| **Q** | \(\Phi_Q^a\), \(\Theta_Q^a\) | Nuclear **EM charge** carrier |
| **L** | \(\Phi_L^a\), \(\Theta_L^a\) | Light opposite-charge sector |

Shared abelian gauge field \(A_\mu\) (lattice links \(U=\mathrm{e}^{i\theta}\), electric field \(E\)) when `complex_gauge=1`.

Potential (each fabric’s private bag product):

\[
s_f = \prod_{a=0}^{2} |\Phi_f^a|^2,\qquad
V_t(s)=\frac{\mu}{2}\frac{s}{1+\kappa s},
\quad m^2=2.25,\;\mu=-41.345,\;\kappa=50.
\]

Standard particle package: `eta=0`, `m_theta=1.6` (θ-drain closed), `g_gauge=0.05`, absorbing BC.

### 1.2 Multi-fabric EM charge (shared \(A\))

Noether density on fabric \(f\) (φ+θ sectors):

\[
\rho_f = \sum_a \Big(
\Phi_f^{a\,\mathrm{Re}}\,\dot\Phi_f^{a\,\mathrm{Im}}
- \Phi_f^{a\,\mathrm{Im}}\,\dot\Phi_f^{a\,\mathrm{Re}}
+ \text{(same for \(\Theta\))}
\Big).
\]

**EM source for Gauss / currents:**

\[
\rho_{\mathrm{em}}
= q_C\,\rho_C + q_Q\,\rho_Q + q_L\,\rho_L.
\]

Config weights (multiples of \(g\) in the covariant derivative / current):

\[
\boxed{q_C=0,\quad q_Q=+1,\quad q_L=-1}
\quad\text{(physical coupling scale \(g=g_{\mathrm{gauge}}=0.05\))}
\]

Covariant transport uses \(U^{q}=\mathrm{e}^{i q\,\theta}\) per fabric (must not drop \(q\) or opposite-charge forces null).

**Gauss law (discrete):** \(\nabla\cdot E - g\,\rho_{\mathrm{em}} \approx 0\) after projection; `gauss_max` must stay \(\sim 10^{-13}\).

**Diagnostics proxies**

| Name | Meaning |
|------|---------|
| `Q_phi` | Integrated Noether on **primary/C** arrays (bag bookkeeping) |
| `Q_flux` | Gauss-cube EM charge proxy \(\propto \int(\nabla\cdot E)/g\) |
| `Q_C,Q_Q,Q_L,Q_em` | CPU multi-fab integrals (when enabled) |
| Track massL / Ql | L-sector survival from SFA columns |

---

## 2. B1 vs B2 dynamics

### 2.1 B1 (Shape β) — isolation + Coulomb proven (F11–F16)

\[
\Phi_Q(x,t)\equiv\Phi_C(x,t)
\quad\text{(copy after each force / half-step)}
\]

- Private bags: \(s_C\) and \(s_L\) never mix → **no C–L merge** at head-on.  
- With lock, effective heavy EM charge is \(q_Q\rho_C\) (Q mirrors C).  
- Seed rule for opposite fabric charges: **same sign** of Noether \(\omega\) on C and L (fabric \(q\) supplies the EM sign flip).

### 2.2 B2 unlock — **required for P/N** (F17+)

\[
\texttt{mf\_lock\_CQ}=0
\quad\Rightarrow\quad
\Phi_Q\text{ evolves independently of }\Phi_C.
\]

Init options:

| Mode | Config | Use |
|------|--------|-----|
| Load Q seed | `init_sfa_Q=...` | Selective charge (protons only) |
| One-shot copy | `mf_seed_Q=1` | Proton-like start without ongoing lock |
| Leave Q zero | default unlock | **Neutron-friendly** (no EM source) |

**Implementation note:** bag \(V_t\) still applied per fabric in the current force path (Q self-binds if seeded). Preferred long-term: bag on C only; not yet required for the P/N claim.

---

## 3. Proton / neutron analogs (P2.0 freeze)

### 3.1 Definition (measured)

| Species | Seed content | EM |
|---------|--------------|-----|
| **Proton analog** | \(\Phi_C\) bag + \(\Phi_Q\) co-located | \(\rho_{\mathrm{em}}\approx q_Q\rho_Q\neq 0\) |
| **Neutron analog** | \(\Phi_C\) bag only (\(\Phi_Q=0\) at that site) | \(\rho_{\mathrm{em}}\approx 0\) |
| **Light / e-analog** | \(\Phi_L\), \(q_L=-1\) | opposite to nuclear \(Q\) |

Bookkeeping:

\[
Z = \text{count of proton seeds (or proxy \(Q_{\mathrm{em}}/Q_{\mathrm{em,1p}}\))},\quad
N = \text{count of neutron seeds},\quad
A = Z+N.
\]

**Isotope control:** vary \(N\) at fixed \(Z\) (and fixed L count \(=Z\)); nuclear \(Q_{\mathrm{em}}\) (and L charge) must stay fixed while bag mass / \(Q_\phi\) change.

### 3.2 What does **not** define P/N

| Approach | Result |
|----------|--------|
| Flavored multiplet \(\Delta\omega\) alone (v71 style) | Stable multiplet, **same-sign EM** — not neutral n |
| Equal-\(f\) cancel multiplet \(\omega=(+w,+w,-w)\) | Residual \(Q\sim Q/3\), not zero |
| B1 lock only | Every nuclear lump sources EM |

### 3.3 Seed tool

`bin/gen_pn_core`:

```text
gen_pn_core N L profNuc omega out_C out_Q out_L \
  nZ  xz yz zz ... \
  nN  xn yn zn ... \
  nL  xl yl zl ... \
  [profL omegaL]
```

- C gets all nuclear balls (Z+N).  
- Q gets **only** Z proton centers.  
- L optional; hierarchy via trailing `profL omegaL` (e.g. ω=1.46 for light).

---

## 4. Experimental record (this version)

### 4.1 Closed gates

| ID | Result |
|----|--------|
| F11–F13 | B1 isolation + Coulomb (U^q) |
| F15 | B4 self-tune packaging (soft θ) |
| F16 | Pair kinematics; Z6+L6 PASS_park (pre-P/N, B1) |
| **F17** | B2 P/N: n \(Q_\mathrm{flux}=0\); p \(Q_\mathrm{flux}\approx Q_\phi\) |
| **F18** | Z2 park PASS; L count=Z; isotope \(Q_\mathrm{flux}\) identical |
| **F19** | Z6: isotope \(Q_\mathrm{flux}\) identical; Z6N6 park PASS; Z6N0 soft; L −12.5% |

### 4.2 Key measured numbers (F17–F19)

**F17 (N=96, single-ball B2)**

- Neutron: \(Q_\mathrm{flux}=0\), \(E_{\mathrm{em}}=0\), bag \(Q_C\approx 210\).  
- Proton: \(Q_\mathrm{flux}\approx Q_\phi\approx 209\).

**F18 (N=128, Z=2)**

- Nuclear: \(Q_\mathrm{flux}(Z2N0)=Q_\mathrm{flux}(Z2N2)=372\) (end).  
- Atom L=2: massL −0.6%; identical L for +N.

**F19 (N=192 L=48, Z=6)**

- Nuclear: \(Q_\mathrm{flux}(Z6N0)=Q_\mathrm{flux}(Z6N6)=990\) (end).  
- Z6N6: \(c_{Q,\mathrm{park}}=0.046\) PASS; Z6N0: \(0.184\) soft.  
- Atom L=6: massL −12.5%; **identical** L evolution for isotope pair.  
- **Visual:** multi-ball nuclear seed **merges to a single droplet** by t∼400 (parked core), not multi-center nucleons.

### 4.3 Grid (boundary safety)

\[
R_L + 2 r_{\mathrm{core}} + \Delta_{\mathrm{buf}} \lesssim R_{\mathrm{damp}}=L-\mathrm{damp\_width}.
\]

F19 freeze: `N=192`, `L=48`, `damp_width=5`, nuclear octa R=8 / n R=5.5, L shell R=22.  
**Do not** put Z6 L@R=22 in F18’s L=32 box (sponge eats L). Details: `PN_GRID.md`.

---

## 5. Scorecard definitions

**Nuclear park (no L):**

\[
c_{Q,\mathrm{park}}=\frac{|Q_{\mathrm{mid}}-Q_{\mathrm{end}}|}{|Q_{\mathrm{mid}}|},\quad
c_{Q\mathrm{em},\mathrm{park}}=\frac{\big||Qf_{\mathrm{mid}}|-|Qf_{\mathrm{end}}|\big|}{\max(|Qf_{\mathrm{mid}}|,1)}.
\]

PASS_nuc: \(c_{Q,\mathrm{park}}\le 0.15\), \(c_{Q\mathrm{em},\mathrm{park}}\le 0.20\), Gauss floor, \(|Q_{\mathrm{end}}|>0.5|Q_{\mathrm{mid}}|\).

**Atom:** nuclear park + L track: \(c_L\le 0.15\), massL_end > 0.5 massL_0.  
**Do not** treat late **net** \(Q_\mathrm{flux}\) drop alone as L death when massL/Ql hold (net = nuclear EM − L).

Tool: `v75/analysis/score_pn_park.py`.

---

## 6. Kernel / config surface (operators)

```
n_fabrics=3
mf_lock_CQ=0          # B2 for P/N
mf_stage=2
q_C=0  q_Q=1  q_L=-1
init=sfa
init_sfa=..._C.sfa
init_sfa_Q=..._Q.sfa  # protons only; empty Q file for pure-n cores
init_sfa_L=..._L.sfa  # zero or light shell
complex_phi=1 complex_gauge=1 g_gauge=0.05 m_theta=1.6 eta=0
```

Code touchpoints (this version):

- `sfa/sim/scp_config.h` — `init_sfa_Q`, `mf_seed_Q`  
- `sfa/sim/scp_sim.c` / `.cu` — B2 force path, Q init, CPU `Q_em` diags; GPU Gauss uses Q as primary when unlocked (`q_C=0`)  
- `sfa/seed/gen_pn_core.c` — Z/N/L seed writer  

**Protected:** do not modify `sfa.h` without explicit auth; sim kernel only when authorized (P/N work was authorized).

---

## 7. Ops notes

- **Runner:** `scp-runner` MCP. **Never** put multi-minute `sleep` inside `sim_exec` (blocks SSH/MCP). Campaign = nohup script; poll with instant `cat`/`ps`/`nvidia-smi`.  
- **Status snippet:** `v75/analysis/remote_status_snippet.sh`  
- **Data:** `/space/scp/v75/pn/` (F17), `.../p2/` (F18), `.../z6/` (F19 SFAs + renders)  
- **Scores:** `v75/results/pn_z6/`  
- **Renders:** `images/C6_atom_sheet.png` · `/space/scp/v75/pn/z6/renders/`  
- Instance used: `v75f16` (V100-16GB; N=192 multi-fab ~8 GB GPU)

---

## 8. Open problems (honest)

1. Multi-ball nuclei **fuse to a droplet** by t∼400 — not a multi-center nuclear crystal.  
2. Z6N0 park softer than Z6N6; L loss larger at Z=6 than Z=2.  
3. Shell-radius diagnostic still needed (COM D≈0 is not multi-L structure).  
4. Long-T time-stable **visual C₁₂** (ideal goal) not yet claimed.  
5. Q self-bag under B2 still on; ε_CQ portal not required yet.

---

## 9. Document map

| File | Role |
|------|------|
| **This file** | Equations + freeze + version snapshot |
| `C12_ATOM_GOAL.md` | Goal phases 1–3 + P/N requirements |
| `MULTIFABRIC_SPEC.md` | Original B1/B2 design |
| `PN_EXPERIMENT.md` / `PN_GRID.md` | P/N matrix + grid sizing |
| `FINDINGS.md` | Chronological F11–F19 |
| `STATUS.md` | Live checklist |
| `FUTURE.md` | Repo pointer to C₁₂ goal |

**Version one-liner:**  
*B2 multi-fabric with \(q_C=0,q_Q=+1,q_L=-1\): proton = C+Q, neutron = C-only; isotopes keep \(Q_{\mathrm{em}}\propto Z\) through Z=6; multi-ball cores park as droplets; C₁₂ long-T package still open.*
