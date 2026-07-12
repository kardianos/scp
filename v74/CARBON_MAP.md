# v74 Stage 0 — Carbon Mapping Spec

**Date**: 2026-07-11  
**Status**: [spec + Q-scale survey] — first fusion runs queued (Stage 2A)  
**Standing goal**: CLAUDE.md § Standing Goal — Carbon from fabric  
**Branch data**: `v69/theory/gscan.tsv` (g=0.05)

This document defines what counts as **carbon** in the SCP fabric theory
before any new nuclear runs. Parallels to real \(^{12}\)C are **structural**,
not quantitative (CONCEPT.md).

---

## 1. What carbon means here

| Real carbon | Fabric analog | Conserved bookkeeping |
|-------------|---------------|------------------------|
| Atomic number \(Z=6\) | Six units of **nucleon inventory** in formation history (or six light free nucleons fused) | Formation history + total \(Q\) |
| Mass number \(A=12\) | Twelve nucleon fusions (when charge budget allows) | \(A\) is history; fused droplet forgets \(A\) (v71) |
| Nucleus | Liquid-drop **fused charge droplet** on the gauged Q-ball branch | Final \(Q\), \(E/Q\), on-branch match |
| Electrons (6) | Light **opposite-charge** solitons (Stage 3) — **not in this stage** | — |
| Atom | Nucleus + 6 light opposite charges Coulomb-bound via \(A\) (Stage 4) | — |

**Stage 2A product name:** *fabric carbon nucleus* (bare), not atom.

Until Stage 3 (electron sector) exists, success stops at a stable or
characterised nuclear droplet with carbon-scale inventory.

---

## 2. Nucleon definition (Stage 1 baseline)

Standard **nucleon** = gauged diagonal-U(1) Q-ball, all three \(\Phi_a\) equal,
\(\Theta=0\), \(\eta=0\), \(m_\theta=1.6\), \(g=0.05\) unless noted.

| label | \(\omega\) | \(Q_N\) | \(E/Q\) | \(r_{1/2}\) | profile | use |
|-------|------------|---------|---------|-------------|---------|-----|
| standard (v71 he2/li3) | 1.42 | 311.46 | 1.4649 | 3.80 | `profiles/f_w142_g005.txt` | He/Li history; too heavy for A≥4 at g=0.05 |
| light (v71 li3s) | 1.46 | 114.13 | 1.5187 | 2.63 | `profiles/f_w146_g005.txt` | default light nucleon for carbon campaign |
| lightest VK-stable | 1.4825 | 90.01 | 1.5321 | 2.34 | (shoot on demand) | minimize \(A\cdot Q_N\) |

VK-stable branch at g=0.05: \(\mathrm{d}Q/\mathrm{d}\omega < 0\) for
\(\omega \in [1.40624, 1.4825]\). Above \(\omega\approx 1.485\), sign flips
(unstable thin-wall edge) — do not seed free nucleons there.

---

## 3. Q-scale survey (Stage 2A prep) — hard bound

At **\(g=0.05\)**, static gauged branch has **\(Q_\mathrm{max} = 921\)**
(\(\omega_\mathrm{min}=1.40624\)). Lightest free nucleon on the stable branch
has **\(Q_N^\mathrm{min} \approx 90\)**.

| free nucleons \(A\) | min inventory \(A\cdot Q_N^\mathrm{min}\) | vs \(Q_\mathrm{max}\) |
|---------------------|------------------------------------------|------------------------|
| 2 (He) | ~180 | mid-branch |
| 3 (Li) | ~270 | mid-branch |
| 6 (“Z-carbon”) | ~540 | mid-branch ✓ |
| 8 | ~720 | high branch ✓ |
| **12 (\(A=12\))** | **~1080** | **always super-critical** |

**Theorem of the charge budget (g=0.05):**  
No choice of free gauged nucleons on the stable branch can fuse twelve of them
into a **subcritical** static droplet. \(12\times Q_N^\mathrm{min} > Q_\mathrm{max}\).

Consequences:

1. **True liquid-drop \(A=12\) at g=0.05** always merges then **evaporates
   toward \(Q_\mathrm{max}\)** (li3 mechanism) — still a valid experiment
   (carbon-scale hot droplet → cap).
2. **Stable parked “carbon” at g=0.05** requires either fewer free nucleons
   (map **\(Z=6\)** as \(A=6\) light fusions) or a **lower \(g\)** (larger
   \(Q_\mathrm{max}\)).
3. At **\(g=0.02\)**, \(Q_\mathrm{max}\approx 5300\) → \(12\times 311 < Q_\mathrm{max}\)
   (standard nucleons allowed). Deferred until g=0.02 profiles are shot.

### Recommended inventory targets (g=0.05)

| run id | recipe | expected seed \(Q\) | expected fate | role |
|--------|--------|---------------------|---------------|------|
| **c6_light** | 6× light (\(\omega=1.46\), \(Q_N=114\)) co-phase | ~750–850 (w/ interference) | park mid-branch \(Q\sim 600\!-\!700\), \(E/Q\sim 1.44\!-\!1.47\) | **primary stable carbon nucleus (Z=6 map)** |
| **c12_light** | 12× light co-phase | ~1400–1600 | merge → evaporate toward 921 | **A=12 supercritical control** |
| c8_light | 8× light | ~1000+ | near-cap / mild super | optional high-A edge |
| c12_g02 | 12× standard at g=0.02 | ~3700+ | park if \(Q_\mathrm{max}\) holds | true A=12 (later) |

Branch prediction helper (g=0.05): interpolate `gscan.tsv` for final \(E/Q(Q)\).

---

## 4. Observables (success criteria)

### Must measure (every nuclear run)

1. **Global \(Q(t)\), \(E(t)\), \(E/Q(t)\)** from diag.tsv  
2. **On-branch test:** final \(E/Q\) vs gscan prediction at measured \(Q\) (he2: Δ0.1%)  
3. **Mass defect:** \(1 - (E/Q)_\mathrm{final}/(E/Q)_N\) relative to free nucleon  
4. **Single droplet:** one connected charge cluster (fragmentation / sfa cluster tools)  
5. **Per-component \(Q_a\)** remain \(Q/3\) (quark bookkeeping through fusion)  
6. **Gauss floor** `gauss_max` ~ 1e-13 (implementation tripwire)

### Carbon-specific labels

| if… | we call it… |
|-----|-------------|
| 6 light nucleons → stable on-branch droplet | **Z-carbon nucleus** (Stage 2A success) |
| 12 light nucleons → single droplet evaporating to \(Q_\mathrm{max}\) | **A=12 hot carbon** (budget-limited; not a parked isotope) |
| multi-center static bound state with \(A\) retained | **multi-center carbon** (Stage 2B; not yet available) |
| nucleus + 6 opposite-charge light partners | **carbon atom** (Stage 4; blocked on Stage 3) |

### Explicit non-goals for Stage 2A

- Chemistry (σ/π bonds, shells, spectroscopy of CBi⁻, etc.)  
- Electrons / orbitals / valence  
- Quantitative match to 12 u or real binding MeV  
- Kernel modifications  

---

## 5. Geometry and numerics (first runs)

- **Tool:** `bin/gen_qball_multi` (MAXBALLS=16)  
- **Phases:** all \(\delta=0\) (co-phase → contact fusion, v71)  
- **Separation:** nearest-neighbour \(D \approx 10\) (he2/li3)  
- **Grid (production):** N=192, L=36 (v70/v71 nuclear resolution, dx≈0.377)  
- **Grid (smoke):** N=128, L=28 if CPU-only  
- **N,L must match the seed in the `.cfg`** — the cfg load path does *not*
  auto-override grid from the SFA header (only `scp_sim file.sfa …` does).  
- **Physics cfg:** `complex_phi=1`, `complex_gauge=1`, `g_gauge=0.05`,
  `m_theta=1.6`, `eta=0`, absorbing BC, init=sfa  
- **Time:** T≥300 (v71); longer if still shedding  
- **Layouts:**  
  - **c6:** regular octahedron, centers on axes at distance \(D/\sqrt{2}\)  
  - **c12:** regular icosahedron, edge \(D\)  

Seed interference: co-phase tails inflate seed \(Q\) by ~7–17% at D=8–10
(v71). Report both seed \(Q\) and asymptotic parked \(Q\).

---

## 6. Order of work (this campaign)

1. **[done]** This mapping + Q-scale survey  
2. **[next]** Seed + run **c6_light** (stable Z-carbon)  
3. Seed + run **c12_light** (A=12 super-critical control)  
4. Analyse on-branch, mass defect, renders; write `v74/RESULTS.md`  
5. Decide: park Z-carbon as the Stage 2A deliverable **or** open g=0.02
   shooter for true A=12  
6. Stage 2B / 3 only after user priority shift  

---

## 7. Files

| path | role |
|------|------|
| `v74/CARBON_MAP.md` | this spec |
| `v74/profiles/f_w146_g005.txt` | light nucleon radial f(r) |
| `v74/profiles/f_w142_g005.txt` | standard nucleon radial f(r) |
| `v74/analysis/make_carbon_seeds.sh` | geometry + gen_qball_multi |
| `v74/cfg/*.cfg` | sim configs |
| `v74/seeds/*.sfa` | multi-ball seeds (local or /space) |

Reference prior art: `v71/NUCLEI.md`, `v69/theory/gscan.tsv`.
