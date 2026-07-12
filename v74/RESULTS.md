# v74 RESULTS — Stage 2A liquid-drop Z-carbon (`c6_light`)

**Date**: 2026-07-11  
**Run**: GPU V100, N=192, L=36, T=300, g=0.05, ω_seed=1.46, 6 co-phase light nucleons (octahedron, D=10)  
**Wall**: 38.1 min (71.8 ms/step)  
**Data**: `/space/scp/v74/results/c6_light.{sfa,diag.tsv,log}` (SFA ~35 GB f32 snaps)

## Verdict

**Z-carbon nucleus forms.** Six light free nucleons fuse into a **single coherent
charge droplet** that relaxes toward the gauged Q-ball branch. Endpoint is
**mid-high branch** (Q≈650 ≈ 0.71 Q_max), not super-critical evaporation to the
cap — despite a seed inflated near Q_max by co-phase interference.

Caveat (same class as v71 he2/li3): at T=300 the droplet is **still slowly
shedding** (dQ/dt ≈ −0.17 / t.u. for t>250). Parked Q is an upper bound; a
longer hold or cooler seed would pin it further.

## Numbers

| quantity | t=0 (seed) | t=300 (end) | notes |
|----------|------------|-------------|--------|
| Q | 907.84 | **649.83** | −28.4% radiated/absorbed |
| E | 1360.7 | 939.0 | radiation + sponge |
| E/Q | 1.4988 | **1.4450** | |
| Branch E/Q(Q) pred | — | **1.4377** | gscan g=0.05 interp |
| On-branch Δ | — | **+0.51%** | he2 was +0.1%; acceptable for still-hot |
| gauss_max | 9.8e-14 | 1.0e-13 | floor held entire run |
| ω_core | 1.460 | 1.365 | drops during merger |
| r_core | 7.61 | 6.94 | |
| Q_core | 591 | 592 | core inventory ~stable late |

**Free-nucleon baseline** (gscan light, ω=1.46): Q_N=114.13, E/Q=1.5187.  
**6× free inventory** = 684.8; seed Q = 907.8 → **+32.6% interference excess**
(stronger than v71’s 7–17% at A=2–3 — six overlapping co-phase tails).

**Mass defect** vs free light nucleon:  
\(1 - (E/Q)_\mathrm{final}/(E/Q)_N =\) **4.85%** of rest mass per unit charge  
(compare v71 li3s ~3.5%, he2 ~1.6%).

## Time history (selected)

| t | Q | E/Q | ω_core | r_core | Q_core |
|---|----|-----|--------|--------|--------|
| 0 | 907.8 | 1.4988 | 1.460 | 7.61 | 591 |
| 30 | 907.8 | 1.4985 | 1.415 | 6.01 | 801 |
| 50 | 907.5 | 1.4980 | 1.358 | 5.54 | 781 |
| 100 | 892.0 | 1.4913 | 1.372 | 12.55 | 607 |
| 150 | 733.0 | 1.4598 | 1.388 | 11.95 | 601 |
| 200 | 668.8 | 1.4482 | 1.334 | 8.34 | 606 |
| 250 | 658.3 | 1.4468 | 1.380 | 7.05 | 605 |
| 300 | 649.8 | 1.4450 | 1.365 | 6.94 | 592 |

Interpretation: contact/merge by t~30–50 (ω drops, Q_core peaks); violent
reconfiguration through t~100–150 (r_core balloons — multipole ringing /
radiation); then settle with gradual charge loss toward branch.

## Mapping outcome (Stage 0 → 2A)

| claim from `CARBON_MAP.md` | result |
|----------------------------|--------|
| Primary stable target = 6 light nucleons at g=0.05 | **Confirmed path works** |
| Seed near interference-inflated Q may sit near Q_max | Seed Q=908 ≈ 0.99 Q_max |
| Final should park mid-branch if subcritical after radiation | **Q→650, on-branch ~0.5%** |
| A=12 free always super-critical | Untested this run (c12_light still queued) |

**Structural name:** fabric **Z-carbon nucleus** (A=6 liquid-drop inventory
in formation history; identity is Q + branch point + fusion history, not retained
nucleon count — v71 doctrine).

## Visual verification — c6 (2026-07-12)

Slice panel `v74/results/c6_light.png` (z=0 plane, frames t≈0/15/30/50/100/200/300):

- **t=0**: four equatorial light nucleons of the octahedron (poles off this plane)
- **t=15–50**: co-phase infall → single central droplet with multipole radiation
- **t=100–300**: single compact charge blob + Coulomb `Emag` ring; no fission
- Component densities stay co-located (rho2_0/1/2 track); θ inert

Matches v71 li3s morphology at higher inventory. **Render OK → ran c12_light.**

---

# c12_light — A=12 supercritical control

**Date**: 2026-07-12  
**Run**: GPU V100, N=192, L=36, T=300, g=0.05, ω_seed=1.46, **12** co-phase light
nucleons (icosahedron, D=10)  
**Wall**: 38.5 min (72.5 ms/step)  
**Data**: `/space/scp/v74/results/c12_light.{sfa,diag.tsv,log}` (SFA ~42 GB)

## Verdict

**A=12 free inventory is super-critical and stays that way through T=300.**  
Twelve light nucleons fuse into **one** droplet (no fission), then the hot droplet
**evaporates charge** toward Q_max — same channel as v71 li3, at larger excess.

| | seed | t=300 |
|--|------|-------|
| Q | **1959.9** (2.13× Q_max) | **1410.5** (1.53× Q_max) |
| E/Q | 1.4968 | 1.4265 |
| Q_max | 921 | still **+489 above cap** |
| Late dQ/dt (t>250) | — | **−0.30 / t.u.** |
| Linear time-to-cap from t=300 | — | ~1600 t.u. more (rough) |
| gauss_max | 2e-13 | floor held |
| Mass defect vs free light N | — | 6.07% (still hot / off-branch above Q_max) |

Seed interference: 12×114.13 free = 1370 → seed Q=1960 (**+43%**).  
Violent merger around t≈50 (s_max peaks ~14, phi_max ~1.56) then single-blob settle
with continued charge loss.

**Confirms CARBON_MAP theorem:** at g=0.05, free A=12 cannot park as a static
branch isotope; the honest A=12 product is a **hot evaporating droplet**, while
stable parked carbon is the **Z=6 (c6_light)** path.

## Visual — c12 (`v74/results/c12_light.png`)

- **t=0–30**: multi-ball icosahedron (z-slice shows ~5–6 of 12 vertices)
- **t=50**: merger violence → one core + strong radiation
- **t=100–300**: single larger droplet + multi-ring radiation + Coulomb halo
- No binary fission; morphology is liquid-drop evaporation, not breakup

## Stage 2A summary

| run | A | seed Q | final Q | vs Q_max | fate |
|-----|---|--------|---------|----------|------|
| **c6_light** | 6 | 908 | **650** | 0.71× | **parks mid-branch** (Z-carbon) |
| **c12_light** | 12 | 1960 | **1411** | 1.53× | **still evaporating** (hot A=12) |

## Not yet done

- Longer T or cooler seed for fully parked dQ/dt≈0 on c6
- True A=12 at g=0.02 (raise Q_max)
- Stage 3 electrons

## Files

| path | role |
|------|------|
| `v74/CARBON_MAP.md` | Stage 0 spec |
| `v74/STATUS.md` | campaign status |
| `v74/cfg/c6_light.cfg`, `c12_light.cfg` | production configs |
| `v74/results/c6_light.png`, `c12_light.png` | slice panels |
| `/space/scp/v74/results/c6_light.*` | Z-carbon data (~35 GB SFA) |
| `/space/scp/v74/results/c12_light.*` | A=12 control (~42 GB SFA) |
