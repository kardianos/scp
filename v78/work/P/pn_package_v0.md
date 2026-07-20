# P/N Package v0 — Proton- and neutron-analogs (B2 multi-fabric)

**Agent:** P (Proton/neutron)  
**Date:** 2026-07-19  
**Campaign:** v78 Phase 2 · Checkpoint **CP-P**  
**Evidence base:** v75 F17–F19 · `v75/STATE.md` · `v75/PN_EXPERIMENT.md` · `v75/C12_ATOM_GOAL.md` §P/N  
**Status:** **CP-P ADOPT** (charge-structure package freeze; residual packaging gaps deferred to C/A)

---

## 0. Claim (one line)

Under B2 unlock with \(q_C=0,\,q_Q=+1,\,q_L=-1\):

\[
\boxed{\text{proton} = C+Q\ ({\rm co\text{-}located});\quad
\text{neutron} = C\text{-only};\quad
Q_{\mathrm{em}}\propto Z,\ N\ \text{adds bag mass only}}
\]

Flavored \(\Delta\omega\) alone is **not** the Z/N knob (F17b/c null).

---

## 1. Field assignment (frozen)

| Symbol | Fabric content | EM source \(\rho_{\mathrm{em}}\) |
|--------|----------------|----------------------------------|
| **p** (proton analog) | \(\Phi_C\) bag + \(\Phi_Q\) charge, co-located seed | \(q_Q\rho_Q \neq 0\) |
| **n** (neutron analog) | \(\Phi_C\) bag only (\(\Phi_Q=0\) at site) | \(\approx 0\) |
| **e** / light | \(\Phi_L\), \(q_L=-1\) | opposite nuclear \(Q\) |

**Config surface (mandatory for P/N):**

```text
n_fabrics=3
mf_lock_CQ=0          # B2 unlock — REQUIRED
mf_stage=2
q_C=0  q_Q=1  q_L=-1
complex_phi=1 complex_gauge=1 g_gauge=0.05 m_theta=1.6 eta=0
init=sfa
init_sfa=..._C.sfa
init_sfa_Q=..._Q.sfa   # protons only; empty/zero Q for pure-n
init_sfa_L=..._L.sfa   # zero (nuclear) or light shell (atom)
```

**Why B2:** B1 lock \(\Phi_Q\equiv\Phi_C\) forces every nuclear lump to source EM — no neutral species.  
**Why \(q_C=0\):** bag fabric does not couple to \(A\); only Q (and L) contribute to Gauss source.

Bookkeeping:

\[
Z=\#\{\text{p seeds}\},\quad N=\#\{\text{n seeds}\},\quad A=Z+N.
\]

Isotope control: vary \(N\) at fixed \(Z\) (and fixed L count \(=Z\)); nuclear \(Q_{\mathrm{em}}\) fixed, bag mass / \(Q_\phi\) change.

---

## 2. What is *not* P/N

| Approach | Result (F17) | Verdict |
|----------|--------------|---------|
| Flavored multiplet \(\Delta\omega\) alone (v71 style) | F17b: still fully charged \(Q_{\mathrm{flux}}\approx 284\) | **not n** |
| Equal-\(f\) cancel \(\omega=(+w,+w,-w)\) | F17c: residual \(Q\sim Q/3\) | **not n** |
| B1 lock only | every nuclear ball sources EM | **no neutral** |

---

## 3. Evidence (F17 primary; F18–F19 scale)

**Data:** `/space/scp/v75/pn/` (F17), `.../p2/` (F18), `.../z6/` (F19)  
**Plan:** `v75/PN_EXPERIMENT.md` · scores: `v75/results/pn_z6/` · scorer: `v75/analysis/score_pn_park.py`

### 3.1 F17 — single-species + small multi (N=96/128)

| ID | Relation | \(Q_\phi\) end | \(Q_{\mathrm{flux}}\) end | Gate |
|----|----------|----------------|---------------------------|------|
| F17d | **B2 neutron** C-only | 209.6 | **0.000** | **S1 PASS** |
| F17e | **B2 proton** C+Q | 209.6 | **209.2** | **S2 PASS** |
| F17f | Z1 N1 | 439.5 | 205.9 ≈1p | S3 |
| F17g | Z2 N0 | 439.5 | 422.7 ≈2p | S3 |
| F17h | Z2 N2 | 1128 | 387.5 (≈g; −8%) | S3≈ |
| F17a/b | flavored baselines | charged | large | **S5 null** |
| F17c | cancel multiplet | 69.9 | 69.8 residual | not n |

Scorecard F17:

| Gate | Criterion | Result |
|------|-----------|--------|
| **S1** | n: \(\lvert Q_{\mathrm{em}}\rvert/\lvert Q_C\rvert < 0.05\) | **PASS** (\(Q_{\mathrm{flux}}=0\)) |
| **S2** | p: \(Q_{\mathrm{em}}\approx Q_C\) | **PASS** (0.998) |
| **S3** | \(Q_{\mathrm{em}}(Z,N)\approx Q_{\mathrm{em}}(Z,0)\) | **PASS≈** (Z2: −8% flux) |
| **S4** | lumps survive T≳80 | **PASS** (Q ret ≥0.94) |
| **S5** | flavor ≠ Z/N | **PASS (null)** |

### 3.2 F18 — Z=2 park + L=Z + isotope (N=128, T=400)

| Metric | Z2N0 | Z2N2 |
|--------|------|------|
| \(Q_{\mathrm{flux}}\) end | **372** | **372** (identical) |
| \(c_{Q,\mathrm{park}}\) | 0.028 | 0.021 |
| PASS_nuc | True | True |
| Atom L=2 massL loss | −0.6% | **same** |

### 3.3 F19 — Z=6 carbon-class (N=192 L=48, T=400)

| Metric | Z6N0 | Z6N6 |
|--------|------|------|
| \(Q_{\mathrm{flux}}\) end | **990.3** | **990.3** (identical) |
| \(Q_{\mathrm{em}}/Z\) | 165 | 165 |
| \(c_{Q,\mathrm{park}}\) | 0.184 **soft** | 0.046 **PASS** |
| L=6 massL loss (atom) | −12.5% | **same** |

**Isotope EM holds through Z=6.** Neutrons add bag mass, not EM flux.

---

## 4. Gates for Z vs N (operational)

Use these for any new (Z,N) package. Scorer: `v75/analysis/score_pn_park.py`.

### 4.1 Species gates (unit p / unit n)

| ID | Gate | Pass bar |
|----|------|----------|
| **PN-S1** | Neutron EM-null | \(\lvert Q_{\mathrm{flux}}\rvert / \max(\lvert Q_\phi\rvert,1) < 0.05\); \(E_{\mathrm{em}}\approx 0\) |
| **PN-S2** | Proton charged | \(\lvert Q_{\mathrm{flux}}/Q_\phi - 1\rvert < 0.05\) (same-seed C+Q) |
| **PN-S4** | Lump existence | Q retention ≥ 0.9 over T≳80; no full dispersal |

### 4.2 Composition / isotope gates

| ID | Gate | Pass bar |
|----|------|----------|
| **PN-Z** | EM tracks Z | \(Q_{\mathrm{flux}} \approx Z\times Q_{\mathrm{em,1p}}\) within ~10% after park |
| **PN-N** | N is isotope knob | \(\lvert Q_{\mathrm{flux}}(Z,N)-Q_{\mathrm{flux}}(Z,0)\rvert / \max(\lvert Q_f\rvert,1) \le 0.10\) |
| **PN-A** | Bag tracks A | \(Q_\phi\) and \(E\) increase with N at fixed Z |
| **PN-L** | L count = Z not A | same L package for isotope pair; massL/Ql identical within noise |

### 4.3 Nuclear park (multi-ball)

\[
c_{Q,\mathrm{park}}=\frac{|Q_{\mathrm{mid}}-Q_{\mathrm{end}}|}{|Q_{\mathrm{mid}}|}\le 0.15,
\quad
c_{Q\mathrm{em},\mathrm{park}}\le 0.20,
\quad
\mathrm{gauss\_max}\sim 10^{-13}.
\]

Do **not** treat late **net** \(Q_{\mathrm{flux}}\) drop alone as L death when massL/Ql hold (net = nuclear EM − L).

### 4.4 Gate → campaign unlock

| Gate set | Unlocks |
|----------|---------|
| PN-S1–S4 + PN-Z/N | **CP-P** (this package) |
| + Z6-class park | C-phase isotope templates |
| + L=Z hold | A-phase Z-count cloud |

---

## 5. Recipes / configs

### 5.1 Seed tool

```text
bin/gen_pn_core N L profNuc omega out_C out_Q out_L \
  nZ  xz yz zz ... \
  nN  xn yn zn ... \
  nL  xl yl zl ... \
  [profL omegaL]
```

Source: `sfa/seed/gen_pn_core.c`  
- C gets all nuclear balls (Z+N)  
- Q gets **only** Z proton centers  
- L optional; hierarchy via trailing `profL omegaL` (e.g. ω=1.46)

### 5.2 Config pointers (`v75/cfg/pn/`)

| Config | Role | Gates |
|--------|------|-------|
| [`f17d_b2n.cfg`](../../v75/cfg/pn/f17d_b2n.cfg) | B2 neutron (C-only) | PN-S1 |
| [`f17e_b2p.cfg`](../../v75/cfg/pn/f17e_b2p.cfg) | B2 proton (C+Q) | PN-S2 |
| [`f17a_sym.cfg`](../../v75/cfg/pn/f17a_sym.cfg) | flavored baseline (null) | S5 |
| [`f17b_flav.cfg`](../../v75/cfg/pn/f17b_flav.cfg) | flavored p-branch (null) | S5 |
| [`f17c_cancel.cfg`](../../v75/cfg/pn/f17c_cancel.cfg) | cancel multiplet (null) | S5 |
| [`f17f_z1n1.cfg`](../../v75/cfg/pn/f17f_z1n1.cfg) | Z1 N1 | PN-Z/N |
| [`f17g_z2n0.cfg`](../../v75/cfg/pn/f17g_z2n0.cfg) | Z2 N0 | PN-Z |
| [`f17h_z2n2.cfg`](../../v75/cfg/pn/f17h_z2n2.cfg) | Z2 N2 isotope | PN-N |
| [`p2_z2n0_nuc.cfg`](../../v75/cfg/pn/p2_z2n0_nuc.cfg) | Z2N0 park T=400 | park |
| [`p2_z2n2_nuc.cfg`](../../v75/cfg/pn/p2_z2n2_nuc.cfg) | Z2N2 park | isotope |
| [`p2_a_z2n0.cfg`](../../v75/cfg/pn/p2_a_z2n0.cfg) | Z2 + L=2 atom | PN-L |
| [`p2_a_z2n2.cfg`](../../v75/cfg/pn/p2_a_z2n2.cfg) | Z2N2 + L=2 | PN-L isotope |
| [`z6_n0_nuc.cfg`](../../v75/cfg/pn/z6_n0_nuc.cfg) | Z6N0 nuclear | C-class |
| [`z6_n6_nuc.cfg`](../../v75/cfg/pn/z6_n6_nuc.cfg) | Z6N6 nuclear | C-class isotope |
| [`z6_a_n0.cfg`](../../v75/cfg/pn/z6_a_n0.cfg) | Z6 + L=6 | atom |
| [`z6_a_n6.cfg`](../../v75/cfg/pn/z6_a_n6.cfg) | Z6N6 + L=6 | atom isotope |
| [`run_z6_remote.sh`](../../v75/cfg/pn/run_z6_remote.sh) | F19 campaign launcher | ops |

### 5.3 Grid freeze (boundary safety)

See `v75/PN_GRID.md`.

| Package | N | L | damp_width | Notes |
|---------|---|---|------------|-------|
| F17 single | 96 | 18 | 3 | unit p/n |
| F18 Z2 | 128 | 32 | 4 | L shell R=18 |
| F19 Z6 | **192** | **48** | **5** | p octa R=8; n R=5.5; L R=22 |

**Do not** put Z6 L@R=22 in F18’s L=32 box (sponge eats L).

### 5.4 Standard physics constants

```text
m^2=2.25  mu=-41.345  kappa=50
m_theta=1.6  eta=0  g_gauge=0.05
absorbing BC  complex_phi=1  complex_gauge=1
```

---

## 6. Diagnostics map

| Proxy | Meaning | Use |
|-------|---------|-----|
| `Q_phi` | C-bag Noether (primary arrays) | bag / A inventory |
| `Q_flux` | Gauss-cube EM \(\propto\int(\nabla\cdot E)/g\) | nuclear \(Q_{\mathrm{em}}\) proxy |
| `Q_C,Q_Q,Q_L,Q_em` | CPU multi-fab integrals | exact fabric charges (when enabled) |
| track massL / Ql | L survival | atom; **not** net flux alone |
| `gauss_max` | discrete Gauss residual | must stay ~1e-13 |

GPU campaigns: use \(Q_{\mathrm{flux}}\) + \(Q_\phi\) as EM/bag proxies (CPU smoke confirmed exact fabric charges for F17).

---

## 7. Gaps (honest; do **not** block CP-P)

| Gap | Severity | Owner phase |
|-----|----------|-------------|
| Multi-ball cores **fuse to single droplet** by t∼400 | packaging / morphology | **C** (nucleus) |
| Z6N0 park soft (\(c_Q=0.184>0.15\)); Z6N6 clean | park retune | **C** |
| L mass loss −12.5% at Z=6 (vs −0.6% at Z=2) | L hold | **L/A** |
| Shell-radius diagnostic (COM D≈0 insufficient) | atom structure | **A** |
| Q self-bag under B2 still on (bag on every fabric) | optional cleanup | kernel later |
| Long-T visual C₁₂ not claimed | product goal | **A/C** |

These are **not** failures of the p vs n *definition*. Charge bookkeeping (S1–S3, PN-Z/N) is closed.

---

## 8. STAMP CP-P

```text
STAMP CP-P: ADOPT
```

**Rationale:** P2.0 charge structure closed at machine-measured level (F17 S1–S5); isotope EM control demonstrated at Z=2 (F18) and Z=6 (F19); recipes + configs + gates packaged. Residual items are C/L/A packaging, not P/N species definition.

**Co-agree:** U must co-stamp on board (`CAMPAIGN_MAP.md`: P+U).  
**Unlocks:** C isotopes (N at fixed Z); A Z-count for L cloud.

**Reject triggers (if reopened):** new data showing \(Q_{\mathrm{em}}\) tracks A not Z under B2; or neutron acquiring finite \(Q_{\mathrm{flux}}\) under \(q_C=0\) dynamics.

---

## 9. Handoff

| To | Action |
|----|--------|
| **U** | Co-stamp CP-P ADOPT; enter package into PARTICLE_LADDER / RECIPES |
| **C** | Use `gen_pn_core` Z/N templates; isotope = +N at fixed Z; retune Z6N0 park |
| **L** | L count always **Z**, never A; F19 L hold is soft at carbon-class |
| **A** | Assemble (Z≈6, N≈6) + L=6 from this package; do not invent alternate P/N |
| **N** | Unit nucleon template feeds profile into `gen_pn_core` (ω, prof) |

**Canonical physics snapshot:** `v75/STATE.md` §3.  
**This file:** v78 freeze of that package for campaign board.
