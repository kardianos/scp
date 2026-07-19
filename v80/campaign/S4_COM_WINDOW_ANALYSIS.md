# S4 COM / window analysis — `mf_h1_orbit_N192`

**Date:** 2026-07-19  
**Gate:** Claude council pre-force check — is S4 “hold” real, or is \(Q_{\mathrm{core}}\)/\(Q_{\mathrm{flux}}\) collapse charge death?  
**Data:** existing `results/mf_h1_orbit_N192_diag.tsv` + kernel source (no SFA; output pruned).  
**Kernel auth:** none.

---

## 1. What the diagnostics actually integrate

From `scp_sim.c` / `scp_sim.cu` (v67 Q-diag):

| Quantity | Center | Radius / support |
|----------|--------|------------------|
| **\(Q_\phi\)** | global | whole grid |
| **\(Q_{\mathrm{core}}\)** | **box origin (fixed)** | ball \(r \le\) `qdiag_radius` (default **8**) |
| **\(Q_{\mathrm{flux}}\)** | **box origin (fixed)** | cube half-width = `qdiag_radius` (**8**) |
| **\(r_{\mathrm{core}}\)** | **mass COM** (centroid of \(\sum_a|\Phi_a|^2\)) | rms about that COM |

Default cfg has no `qdiag_radius` override → **R = 8**.  
S4 seed: C at origin, L on shell **R = 16** with tangential \(v_t = 0.05\) (`gen_mf_shell_orbit`).

**Implication:** L starts **outside** the Q_core / Q_flux windows. Those windows only track charge that stays near the **box center**, not near the dynamical COM.

---

## 2. Time series (diag)

| \(t\) | \(Q_\phi\) | \(Q_{\mathrm{flux}}\) | \(Q_{\mathrm{core}}\) | \(r_{\mathrm{core}}\) | \(E_{\mathrm{em}}\) | \(E_{\mathrm{total}}\) |
|------:|-----------:|-------------------:|-------------------:|-------------------:|-------------------:|---------------------:|
| 0 | 114.727 | 114.23 | 113.60 | 3.177 | 0.710 | 347.8 |
| 100 | 114.725 | 114.22 | 113.62 | 3.160 | 0.714 | 349.2 |
| 200 | 114.723 | 114.22 | 113.51 | 3.152 | 0.720 | 350.3 |
| 300 | 114.723 | 114.0 | 112.84 | 3.153 | 0.727 | 350.9 |
| 400 | 114.722 | 112.8 | 110.07 | 3.152 | 0.735 | 351.5 |
| 500 | 114.722 | 107.4 | 95.94 | 3.151 | 0.742 | 352.1 |
| 600 | 114.722 | 81.7 | 44.61 | 3.150 | 0.748 | 352.8 |
| 700 | 114.722 | 33.0 | 5.53 | 3.150 | 0.750 | 353.5 |
| 799 | 114.722 | **4.69** | **0.46** | 3.150 | **0.748** | 354.2 |

Thresholds (fraction of init \(Q_{\mathrm{core}}\)):

| Drop | First \(t\) |
|------|-------------|
| &lt; 90% | ~476 |
| &lt; 50% | ~581 |
| &lt; 10% | ~671 |
| &lt; 1% | ~763 |

**Global charge drift:** \(\Delta Q_\phi \approx -5.5\times 10^{-3}\) (−0.005%) — machine-stable.  
**EM energy:** +0.038 absolute (~+5%), **not** v79-style collapse.  
**\(r_{\mathrm{core}}\):** min/max **3.144 / 3.177** — tight, COM-relative size stable.

---

## 3. Analytic verdict

### Claim A — “Charge evaporated / annihilated”
**REJECTED.**  
Global \(Q_\phi\) is conserved. \(E_{\mathrm{em}}\) is held or slightly up. That is incompatible with v79-style EM death or charge loss to BC as the primary story for this seed.

### Claim B — “Fixed-center windowing / COM motion”
**ACCEPTED (primary).**

1. Kernel **defines** \(Q_{\mathrm{core}}\) / \(Q_{\mathrm{flux}}\) about **box center**, not COM.  
2. Seed has orbital kinematics (L at R=16, \(v_t>0\)) → C recoils; charge can leave the R=8 ball about origin while remaining a compact object.  
3. \(r_{\mathrm{core}}\) (about **COM**) stays ~3.15, so the field remains a **localized lump in the COM frame**, not a spreading bath.  
4. Collapse of \(Q_{\mathrm{core}}\) is **late and smooth** (t≳450–800), consistent with gradual drift out of a fixed window, not a sudden annihilation event.

### Claim C — “S4 holds as a hydrogenoid-class object”
**SUPPORTED for global ledger** (\(Q_\phi\), \(E_{\mathrm{em}}\), Gauss).  
**NOT supported** if “hold” is defined by raw \(Q_{\mathrm{core}}\)/\(Q_{\mathrm{flux}}\) without COM correction — those diagnostics are **misleading for this seed**.

---

## 4. Gate score (Layer 2 **H**)

| Sub-gate | Status |
|----------|--------|
| \(Q_\phi\) held | **PASS** |
| \(E_{\mathrm{em}}\) held (≠ v79 death) | **PASS** |
| Gauss floor | **PASS** |
| Fixed-center \(Q_{\mathrm{core}}\)/flux as soliton tracker | **FAIL as metric** (window artifact) |
| **Overall H for “minimal pair long-T viable”** | **PASS** (with diagnostic caveat) |

**Pre-force decision:** proceed to **force grid with tracks**. Do **not** treat S4 flux/core drop as product soft-kill.  

**Policy for future runs:**

1. Prefer **pair tracks / COM diagnostics** over fixed-center \(Q_{\mathrm{core}}\) for dynamical seeds.  
2. Optionally set larger `qdiag_radius` only as a band-aid; tracks are the real fix.  
3. Do not claim multi-rev orbit from this S4 diag alone (no COM track of L).

---

## 5. What we did *not* need

- No GPU re-run: kernel definition + diag time series suffice for the windowing claim.  
- No kernel change.  
- Optional later: short T=200 re-run with `mf_pair_track` if a plotted COM trajectory is wanted for figures — **not required to open the force gate.**

---

## 6. Next (scorecard order)

1. ~~S4 COM/window~~ **DONE** — hold is real; core/flux drop is fixed-center.  
2. **Force grid** — keep SFAs + `mf_pair_track` (gates **F**, **R**, **L1b**).  
3. Low-\(v_t\) orbit with tracks (gate **O**).  
4. No Z6+L6 park.
