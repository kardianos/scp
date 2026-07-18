# Operational Rods and Clocks from Free Medium (C1) — v0

**Approach:** C — reverse theoretical  
**Date:** 2026-07-18 (Round 2)  
**Status:** operational sketch (not closed dynamics)  
**Depends on:** NC-1.1–1.3; Ax3, Ax7 (A); PROBLEM.md eligibility  
**Goal:** Constraints on media so “local \(c=\mathrm{const}\)” is medium-made, not external stage instruments.

---

## 0. The circularity problem

If rods and clocks are classical objects *outside* the continuum, then:

- “local \(c=\mathrm{const}\)” is a statement about matter tools on a stage;
- free-medium warp becomes secondary (dualist measurement).

Monism requires: **rods and clocks are free-medium states** (or weakly bound oscillators of the same continuum).

Circularity risk: defining geometry from rods that already assume geometry.

**Escape pattern (operational, Einstein 1905-style but medium-first):**  
Use free-signal *procedures* to *define* simultaneity, length, and time so that free locality is isotropic at value \(c\) **by construction** in each small free patch; then nontrivial global charts are discoveries, not inputs.

---

## 1. Primitive free operations (no chart)

Assume the medium admits:

1. **Free disturbances** — localized free updates that can influence neighbors (Ax3).  
2. **Repeatable free oscillators** — configurations that tick (phase cycles) when free budget is available (clock prototypes).  
3. **Echo / round-trip** — ability to send a free disturbance along a free path and detect return (if the medium allows reflection or a free partner emitter).

No a priori \(\mathbb{R}^{1,3}\) metric.

---

## 2. Local clock (operational)

### 2.1 Definition

A **local free clock** at event-region \(U\) is a free oscillator \(O_U\) whose cycle is a free-medium process (e.g. free-field phase, free standing mode in a small cavity of free medium).

**Tick:** one period \(\tau_U\) of \(O_U\).

### 2.2 Necessary conditions

| ID | Condition |
|----|-----------|
| **RC-C1** | Clock is made of free (or weakly bound) medium; removing free budget changes \(\tau_U\) or kills the oscillator |
| **RC-C2** | In a small free region with no nearby locks, two free clocks of the same construction agree after transport (free isotropy) |
| **RC-C3** | Clock rate is a functional of local free state; not of an independent metric field |

### 2.3 Near locks

If free budget or free path structure changes near a lock, \(\tau_U\) may change when compared via free-signal links to a distant clock.  
**That comparison *is* gravitational redshift bookkeeping** under monism — not a second law — provided the same free medium defines both clocks and links.

---

## 3. Local rod (operational)

### 3.1 Einstein-light rod (medium version)

In a small free patch \(U\):

1. Place free clocks \(A,B\) of identical construction at the ends of a prospective rod.  
2. Send free signals \(A\to B\to A\).  
3. **Define** the rod to be “at rest and of unit length” when:
   - round-trip free time is \(2/c\) in units where the free oscillator defines the second, **and**
   - one-way times match under Einstein synchrony built from free signals,
   - free signal speed is isotropic in the patch.

**Local \(c\):**  
\[
c
\;\equiv\;
\frac{2L_{\mathrm{rod}}}{\tau_{\mathrm{round}}}
\]
is **normalized to a constant** by the rod/clock construction in each free local frame.  
That is Ax7 / NC-1.2 as reverse constraint: media must **admit** such local calibration.

### 3.2 Necessary conditions

| ID | Condition |
|----|-----------|
| **RC-R1** | Rod length is free-signal calibrated, not a fixed lattice edge length as ontology |
| **RC-R2** | Isotropy: rotating the rod in a free patch leaves \(c\) unchanged |
| **RC-R3** | Rod material is free or weakly bound medium; lock formation nearby can change operational length vs distant standard |

### 3.3 Lattice warning (B sandboxes)

Fixed \(\Delta x\) on a numpy grid is **storage**.  
True monist rod length must eventually be **path cost in free medium** (B1 graph) or free-signal calibrated optical length.  
Round-1 B chart rods are dualist residue — allowed only if tagged scaffolding.

---

## 4. Local frame and “\(c=\mathrm{const}\)”

### 4.1 Construction sketch

1. Choose a free event \(p\) with free budget above threshold.  
2. Build free clock \(O_p\).  
3. Emit free signals in several directions; adjust operational rod units so one-way / round-trip free speeds match (isotropy).  
4. Declare \(c_p\) the common free speed in these units.

**Axiom-level reverse demand:** such frames exist wherever free medium exists (NC-1.2).

### 4.2 What is forced on the medium

| Demand | Medium implication |
|--------|-------------------|
| Local isotropy of free speed achievable | Free update law locally isotropic in free state (or frame can absorb anisotropy) |
| \(c\) independent of source motion (usual relativity) | Free update law has no preferred free-emitter rest frame (local Lorentzian or analog) |
| Frames exist near weak locks | Free medium not fully extinguished outside lock core; local free oscillators still work |
| Global charts nontrivial | After local calibration everywhere, free transit times between distant free patches depend on intervening locks (warp) |

---

## 5. Transport and comparison (warp without second sector)

### 5.1 Clock comparison

Two free clocks \(O_1,O_2\) compared by free-signal link:
\[
\frac{\tau_1}{\tau_2}
\;=\;
\text{functional of free path structure between them}.
\]
Under monism this ratio is **not** “time dilation from metric” as external cause; it **is** free-medium path + local free rate structure.

### 5.2 Necessary equality (weak field target)

For consistency with C3 path-cost monopole and constant local \(c\):
\[
\frac{\delta\tau}{\tau}
\sim
\mathcal{O}\!\left(\frac{G_{\mathrm{eff}} M}{c^2 r}\right)
\]
when compared to asymptotic free clocks — same class as \(\delta\ell\) (path_cost_profile §5).

### 5.3 Fail modes

| Fail | Meaning |
|------|---------|
| Clocks immune to free budget | Dualist instruments (S16) |
| Only bound matter clocks, free medium timeless | Split ontology |
| Local \(c\) depends on direction after calibration | Medium lacks free isotropy / frame construction fails |
| Local \(c\) changes with free packet energy arbitrarily | No stable locality bound |

---

## 6. Relation to mass and inertia (C2)

Rods/clocks fix the operational meaning of \(c\) and of energy units (work counted in free-clock ticks × free-rod forces).

Then \(m=E_\star/c^2\) uses **that** \(c\).  
If clocks were external, \(c^2\) would smuggle stage structure into mass.

**NC-link:**  
Unlock work and inertial response must be measured with free-medium rods/clocks (or prove equivalence). Same ledger, same \(c\).

---

## 7. Minimal numeric test (FOR_B)

Without full relativity machinery:

1. Place free “oscillator” proxies: period \(\propto 1/\sqrt{\rho_f}\) or constant free tick if \(\rho_f>\theta\).  
2. Round-trip free signal time between two grid points using free-path cost (\(n\) or kernel \(\Phi\)).  
3. Define local \(c_{\mathrm{op}}=2L_{\mathrm{op}}/\tau_{\mathrm{rt}}\) with \(L_{\mathrm{op}}\) from free calibration.  
4. **Pass:** \(c_{\mathrm{op}}\) constant in free patches after calibration; **varies only as chart artifact** before calibration.  
5. **Near lock:** round-trip excess matches path-cost integral; oscillator rate shift tracks free state if model includes it.

---

## 8. What C1 does *not* solve yet

- Unique free oscillator microphysics.  
- Full Lorentz group from free updates.  
- Machian global inertia from distant free medium.  
- Quantum clocks.

Enough for Round 2: **constraints + procedure schema** so A Ax7 is reverse-backed and B cannot treat fixed \(\Delta x\) as true length forever.

---

## 9. Bottom line

> **Rods** = free-signal calibrated segments of free medium.  
> **Clocks** = free oscillators of free medium.  
> **Local \(c=\mathrm{const}\)** = normalization of free signalling in those units.  
> **Warp** = failure of a single global chart to keep free transit trivial once locks rearrange free medium.

Necessary: RC-C1–C3, RC-R1–R3; lattice edges ≠ ontology; dualist external instruments stripped (S16).
