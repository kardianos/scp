# Shear genesis — simulation design and instrumentation

**Status:** design only (pre-registration). No kernel commitment until
an implementation face is authorized. Ontology: free-cell field sector
(fa1, fa2, pass F) + conversion door; shear ≠ temperature; target =
**genesis**, not stability.

**Success including the “incomplete theory” case:** Em appears while
organized field shear is on, spatial correlation with shear, and Em
**dissipates after shear is removed** (or after a free shear packet
leaves). That is a **positive genesis result**. Stability is a later
arc.

**Kill control (honesty bar):** equal ⟨Ee⟩, **zero organized shear**
(phase-scrambled / isotropic) must birth **strictly less** Em (or
spatially uncorrelated Em). If S0 ≈ S1 at fixed energy, shear is not
generative.

---

## 0. What we are measuring

| is | is not |
|---|---|
| Does **organized field motion** (phase wind, counterflow, shear layer, sustained stream) **create dense Em** via the existing door? | Is Em stable forever? |
| Does birth track **shear meters** (∇φ, J_field, |∇·J|), not only |ψ|²? | Is shear “temperature”? |
| After drive ends, does Em die on a measurable half-life? | Have we found a proton? |

Macro “conveyor” language is only a ladder. At this scale the object is
**structured ψ_field transport**. Gravitational-lensing-like perception
at human scales is a *later* reading of space/field gradients — not a
meter for this face.

---

## 1. Physics sketch (for design, not law)

Field state per cell: `ψ = fa1 + i·fa2`, `Ee = |ψ|²`.

Link field current (same bilinear family as dense J):

```
J_ij ∝ 2 Im( ψ_i* · ψ_j )     (or hop-angle weighted form already in pass F)
```

Shear / organization proxies (per cell or per band):

```
gφ   = |∇φ|           phase gradient (finite-diff on neighbors)
|J|  = mean |J_ij| over contacts
Σ    = |∇·J| or interface indicator (see §4)
Ω    = |ω_field − ω_local|   beat readiness vs dense clocks if Em>0
```

Door already converts on beats + workfn + atoms. Hypothesis:

```
rate_cond = f( Ee, workfn, beats, Σ )
```

with `∂f/∂Σ > 0` in some window — **shear assists genesis**.

Temperature remains an attribute of condensed harmonic stacks **after**
Em exists. Shear is drive on the field book only.

---

## 2. Forcing methods (multiple “conveyors”)

All methods act on **field amplitudes**, not on dragging cell identities.
Prefer zero new knobs first; add knobs only when a method cannot be
seeded with existing apparatus.

### Family A — impulse / free packet (no continuous force)

Uses existing seeds; shear is **transient** as the packet crosses.

| id | method | how (existing or near) | shear character |
|---|---|---|---|
| **A1** | Tilted field pulse | `exp=pulse amp=… sigma=… kx=…` | bulk phase wind; free-propagating |
| **A2** | Spinning / radial tilt | pulse + `spin_m`, `imp_k` | azimuthal + radial phase structure |
| **A3** | Counter-propagating pair | two pulses (or `blob2`-style dual seed on **field**, not Em) with ±kx | **shear layer** on midplane |
| **A4** | Slit beam | `exp=slit` with kx, optional mask | directed beam; shear at edges / after slits |

**Pros:** no continuous apparatus energy; conservation clean; “shear then
gone” maps exactly to the positive transient-genesis criterion.  
**Cons:** shear duration short; duty cycle limited by box crossing time.

### Family B — boundary-driven sustained wind

New apparatus (small): each step (or each beat window), reinject field on
an inlet slab and absorb/open on an outlet slab.

| id | method | pin | shear character |
|---|---|---|---|
| **B1** | Uniform stream | inlet: fixed `|ψ|`, `φ = −k x + ωt` | bulk wind, weak shear (only at walls / structure) |
| **B2** | Couette-like | inlet φ-wind; transverse gradient of k or |ψ| in y | **∇v** across volume |
| **B3** | Shear sandwich | top half +k, bottom half −k (or two inlet bands) | interface shear layer through midplane |
| **B4** | Phase-locked vs random | same |ψ| spectrum; φ organized vs scrambled | the S0/S1 kill pair under continuous drive |

**Pros:** steady shear for long T; clean on/off protocols; map Em vs
shear strength.  
**Cons:** must **book apparatus δE** every reinjection; open system.

### Family C — volumetric drive (no preferred wall)

| id | method | notes |
|---|---|---|
| **C1** | Soft phase-wind force | weak continuous `dφ += α·dt` with spatial mask (e.g. only in a band) — careful: is a gauge drive; keep energy consistent by adjusting |ψ| only via door, not by hand |
| **C2** | Dual-mode beat bath | seed two plane waves (k1, k2); interference fringes = spatial beat map for door |
| **C3** | noise_amp control | existing random light — **this is S0 material**, not shear |

**C1 is last resort** (easy to invent fake energy). Prefer A and B.

### Recommended order of attack

```
1. A1 vs S0 (scrambled equal-Ee pulse)     — cheapest genesis test
2. A3 shear layer vs S0                    — interface nucleation?
3. B3 sustained shear sandwich + off-ramp  — transient after shear dies
4. B1/B4 continuous I(u,k) map             — only if genesis is real
```

Stability runs (R1–R8 style) are **out of scope** for this face unless
genesis is confirmed and someone later asks “does any residue hold.”

---

## 3. Protocol designs

### 3.1 Core suite SHEAR-G (genesis)

Common law surface (start simple; do not mix retention fights):

```
# genesis-first surface (proposal)
wf_on=1 wf_floor=0.01 wf_far=99     # door can open on presence / far demand
k_rad=0 or 0.05                     # two arms: pure birth vs birth+tax
amp_nat=0                           # dense transport additive default first
                                    # (optional amp_nat>0 later — not required for genesis)
bath=1 L=16 or 32  dt=0.02
seed fixed ensemble {111, 20260802, 314159}
```

**Why wf_on=1:** Face B showed dark bath ⇒ cond=0. Genesis needs a door
that can fire when field is present.  
**Why try k_rad=0 first:** tax confounds “born then died because shear
ended” with “born then died because radiance.” Add V3a tax in a second
pass.

#### Arm table

| arm | field seed | organized shear? | purpose |
|---|---|---|---|
| **S0** | pulse/stream, **random φ**, ⟨Ee⟩ matched to S1 | no | kill control |
| **S1** | pulse `kx=k`, same ⟨Ee⟩ | bulk wind | bulk genesis? |
| **S2** | dual ±kx packets / sandwich | **shear layer** | interface genesis? |
| **S3** | S1/S2, **drive off at t=t_off** | on then off | **transient success test** |
| **S4** | S0 with drive off at same t_off | no | off-ramp control |
| **G0** | no field (dark) | no | floor: cond must stay ~0 |

**Matching ⟨Ee⟩:** integrate `Σ Ee` at t=0+ (and for B-family, time-average
inlet power). S0 must not be “weaker light.”

**Positive result (any of):**

1. `Em_peak(S1) > Em_peak(S0)` and/or `Em_peak(S2) > Em_peak(S0)` by a
   pre-registered factor (e.g. ≥2×) at matched ⟨Ee⟩, seed-robust (≥2/3
   seeds).
2. Spatial: new Em correlates with high `gφ` or shear-interface mask
   more than with |ψ|² alone (partial correlation test §4.3).
3. **S3:** Em rises during shear window; after t_off, Em falls with
   measurable τ_half; S4 does not show a shear-timed Em bump.

(3) alone with weak (1) is still **very positive** if the Em bump is
shear-timed and S0-absent — it means genesis tracks the drive even if
net yield is modest.

**Negative result:** S0≈S1≈S2 for Em birth at all matched energies →
shear is not a genesis parameter; only energy density / random door
noise matters. Record and stop; do not “tune until shear wins.”

### 3.2 Off-ramp protocol (the incomplete-theory success)

```
t = 0        seed field shear (A or B)
t ∈ [0, t1]  shear on  — record Em(t), cond rate, shear meters
t = t1       remove organized drive:
               A-family: packet exits / absorb
               B-family: set inlet |ψ|→0 or switch inlet to S0 scramble
t ∈ [t1, t2] free evolution — fit Em(t) ~ Em(t1)·e^{−(t−t1)/τ}
```

**Pass (genesis, incomplete stability):** Em(t1) significantly above S0
floor, and τ finite (Em returns toward floor).  
**Do not require:** τ→∞, monodisperse hoards, phase retention, mobility.

### 3.3 Optional structure-in-flow (secondary)

Only after pure genesis (empty of Em seed) is measured:

- Place a mid-box dense seed (small blob, ring) **or** a Face-C hoard.
- Ask: does shear **increase I** onto existing Em (pickup), not only
  birth in vacuum?
- Meters: ΔEm on tagged support vs new Em outside tag.

This is “structure drinks the wind,” not the primary genesis claim.

---

## 4. Instrumentation (detectors)

Design principle: **separate genesis detectors from stability detectors.**
Reuse existing dumps; add offline analyzers; add only cheap online
RESULT lines if needed.

### 4.1 Existing streams to keep

| source | provides |
|---|---|
| log `RESULT drift_rel` | conservation; open B-family must report apparatus δE too |
| log `cond` / conversion totals (`cond_total`, rough, evap) | door traffic |
| `fcs` / `fcsdump -mode cells` | per-cell x,y,z, Em, Ee, fa1, fa2, th2, tag |
| `analyze_melt.py` | Em_tot, Em_max, PR, rms, comx vs t |
| `analyze_hoard.py` | hoard spectrum if clumps form |
| `p1_meter=1` | field/dense currents if available on arm |

### 4.2 New offline analyzer: `analyze_shear_genesis.py`

**Input:** `fcsdump -mode cells` time series (same as melt).

**Per frame, compute:**

| symbol | definition | role |
|---|---|---|
| `Ee_tot`, `Em_tot` | sums | books |
| `Em_max`, `PR`, `rms` | as melt | clump vs fog |
| `cond_proxy` | cells with Em>ε_em that had Ee>ε_ee recently — or use log cond | birth activity |
| `gφ_mean` | mean |∇φ| on contacts with Ee_i Ee_j > floor | organization |
| `Jf_mean` | mean |2 Im(ψ_i* ψ_j)| on those contacts | field current |
| `shear_iface` | for ±k arms: Em and Ee in midplane slab \|y−L/2\|<w | interface nucleation |
| `C(Em, gφ)` | spatial correlation / binned Em vs gφ | **key genesis detector** |
| `C(Em, Ee)` | same vs energy density | null: pure energy |

**Finite-difference phase gradient on free cells:** for each contact
(i,j) with d < contact rule,

```
dφ_ij = atan2( sin(φj-φi), cos(φj-φi) )
gφ_i  = mean |dφ_ij| / d_ij
```

φ = atan2(fa2, fa1). Empty field cells skipped.

**Output rows:**

```
t  Em_tot Em_max PR rms  Ee_tot gphi_mean Jf_mean  corr_Em_gphi corr_Em_Ee  Em_iface
```

**RESULT summary line (script end):**

```
# RESULT shear_genesis
#   Em_peak=… Em_end=… ratio_S_over_S0=…
#   corr_Em_gphi=… corr_Em_Ee=…
#   tau_half=… (if off-ramp)  pass_genesis=0|1
```

### 4.3 Correlation test (pre-registered)

Bin cells by gφ (or by distance to shear interface). Within each bin,
mean Em. Fit:

```
Em ~ a · Ee + b · gφ + c
```

(or nonparametric: partial correlation Em⊥Ee with gφ).

**Genesis-via-shear signature:** `b` significant and `corr(Em,gφ | Ee) > 0`
on S1/S2; consistent with ~0 on S0.

### 4.4 Online / log meters (cheap, optional kernel)

If implementing a face, add **byte-inert** diag lines (default off):

```
# SHEAR t=… Ee=… Em=… gphi=… Jf=… cond_rate=… app_dE=…
```

`app_dE` = cumulative energy injected by B-family boundaries (0 for A).

Conservation bar for B-family:

```
| (Etot − E0 − app_dE) / E0 | ≤ 1e-12-class
```

### 4.5 Spatial detectors (maps, not only scalars)

Offline 2D histograms (project z or use mid-slab):

| map | shows |
|---|---|
| `Ee(x,y)` | where field is |
| `Em(x,y)` | where matter appeared |
| `gφ(x,y)` | where shear/organization is |
| `Em − â·Ee` residual | matter not explained by energy alone |

For S2/B3, expect Em residual along the **shear interface** if the
hypothesis holds; for S0, residual noise-like.

Snapshot cadence: enough frames during packet transit (A: every
0.5–2 t.u.) and around t_off (S3).

### 4.6 Time detectors (the incomplete-theory pass)

| window | metric |
|---|---|
| shear on | dEm/dt, cond_rate, Em_peak |
| at t_off | Em(t1), gφ(t1)→0 or Jf drop |
| after | τ_half of Em; residual Em_end |

**Pass band (example, sharpen after pilot):**

```
Em(t1) / Em_S0(t1) ≥ 2
τ_half < (t2−t1)      # dissipates within the watch window
Em_end / Em(t1) ≤ 0.5
```

No requirement that Em_end = 0 exactly (numerical floor / tax residue).

### 4.7 What we deliberately do **not** meter as primary

| skip as primary | why |
|---|---|
| winding R2d / W | stability/phase retention (Face B/C territory) |
| blob_drift tagged centroid | misleading (item 2) |
| “temperature” | shear ≠ T; no T meter until condensed H defined |
| galactic lensing analogs | wrong scale; space book is another face |

Secondary optional: if Em clumps, run `analyze_hoard.py` for curiosity
only — not a pass bar.

---

## 5. Implementation path (when authorized)

### Phase 0 — zero kernel change (do first)

1. Build A1/A2 arms with `exp=pulse` + kx / spin_m / imp_k.
2. Build S0 by **phase scramble**: same |fa| pattern, replace φ with
   RNG uniform [0,2π) at t=0 (one-shot script, or small seed hook).
3. Match ⟨Ee⟩ post-seed.
4. Ship `report/analyze_shear_genesis.py` + shell `runs/shear_g0.sh`.
5. Run ensemble; write `SHEAR_G0_RESULTS.md`.

If Phase 0 already shows S1>S0 birth or shear-timed S3 transient →
**genesis evidence without new forcing code.**

### Phase 1 — dual packet / shear layer seed

- Extend pulse seed or add `exp=pulse2` (mirror blob2 but field-only):
  two Gaussians, ±kx, separation along y or x.
- Interface meters in analyzer (midplane slab).

### Phase 2 — sustained B-family (only if Phase 0–1 warrant)

- Byte-inert knobs e.g. `field_wind=0`, `wind_k`, `wind_amp`,
  `wind_ygrad`, `wind_app` (book δE).
- Inlet slab x < x_in: set ψ each step; outlet x > x_out: damp Ee.
- Full battery at knobs=0 must stay ALL GREEN 87.

### Phase 3 — not in this face

Metabolic stability, monodispersity, amp_door coherent refill, pass S,
DM/lensing narratives.

---

## 6. Config sketches (illustrative)

### 6.1 Free pulse bulk wind (A1) vs scramble (S0)

```
# S1 — organized
exp=pulse bath=1 L=16 dt=0.02 T=40
amp=2.0 sigma=3.0 kx=1.2
wf_on=1 k_rad=0
diag_every=50 snap_every=100 snap_file=runs/shear/S1.fcs
seed=20260802

# S0 — equal energy, random phase (seed hook or preprocess)
# same amp,sigma; kx ignored; fa phase scrambled; ⟨Ee⟩ matched
```

### 6.2 Off-ramp (S3) conceptual

A-family off-ramp is automatic when the pulse leaves the box if
absorbing/open edges exist; on periodic box, define t_off by
**zeroing fa when packet COM exits a window**, or switch to short T
and compare late-time Em after constructive interference dies.

B-family off-ramp is cleaner: `wind_amp→0` at t1.

### 6.3 Sustained shear sandwich (B3) — future knobs

```
field_wind=1 wind_amp=… wind_k=1.2 wind_mode=sandwich
# top y>L/2: +k; bottom: −k; book app_dE
```

---

## 7. Decision tree

```
Phase 0 A1 vs S0
    │
    ├─ S1 ≲ S0  ──► try A3 interface; if still null, SHEAR-G negative
    │                 (record: genesis not shear-organized under this door)
    │
    └─ S1 > S0 or interface hit
            │
            ├─ S3: Em dies after shear  ──► POSITIVE genesis, incomplete
            │                              theory (expected best case)
            │
            └─ Em holds after shear     ──► genesis + accidental residue;
                                           then (and only then) hand to
                                           stability reckoning R1–R2
```

---

## 8. Risks and pushbacks (design-time)

| risk | mitigation |
|---|---|
| S0 under-powered vs S1 | match ⟨Ee⟩; also match max Ee if door is local |
| Periodic box re-shears forever | use T shorter than wrap, or absorb BC, or explicit t_off |
| Workfn blocks all cond | confirm G0 dark; if all arms cond=0, lower wf or light structure seed secondary |
| Radiance erases birth | k_rad=0 first |
| Apparatus δE hidden | B-family books app_dE; fail run if unbooked |
| Overclaim stability | pass bars are genesis-only; dissipation after shear is success |
| Intuition lock-in | if S0 wins, stop; do not add shear meters until energy matched |

---

## 9. Deliverables checklist

- [ ] `report/analyze_shear_genesis.py` — correlations + off-ramp τ
- [ ] `runs/shear_g0.sh` — S0/S1/S2/G0 ensemble
- [ ] `runs/shear_g0/` logs + fcs
- [ ] `v93/theory/SHEAR_G0_RESULTS.md` — ratios, maps, pass/fail
- [ ] (optional) Phase 1 dual-pulse seed
- [ ] (optional) Phase 2 field_wind knobs + battery green at 0

---

## 10. One-sentence design

**Drive organized field transport (pulse wind, counterflow, or booked
boundary stream), hold energy-matched scrambled controls, meter Em
birth against phase-gradient and interface maps, and treat Em that
appears only while shear lasts — then dies — as a successful genesis
detection, not a failed particle.**
