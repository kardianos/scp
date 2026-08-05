# CANTUS — the H campaign: the superimposed harmonic-lock field
## (coherent-channel candidate B — "atoms are NOT cells")

Opened 2026-08-04 by the v91 session following the R campaign, the pad
campaign, and INTEGRATION. **User directive (2026-08-04, verbatim):**

> "Atoms are NOT cells. I may suggest superimposing a field fitted to
> the cells which can hold it's own gauges, to simulate the overall
> harmonic lock. Try that along with or in addition to the other
> possible ways forward to the goal of stable atomic mass and other
> endevors to make the model corraspond more closely to reality."

This authorizes the S2-lane law-candidate implementation (kernel change
under the ratchet). Format precedent: `RADIANCE.md`. **§1–§3 of this
document are written to disk BEFORE the first k_cant>0 experiment
runs.** §4+ records what happened. The candidate is implemented INERT
(defaults k_cant=0, k_tune=0 ⇒ byte-identical to the current kernel);
the battery stays green at defaults throughout.

## 0. The thesis, and what is being repaired

**Atoms are not cells.** A cell is substrate. The atom-analog is a
GROUP of cells sharing one mode — and the standing law already says
where mode-sharing must live: *amplitudes within a mode; atoms at mode
boundaries*. What the dense sector lacks (measured, S2_CHANNEL.md) is
any way for N bonded cells to BE one mode: transfers carry magnitude
only, phases only gate, so every voice jitters alone and the ensemble
is N doors, not one body.

The measured killers of small/curved matter at the radiance balance
(RADIANCE.md §4.7–4.8, INTEGRATION):

* **K1 — gate-break in vacuum at x̂\*:** each radiance atom kicks
  holdings ~10–20%, the pitch w2e jumps ~4–9%, the bond phase-closure
  ψ random-walks, and the cos^8 gates die (ring6@0.62 vacuum: ret
  0.91, ALL SIX channels dead).
* **K2 — pile-pump-out in bath:** Em piles unevenly (Jensen on x⁴ +
  flight-loaded det + roughness), the piled cells radiate ~2.3× the
  single-cell law, and the object boils off (ring6/comp6 t_half
  80–140 vs k0's 260–510).
* **K3 — the negative template:** the rate-level correction
  (kappa_reac=1) retires the fifth wherever the fifth is load-bearing
  (FCQ ggm→0, D bloats to x̄ 0.83, composite boundary dies). Any
  completion that flattens p:q structure toward unison is dead on
  arrival.

**Candidate B (CANTUS):** superimpose an order-parameter field FITTED
TO THE CELLS — per-cell state, living and dying with the cells' own
bonds, no background — which (i) holds the object's gauge frame (the
p:q chord chart on its live links + a per-voice pitch memory), and
(ii) enforces the overall harmonic lock two ways: a Kuramoto-style
phase consensus on the matter clocks (not the transfer rates — the
anti-K3 decision), and a within-mode retuning current that flattens
fast holdings fluctuations while preserving the slow chord structure
(flavor). The name: cantus firmus — the held line the voices sing
around.

## 1. The design, precisely (implemented in both kernels, inert)

### 1.1 New per-cell state (allocated always, zeroed; no energy content)

* `ca[i]` ∈ [0,1] — cantus amplitude, the local order parameter. A
  pure gauge/plasticity object (kernel precedent: P15 moves ld with no
  ledger entry; geometry is not a ledger). It carries NO energy.
* `cxl[i]` — holdings memory: low-passed x_hold = Em/cap — the same
  base the radiance law reads (R-D1: taxes read holdings). This is
  each voice's remembered "part" in the chord. Flight is deliberately
  NOT in the memory: the holdings/pitch split is the built-in
  regulator INTEGRATION measured, and the cantus must serve it, not
  erase it.

Initialization: `ca = 0`, `cxl = Em/cap` at seed time. Apparatus key
`cant_seed=1` additionally sets `ca = 1` on tagged object voices at
init (constructed-coherent start); default 0 (self-growth only).

### 1.2 Pass H (new; after pass 5 deliveries, before pass D motion)

Entirely guarded by `(k_cant > 0 || k_tune > 0)`; when off, zero FP
work, zero RNG, byte-identical step.

Eligibility: live channel s (`sst != FREE`, `sA > 0`) with
`Em_i > 1e-15 && Em_j > 1e-15`, using the slot's standing comb chart
(p, q) = (slp, slq). Per eligible link, with the LIVE ladder phases:

    psi_f = wrap_pi(q·th2_i − q·w2e_i·d/C − p·th2_j)
    psi_b = wrap_pi(p·th2_j − p·w2e_j·d/C − q·th2_i)
    gg    = gate_of(psi_f) · gate_of(psi_b)      (same LUT gate)

**(a) Support (slot loop 1, read-only on state):**
`sup_i = max(sup_i, gg)` and `sup_j = max(sup_j, gg)` — support is
the best two-sided gate quality among the cell's live dense links.

**(b) The harmonic lock (slot loop 2, sequential in slot order, both
kernels identical order):** amplitude-gated, GATE-INDEPENDENT (the
flywheel decision C-D5):

    wl   = k_cant · dt · sqrt(ca_i · ca_j)
    e    = 0.5 · wrap_pi(psi_f − psi_b)     (differential error)
    th2_i −= wl · (q/(p²+q²)) · e
    th2_j += wl · (p/(p²+q²)) · e

so that de = −wl·e exactly (the (q,p)-metric normalization). Note the
phase sum psi_f + psi_b contains NO phase at all — it is the ladder
misfit −(q·w_i + p·w_j)·d/C (mod 2π), pure geometry: the phase DOF can
only ever fix the differential error, and the bond's plastic d-walk
can only ever fix the common one (decision C-D4 is a theorem of the
closure algebra, not a choice — phases keep the gates closed, geometry
keeps the rung; the two channels cannot compete).

**(c) The within-mode retuning current (same slot loop 2):** the
locked ensemble is ONE extended dense mode (the user's "atoms are not
cells"); within a mode, transport is never quantized (the standing
law), and the field sector's pass-F local hops are the precedent for
un-retarded within-mode spatial redistribution. Drive on the
pitch-memory ERROR — not the level — so flavor structure (D's heavier
part) is preserved while fast piles flatten:

    u_i = Em_i/cap − cxl_i              (fast deviation from own part)
    J   = k_tune · dt · sqrt(ca_i·ca_j) · cap · (u_i − u_j)
    clamp: J ≤ 0.98·Em_src ; J ≤ max(0, cap − (Em_dst + Ee_dst))
    Em_i −= J ; Em_j += J               (pairwise-conserving, exact)
    tune_total += |J| ; p1 site H: p1fl += J · d · u_hat (COE meter)

**(d) State update (cell loop, canonical order):**

    cxl_i += (dt/cant_tau) · (Em_i/cap − cxl_i)
    ca_i  += (dt/cant_tau) · (sup_i − ca_i),  clamped to [0,1]

ca is a low-passed record of the best gate the cell actually holds:
it grows only where locked exchange PERSISTS (transient bath contacts
< cant_tau build nothing), and it decays in ~cant_tau when the bonds
die (no zombie locks — the field dies with the object; no background).

### 1.3 New keys (one header line; battery law-purity pins defaults)

`k_cant` (lock strength, default 0), `k_tune` (retune rate, default
0), `cant_tau` (memory time, default 50), `cant_seed` (apparatus,
default 0). Header: `# v91 cantus (coherent-channel candidate B): ...`
printed always; `# CANT` diag line and `# RESULT cantus` only when
enabled (logs at defaults stay byte-identical).

### 1.4 Decision points (recorded before running)

* **C-D1 (level):** the lock acts on the matter CLOCKS (th2), never on
  the wants. kappa_reac acted on rates and killed the fifth (K3); a
  clock consensus leaves the want formula untouched — transport
  through an open gate is still magnitude/headroom-driven, so U↔D
  exchange through a locked fifth should INCREASE, not die.
* **C-D2 (amplitude source):** ca = low-passed best-gate-quality.
  Alternatives (transfer-weighted, per-link amplitudes) recorded as
  refinements; max-over-links is the minimal object that
  distinguishes "bonded member" from "bath transient".
* **C-D3 (within-mode current legality):** justified by the mode
  argument + pass-F precedent; un-quantized by the standing law
  (transport within a mode). It enters the COE meter (p1 site H).
  If R2-class runs show the current carrying net momentum artifacts,
  revisit with a flight-routed variant (the ablation).
* **C-D4 (differential/common split):** phases correct only the
  gate-killing differential error; the ladder misfit stays with the
  bond walk. No competition between the two channels.
* **C-D5 (flywheel gating):** lock and current scale with ca (memory),
  NOT with the instantaneous gate — a gate-gated lock releases exactly
  during the dips it exists to bridge.
* **C-D6 (what must not change):** atoms machinery, conversion door,
  space sector, contact rule, radiance candidate A — untouched. The
  cantus adds no energy mode; conservation stays at the FP floor.

## 2. Closed-form expectations (from the implemented term, before any run)

* **Jitter budget at x̂\* = 0.62** (k_rad=0.05, p_rad=4): radiance atom
  ε(x̂\*) = 0.304 ⇒ holdings kick δx = 0.122; pitch jump δw2e =
  w2·q_detune·δx/det² ≈ 0.14 ⇒ closure drift rate |ψ̇| ≈ q·δw ≈ 0.14
  rad/t.u. between fire and recondensation.
* **Lock steady-state error** e_ss ≈ ψ̇/(k_cant·ca) (first-order PLL):
  k_cant=0.3 ⇒ e_ss ≈ 0.47 rad ⇒ gg ≈ 0.40 (marginal); k_cant=1 ⇒
  0.14 rad ⇒ gg ≈ 0.93; k_cant=3 ⇒ 0.047 ⇒ gg ≈ 0.99. **The gate-
  survival threshold should sit between k_cant 0.3 and 1** — the H1
  sweep {0, 0.3, 1, 3} brackets it.
* **Retune flattening time** τ_flat ≈ 1/(2·k_tune): k_tune=0.2 ⇒ 2.5
  t.u., well under the pile-buildup scale (tens of t.u.). Sweep
  {0, 0.05, 0.2}.
* **Memory time** cant_tau must exceed the gate-dip/recovery scale
  (~10–20 t.u. at x̂\*): default 50; sweep {20, 150} at the selected
  k_cant.
* **Stability:** per-step phase correction ≤ k_cant·dt·(degree) ≈
  0.14·err at k_cant=1, deg 7 — far from the ~2/step instability; the
  scheme is safe to k_cant ~ 10.
* **The fifth keeps its chart:** for a U–D fifth (p:q = 3:2) the lock
  corrects 2·th_U − 3·th_D differentially — it enforces the 3:2
  closure, not unison; and the retune current drives on (u_U − u_D)
  with u the DEVIATION from each voice's own memory, so the D voice's
  standing heavier load is not a drive term. K3 should not reproduce.

## 3. Pre-registered predictions (BEFORE first k_cant>0 run)

* **P-H1 (byte-inertness):** at defaults the kernels are byte-identical
  to the pre-cantus binaries on every arm (bath diff empty, battery
  93 ALL GREEN, abx green), and the committed R2/pad logs remain the
  valid k_cant=0 controls for every H2 arm.
* **P-H2 (wiring parity):** a k_cant>0, k_tune>0 pair run prints
  identical `# RESULT cantus` numbers in C and Go (the rad=4.584987
  precedent), drift at the FP floor.
* **P-H3 (the lock holds gates under thermal fire — K1 repaired):**
  ring6@0.62 VACUUM k005 with cantus at the selected point: ≥5 of 6
  ring channels alive (gg > 0.5) at T=20000, vs the committed control
  0 of 6 — with ret ≥ 0.85. Threshold behaviour per §2: k_cant=0.3
  marginal, 1 and 3 alive; failure of all three falsifies the clock-
  consensus mechanism.
* **P-H4 (the current flattens the pile — K2 repaired):** ring6@0.62
  BATH k005 + cantus: t_half(ret 0.5) ≥ 10× the measured k005 hoard
  ceiling 140 ⇒ ≥ 1400, target alive at T=5000; the tagged outtake
  tracks the single-cell law (ledger out within ~1.3× of 6× the
  per-cell D(x̂\*) + intake balance) instead of the 2.3× Jensen
  amplification; windowed in ≈ out (flux-machine bar).
* **P-H5 (small matter lives at the sweet spot):** the i5 arm (nv=6
  ring, x=0.47, amp=0.15, k005) with cantus: t_half ≥ 10× the
  measured 230 ⇒ ≥ 2300, target alive at T=5000 with ≥ 4/6 bonds.
  This is the acceptance-surface row the nv<36 zoo lives or dies by.
* **P-H6 (the fifth survives coherence — K3 does NOT reproduce):**
  UUD triad (tri_xU=0.28, bath, k005) + cantus: ggm ≥ 0.9 at T=2000
  with no D bloat (x_D within ±0.1 of its chart value); composite
  nv=6 (rings_kind=1): boundary fifth gg ≥ 0.9 at T≥100 and D-ring
  x̄ ≤ 0.35 (the S2_CHANNEL.md quantified faces). If instead the
  fifth dies or D bloats, candidate B inherits kappa_reac's disease
  at clock level and is rejected as-is.
* **P-H7 (self-growth, no hand-seeding needed):** with cant_seed=0,
  a bonded object's ca grows from 0 to ≥0.7 within ~3·cant_tau of
  bond formation (measured on the H1 pair); with cant_seed=1 the
  early transient shortens but the t ≥ 1000 behaviour is the same
  (the field is self-sustaining, not apparatus-sustained).
* **P-H8 (the bath is not rewritten):** a no-seed BATH run at the
  selected (k_cant, k_tune): the glow economy stays within the
  seed-to-seed scatter of the k005 baseline (global cond/rad within
  ±15%), and system-spanning lock does not form (cells with ca > 0.9
  < 5% of NC). The cantus must stabilize OBJECTS, not crystallize the
  vacuum.
* **P-H9 (conservation + purity):** |drift_rel| ≤ ~1e-13 every run
  with the term firing; battery at defaults ALL GREEN after
  implementation (P-H1); no bar moves at defaults.
* **P-H10 (size selection gets its measurand):** if P-H3/P-H4 pass,
  surviving ring bond lengths land on the live pitch ladder rung
  (d\* = πC/w at the measured holdings+flight x, the R3 zoo table) —
  the R3 "no measurand" verdict is repaired and sizes become derived
  quantities under the coupled law.

**Pre-registered selection rule for (k_cant, k_tune, cant_tau):** the
smallest k_cant (then smallest k_tune) that passes P-H3 AND P-H4 with
margin (gg ≥ 0.7, t_half ≥ 10×), with cant_tau kept at 50 unless the
{20,150} probes show a measured preference; ties broken toward fewer
active knobs (k_tune=0 preferred if the lock alone passes both).
Algebraic expectation: (1.0, 0.2, 50). The selection is made by the
measured curves, not by this expectation.

### 3.1 The H arms (protocols pinned to committed controls)

| arm | protocol (all + `k_rad=0.05 p_rad=4 rad_clock=0` unless noted) | control (committed) |
|---|---|---|
| H0 | build both kernels; battery; bath byte-diff C≡Go; wiring pair | `runs/BATTERY.log`, VERIFY protocol |
| H1a | pair unison x=0.62 vacuum (bath=0) noise_amp=0.5, T=2000, k_cant ∈ {0,0.3,1,3} × k_tune ∈ {0,0.2}, cant_seed ∈ {0,1} probe | k_cant=0 arm |
| H1b | pair fifth (pair_pp=3 pair_qq=2, x_U=0.28, x_D chart 0.92) bath=1, T=1200, same grid | k_cant=0 arm |
| H2a | ring6 x=0.62 VACUUM T=20000 at selected + one-below | R-campaign §4.8 (0/6 alive) |
| H2b | ring6 x=0.62 BATH noise 0.5 T=5000 at selected | §4.7 (t_half ≈ 140) |
| H2c | comp6 rings_kind=1 xU=0.28 bath T=5000 at selected | `runs/r2/comp6_bath_k005.log` (t_half ≈ 85–90) |
| H2d | UUD tri_xU=0.28 bath T=2000 at selected | `runs/r2/uud_bath_k005.log` (t_half ≈ 140) |
| H2e | i5 rerun: ring6 x=0.47 amp=0.15 bath T=5000 at selected | `runs/pad/i5_nv6_x047_a015.log` (t_half ≈ 230) |
| H4 | bath no-seed T=120 + T=480 at selected (P-H8); ca census | R1b k005 bath numbers |
| H5 | full battery `-extra "k_cant=… k_tune=… "` — the reckoning | battery at defaults |

Logs land in `runs/cantus/`, named `h<arm>_<config>.log`.

### 3.2 Addendum (2026-08-04, post-wiring-probe, pre-campaign — dated; nothing above rewritten)

The H0 wiring probes (pair vacuum x=0.62 and x=0.325, k_cant=1
k_tune=0.2) and a re-read of the committed R-campaign log they pointed
at expose a measured fact that re-aims the campaign:

**The R-campaign's ring6@0.62 VACUUM arm was dead at birth.**
`runs/r2/ring6x62_vac_k005.log` shows `slots=0, live=0` from t=0 to
t=20000 — the six voices were NEVER connected. Mechanism, now derived:
seeding a voice at x draws s_pull·x·cap from its own space store, so a
vacuum-seeded voice has Es = e_s0 − s_pull·x·cap and radius
r = r0·Es^{1/3}. At x=0.62: Es=0.7675, r≈0.778, candidate cutoff
cfac·2r ≈ 1.790 < d*(0.62)=1.889. **The unison rung at the fixed
point exceeds the candidate cutoff: unison matter at x̂\* cannot even
CONNECT, let alone bond.** RADIANCE §4.8's "thermal atom-kick jitter
breaks the gates" was the wrong mechanism attribution for that arm
(nothing existed to break); OUTLOOK's empirical bond ceilings
("exist only below x≈0.52") are this inequality, softened by Es
breathing in live environments. Thresholds at vacuum-seed Es: lens
contact 2r(x) crosses d*(x) near x≈0.30; the candidate cutoff
cfac·2r(x) crosses near x≈0.52.

**The chord consequence (new, pre-registered):** interval rungs
shrink as p+q grows — d* = 2πmC/((p+q)·w) at matched pitch scale. At
x̂\*=0.62 the unison m=1 rung (1.889) does not fit contact, but the
FIFTH m=2 line (≈1.51 at chart loads) and every deeper consonance DO.
**The comb admits only multi-pitch matter at the balance — atoms at
x̂\* are forced to be chords, not unison choirs** ("atoms are not
cells", now as a selection rule). Note the FCQ UUD chart at
tri_xU=0.28 has mean load (0.28+0.28+0.837)/3 = 0.466 ≈ the
INTEGRATION sweet spot, and a fifth pair at (0.35, 0.94) has mean
0.645 ≈ x̂\* — the standing flavored objects already sit near the
balance in the mean. The flux-machine picture sharpens: light voices
absorb (headroom), the held fifth conducts, the heavy voice radiates
(x⁴) — a two-temperature metabolizing atom.

**Wiring-probe results (recorded):** (i) pair@0.62 vacuum: no
eligible channel ever (support 0, ca decays exp(−t/τ) exactly, tune
0) — the contact rule is respected; the lock cannot and does not
reach out. (ii) pair@0.325 vacuum + noise + k005, cant_seed=0: ca
self-grows 0→0.98 (P-H7 behaviour confirmed at probe level), tune
fires (0.94 units over T=400), drift 1.0e-13, and the FIRST-SIGNAL
numbers move the pre-registered way: holdings symmetric 0.940/0.941
(control 0.810/0.972) and rung offset +0.026 vs control +0.096
(3.7× tighter). Not yet a claim — single seed, probe length.

**Amended arms (replacing the §3.1 rows where they conflict):**

| arm | protocol (all + k005 radiance) | control |
|---|---|---|
| H1a | pair OF the bath (bath=1) x=0.47, noisy ambient, T=1200, k_cant ∈ {0, 0.3, 1, 3} × k_tune ∈ {0, 0.2} | k_cant=0 arm, same seed |
| H1b | fifth pair vacuum+noise: pair_pp=3 pair_qq=2, x_U=0.28, x_D=0.8367 (chart), m=2, T=1200, grid as H1a | k_cant=0 arm |
| H2a | ring6 x=0.47 VACUUM noise_amp=0.5 T=20000, selected point + k_cant=0 control (no committed control exists at 0.47 vacuum) | in-arm |
| H2b | ring6x62 bath (as committed protocol) at selected — SECONDARY (skin-lock observation; unison interior cannot connect) | committed log |
| H2c/d/e | unchanged (comp6, UUD, i5 — the PRIMARY life/death bars) | committed logs |
| H3 | the chord at the balance: fifth pair (3:2, m=2) at (x_U, x_D) = (0.35, 0.94), bath, T=5000, selected point + k_cant=0 control — pre-registered: with cantus the chart holds (pitch ratio stays 3:2±2%, both voices alive, windowed in≈out); without, D drains/detunes | in-arm |

P-H3 is restated onto H2a (x=0.47 vacuum): ≥5/6 ring channels alive
(gg>0.5) at T=20000 at the selected point. Its control may also
survive (the 0.325 vacuum k005 ring did) — if both live, the K1 bar
moves entirely to the bath arms (H2c/d/e), where death is measured.
P-H4's flux-machine ledger bar moves to H2e (i5) and H2c (comp6).
P-H6 and all other predictions stand as written.

---

## 4. Results

(§4 is filled ONLY after the corresponding runs complete; nothing
above this line is edited after the first k_cant>0 run.)
