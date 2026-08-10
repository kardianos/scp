# v93 — the retention face: closed dense ring (route A) + DNLS map (route B)

2026-08-07. Executes RESUME §7 / ITEMS_OUTCOMES "the clear next face": item 4
diagnosed that the linear unitary channel diffracts a *localized* winding
(blob R2d 1.0→0.21 in 4 t.u.) and named two candidate rescues — topological
support (a closed cycle with |ψ|>0 everywhere) and DNLS nonlinear binding
(the load-detuning w₂e(Em) already in the law). Both routes are executed
here, to completion. Apparatus: `seed_mw` extended to `exp=ring` (th2 = m·φ
on the cycle, overriding the chain lock; BOTH kernels — the Go mirror also
gains seed_mw on the blob and the previously-unmirrored MOTION-#31 ring kx
tilt; byte-inert at 0 — full battery **ALL GREEN 87** after the edit), plus
two offline meters (`report/analyze_ring.py`: cycle winding W, seeded-m
coherence R_m, Em_min on the cycle, keep = Em_tag/Em_tot;
`report/analyze_melt.py`: Em_max, participation number PR, Em-weighted rms
about the packet centroid, all-cell COM x).

Common config (item-4 comparable): `L=16 dt=0.02 T=80 k_rad=0 wf_on=0 kx=0
seed=20260802`, ring nv=6: edge d*=1.5058 < link cutoff 1.955 < next-neighbor
2.61 ⇒ the vacuum ring is an exact C6 cycle graph; x=0.325 ⇒ Em=0.8125/voice,
w₂e=2.086 (empty-cell pitch 2.9). Unitary arms are the LINEAR channel
(amp_nat=2 amp_logate=1) unless noted. On nv=6 the m-spectrum is defined mod
6 (pk=−4 ≡ +2 is aliasing, not decay); with n=6 cells the R2 noise floor is
≈1/√6≈0.41 (measured 0.28–0.52 on the m0 control) — only R2 near 1 WITH W
locked counts as retention. Conservation (L1-B) held at the FP floor on
every run in this file (|drift_rel| ≤ 1.7e-14, most ~1e-15–1e-16).

## A.1 Vacuum closed cycle: the winding HOLDS — the existence proof

| arm | law | W | R2@80 | Em_min@80 |
|---|---|---|---|---|
| rv6_add_m2 | additive | **+2 exact to t=80** | 1.000 | 0.8125 |
| rv6_uni_m2 | unitary linear | **+2 exact to t=80** | 0.988 | 0.8089 |
| rv6_gat_m2 | unitary gated | **+2 exact to t=80** | 1.000 | 0.8124 |
| rv6_uni_m1 | unitary linear, m=1 | **+1 exact to t=80** | 0.955 | 0.8053 |

**The closed dense cycle retains the winding the blob could not hold** —
under the SAME unitary linear channel that scrambled the blob vortex to
R2d≈0.21 in 4 t.u., the C6 ring carries W=+2 exactly for ≥80 t.u. with
R2=0.988 (drift ~1.5e-4/t.u., consistent with the Em-breathing detuning
spread ±0.4% → δw₂e≈4e-3). The m=+2 state is an (approximate) eigenmode of
the cycle's hop graph; |ψ|>0 everywhere on the cycle ⇒ no phase-slip site ⇒
**topology beats diffraction**. This is the programme's first retained
matter winding, and it confirms item 4's diagnosis in its clean form. (The
gated arm holds near-frozen: the wound state's closure phase puts the gate
at cos⁸≈0.10, throttling exchange ~10× — consistent with item 1's
gate-suppression story. The additive arm holds because with only ring cells
present, κ_lock has nothing foreign to entrain against — contrast A.2.)

## A.2 Bath-embedded: the winding dies by t≈8–16 in BOTH laws

| arm | keep@80 | Em_min@80 | winding | note |
|---|---|---|---|---|
| rb6_add_m2 | 0.925 | ~0.003 | dead t≈8 (W→−1, R2 0.24) | additive |
| rb6_uni_m2 | 0.711 | ~0.001 | dead t≈8 (W→0, R2 0.22) | unitary linear |
| rb12_uni_m2 | 0.785 | ~0.01 | dead t≈8 (W→−3) | nv=12: finer cycle, 2× contacts — no rescue |
| rb6_uni_m0 (ctl) | 0.605 | — | R2 0.28–0.52 | the 6-cell noise floor |

Embedded in the jammed 3D bath (cavity carve clear=2.71; first shell ≥1.2
from the voices ⇒ linked), the winding is destroyed within ~8 t.u. under
BOTH transport laws, while most of the ring's energy is still in place.
Individual voices drain toward Em≈10⁻³ (phase-slip sites open) and the phase
order drops to the noise floor.

## A.3 Mechanism split (six diagnostic arms; predictions registered before each)

| arm | keep@80 | winding | lesson |
|---|---|---|---|
| freeze_geo=1 | **0.993** | dead t≈16, Em_min~0.1 | leak ≈ all NEW-contact accretion; winding dies anyway |
| q_detune=0 | **0.319** | dead | resonance barrier removed ⇒ ~10× leak (prediction ✓) |
| q_detune=3.6 | 0.578 | dead | NON-monotone vs law 1.2 (0.711): steeper drained-cell runaway (prediction ✗, recorded) |
| m=1 (margin 2.09 rad) | 0.706 | dead t≈16 | doubling the slip margin does not save it |
| additive m=1 | 0.984 | dead t≈16 | the additive killer is κ_lock entrainment toward ITS chain state, not leak |
| amp_nat=0.5 | **0.991** | dead t≈16, Em_min≥0.1 | leak ≈ eliminated; winding STILL dephases away |

Three separated mechanisms:

1. **Bulk leak = accretion crowding.** The carved cavity is a density
   depression; the bath presses in, new contacts form, and the ring exports
   Em through them (freezing topology keeps 99.3%). The leak through a
   ring→empty-bath contact is throttled by the **q_detune resonance barrier**
   (loaded voice w₂e=2.086 vs empty bath 2.9, res≈0.087): removing it (qd=0,
   everything resonant at 2.9) collapses keep to 0.32. But the barrier is
   **non-monotone**: a draining voice retunes toward the bath pitch, and
   steeper detuning makes that runaway steeper once started (qd=3.6 keep
   0.578 < law's 0.711). The law value sits near the benign range
   ("resonant melt" is self-limiting there).
2. **Phase slip without support loss.** With the leak eliminated
   (freeze_geo, or amp_nat=0.5) Em_min stays ≥0.1 — |ψ|>0 on the whole
   cycle — and the winding is STILL destroyed by t≈16. Bath contacts are a
   random-phase reservoir; on a 6-point discrete cycle an increment only
   has to wander through π on ONE edge to slip — no zero crossing needed.
   Doubling the margin (m=1) does not help, nor does halving the increment
   with 2× the contacts (nv=12). **Discrete small cycles have too little
   phase stiffness for jammed-bath contact noise.**
3. **The additive path loses the winding its own way** — κ_lock entrains
   the voices toward the bond-locked chain state (not the seeded winding),
   even with keep 0.98. Holding ≠ its preferred state.

**Refined requirement (supersedes item 4's clause):** |ψ|>0 on a closed
cycle is *necessary and sufficient in vacuum* (A.1) but **not sufficient in
medium** — retention there additionally needs phase stiffness ≫ contact
noise, or isolation, or an active mechanism that re-pumps phase order.

## A.4 Endurance + the vacuum ceiling: the kernel's own sweep order

T=1000 vacuum arms, unitary linear, all Em_min > 0 throughout (no support
holes anywhere in this section — every slip below is a pure phase slip):

| arm | W=+2 exact until | then |
|---|---|---|
| live geometry, tumble on (T1000) | **t≈200–250** | one slip → W=+1, slow walk; ring expanded 8.7% |
| live geometry, tumble OFF | t≈300 | slip; Em_min 0.69 |
| frozen scaffold (freeze_geo=1), tumble on | t≈50 | walk; INTERNAL condensation Em_min→0.07 at keep=1.000 |
| frozen scaffold, tumble OFF | t≈100 | walk; Em_min→0.01 at keep=1.000 |

Predictions registered: "if geometric noise is the vacuum slip driver,
freeze_geo holds to t=1000" — **REFUTED**, and the refutation is the
finding. With positions AND normals frozen and no stochastic term left, the
only C6-breaking perturbation remaining is the **canonical sequential order
of the pairwise Givens hops** (pass U applies slots in i<j sweep order — a
deterministic Trotter-ordering asymmetry). That seed asymmetry is then
AMPLIFIED by the q_detune resonance mechanism (a voice that gains Em detunes
away from its drained neighbour and keeps its gain — route B's condensation
operating *inside* the cycle): the frozen scaffold internally condenses
(Em_min→0.01 with zero leak) and the winding slips by t≈50–100. **Live
geometry is a symmetrizer** — bond feedback re-equalizes the edges and
extends the hold 4–6× (to t≈250–300); σ_tumble contributes only ~20%.

So the vacuum retention ceiling (~200–300 t.u., still 50–75× the blob's
scramble) is set by the kernel's own hop-ordering artifact — the same
Trotter family as item 1's dt-invariance failure. A symmetrized hop
schedule (Strang split / paired sweep / randomized order) is the identified
kernel-face candidate that would address BOTH open Trotter issues at once.

## A.5 Live-law arms (V3a) — refill never engaged; verdict unchanged

Full laws_V3a (k_rad=0.05 p_rad=4 wf_on=1), same ring, T=80: winding dead by
t≈16 in all three arms (additive; unitary; unitary+amp_door=1). **cond=0.000
in every arm** — the bath is dark and the workfn gate never opens, so the
condensation door NEVER FIRED: the amp_door=1 arm is byte-identical to
amp_door=0 (fcs `cmp` equal), and the hoped-for metabolic-refill rescue went
**untested**, not refuted. What was measured: the radiance tax alone (rad
0.27–1.38 over the run) does not change the retention verdict; the additive
V3a ring's Em decays 0.81→0.48 (taxed hoard below the x̂*=0.62 fixed point);
the bath compresses the ring (edges to d≈1.0–1.38, dev −0.5). A real refill
test needs a LIT bath (QUENCH-style warm medium) with door traffic on the
cycle — registered as the open arm.

---

# Route B — the DNLS map: the law's detuning binds ENERGY, not PHASE

`exp=blob bath=1 L=16 dt=0.02 T=80 sigma=2.5 kx=0 seed=20260802`, unitary
LINEAR (amp_logate=1), uniform-phase seed (seed_mw=0) unless noted. Meters
@t=80 unless noted; seed values Em_max=0.5634·(amp/0.5), PR=456, rms=4.30.
q_detune≠1.2 arms are off-law diagnostics (battery precedent: the xs bars
themselves run q_detune=0). For the two ADDITIVE controls the cell-sum
Em_tot undercounts by the in-flight parcel ledger (~⅓ at amp=0.5!) —
conv channels all zero, drift_rel at floor; cell-side meters remain usable.

## B.1 The discriminator: q_detune ON/OFF (amp=0.5, amp_nat=2)

| arm | Em_max@80 | PR@80 | rms@80 | verdict |
|---|---|---|---|---|
| **qd0** | 0.55 (transient 2.08@t4) | **853** | **6.79** | **MELTS** — pure linear diffraction spreads the packet |
| qd06 | 2.56 (pk 2.85@40) | 213 | 5.75 | condenses |
| qd1.2 (law) | 1.97 | 193 | 5.31 | condenses |
| qd3.6 | 2.12 | 224 | 4.83 | condenses, tighter |
| qd12 | 2.00 | 253 | **4.48** | condenses, tightest envelope |

**The load-detuning w₂e(Em) is a real self-trapping nonlinearity for the
unitary channel.** Remove it (qd=0) and the packet disperses (PR 456→853,
rms→6.8). At any qd ≥ 0.6 the packet CONDENSES: peak grows 3.5–3.6×
(amplitude-ratio roughly amp-independent: am025 0.28→1.01; am05 0.56→1.97;
am1 1.13→3.54), participation halves, and envelope confinement strengthens
monotonically with qd. Mechanism (same as A.4): **resonant rich-get-richer**
— a cell that gains Em detunes away from its poorer neighbours (res
collapses), decouples, and keeps its hoard; a draining cell retunes toward
the bath and drains resonantly. Phase-blind condensation, not coherent
focusing. Condensation appears at every hop strength (amp_nat 0.5–4:
Em_max@80 1.45/1.84/1.97/1.59, PR 248/222/193/205).

## B.2 The deep corner: full self-trap (and its cost)

| arm | Em_max | PR@80 | rms 0→80 | note |
|---|---|---|---|---|
| am2 (law qd) | 5.40@16 → 4.24 | 295 | 4.30→4.70 | evap 16.6 (cap bleeds the over-cap hoards) |
| qd36am2 | **7.75@40** | 260 | 4.30→**4.32 (flat)** | fully self-trapped; −0.75% Em_tot/80 t.u. |
| qd12am2 | 7.19 | 213 | 4.30→**4.25 (flat)** | frozen hoard (−0.06%; w₂e≈0.08, all channels quiet) |

At (qd≥3.6, amp=2) the packet is **fully self-trapped — zero envelope
spread over 80 t.u.** The existence region of non-melting packets is
mapped: it opens between qd 0 and 0.6 and deepens monotonically; depth
(amp) mainly raises the hoard peaks. The hoards run to 3× cap — pass U
conserves but does not cap; the law defends cap by EVAPORATING overflow to
the field sector (am2: −4.4%/80 t.u.), except at extreme detuning where
even that channel quiets (qd12am2: frozen).

## B.3 Mobility: self-trap and translation are antagonists

| arm (amp=0.5, law qd) | Em_max@80 | PR_min | rms@80 | COM x (start 7.99) |
|---|---|---|---|---|
| kx=1.1 | 2.04 (pk 2.54@16) | 151@16 | 6.15 | **9.31@t4, then back to 8.80** |
| kx=2.0 | 2.70 | 165 | 5.99 | 8.93@t4 → 8.72 |
| kx=2.6 | 1.65 | 231 | 6.39 | 8.41@t4 → 8.31 |

A moving packet condenses HARDER (kx=1.1 is the strongest condenser of the
amp=0.5 family) and then **stalls**: the all-cell COM advances ~1.3 cells in
the first ~4 t.u. and then anchors/reverts as the hoards decouple. No
band-edge melt suppression distinct from this: the kx axis shows the
textbook DNLS mobility/trapping tradeoff — **the corner that binds is a
corner that does not translate** (and vice versa: item 2's coherent current
lives in the un-bound regime).

## B.4 Retention at the self-trap corners — the direct probe

seed_mw=2 blob vortex at the condensing corners (predictions registered:
condensation fragments phase → dies like the law corner):

| corner | R2d@t8 | R2d@t80 |
|---|---|---|
| item-4 law corner (reference) | 0.07 | 0.07 |
| qd3.6, amp=0.5 | 0.044 | 0.05 |
| qd12, amp=0.5 | 0.137 | 0.04 |
| qd12, amp=2 (frozen hoard) | 0.134 | 0.01 |

**Dead everywhere, prediction ✓** — the winding scrambles during the
condensation transient itself, before any kinetic freezing could preserve
it. Item 4's "possibly DNLS binding" clause is now CLOSED by measurement:
**the law's nonlinearity binds energy, not phase.**

## B.5 Additive controls

No condensation under the additive law (Em_max stays ≤ cap: 0.85 / 2.33;
PR 298/340; rough=9.0 at amp=2, plus the in-flight ledger noted above).
Condensation is a unitary-channel phenomenon — conservative resonant
decoupling, unavailable to the dissipative want path.

---

# Synthesis — what the retention face measured

1. **Topological protection is real and now demonstrated** (first retained
   matter winding: W=+2 exact, vacuum C6, ≥200 t.u.) — but its medium
   version fails: contact noise dephases small discrete cycles in ~10 t.u.
   regardless of law, margin, hop strength, or cycle size at this scale.
2. **The DNLS route binds the wrong thing.** The q_detune detuning is a
   genuine, mapped self-trapping nonlinearity (qd0 melts / qd≥0.6
   condenses / deep corner fully trapped) — but it traps by incoherent
   fragmentation: phase retention dead at every corner, and mobility dies
   with it.
3. **Transport, binding, and phase retention are pairwise antagonists in
   the current channel**: linear ⇒ transports, diffracts, doesn't bind;
   detuned ⇒ binds, stalls, fragments phase; closed cycle ⇒ holds phase,
   only in isolation.
4. **The vacuum retention ceiling is the kernel's own hop-sweep order**
   (A.4) — the same Trotter family as item 1's dt artifact. One kernel face
   (symmetrized/Strang hop schedule) targets both.
5. **Spontaneous condensation is a creation-adjacent discovery in its own
   right**: under the unitary channel + the law's own detuning, a spread
   packet spontaneously condenses into long-lived, envelope-frozen dense
   hoards — conservative clumping with no lock, gate, or door (the additive
   law cannot do this). This echoes QUENCH ("dynamics created what gates
   could not") and is the natural bridge to the programme's mass-formation
   story, independent of the spin-retention goal it was aimed at.
6. **Open, registered:** the lit-bath metabolic-refill arm (A.5 — the door
   never fired; amp_door untested where it matters); the symmetrized hop
   schedule; tolerance-C≡Go at amp_nat>0 (all results here are C-only, per
   the standing note).

## Files

- Kernel: `kernel/freecell.c` (ring seed_mw block, marker `v93
  ring-retention face`); Go mirror `fab/run.go` (ring kx tilt + seed_mw,
  blob seed_mw), `fab/params.go` (SeedMw). Byte-inert at 0; battery ALL
  GREEN 87.
- Meters: `report/analyze_ring.py`, `report/analyze_melt.py`.
- Runs: `runs/ringret.sh` → `runs/quench/rv6_*, rb6_*, rb12_*` (+ diagnostic
  arms `rb6_uni_m2_{fg,qd0,qd36,an05}`, `rb6_{uni,add}_m1`, endurance
  `rv6_uni_m2_{T1000,fgT1000,fgntT1000,ntT1000}`, live-law `rbV3a_*`);
  `runs/dnls_map.sh` → `runs/dnls/*.{log,melt,fcs}` + retention probes
  `runs/dnls/ret_*`. Console summaries: `runs/ringret.out`,
  `runs/dnls_map.out`.
