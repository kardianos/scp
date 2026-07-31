# v89 — clean start

## The isolation rule

**Do not consult any version before v89 for concepts, models, code, methods, or
experimental design.** Not v88, not v87, not the census, not the Q-ball era.

This is not tidiness. Every prior version was built on a **background** — a
permanent index set with evolving contents — and that assumption was
reintroduced repeatedly, by different people, in different formalisms, *after*
being explicitly identified and rejected each time. It comes back through
inherited habits of construction, not through reasoned choice. The only reliable
defence is to not read the prior material.

Concretely, do not look at earlier versions for:

* how to represent state (they all use arrays indexed by immortal coordinates)
* how to write dynamics (they all use differential equations on a fixed index set)
* what a particle is (they all use profiles — shapes *in* something)
* instruments, seeds, or analysis tools (all assume a stage)
* parameter values, tuned constants, or "what worked before"

If a question seems to need prior context, the answer is to derive it from
`PRINCIPLE.md`, not to go looking.

## The standing constraints

1. **Energy is never destroyed.** It only changes mode. This is the only law;
   there is no separate equation of motion for it to constrain.

2. **No background.** Space is a mode, not a stage. A construction has
   reintroduced a stage if it contains: labels that persist while values change;
   a permanent reference geometry from which physical geometry deviates; fixed
   connectivity not itself produced by the current energy structure; or a
   quantity whose conservation requires fixed ambient space (ordinary momentum).

3. **No imported field or species.** Nothing may be added to carry an effect.
   If an effect requires a new field, that is a failure of the construction, not
   a licence to add one.

4. **Emergent from the fabric alone.** Particles, species, masses and
   interactions must arise from the fabric's own degrees of freedom.

## The law, as it stands (2026-07-28)

What the constraints have produced so far, in one place. One table of
constants — `battery/laws_V2g.cfg` — passes every experiment in the program
(20/20), with only apparatus differing between them: conservation, field
packets, heavy and light blobs, curvature linearity, CHSH, the rung tongue,
the pair tuning curve, the comma, the fifth, double-slit fringes,
single-quantum clicks, the eraser, HOM ordering, sub- and above-threshold
conversion, ħ-linearity, momentum of light, the gravitational footprint,
and occultation. The table is not a list of fixes; it states a short set
of laws, and it is guarded by the **ratchet rule** (`battery/README.md`):
every modification runs the full battery, and experiments that pass join
the gate.

1. **Amplitudes within a mode; atoms at the boundaries; integers in
   closure.** Inside a mode, energy moves as continuous amplitude —
   transport is never quantized. Every conversion between modes
   (space ↔ dense ↔ field) fires in whole atoms of action, ε = (A₀/2π)·ω
   at the source's pitch, carried by a credit that lapses at two atoms.
   This is measured, not assumed: every fired grain in every log sits on
   the ε(ω) grid to ~10⁻⁸ (E = ħ_eff·ω as a battery invariant), and a
   pool that cannot afford one atom never fires at all. The credit
   memory is physics, not bookkeeping — without it the comma cannot be
   paid.

2. **Load flattens pitch, and flight is load.** A voice's detune is
   x = (E_store + E_flight/2)/cap: bound energy re-pitches a cell, and a
   channel is a joint process of its two ends, so each end carries half
   of the channel's in-flight energy. Passing field amplitude does not
   re-pitch: vacuum optics is linear automatically, and Kerr
   nonlinearity is a property of loaded matter. Inertia emerges from the
   same law — light blobs translate, heavy blobs freeze.

3. **The choir's correction — derived, then retired.** A mistuned pair
   restores its own rung: the interference cross-flow of coherent
   exchange at the locked phase offset feeds the sharp voice and
   flattens it — restoring sign, ~2Γ acceptance window, and
   mutual-gate riding all *derived* in a two-voice amplitude micro-model
   with no bias law (`s2/`, CONSONANCE Part VII). At rate level the
   corresponding constant proved unnecessary: after flight-loads-pitch,
   the battery passes 17/17 with the bias at zero, so the law table
   carries no bias constant at all. Rungs are two-sided attractors
   without it — the fifth locks and holds. The raw derived term breaks
   tilted transport at rate level (16/17), which is now the acceptance
   criterion for the dense amplitude completion: it passes when
   translation IS the current.

4. **The vacuum skirt** — a scored prediction, confirmed. A voice loaded
   to within ~2Γ of the room's pitch dissolves into it; pair deaths
   retreat inside the computed boundary x_skirt = 2Γ/(q·(w₂ − 2Γ)), and
   the battery guards it.

5. **Momentum is the first moment of conversion.** Not fundamental — it
   is the energy centroid's velocity, and its conservation is the
   center-of-energy theorem: reciprocity of paired two-ended ledger
   moves plus statistical isotropy of the foam, replacing Noether.
   Measured: a field packet is ballistic through the foam (cos 0.999 at
   0.40 C) with its momentum current forward in the books; a sealed
   blob exerts no self-force. Where it fails is load-bearing: absorbed
   light's momentum does not survive the conversion (recoil deficit
   ~100×) — radiation pressure joins kappa_reac = 1 and moving-clock
   dilation as acceptance criteria for the amplitude completion.

6. **Space flows — pressure pushes, nothing reaches out.** Space is a
   mode like the others, so it moves through the same channels: each
   link carries space from its higher-pressure end at the channel's own
   conductance, outflow-limited, never drawing on the empty side (there
   is no suction). Pressure carries a displacement term from load —
   matter is converted space, so a mass pushes space out and MAINTAINS
   an extended graded depression: the gravitational footprint (core
   ~0.5 of vacuum, halo to r ≈ 6, a faint displaced-space ridge, no far
   tail). The throughput measurement closed the 1/r question for what
   we have: every far-shell space flux of a leaking blob is bookkeeping
   of its mass rate — outward wind while it evaporates, inward refill
   as its depression relaxes — decaying with the leak and subdominant
   to the radiative channel, with no steady accretion monopole. A
   sealed static mass sources no 1/r field; the far field requires
   steady space cycling at constant mass — a property of the stable
   particle the program does not yet have. The loaded core is a mirror
   by the same pitch law that mutes the slit wall — matter is
   emergently opaque, and a passing beam exits shifted away
   (occultation); the halo's true lensing pull is measured at the right
   sign, two orders below the foam floor (recorded).

Four placements of the atom competed as whole law tables (source- vs
destination-sized; per-cycle floor vs credit). V2 — source-sized atoms
with the two-atom credit — is the unique full pass (V1 15/17, V3 16/17,
V4 16/17), and the failures are discriminations, not accidents: a floor
cannot pay the comma, and destination-sizing taxes a translating front.
V2z retired the bias constant (§6.10); V2g added the space-transport
law (s_k, s_disp — both load-bearing) and is the standing table. Full
record: `ROADMAP.md` §6 and §7; harness and verdicts: `battery/`.

## What is inherited

Exactly two documents, and nothing else:

| file | role |
|---|---|
| `PRINCIPLE.md` | **foundational.** Energy is never destroyed; space is a *mode* of energy; matter is converted space; curvature follows from conservation. Everything else is subordinate. |
| `CONSTRUCTION.md` | the first background-free construction: state as a combinatorial energy complex, dynamics as rewrites, time as rewrite order, plus a background audit. |

Two programs are carried as working starting points, not as results:

| file | status |
|---|---|
| `construct_species.c` | enumerates closed integer locks. Produces a **finite** species count, which is the structural claim. The resonance rule it uses is **schematic**, so the count is not a prediction. |
| `construct_gauss.c` | combinatorial curvature: K − K₀ linear in converted energy, r² = 0.903, no ambient metric. The **linearity** is the claim; the coefficient is foam-dependent and not universal. |

Measured facts worth keeping are already embedded in those documents. Nothing
else needs to be retrieved.

## v89-native work

| file | role |
|---|---|
| `ROADMAP.md` | **program state, open areas, and the reality ladder** — the standing execution document. Tiers 0–3 of the double-slit ladder green (single-quantum clicks, duality equality, eraser, delayed choice); tier 4 (HOM) **partial**: exchange registry works, boson dip / fermion peak with correct signs and delay recovery (g_b 0.42 < 0.5 < g_f 0.58), depth limited by mode matching and coupler chromaticity — the same limits as real experiments. Forks F-ħ and F-NS resolved. **§6 = the standing first priority: unification** — one law table for all experiments; the parameter audit, the thesis (amplitudes within a mode, atoms at mode boundaries, integers in closure), construction U1–U6. **S1 executed (§6.8) then recursively annealed (§6.9): 17/17 under one law table** — variant V2 (source-sized atoms + credit) is the unique full pass; ħ-linearity a measured battery invariant; emergent inertia; A1 closed as a law (the choir's correction: dispersive bias × acceptance window × mutual gate closure) plus the zero-constant flight-loads-pitch law; the vacuum skirt a confirmed scored prediction. **§6.10: S2 first execution — the choir's correction derived at amplitude level (`s2/`, CONSONANCE Part VII) and its constant retired: 17/17 with zero bias (`laws_V2z.cfg` standing); kappa_reac=1 full-pass = the S2-full acceptance criterion.** |
| `battery/` | **the unification battery** — one `laws_V2z.cfg` shared byte-identically by 16 experiments (purity-checked), apparatus-only per-experiment cfgs, physics acceptance, variant cross-table. `python3 battery.py --laws laws_V2z.cfg`. |
| `s2/` | **the S2 derivation micro-model** (`choir_pull.c`): two-voice amplitude-complete channel; twin A (rate + posit) vs twin B (coherent, no posit); curve/heal/dscan sweeps in `s2/results/`. |
| `MASS.md` | **living lab notebook: the stable particle** — **first goal = exact mass (M-R1/P19)** (§5c); C1–C5; PRESTRESS waves; MASS↔EMF coupling; routes, instruments, ladder. |
| `EMF.md` | **living lab notebook: the field mode vs real EMF** — gaps, EM1–EM6; **§5 MASS↔EMF workflow** (parallel cheap lane, alternating law gates, joint C4/P2/g5). |
| `CHARGE.md` | **the charge construction (effects-first design, 2026-07-31)** — the QFT effect catalog E1–E11 (what charge DOES, quantified), the v89 mapping (Gauss = index theorem; quantization/universality = winding), ontology options (Poincaré-texture hedgehog / helicity-throughput imbalance of the flux-machine books / fractional shares on odd 3-cycles = the quark pattern), numerical means + ladder Q1–Q8. **Nucleon-class MASS is blocked on this** (user, 2026-07-31). |
| `INT_LEDGER.md` | **tentative integer energy-ledger simulator design** (int64 accounts + FP kinematics); does not block waves; FP64 remains production. |
| `bench/` | **the speedup ladder** — standard benchmark boxes + `BENCH.md` recording measured per-rung speedups (baseline: 17.9 ms/step single-thread, L=24 blob box). |
| `HBAR.md` | **the F-ħ situation, written out**: the cycle gate quantizes *when*, nothing quantizes *how much*. The reframe — the atom of transfer is action (fixed tail-call frame ⇒ E = (A₀/2π)ω) — plus ranked candidate origins (integer harmonic content restored; phase-space resolution from tongue widths; the space-cell grain e_s0·d̄/C ≈ 1.15 vs measured median act 0.617; zero-point self-consistency; winding), a rejected branch, and five discriminating tests. |
| `CELLFAB.md` + `cellfab.c` | the **cell-fabric kernel** (**default = integer ledger mode 3**, 2026-07-29): no-grid cells, two-plane harmonics, quant at mode boundaries. FP64 reference: **`cellfabf.c`**. Design/results: `INT_LEDGER.md`, `int_ledger/RESULTS.md`. |
| `CONSONANCE.md` | **locks without binding**: the musical resonance theory (notes as tail-calling processes, consonance as closure, the comma, the computed harmonic limit) back-ported to cell harmonics. Derives the pair separation ladder (ω₁+ω₂)d/C = 2πm from the existing gate law — no new mechanism — and measures it (E6/E7: tongue, tempered comma, defect payment, Huygens acquisition). Round 2 (Part VI) promotes the theory to kernel law (partial comb, roughness radiates, mutual coupling) and measures the comma paid (E8) and the interval lifetime hierarchy (E9). "Binding" is retired. |
| `DOUBLESLIT.md` | the double slit, postulated and simulated with instruments (wall, slits, screen, shutter, which-path recorders). Round 1 falsified the single-phase field sector exactly (I_AB = I_A + I_B, no fringes — refusal is delay, not cancellation). Round 2 **made the repair**: the field sector is now a signed two-plane amplitude (CONSTRUCTION §2.2's chiral pair) with symmetric-normalized unitary hops — **fringes at the parameter-free loci** (V_norm 0.32), complementarity as a measured dial (one strong meter: V → 0.025), 223 click grains, conservation ≤1.1e−15 through all detection. Dense-sector results unaffected. |

## What is open

From `PRINCIPLE.md` §7 and `CONSTRUCTION.md`:

* **The resonance rule.** Schematic at present. This is what stands between a
  combinatorial species count and a species *prediction*. Highest priority.
* **Direction and rate of conversion.** Conservation says energy changes mode;
  it does not say which way a transformation runs, or when.
* **Mode structure.** What distinguishes the patterning of mass from that of
  EMF; whether the mode list is complete.
* **Binding.** Whether an attractive channel exists that vanishes with
  separation. Not yet demonstrated in any form.

## Working rule

Before writing any model, state where its state lives and check it against
constraint 2. If anything in it persists and is merely re-valued, stop — the
stage is back.
