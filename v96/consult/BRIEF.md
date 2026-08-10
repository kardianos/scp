# v96 — claude developmental session: the ground state (size selection + quark ontology)

You are the **developmental session** for SCP v96. This brief is **self-contained** but assumes the v89–v95 record. You have read access (`--add-dir` to v94, v95, v93, v89).

Produce `v96/consult/DESIGN_claude.md` — theory, derivations, numerical simulation method, bars. Print your session ID at the head.

---

## THE USER'S TWO QUESTIONS

### Question 1: why does matter have a specific size?

> The size must be correct. Think about why a proton and neutron are a specific size. Everything else is meta stable.

The proton is stable at ~0.84 fm because confinement (volume) balances kinetic pressure (surface). Everything else — the Delta, excited baryons, resonances — is meta-stable because it's NOT at the balance. The programme has created locally-stable configurations (q_detune self-traps, metabolic hoards) but has never identified the balance that selects a UNIQUE size. The v93 Face C hoards showed a hierarchy (0.5–7.2); RECKONING correctly identified them as transients. The M6 metabolic park is a load fixed point (x* ≈ 0.62), not a size fixed point.

**v96 must derive the liquid-drop balance:** volume intake (∝ N, from the metabolic cycle / q_detune trap) vs. surface leak (∝ N^(2/3), from contact drains). Find the N where total energy is minimized. Show objects at other N are meta-stable. Show N emerges from law constants. The J-metric (v95) modulates the surface term (starves contacts) — derive how this shifts the selected N.

### Question 2: are quarks distinct particles, or emergent attributes?

> Are quarks distinct particles or are they 3 attributes that emerge after a proton/neutron is formed, and we see them emerge transiently when they collide and we break them apart?

This is the deepest ontological question the programme faces. The pre-v89 "quark" model (multi-component fields on a fixed lattice) was **explicitly rejected** — it was a background-dependent, separately-ledgered, identity-carrying fiber. But the QUESTION remains: when a proton is probed at high energy (deep inelastic scattering), three sub-structures appear. Are they:

**(A) Distinct pre-existing components** — the proton is made OF three things, which exist before and after the proton. This is the QCD quark picture. In the SCP substrate, this would mean three separately-evolving dense amplitudes per cell (the rejected quark road).

**(B) Emergent attributes of the whole** — the proton is ONE object. When you break it (collision), the energy of the break distributes into three channels (the three lowest modes of the broken object). The "three-ness" is a property of the DECOMPOSITION, not of the ORIGINAL. Like the harmonics of a struck bell: the bell is one object, but its vibration decomposes into modes. The "quarks" are what the fragments look like when you smash the bell.

**(C) Something else** — the substrate's own answer, which may not be either (A) or (B).

The user's intuition leans toward (B): quarks are what EMERGE when you break a proton, not what it's MADE of. The proton forms first (from the fabric), and the three-ness is a property of its internal structure (three phase sectors? three winding modes? three beat-frequency channels?) that becomes visible only under destructive probing.

**v96 must address this ontologically and physically:**
- What is the SCP substrate's answer to the quark question?
- Can the "three-ness" emerge from the substrate's own structure (three sectors: field/matter/space; or three voices on a ring; or three phase-winding modes)?
- Is there a falsifiable prediction that distinguishes (A) from (B) in the substrate?
- How does this connect to the CLAUDE.md prohibition on the "quark road" (separately-ledgered components)?

---

## WHAT TO PRODUCE

### I. Theory: the ground state and the size selection

1. **Derive the liquid-drop balance** from the SCP law constants:
   - Volume term: the metabolic cycle's intake rate per interior cell (∝ x* from the M6 park).
   - Surface term: the contact leak per boundary cell (∝ the contact drain rate, modulated by the J-metric).
   - The binding energy curve: E(N) = -a_vol·N + a_surf·N^(2/3) + ... (derive the coefficients from the law).
   - The ground state N*: where dE/dN = 0. Show N* is unique and emergent.

2. **Show meta-stability of other sizes:** objects at N < N* evaporate (surface dominates); objects at N > N* are the "box-fill" runaway (volume dominates). Only N* is stable.

3. **The J-metric's role:** the surface term is modulated by jbed (contact starvation). Derive how N* shifts with j_metric. This is the SIZE TUNING prediction — falsifiable.

4. **Connect to the proton/neutron:** in the SCP substrate, the "proton" is the ground state at N*. The "neutron" might be the same N* with a different phase configuration (seed_mw = 0 vs. 1, or a different winding). The size is the SAME (both proton and neutron are ~0.84 fm) but the phase structure differs.

### II. Theory: the quark ontology

1. **State the substrate's answer** to (A) vs. (B) vs. (C). Argue from the programme's principles (one law, no background, no imported species, emergent from fabric).

2. **The "three-ness" derivation:** if the answer is (B) (emergent attributes), derive WHERE the three comes from. Candidates:
   - The three sectors (field/matter/space) — but these are MODES, not sub-particles.
   - The three lowest winding modes of a bound ring (m=0, ±1) — the proton has spin 1/2 (two states) and the third might be a radial mode.
   - The beat-frequency decomposition: the proton's internal cycle has three phase relationships (field-register, register-matter, field-matter) — the "three quarks" are these three phase channels.
   - A topological constraint: the bound state's winding number on its support graph gives a three-fold degeneracy (the SU(3) analogy, but emergent not imposed).

3. **The falsifiable prediction:** under "deep inelastic scattering" (a high-energy probe that breaks the bound state), what emerges? If (A), three pre-formed sub-objects fly apart. If (B), the energy distributes into the three lowest decomposition modes of the shattered object. The SUBSTRATE predicts one of these — derive which.

4. **The CLAUDE.md prohibition:** the "quark road" (separately-ledgered, identity-carrying, separately-clocked fiber components) is rejected. But emergent attributes of the WHOLE are not the quark road. Clarify the boundary: what IS allowed (mode structure of a bound state) vs. what is NOT (separate particles with their own ledgers).

### III. Numerical simulation method

1. **The size-selection experiment:** seed objects of different sizes (N=3, 6, 12, 24, 48 cells) in a uniform bath. Run each to T=3000. Measure: does each converge to the same N* (the ground state), or do they maintain their seed size? The ground state is where ALL initial sizes converge.

2. **The binding energy curve:** for each N, measure the equilibrium Em (the metabolic park level) and the leak rate. Plot E(N) and find N*.

3. **The quark-emergence experiment:** create a stable bound state at N*. Then "smash" it (inject a high-energy pulse). Measure the decomposition: how many fragments? What phase structure? Does the "three-ness" appear?

### IV. Acceptance bars (pre-registered)

- **S1 (size selection):** objects seeded at N ∈ {3, 6, 12, 24, 48} all converge to the same N* at T=3000. N* emerges from law constants.
- **S2 (meta-stability):** objects at N ≠ N* either grow toward N* (if too small) or evaporate (if surface-dominated). No stable size other than N*.
- **S3 (J-metric tuning):** N* shifts with j_metric (the surface modulation changes the balance point).
- **Q1 (quark emergence):** a bound state at N*, when smashed, decomposes into a PREDICTED number of fragments with a PREDICTED phase structure. The "three-ness" (or whatever the substrate predicts) emerges from the decomposition, not from pre-existing components.
- **Conservation, anti-ignition, invariant surface:** as always.

### V. What NOT to do

- Do NOT impose a size (no "nv=6" seed forcing).
- Do NOT impose the three-ness (no "three-component fiber").
- Do NOT import QCD or the Standard Model.
- Do NOT violate the CLAUDE.md prohibitions (no background, no separate ledgers, no quark road).

---

## CONTEXT FROM v89–v95

Read these in order:
1. `v89/PRINCIPLE.md` — the one law and four constraints.
2. `v93/RECKONING.md` — Face C transients are not matter.
3. `v94/consult/DESIGN_claude_face3.md` — the unification theorem (fiber minimality).
4. `v94/consult/NOCANCELLATION_claude.md` — the no-cancellation primitive.
5. `v94/FINAL_STATUS.md` — the v94 scoreboard.
6. `v95/consult/DESIGN_claude.md` — the J-metric design.
7. `v95/kernel/freecell.c` — the kernel (all knobs: door_g, door_avg, chi_load, j_metric).

The programme's measured record:
- DS-G3: r(Em,Ee)=+0.684 (coherent conversion confirmed).
- J1: R2=0.774 (best phase retention, via J-metric).
- M6: metabolic park approaching (chi_load=1, 10× slope reduction).
- M1: door pair confines (Em=17 vs legacy 417).
- The door pair's stimulated outtake is the single pathology blocking compositions.

The J-metric is an amplifier, not a generator. The metabolic cycle is a candidate generator. The SIZE SELECTION is the thermodynamic/structural question that makes the ground state a proton, not a meta-stable hoard.
