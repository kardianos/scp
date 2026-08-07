I'll ground this review in the project's actual law, free-cell substrate, and failed amplitude/feedback attempts so the critique is mechanism-level rather than metaphorical.# Critical design review: 2-component dense amplitude (uptake/outtake)

Grounding: S2 design memo, Phase-M map, L-1 feedback closure (3 forms + Path-A review), π-flatness / open-books programme, and the field hop algebra already in `freecell.c`. Verdict first, then the four asked aspects.

---

## Executive verdict

**Two equal planes `(m1, m2)` with U(1) rotations *like the field* is the right *algebra* for dense transport — but “uptake/outtake” as labels on those planes is mostly storytelling unless the *couplings* are asymmetric.**

What failed is not “we forgot complex numbers.” What failed is:

1. **Transport primitive:** dense energy moves by real debit/credit (`Em[src]−=…`, `Em[rcv]+=…`) while phase only gates. Feedback cannot invent a first moment the ledger never carries.  
2. **Composition-as-energy:** any design that sets energy from `|Σ A|` collides with `|ΣA|² ≠ Σ|A|²`.  
3. **Open books are spatial:** gravity needs *separated* intake and exhaust through the medium (π-flatness theorem), not two co-located numbers.

So: **complexify Em the field’s way (pairwise norm-preserving hops)** is real physics. **Uptake/outtake, projective line, and DEC-as-cell-interior** are either notation or point at a *different* design (asymmetric door dual / edge forms / multi-cell rings).

---

## 1. Uptake and outtake as the two components

### Is it sound?

**As phenomenology: yes.** The programme already measures and names this:

- Door economy: condensation vs evaporation; radiance fixed point \(\hat x^* \approx 0.62\) is *intake = outtake* in energy books.  
- Open books: “a chain that eats at one end and shines at the other” is the gravity acceptance in the v92 charter.  
- FLOW: only assertion apparatus with open books digs beds; closed-book matter is bed-free (π-flat).

**As a map from two *cell-local* numbers → “intake face / exhaust face”: mostly not.**

A free cell is a point (plus radius). `(m1, m2)` lives at one roster index. “Movement = throughput from uptake to outtake” needs either:

- **space** (different cells / ends of a chain), or  
- **oriented edges** (link-directed flux), or  
- **asymmetric process ports** (different laws for the two components).

Two co-located equal planes do not automatically define faces. They define a **phase and a magnitude** — which the dense sector already has as `(Em, th2)`.

### What could uptake/outtake *actually* be?

| Candidate | Already exists? | Matches cell-local `(m1,m2)`? |
|---|---|---|
| Field↔matter doors (F→D / D→F) | Yes, pass 6 | No — events, not sub-state of Em |
| `s_pull` / backsplash (space books) | Yes, door routing | No — π bookkeeping, not transport planes |
| Forward/back want on a slot | Yes, `swant[2s+dir]` | **Yes as link duality**, not cell duality |
| Phase-M shadow per slot-dir | Yes, `sre_/sim_` | Same — **edge-borne**, not cell-borne |
| Asymmetric dual: “in-mode” couples only to condensation, “out-mode” only to emission, with an internal transfer in→out | **Not built** | **Yes — this is the only cell-local reading that is not relabeling** |

The honest geometric match for “uptake / outtake of a parcel” in *this* kernel is already **oriented slot flux**, not two halves of Em. Phase M put the complex shadow exactly there for a reason: directionality grain and chart order live on links.

### Does the framing change the math?

**Re/Im vs uptake/outtake as *names* for equal planes: interpretation only.**  
Dynamics `U(1)` on `(m1,m2)` with `Em = m1²+m2²` is isomorphic to \(\psi = \sqrt{\mathrm{Em}}\,e^{i\,\mathrm{th2}}\).

**It becomes new math only if the planes are *not* dual under the same group action**, e.g.:

- different door couplings (in vs out),  
- different mobility / chart roles,  
- irreversible atom routing that treats them unequally,  
- or they are not co-located (dipole / two poles).

Without that, “uptake/outtake” is a story painted on complex Em.

**Critical flag:** the v92 gravity sentence (“separated intake and exhaust”) is a **multi-site open-book** design. Collapsing it into two components of one cell is a category error that will not break π-flatness (Wall 3 is about *event routing proportionality*, not missing an Im part of Em).

---

## 2. Decay rates and translation without stable uptake–outtake modes

### Is the diagnosis coherent?

**Partly. The chicken-and-egg is real; the “particle that can’t eat/shine decays” slogan overfits.**

**Solid, measured chain (L1 addendum 3 / S2):**

- Fifth dies → no transport on those links → no deposits → no shadow → feedback cannot revive.  
- Same root for e3b: magnitude want has no coherent phase cargo; feedback amplifies residual noise.  
- Foam phase memory \(\tau_\phi \sim 3\) t.u. for flying packets; door can imprint phase (qp_phase 3.5σ) but **retention is zero** if phase lives only on naked `th2`.

That is: **phase structure dies; energy books can still be full.** “Decay” here is dephasing / channel death, not Em → 0.

**Weaker part of the diagnosis:**

- Many dense cells *do* metabolize (radiance drive, door events). They “eat and shine” in the scalar sense and still do not translate (B2, v_COE ceilings).  
- Radiating equilibrium at cap is **scalar book balance**, not proof of missing sub-modes.  
- The dead fifth is a **chart-order / gate / circulation** failure on *links*, not “this cell lacks m2.”

So: “no stable separated modes ⇒ no coherent current” is right if “modes” means **a protected, directed phase current (edge- or dual-borne)**. It is wrong if it means “each cell needs two named faces or it evaporates.”

### How would a decay rate emerge from a 2-component cell?

Distinguish carefully:

1. **If transport is unitary on `(m1,m2)` (field-like hops)**  
   - `|ψ|²` does **not** decay under hops (norm-preserving).  
   - What can “decay” is **relative phase** between neighbors, chart frames, or dual ports → directed current \(J \sim \mathrm{Im}(\psi_i^*\psi_j)\) dies while Em remains.  
   - That matches: blob can sit as a hot mass; structure/current dies; fifth gg → 0 with mass piled on the wrong voice.

2. **If you keep real wants + phase cargo**  
   - “Decay rate” is just the old dephasing of `th2` / shadow window — already measured; not new.

3. **If uptake/outtake are asymmetric duals with internal throughput**  
   - A *rate* appears as imbalance: \(r_\mathrm{in} - r_\mathrm{out}\).  
   - Stationary particle: \(r_\mathrm{in} = r_\mathrm{out}\) (radiance fixed point language, now *internal*).  
   - Moving particle: spatial offset of in- vs out-dominant regions → net medium current (open books).  
   - “Sub-atomic decay” ≈ loss of the **phase lock that keeps in→out coherent** → throughput randomizes → looks like burn-down / no translation.  
   - This is the only version that ties **decay rate**, **radiating equilibrium**, and **dead fifth** into one mechanism without handwave.

**Chicken-and-egg under native unitary hops:** a dead fifth still has a problem, but a different one. Feedback needs deposits; unitary hops need **finite hop angle on those edges**. If hop weight is gated by the same collapse that kills gg, you still get death. You need either:

- hop weights not solely owned by the dying gate, or  
- a maintenance bar framed as *sustain live fifths* (Path-A review was right: sustain ≠ ignite).

---

## 3. Geometry — three framings

### (a) Projective line `[m1 : m2]`

**Assessment: almost pure notation.**

- Phase \(\arctan(m2/m1)\) + scale \(Em = m1^2+m2^2\) **is** polar form of what you already store as `(Em, th2)`.  
- “Care about ratio, not magnitude” = “phase is shape, Em is energy” — already the sector split.  
- Projectivizing does **not** fix conservation: conservation lives in the affine scale; the RP¹ piece carries no energy. Path A’s failure mode was exactly “energy ledger real + phase informational.” Projective language *encourages* that split.

**Useful only as:** a reminder that identity/coherence lives in relative phase (chart-folded), not in Em. Not a design of transport.

**Verdict: notation.** Do not build a campaign around projectivizing matter.

---

### (b) Toroidal geometry (two cycles)

**Assessment: strong for *composites and spin*, weak as *cell-local 2-component*.**

What the project actually measured:

- Winding can be injected in the field; door is magnitude-blind unless instrumented; matter does not retain winding (QUENCH).  
- Ring-lock needs winding-sector bookkeeping (CANTUS).  
- Flux-machine interior needs **persistent internal circulation at zero net transport** (S2 item 4).

A \(T^2\) with (translation-like phase, spin-like phase) is a natural **configuration space for a multi-cell ring / chord machine**, not for one cell’s Em. “Uptake cycle / outtake cycle” on a torus is almost the same as **two independent U(1) currents on a physical loop** — which is L2 territory (charts + circulation), not “split Em into m1,m2.”

**Risk:** permanent dual windings as fixed cell attributes reeks of the rejected multi-component lattice (background substructure). Safe form: **emergent** windings of dynamical rings (roster of cells + links), metered as \(W\), not hard-coded dual fields per cell.

**Verdict: best as multi-cell bookkeeping / L2–spin acceptance language; not as the definition of a 2-component dense cell.**

---

### (c) Discrete exterior calculus

**Assessment: excellent language for the *slot layer*; wrong motivation for *splitting Em*.**

What DEC buys:

- 0-cochains ≈ cell scalars (Em, Es, …)  
- 1-forms ≈ oriented slot fluxes (`swant`, `slem`, Phase-M `sre_/sim_`)  
- \(d\) / \(\delta\) ≈ discrete Stokes → telescoping conservation on the *current* complex  
- Currents and “uptake − outtake” as **coboundary balance at a node** — already what deposit pass does

What DEC does **not** buy:

- A background-free dynamical roster is a **time-varying complex**. Instantaneous Stokes still holds; you must rebind forms when cells/links appear or die (free-cell already does bookkeeping surgery). DEC does not remove that cost.  
- Interpreting the two components of a cell as “\(d\) and \(\delta\)” is category-wrong: \(d\) maps form degree, it is not a pair of 0-forms.  
- Hodge star needs a dual mesh / metric — free-cell metric is emergent packing; you can *speak* DEC there, but implementing full DEC machinery is a research tax with little new physics beyond “put amplitude on oriented edges.”

**Verdict: use DEC as vocabulary for “amplitude lives on 1-chains; conservation is Stokes.” Do not use it to justify a 2-plane cell interior.** Phase M already started the right DEC move.

---

## 4. The big question: what is more than “just complexify Em”?

### What “just complexify” usually means (and why it fails)

\[
\psi_i = \sqrt{\mathrm{Em}_i}\, e^{i\,\mathrm{th2}_i},\quad
R_j = \sum_k \sqrt{w_{kj}}\, e^{i\phi_{kj}},\quad
\text{drive wants from }|R|\text{ or }\rho_\mathrm{coh}.
\]

That is **Path A / L-1 form-3** in different clothes:

- energy still moves as real \(w\);  
- phase is cargo or gate garnish;  
- composition either does nothing useful (feedback) or **creates/destroys energy** if you promote `|R|²` to Em;  
- chicken-and-egg when \(w\approx 0\).

Grok’s Path-A review already scored this: not S2-native until the **in-flight object is one complex quantity that energy and phase both ride**, with transport constitutive of composition — or, more sharply, **until the dense hop algebra itself is unitary.**

### What *is* genuinely new (strongest framing)

**Mirror the field’s transport algebra, not only its state shape.**

Field (pass F, already law):

1. Local precession: 2×2 rotation of `(fa1, fa2)` — preserves `Ee`.  
2. Pairwise Givens on each link: mixes `(fa_i, fa_j)` so **`Ee_i + Ee_j` is exact**.  
3. First moment emerges as energy the rotation moved (P1 site 5).  
4. Interference is free; no `|ΣA|²` ledger.

**Dense completion of the same shape:**

- State: `(m1, m2)` per cell with `Em = m1² + m2²` (or keep Em and store a protected phase on an account — but then you must not call it native).  
- **Transport primitive:** pairwise rotations (or Givens-like) on dense planes between neighbors, with angles from chart order, closure, mobility, headroom — **not** `Em−=` / `Em+=` of a real want as the primary mover.  
- Conservation: automatic at FP, same class as field.  
- Translation current: first moment of energy moved by the hop (dense analog of P1 site 5) — **this is “translation IS the current” as mechanism, not slogan.**  
- Anti-ignition: hop *angles* must be zero-sum / balanced in the bath the way FLOW beds and L-1 bias were starved by isotropy — **structural**, not “gate on ρ_coh” (bath ρ_coh ≈ 0.77; that path is ignition-shaped).  
- Retention: phase rides the dense planes (or slot-borne rotated state), not naked `th2` churn.  
- Door remains special: atoms, cap, workfn, radiance stay on the **conversion face**; unitary hops are the **within-mode dense transport face**. Do not unitarize the door into a continuous rotation or you lose ε-atoms and XSEC/PAULI structure.

That design is **not a relabeling**. It changes the dense sector from an irreversible real transport PDE-on-graph to a **unitary (or unitary-plus-door) channel**, which is exactly the field’s proven success mode.

### Where uptake/outtake earns its keep (second-strongest, optional face)

Keep unitary (or near-unitary) dense hops for **within-mode** motion, and add an **asymmetric dual only if open books are the target**:

- `m_in` preferentially fed by F→D, `m_out` preferentially drained by D→F;  
- internal rotation/transfer `m_in → m_out` is the particle’s metabolism current;  
- stationary: balanced throughput, \(\hat x^*\)-class equilibrium;  
- lawful gravity: *spatial* polarization of in- vs out-dominance along a chain so medium π sees net flow — **breaks π-flatness by breaking the proportional closed-book routing**, not by “having an Im part.”

That is **new physics** and matches the blackhole/flow programme. It is **not** “Re and Im.” Calling it uptake/outtake is then justified.

### How framings score on conservation / anti-ignition / retention

| Framing | Conservation | Anti-ignition | Retention | New physics? |
|---|---|---|---|---|
| Projective line | Avoids problem by not moving energy in phase | N/A | Phase-only | **No — notation** |
| Feedback / compose-then-bias | OK if zero-sum; useless for L1-A | Measured OK | Shadow helps; not enough | **Closed experimentally** |
| Compose-to-energy `|R|²` | **Broken** (`|ΣA|²≠Σ|A|²`) | Dangerous | — | **No-go** |
| Cell-local Re/Im = uptake/outtake (equal) | OK if hops unitary | Depends on hop weights | Possible | **Only if hops change** — labels still empty |
| **Unitary dense hops (field-symmetric)** | **Best** (pairwise norm) | Must design zero-sum / isotropic starve | **Best match to field** | **Yes — strongest** |
| **Asymmetric dual + spatial open books** | OK if dual transfer conserves Em | Dual must not reward bath unison doors | Internal cycle can outlive `th2` | **Yes — gravity face** |
| Torus as multi-cell windings | Bookkeeping | — | Spin/L2 | **Yes for L2, not L1 cell** |
| DEC | Stokes on links | — | Edge amplitude | **Language for slots, not cell split** |

---

## Concrete recommendation (mechanism, not analogy)

1. **Do not** implement “ψ_matter = m1 + i m2” as a rename of `(Em, th2)` with composition feedback. That path is closed.  

2. **Do** treat L1-A as: **dense within-mode transport becomes pairwise plane rotations** (field pass F’s cousin), with chart-folded angles and structural anti-ignition on hop weights. Acceptance stays e3b translation + bath byte-stable + drift floor.  

3. **Use “uptake/outtake” only for:**  
   - oriented **link** currents (already real), and/or  
   - a later **asymmetric dual + door polarity** face aimed at open books / beds — after translation exists.  

4. **Park** projective language; **use** torus language for composite circulation / spin retention (L2, QUENCH rerun); **use** DEC language for “amplitude on oriented slots,” which Phase M already partially is.  

5. **Historical caution:** any permanent multi-component scaffold per cell that does not emerge from energy and contacts is a background. Unitary dense planes are acceptable only as **modes of the same free cell** (like fa1/fa2 already are for the field sector — the programme already has two field planes without calling that a lattice of quarks). The red line is **fixed multi-site substructure** or imported species, not “a second plane of an existing mode.”

### Strongest framing in one line

**Field-symmetric unitary dense hops** make the 2-component cell real physics; **uptake/outtake** is either edge-oriented current (already the right place) or a later asymmetric dual for open books — not the default meaning of Re/Im.

Anything that keeps real `swant` as the energy mover and only dresses it with phase remains feedback, however elegantly narrated.
