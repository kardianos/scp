# Metabolic matter — temperature, size, and why condensed mass still melts

**Status:** analytical / symbolic only. No kernel changes. Segregated under
`v93/theory/`. Written 2026-08-08 against the Face A/B/C + RING_DNLS record
and the long-standing MASS/PRINCIPLE framing of intake–outtake balance.

**Purpose:** name the structural concept the programme is circling but has
not yet closed: **stable matter is not a sealed lump; it is a continuous
intake/outtake cycle whose rate is temperature, whose fixed point is mass,
and whose far-field signature is the imbalance (or the cycling) of space
and field through that cycle.** The recent findings (hoards condense, then
radiate toward Em\* or resist and sit frozen) are re-read as partial
realizations of one side of that cycle — outtake without matched intake.

---

## 0. Where we stand (the melt that will not stop)

### 0.1 Chronology of attempts (last few histories only)

| era | what was tried | what it bought | what still melted |
|---|---|---|---|
| **v89 MASS** | sealing, pressure annealing, rings, flux balance `dM/dt = I−O` | doctrine: *stability ≠ leak=0*; needs continuous cycle | every free blob under V2g: `I < O` everywhere above the skirt |
| **v91 RADIANCE** | graded outtake `k_rad·cap·x^{p_rad}` (`p_rad=4`) | interior fixed point of *tax* (x\*, Em\*) instead of cap-wall only | tax alone: drained ring Em 0.81→0.48; no refill engine |
| **v91 WORKFN** | condensation demand lives in bound matter, not empty space | door does not eat vacuum; optics-grade far demand | dark bath ⇒ door never fires (Face B dark control) |
| **v92 Phase-M / L1 feedback** | shadow amplitude, three feedback drivers | proved the additive Em-ledger *rejects* interference | coherent dense transport never arrived via feedback |
| **v93 pass U** | unitary dense hops (Givens); linearize; hop_order Strang | **transport works** (fd seed-robust +x); conservation theorem | localized winding diffracts; medium dephases ~10 t.u. |
| **v93 RING / DNLS** | closed cycle + q_detune self-trap | vacuum W=+2 holds; energy condenses (qd≥0.6) | **binds energy, not phase**; mobility dies with trap |
| **v93 Face A/B/C** | Strang; lit-bath door; condensation lane | stable *hierarchical hoards* to t=300; radiance truncates spectrum | deep traps resist tax but are not monodisperse particles; door traffic scrambles phase; **no metabolic refill** |

**One-line diagnosis of the present:** we have built a **conservative clumper
(outtake-shaped filter on a self-trapped hoard)** and a **coherent transporter**.
We have **not** built a **closed intake↔outtake machine**. That is why condensed
mass still “melts away” under the adopted law, or freezes into a non-radiating
over-cap lump when detuning is extreme: both are failures of *balance*, not of
clumping.

### 0.2 The measured objects (Face C, re-read)

- **qd=0:** packet melts (diffraction). No trap ⇒ no condensed matter.
- **qd≥0.6, no radiance:** hierarchical hoards form, hold Em_max to t=300.
  Spectrum broad / long-tailed — **not a single exact size**.
- **V3a radiance on mild hoards:** peak walks to Em\*≈1.55 (x\*≈0.62). Radiance
  *selects*, but by **erosion from above**, not by minting a unique species.
- **Deep self-trap + radiance:** Em_max stays ≫ cap; tax runs but the trap
  decouples — a **frozen over-cap object**, not a proton analog.
- **π-flatness / g4 (standing):** closed-book matter loads the medium with
  identically zero far field; only *dying* matter (dM/dt≠0 bookkeeping)
  disturbs space at range.

So: **condensation is real; monodisperse stable matter is not yet.** The gap
is structural, not a knob.

---

## 1. Postulate (user’s core claim, formalized)

### P1 — Existence is process

> All matter exists only as a continuous **intake and outtake** of convertible
> energy (field ↔ dense ↔ space). A sealed Em store with zero throughput is
> not “cold matter”; it is either **absent** (no cycle) or **unstable**
> (will be eroded by any open channel).

This is already programme language (PRINCIPLE: identity is succession of
conversions; MASS user correction 2026-07-28: *“stability does not require
leak = 0; it requires continuous movement and a structure that FORCES field
intake in movement to equal field outtake”*). Face C’s frozen deep trap is
the negative example: mass without a live cycle.

### P2 — Temperature is the speed of the cycle

> **Temperature is not a property of vacuum energy density.** It is the
> **rate** of the matter-bound intake/outtake cycle — how fast the machine
> turns. Only condensed matter can have temperature, because only condensed
> matter can sustain a closed conversion cycle against a continuum of modes.

Symbolic bookkeeping (local, per-particle or per-hoard):

```
M(t)     = dense inventory (Em integrated over the support)
I(t)     = intake rate   (field→dense + space→dense pathways that *enter* M)
O(t)     = outtake rate  (dense→field radiance + dense→space return)
Φ(t)     = throughput    = min(I,O)  when balanced;  the “breathing” rate
T        ∝  Φ / M        (cycle frequency per unit mass)
         or T ∝ Φ^{1/4}  if outtake is Stefan-like in T  (see §2)
```

Balance:

```
dM/dt = I − O.
```

- **Fixed point of mass:** `I(x*) = O(x*)` with restoring slope `I' < O'`
  (excess mass radiates faster than it refills; deficit refills faster than
  it radiates). That fixed point *is* the species mass window.
- **Temperature:** nonzero exactly when Φ > 0 at that fixed point.
  **T = 0 is only asymptotic** — approachable by Φ→0, never a stable resting
  state of matter (MASS doctrine mapping already recorded).

### P3 — Only matter has temperature

Vacuum / pure field / pure space can carry energy and can be *hot* in the
sense of mode occupation — but in this programme’s no-background ontology:

- **Field** has amplitudes and a unitary transport algebra; it does not host
  a *condensed* intake/outtake fixed point.
- **Space** is a convertible mode (Es); pressure pushes; it is not a
  thermometer by itself.
- **Dense matter** is the only sector that can *close a book*: condense at a
  work-function door, radiate through a graded outtake, load-detune its own
  pitch, and thereby select a load x.

Hence: **thermometry is a matter instrument.** A galactic “temperature” of
field density without condensed cycles is a different object (see §6).

---

## 2. Energy transfer above absolute zero: why T⁴

### 2.1 The law already encodes a Stefan-like outtake

Adopted V3a radiance (kernel comment, pass 6; `laws_V3a.cfg`):

```
demand_rad = k_rad · cap · x^{p_rad} · det_clock,
x = Em / cap,
p_rad = 4,   k_rad = 0.05.
```

Per-time outtake rate (with `rad_clock=0` compensating beat stretch):

```
O_rad(x) ∝ x^4.
```

This is **not** an accident of notation. In continuum thermodynamics the
Stefan–Boltzmann law is

```
j = σ T⁴
```

for the energy flux of a black body. The structural analogy the programme
needs is:

| continuum thermo | fabric metabolic picture |
|---|---|
| black body | condensed dense support with open conversion doors |
| T | cycle rate / load proxy (see below) |
| σ T⁴ | radiance channel O_rad |
| absorption cross-section · u_rad | intake I from ambient field |
| detailed balance | I(x,T_bath) = O(x,T_self) at equilibrium |

So **“energy transfer in anything above absolute zero by T⁴”** is, in-fabric:

> any live condensed object with Φ>0 *must* exchange energy with its
> environment at a rate that scales as the fourth power of the load (or of
> the temperature assigned to that load), because that is how the outtake
> channel is written — and a matched intake must exist or M walks.

### 2.2 Relating load x to temperature T (symbolic)

Two consistent assignments (not yet decided by measurement):

**(A) Load-is-temperature (simplest).**
Identify `x ∝ T` (or `x = T/T_ref` in dimensionless form). Then O ∝ T⁴
*identically* with p_rad=4. The radiance fixed point Em\* = x\*·cap becomes
a preferred **operating temperature** of mild hoards — Face C’s walk to
Em\*≈1.55 is then “thermal equilibration of load,” not particle mass.

**(B) Cycle-rate-is-temperature (metabolic, preferred by P2).**
Let the local dense pitch be `ω₂(x) = w₂ / (1 + q·x)` (load-detuned). Cycle
throughput Φ is some function of door traffic, radiance clicks, and
coherent recapture. Define

```
T  :=  ħ_eff · Ω_cycle,     Ω_cycle = Φ / (action per cycle),
```

with action atoms `ε = A₀ ω / 2π` already in the door. Then Stefan form
requires the *macroscopic* outtake to satisfy

```
O = κ T⁴
```

while the *microscopic* law is still `O ∝ x⁴`. Closing the dictionary:

```
x = x(T)  such that  k_rad cap x⁴ = κ T⁴
  ⇒  x ∝ T    (same as A, if κ is load-independent)
```

or, if κ itself depends on headroom / geometry / Es depression,

```
x = f(T; geometry),
```

and monodispersity becomes a **geometry-selected** fixed point of
`I(x, geometry) = κ[T(x)]⁴` — the missing structural concept of §5.

**Either way:** p_rad=4 is the programme’s already-adopted T⁴ skeleton.
What is missing is not the power; it is **the intake side of the ledger
and the structure that forces I=O at a unique x\***.

### 2.3 Why T=0 is unreachable for matter

If T∝Φ/M and matter *is* the cycle, then:

```
T → 0  ⇔  Φ → 0  ⇔  I,O → 0.
```

But with any open radiance channel O_rad(x) = c x⁴ > 0 for x>0, holding
Φ=0 requires either:

1. **x=0** (no matter), or
2. **I structured to cancel O exactly** at some x\*>0 with Φ = I = O > 0
   (finite T), or
3. **decoupling** that shuts all channels (deep q_detune freeze) — which is
   not T=0 matter, it is a **metastable sealed store** outside the metabolic
   definition (Face C deep corner).

Thus: **protons that “don’t melt” are not cold sealed bags; they are warm
machines at a restoring I=O fixed point.** Absolute zero is the absence of
the machine.

---

## 3. Why protons and neutrons have one exact size

### 3.1 What Face C shows we do *not* yet have

Radiance truncates a **broad hierarchy** at Em\*. That is a **soft cutoff**,
not a **sharp mass eigenvalue**. Real nucleons have:

- a single charge radius (up to isospin),
- a discrete spectrum of excitations,
- no continuous family of stable “half-protons.”

So monodispersity is **not** automatic from `O ∝ x⁴`. You need a fixed
point that is:

1. **isolated** (not a valley segment — MASS M-R1 failure mode),
2. **restoring in mass** (`∂(I−O)/∂M < 0`),
3. **selected by dimensionless ratios only** (no-background forbids an
   absolute imported scale — RECONCILIATION),
4. **stable under isolation** (R2 / PATH P2: one hoard in vacuum).

### 3.2 Symbolic condition for a unique size

Let state variables of a candidate particle be

```
x = Em/cap,     s = local space contrast (Es depression),
γ = phase-order / winding support,   geometry encoded in contacts.
```

Rates:

```
I = I(x, s, γ; bath),
O = O_rad(x) + O_evap(x,cap) + O_space(x,s) + …
  = k_rad cap x⁴ + … .
```

A **species** is an isolated root

```
I(x*, s*, γ*) = O(x*, s*, γ*)
```

with negative Jacobian eigenvalue along the mass direction and a basin of
attraction from a range of seeds. **Exact size** means the root’s mass

```
M* = ∫ Em* = x* · cap · N_eff(geometry*)
```

is discrete because **geometry\* is discrete** (topological sector, shell
count, winding m, contact graph) — not because Em is continuous but
“somehow prefers 1.55.”

Face C’s hierarchy is what you get when **geometry is not quantized**: many
clusters, many N_eff, one soft tax ceiling.

### 3.3 Electron orbits and stability (reframe without importing Bohr)

Standard language: “electrons in orbits keep atoms stable; bare nuclei
would be different objects.” In fabric language (PATH_claude shells-without-
charge path, RECONCILIATION §2):

1. **The dense hoard depresses Es** — a graded well (space is a mode).
2. **Field modes refract / bind in that well** — “orbits” are **standing
   field (or light-sector) patterns** on the depression, not point particles
   on rails.
3. **Cap refusal (PAULI-0)** — overflow dense energy is refused or evaporated;
   an occupation limit per well is an exclusion *seed* without fermions yet.
4. **Stability of the atom** = mutual completion of books:
   - nucleus cycle: I_n = O_n at M_n\*,
   - shell cycle: field pattern continuously converts against the well,
   - the two couple so neither can run down without the other.

**What electron orbits have to do with stability:** they are not decorative
quantization; they are **the intake/outtake completion of the outer book** —
the part of the machine that recaptures and re-emits at the soft scale,
while the nuclear book runs the hard scale. An isolated bare hoard with only
O_rad and no resonant recapture **must** melt or freeze (measured). An
atom-like composite is the minimal structure that can hold **both** books
open and balanced.

This does **not** require importing U(1) charge yet. The first shell result
the path names is **charge-free**: well + field modes + cap refusal.

### 3.4 Why a proton does not melt (target statement)

A proton analog, in this theory, is a configuration such that:

```
I_struct(x*, s*, γ*) = O_rad(x*) + O_other(x*, …) > 0,
dM/dt|_* = 0,   ∂(I−O)/∂M|_* < 0,
γ* topologically or nonlinearly protected against bath dephasing,
mobility compatible with binding (the present antagonism resolved).
```

**It radiates continuously** (O>0, T>0) **and does not evaporate** because
**every outtake quantum is matched by intake** — not because radiation is
off. Cosmic “cold” matter is low-T (small Φ/M), not zero-cycle.

The Face C deep trap fails clause `I_struct` (decoupled, Φ≈0 channels).
The mild radiance-truncated hoard fails clause `I_struct` in the other
direction (O runs, I does not keep up ⇒ walk to Em\* then continued slow
death if I remains under-supplied — the V2g no-particle theorem).

---

## 4. Why matter in space keeps emitting without evaporating

### 4.1 The paradox in ordinary language

Interstellar dust, stars, planets: continuous electromagnetic emission, yet
they persist for geological / stellar times. Naive sealed-bag models predict
either (a) radiation ⇒ mass loss until gone, or (b) stop radiating ⇒ T=0.

### 4.2 Resolution in metabolic language

```
Emit continuously  ⇔  O > 0  ⇔  T > 0.
Not evaporate      ⇔  I = O  (or I ≈ O with slow secular drift).
```

Sources of **I** for real matter (fabric analogs in parentheses):

| continuum process | fabric analog |
|---|---|
| absorption of CMB / starlight | condensation door from ambient Ee (workfn-gated) |
| nuclear / chemical power | internal conversion cycles (dense internal process) |
| gravitational / compressional work | Es pressure / piston intake (MASS R1) |
| stimulated recapture of own radiation | resonant recondensation of self-shed field (MASS “choir at the door”) |

**g4 / π-flatness constraint:** a *stable* object’s *space* book may close
locally (no steady monopole from stable mass — measured). That does **not**
forbid a steady **field** book (photon-analog radiance) if intake matches
outtake in the dense↔field door. The far-field gravity problem and the
thermal radiation problem are **different books**:

- **Space book (G):** stable matter may be π-flat; dying matter shows
  dM/dt wind.
- **Field book (P/T):** stable matter *must* run O_field = I_field at
  finite T, or it is not warm matter.

Stars are extreme: internal I ≫ environmental I, high T, O∝T⁴ huge, mass
changes on nuclear timescales. Planets/dust: low T, environmental I
dominates, slow balance. Protons: internal cycle at nuclear scale, T small
in macroscopic units, lifetime ≫ age of universe if |I−O|/M is tiny.

### 4.3 What the sim is missing (honest)

Present V3a + unitary channel implements **O_rad** and **cap evaporation**
and a **door** that can condense when the bath is lit — but:

1. **Dark bath ⇒ I=0** (Face B): workfn keeps vacuum from condensing;
   no metabolic refill.
2. **Lit random-phase bath ⇒ I phase-hostile** (Face B): door traffic
   *destroys* winding; amp_door cannot fix random-phase injection.
3. **No structural I(x)** that *forces* recapture of self-radiation
   (stimulated, phase-matched) — the standing MASS proposal, unbuilt.
4. **Deep trap shuts O and I together** — sealed bag, not proton.

So continuous emission without melt requires **coherent or structurally
gated intake**, not merely more noise_amp.

---

## 5. The core missing structural concept

### 5.1 What is *not* missing (already present)

| piece | status |
|---|---|
| Outtake ∝ x⁴ (T⁴ skeleton) | adopted (radiance) |
| Cap / workfn doors | adopted |
| Conservative dense transport | pass U (when linearized) |
| Energy condensation / self-trap | q_detune DNLS (binds energy) |
| Topological phase support (vacuum) | C6 ring + Strang |
| Doctrine I−O balance | MASS, measured no-particle under V2g |

### 5.2 What *is* missing (the hinge)

**A self-consistent metabolic fixed point:** a structure that

1. **hosts a nonzero throughput Φ\*** (T\*>0),
2. **forces I = O at an isolated x\*** (unique mass / size),
3. **protects the phase/order parameter γ** that makes I *structural*
   rather than random-bath (Face B: random I kills order),
4. **remains mobile or composite-capable** (present transport/binding
   antagonism is a wall — RING_DNLS synthesis #3),
5. **closes or cycles the space book** in a way compatible with π-flatness
   *and* galactic phenomenology (§6).

Name for the concept (working):

> **The metabolic attractor** — matter as a **limit cycle** (or standing
> conversion cycle) in (Em, Ee, Es, phase) with finite flux Φ, not as a
> static Em lump.

Uptake/outtake are **moments of the current and the doors** (v93 README
§II.6) — not two state components. The missing piece is the **geometry +
phase lock that makes those moments self-restoring**, not another float on
the cell.

### 5.3 Why recent routes almost-but-not

| route | almost | not |
|---|---|---|
| Radiance fixed point | selects Em\* under tax | tax without matched I ⇒ erosion, not life |
| DNLS hoard | holds energy | phase-blind; freezes or hierarchy; no Φ\* |
| Closed ring | holds γ in vacuum | medium dephases; no I engine |
| amp_door | correct merge algebra for coherent field | thermal field is random ⇒ same as scramble |
| Unitary transport | moves M coherently | does not bind γ; not a cycle by itself |

**The core gap is not transport algebra (v93 closed that gap partially).
The core gap is the closed metabolic loop with phase-safe intake.**

### 5.4 Symbolic minimal model (for later measurement, not code)

Let a single degree of freedom x(t) and an order parameter r∈[0,1]
(coherence / winding fidelity):

```
ẋ = I₀ · r · B(x; bath) − k x⁴ − λ (x − x_cap)_+
ṙ = −Γ_bath · r + G(x) · r · (1−r) − D_door · (1−r_field)
```

- `I₀ r B` — intake only when order is present and bath/structure admits
  (workfn, resonance).
- `k x⁴` — radiance outtake.
- `Γ_bath` — contact dephasing (~10 t.u. measured scale).
- `G(x)` — nonlinear self-phase regeneration (missing; DNLS did energy not
  phase).
- `D_door` — door traffic cost unless field is phase-matched (`r_field≈1`).

Metabolic matter = attracting fixed point `(x*, r*)` with `r*>0`, `x*>0`,
`Φ* = k(x*)⁴ > 0`. Face C deep trap ≈ large x, r≈0, Φ≈0. Face C melt
control ≈ x→0. V3a mild hoard ≈ x→x_rad, r≈0, slow leak if I₀ r B under-
supplied.

---

## 6. Dark matter as field-density contrast — implications for stable matter

### 6.1 The postulate (user)

> Dark matter is not a particle species. What astronomy attributes to dark
> matter is **field (or space) density structure** that acts as relative
> matter for the laws of motion at galactic scales.

### 6.2 Compatibility with the programme’s measured far-field wall

Standing result: **stable condensed matter is π-flat** — it does not source
a steady space monopole; g4’s 1/r from *stable* mass is absent; only
evaporating/refilling bookkeeping shows at range.

If galaxies show “extra gravity” from **field density gradients** (or Es
gradients) rather than from extra Em particles, then:

1. **Galactic mass budgets** are mis-attributed Em; the true dynamical mass
   includes **convertible mode gradients** (Ee, Es, or their fluxes).
2. **Stable baryonic matter** can remain **locally π-flat** while large-scale
   **field/space texture** supplies the rotation-curve phenomenology —
   no contradiction with pad-11/g4 if the texture is *not* a stable Em
   monopole.
3. **“Relative matter”** means: equations of motion (conversion-driven
   translation of test masses) respond to **gradients of the fabric’s energy
   partition**, not only to integrated Em.

Symbolic schematic:

```
a_test ~ ∇( effective potential of {Ee, Es, Em} ),
```

with Em-stable objects contributing negligibly to long-range Es monopole,
while **cosmic web field/space contrast** dominates at kpc–Mpc scales.

### 6.3 What that forces on stable matter (the sharp implications)

If dark matter is **not** a second stable Em species, then:

| implication | content |
|---|---|
| **1. Stable matter is rare and special** | Most of the energy in a galaxy may be *non-condensed* modes. Condensed Em (stars, gas, people) is the **metabolic minority** that can have T and spectra. |
| **2. Stable matter must be *harder* to form, not easier** | Face C’s easy hierarchical hoards are *not* yet “baryons.” Real monodisperse M\* should be a **small measure** of configuration space — a narrow metabolic attractor — consistent with most energy remaining diffuse field/space. |
| **3. No Em halo of slow particles** | Looking for stable Em lumps as DM is the wrong search. Instruments should look for **space/field gradient phenomenology** (π-structure at large scale, anisotropic conversion, lensing-without-Em). |
| **4. Temperature is reserved for the condensed minority** | Diffuse field can have mode temperature in a QFT sense, but **thermodynamic T as cycle rate** belongs to metabolic Em. CMB-like baths are *boundary conditions for I*, not dark matter. |
| **5. Far-field gravity from baryons is weak by construction** | Matches π-flatness: planets/stars don’t need to source strong monopoles if galactic dynamics are mostly texture-driven. The **missing 1/r from stable mass** is then not a bug to brute-force for galaxies — it may be **correct physics**, with galactic “extra mass” living in non-Em modes. |
| **6. Still need *local* metabolic closure for chemistry** | DM-as-texture does **not** excuse missing I=O for protons. Lab-scale matter still needs the metabolic attractor; it only changes what we demand of **long-range Em gravity**. |

### 6.4 Tension to keep honest

- If *all* long-range attraction is field/space texture, then **C5** (space
  throughput at constant mass as 1/r source) may be **reclassified**: from
  “required for gravity” to “required for inertia / equivalence / local
  space cycle,” or delayed.
- If *some* local gravity still needs Em-sourced curvature, then stable
  matter needs a **higher multipole or cycle-averaged** space signature that
  evades π-flat nulls — an open design (pass S / coherent space channel in
  RECONCILIATION long arc).

**Recommendation:** treat DM-as-texture as a **galactic boundary condition**
that *relaxes* pressure to force monopole 1/r out of every proton, while
**keeping** metabolic I=O as the definition of stable *local* matter.
Do not use DM-as-texture to justify sealed Em bags.

---

## 7. Synthesis — one picture

```
                    field bath (Ee)
                         │
              intake I   │   outtake O ∝ T⁴ ~ x⁴
                         ▼
              ┌─────────────────────┐
              │  dense cycle (Em)   │  ← only here is Temperature
              │  Φ = I = O  at M*   │
              │  phase γ protected  │
              └─────────────────────┘
                         │
                    space book (Es)
                 (local cycle; π-flat
                  if stable monopole)
```

**Temperature** = speed of the cycle (Φ per mass), only for condensed matter.

**T⁴ transfer** = the radiance law we already run; needs matched I.

**Electron orbits** = outer book (field patterns on Es well) completing
atomic stability — first probes can be charge-free.

**Unique proton size** = isolated metabolic root (x\*, geometry\*, γ\*), not
a radiance ceiling on a hierarchy.

**Protons don’t melt** because I=O at finite Φ, not because O=0.

**Matter radiates in space without evaporating** for the same reason.

**Dark matter as field density** ⇒ stable Em is the metabolic minority;
galactic dynamics need not be Em halos; π-flat stable matter may be right;
local metabolic closure remains mandatory.

### 7.1 Against the melt: what “close” means

We are close because:

- O∝x⁴ exists,
- condensation exists,
- transport exists,
- doctrine I−O exists,
- topology can hold γ in vacuum.

We lack the **core structural concept**:

> **phase-safe, self-restoring intake locked to outtake at an isolated
> load — the metabolic attractor.**

Until that closes, condensed mass will keep melting, freezing, or
hierarchically radiating toward a soft ceiling — every outcome Face C
already measured.

---

## 8. What this document is *not*

- Not a kernel proposal and not a knob list.
- Not a claim that Face C hoards are nucleons.
- Not a license to import pre-v89 Q-ball shells or U(1) species.
- Not a replacement for the R1–R8 reckoning (RECONCILIATION): those
  experiments decide whether the present hoard is even particulate. This
  note says what **particle** must mean if the overall theory is right.

## 9. Reading order / anchors

| file | role |
|---|---|
| `v93/RESUME.md` | campaign status |
| `v93/FACE_C.md` | condensation + radiance selection (the melt/tax data) |
| `v93/FACE_B.md` | door traffic vs phase-random intake |
| `v93/RING_DNLS.md` | transport/binding/phase antagonism |
| `v93/consult/RECONCILIATION.md` | next reckoning R1–R8 |
| `v89/MASS.md` (§ flux balance, T>0 doctrine) | I−O as stability |
| `v93/laws/laws_V3a.cfg` | k_rad, p_rad=4, workfn |
| `v93/README.md` §II.6 | uptake/outtake as current moments |

---

*End of analytical note. Next work, if desired: (i) map Face C meters onto
(I,O,Φ) time series without new physics; (ii) design a phase-matched intake
probe (coherent driver) as the first empirical attack on §5; (iii) keep R2
isolation as the lynchpin — if isolated hoards dissolve, the metabolic
attractor is not yet in the law.*
