# CELLFAB — the cell-fabric kernel (v89-native)

**Proponent numerical construction.** Subordinate to `PRINCIPLE.md`. First
kernel written inside v89: no prior kernel, no lattice, no SFA, no imported
field. Implements the 2026-07-27 framework direction:

> Fill a field with cells. Each cell has internal bounds; internally it aligns
> on two planes with different harmonic frequencies. Frequency interaction
> gives rise to Bell statistics. Energy transfers at rate C between cells
> through harmonic interaction (resonant joining), but only if the tail phases
> match. Cells are not on a grid: each cell has an area of effect, harmonics
> operate geometrically, and each cell operates through just a list of cells.
> Interaction occurs at cell boundaries according to distance and surface area
> of the adjacent cell — modified by how the harmonic planes overlap. Mass
> doesn't move: cells re-align and transfer energy in stable configurations.
> Cell boundaries (not real, but a useful approximation) prevent collapse to a
> point, because the harmonics will not admit energy that way, and energy
> cannot be destroyed. Energy is conserved across space, dense mass, and EMF
> patterns.

Tags as in `CONSTRUCTION.md`: **[D]** derived from the principle, **[P]**
postulated, **[G]** guess, **[M]** measured (below, by this kernel).

---

## 0. Mapping to the principle

| framework element | kernel object | PRINCIPLE / CONSTRUCTION |
|---|---|---|
| cell | parcel of energy with internal two-plane harmonic structure | parcel `v`, harmonic content `h(v)` |
| cell's space content `E_s` | the cell *being* space | space mode `S` |
| dense content `E_m` | locked pattern energy | dense mode `D` (mass) |
| field content `E_e` | open traveling pattern energy | field mode `F` (EMF) |
| in-flight channel energy | `e_mid` on a link, cycle-gated | transfer mode `T` |
| area-of-effect overlap | channel (link) between two cells | channel `ε` |
| two internal planes, two frequencies | chiral transverse pair (ω₁, n̂₁), (ω₂, n̂₂) | "a link actuates in the two planes transverse to it" (PRINCIPLE §3) |
| resonant joining | frequency-match × plane-overlap × tail-phase gate | resonance condition §2.2, cycle gate R1–R3 |
| rate C | one link crossing per cycle time d/C | c = conversion rate (PRINCIPLE §4.5) |
| mass movement | phase-tilt-gated transfer + re-alignment; cells never move | motion is regular conversion (PRINCIPLE §4.4) |
| collapse prevention | saturation capacity + detuning choke | saturation §2.3(A) |
| curvature | space energy consumed by locking → area-of-effect shrinks → contact defect | §4.3 combinatorial Gauss |

**Honesty note on the approximation.** The cell set and its candidate
adjacency are generated once and persist — that is a scaffold, and under the
strict audit of `CONSTRUCTION.md` §0 a scaffold is a stage. The framework
direction sanctions it explicitly: *cell boundaries are not real, but a useful
approximation*. The defence is: (1) every cell IS energy (its space ledger),
not an empty site — a cell whose energy is fully converted has physically
ceased even though its handle remains; (2) dynamics uses only **relational**
geometry (distances, contact areas, orientations relative to the link) — no
absolute coordinate enters any law; (3) **live** connectivity is dynamical:
a channel exists only while the areas of effect actually overlap, and overlap
depends on the current energy state. Positions are used for initialization
and for *diagnostic reconstruction only*. This is a v89 construction with one
flagged, sanctioned residue — see audit, §8.

**No SFA.** SFA and its viewers are lattice instruments (see CLAUDE.md v89
preamble). Output is self-describing TSV: diagnostics on stdout, optional
per-cell snapshots.

---

## 1. State

### Cell (parcel)

| field | meaning |
|---|---|
| `p` | relational scaffold position (init-time; never enters dynamics except through link geometry) |
| `r` | area-of-effect radius, **dynamical**: `r = r0 · (E_s/e_s0)^(1/3)` [P] |
| `n̂₁, n̂₂` | the two internal plane normals (headless: n̂ ≡ −n̂) |
| `θ₁, θ₂` | harmonic phases of the two planes |
| `ω₁, ω₂` | base frequencies (ω₁ ≠ ω₂ — *different* harmonic frequencies) |
| `β` | beat phase, `β̇ = ω₁ᵉ − ω₂ᵉ` — the frequencies' interaction clock |
| `E_s, E_m, E_e` | energy ledger: space, dense (mass), field (EM) |

Plane 1 (faster ω₁) is the **field plane**: `E_e` transfers through it.
Plane 2 (slower ω₂) is the **dense plane**: `E_m` transfers through it. [P —
the two planes are the cell's two transfer sectors; their beat is the
conversion clock between the sectors.]

### Channel (link)

Two cells are channel candidates iff their areas of effect overlap at init.
Each channel carries, per plane-sector and per direction: in-flight energy
`e_mid` and cycle progress `φ ∈ [0,1)`. In-flight energy is **transfer mode
T**: it is nowhere else, and it is counted in the ledger.

### Effective frequency (load detuning) [P]

```
x = (E_m + E_e)/cap            occupancy
ωᵃᵉ = ωᵃ / (1 + q·x)           a = 1,2
```

Loaded cells run slow. This one rule produces three phenomena at once:
saturation choke (a full cell detunes off resonance from its feeders — *the
harmonics will not admit more energy*), spectral self-sealing of dense blobs
(interior cells share a common detuned frequency, edges mismatch → the
boundary is a failure of resonance, exactly CONSTRUCTION §5.2.3), and
gravitational time dilation (clocks near mass tick slow, §5).

---

## 2. The transfer law (resonant joining)

Per channel {i,j}, per plane-sector a, per direction i→j, the deposit rate is
a product of independent gates:

```
F(i→j) = k_dep · G_geo · G_plane · G_res · G_tail · H_j · E_mobile,i
```

### G_geo — boundary geometry [P, per framework: distance & surface area]

Sphere–sphere lens contact area with **current** radii:

```
a² = [4d²rᵢ² − (d² − rⱼ² + rᵢ²)²] / 4d²      A = πa²   (0 if no overlap)
G_geo = (A/πr0²) · (2r0/d)
```

No overlap → channel dormant (channel death); overlap restored → channel
reopens. Live connectivity is a function of the current energy state.

### G_plane — harmonic overlap along planes [P]

A plane transmits in-plane (best perpendicular to its normal). The two
planes of a cell intersect in a line — **the cell's axis** — and a channel
carries only along directions lying in *both* planes, which collimates
transfer onto the axis the planes jointly define. Per channel a:

```
axis(cell) = [1−(n̂₁·û)²] · [1−(n̂₂·û)²]
G_plane    = axis(i) · axis(j) · (n̂ₐᵢ·n̂ₐⱼ)²        (headless normals)
```

Without the two-plane axis condition (a single plane per channel), oblique
links carry half-open flux in a wide cone and a directed packet diffuses
instead of propagating — measured, see §10 tuning notes.

> **Round 3 (the repair — see `DOUBLESLIT.md` §8).** The field sector no
> longer uses gated transport at all. It is a **two-component signed
> amplitude** (a₁, a₂) on the cell's plane pair (ψ = a₁+i·a₂, E = |ψ|²),
> evolved by exact unitary steps: onsite rotation at ω¹ᵉ (load and wall
> detuning enter here) plus exact pairwise hop rotations over live
> channels with **symmetrically normalized** weights
> ŵᵢⱼ = wᵢⱼ/√(sᵢsⱼ) — a cell's total joining bandwidth is a property of
> the cell, not of its accidental contact geometry (without this the foam
> speckles waves within a few cells). Superposition, diffraction, and
> interference are native; c is emergent group velocity (∝ `field_J`).
> All exchanges with other modes scale |ψ| by the exact energy moved.
> The tail gate, entrainment, flight, and everything below now apply to
> the **dense sector only** — harmony stays gated transport; melody
> became amplitude mechanics.

### G_res — resonant joining [P]

Field sector (historical, rounds 1–2; superseded by the amplitude
dynamics above): a single Lorentzian in the fundamentals,

```
G_res = Γ² / (Γ² + (ω¹ᵉᵢ − ω¹ᵉⱼ)²)
```

Dense sector (round 2, C1 of `CONSONANCE.md` Part VI): **the partial
comb** — sympathetic joining through any coincident partial, best over
coprime ratios p:q with p·q ≤ `comb_limit`:

```
G_res = max over p:q of  [Γₘ/(pq)]² / ([Γₘ/(pq)]² + (q·ωᵢ − p·ωⱼ)²) / (pq)
```

Tongues narrow and weaken with complexity (Tenney height); ties resolve to
the simpler ratio. The winning (p,q) also multiplies the tail-gate and
entrainment phases (Δ = q·θᵢ − q·ωᵢd/C − p·θⱼ), so interval species (the
fifth, the octave) are dynamical objects. Detuned cells do not join; with
load detuning this is the saturation choke and the blob-sealing mechanism —
a **narrow dense width is the rim seal** (spectral boundary of CONSTRUCTION
§5.2; a kinetic freeze, not an attractive channel; see `CONSONANCE.md`).

### G_tail — the tail-phase gate [P, load-bearing]

The sender's harmonic tail, extrapolated across the gap at rate C, must
arrive in phase with the receiver's clock:

```
Δ(i→j) = wrap( θₐᵢ − ωᵃᵉᵢ·d/C − θₐⱼ )
G_tail = [ (1 + cos Δ)/2 ]^p            p = p_gate (default 4)
```

**What this buys — motion without movement [D from the gate]:** write a phase
pattern θ(p) = θ₀ − k·p (downstream clocks lag — the retarded convention).
Then for a link of length d along k̂:

```
Δ(forward)  = d·(k − ωᵉ/C)   → open exactly when k = ωᵉ/C
Δ(backward) = −d·(k + ωᵉ/C)  → strongly suppressed
Δ(static, k=0) = −ωᵉd/C both ways → generically sealed (mod 2π accidents)
```

A static in-phase configuration is *stable* — it does not spread. A phase
**tilt** opens the forward gates only: energy conveys down-tilt at up to one
link length per cycle time d/C, i.e. at up to C. Momentum is not a conserved
fundamental (PRINCIPLE §5); **tilt is its ghost**: a re-alignment of internal
clocks, not a property of moving stuff.

**Frequency window [M].** The backward seal is Δ(backward) ≈ −2ωᵉd/C *after
wrapping*: near ωd/C ≈ π it wraps back to ~0 and the backward gate reopens —
a foam-geometry resonance that melts phase patterns (measured: at ωd ≈ 2.5
the backward gate sat at 0.2 and a tilted packet thermalized into a phase
glass within ~2 t.u.). Keep ωᵉd/C near π/2 for typical link length d: with
d ≈ 1.0–1.45 and ω = 1.5 the backward gate stays ≤ 0.01 across the spread.
Hence the defaults w1 = 1.5, w2 = 0.93.

### H — internal bounds (capacity backstop) [P]

```
H_j = clamp(1 − (E_m,j + E_e,j)/cap, 0, 1)
```

plus a hard cap at delivery. Together with the detuning choke this is why
everything cannot collapse into one cell: the harmonics refuse the energy,
and energy cannot be destroyed, so it stays distributed.

### Cycle-gated flight (transfer takes time) [D: PRINCIPLE §3]

Deposits enter `e_mid`; `φ` advances by `dt·C/d`; **delivery only at φ = 1**
(a partial cycle converts nothing). Delivery is hard-capped by the receiver's
remaining capacity; the blocked remainder stays in flight (stalled
conversion, never destroyed). Latency per hop is exactly d/C → the maximum
pattern advance rate is C by construction; measured front speeds fall short
of C only by foam tortuosity.

Round-2 additions at dense delivery (CONSONANCE Part VI):

* **C2 — dissonance radiates** [P]: the rough fraction of an off-comb
  delivery, R = 2|det|·Γᵣ/(det² + Γᵣ²) · `rough_k` (Plomp–Levelt bump on
  the winning ratio's detuning), lands as field energy instead of dense —
  and returns its space share (s_pull/(1+s_pull) → E_s), so matter leaving
  the dense mode un-converts its space and **curvature tracks conversion in
  both directions**. Ledgered globally as `roughness_radiated`.
* **C3 — harmony is mutual** [P]: dense exchange scales with
  √(mob_src · max(mob_rcv, `mob_floor`·cap)) — it takes two to sound; a
  silent partner answers only at its sympathetic readiness. The field
  sector stays source-driven (melody flows from the singer).

### Resonant joining re-aligns the receiver [P]

Joining begins at the boundary: **every deposit entrains the receiver's
clock** toward the retarded tail (and delivery reinforces it), weighted by
how much is arriving against what the receiver already coherently holds:

```
mix = f / (f + E_mobile,recv + lock_floor·cap)
θₐ,recv += κ_lock · mix · wrap(θ_tail − θₐ,recv)
```

This **regenerates the tilt downstream**: an advancing packet rewrites the
clocks ahead of it to the exact retarded pattern, so motion self-sustains.
(Locking only at delivery — once per d/C — was measured to bottleneck the
front at ~0.1 C; entrainment at deposit time is what makes joining fast
enough to carry a wave.) Heavily loaded receivers re-align less (mix small)
— an inertia-flavored resistance [G, noted not claimed]. The same channel is
the kernel's decoherence mechanism: oblique and backward trickles write
*their own* retarded phases, which is what the axis collimation, the sharp
gate exponent, and the frequency window all exist to starve.

Plane normals re-align the same way: each plane relaxes toward the
flux-weighted orientation of the planes it exchanged energy with (headless
sign-corrected) **and toward containing the flux direction** (n̂ pulled
perpendicular to û, so the planes' intersection line rotates onto the
direction of transfer), plus a small thermal tumble. *Cells re-align; energy
transfers; nothing moves.*

---

## 3. Conversion between modes (the beat) [P]

The two frequencies interact through their beat `β̇ = ω₁ᵉ − ω₂ᵉ`. Once per
completed beat cycle (β wraps 2π) — a *complete cycle*, per PRINCIPLE §3 —
the cell converts:

- **Condensation (F → D, consuming space):** if `E_e > e_cond`:
  `δ = f_conv·(E_e − e_cond)` moves E_e → E_m, and additionally
  `δ_s = min(s_pull·δ, E_s − e_floor)` moves **E_s → E_m**.
  *Matter is converted space* (PRINCIPLE §4.2): locking pattern consumes the
  cell's own space energy, shrinking its area of effect.
- **Evaporation (D → F,S):** if `E_m + E_e > cap` (over-full):
  `δ = f_evap·excess` leaves E_m, split back in the ratio it was built
  (1 : s_pull) to E_e and E_s.

Loaded cells beat slower (both ωᵉ share the detune factor) → conversion near
mass is time-dilated, consistently with §5.

---

## 4. The only law, mechanically

Every operation in the kernel is a **paired ledger move**: the same double
subtracted from one account and added to another (cell↔channel, or
intra-cell mode↔mode). There is no source term anywhere. Total energy

```
E_tot = Σ_cells (E_s + E_m + E_e) + Σ_channels Σ e_mid
```

is conserved to floating-point roundoff, and the kernel reports the relative
drift every diagnostic row. There is no equation of motion for conservation
to constrain — transfer gates and conversion gates are the entire dynamics.

---

## 5. Curvature and time dilation [P → measured]

Condensation drains `E_s` at the locus of the pattern; `r = r0·(E_s/e_s0)^⅓`
shrinks; contact areas with neighbours drop; some channels pinch off
entirely. Around a dense region the space fabric is thinner and less
connected — *where there is matter, there is less space*, and the
surrounding structure must account for it (PRINCIPLE §4.3).

Measured defect: `ΔA = Σ_links (A₀ − A)` (contact-area deficit) and the
count of dead channels. Prediction in the spirit of the combinatorial Gauss
law: **ΔA grows linearly with converted dense energy** while conversion
proceeds. Time dilation: mean ωᵉ in the core vs far field < 1.

---

## 6. Two planes, two frequencies → Bell [P]

A cell's transverse pair (two planes, two harmonic frequencies, relative
phase) is structurally a polarization qubit; the plane-overlap factor
`(n̂ᵢ·n̂ⱼ)²` in every transfer is a Malus cos² law. A **pair conversion**
(one parcel of field energy splitting into two packets) leaves the two
packets sharing **one unfinished joint harmonic** — a single conversion
process spanning both (transfer mode T is one object with two ends), not two
cells each owning a private phase.

Measurement = transfer through an analyzer cell whose plane is set at angle
a: passes with the kernel's own cos² law. First passage anywhere completes
the joint cycle and collapses the shared phase; the partner's subsequent
passage follows cos² against the collapsed phase. Coincidence statistics:

```
E(a,b) = cos 2(a−b)     →  CHSH S = 2√2 ≈ 2.828  at standard angles
```

Control: the same source with a *pre-assigned private* phase λ per packet
(a local-hidden-variable model with the identical cos² responder) yields
E = ½cos2(a−b) → S = √2 < 2. **The violation is bought exactly by the joint
unfinished harmonic** — the frequencies' interaction as one process — not by
the response law. This is not an LHV evasion of Bell's theorem; it is the
model agreeing with quantum statistics because a pair is one conversion, not
two objects with properties.

`mode=bell` runs both as event-level Monte Carlo on the kernel's transfer
law. Field-level entangled transport (registry of joint phases across the
cell fabric) is open — §11.

---

## 7. Experiments

| id | cfg | question | success |
|---|---|---|---|
| E1 | `e1_conserve.cfg` | ledger closure with all mechanisms firing | \|drift\| ≲ 1e−11 |
| E2 | `e2_pulse.cfg` | tilted field packet propagates; front speed vs C | coherent front, v ≲ C |
| E3a | `e3_blob.cfg` | untilted dense blob is static and self-sealed | containment ≈ 1, centroid still |
| E3b | `e3_blob_tilt.cfg` | tilted blob translates while no cell moves | centroid drifts along k̂ |
| E4 | `e4_curve.cfg` | condensing packet curves space | ΔA ∝ E_m (linear), ωᵉ dip in core |
| E5 | `e5_bell.cfg` | joint harmonic vs LHV control | S_joint ≈ 2.83, S_LHV ≈ 1.41 |

---

## 8. Background audit

| object | persists? | merely re-valued? | verdict |
|---|---|---|---|
| cell scaffold (positions, candidate links) | yes | — | **flagged**: sanctioned approximation ("boundaries not real, but useful"); dynamics touches only relational geometry; audit-honest residue |
| cell energy ledgers | cell = its energy; empty ⇒ physically ceased | — | OK |
| live connectivity (A > 0) | no — overlap is energy-dependent | no | OK |
| plane normals, phases | internal harmonic content h(v) | — | OK (algebra, not place) |
| in-flight e_mid | born/dies with its cycle | no | OK (mode T) |
| total energy | scalar | no location | OK — the law |
| positions in diagnostics | reconstruction only | — | OK if never fed back (they are not) |
| momentum | not tracked | — | correctly absent; tilt is its ghost |

---

## 9. Build & run

```
gcc -O2 -o v89/cellfab v89/cellfab.c -lm
./v89/cellfab v89/cellfab_runs/e2_pulse.cfg          # config file
./v89/cellfab e2_pulse.cfg k_dep=2.0 T=30            # overrides
```

Key parameters (defaults in `cellfab.c`, all overridable):

| knob | default | meaning |
|---|---|---|
| `L, dmin, r0, rjit` | 24, 1.0, 0.85, 0.06 | box, dart min separation, area-of-effect radius, jitter |
| `C` | 1 | conversion rate (cycle time per link = d/C) |
| `w1, w2` | 1.5, 0.93 | the two plane frequencies (field / dense); keep ωd/C ≈ π/2 (§2 window) |
| `q_detune` | 0.35 | load detuning strength |
| `gamma_res, gamma_res_m` | 0.25, = | resonance acceptance Γ per sector (narrow `gamma_res_m` = rim seal) |
| `p_gate` | 4 | tail-gate sharpness (8 for long-lived dense tilts) |
| `k_dep, k_dep_m` | 1.2, 1.0 | deposit rate scale (field; dense multiplier) |
| `cap` | 2.5 | internal bound (E_m+E_e capacity) |
| `e_s0, es_floor` | 1.0, 0.05 | space energy per cell; conversion floor |
| `e_cond, f_conv, f_evap, s_pull` | 0.3, 0.25, 0.5, 0.5 | beat-conversion knobs; space pulled per locked unit |
| `kappa_lock, lock_floor` | 0.9, 0.005 | clock entrainment strength; imprint inertia floor (×cap) |
| `kappa_align, sigma_tumble` | 0.5, 0.01 | plane alignment rate, tumble |
| `comb_limit` | 6 | C1: max p·q in the partial comb (1 = round-1 unison-only) |
| `rough_k, gamma_rough` | 0.35, 0.5 | C2: roughness radiated fraction; Plomp–Levelt peak detuning |
| `mob_sym, mob_floor` | 1, 0.01 | C3: mutual dense coupling; sympathetic readiness floor (×cap) |
| `npairs, pair_x0, pair_x1, pair_p, pair_q` | 48, 0.4, −1, 1, 1 | standing-pair seeds: occupancies (x0<0: on-curve; x1<0: same) and interval |
| `field_J` | 0.06 | round 3: field hop coupling (group velocity ∝ J; slit runs use 1.8) |
| `dt, T` | 0.02, 40 | integrator step, span |
| `debug` | 0 | per-diag gate/flux statistics (`# DBG0/1` rows) |

Output: `#`-prefixed header + TSV diagnostics on stdout; `# RESULT` summary
lines at end; optional snapshots `snap_dir/cells_NNNNNN.tsv` when
`snap_every > 0`.

---

## 10. Results — first campaign (2026-07-27, seed 20260727)

Logs in `cellfab_runs/*.log`. Typical fabric: 10–32k cells, 75–282k
channels, mean degree ~18.

| id | claim | measured |
|---|---|---|
| E1 | ledger closure, all mechanisms firing | **rel_drift 3.7e−16** over 2000 t.u.-steps (E0 = 7365.28402937 unchanged to the last digit printed) |
| E2 | tilted field packet propagates | front speed **v = 0.48 C** (r² = 0.93), launch asymmetry ~300:1 (gate+x 0.43 vs gate−x 0.006); the shortfall from C is redeposit lag + cone obliquity, not a posted limit |
| E3a | untilted blob static & sealed | centroid speed **0.0026** (12× below E3b), \|Δcm\| < 0.06 over 80 t.u.; containment 0.91 → 0.64 (slow halo — kinetic freeze, not binding) |
| E3b | tilted blob translates, cells never move | centroid **7.0 → 12.3 (+5.3 units) over 80 t.u.**, direction cos to k̂ = **0.9996**; late-half fitted speed 0.032 (early ~0.15, decaying as the tilt decoheres); containment 0.90 → 0.48 |
| E4 | condensation curves space | contact-area defect **ΔA = 4.79·E_dense, r² = 0.99** (E1's independent foam: 4.10, r² = 0.998 — slope foam-dependent, linearity the claim); core radius −5.4% vs far field; core tick rate ratio 0.9961 (time dilation, small because converted mass is modest) |
| E5 | joint harmonic beats the classical bound | **S_joint = 2.826** ≈ 2√2 = 2.828; identical cos² responder with private phases: **S_LHV = 1.414 < 2** (500k trials/combo) |

Conservation held at the floating-point floor in every run, including the
ones with heavy churn (E3b cycles ~40% of the blob through transfer mode at
any moment).

### Round 2 (music mechanisms C1–C3 active; see CONSONANCE.md Part VI)

Same seeds, kernel now carrying the comb, roughness radiation, and mutual
coupling. Conservation at the floor everywhere. E2/E5 byte-identical (field
and bell untouched). E1 curvature fit sharpened to r² = 0.9997. E3a
containment 0.64 → **0.84**; E3b trades reach for cohesion: +2.3 units at
containment **0.81** (was +5.3 at 0.48), radiating 30.9 units of roughness
en route. E4's defect now relaxes when its mass radiates (no orphaned
curvature). E6/E7 reproduce with the shed-ratchet cleaner; E8 weighs the
comma; E9 measures the interval hierarchy (unison ≫ fifth ~20 t.u. >
octave). Details and tables: `CONSONANCE.md` Part VI.

### Tuning notes (what it took, kept honest)

1. **Locking only at delivery** (once per d/C) bottlenecks the front at
   ~0.1 C. Entrainment must begin at deposit — resonant joining is a
   boundary process, not an arrival event.
2. **Single-plane conductance** lets oblique links carry half-open flux; the
   packet diffuses (~0.2 C, melting tilt). The two-plane axis condition is
   what collimates propagation.
3. **The wrap window**: at w1 = 2.2 (ωd ≈ 2.5) the backward gate wraps to
   0.2 open and its deliveries imprint wrongly-retarded phases; the tilt
   thermalized in ~2 t.u. Moving to ωd ≈ π/2 sealed it (0.006).
4. **Dense sector wants a narrow Γ and sharp gates**: with Γ_m = 0.12 the
   rim seal never engages (res ≈ 0.98 in the halo) and the blob dissolves
   isotropically; Γ_m = 0.02 + p_gate = 8 gives the sustained directed
   drift of E3b.

---

## 11. Open

* ~~**Binding.**~~ **Retired and answered in `CONSONANCE.md`**: what was
  sought under this name is consonance — closure of mutual conversion
  cycles. The separation ladder (ωᵢᵉ+ωⱼᵉ)d/C = 2πm, the tuning curve
  x*(d), the tempered comma, entrainment acquisition, and the one-sided
  relaxation ratchet are all derived from this kernel's existing gate laws
  and measured (E6/E7). What remains open there: extended *many-cell*
  patterns still shed halos (blob containment 0.64/0.48 at 80 t.u. in
  E3a/E3b) — rings and composite locks are the path, not an attractive
  channel.
* **Tilt decoherence.** The imprint channel that carries a wave is also what
  erodes it (oblique writes). E3b's drift decays on ~30 t.u. A closed
  integer lock (CONSTRUCTION §2.2) rather than a continuous phase may be the
  real answer; the continuous kernel can only starve the channel.
* Field-level Bell: carry the joint-phase registry through actual fabric
  transport and analyzer cells (event-level version is §6).
* The tilt→speed law v(k) and whether entrainment pumps tilts toward
  k = ωᵉ/C (a c-attractor for massless patterns; what resists it for dense
  ones — inertia [G]).
* Cell birth/death: currently cells persist as handles when drained; a true
  rewrite version would recycle them (remove the flagged audit residue).
* The two base frequencies are uniform constants [P]; species structure
  would make ω-words per-cell integer content (CONSTRUCTION §2.2).
