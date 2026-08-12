# v97 — adopted metrics (the conversion-dipole hunt)

**Status:** adopted (user, 2026-08-12). These replace lump-persistence
as the primary scoreboard. Config-only on the carried kernel until a
principle-face is authorized.

The last two years scored hoards: `Em_max`, participation ratio, tagged
centroid, R2d of a painted vortex. Those meters answer “does a store
persist / drift / hold a phase we wrote on it?” The principles say
that is the wrong object. Properties are moments of a **live conversion
cycle**. This file is the hunt for that cycle.

Do not retire the old meters. They remain diagnostic. They no longer
decide the claim.

---

## 0. What we are looking for

A **persistent conversion dipole**:

| clause | meaning | why |
|---|---|---|
| **I and O spatially separated** | intake and outtake do not fire at the same cell | otherwise the book is a point and π-flatness is a theorem |
| **Φ > 0** | throughput `min(I, O)` (or `I` and `O` separately when they match) is nonzero | a sealed store has no temperature and no inertia |
| **dM/dt ≈ 0** | inventory of the motif is stationary (restoring, not a walk) | otherwise it is a leak or a refill, not a species |
| **dipole of J steady** | first moment of the link current does not average to zero | translation / inertia **is** that dipole |
| **Es throughput at constant M** | space book cycles while mass is parked | the g4 sentence: the far field awaits this, not a leak |

If no existing stream satisfies these five clauses, that is the
finding, and it tells you the slot-door (P3) is load-bearing rather
than decorative.

The vacuum C6 ring is the existence proof of a *closed* current. It is
the wrong topology for translation and for a far field (a loop is
monopole-neutral). The object is a **polarized** cycle: condensation
preferring one side of the current, evaporation the other — a helix on
the torus v93 already named, not a tighter Gaussian.

---

## 1. Property table (what to score, not what to hope)

| wanted | look for | do not score as |
|---|---|---|
| **mass** | isolated restoring root of `I − O`, with discrete geometry\* so `M*` is not a continuum | `Em_max` of a hoard; radiance ceiling `Em* ≈ 1.55` |
| **temperature** | attribute of the occupied harmonic stack; fluxes `I, O` are what the open world sees (`theory/T_HARMONICS_AND_WIND.md`) | `noise_amp`; vacuum `Ee`; `Φ` itself |
| **momentum / translation** | dipole of `J` with `Φ > 0` | tagged-centroid `blob_drift` (v93 item 2: lies about the unitary channel); COM of a draining tag |
| **spin** | winding of a closed high-`|ψ|` cycle that is *part of the machine* | a painted `seed_mw` vortex on a blob; R2d of a sparse cloud |
| **charge** | incomplete shell `Q = Z − N` on a cyclic contact graph | renamed `Em`; a painted integer |
| **well / gravity** | Es throughput at constant `M`, or a depression that actually binds a field mode | any far monopole of a sealed lump (P5 says that should be ~0) |
| **shell** | standing field mode trapped in that well (no imported opposite charge) | imported orbitals; kernel-painted degeneracy labels |

Charge remains a **label** until P6 produces a joint-graph energy drop.
Do not grade Coulomb.

---

## 2. How to read the carried meters

The carried kernel already prints pieces. Use them. Do not invent a
parallel ledger.

| already there | use as | warning |
|---|---|---|
| `p1_meter=1` → `fd` | dense-hop current (the right translation signal under pass U) | tagged `blob_drift` is the wrong direction for the unitary channel |
| `p1_meter=1` → `sp` | space-transport current (pass S) | gradient flow of `Es`; zero on a static depression |
| `cond_total`, `rad`, `evap`, `door_in` | ingredients of `I` and `O` | they are **global** sums; the hunt needs them *per motif* and *per side* |
| `drift_rel` | conservation (must stay at the FP floor) | never grade this on `exp=bath` (Em ≡ 0; vacuous) |
| `analyze_winding.py` | field / matter winding `W`, `R2d`, m-spectrum | on a dipole, compute about the motif, not the box centre |
| `analyze_ring.py` | cycle `W`, `R_m`, `Em_min`, keep | vacuum hold ≠ medium hold |
| `analyze_melt.py` | `Em_max`, PR, rms, all-cell COM | diagnostic of stores; not acceptance |
| `analyze_hoard.py` / `_phase.py` | occupancy histogram, within-hoard `R2` | Face C's 189 “hoards” were a cell histogram, not a particle spectrum |
| `scpan` | column-by-name parse of RESULT lines | do not `awk` column offsets |

**Instruments (Go-only, offline on FCS; 2026-08-12):** package
`hunt/`, CLI `./dipole`. H0 honesty suite is green — `H0.md`,
`runs/H0.out`. C is unused.

1. **Motif census** — union-find on live-Em cells via LINK contacts;
   tracked by wrapped centroid. Face C reports a handful of connected
   components, not a 189-hoard histogram.
2. **Per-motif I / O** — frame-to-frame mode-crossing heuristic
   (cells matched by position, not roster index); split fore/aft
   along J. Hops contaminate; sealed/vacuum streams correctly read
   Φ ≈ 0.
3. **Dipole of J** — `J = 2 Im(ψ_i* ψ_j)` on links. Reports Jnet,
   Jdip (monopole first moment), Circ, plus global and tagged sums.
   Vacuum ring: Jnet ≪ Circ (loop, not dipole).
4. **Es throughput** — pass-S pressure law on the motif surface and
   on the Em>0.5 core surface. Leak streams: core `|Thru|` tracks
   `|dEmTot|`.
5. **Separation of I and O** — `(c_I − c_O) · ĵ`. Zero on the ring
   (point book / no I/O).

```
./dipole FILE.fcs                         # per-frame TSV
./dipole -h0 sealed|leak|ring|bathring|current FILE.fcs
./dipole -h1 FILE.fcs [...]               # inhabited / empty
./dipole -tag FILE.fcs                    # census on painted cells
```

**Limit (H0 current SKIP):** snapshot J is not `fd`. `fd` is a
time-integral of hops and is not in the FCS. Instantaneous TJnet
sloshes on e3b; that is reported, not forced into a direction bar.

---

## 3. Acceptance bars (pre-registered; sharpen, never soften)

Staged. Nothing here is a law-table bar. The 87 V3a bars stay the
inertness gate.

### H0 — honesty of the meter (do this first)

On one known stream of each class below, the new meters must report
the *already-measured* verdict. If they do not, the meter is wrong.

| class | stream (examples) | H0 must say |
|---|---|---|
| sealed / frozen hoard | Face C deep corner; R1 quiet HOLD | Φ ≈ 0 or I, O co-located; dipole of J not steady; Es throughput tracks nothing or tracks melt |
| leaking blob | v91 g4 / ASTRO leaking mass | Es throughput **tracks** `dM/dt` (bookkeeping wind) |
| vacuum ring | RING_DNLS A.1 / Face A Strang | closed current, monopole-neutral, I/O not split; winding held |
| bath-killed ring | RING_DNLS A.2 | winding dead ~10 t.u.; no dipole rescue |
| unitary current | v93 e3b `amp_nat=2 amp_logate=1 p1_meter=1` | `fd` seed-robust +x; tagged centroid may still lie |

### H1 — does any existing stream host a dipole?

Offline / config-only on carried knobs. Hunt; do not retune.

**Bar:** at least one motif with all five §0 clauses held over a window
`≥ 100` t.u. (stretch: `≥ 1000`). Seed-robust on ≥3 seeds or recorded
as a singleton.

**Decides:** inhabited or empty. Empty ⇒ P3 (slot-local promotion) is
the first kernel face, not another hop.

### H2 — restoring mass, not a walk

Given an H1 motif (or a candidate), perturb inventory ±10–20% (seed
amp, not a new law).

**Bar:** `I − O` changes sign across the unperturbed `M`, and the
motif returns toward it. A one-sided slope is SLOW_DRAIN / SLOW_REFILL
(already measured for radiance and `rad_dual`) and does **not** pass.

### H3 — topology of the machine

**Bar:** `|ψ| > 0` on a closed cycle that carries the winding, *and*
I/O are polarized along that cycle (helix), not a painted blob vortex.

Vacuum ring without polarization is H0, not H3.

### H4 — space book at constant mass

**Bar:** Es throughput through a surface around the motif is nonzero
while `dM/dt ≈ 0`, and is *not* the g4 leak-tracking signal.

This is the P5 discriminator. A pass does not yet mean 1/r. It means
the space cycle exists. Range coupling of two such motifs is a later
face (v96 W-A, reinterpreted).

### H5 — contact interaction before Coulomb

Two H1/H2 motifs at contact, complementary vs like occupancy deficit
(once a graph can be assigned).

**Bar:** complementary contact lowers joint energy / completes a shell
more than like contact, seed-robust. Null at range is the expected
P5/P6 default, not a failure.

---

## 4. What a green H* is not

- Not a proton.
- Not a licence to adopt `amp_nat`, `door_g`, or `j_metric`.
- Not SI. Pure numbers only (`METRICS` §1 ratios, `PRINCIPLES` P4).
- Not a reason to soften a V3a bar.

A green H1 is the first time the programme has seen the object the
principle describes. That is enough to open a principle-face. It is
not enough to claim correspondence.

---

## 5. First actions (no kernel change)

1. ~~Build the five missing instruments~~ **done** (`hunt/`, `./dipole`).
2. ~~Run H0~~ **done** — `H0.md` / `runs/H0.out`. All five classes
   PASS; current-direction SKIP (snapshot J ≠ fd).
3. ~~Run H1~~ **done — EMPTY.** `H1.md` / `runs/H1.out`. 0 / 203
   streams inhabited. Closest 3/5. Polarized I/O and Es-at-constant-M
   never co-occur. P3 (slot door) is the named first kernel face;
   not authorized by the hunt.
4. Stop. H2–H5 have no motif. Kernel faces wait on explicit
   user authorization (`README.md` §V).
