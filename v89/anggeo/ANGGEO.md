# ANGGEO — alternative kernel A: natural angular geometry

**Opened 2026-07-31 on user proposition:** vary the angular
interference patterns to be slightly more natural — if the two planes
are perpendicular, they should FORM that way, not be rigidly enforced.
Experiment on a **copy** of the kernel (`cellfab_ag.c`; production
`cellfab.c` untouched, per kernel policy).

## 1. The facts pass — where right angles actually live today

Map of every angular touchpoint in the production kernel:

| site | behavior | square enforced? |
|---|---|---|
| vacuum cells at birth (cellfab.c:850) | n1, n2 = two INDEPENDENT random unit vectors | **no** — β free from birth |
| pass-7 dynamics (2866) | per plane, independently: pull toward flux-partners' plane (parallelism, sign-folded), pull to CONTAIN the flux direction (the axis rotates onto active channels), tumble noise | **no** — no intra-cell n1↔n2 coupling anywhere in law |
| structure seeders — pairs (1103), ring (1190), shell (1324), net (1492) | exact orthonormal frame: n1 = t1 ⊥ û, n2 = û × t1 | **YES — the only rigid squares in the system, imposed by hand** |
| prealign (1034) | beam-prep apparatus frame | yes, kept deliberately (e3b/g3 beam tool) |
| cbeta | the two-plane BEAT clock (cbeta += (w1e−w2e)·dt), not a geometric angle | n/a |

**So the fabric is already "natural" everywhere except at our seeded
structures.** The proposition becomes a sharp A/B: seed structures
with FREE frames and watch what the dynamics does.

## 2. Kernel A (`ang_free`, `ang_diag`; defaults 0 = production-exact)

`ang_free=1`: the four structure seeders draw BOTH normals at
independent random azimuths in the transverse plane span{t1, û×t1} —
each plane still CONTAINS the channel (the functional requirement:
axi = 1 exactly as before); only the square is released (β uniform).
`ang_diag=1`: `# ANG` rows — mean |cos β| over free cells and over
loaded cells (Em > 0.05). β = 90° ⇒ |cos β| = 0; isotropic free ⇒ 0.5.

## 3. Pre-registered questions and predictions

| id | question | prediction (to be scored honestly) |
|---|---|---|
| **AG-1** | vacuum: what β distribution does the free dynamics maintain? | isotropic (mean\|cosβ\| ≈ 0.5), tumble-dominated; any structure = a hidden angular channel in the alignment term |
| **AG-2** | pairs seeded square vs free: does lock quality care about β? does β(t) of members drift toward 90°? | the law's angular content is channel-containment + PARTNER parallelism, not intra-cell squareness: free-frame pairs start conductance-reduced (plane-2 azimuth mismatch between the two ends, gpl ∝ cos²Δφ), kappa_align pulls the partners parallel in ~10 t.u., locks then equivalent to square-seeded; **β within a cell stays free — perpendicularity does NOT form; inter-cell frame SYNC does**. If instead locks stay degraded, β is load-bearing and a missing intra-cell law (frame stiffness) is exposed; if β → 90°, the user's hypothesis wins and the square is a genuine attractor |
| **AG-3** | full battery on kernel A, ang_free=1 | pair/ring/net experiments transiently slower to acquire; end-state bars hold |

Doctrinal stake: if AG-2 shows locks indifferent to β, the "square
frame" was scaffolding, never law — the real angular order of the
fabric is (i) axes onto channels, (ii) partner planes parallel. If β
matters and nothing restores it, "natural formation" requires a
derived intra-cell frame law — a finding either way.

## 4. Runs

On the S1 substrate (V2s + anneal extras), kernel `anggeo/cellfab_ag`:
AG-1: vacuum, T=200, ang_diag=1 (β drift control).
AG-2: e6-pattern pairs ×2 arms (ang_free 0 vs 1), T=60, ang_diag=1,
score: gg(t), # ANG member trace, acquisition time.
Verdicts appended below.

## 5. Verdicts

(appended as measured)
