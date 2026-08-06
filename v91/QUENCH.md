# QUENCH — genesis by dynamics: compress fast, release fast (user-directed 2026-08-06)

User question (2026-08-06, verbatim): *"Have we tried compressing the
field quickly in the first N frames (eg 10 or 50) to a high degree
then quickly releasing the bounds?"* Answer measured in the record:
never. This campaign runs the config-only forms of that protocol on
the standing binary (kernel and laws untouched; the scheduled-anneal
/ breathing-box knobs are kernel work, conditionally authorized by
the user ONLY if these arms show a pulse). Registered before any run;
results land in §4 per arm; nothing is adopted.

## 0. Thesis

Every creation-class mechanism tried so far is a GATE, and four
rounds across three lanes closed that door: gates preserve, they do
not create. A quench creates by DYNAMICS: drive the medium through a
transient over-compressed state and release; nucleation happens in
the sweep, not through a gate. Two measured facts aim it: (i) the
law-regime vacuum CONDENSES above-threshold beams into matter trails
(XSEC/B7 — the condensation door exists and is strong: cond ≈ 77 of
a ~180-unit beam); (ii) during release the local fullness sweeps DOWN
THROUGH the balance locus x̂* — and CANTUS says chords are the only
stable fits there. The registered question: does sweeping the balance
locus ever nucleate persistent bonded matter — episodes + mutual
stamps that outlive every null?

The DETECTOR is the identity meter (IDENTITY.md, physics-silent):
nucleation is claimed only for a connected cluster of ≥3 cells whose
episode ages exceed 3× the bath mean (>1600 t.u.) with ≥2 mutual
stamps older than 1500 t.u. (bath stamp tail: P < 0.3%), alive at
run end with roughly balanced books. No eyeballing.

## 1. Arms (exact commands, standing binary, from `v91/`; seed 20260802)

The beam is a one-shot packet: its transit IS the fast compression
(~σ/C t.u. locally — the user's "first N frames"), and its passage IS
the release. No knob is needed for the vacuum forms.

```
# Q1  beam dump into quiet foam, law regime (XSEC class + long watch)
./freecell exp=slit slit_mask=3 L=64 T=1500 sigma=4 slit_sy=3 kx=2.0 amp=4 slit_srcx=14 slit_screenx=40 slit_t0=0 slit_t1=30 diag_every=200 seed=20260802 k_rad=0.05 p_rad=4 rad_clock=0 par_tau=10 qatom_every=200 snap_every=2500 snap_file=runs/quench/q1_beam4.fcs > runs/quench/q1_beam4.log 2>&1

# Q2  deep quench: amp 4 -> 8 (only knob changed)
./freecell exp=slit slit_mask=3 L=64 T=1500 sigma=4 slit_sy=3 kx=2.0 amp=8 slit_srcx=14 slit_screenx=40 slit_t0=0 slit_t1=30 diag_every=200 seed=20260802 k_rad=0.05 p_rad=4 rad_clock=0 par_tau=10 qatom_every=200 snap_every=2500 snap_file=runs/quench/q2_beam8.fcs > runs/quench/q2_beam8.log 2>&1

# Q3  packet through the LIVE warm bath (exp=pulse at quench amplitude)
./freecell exp=pulse bath=1 noise_amp=0.5 L=24 T=1500 amp=3 sigma=3 kx=1.1 diag_every=200 seed=20260802 k_rad=0.05 p_rad=4 rad_clock=0 par_tau=10 qatom_every=200 snap_every=2500 snap_file=runs/quench/q3_pulse_bath.fcs > runs/quench/q3_pulse_bath.log 2>&1

# Q4  dense-grain quench: driven two-blob collision (MOTION #28 seed)
./freecell exp=blob2 amp=1.6 sigma=2.2 blob2_sep=8 blob2_kx=0.9 L=24 T=1500 diag_every=200 seed=20260802 k_rad=0.05 p_rad=4 rad_clock=0 par_tau=10 qatom_every=200 snap_every=2500 snap_file=runs/quench/q4_blob2.fcs > runs/quench/q4_blob2.log 2>&1
```

Controls, standing (no new runs): the warm-bath identity baselines
are M-I1 (episode mean 530 t.u.; stamp ages memoryless mean 261;
`runs/identity/mi1_bath_warm.log`); the quiet foam does not densify
spontaneously (FC0, gated). 2D for Q1/Q2 (the slit substrate — B9's
dimensionality caveat applies and is recorded).

## 2. Closed-form expectations (before any run)

1. **The nulls to beat.** RADIANCE: released dense hoards decay with
   t_half 80–140 → the condensed trail should be substantially gone
   by t ≈ 500. CANTUS geometry wall: quench condensate is
   unison-grade matter; unison cannot bond at the balance → no
   chords expected. FC0: nothing condenses without the beam.
2. **The hypothesis's one opening.** During trail decay each dense
   cell's x sweeps down THROUGH the chord operating band (x_D ≈
   0.55–0.65, x_U ≈ 0.28); adjacent cells crossing at m=2-compatible
   loads could, in principle, latch a fifth. Nothing measured says
   they will; everything measured (unison wall) says they won't —
   the arm exists because the sweep-through-balance configuration
   has never actually been produced.
3. **Q3 caution registered:** the pulse deposits into a LIVE bath —
   any "nucleus" must beat the bath's own 530/261 baselines, and
   spatial attribution needs the snapshots (global PAR is
   bath-dominated in Q3).
4. **Q4 expectation:** blobs profile-merge (COE: v_COM ≤ 1e-3 of
   closing) — the "collision" is slow; the compressed interface is
   the quench zone; RADIANCE null applies to the merged hoard.

## 3. Bars

* **Q-DET (the only claim bar):** a cluster per §0's detector at run
  end in any arm → NUCLEATION CLAIMED, arm replicated ×2 seeds
  before anything further; then and only then the scheduled-anneal
  knob (kernel) is authorized per the user's condition.
* **Q-NULL (per arm):** trail/hoard decay curve reported against the
  RADIANCE t_half band; episode/stamp demographics reported against
  M-I1; "no chords" checked via locks=0 (k_cant=0 here — no lock
  apparatus: whatever persists must persist by CONTACT alone, the
  nv48-vacuum class) and via the QATOM spectrum (any D-band emission
  line standing over the trail's own band would be a chord
  signature).
* Drift ≤ 2e-15 as always; battery untouched (config-only).

## 4. Results

(registered empty; filled per arm)

## 5. Verdict

(after all arms; decisions are the user's)
