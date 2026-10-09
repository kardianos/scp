# v102 kernel step, K2e — decision (singlet core in the formation event)

**Written:** 2026-10-09, after the runs in `cas/out/k2e_*`. Gates:
`K2e_PREREG.md`. Code:
- `v97/kernel/lockq.c`: `lockq_proj` (`lq_proj_link`, `lq_proj`, the
  call in `lq_map`), the second seeded jump;
- `v97/kernel/freecell.c`: the config fields `lockq_proj`,
  `lockq_force2`, `lockq_fmode2`;
- the tools `cas/k2d_replay.py` (the `proj` argument), `cas/k2e_ab.py`
  (the E2 reference), `cas/analyze_k2e.py` (E3) and `cas/k2e_conv.py`
  (the projection check, §E4).
- `cas/k2d_replay.py` also gained K3's readouts in this commit: seeded
  jumps are replayed, and it prints the profile after the last event and
  the pair's distances. They add output lines and change no existing
  number.

## Ending: **Go**

- E1–E3's bars all pass.
- R12 is now exact in the kernel: under `lockq_p11=2 lockq_proj=1`, each
  jump state is projected onto the allowed subspace before its window
  and rate are evaluated.
- **The physics changes sharply from K2d.** Under bright light, holes
  now form beside holes in the allowed channels, and the lock fills
  with blocks of adjacent holes that are far from bound. K2d's jam at
  three holes, with born-bound pairs, came from the mean-energy window,
  not from R12 (§E3).

## E1. Inert at default: **PASS**

- **The full battery, both kernels** (`cas/out/k2e_battery.out`):
  **ALL GREEN (87 bars)**, bar table identical to K2d's.
  - 70 of 71 run logs are byte-identical to K2d's.
  - `p2lc` (the Go kernel) differs only in its wall-clock timings.
- **Kernel stdout against the K2d binary** (`cas/out/k2e_e1_stdout.out`)
  is byte-identical for:
  - `exp=bath`, `blob`, `slit`, `pulse` (T = 20);
  - K2d's D3 configurations `M2_b_p2_s1` and `M2_b_p0_s1`.
- There are no new compiler warnings (9 before, 9 after). The gate binary
  was rebuilt from the working tree and is byte-identical.

## E2. A/B of the projected jump: **PASS**

The setup:
- `lockq_p11=2 lockq_proj=1 lockq_rate=0`, 2000 steps (tL = 40);
- the reference is `cas/k2e_ab.py`, which projects each jump by the
  filter exp(−τ H_c) with H_c = Σ_l Π_l and τ = 1000;
  - the pre-registration said τ = 200; the larger τ only filters harder;
  - after the filter, |H_c x|/|x| ≤ 1.3e-16 in every case;
- outputs in `cas/out/k2e_e2_*`.

| case | projected weight, kernel | reference | ⟨H⟩/J0 at tL (both) |
|---|---|---|---:|
| 3+2 n = 8, (a) site 1 | 0.258305709950745 | 0.258305709950748 | −14.9328964406 |
| 3+2 n = 8, (b) site 2, then site 3 (mode 1) | 1.0, then 0.615265225177836 | 1.0, then 0.615265225177866 | −15.3530902467 |
| W n = 6, (a) site 1 | 0.190113749576395 | 0.190113749576391 | −16.5834550989 |
| W n = 6, (b) site 2, then site 3 (mode 3) | 0.5, then 0.709630450838198 | 0.5, then 0.709630450838202 | −19.5534380319 |

- The nh, S^z and I^z profiles agree at every printed digit (10). The
  bar is 1e-6.
- A jump with no (1,1) link keeps its unprojected weight: 1 in 3+2, ½ in
  W, where the removed quantum has a definite mode.
- Every E2 projection converged in 2 sweeps.

## E3. Singlet core in light: **PASS** (bars), and the readouts

The setup: K2d's D3 configurations at `lockq_p11=2` with `lockq_proj=1`.
- 3+2 n = 10 and W n = 7;
- bright light (amp 4, σ 3) and dim light (amp 1, σ 2.2);
- seeds 1–3;
- outputs in `cas/out/k2e_e3_*`, and the analysis in
  `cas/out/k2e_e3_analysis.out`.

**The bars, every run:**
- |drift_rel| ≤ 3.5e-15;
- `check_violations` = 0 and `window_viol` = 0;
- the replay with projection matches the kernel's final lock line to
  ≤ 5.0e-11 (bar 1e-8).

**Readouts.** Dim light made one hole in 3+2 (at site 1, beside the end
cell) and none in W, in every seed, as in K2d.

Bright light:

The bound weight is the whole lock's weight on levels of its sector
with B_k < −0.05 J0.

| fabric | seed | events | beside a hole | final holes (sites) | adjacency | mean B / J0 | bound weight |
|---|---|---:|---:|---|---:|---:|---:|
| 3+2 | 1 | 6 | 2 | 1, 2, 3, 4, 6, 8 | 3.0 | +16.56 | 0 |
| 3+2 | 2 | 6 | 4 | 1, 2, 3, 4, 5, 7 | 4.0 | +15.00 | 0 |
| 3+2 | 3 | 6 | 3 | 1, 2, 3, 4, 6, 7 | 4.0 | +16.22 | 0 |
| W | 1 | 2 | 0 | 1, 3 | 0 | +4.46 | 0 |
| W | 2 | 2 | 1 | 1, 2 | 1.0 | +8.33 | 0 |
| W | 3 | 3 | 1 | 1, 2, 3 | 2.0 | +12.30 | 0 |

Against K2d's D3 at the same `lockq_p11=2`, without the projection:

| | 3+2 bright, final holes | events beside a hole | bound weight |
|---|---|---:|---|
| K2d (mean window) | 3, always at sites 2, 4, 7 | 0 | 0.51, 0, 0.51 |
| **K2e (R12 exact)** | **6 of 8** | **9 of 18** | **0, 0, 0** |

**What this shows:**

1. **Under R12, holes form beside holes.**
   - A jump beside an N = 1 cell loses only its singlet part. The rest
     fires.
   - In 3+2 (P = ε) the allowed state of two adjacent N = 1 cells is
     their triplet. A run of adjacent N = 1 cells with no singlet on any
     link is symmetric under every exchange, so it has the largest
     spin. That follows from the rule; it was not measured separately.
2. **The blocks are not bound.**
   - Bright 3+2 ends with six holes and mean B = +15 to +17 J0.
   - Bright W ends with B = +4.5 to +12 J0.
   - No run has any weight on bound levels.
   - The lock cannot shed this energy (K1-2). The blocks are held by
     the closed lock, not by binding.
3. **The window still costs adjacency, but only a little.** Each N = 1
   neighbour of the new hole (another hole, or an end cell) moves the
   event's δ by about −0.035 in rate units (≈ 3.6 J0). That is the
   lock's measured energy cost of the neighbour; it is not decomposed
   here:

   | N = 1 neighbours | δ, 3+2 (events) | δ, W (events) |
   |---|---|---|
   | 0 | −0.036 to −0.038 (6) | −0.036 (2) |
   | 1 | −0.057 to −0.072 (9) | −0.073 to −0.076 (4) |
   | 2 | −0.087 to −0.100 (3) | −0.110 (1) |

   - The band edge is −0.12, so a hole can still form between two N = 1
     cells.
   - An isolated hole already forms at δ ≈ −0.036, because the vacuum
     pitch (1.65) sits below the resonance. Adjacency pushes events
     toward the band's edge, not out of it.
   - The rate's band factor √(band² − δ²) falls only by about 40% from
     δ = −0.036 to −0.10.
4. **K2d's born-bound pairs were an artifact of the mean window.**
   - Under score without the projection, the singlet part priced at
     104 J0 pushed the whole adjacent jump outside the window. Adjacency
     was then blocked outright.
   - The lock jammed at three separated holes, and when one had the
     complementary mode, two of them sat at d = 3, bound.
   - R12 exact removes that block. The lock fills with adjacent holes,
     and no seed is left with two holes alone at d = 3 in complementary
     modes.
5. **The singlet weight of the final states** is 5e-5 to 2e-3. It is
   the virtual admixture the penalty allows during evolution. The
   pre-registration estimated (J0/104 J0)² ≈ 1e-4.

## E4. The projection's convergence (reported, not gated)

The kernel caps the alternating projections at 200 sweeps, with a
tolerance of 1e-15 of the norm. It reports the most sweeps used. Over
all candidate jumps of a run:

| run | most sweeps |
|---|---:|
| 3+2 bright, seeds 1, 2, 3 | 187, **201**, **201** |
| 3+2 dim | 27 |
| W bright | 18, 19, 21 |
| W dim | 0 |

201 means that in two runs at least one candidate jump stopped at the
cap before reaching the tolerance.

`cas/k2e_conv.py` replays each bright run on the kernel's own path.
Before each event it projects **every** candidate jump twice:
- by the 200-sweep rule;
- where that rule stopped at its cap, onward to deep convergence (a
  removal per sweep ≤ 1e-34 of the norm, at most 100 000 sweeps).

Any result of alternating projections is x minus a sum of vectors in
the ranges of the Π_l. Once its singlet weight is zero, it is the
orthogonal projection. The check reports, for each candidate:
- that weight;
- the overlap ⟨x₂₀₀|x − x₂₀₀⟩;
- the change in the rate weight when the projection is continued.

Outputs: `cas/out/k2e_e4_events_*` (every event moment) and
`cas/out/k2e_e4_sampled_*` (every 0.5 t.u. as well).

| run | moments checked | most sweeps | stopped at the cap | rate-weight change when continued | residual singlet weight | ⟨x₂₀₀\|x − x₂₀₀⟩/\|x\|² |
|---|---|---:|---:|---:|---:|---:|
| 3+2 bright s1 | 6 events | 120 | 0 | 0 | 2e-28 | 8e-16 |
| 3+2 bright s2 | 6 events, and every 0.5 t.u. to t = 15.0 | 201 | 13 | 3.5e-16 | 4e-27 | 1e-14 |
| 3+2 bright s3 | 6 events, and every 0.5 t.u. to t = 18.0 | 201 | 10 | 3.7e-16 | 8e-27 | 9e-15 |
| W bright s3 | 3 events | 18 | 0 | 0 | 4e-32 | 6e-17 |

- The cap is reached only in the six-hole sectors, where five or six
  adjacent N = 1 cells (holes and the end cell) make the alternating
  projections slow.
- There, continuing to deep convergence changes a candidate's rate
  weight by at most 3.7e-16, and the fired states not at all.
- **So the 200-sweep cap costs nothing measurable.** It is kept.
- The sampled check was stopped by hand at t = 15.0 (s2) and t = 18.0
  (s3). The remaining six-hole samples took several minutes each. They
  are the same sector as the last checked ones.

**A rejected reference.** The first version of the check used a MINRES
solve of H_c z = H_c x as the "exact" projection. In floating point
its z picked up a kernel component: ⟨x_d|x − x_d⟩ = 5e-4, against
2e-17 for the 200-sweep result. The MINRES reference was dropped. The
200-sweep result was the correct one.

## Limits stated

- **Formation does not see binding.** The window (±12.5 J0) and the
  band factor discriminate only at the scale of the lost (2,1) links
  (≈ 3.6 J0 per neighbour). Binding (0.13–0.50 J0) is far below that.
  What light builds in the lock is set by where and in which mode the
  light lands, not by the lock's ground state.
- The exact sectors cap the lock at 3+2 n ≤ 12 and W n ≤ 7. W n = 7
  has only five bulk sites.
- T = 20. Nothing here shows the long-time fate, but a closed lock
  holds its energy (K1-2), so the blocks cannot relax.

## Then

K3, under R12, with its own pre-registration: a lock in light over
long times. These results change what K3 should expect:
- bright light builds unbound blocks of adjacent holes;
- dim light makes one hole per pulse;
- a born-bound pair now needs its two holes to land apart (d ≥ 2), in
  complementary modes, before the lock fills.
