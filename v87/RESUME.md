# Resume after restart — 2026-07-26

## The one thing that matters

**EX-1 is running on a remote Vast.ai GPU and will keep running through a local
restart.** The simulation process lives on the remote instance, detached; the
local machine only monitors it. Restarting here does not touch it.

Billing continues while it runs: **$0.169/hr**, so a full 20 h run is ~$3.40.

## Reconnect

```
sim_setup(name="gpu1", executor="remote", host="ssh9.vast.ai", port=13566)
sim_run_status(name="gpu1", id="ex1_boost_v002")
```

`sim_setup` with an explicit host/port attaches to the existing instance rather
than provisioning a new one. Identifiers:

| field | value |
|---|---|
| instance name | `gpu1` |
| run id | `ex1_boost_v002` |
| host / port | ssh9.vast.ai : 13566 |
| Vast instance id | 45933567 |
| machine id | 58241 |
| GPU | Tesla V100-32GB |
| rate | $0.169/hr |

Runner state is file-based and survives restart: `~/.scp-runner/state/gpu1/`
holds `instance.json` and `ex1_boost_v002.status.json`. The `.done` / `.failed`
marker files appear there at terminal state.

## EX-1 status at restart

Init **complete** (the ~67 min serial Gauss projection), physics integrating.

```
Gauss projection: 1377 CG iters, gauss_max(0) = 1.84e-13     <- exact floor, clean
GPU allocated 20.38 GB (45 N^3 arrays)
INIT: E_total = 3.9345e+02   Q_phi = 266.33   E_em = 1.9258
      gauss_max = 1.837e-13  gauss_l2 = 1.336e-14
```

Seed validation held: the kernel's E_matter = 391.52 against the 391.53 measured
locally before upload (0.01%), and Q = 266.33 against 266.29. The boosted seed
loaded correctly and the Gauss projection rebuilt E properly.

Config: N=384, L=55, T=2000, g=0.05, absorbing BC (width 3.0, rate 0.01),
snap_dt=100, diag_dt=0.2. Seed `/root/ex1_seed.sfa`, output `/root/ex1.sfa`.

**On completion:** download `ex1.sfa` and `ex1_diag.tsv`, then
`sim_teardown(name="gpu1")` — teardown verifies downloads before destroying.
Do not tear down before the download completes.

## What to do with the result

Pre-registered from EX-2: radiative losses should be dominated by above-gap
matter waves, not gauge radiation. The analysis chain already exists:

- `v86/n_battery/sfa_radial` — radial energy/charge decomposition
- `v86/n_battery/sfa_momentum` — P, COM, stress flux (validated on x10c)
- `v86/n_battery/ex2_spongeflux.py` — matter/gauge channel split

Measure: co-motion vs stripping of the cloud, cloud binding and radius in
flight, and the ball clock's transport (EX-4). The seed carries
**v_eff = P/E = 0.01998**, 38% of the adiabatic threshold α_f = 0.053.

## Local work

Completed:

**v86 Part 0 + census** (`v86/PART0_RESULTS.md`): N1–N6, N9, N10, HC-1, HC-2,
HC-3, EX-2, N7, D7-lite. D5′ closed — M = E/c² to 1–4 parts in 10⁴, dx-converged
across a 1.43× resolution change.

**v87 Bell rung** (`v87/`): three seats (FABLE, GROK, GEOM) plus an independent
review; `CRANK2_RESULTS.md` corrected against it (`crank2_review.md`);
`NEXT_PROGRAM.md` for the census triage.

**v87 B1 CHSH rung — DONE** (`v87/B1_CHSH_RESULTS.md`). Both instruments carry
all three protocol requirements. `work/geom/bell_grid.c` (analytic layer,
N = 10¹⁰) and `work/kernel/bell_kernel.c` (in-kernel, reads the pair SFA).
Headlines:
- Search bias is a property of the **waveform**: +3.95/√N on the triangle
  (12/12 replicas positive), ~0 on the cosine. The triangle is what the fabric
  produces, so this control was load-bearing.
- **The in-kernel CHSH with an offline per-object readout is |S| ≤ 2 by an
  algebraic identity of the readout**, before the kernel runs. It cannot be
  evidence about fabric locality. CRANK2 §5.4 gap 1 (settings as fabric DOF) is
  therefore a *prerequisite* for that experiment, not a parallel gap. B2 is
  designed in §4.2 (dynamical analyzers, one run per setting pair).
- Transport measured: pair clock at D = 13 stays locked, **|R| = 1.000000**,
  drift −0.0028°/unit over 9.4 clock periods. Favourable against G14 pessimism.
- Fabric phases reach 2√2 under the p = 1 tilt **only** at η = 0.6366, i.e.
  discarding 36.3% — CRANK2 §4 Reading 2 demonstrated on real data.

**Open, ready to pick up:**
- `v87/CRANK2_RESULTS.md` §5.3 — the derivation checklist for the sampling weight
- **B2** — the dynamical-analyzer CHSH, spec'd in `B1_CHSH_RESULTS.md` §4.2;
  the only version of the in-kernel test that is a real locality measurement
- `v86/NEXT_PROGRAM.md` — HC-1-gauged, HC-3-volume, D8b/D8c, N8, HC-4-lite

**Geometry convention that cost two runs on B1:** in both `scp_sim` and the
seed generators **`L` is the HALF-domain**, `dx = 2L/(N−1)`. `N=128, L=19` is
the box [−19, 19] at dx = 0.29921.
