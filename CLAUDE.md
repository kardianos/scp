# SCP Project Conventions

---

# CURRENT VERSION: v92 — READ THIS FIRST

**v92 IS THE ACTIVE PROGRAMME (opened 2026-08-07, user-authorized: "Open
v92 with the Phase L charter and combined reckoning"): `v92/README.md` is
the single entry point — read it top to bottom, then `v92/RECKONING.md`.
v92 is the AMPLITUDE SUBSTRATE (Phase L: dense transport promoted from
magnitude to complex amplitude — the carrier for momentum, spin, exclusion,
and lawful gravity's open books). STATUS: DOCS-ONLY OPENING per the user's
instruction — the tree is NOT yet carried, nothing has run; the first v92
agent executes STEP ZERO (mechanical carry + baseline battery), then ACT
ONE = the combined adoption reckoning (radiance / workfn / flow-bed) WITH
THE USER before any Phase-L design. v91 is the frozen evidence substrate:
its record (RADIANCE fixed point x̂*=0.62; CANTUS/REGISTRY/IDENTITY —
identity is ontological; ASTRO/COMBINE — species spectra, π-flat far-field
null; QUENCH 1-3 — first creation, the door-phase instrument, foam phase
memory τ_φ≈3 t.u.; BLACKHOLE/HORIZON — the π-flatness theorem, the five-wall
no-go, the forced-hole exterior, the prover's ROSTER_DEATH derivation;
ASYM/FLOW — the bed-digging channel law, six bars green, anti-ignition
held) lives in `v91/*.md`; its battery stands ALL GREEN 93 at inert
defaults (commit f7f5486). All law-candidate knobs remain default-inert
everywhere (= laws_V2g byte-exactly) until the v92 reckoning adopts a
table; v90 remains the evidence substrate below.**

**Start at `v90/README.md` (the charter). The carried foundation is v89:
`v89/README.md`, `v89/PRINCIPLE.md`, `v89/FREECELL.md` — consult v89 freely
for the law and the evidence; it is the substrate v90 stands on.**

**v90 (opened 2026-08-03)** carries forward from v89 exactly three working
methods: the **free-cell substrate** (`v90/kernel/freecell.c` — **C is the
kernel of record**, GPU path preserved; `v90/fab/` is the sanctioned Go-kernel
*experiment*, A/B-verified against it), **local-clock execution** (FREECELL
§2's four measured conditions are the implementation spec), and the **battery
discipline** (`v90/cmd/battery` — **Go replaces python** as the harness
language). The law is unchanged: laws_V2g constants VERBATIM,
`v89/battery/laws_V2g.cfg` canonical. Goals in order: re-verify the standing
experiments on the free substrate (double slit first), larger simulation
areas, emergent electron shells, the Go volumetric viewer, then numeric
speedups; stretch (user, 2026-08-03): a Pauli-exclusion demonstration and
shell visualization, if reachable. The ratchet rule governs v90 exactly as it
governed v89.

## Do not consult any version before v89

Not v88, not v87, not the census, not the Q-ball era. **Not for concepts,
models, code, methods, parameters, or experimental design.**

Every version before v89 was built on a **background** — a permanent index set
with evolving contents. That assumption was identified and rejected explicitly,
and then reintroduced anyway, repeatedly, in different formalisms. It returns
through inherited habits of construction rather than through reasoned choice,
which is why the defence has to be *not reading the prior material* rather than
reading it carefully.

Specifically, do not go to earlier versions for how to represent state (arrays
on immortal coordinates), how to write dynamics (PDEs on a fixed index set),
what a particle is (profiles — shapes *in* something), instruments and seeds
(all assume a stage), or "what worked before" (it worked in a background).

If something seems to need prior context, derive it from `v89/PRINCIPLE.md`
instead of looking it up.

## The principle, in one line

**Energy is never destroyed. It only changes mode. Space is one of the modes.**

Therefore: there is no background; matter is converted space; space curves
because energy is conserved and space is convertible; motion is regular
conversion rather than displacement; and there are no equations of motion in the
usual sense, because nothing has a trajectory.

## Standing constraints

1. Energy is never destroyed — the only law.
2. No background. Nothing may persist and be merely re-valued.
3. No imported field or species.
4. Emergent from the fabric alone.

## The unified law (synthesis, 2026-07-28)

One law table — `v89/battery/laws_V2g.cfg` — passes all 20 experiments with
only apparatus varying (`python3 v89/battery/battery.py --laws laws_V2g.cfg`).
The synthesis, in full in `v89/README.md` ("The law, as it stands"):
**amplitudes within a mode** (transport is never quantized); **atoms at mode
boundaries** (every conversion fires in whole action atoms ε = A₀ω/2π at the
source's pitch, two-atom credit; ħ-linearity measured to ~1e−8); **load
flattens pitch and flight is load** (x = (Em + flload)/cap; vacuum optics
linear, Kerr = loaded matter); **the choir's correction derived then
retired** (the restoring bias is the interference cross-flow of coherent
exchange — derived in `v89/s2/`, unnecessary at rate level post-flload;
kappa_reac=1 full-pass is an S2-full acceptance criterion); **the vacuum
skirt** (confirmed scored prediction); **momentum as the first moment of
conversion** (center-of-energy theorem replaces Noether; p1 gated; the ~100×
radiation-pressure deficit at the conversion door is recorded — another
S2-full criterion); **space flows: pressure pushes, nothing reaches out**
(s_k/s_disp: a mass maintains its extended graded depression — the
gravitational footprint g1; matter is emergently opaque — occultation g3;
g4 measured: a leaking blob's space flux is all mass-rate bookkeeping, no
steady monopole — the 1/r far field awaits a stable particle's internal
space cycle). Record: `v89/ROADMAP.md` §6–§7.

**The ratchet rule (v89/battery/README.md):** every kernel or law-table
modification runs the FULL battery before commit; experiments that pass and
defend a claim are added to the suite and gate all future modifications; bars
encode claims (sharpen by measurement, never soften to pass); a green test
leaves the suite only by explicit user decision.

## Everything below this line is HISTORICAL

The Standing Goal, stage table, Q-ball tooling, simulation-kernel policy,
parameter conventions and analysis instruments below describe the pre-v89
program. They are retained as a record of what was done and are **superseded as
a guide to work**. `scp_sim`, the seed generators, and the SFA analysis tools
all presuppose a fixed lattice and are not instruments for v89.

Use them only if explicitly asked to work on historical material.

---


## Standing Goal (HISTORICAL, pre-v89) — Carbon from fabric

**SUPERSEDED by v89.** Retained as a record. Primary program goal was: produce a **Carbon atom** (structural
analog) from the space-time fabric alone — the gauged complex Cosserat field
content — not by importing chemistry or external particle species.

Carbon is a **scale stack**. Work the stages in order; do not skip to
chemistry geometry or hand-waved shells. Until Stage 3 exists, the honest
product is fabric **carbon nuclei**, not atoms.

| Stage | Target | Status / notes |
|-------|--------|----------------|
| **0** | **Carbon mapping spec** — what \(Z\), \(A\), electron count, and observables count as success (binding, multipoles, shell radii). Map via conserved bookkeeping (\(Q\), \(Q_a\), \(\omega_a\), Gauss flux), not real-world chemistry. | **Done** — `v74/CARBON_MAP.md` |
| **1** | **Nucleon template** — stationary gauged Q-ball (fixed-\(Q\) / \(\eta\)-state when \(\eta>0\); ledger-closed; per-component clocks). Formation via Affleck–Dine condensate is the fabric-only path. | Largely done (v66–v73) |
| **2A** | **Liquid-drop carbon nucleus** — `v74/CARBON_MAP.md` + `v74/RESULTS.md`. **c6_light** Z-carbon parks (Q→650); **c12_light** A=12 stays super-critical (Q→1411>921). Both rendered. | **Done (primary)** |
| **2B** | **True multi-center carbon** (optional research) — retain \(A\) substructure: gauged interlock molecule, cold fission channel, flavored GDR. | After 2A; may need new binding |
| **3** | **Electron / multi-fabric sector** — light **opposite-charge** stable solitons; prove **positronium analog** first. Architecture: **v75 three-fabric (C/Q/L)** — see `v75/PROPOSAL.md`. Kernel change only with explicit user authorization. | **Blocking for atoms** — v75 design |
| **4** | **Atom** — C-nucleus + 6 light opposite charges, Coulomb-bound via \(A\); large box; force \(n=2\); multi-clock stability. | After 2 + 3 |
| **5** | **Spontaneous production** — condensate → fusion tree → recombination; abundance peak near carbon inventory. Fabric-only, not hand-placed. | After engineered seeds work |

**Success criteria (structural):** conserved-quantity bookkeeping at machine
precision; objects on the measured Q-ball / fusion branch; closed fabric
ledger where process-form applies; no fatal radiation. Parallels to real
carbon are structural, not quantitative (see CONCEPT.md).

**Default next work when the user does not specify otherwise:** **v86**
(`v86/THEORY_v86.md` is the canonical clean-state entry point, then
`v86/PART0_RESULTS.md` for what has actually been measured; program order in
Part C). The entire **zero-spend** program is done (2026-07-26): N1–N5, N9,
N10, HC-1, HC-2, HC-3 and EX-2 — instruments in `v86/n_battery/`. **Next rung:
N7, the inertia lock** (decides D5′; must precede any EX-1 boost work); cheap
CPU debt alongside it is **D7-lite** and a **gauged BdG**. v85 closed with the
atom-arc campaigns; the v86 council corrections govern. The stage table above
predates the v85 finding that 2B multi-center is load-bearing for atoms
(change pending explicit user sign-off).

**Tools (existing):** `radial_qball` / gauged shooter, `gen_qball_multi`,
`gen_qball_flavored`, `eta_qflow`, collision/pair seeds, SFA charge/flavor
views, `scp-runner` for local/remote runs. Do **not** modify `scp_sim` /
`sfa.h` unless the user explicitly requests kernel/format changes (Stage 3
may eventually need a design discussion).

---

## Document Roles — CRITICAL
- **CONCEPT.md** is the THEORY DOCUMENT — it presents the current best understanding
  of the physics as a cohesive, replicable description. Write it like a textbook chapter,
  NOT a lab notebook. Someone reading only CONCEPT.md should be able to understand and
  replicate the theory without knowing about the trial and error that got there.
  - Do NOT include chronological experiment narratives or failed approaches
  - DO include null results that clarify the theory ("it is NOT X, it IS Y")
  - DO include key numerical data that confirms claims
  - Present the physics as understood NOW, not as it was discovered over time
  - Use correct terminology (U(1) era, v66+): "ball"/"baryon" for the stable
    charged Q-ball, "component"/"quark" for the Φ_a sub-fields, "flavor" for the
    frequency/charge partition (ω_a, Q_a), "nucleus" for fused multi-ball droplets.
    "Braid", "oscillon", "UUD/UDD proton/neutron" are real-field-era (≤v53) terms —
    use only in historical context.
- **DISCOVERIES.md** is the chronological lab notebook — records what was found, when,
  including failed approaches and process. This is where the history lives.
- **EM_THEORY.md** is the detailed electromagnetic sector theory document (same
  standards as CONCEPT.md — cohesive, not chronological). NOTE: the photon analog is
  the gauge field A (v68+); the old θ-polariton mapping is historical.
- **FUTURE.md** tracks open questions and proposed experiments.
- **Version directories** (v28/ … v73/) contain per-experiment plans, results,
  analysis, and generated data. Current era: v66 (complexified U(1) Q-balls),
  v67–v68 (characterization + gauge design), v69 (gauged kernel-v3), v70 (skeptical
  verification + existence dynamics), v71 (quark substructure, nuclei, flavor),
  v72 (stationary η-soliton via fixed-Q relaxation — η-drain closed),
  v73 (process-form stability: fabric uptake/layment ledger; the spinning
  Q-ring — spin = winding of real-space circulation, L_z = nQ),
  v74 (carbon Stage 0–2A: Z-carbon parks; A=12 evaporates),
  v75 (multi-fabric proposal toward atoms — three fabrics C/Q/L first).

## Data Format Policy
- **ALL simulation output MUST use SFA format** (`.sfa` files via `sfa.h`)
- **NO `.bin` file output** — remove `save_field()` functions that write raw binary
- SFA header: `/home/d/code/scp/sfa/format/sfa.h` — single-header C library
- Usage: `#define SFA_IMPLEMENTATION` in exactly ONE .c file before including
- Link with: `-lzstd`
- Column convention: phi_x/phi_y/phi_z (SFA_POSITION), theta_x/theta_y/theta_z (SFA_ANGLE)
- All analysis tools must READ SFA files using `sfa_open()` / `sfa_read_frame()`
- SFA viewer: `sfa/volview` (binary `bin/volview`) — handles 12/24/30-column files;
  views: 4 field / 5 velocity / 6 accel / 8 U(1) gauge (|E|,|A|) / 9 charge (±ρ_Q) /
  0 flavor (per-component RGB + inline clock analysis) / C local clock;
  headless export: `volview -snapshot N -view K -out f.webp file.sfa`

## Simulation Parameters (standard)
- **Current standard (U(1) era, v66+)**: complexified 12-field Cosserat, optionally
  gauged. Config: `complex_phi=1`, `complex_gauge=1`, `g_gauge=0.05`, `m_theta=1.6`
  (closes the θ-drain), `eta=0` for particle experiments, absorbing BC.
  Vt(s) = (mu/2) s/(1+kappa s), s = Π|Φ_a|²; m²=2.25, mu=-41.345, kappa=50.
  Q-ball window ω ∈ (1.3087, 1.5) ungauged; (1.406, 1.5) and Q_max=921 at g=0.05.
  Seed via radial profiles (init=qball on CPU; gen_qball_* + init=sfa on GPU —
  the GPU kernel has no init=qball path; the init Gauss projection builds E).
- Real-kernel legacy defaults (≤v65): m^2=2.25, mu=-41.345, kappa=50, eta=0.5,
  A_bg=0.1, delta = {0, 3.0005, 4.4325} (from v28 optimization)

## Running Simulations — Use the MCP Runner

The `scp-runner` MCP server (source: `sfa/runner/`, binary: `bin/scp-runner`) manages ALL simulation
execution, both local and remote. It exposes `sim_*` tools via MCP that handle
instance management, building, running, monitoring, and downloading — with NO
sleep polling required. Use these tools instead of manual SSH/SCP/rsync.

**Instances are named and independent.** Every `sim_*` tool takes a `name`
parameter identifying the instance (chosen at `sim_setup`, e.g. "local1",
"gpu1"). Multiple named instances — local and remote, or several GPUs —
coexist and run concurrently; `sim_teardown(name=...)` destroys one without
touching the others. `sim_status` with no name lists all instances.

**Local execution** (CPU, for quick tests):
```
sim_setup(name="local1", executor="local")
sim_build(name="local1", sources=["sfa/sim/scp_sim.c"])   ← auto-detects gcc; quoted #includes auto-uploaded
sim_run(name="local1", config="N=64\nL=15\nT=5\n...", id="test_001")  ← config content, uses last-built binary
sim_run_status(name="local1", id="test_001")   ← instant, no polling
```

**Remote execution** (GPU, for production runs):
```
sim_setup(name="gpu1", executor="remote", disk_gb=40)  ← auto-provisions Vast.ai (reuses existing if running)
sim_build(name="gpu1", sources=["sfa/sim/scp_sim.cu"])  ← auto-detects nvcc; sfa.h etc. discovered from #include "..."
sim_run(name="gpu1", config="N=384\nL=100\nT=200\n...", id="gradient_test")
sim_run_status(name="gpu1", id="gradient_test")
sim_download(name="gpu1", remote_path="output.sfa", local_path="/space/scp/results/", wait=true)
sim_teardown(name="gpu1")        ← verifies downloads, destroys instance
```

**Remote with existing instance** (connect to already-running GPU):
```
sim_setup(name="gpu1", executor="remote", host="ssh5.vast.ai", port=12345)
```

**Custom GPU filter** (default is V100 16GB):
```
sim_setup(name="gpu1", executor="remote", gpu_filter={"gpu_name": "Tesla_V100"})
```
The `sim_setup` result includes `machine_id` and `dph_total` ($/hr). Machines
that fail provisioning are destroyed automatically and blacklisted for the
session.

**sim_run config**: Pass config file content (key=value lines). The runner writes it to
a `.cfg` file and invokes the last-built binary automatically. Alternatively, pass a
command string (e.g. `/path/to/binary /path/to/config.cfg`) for direct execution.

**sim_run extras**: `on_complete="<shell cmd>"` runs a local command (bash -c,
env `SCP_RUN_ID`/`SCP_STATE`/`SCP_INSTANCE`) when the run reaches a terminal
state. Local runs are queued nproc-aware (max concurrent = max(1, nproc/16),
override with env `SCP_RUNNER_MAX_LOCAL`; bypass per-run with `no_queue=true`);
queued runs show state/phase "queued" and start as slots free.

**sim_run_status** reports `phase` alongside `status`: `init` (process alive,
zero diag rows yet — GPU runs have a long silent CPU projection at startup, do
NOT kill an "init" run), `running`, `stalled` (alive but diag output older
than 10 min), `queued`, `complete`, `failed`. It also reports `proc_cpu_pct`,
`last_diag_age_s`, and `log_bytes`.

**sim_download**: default is async (poll `sim_download_status`); pass
`wait=true` to block until the transfer finishes and get
`{remote_bytes, local_bytes, verified}` back.

**File-based signals** (MCP push notifications can be missed — these files are
the reliable anchor): `~/.scp-runner/state/` holds per-instance dirs with
`<run_id>.status.json` (state/sim_time/total_time/last_update/wall_seconds,
refreshed every ~5 s), an empty `<run_id>.done` or `<run_id>.failed` marker
created atomically at terminal state (watch these with a file-watcher),
and `instance.json` ({reachable, degraded, last_contact, host, port,
machine_id}). `heartbeat.log` (one line/30 s) and `heartbeat.json` prove the
runner itself is alive — stale age means the runner died. `sim_status` also
reports `reachable`/`last_contact`/`degraded` per instance; a degraded
instance (3+ failed liveness probes) skips the SSH stat queries instead of
returning zeros.

**sim_build**: Omit `cmd` for auto-detection (gcc for .c, nvcc for .cu). If you need
a custom build command, use `${OUTPUT}` as the output path placeholder. Local
`#include "..."` headers are discovered recursively and added to the upload/hash
set automatically — passing just the main source file is enough.

**DO NOT** use manual SSH, SCP, rsync, or `sleep N` polling for simulation work.
The runner handles all of this internally with goroutines — every tool call
returns instantly with cached state.

**Key tools**: `sim_setup`, `sim_status`, `sim_build`, `sim_run`, `sim_run_status`,
`sim_run_cancel`, `sim_upload`, `sim_download`, `sim_download_status`,
`sim_list_files`, `sim_exec`, `sim_teardown`

## GPU Notes
- **V100-16GB**: Fits N=384 (10.3 GB). Use for standard proton/gradient tests.
- **V100-32GB**: Fits N=512 (19.3 GB). Use for large-grid or high-resolution runs.
- **RTX 4090**: 24 GB, sm_89, faster single-GPU but no multi-GPU.
- The CUDA kernel (`sfa/sim/scp_sim.cu`) uses async hooks for I/O — snapshots
  and diagnostics overlap with physics compute. GPU utilization should be ~100%.
- The `init=template` mode loads a small proton template (5 MB) and generates
  the background analytically — no large seed files needed.
- For gradient tests, use `bc_type=1` with `gradient_A_high`/`gradient_A_low`.

## Build Convention
- **`make -C sfa install`** builds everything and installs to `bin/` (18 binaries)
- **`make -C sfa runner`** builds just the MCP runner
- **`make -C sfa analysis`** builds just the analysis tools
- **`make -C sfa cuda`** builds CUDA kernels (requires nvcc)
- C with OpenMP: `gcc -O3 -march=native -fopenmp -o <binary> src/<file>.c -lzstd -lm`
- CUDA: `nvcc -O3 -arch=sm_70 -o scp_sim_cuda scp_sim.cu -lzstd -lm -lpthread`
- **ONE copy of sfa.h** at `/home/d/code/scp/sfa/format/sfa.h` — do NOT copy to other directories
- Include via relative path: `#include "../../sfa/format/sfa.h"` (adjust depth as needed)
- All binaries go to `bin/` via `make install` — do NOT commit binaries

## Physics Requirements — CRITICAL
- **ALL simulations MUST use the unified kernel with the full field content** —
  never a reduced/3-field-only simulator. The theta sector and (in complex mode)
  the imaginary sector are essential physics even when seeded to zero.
- Real mode: d²φ/dt² = ∇²φ - m²φ - V'(P) + η×curl(θ), d²θ/dt² = ∇²θ - m_θ²θ + η×curl(φ)
- Complex/gauged mode (v66/THEORY.md, v68/GAUGE_DESIGN.md, v69/SPEC.md): all matter
  derivatives link-covariant; charge conservation and the discrete Gauss law are
  exact by construction — verify `gauss_max` stays at the 1e-13 floor (a drift is
  implementation-bug tripwire #1).
- Verify (complex): 36 arrays (+ gauge blocks when complex_gauge=1), Q_phi conserved
  to the integrator floor at eta=0, g=0 byte-identical to the ungauged path.
- Reference docs: `v66/THEORY.md` (complexification), `v69/SPEC.md` (gauged lattice);
  historical real-field reference: `v34/torsion_coupling/src/v33_cosserat.c`

## Simulation Kernel Policy — CRITICAL
- **The unified sim kernel (`sfa/sim/scp_sim.c` and `sfa/sim/scp_sim.cu`) MUST NOT be modified unless the user EXPLICITLY requests changes to the simulation code.**
- "Explicitly" means the user says something like "modify scp_sim", "change the kernel", "update the simulator code", "add X to scp_sim". General requests like "run a simulation", "test this parameter", "try a new experiment" do NOT authorize modifying the kernel.
- To run experiments: write a config file and/or a seed generator. The kernel reads config — never modify it for a specific experiment.
- To add new initialization patterns: write a new seed generator in `sfa/seed/`, NOT by editing the kernel.
- The same protection applies to `sfa/format/sfa.h` — do NOT modify the SFA format without explicit user authorization.
- If a simulation requires physics not supported by the current kernel, ASK the user before modifying it. Describe what change is needed and why.

## Simulation Kernel Location
- CPU kernel: `sfa/sim/scp_sim.c` — build with `gcc -O3 -march=native -fopenmp -o scp_sim scp_sim.c -lzstd -lm`
- GPU kernel: `sfa/sim/scp_sim.cu` — build with `nvcc -O3 -arch=sm_70 -o scp_sim_cuda scp_sim.cu -lzstd -lm`
- Seed generators (U(1) era): `gen_qball_pair.c` (two balls, phases/signs),
  `gen_qball_boost.c` (boosted single ball; v=0 → single static ball),
  `gen_qball_bath.c` (ball + θ bath), `gen_qball_multi.c` (N balls — nuclei),
  `gen_qball_quark.c` (single-COMPONENT lumps, mask + per-component centers),
  `gen_qball_flavored.c` (per-component profiles/frequencies; multi-ball;
  accepts 2-col symmetric or 4-col flavored profiles),
  `eta_qflow.c` (fixed-Q relaxer: the stationary η-coupled (Φ,Θ,ω) state —
  REQUIRED for any η>0 particle experiment; Θ=0 seeds shed 1.4–5.6%),
  `radial_eta_soliton.c` (linear-BVP Θ partner). Profiles from
  `radial_qball` (ungauged) or the v69 gauged shooter; profile f(r) is the
  PER-COMPONENT amplitude (v72 convention fix). NOTE: seed writers MUST
  call `sfa_finalize_header()` before `sfa_write_frame()`, and SFA semantics
  are SFA_POSITION=0/ANGLE=1/VELOCITY=2 (numeric {1,3,2} was the v72-found
  scrambled-seed bug).
- Seed generators (real-field era, historical): `gen_braid.c`, `gen_oscillon.c`,
  `gen_phase_confined.c` (UUD/UDD), `gen_proton_analytical.c`, `gen_composite.c`
- **Analytical seed warning**: `gen_proton_analytical` (Level 2) produces baryons
  with ~4× the equilibrium binding density (P_int ~1270 vs ~320 per baryon in V42).
  These seeds need T=2000+ to relax and are NOT suitable for energy comparisons.
  For energy-sensitive experiments (mass defect, binding), use template seeding
  (`init=template` with `proton_template.sfa`) or the `gen_deuterium.c` generator
  which produces better-equilibrated initial conditions.
- Pre-converged templates: `v43/proton_formation/proton_template.sfa` (64³, UUD),
  `v43/proton_formation/neutron_template.sfa` (192³, UDD from V41)
- Config files: `sfa/sim/*.cfg`
- SFA header: `sfa/format/sfa.h` (single copy, include via relative path `../format/sfa.h`)
- Reference implementation (historical): `v34/torsion_coupling/src/v33_cosserat.c`

## Data Transfer and Storage
- The `scp-runner` handles all remote file transfers via `sim_download`.
  It uses rsync with `--append-verify` internally — no manual rsync needed.
- `sim_teardown` verifies downloads before destroying instances.
- Large SFA files should go to `/space/scp/` (separate disk, 600+ GB free).
- Local working disk (`/home/d/code/scp/`) has limited space — keep only
  small files (diag.tsv, templates, analysis results) on the local disk.

## SFA Archival (rclone)
- Completed SFA files can be archived to cloud storage via rclone:
  ```bash
  rclone copy local_file.sfa scpsfa:scpsfa/v42/
  rclone ls scpsfa:scpsfa/           # list archived files
  rclone copy scpsfa:scpsfa/v42/file.sfa ./  # retrieve when needed
  ```
- The `scpsfa` rclone remote is pre-configured locally for B2 cloud storage.
- Archive SFA files after analysis is complete to free local disk space.
- Keep diag.tsv, analysis.json, and freq.json locally (small files, always needed).
- Large f32 output SFAs (10-30 GB) should be archived and deleted locally.
  The f16 viewing copies (1-9 GB) can be kept locally or also archived.
- To check if a file is already archived: `rclone ls scpsfa:scpsfa/`

## Diagnostics Requirements
- Every simulation MUST include fragmentation detection (connected component analysis)
- Track per-cluster mass, centroid, aspect ratio — not just global totals
- Time-averaged death check (rolling window, not instantaneous)
- Absorbing boundary damping for compact objects (not periodic-only)
