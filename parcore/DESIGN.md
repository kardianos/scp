# parcore — a generic, pluggable, accelerated parallel core for free-cell simulations

`parcore/` generalizes the structurally-isolated `pfc/` experiment into a
reusable Go core. `pfc` proved the structure (contiguous AoS cells with no
pointers, int32 index topology, two-phase publish/compute/apply barriers,
edge-color rounds for dual-endpoint hops, z-slab ownership, per-shard ghost
halos, one-goroutine-per-shard procedural runtime, zero-alloc step). `parcore`
makes the **cell type, slot type, and pass list generic** so a new simulation
plugs in without touching the harness.

This is a self-contained engineering prototype. It does **not** touch `pfc/`,
`v93/fab/`, `v93/kernel/`, or `v93/runs/` (those trees are read-only
references). It is a standalone module (`GOWORK=off go build|test|run ./...`),
not added to the repo `go.work`.

---

## 1. The driver / core contract

### 1.1 The split

| Owned by | What |
|----------|------|
| **Driver** (user) | the per-cell value type `C`, the per-slot type `S`, the immutable law, and the ordered `[]Phase` that composes one timestep |
| **Core** (`parcore`) | sharding into contiguous AoS shards, the publish/compute/apply barriers, ghost-halo exchange, edge-coloring, the CSR topology, the procedural worker pool, timers |

The driver never sees a lock, a goroutine, a shard boundary, or another shard's
state. It writes directly on its shard's `[]C` via `&v.Cells[i]`.

### 1.2 The types (the public contract)

```go
// A Phase is one step in the timestep graph. Consecutive phases are separated
// by a barrier. The driver composes the graph; the core enforces the sync.
type Phase[C any, S any] struct {
    Name string
    Kind PhaseKind   // KindParallel (per-shard) or KindBoundary (serial, once)
    Run  func(v *View[C, S])        // KindParallel: executes on each shard
    Bnd  func(bv *BoundaryView[C])  // KindBoundary: executes once on caller
}

type Driver[C any, S any] interface {
    Step() []Phase[C, S]   // called ONCE at build; replayed every Step
}

// Optional: size per-shard scratch at build (keeps Step zero-alloc).
type ScratchDriver[C any, S any] interface {
    Driver[C, S]
    InitShard(nOwn, nGhost, nSlot int) any
}
```

`View[C,S]` is what a `KindParallel` closure operates on. It carries the
shard's owned+ghost `[]C`, the `[]S` channel data, the core-owned endpoint /
CSR / color topology (`SI`, `SJ`, `Cls`, `Inc`, `Colors`, `Live`), the per-slot
hop-angle scratch (`SlotParam`), the boundary-touch list (`BndTouch`), and the
driver's `Scratch`. Exchange is two methods: `v.Publish()` (face → global
buffer) and `v.Pull()` (global buffer → ghosts).

### 1.3 The value-struct driver pattern

The driver is a plain struct holding the immutable law (e.g.
`type Driver struct { Law Law }`). It is passed to `FromWorld` as a
`Driver[C,S]` interface value; the core calls `Step()` once, caches the
`[]Phase`, and replays it. Because the phase closures capture only the
immutable driver, sharing them across shards is safe — **there is no shared
mutable state on the driver**. All per-shard mutable state is either inside the
cell value `C` or in `View.Scratch` (sized per-shard by `InitShard`).

To plug in a new simulation:

1. Define `C` (cell value) and `S` (slot value) — no pointers inside either.
2. Write pure pass functions taking plain slices (`func PassPitch(law, []C, a, b int)`).
3. Write a `Driver` whose `Step()` returns those passes wired into `[]Phase`,
   with `PublishFaces`/`PullGhosts`/`BoundaryHops` placed at the halo-refresh
   points the physics requires.
4. Build a `World[C,S]` (the `BuildLattice` helper + a cell-seed callback does
   the common case) and call `parcore.FromWorld(w, d, nWorkers)`.

The bundled `example/` driver does exactly this. Its pass functions are PURE —
the same functions are called (a) directly on flat arrays by the serial
baseline (`SerialSim.Step`) and (b) from the parallel `Phase` closures. This is
the proof that the abstraction imposes no physics-side overhead and that the
driver is reusable without the core.

---

## 2. The parallel model

### 2.1 Sharding (z-slabs, contiguous AoS)

`FromWorld` partitions cells into `Workers` z-slabs (whole xy-planes when the
world is a cubic lattice — minimal cut-plane boundary; falls back to contiguous
id ranges otherwise). Each shard gets its own `make([]C)` and `make([]S)` —
**distinct backing arrays** (verified by `TestDistinctSlabs`). Layout per shard:

```
Cells [0:NOwn)            owned (global G0..G0+NOwn-1)
Cells [NOwn:NOwn+NGhost)  ghosts (frozen copies of remote owned cells)
```

Slot endpoints `SI/SJ` are **local** cell indices. The core owns them (they are
NOT inside the user's `S`), which decouples topology from channel physics and
makes the generic core's CSR/color/build job trivial.

### 2.2 The two-phase publish / compute / apply cycle

Every consecutive `Phase` pair is separated by a barrier (`runPhase` returns
only when every shard finished the phase). A halo refresh is just a
`PublishFaces` half-phase followed (across the barrier) by a `PullGhosts`
half-phase — but the bundled driver **fuses** the publishes into the compute
phases (each compute phase ends with `v.Publish()`), so a full timestep is:

```
phase 0  start     Pull + Pitch + Geom + Wants + Scale + Publish(Sscl)  → barrier
phase 1  mid       Pull + Apply + Pitch + Precess + Sum + HopWeight +
                   InteriorHops + Publish(Fa)                          → barrier
phase 2  boundary  BoundaryHops over pub (serial, caller thread)
phase 3  finish    MergeFaBack + Commit + Publish(end)                  → barrier
```

3 barriers/step — pfc's fused-wave target. The `Pull` at the top of phase 0
reads the buffer phase 3 of the previous step published (settled by the
step-boundary barrier); the `Pull` at the top of phase 1 reads the `Sscl` phase
0 just published. Both correct because a full-cell `Publish` carries every
field.

### 2.3 Edge-coloring for dual-endpoint updates

Dual-endpoint slot updates (unitary hops) read AND write both endpoints, so two
slots sharing a cell cannot update concurrently. The core **edge-colors the
interior slots** (greedy, flat `used` table — not one `[]bool` per cell, which
was pfc's catastrophic-allocation bug) so each color bucket is
endpoint-disjoint. A `KindParallel` phase iterates `for _, bucket := range
v.Colors { for _, s := range bucket { ... } }`; within a bucket the writes are
race-free and the pairwise conservation holds exactly.

### 2.4 Cross-shard boundary hops (the thin cut-plane)

Slots with endpoints on two different shards cannot be edge-colored locally.
The core collects them into a global boundary table `(gi, gj, shard, slot)`.
After the interior rounds, a single `KindBoundary` phase runs **serially on the
caller thread** over the global publish buffer, applying the same dual-endpoint
update the interior rounds used. Running it single-threaded avoids a per-edge
region mutex on the dense contact graph; boundary edges are a thin set (the
cut-plane: ~1 xy-plane of cells per slab boundary). A following `KindParallel`
phase merges the published result back into owned face cells (`v.MergeBndTouch()`).

### 2.5 The procedural worker runtime

At `FromWorld` (when `Workers>1`) the core starts **one long-lived goroutine
per shard**, bound for life. `Step` never spawns; it posts a phase index
(`curPhase`) to a condition variable and barriers. There is no work queue, no
`go` per phase, no pool abstraction. The worker reads `c.curPhase` under the
mutex and calls `c.phases[curPhase].Run(&c.shards[id].view)` — **no closure is
allocated per step** (this is why `TestStepZeroAlloc` passes at 0 allocs).

---

## 3. What is abstracted vs exposed

| Abstracted away (driver does not see) | Exposed to the driver (on `View`) |
|---|---|
| sharding, z-slab partition | `G0`, `NOwn`, `NGhost` |
| per-shard backing arrays | `Cells []C`, `Slots []S` |
| ghost halo internals | `Publish()`, `Pull()`, `MergeBndTouch()`, `PublishedCell(li)` |
| CSR build/rebuild | `Cls`, `Inc`, `SI`, `SJ` (read-only) |
| edge coloring | `Colors [][]int32` (interior buckets, read-only) |
| live-slot management | `Live []int32`, `SetLive(s)`, `ResetLive()` |
| boundary edge table | `BndTouch []int32` + `SlotParam` (driver writes hop angles) |
| worker pool, barriers, cond-var | (none — invisible) |
| phase timers | `EnablePhaseTimers`, `PhaseTimerReport` |

The driver controls **granularity** (how many phases, how much fused into each)
and **where halos refresh** (by inserting `Publish`/`Pull`). The core controls
**how** the phases run in parallel and stay consistent.

---

## 4. The dynamic-topology seam (named, deferred)

The free-cell programme requires that cells be born and die (no background);
the roster/graph is not permanent. `parcore` v1 **freezes the graph** at
`FromWorld` (the `frozen` topology policy). The seam is named explicitly:

```go
type Topology interface {
    Dynamic() bool   // frozen returns false; the core skips refresh work
}
```

A future `DynamicTopology` would implement refresh: (a) spatial re-neighboring
from current positions, (b) birth/death of cells into a freelist, (c)
open-addressing pair→slot hash rebuild, (d) canonical CSR rebuild, (e)
re-edge-coloring — all per shard, then a cross-shard ghost-set reconciliation
at a barrier. This is `pfc/README` pushback item 1 verbatim: dynamic
`topo_refresh` across shards is the hard parallel-systems problem and is
deferred on purpose. v1 still lets **geometry** (positions, lens areas, hop
weights) vary every step — the driver's geom pass recomputes those; only the
**roster** (which cells exist, which slots connect them) is frozen.

Note: even with a frozen roster, slot liveness (lens area > 0) varies every
step as cells move; the `Live` list + `ResetLive` capture that.

---

## 5. The conservation contract

The core makes the free-cell invariants **hard to break by accident**:

1. **Contiguous value cells, no `*C` fields.** All mutation is `&Cells[i]` on a
   slice backing array; there is no background of persistent re-valued
   pointers. The driver literally cannot hold a stale cross-shard reference.

2. **Edge-coloring of dual-endpoint slot ops.** A colored round guarantees
   endpoint-disjoint slots, so a pairwise rotation is race-free and its
   two-cell norm preservation is exact (not probabilistic).

3. **`Givens2` / `Givens2Buf` — the blessed pairwise unitary primitive.** This
   is the v93/freecell pass-F (field) and pass-U (unitary dense channel)
   rotation. Drivers SHOULD compose dual-endpoint updates from it so
   conservation is a theorem of the update:
   ```
   |a'i|² + |a'j|² == |ai|² + |aj|²   (exact, by construction)
   ```
   `TestGivens2NormConservation` checks this to 1e-15. The bundled driver's
   `ApplyHop` and `ApplyHopPub` both build on `Givens2`.

4. **Full ledger check.** `TestParallelEnergyConservation` confirms total
   `E = Es+Em+Ee` drift is **0.000e+00 over 30 steps** at 1/2/4/8 workers.
   Space transport conserves `Es` pairwise (src debits exactly what dst
   credits, scaled by `Sscl`); the Givens hops conserve `Ee = Fa1²+Fa2²`
   pairwise; so the sum is conserved by construction. Multi-worker runs need
   not be bit-identical to 1-worker (the cross-shard hop ordering differs),
   but the ledger is conserved at every worker count.

---

## 6. What is NOT supported (v1 seams)

| Not supported | Why / seam |
|---|---|
| **Dynamic topology** | roster frozen at `FromWorld`; the `Topology` interface names the seam (§4) |
| **Selective (per-field) ghost exchange** | v1 publishes **whole cells**. pfc's `Fa1Next`-only and `Cr`-only exchanges are the bandwidth optimization; the hook is a future `ExchangeField(name, get, set)`. Full-cell exchange is correct, just bandwidth-heavier (see §8 bottlenecks). |
| **Hot/cold cell split** | v1 uses a single AoS slab (the demo cell is 88 B). pfc's `CellLo`(112B hot)/`CellHi`(96B cold) split lets Lo-only passes stream less. The generic core can support a `[2]C` split later; documented as the #1 perf lever. |
| **Full v93 law parity / C≡Go** | the bundled driver is a freecell-**shaped** demo (pitch + space transport + field precession + unitary hops), NOT a port of the law table. It exists to exercise the machinery. |
| **Bit-identical serial Givens ordering** | parallel colored application is physically related to the serial canonical sweep but not float-identical; the ledger is conserved either way (proven by the tests). |
| **GPU / SIMD** | CPU-first; Go 1.27 `simd` not wired. |

---

## 7. File layout

```
parcore/
  go.mod                  module scp/parcore (standalone; GOWORK=off)
  core.go                 Phase/PhaseKind/Driver/World/Cluster types (the contract)
  view.go                 View, BoundaryView, ShardView; Publish/Pull/MergeBndTouch
  topo.go                 Topology seam, BuildLattice, CSR, edge-color, partition
  shard.go                shard struct, FromWorld partitioner, ghost/boundary tables
  workers.go              procedural worker runtime (one goroutine/shard, cond-var barrier)
  step.go                 Step: replay phase graph with barriers + boundary round
  phases.go               PublishFaces / PullGhosts / BoundaryHops constructors
  hop.go                  Givens2 / Givens2Buf (conservation primitive) + LensArea
  timers.go               coarse phase timers
  topo_test.go            edge-color disjointness, CSR, partition, Givens2 conservation
  cluster_test.go         parallel conservation, distinct slabs, zero-alloc, no-deadlock
  example/                the bundled driver (a real free-cell-shaped sim)
    cell.go               Cell, Slot, Law (AoS, no pointers)
    passes.go             PURE pass functions (shared by serial + parallel)
    driver.go             Driver.Step() — the 4-phase timestep graph
    serial.go             SerialSim baseline (same passes, flat arrays)
    driver_test.go        serial conservation, zero-alloc, pitch monotonicity
  cmd/bench/
    main.go               scaling sweep: parcore vs serial, strong + weak, phase timers
    profile.go            -cpuprofile flag helper
```

The core never imports the driver; the driver imports the core. `example/`
proves the public API is sufficient to express a real sim.

---

## 8. Benchmark numbers (Ryzen 5700G, 16 logical cores, GOMAXPROCS=workers)

Raw sweep (`GOWORK=off go run ./cmd/bench -n 8000,64000,216000,512000 -steps 30
-workers 1,2,4,8,16`). Cell = 88 B single AoS slab; Slot = 48 B.

**Speedup vs 1 worker (strong scaling):**

| NC \ W | 1 | 2 | 4 | 8 | 16 |
|---:|---:|---:|---:|---:|---:|
| 8 000 | 1.00× | 1.49× | 1.84× | **2.04×** | 1.86× |
| 64 000 | 1.00× | 1.43× | 2.18× | **2.20×** | 2.01× |
| 216 000 | 1.00× | 1.39× | **1.56×** | 1.50× | 1.52× |
| 512 000 | 1.00× | 1.34× | **1.48×** | 1.44× | 1.37× |

**Cell throughput per proc (Mcell/s/proc) — the weak-scaling / bandwidth signal:**

| NC \ W | 1 | 2 | 4 | 8 | 16 |
|---:|---:|---:|---:|---:|---:|
| 8 000 | 5.32 | 3.97 | 2.45 | 1.36 | 0.62 |
| 64 000 | 3.52 | 2.52 | 1.92 | 0.97 | 0.44 |
| 216 000 | 2.99 | 2.08 | 1.17 | 0.56 | 0.28 |
| 512 000 | 2.97 | 2.00 | 1.10 | 0.54 | 0.25 |

**Framework overhead:** parcore at 1 worker matches the serial baseline within
~7% (e.g. 64 000 cells: serial 53.1 steps/s, par 1W 55.1 steps/s — the 1W path
has no ghosts and no barriers, so the only cost is phase indirection). The
parallel path adds real speedup with no correctness tax.

### Reading the numbers

- **At the cache-friendly 64 000 cells, ~2.2× at 8 workers** — this matches
  `pfc`'s documented local-Ryzen scaling (pfc 64k/8w ≈ 2.8×; the gap is the
  hot/cold split, see §6).
- **Past ~216 000 cells, scaling flattens to ~1.5×** — the **memory-bandwidth
  wall**. A CPU profile at 512k/8w shows 581% CPU utilization (5.8 of 8 cores
  busy) with **no lock contention**: the work IS distributed, but per-core
  throughput drops as shards share the memory bus (Mcell/s/proc falls from
  2.97 at 1W to 0.54 at 8W). This is exactly the regime `pfc/README`
  documents: "After ~0.5M cells, aggregate performance declines gently — no
  cliff." A machine with more memory channels (EPYC) would scale further, as
  it does for pfc (~3.3× at 16W on large NC).

### Where the bottlenecks are (ranked)

1. **Memory bandwidth at large NC** (#1). The single 88 B slab streams the
   whole cell — including cold fields (`Sscl`, `Fa1`/`Fa2`) — through every
   pass. The hot/cold split (pfc's `CellLo`/`CellHi`) is the largest available
   win; a `[2]C` generic split is the natural extension.
2. **Full-cell ghost exchange.** v1 publishes whole cells where pfc publishes
   only the changed field (`Cr`, or `Fa1Next`/`Fa2Next`). Selective exchange
   (a driver-supplied pack/unpack) would cut halo bandwidth several-fold.
3. **Barrier count (3/step).** Already at pfc's fused target; further fusion
   would merge the `finish` Publish into the next step's `start`, but that
   crosses the step boundary (a different sync structure).
4. **`PassGeom` full-slot scan.** Refreshes geometry over ALL slots each step
   (Sqrt + lens math). pfc's `Live` list skips dead slots in the downstream
   passes, which the demo already does; the geom scan itself is inherent.
5. **Allocator (init only).** `FromWorld` does the sharded `make`s, ghost
   maps, and per-shard `g2l`/`seen` maps (build-only, same class as pfc's
   `FromWorld` jagged allocs). `Step` itself is verified 0-alloc.

---

## 9. C / OpenMP note (assessment, not implemented)

The same sharded two-phase model maps cleanly to OpenMP for the C kernel
(`v93/kernel/freecell.c`). The question is which of the C step's passes are
**parallelizable** vs **serial-by-design**, and the expected C-vs-Go gap.

### 9.1 Which passes parallelize

The C step (passes 0, S, 1, G0, 2, U, 3, 4, 5, H0, H0b, H0c, H, D, F, 6, 7 —
see `v93/fab/step.go` for the Go mirror) splits into three classes:

| Class | Passes | Parallelizable? | How |
|---|---|---|---|
| **Local per-cell** | 0 (flload), 1 (pitch/cr/w), F-precess, 6-clock, 7 (plane tumble) | **Yes — embarrassingly** | `#pragma omp parallel for` over cells. Pure local writes. |
| **Slot / gather (CSR)** | S (space wants), 2 (lens+dense wants), 3 (outflow limiter), F-sum, F-hopweight | **Yes — with two-phase buffering** | Each pass writes a Jacobi "want"/"next" buffer, applied in a separate per-cell pass. This is exactly what freecell already does (it buffers wants/space-flux/bond-kicks). |
| **Dual-endpoint unitary** | U (unitary dense hops), F-hops, D-motion bond kicks | **Yes — with edge-coloring** | Color slots so endpoints are disjoint within a color; `#pragma omp parallel for` over a color's slots, color-sequential. (freecell applies these in a canonical serial Givens order today; the parallel colored form is *physically related* but not FP-identical — same caveat as parcore.) |
| **Serial-by-design** | **pass 6 (the door / beat-gated conversion)**, identity ledger (H0b gid hysteresis), atoms credit accrual | **No — and should not be** | The door fires whole action atoms `ε = A₀ω/2π` at the source's pitch with a two-atom credit; it is a per-cell sequential decision that reads accumulated credit and writes `Em`/`Ee`/`Es`/`qcnv` atomically. Parallelizing it would race the credit. It is also where the conservation *ledger* is closed per cell. This is the irreducible serial core. |

So the parallelizable fraction is high: the per-cell and slot-gather passes
dominate the wall time (pitch, geom, wants, hops), and they all admit the
sharded two-phase form. The serial-by-design fraction (the door + identity
ledger) is a per-cell loop with no cross-cell dependency, so even it can run
as `#pragma omp parallel for` over cells — it is "serial" only in the sense
that each cell's door decision is sequential, not that the loop is serial.

### 9.2 The OpenMP mapping

```c
// Pseudocode: the sharded two-phase step in OpenMP
step() {
    // shard = thread owns Cells[g0:g0+nOwn] + ghosts (domain decomposition)
    PHASE(pull_ghosts);          // halo exchange (MPI or shmem between threads)
    #pragma omp parallel
    {
        // local per-cell passes — fully parallel
        pass_pitch();            // per owned cell
        pass_geom();             // per slot (read-only ghost)
        pass_space_wants();      // Jacobi-buffered
        pass_space_scale_apply();// per cell, reads buffered wants
        pass_field_precess();    // per cell
        // dual-endpoint hops — colored, color-sequential, parallel within a color
        for (c = 0; c < nColor; c++)
            #pragma omp for
            for (s in color[c]) apply_uhop(s);   // disjoint endpoints
        // the door — per cell, parallel (each cell's door is independent)
        #pragma omp for
        for (i in owned) pass_door(i);           // serial-by-design PER CELL
    }
    PHASE(boundary_hops);        // cut-plane edges, single thread or reduced
    PHASE(publish_faces);        // halo for next step
}
```

The hard part is **not** the OpenMP pragmas — it is the **domain decomposition
+ ghost halo** machinery, which is what `parcore` (and `pfc`) already solved in
Go. Porting that machinery to C with MPI or with OpenMP shmem shards is the
engineering work; the step-internal parallelism falls out.

### 9.3 Expected C-vs-Go speed gap

| Factor | Go (`parcore`) | C (`freecell.c`, `-O3 -fopenmp`) | Ratio |
|---|---|---|---|
| Per-cell arithmetic (pitch `Cbrt`, `Sincos`, etc.) | Go `math` (libm calls, no inlining across boundaries) | `-O3` inlines + vectorizes libm | ~2–4× faster in C |
| Memory bandwidth | same hardware limit | same hardware limit | ~1× (bandwidth-bound at large NC anyway) |
| Barrier cost | Go cond-var (~1–2 µs) | `#pragma omp barrier` (~0.1–0.3 µs) | C cheaper, but not the bottleneck |
| AoS locality | Go GC layout is dense | C struct is dense | ~1× |

**Realistic expectation: a C/OpenMP port of the same sharded two-phase model
would be ~2–4× faster than the Go core on the per-cell arithmetic**, but would
hit the **same memory-bandwidth wall** at large NC (the bandwitch limit is
hardware, not language). At small NC (cache-resident), the C advantage is the
full ~2–4×; at large NC (bandwidth-bound), the advantage shrinks toward ~1.5×
because both are waiting on DRAM.

The C kernel is currently **serial by design at the step level** (`-O3
-fopenmp` is linked but the freecell step is a single-threaded canonical
sweep — the `#pragma omp` would need the pass restructuring above to do
anything). So the comparison today is **parallel-Go vs serial-C**: parcore at
8 workers (~2× over its own 1-worker) vs serial C at 1 worker (but ~2–4×
faster per worker). The net is that **parallel Go and serial C are in the same
ballpark wall-clock**, and a parallel C port would beat both by the product
(~2× parallel × ~2–4× per-core).

### 9.4 Recommendation

The Go core (`parcore`/`pfc`) is the right vehicle to **validate the parallel
model and the driver/core contract** — generics make the pluggability clean
and the iteration fast. Once the model is proven (this prototype), porting the
**machinery** (sharding, ghost halo, edge-color, two-phase barriers) to C with
OpenMP/MPI is mechanical, and the per-cell passes port 1:1 from `freecell.c`'s
existing serial forms. The serial-by-design door stays per-cell-parallel (it is
already structured that way). Don't port the Go to C line-for-line; port the
**structure** and reuse the C step's pass bodies.

---

## 10. Verification

```
$ cd parcore && GOWORK=off go build ./...    # clean
$ cd parcore && GOWORK=off go vet ./...      # clean
$ cd parcore && GOWORK=off go test ./...     # all pass
```

Key test results:
- `TestParallelEnergyConservation`: **0.000e+00 drift** over 30 steps at 1/2/4/8 workers.
- `TestStepZeroAlloc`: **0 allocs/run** at 1 and 4 workers.
- `TestDistinctSlabs`: distinct `[]C` backing arrays per shard.
- `TestGivens2NormConservation`: pairwise norm conserved to 1e-15.
- `TestParallelVsSerialEnergy`: 1-worker parallel matches serial-baseline conservation.

---

## 11. Verdict — is this acceptable to adopt?

**Yes, as a v1 generic core / parallel-model validator. No, as a drop-in
replacement for the serial C kernel today.**

What it proves:
- The **driver/core contract is clean and sufficient**: a real free-cell-shaped
  sim (pitch, space transport, field precession, unitary hops, cross-shard
  boundary hops) plugs in as pure pass functions + a 4-phase graph, with zero
  physics-side overhead vs the serial baseline of the same passes.
- The **parallel model is real**: 2.0–2.2× at 8 workers on cache-friendly sizes,
  zero allocations on the step path, conservation to machine precision at every
  worker count, no lock contention. This generalizes pfc's proven structure
  into a reusable, type-parameterized core.
- The **seams are honest**: dynamic topology, selective exchange, hot/cold
  split, law parity, and bit-identical ordering are all named with concrete
  hooks, not papered over.

What's missing before it could carry production campaigns:
1. **Hot/cold cell split** (the #1 perf lever — closes most of the gap to pfc).
2. **Selective ghost exchange** (the #2 lever — cuts halo bandwidth).
3. **A real driver**: port the v93 law table's passes into a `parcore` driver
   and A/B against `freecell.c` (the C≡Go ratchet test).
4. **Dynamic topology** (the hard parallel-systems problem; the seam is named).
5. **C/OpenMP port** of the machinery (§9) for the ~2–4× per-core headroom —
   only worth it once the Go model is proven against the law.

The recommended adoption path: keep `parcore` as the **parallel-model
reference and driver sandbox**, port the v93 law into an `parcore` driver to
run the A/B, and treat the C/OpenMP port as the production target once the
model clears the ratchet.
