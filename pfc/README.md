# pfc — parallel free-cell experiment

Local experiment (not on the v93 ratchet): can a **sharded, two-phase**
free-cell step with **contiguous AoS cell storage** scale past the serial
C/Go freecell kernels?

Base physics is **shaped like** the most recent freecell (`v93/kernel/freecell.c`,
`v93/fab/`), not a line-for-line port and not bit-identical to the law table.

## What this is testing

| Idea | Status in this tree |
|------|---------------------|
| `[]Cell` contiguous AoS, **no pointers inside `Cell`** | yes |
| Slots as index pairs (`int32`), CSR of indices | yes |
| Core physics **structurally isolated** from the parallel harness | `core/` vs `par/` |
| Two-phase (publish / compute / apply) with phase barriers | yes |
| Per-proc continuous owned slices + ghost halo | yes |
| Cross-proc dual-endpoint hops without full-graph mutex | edge-color rounds |
| Full v93 law (all passes, atoms, cantus, amp_nat, …) | **not yet** |
| Dynamic topology every step (birth/death freelist) | **frozen graph** for now |
| C≡Go / FP-ledger A/B vs freecell | **out of scope** |

## Pushback (what is / is not realistic)

**Realistic and worth doing**

- Contiguous AoS + index-only topology: good locality when a pass touches
  many fields of the same cell (the freecell pattern).
- Two-phase / Jacobi: freecell already buffers wants, space flux, bond
  kicks. Making *published* state read-only for a phase removes most locks.
- Sharding by slab (or by contiguous global-id ranges) with **ghost copies**
  refreshed at barriers: standard domain decomposition; memory stays local.
- Edge-coloring dual-endpoint unitary hops: degree‑bounded graph ⇒ few
  colors; within a color, endpoints are disjoint ⇒ no cell lock.

**Not realistic as a first cut (deferred on purpose)**

1. **Dynamic `topo_refresh` across shards** — freelist, open-address pair
   hash, cross-proc birth/death, canonical CSR rebuild. This is the hard
   parallel systems problem. Start with `freeze_geo`-style static links.
2. **Bit-identical serial Givens order** — pass F/U apply rotations in a
   global canonical order. Parallel Jacobi double-buffer is *physically
   related* but not FP-identical. Conserving the ledger to 1e‑13 under
   reordering needs a separate proof/port.
3. **Per-hop `sync.RWMutex` / CAS on cells** — fine for correctness toys,
   kills speedup under dense degree‑15 contact graphs. Prefer barriers +
   coloring; use a **region mutex only** if a boundary edge cannot be
   colored (should not happen if coloring is global).
4. **GPU drop-in** — same objections as the earlier hardware note; this
   experiment is CPU-first.
5. **Generics/`interface` cell processors** — deferred; packages isolate
   core vs harness without them.

If the only success criterion were “full v93 battery, bit-identical, on
1000× cells tomorrow,” this design would be the wrong vehicle. If the
criterion is “can two-phase + contiguous shards + edge-color hops show
multi-core speedup on freecell-*shaped* work?” — that is what we measure.

## Layout

```
pfc/
  core/     CellLo + CellHi (split), Slot, World, serial Step
  par/      Shard (owned Lo/Hi), one goroutine per shard, Step
  cmd/bench scaling microbench
  cmd/smoke correctness + prove distinct Lo pointers per shard
  cmd/profile  phase timers + pprof
```

`core` never imports `par`. **No FCS writer in this tree yet** — see
[FCS / snapshot integration](#fcs--snapshot-integration-future) below.

### Ownership model (procedural workers, not a pool type)

- Worker `id` is bound for life to `shards[id]`.
- Each shard has its own `make([]CellLo)` and `make([]CellHi)` — **different
  backing arrays** (smoke prints distinct `Lo=` pointers).
- Layout per shard: `[0:NOwn)` owned, `[NOwn:NOwn+NGhost)` ghosts.
- `FromWorld` starts **one goroutine per shard**; `Step` only posts a phase
  id and barriers (`par/workers.go`). No work queue, no `go` per phase.

### CellLo / CellHi split

| Slab | Size | Typical passes |
|------|-----:|----------------|
| `CellLo` | **112 B** | pitch, space, field fa, clocks (hot) |
| `CellHi` | **96 B** | planes, FaNext/EsNext, forces (scratch / cold) |

Index `i` names the same site in both slabs. A Lo-only pass does not stream Hi.

## Build / run

```bash
cd pfc
GOWORK=off go test ./...          # or add ./pfc to the repo go.work
GOWORK=off go run ./cmd/smoke
GOWORK=off go run ./cmd/bench -n 64000 -steps 20 -workers 1,2,4,8,16
```

## Scaling findings (Vast EPYC campaigns, 2026-08)

Raw CSVs: `runs/scale_epyc7c13.csv` (NC 8k–216k), `runs/scale_large_2x16x.csv`
(NC ~439k–3.4M). Workers 1…32. Normalized work = `(10k)·steps/s = steps/s × NC/1e4`.

**For this design, ~16 workers is generally the fastest wall-clock config.**
Across sizes, the best steps/s is almost always at **8–16 workers**; **32 is
rarely better** (only once, NC≈439k on a quiet 7742) and often slightly slower
than 16. Small problems (NC ≲ 64k) barely gain past 1–4 workers.

**After ~0.5M cells, aggregate performance declines gently — no cliff.**
Best-case size-normalized throughput peaks near **NC ≈ 0.2–0.5M**
(~650–760 (10k)·steps/s, ~1.3–1.5 ms per 10k cells per step). From **~0.9M →
3.4M**, that rate falls only modestly (~630 → ~550 (10k)·s/s at 16W). Absolute
steps/s drops ~∝ 1/NC as expected; Mcell·s/s at fixed W stays roughly flat
(weak-scaling-ish), with a slow roll-off rather than a break.

| Regime | Typical best W | Notes |
|--------|----------------|-------|
| NC ≲ 64k | 1–8 | Multi-W often ≤1W; boundary tax dominates |
| NC ~0.1–0.5M | **16** (sometimes 32) | Best normalized rate |
| NC ≳ 0.9M | **16** | ~3.3× vs 1W; 32W ≤ 16W |

**Wall-clock 1 worker → 16 workers (EPYC 7742 large suite)** — not stuck at ~2×:

| NC | 1W steps/s | 16W steps/s | **speedup** | efficiency |
|---:|-----------:|------------:|------------:|-----------:|
| ~439k | 3.36 | 15.4 | **4.57×** | 29% |
| ~857k | 1.81 | 7.40 | **4.08×** | 26% |
| ~1.73M | 0.99 | 3.32 | **3.35×** | 21% |
| ~3.44M | 0.49 | 1.61 | **3.31×** | 21% |

On large NC the gain is **~3.3×** (mid band ~4–4.6×). Soft ceiling: 8W ≈ 16W ≈ 32W
at multi‑M cells; still well above 2×, far from ideal 16×.

2-worker efficiency is excellent on large NC (~90%+); it is poor only on tiny
fields. Per-core work at 16W is ~20–25% of solo 1W — barrier/ghost bound, not
DRAM-size cliff.

### Local baseline (Ryzen 5700G, GOMAXPROCS=16) — older v5

| NC | workers | speedup | with `-B` |
|---:|--------:|--------:|----------:|
| 64k | 8 | **~2.8×** | **~2.7×** |
| 64k | 16 | **~3.1×** | **~2.8×** |
| 216k | 8 | **~1.85×** | **~1.8×** |
| 216k | 16 | **~2.0×** | **~1.9×** |

Profile throughput @ 216k/w8: **~33 steps/s** (was ~15 earlier with real contacts).

`perf` blocked (`paranoid=4`); use `go run ./cmd/profile -n 216000 -workers 8`.

### What moved the needle (v5)

| Lever | Idea |
|-------|------|
| **Face `Publish` set** | Ghost sync only cut-plane cells, not all owned |
| **Publish-at-end / pull-at-start** | Drop mid-step full `syncGhostLo` |
| **Cr face exchange** | Mid-step only radii for boundary forces |
| **Fused phases** | ~4 pool waves/step instead of ~10 |
| **`Live` / `Active` lists** | Space wants/scale/apply + hops on contacts only |
| **No second area scan** | Keep pre-space A for hops; forces use Cr/d |
| **√ reject** | Skip lens math when `d² ≥ (ri+rj)²` |

Remaining #1 hotspot: **`refreshGeom`** (~25% CPU) — one full slot pass/step.
Sparse freecell (blob + bath) is where `Live` will cut that further.

### Allocations / `[][]` (init vs step)

**`Cluster.Step` / `World.Step`:** intended **zero heap** after init
(`par/alloc_test.go`). Live/Active lists reset with `[:0]` + bit clear; no
`make` on the hot path if FromWorld pre-sized.

**Fixed (was catastrophic at large NC):**
- `core.EdgeColor` — was `[][]bool` with **one `make` per cell**; now flat
  `[]bool` of length `nc*stride`.

**Still slice-of-slice / many small allocs (build / FromWorld only, not Step):**

| Site | Shape | When | Severity |
|------|--------|------|----------|
| `par.buildInteriorColors` | **`[][]bool`** `used[nOwn][8+]` + grow via `append` per row | every `FromWorld` | **Same class as fixed EdgeColor** — NOwn small allocs per shard |
| `par.buildInteriorColors` | **`buckets [][]int32`** (color → interior slots) | FromWorld | Fine if nColor tiny; still jagged |
| `par.FromWorld` | `[]map[int32]struct{}` ghostNeed (one map/worker) | FromWorld | OK for W≤32 |
| `par.FromWorld` | `[][]gslot` per worker, `[][]int32` needPub | FromWorld | Build-only append growth |
| `par.FromWorld` | `map[int32]int32` **g2l per shard** | FromWorld | Large: ~NOwn map entries |
| `par.FromWorld` | `map[[2]int32]struct{}` slot dedup per shard | FromWorld | O(slots on shard) |
| `core.BuildLattice` | `map[pair]struct{}` edge dedup | BuildLattice | Unnecessary on lattice (edges unique); O(E) map |
| `core` / `par` CSR | flat `[]int32` Cls/Inc + temporary `fill` | rebuildCSR | One-shot, flat — fine |

**Maps (`g2l`, dedup)** dominate FromWorld memory structure more than jagged
bool tables once EdgeColor is flat; **`buildInteriorColors` is the remaining
twin of the old EdgeColor bug** (per-owned-cell `make([]bool, 8)`).

### Structure notes

- z-slab ownership, exclusive per-worker `[]CellLo` / `[]CellHi`
- `Publish` / `BndTouch` / `Live` / `Active` auxiliary index lists
- R0=0.58 so unit lattice has real contacts (R0=0.5 ⇒ A=0 vacuum)

## Memory model (target)

```
Proc P:
  Cells []Cell     // [0:NOwn) owned, [NOwn:NOwn+NGhost) published halo
  Slots []Slot     // links with at least one owned endpoint
  Inc   []int32    // CSR of local slot indices into owned cells
```

No `*Cell` fields. Temporary `*Cell` from `&Cells[i]` is fine for local
mutation; the storage lives in the slice backing array.

## simd

Go 1.27+ has an experimental `simd` package in this environment. Not wired
yet — AoS layout and pass isolation come first; SIMD is a later pass on
contiguous field loops (or a SoA sidecar for hot meters).

---

## FCS / snapshot integration (future)

**Status:** design only — not implemented in `pfc/`. Reference writers:
`v93/fab/fcs.go`, `v93/kernel/freecell.c` (FCS v3, see `v90/FCS.md`).
Readers: `v93/cmd/volview`.

### What the production path does today

| Piece | Go (`fab/fcs.go`) | C (`freecell.c`) |
|-------|-------------------|------------------|
| User buffer | `bufio.NewWriterSize(f, **1<<20**)` (1 MiB) | stdio + `fflush` per chunk |
| CELL/LINK body | **`make` every frame** | **`malloc`/`free` every frame** |
| Codec 1 | new slices for transpose + shuffle + `zstd.EncodeAll` | **static** `fcs_s1/s2/s3` grown once |
| Flush | **`Flush()` every chunk** | `fflush` every chunk |
| Inner loops | linear pack over cells/slots; nested transpose | same |

Largest wins are **not** more clever format tweaks — they are: **pre-size
frame buffers**, **stop allocating on the snap path**, **pack in parallel
from shards**, **overlap compress/write with the next physics steps**.

---

### Pattern A — double-buffered FCS frames (primary sketch)

Allocate **two complete frame workspaces at init**. Each snap: **fill the
back buffer**, hand it to I/O, **swap**, physics continues into the other.

```
  INIT (once, sized from NC, max live links, max compress bound)
  ─────────────────────────────────────────────────────────────
  type FrameBuf struct {
      cellRaw  []byte   // 20 + NC * 12 * 4   (CELL payload)
      linkRaw  []byte   // 12 + NL_max * (8 + 4*4)
      nLink    int
      // codec scratch (optional; can live outside the pair if serial compress)
      s1, s2   []byte   // transpose + shuffle, len = max(len(cell), len(link))
      zdst     []byte   // zstd dest, cap = CompressBound(max raw)
      t        float64
      ready    bool
  }
  var frames [2]FrameBuf
  var fill, drain int   // indices 0|1
  // open snap_file once; bufio 1 MiB; Flush once per frame (not per chunk)

  SNAP DUE (after physics barrier — all shards quiescent)
  ────────────────────────────────────────────────────────
  F := &frames[fill]

  1) PARALLEL PACK  (same N shard goroutines, new phase phPackFCS)
        each shard id writes OWNED cells into F.cellRaw at global offsets:
            base := 20 + G0*12*4
            for i := 0..NOwn-1:
                pack f32 row from Lo[i] / Hi[i]  →  F.cellRaw[base + i*48 :]
        each shard packs INTERIOR links (both ends owned) into a
            pre-reserved link region OR a per-shard staging slice of fixed
            capacity (allocated at init), then a serial concat into F.linkRaw
        BOUNDARY links: one owner only (e.g. min global endpoint), global i,j

  2) BARRIER — pack complete; F is immutable until drain finishes

  3) SWAP
        fill, drain = drain, fill
        // frames[drain] is now the completed frame

  4) ASYNC DRAIN  (I/O goroutine, or main if snap is rare)
        optionally: colT → shuffle → zstd into frames[drain].zdst
        write FCS chunks (CELL, LINK, optional ANLZ)
        Flush() once
        mark frames[drain].ready = false

  5) PHYSICS continues using frames[fill] on the *next* snap
        — must not overwrite frames[drain] until I/O done
        — if I/O still busy when next snap is due: wait, drop snap, or
          triple-buffer
```

**Why this is the largest win**

| Before (current fab) | After (double frame) |
|----------------------|----------------------|
| alloc pack + transform + compress every snap | **0 heap** on snap after init |
| physics blocked on compress + disk | physics packs, then **overlaps** I/O |
| Flush per chunk | Flush **per frame** (or rarer) |
| single global SoA walk | **parallel pack** by shard into global index layout |

**Sizing (order of magnitude)**

- NC = 216k → CELL raw ≈ 10 MB  
- links ~ 3×NC undirected → LINK raw ≈ same order  
- two frames × (cell+link+scratch) ≈ tens of MB — fine next to the sim state  

**Invariants**

- Pack only **reads** physics state (after `Step`); does not disturb zero-alloc physics.  
- Global cell index for FCS = `G0 + local` (stable under fixed z-slab ownership).  
- Do not call into `fab.fcsCellFrame` as-is — it assumes one SoA `Sim` and allocates.

---

### Pattern B — other useful patterns (from the same insights)

**B1. Serial “C-style” scratch without double-buffer**  
One set of `cellRaw/linkRaw/s1/s2/zdst` reused every snap; physics **blocks**
until write returns. Still huge vs `make` every frame; simpler to ship first.
Double-buffer is B0 once snaps are frequent enough to matter.

**B2. Triple buffer**  
If compress+disk can exceed one `snap_every` interval: fill / compress / write
as a 3-stage pipeline so physics never waits on disk.

**B3. Face / active-only pack (sparse)**  
Reuse `Live` / `Active` / `Publish`: pack only non-quiet cells, or only a
ROI. Requires either dense zero-fill for volview compatibility or a schema
extension (column or ANLZ “index list”). Good when bath ≫ blob.

**B4. SoA pack buffer**  
Pack columns already transposed (12 planes of f32) so codec-1 skips
`fcsColT`’s nested loop — write straight into columnar layout. Shuffle+zstd
remain. Aligns with a future SoA physics sidecar.

**B5. Dedicated I/O thread + channel of frame indices**  
`chan int` of buffer index (0/1), capacity 1 — never allocate frames on the
channel. Physics sends `drain` index; I/O owns compress+write. Matches the
procedural worker style (no pool abstraction).

**B6. Memory-map growing `snap_file`**  
Optional: `mmap` append region for uncompressed debug streams; skip zstd
when interactive follow (`volview -follow`) wants low latency over size.

**B7. Keep ANLZ off the hot pack**  
Meters/tables are tiny; write from main after pack. Do not pull full-shard
scans for ANLZ unless needed — compute meters in a cheap parallel reduce
into pre-sized `float64` rows.

**B8. Flush policy**  
- Debug / follow: flush each frame.  
- Throughput campaign: buffer several frames in userspace (or rely on 1 MiB
  bufio) and flush every K snaps; accept truncated-tail risk on kill (FCS v3
  already parses up to last complete chunk).

**B9. Don’t put FCS inside physics phases**  
Never pack mid-`Step`. Extra traffic on `Lo`/`Hi` fights the cache plan.
Snap = **post-barrier phase** only.

**B10. Zero-alloc contract**  
Physics `Step`: already measured 0 allocs after init. Snap path: separate
budget — either also 0 (pre-sized + zstd into `zdst[:cap]`) or explicitly
“I/O may alloc” isolated from `par` tests.

---

### Integration guide (when we want this in `pfc`)

Do **not** bolt `v93/fab` into `Step`. Add a thin output package and a pack
phase.

#### Suggested layout

```
pfc/
  out/           # or fcs/ — format + double-buffer writer only
    fcs.go       # FCS1 v3 chunk headers, SCHM/CFG once
    frame.go     # FrameBuf, pair, swap, sizes from NC/NLmax
    pack.go      # pack cell row / link row (f32 LE) — no Sim type
  par/
    phases.go    # + phPackFCS
    workers.go   # unchanged barrier style
    step.go      # Step(); SnapIfDue() called by driver, not inside hot loop
  cmd/fabrun-like driver
    for step {
      cl.Step()
      if due { cl.PackFCS(frames[fill]); swap; kick I/O }
    }
```

`out` must not import physics beyond reading plain slices / a small
`PackView` struct. `par` may import `out` for `PackFCS`, or the driver
wires them (cleaner for testing).

#### Steps to implement (order)

1. **Port FCS v3 write primitives** without `fab.Sim` — take `[]byte` + counts.  
2. **Size and allocate `frames[2]`** at cluster build from `NC` and max links.  
3. **`phPackFCS`**: each shard fills its cell range; link staging as above.  
4. **Serial write path** with reused codec scratch (match C’s `fcs_s*`).  
5. **One Flush per frame**; open `snap_file` once for the run.  
6. **Double-buffer + async I/O** when profiling shows snap on the critical path.  
7. **Alloc test**: `Pack` + `Write` after warmup → 0 heap (or document exception for zstd if the encoder forces it).  
8. **volview smoke** on a short stream for byte compatibility with C/Go v93.

#### Wire-up sketch (API shape, not code in tree)

```text
// build
cl := par.FromWorld(w, workers)
fw := out.OpenStream(path, out.Opts{Compress: true, BufSize: 1<<20})
fw.WriteHeader(cfg, schema)
fb := out.NewDoubleFrame(cl.GlobalNC(), cl.MaxLinks())

// run loop
for step := 0; step < nsteps; step++ {
    cl.Step()
    if snapEvery > 0 && step % snapEvery == 0 {
        fb.WaitDrainIfNeeded()          // if previous async write still holds drain buf
        cl.PackInto(fb.Fill())          // barrier + parallel pack
        fb.Swap()
        fw.WriteFrameAsync(fb.Drain())  // or sync WriteFrame
    }
}
fw.Close()
cl.Stop()
```

`PackInto` maps to: `runPhase(phPackFCS)` then main finalizes LINK header
counts / boundary links if not fully parallel.

#### Mapping pfc state → FCS columns

| FCS CELL col | pfc source (owned cell `i`, global `g=G0+i`) |
|--------------|-----------------------------------------------|
| x,y,z,r | `Lo[i].X,Y,Z,Cr` |
| es,em,ee | `Lo[i].Es,Em,Ee` |
| xload | `(Em+flload)/cap` — flload may be 0 or derived |
| tag | apparatus (0 if unused) |
| fa1,fa2,th2 | `Lo[i].Fa1,Fa2,Th2` |

| FCS LINK col | pfc source |
|--------------|------------|
| i,j | **global** endpoint ids |
| d, A | `Slot.D`, `Slot.A` |
| lem | in-flight dense if present (0 in current pfc) |
| gg | gate product if clocks/gates exist (0 or compute from th2) |

Early pfc can write **subset columns** if SCHM lists only what exists —
volview maps by name. Prefer full 12/4 for drop-in volview of campaign runs.

#### Barriers and determinism

- Pack after a full `Step` barrier so all shards see a consistent time.  
- Byte-identical to v93 is **not** required for pfc experiments unless
  you need A/B against freecell streams; document float32 cast order.  
- Parallel pack into disjoint cell ranges is race-free; links need a
  clear ownership rule (interior vs boundary).

#### What to leave alone

- Do not add FCS to the zero-alloc physics `Step` test without separating
  packages — keep `TestStepZeroAlloc` about dynamics only.  
- Do not use roaring/bitmaps for FCS unless doing sparse/ROI snaps (B3).  
- Do not compress on every shard into one file without a single writer.

---

### Summary for future self

| Priority | Action |
|---------:|--------|
| 1 | Pre-allocate **two** frame buffers (+ codec scratch) at init |
| 2 | Parallel **pack** by shard → global index layout |
| 3 | **Swap**; single writer **flush once per frame** |
| 4 | Async drain so physics overlaps I/O |
| 5 | Optional sparse/ROI pack when activity masks are real |

That path reuses the multi-proc structure we already trust (fixed shard
goroutines, phase barriers, face publish) without reintroducing a job pool
or per-frame heap traffic on the hot path.
