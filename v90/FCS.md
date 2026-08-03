# FCS — the free-cell snapshot format

`.fcs` is the v90 snapshot format for free-cell state: **a point cloud of
cells whose geometry is state, plus the live channel graph**. Written by
both kernels (`kernel/freecell.c` `write_fcs()`, `fab/fcs.go`) when
`snap_every > 0` and `snap_dir` (per-frame files) or `snap_file` (one
stream) is set; read by `cmd/volview`. Pure
output — writers consume no RNG and cannot perturb determinism.

## Why not SFA

SFA (`sfa/format/sfa.h`) is the pre-v89 archive format, and the
difference is not cosmetic — it is the background assumption itself:

| | **SFA** | **FCS** |
|---|---|---|
| what a file describes | a **fixed N³ lattice** — a permanent index set whose *contents* evolve | a **set of cells** — position, size, and connectivity are themselves dynamical state |
| the index | grid coordinate (i,j,k), immortal | cell index within one run; the *geometry* carries the physics |
| columns | 12/24/30 fixed columns (φ, θ vectors, velocities, gauge blocks, ρ_Q) tied to the lattice field content | per-cell energies by **mode** (space / dense / field), the load, the phases — laws_V2g's state |
| connectivity | implicit (grid neighbours, forever) | **explicit link records** — channels exist, are born, and die (rule α) |
| file shape | one large archive, many frames, zstd-compressed, seek table | one small file per snapshot (a run is a file sequence); no compression needed (~50 B/cell) |
| scale | GB (N³ float fields) | KB–MB (thousands of cells) |

v89 rejected the background — nothing may persist and be merely
re-valued — so it explicitly excludes SFA and its instruments. FCS is
what remains when the grid goes: the cell list *is* the state.

## Versions at a glance

* **v1/v2** — one snapshot per file, fixed record layout (legacy;
  readers keep them working). A run is a directory of files.
* **v3 (current)** — a **chunked stream**: one file carries a whole run —
  config provenance, a self-describing column schema, CELL/LINK field
  frames, and typed **ANLZ instrumentation frames interleaved freely**
  between (or instead of) field frames. Extending the format means adding
  a column name or a chunk type; readers skip what they do not know.
  Both kernels write v3 (`snap_dir` = per-frame v3 files; **`snap_file` =
  one appended stream**). Appended chunks are flushed per frame, so a
  truncated file parses up to the last complete chunk.

## v3 layout (little-endian)

**File header — 8 bytes:** magic `"FCS1"`, `uint32` version = 3.
Then chunks until EOF:

| type | field |
|---|---|
| `[4]byte` | chunk fourcc |
| `uint64` | payload length |
| … | payload |

**Chunk types:**

* `CFG ` — UTF-8 text: the configuration the kernel actually ran
  (laws + apparatus + seed + exp), written once at stream start. The
  format carries its own provenance — "print the configuration beside
  every result" applies to files too.
* `SCHM` — the column schema, once, before the first frame:
  `uint32 ncellcols`, then per column `uint8 len + name`; then the same
  for link float-columns. **Readers map columns by name** — a kernel that
  appends a new per-cell quantity changes nothing for old readers, and
  the viewer offers the new column as a channel automatically.
* `CELL` — a field frame: `float64 t`, `float64 L`, `uint32 ncells`,
  then ncells × ncellcols × `float32`.
* `LINK` — the channel graph at the same t (written right after its
  CELL): `float64 t`, `uint32 nlinks`, then per link `uint32 i`,
  `uint32 j`, nlinkcols × `float32`.
* `ANLZ` — an instrumentation frame, either form:
  * kind 0, **text**: `uint8 0`, UTF-8 payload (e.g. `# RESULT …` lines);
  * kind 1, **table**: `uint8 1`, `uint8 len + name`, `float64 t`,
    `uint32 ncols`, per column `uint8 len + name`, `uint32 nrows`,
    nrows × ncols × `float64`. This is the interleaved-analysis carrier:
    e.g. the double slit writes its accumulating screen profile
    `ds_screen(y, I)` at every snapshot time — the fringe build-up is in
    the stream — plus a per-frame `meters` table (drift, φ, z, births,
    deaths, …). An analysis-only stream with a few field frames woven in
    is equally legal: chunk order is unconstrained after SCHM.

* `CMPD` — **lossless compression wrapper** (default on; `snap_comp=0`
  writes raw): payload = inner fourcc, `uint64` raw length, `uint8` codec,
  compressed bytes. **Codec 1** = columnar transpose (4-byte units per
  record) + byte shuffle + zstd, chosen by measurement on real streams
  (double-slit run, 3.05 MB of CELL/LINK payloads):

  | codec | ratio | encode |
  |---|---|---|
  | zstd raw | 0.668 | 84 MB/s |
  | flate-6 | 0.674 | 52 MB/s |
  | **columnar+shuffle+zstd (codec 1)** | **0.547** (0.464 on quasi-static streams) | **174 MB/s** |
  | + XOR frame delta | 0.524 | 162 MB/s — rejected: only ~2-4% more for losing random access and per-chunk self-containment |

  Per-chunk self-contained: random access, truncation safety, and live
  following all survive compression. Verified lossless end-to-end
  (renders from compressed and raw streams are byte-identical). Note:
  the C kernel (libzstd) and Go kernel (klauspost) emit different valid
  zstd bytes — parity across kernels is of the *decoded* chunks.

**Current cell columns** (v3 schema names): `x y z r es em ee xload tag
fa1 fa2 th2`. **Link columns:** `d A lem gg`.

## Live streaming

Kernels flush whole chunks per frame, so a stream is followable while
the run is writing it: `volview -i -follow run.fcs` waits for the file,
attaches on the first frame, and appends frames live (pinned to newest;
scrubbing back unpins, `E` re-pins). The battery's heavier experiments
(blob, pulse, ds_m0) write streams to `runs/streams/` as standard output
(regenerable, gitignored).

## Legacy v1/v2 layout

**Header — 28 bytes:** magic `"FCS1"`, `uint32` version (1 or 2),
`uint32` ncells, `float64` t, `float64` L.

**Cell records — ncells × (9 floats v1 / 12 floats v2), float32 each:**

| # | field | meaning |
|---|---|---|
| 0–2 | x, y, z | position (dynamical — this is data, not an index) |
| 3 | r | live radius r = r0·(Es/e_s0)^{1/3} |
| 4 | Es | space-mode store |
| 5 | Em | dense-mode store |
| 6 | Ee | field-mode energy \|ψ\|² |
| 7 | xload | (Em + flload)/cap — the pitch-detune load |
| 8 | tag | seeded-object flag |
| 9 | fa1 | *(v2)* field amplitude, plane 1 (Re ψ) |
| 10 | fa2 | *(v2)* field amplitude, plane 2 (Im ψ) |
| 11 | th2 | *(v2)* dense clock phase |

v2 exists so a viewer can show **aspects of the field**, not just its
energy: phase arg ψ = atan2(fa2, fa1), the signed components, and the
dense phase — the two-plane amplitude is the repaired field sector
(DOUBLESLIT §8.1), and \|ψ\|² alone cannot show interference *mechanism*,
only its result.

**Link records — v2 only, after the cells:**

`uint32 nlinks`, then per link (every non-free slot):

| type | field | meaning |
|---|---|---|
| `uint32` | i, j | endpoint cell indices |
| `float32` | d | live separation |
| `float32` | A | lens (contact) area — A > 0 is a live channel |
| `float32` | lem | in-flight dense energy (both directions) |
| `float32` | gg | gate product g_f·g_b at the locked partial (bond visibility) |

The link block is what SFA structurally cannot express: channels with
identity. A dying channel (out of the candidate set, still flushing) has
A = 0 with lem > 0.

## Versioning

Readers accept v1, v2, and v3 (`cmd/volview` does; `volview -info`
prints a stream's chunk inventory and its ANLZ tables). Writers write
v3. Within v3, evolution is additive by construction: new cell/link
columns are announced in SCHM and matched by name; new chunk types are
skipped by length; nothing re-breaks the version.
