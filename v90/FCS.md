# FCS — the free-cell snapshot format

`.fcs` is the v90 snapshot format for free-cell state: **a point cloud of
cells whose geometry is state, plus the live channel graph**. Written by
both kernels (`kernel/freecell.c` `write_fcs()`, `fab/fcs.go`) when
`snap_every > 0` and `snap_dir` is set; read by `cmd/volview`. Pure
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

## Layout (little-endian)

**Header — 28 bytes:**

| offset | type | field |
|---|---|---|
| 0 | `[4]byte` | magic `"FCS1"` (file type; the version is the next field) |
| 4 | `uint32` | version (1 or 2) |
| 8 | `uint32` | ncells |
| 12 | `float64` | t (simulation time) |
| 20 | `float64` | L (periodic box edge) |

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

Readers accept both versions (`cmd/volview` does). Writers write v2.
Fields are only ever appended; the version field gates the record size.
