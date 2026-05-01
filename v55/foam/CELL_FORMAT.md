# SFA cell-native frame format (FMSH + FCEL)

Extension of `sfa/format/sfa.h` that lets foam_sim write per-cell field
data directly into an SFA file, without intermediate voxel resampling.
Backwards compatible — existing voxel SFAs still load.

## Concept

The mesh is stored as a frame, not a header chunk. The reader iterates
frames in time order and maintains "current mesh" state:

- See `FMSH` frame → parse and replace current mesh
- See `FCEL` frame → decode using current mesh
- See `FRMD` / `FRVD` frame → as before (voxel / vec)

This allows multiple resolutions in one file:

```
t=0   FMSH        ← high-res mesh (e.g. 3M cells)
t=0   FCEL        ← cell data on high-res
t=1   FCEL
t=2   FCEL
...
t=10  FMSH        ← coarsened mesh (e.g. 500k cells)
t=10  FCEL        ← cell data on coarse mesh
t=11  FCEL
...
```

A reader without FMSH/FCEL support skips them via the JMPF index
(frame_type tells it to skip), so old tools degrade gracefully.

## New chunk types

```c
#define SFA_CHUNK_FMSH  0x48534D46  /* "FMSH" — Voronoi mesh   */
#define SFA_CHUNK_FCEL  0x4C454346  /* "FCEL" — cell-data frame */

#define SFA_FRAME_MESH  4   /* mesh frame (FMSH chunk) */
#define SFA_FRAME_CELL  5   /* cell-data frame (FCEL chunk) */
```

## FMSH payload (zstd-compressed)

```
struct FMSH_payload {
    uint32_t magic;          // 'FMSH'
    uint32_t version;        // 1
    uint32_t N_cells;
    uint32_t N_faces;        // 0 if positions-only
    uint8_t  flags;          // bit0:has_faces  bit1:has_csr  bit2:pos_f64
    uint8_t  reserved[3];
    double   L;              // box half-extent (cubic box)

    // cell positions: 3*N_cells values, f32 (default) or f64
    float|double  pos[3 * N_cells];
    float|double  volume[N_cells];

    // optional face records (40 bytes each) iff flags & 0x1
    struct { uint32_t a, b; double area; double dx, dy, dz; }
        faces[N_faces];

    // optional CSR (cell→face) iff flags & 0x2
    uint32_t cell_face_off[N_cells + 1];
    uint32_t cell_face_idx[total];   // total = cell_face_off[N_cells]
};
```

The face data and CSR are optional. For visualization-only SFAs we can
skip them and store ~10× less mesh data. For restart-capable SFAs we
need them.

## FCEL payload (zstd-compressed)

```
struct FCEL_payload {
    uint32_t magic;          // 'FCEL'
    uint32_t version;        // 1
    uint32_t N_cells;        // must match current mesh
    uint32_t n_columns;
    uint8_t  dtype;          // SFA_F16/F32/F64
    uint8_t  flags;          // bit0:bss_encoded  bit1:delta-from-prev
    uint8_t  reserved[2];
    // column-major: column 0 (N_cells × dtype), column 1, ...
    // or BSS-interleaved if flags & 0x1
};
```

## API additions to sfa.h

```c
// Write a mesh frame. cell_pos is interleaved x,y,z.
int sfa_write_mesh_frame(SFA *s, double time,
                         uint32_t N_cells, uint32_t N_faces,
                         double L,
                         const double *cell_pos,    // 3*N_cells doubles
                         const double *cell_vol,    // N_cells doubles
                         const void *face_records,  // optional
                         const uint32_t *csr_off,
                         const uint32_t *csr_idx,
                         uint8_t flags);

// Write a cell-data frame. column_data[c] is N_cells × dtype_size bytes.
int sfa_write_cell_frame(SFA *s, double time,
                         uint32_t N_cells, uint32_t n_columns,
                         uint8_t dtype, void **column_data);

// Read: extends sfa_read_frame to handle FCEL when current mesh is loaded.
// New helper to read the mesh frame at index frame_idx into a struct:
typedef struct {
    uint32_t N_cells, N_faces;
    double L;
    double *cell_pos;       // 3*N_cells, malloc'd
    double *cell_vol;       // N_cells
    void   *face_records;   // optional
    uint32_t *csr_off, *csr_idx;
    uint8_t flags;
} SFA_Mesh;

int sfa_read_mesh_frame(SFA *s, uint32_t frame_idx, SFA_Mesh *out);
void sfa_mesh_free(SFA_Mesh *m);
```

The existing `sfa_read_frame` is voxel-only and continues to work for
FRMD frames. To read FCEL data, callers use a new function:

```c
int sfa_read_cell_frame(SFA *s, uint32_t frame_idx,
                        uint32_t *out_N_cells, uint32_t *out_n_cols,
                        uint8_t *out_dtype, void **out_columns);
```

## volview integration

volview reads SFA in its existing pipeline, which builds a 3D voxel
texture for the ray-marcher. New behavior:

1. Iterate JMPF entries in time order.
2. On `SFA_FRAME_MESH`: parse and store current mesh + spatial bin index.
3. On `SFA_FRAME_CELL`: load cell values, resample to a voxel buffer
   using the current mesh's spatial index, upload as a voxel texture.

The voxel resolution is a volview parameter (CLI flag `-voxel-N` default
192 for L=40 boxes). This effectively moves the foam→voxel resample from
simulation time to view time, where it amortizes across cached frames.

## Migration plan

| Step | Files touched | Outcome |
|------|---------------|---------|
| 1. sfa.h additions | `sfa/format/sfa.h` | Reader+writer for FMSH/FCEL |
| 2. foam_sim cell output | `v55/foam/foam_sim.c` | Writes FMSH+FCEL instead of FRMD |
| 3. foam_to_voxel reader | `v55/foam/foam_to_voxel.c` | Reads cell SFA, writes voxel SFA |
| 4. volview support | `sfa/volview/main.go` | Renders cell SFA directly |
| 5. Validation | re-run L=40, compare files | Verify size + visual match |

After step 1-2, foam_sim writes a 50 MB cell SFA in place of the 3.4 GB
voxel SFA. After step 4, volview displays it directly without an
intermediate conversion step.
