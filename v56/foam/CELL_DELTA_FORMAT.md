# SFA cell-native delta frames (FCEP) — full resolution + sparse temporal deltas

Extension of the FMSH/FCEL format for dense temporal capture with
threshold-gated sparse delta encoding.

## Goal

Capture the simulation at fine `snap_dt` (e.g. 0.1, fully resolving the
T_breathe ≈ 2.2 oscillation cycle) without 30 GB SFA files. Achieved by
storing a per-cell temporal model (`mean + amp·cos(ω·t + φ)` per column,
matching the FRVD pattern in `scp_sim.c`) plus sparse FCEP frames that
list only cells whose actual value departs from the model prediction by
more than `delta_tol`.

For an L=40 box where ~95% of cells follow the carrier wave (which is
exactly `mean + amp·cos`), the per-frame residual is near-zero on most
cells and the FCEP payload is tiny.

## Frame types

```
FMSH   mesh frame (existing)            — Voronoi cells, written once or on remesh
FCEL   cell-data I-frame (existing/ext) — full state, optionally embeds temporal model
FCEP   cell-data P-frame (new)          — sparse delta from model prediction
```

Reader state machine:

- See `FMSH`           → load mesh, allocate state buffer
- See `FCEL`           → replace state buffer; if temporal-flag set, also load model
- See `FCEP`           → apply sparse deltas: `state[c] = predicted(t) + delta[c]` where
                          delta[c] = stored value if c is in the change-list, else 0

The state buffer is `N_cells × n_columns × dtype` floats; one of them
per file at any given time. After applying any frame the buffer holds
the full per-cell field at that frame's `t`.

## FCEL extended (backwards-compatible)

```
struct FCEL_payload {
    uint32_t magic;          // 'FCEL'
    uint32_t version;        // 1 (existing) or 2 (with model)
    uint32_t N_cells;
    uint32_t n_columns;
    uint8_t  dtype;
    uint8_t  flags;          // bit0=bss  bit1=delta_from_prev
                             // bit2=has_temporal_model
    uint8_t  reserved[2];

    // column-major cell values (existing)
    <dtype>[n_columns × N_cells]

    // present iff flags & bit2:
    float omega;
    float mean[n_columns × N_cells]
    float amp[n_columns × N_cells]
    float phase[n_columns × N_cells]
};
```

Old readers (version=1, bit2=0) parse exactly as before. New readers see
bit2 and consume the trailing model fields.

## FCEP payload (new)

```
struct FCEP_payload {
    uint32_t magic;          // 'FCEP'
    uint32_t version;        // 1
    uint32_t N_cells;        // total cell count of current mesh
    uint32_t n_columns;
    uint32_t n_changed;      // count of cells with |delta| > tol
    uint8_t  dtype;          // SFA_F32 expected
    uint8_t  flags;          // reserved=0 for now
    uint8_t  reserved[2];

    // change list:
    uint32_t cell_ids[n_changed];
    float    delta[n_changed × n_columns];   // residual: actual − predicted(t)
};
```

Wire size: `28 + n_changed × (4 + n_columns × 4)` bytes. For
n_columns=6 and 5% of 3.26M cells changing: ~28 + 163k × 28 = 4.5 MB
before zstd, ~1.5 MB compressed.

## Reader algorithm (volview, cellsfa_to_voxel)

Maintain:
- `state[N_cells × n_columns]` — current per-cell field values
- `mean, amp, phase[N_cells × n_columns]` — temporal model parameters
- `omega` — temporal-model angular frequency
- `model_valid` — bool

For each frame in time order:

| Frame type | Action |
|------------|--------|
| FMSH       | reset state/model, build voxel→cell index for rendering |
| FCEL v1    | copy values into `state`; mark model as invalid |
| FCEL v2 w/ temporal | copy values into `state`; load (omega, mean, amp, phase); mark valid |
| FCEP       | for each c in cell_ids: state[c, k] = mean[c, k] + amp[c, k]·cos(omega·t + phase[c, k]) + delta[c, k]<br>for c NOT in cell_ids: state[c, k] = mean[c, k] + amp[c, k]·cos(omega·t + phase[c, k]) |

After the frame is applied, render `state` via voxel→cell resampling.

## Writer algorithm (foam_sim)

Per-cell sliding accumulator (same pattern as scp_sim's FRVD writer):

```
for each frame at time t:
    accumulate sums:
        sum_cos[c, k]  += value[c, k] · cos(ω·t)
        sum_sin[c, k]  += value[c, k] · sin(ω·t)
        sum_mean[c, k] += value[c, k]
        count++

    is_iframe = (frame_index % cell_iframe_interval == 0)

    if is_iframe and count >= refit_min:
        # Refit model from accumulated sums (Fourier projection)
        mean[c, k]  = sum_mean[c, k] / count
        sc          = sum_cos[c, k] / count − mean × cos(ω·t)
        ss          = sum_sin[c, k] / count − mean × sin(ω·t)
        amp[c, k]   = 2 · sqrt(sc² + ss²)
        phase[c, k] = atan2(−ss, sc)
        zero the sums; count = 0

    if is_iframe:
        write FCEL with embedded temporal model
    else:
        # P-frame: compute deltas, threshold
        for c in 0..N_cells:
            for k in 0..n_columns:
                pred = mean[c, k] + amp[c, k] · cos(ω·t + phase[c, k])
                delta[c, k] = value[c, k] − pred
            if max_k |delta[c, k]| > cell_delta_tol:
                add c to change_list with deltas[k]
        write FCEP with change_list
```

## Bootstrap

First FCEL writes the actual values without a meaningful model
(mean = actual, amp=0, phase=0, ω=initial guess). Reader applies:
predicted = mean = actual, delta = 0, so first I-frame round-trips
exactly. The first refit happens after `refit_min` frames have
accumulated. Until then, P-frames degrade to "actual − mean" with
threshold — still useful for vacuum cells where mean is stable.

ω initial guess: `2π / T_breathe`, with `T_breathe` from the config or
inferred from the parameters. Default 2π/2.2 = 2.857 rad/time.

## Estimated compression (L=40, T=50, snap_dt=0.1)

| Format | Frames | Size |
|--------|--------|------|
| Voxel f16 BSS+zstd | 500 | ~30 GB |
| FCEL only (no deltas) | 500 | ~11 GB |
| FCEL + FCEP, no model | 50 + 450 | ~4 GB |
| **FCEL with model + FCEP, threshold 0.01** | **50 + 450** | **~3 GB** |

The model captures the carrier wave + breathing oscillator exactly for
~95% of cells. Only the soliton core and its near-field θ shell need
non-trivial delta encoding.

## API additions

```c
// Variant of sfa_write_cell_frame that embeds a temporal model.
int sfa_write_cell_iframe_temporal(SFA *s, double time,
                                    uint32_t N_cells, uint32_t n_columns,
                                    uint8_t dtype, void * const *column_data,
                                    float omega,
                                    const float *mean,
                                    const float *amp,
                                    const float *phase);

// Sparse cell P-frame.
int sfa_write_cell_pframe(SFA *s, double time,
                          uint32_t N_cells, uint32_t n_columns, uint8_t dtype,
                          uint32_t n_changed,
                          const uint32_t *cell_ids,
                          const float *delta_values);

// Reader: returns model-and-state expanded into the caller's buffer.
int sfa_read_cell_pframe(SFA *s, uint32_t frame_idx,
                         uint32_t *out_N_cells, uint32_t *out_n_columns,
                         uint32_t *out_n_changed,
                         uint32_t **out_cell_ids,
                         float **out_delta_values);
```
