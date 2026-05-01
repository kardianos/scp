# Field plug-in interface

The foam kernel (`foam_sim.c`) is generic over the per-cell field type.
A *field plug-in header* provides the field-specific bits:

- field count `NCOMP`
- bulk-norm metric `g_metric[NCOMP]`
- SFA column metadata (names, semantic, component)
- initial-condition functions (one per `init` string)
- self-coupling step that combines the kernel-computed Laplacian with
  the field's potential / interaction terms

To swap field theories, replace the `#include "field_pga.h"` line near
the top of `foam_sim.c` with a different header and rebuild. Examples:

- `field_pga.h` — V56 stage A: 8-component multivector in ℝ(3,1,1)
  with a Higgs quartic potential.
- `field_cosserat.h` (TODO) — V55-equivalent: 6 fields (3 φ + 3 θ),
  V(P) saturation potential, η × curl coupling.
- `field_dirac.h` (TODO) — V56 stage B: Dirac equation on the
  relativistic-quaternion field (first-order; would also need a new
  time integrator).

## Contract — a field header must define

### Always (constant block, before `#define FIELD_IMPL`)

| Symbol | Type | Description |
|--------|------|-------------|
| `NCOMP` | `#define` | number of field components per cell |
| `FIELD_NAME` | `#define` (string) | label printed at startup |
| `g_metric` | `static const double[NCOMP]` | diagonal metric for the bulk norm |
| `field_names` | `static const char *const[NCOMP]` | SFA column names (≤11 chars) |
| `field_semantic` | `static const uint8_t[NCOMP]` | SFA semantic codes (`SFA_POSITION`/`SFA_ANGLE`/`SFA_VELOCITY`/`SFA_CUSTOM`) |
| `field_component` | `static const uint8_t[NCOMP]` | 0..3 component index in the semantic group |

### Implementation (when `FIELD_IMPL` is defined)

```c
static void field_init(Sim *s);
static void field_forces(Sim *s,
                         double *const lap[NCOMP],
                         double *acc[NCOMP]);
```

- `field_init(s)` populates `s->M[k][i]` and `s->M_vel[k][i]` based on
  `s->c.init` (one of "vacuum", "qball", "skyrme", or whatever this
  field supports). Standard scalar config knobs available: `s->c.A`,
  `s->c.R_seed`, `s->c.omega_seed`, `s->c.v2`, `s->c.mass2`,
  `s->c.lambda`. Add field-specific knobs by extending `Config` and
  the `config_load` parser in `foam_sim.c`.
- `field_forces(s, lap, acc)` is called every step after the kernel
  has computed `lap[k][i] = ∇² M[k]` (per-cell Laplacian). The plug-in
  combines this with its own potential / coupling and writes
  `acc[k][i]`. Optionally updates `s->bulk_norm[i]` for diagnostics.

## Sim struct fields the plug-in may use

```c
struct Sim {
    Mesh   *m;             /* cell mesh, has cell_x, cell_y, cell_z, cell_vol */
    Config  c;             /* config struct (parsed from .cfg file) */
    double *M[NCOMP];      /* current field values per cell */
    double *M_vel[NCOMP];  /* time derivatives */
    double *M_acc[NCOMP];  /* accelerations (output of field_forces) */
    double *grad_M[NCOMP][3]; /* spatial gradient components (kernel fills) */
    double *bulk_norm;     /* per-cell bulk norm (plug-in fills) */
};
```

## Standard config fields

The kernel parses these out of `.cfg` files and exposes them via
`s->c`:

| Key | Field | Default | Meaning |
|-----|-------|---------|---------|
| `m` (or `mass2`) | `mass2` | 0 | bare-mass squared term |
| `lambda` | `lambda` | 5 | quartic coupling |
| `v` (or `v2`) | `v2` | 0.25 | vacuum norm-squared (\|M\|² = v²) |
| `A` | `A` | 0.3 | seed amplitude |
| `R_seed` | `R_seed` | 4 | seed localisation radius |
| `omega_seed` | `omega_seed` | 1 | seed time-frequency |
| `T` | `T` | 50 | total simulation time |
| `dt_factor` | `dt_factor` | 0.025 | dt = dt_factor × min cell distance |
| `init` | `init` | "vacuum" | seed selector — field-specific |
| `cell_native`, `sfa_output`, `cell_iframe_interval`, `cell_model_interval`, `cell_delta_tol`, `cell_omega`, `voxel_N`, `output`, `diag_file`, `snap_dt`, `diag_dt` | — | — | I/O knobs (see `CELL_FORMAT.md`, `CELL_DELTA_FORMAT.md`) |

Field-specific config keys can be added by patching `Config` and
`config_load` in `foam_sim.c`.

## Why this split

- The mesh, finite-volume operators, time integration, SFA I/O, and
  diagnostics framework don't depend on the field type. They live in
  the kernel.
- Whether the field is 6-component Cosserat, 8-component multivector,
  16-component bivector + spinor, or even an integer-quantised field
  for stage 2 — only the plug-in needs to change.
- The interface is small enough that a new field header is ~100–200
  LOC. Adding a new field theory does not touch any other file.

## Adding a new field plug-in

1. Copy `field_pga.h` to `field_<yourname>.h`.
2. Update `NCOMP`, `FIELD_NAME`, `g_metric[]`, `field_names[]`,
   `field_semantic[]`, `field_component[]`.
3. Replace the `field_init_*` and `field_forces` bodies with your
   theory's seeds and force-update.
4. Add any new config keys: extend the `Config` struct in `foam_sim.c`
   and the `config_load` parser. (TODO: move these into the field
   header too.)
5. Edit the `#include` in `foam_sim.c` to point at your new header.
6. Rebuild.

## Future cleanup

The next step for full plug-in independence is to move the seed-knob
fields (`A`, `R_seed`, `omega_seed`, `lambda`, `v2`, `mass2`) into the
field header so plug-ins declare their own config schema. The kernel
would expose only generic I/O + integration knobs.
