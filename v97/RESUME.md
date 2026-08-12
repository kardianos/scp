# v97 RESUME — read this first after compaction

Working dir: `/home/d/code/scp`. Active programme: **v97**.

## 0. What v97 is

**THE REGULAR METHOD.** Conversion is promotions of one energy on
oriented links, allowed only where a cycle can close. Properties are
moments of a live conversion dipole, not attributes of an Em store.
Correspondence is identity or it is nothing. The six principles and
the dipole hunt are adopted. The kernel is the v95 carry, byte-inert
at defaults. No principle-face is implemented yet.

## 1. Reading order

1. This file.
2. `v97/README.md` — charter (self-contained).
3. `v97/PRINCIPLES.md` — the six, in full.
4. `v97/METRICS.md` — the adopted hunt (H0–H5).
5. `CLAUDE.md` — project conventions (banner now points here).
6. Evidence as needed, newest first: `v96/RESUME.md` (no-consumption /
   no-pricing diagnosis), `v95/consult/DESIGN_claude.md` (J-metric),
   `v94/FINAL_STATUS.md` (door pair), `v93/RECKONING.md` +
   `v93/RING_DNLS.md` + `v93/theory/FINAL_STATUS_ADOPTION.md`,
   `v93/theory/qgraph/docs/UNIFIED_LAW.md` (promotion algebra,
   symbolic only).
7. Foundation: `v89/PRINCIPLE.md`, `v89/FREECELL.md`.

Do not read anything before v89.

## 2. Status

| item | status |
|---|---|
| Charter / principles / metrics | **adopted** (2026-08-12) |
| STEP ZERO carry v95 → v97 | **done** (this opening) |
| Battery 87 at defaults | **ALL GREEN 87** (`runs/BATTERY.log`) |
| H0 meter honesty | **PASS** (5/5 classes; current-direction SKIP) — `H0.md` |
| H1 dipole hunt | **EMPTY** (0/203 inhabited) — `H1.md` |
| Promo leads (Go `promo=1/2/3`) | **tried** — sep+Es now co-occur; H1 still EMPTY — `LEADS.md` |
| Kernel faces (slot door, promotion write, …) | leads only; **not adopted** |
| Law-table change | **none** |

## 3. The code (carried, all default-0)

C kernel of record: `v97/kernel/freecell.c`.
Go A/B mirror: `v97/fab/`.
Build: `cd v97 && make all` (includes `dipole`).
Hunt: `./dipole -h0 ring ../v93/runs/quench/rv6_uni_m2.fcs`.

Carried knobs that exist and must stay inert at 0:

| knob | from | what |
|---|---|---|
| `amp_nat` / `amp_logate` / `hop_order` | v93 | unitary dense hop, linearize, Strang |
| `amp_door` | v93 | arg(ψ) door + reverse evap compose |
| `seed_mw` | v93 | painted matter winding |
| `door_g` / `door_avg` / `chi_load` / `door_out_enable` | v94 | door pair + register `χ` |
| `j_metric` | v95 | FLOW bed on dense J |
| `rad_dual` / `bath_xb` | v93 theory | dual intake (not adopted) |
| `promo` | v97 leads | oriented door (0=off; not adopted) |

Do not flip any of these as “the v97 law.” They are instruments.

## 4. Next work (when the user wants to continue)

Instruments are built (Go-only, C ignored): `hunt/` + `./dipole`.
H0 is green (`H0.md`). H1 is **EMPTY** (`H1.md`).

1. H1 is still EMPTY after the promo leads (`LEADS.md`). The
   missing pair (sep + Es cycle) now exists; Φ/Em and J-lock do not.
2. H2–H5 still have no motif. Do not adopt `promo`.

## 5. Conventions (unchanged)

- Ratchet: full battery before any kernel/law commit; bars sharpen
  by measurement, never soften; 87 GREEN at defaults is the
  invariant surface.
- Byte-inertness of new knobs at 0.
- Both kernels (C record, Go mirror). `math.Sincos` returns
  `(sin, cos)`.
- Pre-v89 ban.
- No kernel change without explicit user sign-off. Opening v97 is
  **not** that sign-off for P3/P1/P5 faces.
- Detached long runs: `setsid bash -c '…; touch out.done' </dev/null >/dev/null 2>&1 & disown`.

## 6. Gotchas carried forward (do not re-learn)

- `exp=bath` is matterless — never grade conservation or ignition on it.
- Tagged-centroid `blob_drift` lies about the unitary channel; use `fd`.
- Empty-cell `th2 → 0` was a real bug (fixed in v93); do not reintroduce.
- `p_gate=8` random-phase mean is ~0.1–0.2, not ≈0.
- Face C t=300 was a transient peak (R1/R2/R3).
- v94 stimulated outtake (`door_out_enable=1`) is the named pathology
  across M4/M5; do not “just turn it on” for the hunt.
- v95 J-metric amplifies a current; it does not generate one.
