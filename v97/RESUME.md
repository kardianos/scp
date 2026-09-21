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
| Battery 87 at defaults | **ALL GREEN 87** (`runs/regrade/BATTERY_helixseed.log`; the leads commit had silently broken 6 C≡Go parity bars — found, proven pre-existing, fixed — `REGRADE.md` §4b) |
| H0 meter honesty | **PASS** (now scripted `runs/h0.sh`; + `exact` class) — `H0.md`, `runs/regrade/H0_recheck.out` |
| H1 dipole hunt | **EMPTY** — archive 0/259 on honest meters (`runs/regrade/`), leads 0/10 (`REGRADE.md`) |
| PLAN (2026-08-13, grok-consulted) | **executed** — `PLAN.md`; consult `consult/REVIEW_grok_next.md` |
| Exact-ledger meters | **adopted in hunt/** — heuristic was inflating I ~14–39×; pI/pO now grade whenever present; `-dilate` book-scope; regression tests pin kernel RESULT totals |
| A5 ring campaign (7 arms + T=1000) | **done** — `A5_RING.md` + dated correction: the apparent Φ pass was dilate-halo bookkeeping (**retracted** — see §5b row and `REPLICATE.md`); what survives: orientation splits I/O wherever the door fires; J and park failures localized |
| Theory | `theory/PHI_BAR.md` (Φ bar passable at g=0; g-term is H2), `theory/DOOR_RULE.md` (selection rule = DFT of door profile; oriented edge door is ANTI-selective — pre-registered Δm=±1 prediction REFUTED), `theory/calc/` tests green |
| Kernel edits (Go) | **two, user-authorized 2026-08-13**: composed helix seed (kx now composes with seed_mw) + promo-header parity gate; battery ALL GREEN 87 after |
| §5b queue (grok state-check, user-directed L=24 dt=0.01) | **executed** — `REPLICATE.md`: Φ never passed on honest books (dest kill band + 0/6 seeds); drive window starve↔flood, no balance; promo=3 closed; `pJ` meter H0 GREEN (snapshot J can invert sign — J clause now grades integrated channel); rf_pj = first trustworthy J verdict (net directed current, direction unpinned) |
| Kernel faces (steer/door-locus signed credit; g-term/H2 after a 0.2–0.3 drive sweep; source-booked ledgers) | named with honest exhibits in `REPLICATE.md`; **not opened, need sign-off** |
| `CHARTEXP.md` (other agent, fork) | noted; 3+1 chart projection = candidate meter-only instrument; chart write-paths refuted — do not adopt |
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

Read `PLAN.md` → `REGRADE.md` → `A5_RING.md` → `REPLICATE.md` →
`theory/` first.
H1 is EMPTY on honest meters at both scales (`REPLICATE.md` campaign
verdict). The measured blockers: Φ has never passed honest books
(asynchrony + one-sided intake + seed-contingent feeding + a
starve↔flood drive window with no observed balance); the J direction
is genuinely unpinned (rf_pj, integrated channel — the first
trustworthy J verdict). Candidate next steps, each needing a brief +
user sign-off where they touch the kernel:

1. **Door-locus signed credit** (steer face): stop `promoPick`
   re-picking from the snapshot compass; per-link signed promotion
   integral with FLOW-τ decay. Now backed by an honest exhibit
   (rf_pj) — brief it as "the door listens to a bad clock."
2. **Drive sweep first, then the g-term (H2)**: a 0.2–0.3
   `noise_amp` sweep decides whether any g=0 balance exists between
   the measured starve (0.1) and flood (0.4) ends. Only then does
   the stimulated-intake brief have a complete exhibit. v92
   feedback-closure priors must be addressed in it.
3. **Instrument debts**: source-booked ledgers (dest replay is a
   lower bound); unconditional ledger columns (promo=0 controls
   can't be exact-graded); optionally the CHARTEXP 3+1 projection
   as a meter-only hunt/ANLZ column.

Do not adopt `promo`. Do not move 0.015/0.80. H2–H5 still wait on
an inhabited H1.

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
- (2026-08-13) ring seeds: `kx` now composes ON TOP of `seed_mw`
  (helix); previously kx+seed_mw silently discarded the tilt.
- (2026-08-13) the Go `# v97 promo leads` header prints only at
  promo≠0 (C≡Go parity at defaults depends on it).
- promo>0 books door events at the link's FAR cell — grade motifs
  with `dipole -dest` (replayed picks; `-dilate` is diagnostic only —
  its halo faked the retracted Φ pass).
- The mode-crossing heuristic INFLATES I/O (~14–39×); trust it only
  on streams with no `promo_*` columns, and prefer small tagged
  motifs over foam censuses (census dilution, PHI_BAR §C).
