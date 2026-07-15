# /space/scp disk plan — 2026-07-13

## Snapshot (start → mid-cleanup)

| | Start | After deletes (pre-archive complete) |
|--|------:|--------------------------------------:|
| /space used | 821G (95%) | **516G (60%)** |
| free | 49G | **354G** |
| Local SFAs | 357 / 583 GB | ~ fewer; archive in flight |
| B2 `scpsfa:` | 138 objs / 796 GB | +v74/v75/v70/v52/v66/v67 in progress |

## Actions executed

### DELETE (corrupt / 0-frame / bad download) — no archive

| Path | ~Size | Why |
|------|------:|-----|
| `v75/results/pm_orbit_D20_{vc,v05,v15}.sfa` | 48G | frames=0 |
| `v75/results/pm_headon_D20_vr01.sfa` | 3.2G | frames=0 partial |
| `v75/results/pp_force_D20.sfa` | 16G | frames=1 |
| `v52/proton_vecstream/proton.sfa` | 27G | frames=0 |

**Kept:** all tracks, diags, logs, seeds, renders. Step 2 science is in TSV/logs.

### DELETE (already on B2, size-matched)

~**157G**: `v52_expA_uu`, `v53/vacuum_full` (+dup), v55 foams, v72 qfi×4, v73 rings/fx,
v48/v50 previews, v44 seeds, v70 copies of v71 archive (mc2/im2/fc*).

### DELETE (superseded intermediates)

| Path | ~Size | Why |
|------|------:|-----|
| v53/grid_tests/* | 21G | grid scouts, not production physics |
| v52 proton_vecstream intermediates | ~44G | repaired/vec/vec_final/v2; keep only `proton_unified` for archive |
| v56 `*_OLD` | 3G | superseded |
| `results/v100_medium.sfa` | 2.5G | throwaway GPU test |

### ARCHIVE then delete (high value) — rclone batch

| Path | Size | Frames | Remote |
|------|-----:|-------:|--------|
| v74 `c6_light.sfa` | 36G | 122 | `scpsfa:scpsfa/v74/` |
| v74 `c12_light.sfa` | 40G | 122 | `scpsfa:scpsfa/v74/` |
| v75 `pm_force_D20.sfa` | 16G | 52 | `scpsfa:scpsfa/v75/` |
| v75 `pp_force_D16.sfa` | 16G | 52 | `scpsfa:scpsfa/v75/` |
| v75 `pm_force_D16.sfa` | 16G | 51 | `scpsfa:scpsfa/v75/` |
| v70 mc1 / fe1 / im1 | ~11G | ok | `scpsfa:scpsfa/v70/` |
| v52 `proton_unified.sfa` | 27G | 221 | `scpsfa:scpsfa/v52/` |
| v66/v67/*.sfa, v69_verify_r3 | ~15G | char. | matching prefixes |

**Note:** B2 SHA1 hashes full object before/during large upload — 36G ≈ few min
hash + transfer. Log: `/tmp/scp_rclone_logs/archive_20260713b.log`.

After remote size matches local: **delete local copies** of archived result SFAs.
Keep seeds + tracks + diags + pngs forever on disk.

## Still optional later (~50–80G)

- v50 em_wave/polariton (historical EM sector)
- v56 braid/skyrme L40 (real-field era)
- v51 collision_zoom (partial; r2 already on B2)
- v67/v66 after archive verify

## Keep local always

- `**/_diag.tsv`, `*_track.tsv`, `*.log`, `*.png`, `*.webp`
- Seeds under `/space/scp/v7{4,5}/seeds/` (regenerable but cheap)
- Templates referenced by configs

## Recovered headroom target

| Step | Free Δ |
|------|-------:|
| Done (corrupt + already-arch + intermediates) | ~**+300G** |
| After archive+delete of v74/v75/v52/v70 batch | ~**+150G more** |
| Stretch historical purge | +50–80G |
