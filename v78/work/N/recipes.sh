#!/usr/bin/env bash
# Phase N — nucleon profile + seed recipes (existing tools/profiles only).
# Usage:
#   ./recipes.sh              # print recipes (no write if OUTDIR missing tools)
#   ./recipes.sh run          # generate single-nucleon seeds under SEED_OUT
#   ./recipes.sh profile-note # print how profiles were / would be (re)shot
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/../../.." && pwd)"
V74="$ROOT/v74"
V69="$ROOT/v69"
BIN="$ROOT/bin"
PROF_LIGHT="$V74/profiles/f_w146_g005.txt"
PROF_STD="$V74/profiles/f_w142_g005.txt"
PROF_EDGE="$V74/profiles/f_w1485_g005.txt"
GSCAN="$V69/theory/gscan.tsv"

# Default single-nucleon grid (smoke / local). Production nuclear often N=192 L=36.
N="${N:-128}"
L="${L:-24}"
SEED_OUT="${SEED_OUT:-/space/scp/v78/seeds}"
MODE="${1:-help}"

echo "=== v78 Phase N nucleon recipes ==="
echo "ROOT=$ROOT"
echo "profiles: light=$PROF_LIGHT"
echo "          std  =$PROF_STD"
echo "gscan   : $GSCAN"
echo "N=$N L=$L SEED_OUT=$SEED_OUT"
echo

# ---------------------------------------------------------------------------
# A) Use existing gauged profiles (RECOMMENDED — already match gscan g=0.05)
# ---------------------------------------------------------------------------
# Files are 2-col "r f" (per-component amplitude). Do not re-normalize by √3.
#
# light  ω=1.46  Q≈114.13  E/Q≈1.5187  r½≈2.63
# std    ω=1.42  Q≈311.46  E/Q≈1.4649  r½≈3.80
# edge   ω=1.485 (VK sign flip at g=0.05 — not free-nucleon default)

# ---------------------------------------------------------------------------
# B) Optional: re-shoot gauged branch (only if regenerating profiles)
# ---------------------------------------------------------------------------
# cd "$V69/theory" && python3 gauged_shooter_fast.py
# → writes gscan.tsv and multi-col gprofile_*.txt; extract "r f" columns to
#   match v74/profiles/ layout for gen_qball_* loaders.

# ---------------------------------------------------------------------------
# C) Ungauged radial_qball (NOT for g=0.05 particle seeds — no Coulomb)
# ---------------------------------------------------------------------------
# "$BIN/radial_qball" -omega 1.46 -profile /tmp/f_w146_ungauged.txt
# Window: ω ∈ (1.3087, 1.5). Prefer gauged profiles above for all campaign work.

# ---------------------------------------------------------------------------
# D) Single static nucleon seed — gen_qball_boost vx=0
# ---------------------------------------------------------------------------
# Usage: gen_qball_boost N L profile omega vx out.sfa

seed_light() {
  local out="$SEED_OUT/nucleon_light_w146_N${N}.sfa"
  mkdir -p "$SEED_OUT"
  "$BIN/gen_qball_boost" "$N" "$L" "$PROF_LIGHT" 1.46 0.0 "$out"
  echo "wrote $out"
}

seed_std() {
  local out="$SEED_OUT/nucleon_std_w142_N${N}.sfa"
  mkdir -p "$SEED_OUT"
  "$BIN/gen_qball_boost" "$N" "$L" "$PROF_STD" 1.42 0.0 "$out"
  echo "wrote $out"
}

# Equivalent one-ball multi (same physics, origin-centered):
# "$BIN/gen_qball_multi" "$N" "$L" out.sfa "$PROF_LIGHT" 1.46 0.0  0 0 0

# ---------------------------------------------------------------------------
# E) Boosted single nucleon (optional kinematics)
# ---------------------------------------------------------------------------
# "$BIN/gen_qball_boost" "$N" "$L" "$PROF_LIGHT" 1.46 0.1 \
#   "$SEED_OUT/nucleon_light_vx0.1_N${N}.sfa"

# ---------------------------------------------------------------------------
# F) Multi-nucleon fusion seeds (Phase C — reference only)
# ---------------------------------------------------------------------------
# Light co-phase balls for carbon inventory:
#   "$V74/analysis/make_carbon_seeds.sh" 192 36
# Manual single pair (D=10 along x):
#   "$BIN/gen_qball_multi" 192 36 out.sfa \
#     "$PROF_LIGHT" 1.46 0.0   5 0 0 \
#     "$PROF_LIGHT" 1.46 0.0  -5 0 0
# Opposite charge (Phase L prep): second omega negative
#   "$BIN/gen_qball_multi" 192 48 out.sfa \
#     "$PROF_LIGHT"  1.46 0.0   8 0 0 \
#     "$PROF_LIGHT" -1.46 0.0  -8 0 0

# ---------------------------------------------------------------------------
# G) Sim config skeleton (init=sfa) — N,L MUST match seed
# ---------------------------------------------------------------------------
# N = $N
# L = $L
# T = 100
# dt_factor = 0.025
# m = 1.5
# m_theta = 1.6
# eta = 0
# mu = -41.345
# kappa = 50
# complex_phi = 1
# complex_gauge = 1
# g_gauge = 0.05
# bc_type = 0
# damp_width = 3.0
# damp_rate = 0.01
# init = sfa
# init_sfa = <seed path>
# init_frame = 0
# output = ...
# diag_file = ...
# precision = 0
# snap_dt = 2.5
# diag_dt = 0.25
#
# Run via scp-runner (preferred) or local:
#   bin/scp_sim path/to.cfg
# Do NOT edit scp_sim for this recipe.

# ---------------------------------------------------------------------------
# H) Stability / branch checks (analysis, not seed)
# ---------------------------------------------------------------------------
# - Retain Q over long T; expect negligible dQ/dt for isolated free ball
# - On-branch: E/Q vs gscan at measured Q (awk/python on gscan.tsv)
# - gauss_max ~ 1e-13 floor entire run
# Example gscan lookup (light):
#   awk -F'\t' '$1==0.05 && $2==1.46 {print "Q="$6,"Etot="$9,"rhalf="$12,"dQsign="$11}' "$GSCAN"

case "$MODE" in
  help|print|"")
    cat <<'EOF'
Commands documented above (A–H). Existing profiles are the campaign freeze.

  ./recipes.sh run           # write light + standard static seeds
  ./recipes.sh run-light     # light only
  ./recipes.sh run-std       # standard only
  ./recipes.sh profile-note  # profile provenance
  ./recipes.sh gscan-light   # print gscan row ω=1.46 g=0.05

Env: N L SEED_OUT (defaults 128 24 /space/scp/v78/seeds)
EOF
    # Still print the gscan fingerprint when available
    if [[ -f "$GSCAN" ]]; then
      echo
      echo "--- gscan fingerprint g=0.05 (from file) ---"
      awk -F'\t' 'NR==1 || ($1+0==0.05 && ($2+0==1.42 || $2+0==1.46 || $2+0==1.4825 || $2+0==1.485)) {
        if(NR==1) print; else printf "ω=%s Q=%s Etot=%s E/Q=%.4f r½=%s dQsign=%s\n",$2,$6,$9,$9/$6,$12,$11
      }' "$GSCAN"
    fi
    ;;
  run)
    [[ -x "$BIN/gen_qball_boost" ]] || { echo "FATAL: $BIN/gen_qball_boost missing — make -C sfa install" >&2; exit 1; }
    [[ -f "$PROF_LIGHT" && -f "$PROF_STD" ]] || { echo "FATAL: profiles missing under $V74/profiles" >&2; exit 1; }
    seed_light
    seed_std
    ls -lh "$SEED_OUT"/nucleon_*_N${N}.sfa
    ;;
  run-light)
    [[ -x "$BIN/gen_qball_boost" ]] || { echo "FATAL: gen_qball_boost missing" >&2; exit 1; }
    seed_light
    ;;
  run-std)
    [[ -x "$BIN/gen_qball_boost" ]] || { echo "FATAL: gen_qball_boost missing" >&2; exit 1; }
    seed_std
    ;;
  profile-note)
    cat <<EOF
Provenance:
  v74/profiles/f_w*_g005.txt  — 2-col r,f extracts of gauged radial solutions
  at g=0.05 matching v69/theory/gscan.tsv (f0,Q,r_half).

  Full multi-col shooter output (r f Er weff) example:
    $V69/theory/gprofile_w142_g005.txt

  Re-shoot: python3 $V69/theory/gauged_shooter_fast.py
  Ungauged (wrong for g=0.05 seeds): $BIN/radial_qball -omega W -profile out.txt

  Edge profile f_w1485_g005: ω≈1.485, dQ/dω > 0 on gscan — not Phase N default.
EOF
    ;;
  gscan-light)
    awk -F'\t' 'NR==1 || ($1+0==0.05 && $2+0==1.46){print}' "$GSCAN"
    ;;
  *)
    echo "unknown mode: $MODE" >&2
    exit 2
    ;;
esac
