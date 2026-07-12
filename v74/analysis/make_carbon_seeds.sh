#!/usr/bin/env bash
# Build Stage 2A multi-ball carbon seeds (liquid-drop fusion).
# Usage: ./make_carbon_seeds.sh [N] [L]
# Defaults: N=192 L=36 (v70/v71 nuclear resolution)
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
V74="$ROOT/v74"
GEN="${ROOT}/bin/gen_qball_multi"
PROF_LIGHT="$V74/profiles/f_w146_g005.txt"
N="${1:-192}"
L="${2:-36}"
OUTDIR="$V74/seeds"
mkdir -p "$OUTDIR"

if [[ ! -x "$GEN" ]]; then
  echo "FATAL: $GEN missing — run: make -C sfa install" >&2
  exit 1
fi

D=10.0
W=1.46
DELTA=0.0

# Emit center lists via Python (octahedron + icosahedron, edge D)
eval "$(python3 - <<PY
import math, itertools
D = $D
# octahedron: edge D → axial a = D/sqrt(2)
a = D / math.sqrt(2)
octa = [
    ( a, 0, 0), (-a, 0, 0),
    (0,  a, 0), (0, -a, 0),
    (0, 0,  a), (0, 0, -a),
]
# icosahedron regular, edge D
phi = (1 + math.sqrt(5)) / 2
verts = []
for s1 in (-1, 1):
    for s2 in (-1, 1):
        verts += [(0.0, s1, s2*phi), (s1, s2*phi, 0.0), (s2*phi, 0.0, s1)]
md = min(math.dist(u, v) for u, v in itertools.combinations(verts, 2))
scale = D / md
ico = [(scale*x, scale*y, scale*z) for x, y, z in verts]
assert len(ico) == 12
def fmt(name, pts):
    parts = []
    for (x,y,z) in pts:
        parts.append(f"{x:.10f} {y:.10f} {z:.10f}")
    print(f'{name}=({" ".join(repr(p) for p in parts)})')
fmt("OCTA", octa)
fmt("ICO", ico)
print(f"A_AXIAL={a:.10f}")
PY
)"

echo "=== c6_light N=$N L=$L D=$D a=$A_AXIAL ==="
args6=()
for c in "${OCTA[@]}"; do
  # shellcheck disable=SC2086
  set -- $c
  args6+=("$PROF_LIGHT" "$W" "$DELTA" "$1" "$2" "$3")
done
"$GEN" "$N" "$L" "$OUTDIR/c6_light_N${N}.sfa" "${args6[@]}"

echo "=== c12_light N=$N L=$L D=$D ==="
args12=()
for c in "${ICO[@]}"; do
  # shellcheck disable=SC2086
  set -- $c
  args12+=("$PROF_LIGHT" "$W" "$DELTA" "$1" "$2" "$3")
done
"$GEN" "$N" "$L" "$OUTDIR/c12_light_N${N}.sfa" "${args12[@]}"

ls -lh "$OUTDIR"/c*_N${N}.sfa
echo "done."
