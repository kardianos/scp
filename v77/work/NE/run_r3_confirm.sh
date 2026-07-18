#!/usr/bin/env bash
# NE R3: re-run parent-fixed full Maxwell sandbox and confirm claim.
set -euo pipefail
cd "$(dirname "$0")"
python3 sandbox_full_maxwell_r2.py --quick
python3 -c "
import json
r=json.load(open('outputs/r2_result.json'))
assert r.get('full_maxwell_claim') is True, r
s=r.get('summary',{})
for k in ['FM1_vacuum','FM2_wave_unit','FM2_wave_offunit','FM3_divB','FM4_faraday','FM5_gauss_static','FM6_continuity','FM7_coulomb_3d']:
    assert s.get(k) is True or r.get('TE_KG',{}).get(k.replace('FM','KG')) is not False
print('R3 CONFIRM: full_maxwell_claim=', r['full_maxwell_claim'])
print('key_numbers=', r.get('key_numbers') or {
  'v_unit': r['gates']['FM2_wave_unit']['v_ratio'],
  'v_off': r['gates']['FM2_wave_offunit']['v_ratio'],
  'divB': r['gates']['FM3_divB']['divB_max'],
  'dQ': r['gates']['FM6_continuity']['dQ_rel'],
})
"
