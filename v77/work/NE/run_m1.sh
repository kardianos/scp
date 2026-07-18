#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
python3 sandbox_m1_2d.py --quick
python3 -c "
import json
r=json.load(open('outputs/m1_result.json'))
print('m1_claim=', r.get('m1_claim'))
for k,v in r.get('mandatory_pass',{}).items():
    print(f'  {k}: {v}')
assert 'm1_claim' in r
"
