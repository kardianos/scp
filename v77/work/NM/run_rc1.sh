#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
python3 run_rc1.py
python3 -c "
import json
r=json.load(open('outputs/rc1_result.json'))
print('rc1_claim=', r.get('rc1_claim'))
print('stamps=', r.get('stamps'))
print('RG=', r.get('gates'))
assert 'rc1_claim' in r
"
