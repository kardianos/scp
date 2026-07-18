#!/usr/bin/env python3
"""Run RC1 co-field and exit 0 iff rc1_claim."""
import os
import sys

os.chdir(os.path.dirname(os.path.abspath(__file__)))
sys.argv = [sys.argv[0], "--nx", "48", "--steps", "24"]
from sandbox_rc1_cofield import main

r = main()
print("RC1 DONE rc1_claim=", r["rc1_claim"], "stamps=", r["stamps"])
raise SystemExit(0 if r["rc1_claim"] else 1)
