#!/usr/bin/env python3
"""
Run full M1 suite (same as sandbox_m1_2d.main) and write outputs.
Prefer: python3 sandbox_m1_2d.py --quick
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
os.chdir(os.path.dirname(os.path.abspath(__file__)))
sys.argv = ["offline_m1_emit.py", "--quick"]
from sandbox_m1_2d import main

if __name__ == "__main__":
    r = main()
    print("m1_claim=", r["m1_claim"])
    sys.exit(0 if r["m1_claim"] else 1)
