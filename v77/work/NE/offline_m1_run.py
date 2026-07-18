#!/usr/bin/env python3
"""Offline/local M1 runner — identical to sandbox_m1_2d.main --quick."""
import sys
import os

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.argv = ["offline_m1_run.py", "--quick"]
from sandbox_m1_2d import main

if __name__ == "__main__":
    r = main()
    sys.exit(0 if r.get("m1_claim") else 2)
