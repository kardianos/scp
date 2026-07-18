#!/usr/bin/env python3
"""Offline fallback: import sandbox_rc1_cofield.main if full run desired.

Prefer: python3 run_rc1.py
This file exists so parent/agent can `python3 _emit_rc1_offline.py`.
"""
from sandbox_rc1_cofield import main

if __name__ == "__main__":
    main()
