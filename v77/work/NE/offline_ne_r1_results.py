#!/usr/bin/env python3
"""
Offline / always-runnable result package for NE R1.

Runs the same physics as sandbox_ne_r1_em.py with defaults tuned for
stdlib pure-Python on a laptop (N=20 SOR, modest wave). Prefer:

  python3 sandbox_ne_r1_em.py

when a shell is available; this file is the portable twin and is
executed by `run_and_export.py` pattern: import + main.

If imported as module, call main().
"""
from __future__ import annotations

import os
import sys

# Re-use full sandbox as single source of truth
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from sandbox_ne_r1_em import main  # noqa: E402

if __name__ == "__main__":
    # Slightly smaller defaults for quick offline
    sys.argv = [
        sys.argv[0],
        "--N", "20",
        "--iters", "450",
        "--wave-nx", "301",
        "--wave-steps", "500",
    ]
    main()
