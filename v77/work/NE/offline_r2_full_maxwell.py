#!/usr/bin/env python3
"""NE R2 offline runner — invokes sandbox_full_maxwell_r2 with --quick."""
from __future__ import annotations

import sys
import os

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

if __name__ == "__main__":
    sys.argv = [sys.argv[0], "--quick"]
    from sandbox_full_maxwell_r2 import main

    main()
