#!/usr/bin/env python3
"""Execute M1 suite and print m1_claim. Parent/local: python3 _run_m1_once.py"""
import sandbox_m1_2d as m1

if __name__ == "__main__":
    import sys

    sys.argv = [sys.argv[0], "--quick"]
    r = m1.main()
    print("DONE m1_claim=", r["m1_claim"])
    raise SystemExit(0 if r["m1_claim"] else 1)
