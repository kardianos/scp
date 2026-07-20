#!/usr/bin/env python3
"""CLI: run T0–T2 for v81 path3_token."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

# allow `python -m src.run_all` or `python src/run_all.py`
_ROOT = Path(__file__).resolve().parent.parent
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))

from src.experiments import run_all, run_t0, run_t1, run_t2  # noqa: E402


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description="v81 path3_token T0–T2 runner")
    p.add_argument(
        "--out",
        type=Path,
        default=_ROOT / "results",
        help="results directory (default: path3_token/results)",
    )
    p.add_argument(
        "--only",
        choices=("T0", "T1", "T2", "all"),
        default="all",
        help="which experiment block to run",
    )
    args = p.parse_args(argv)
    out = args.out
    out.mkdir(parents=True, exist_ok=True)

    if args.only == "all":
        summary = run_all(out)
        return 0 if summary["all_passed"] else 1
    if args.only == "T0":
        r = run_t0(out / "T0")
        print(r)
        return 0 if r.passed else 1
    if args.only == "T1":
        r = run_t1(out / "T1")
        print(r)
        return 0 if r.passed else 1
    if args.only == "T2":
        r = run_t2(out / "T2")
        print(r)
        return 0 if r.passed else 1
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
