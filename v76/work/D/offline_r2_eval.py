#!/usr/bin/env python3
"""
Offline Round-2 evaluation with dense multi-start (stdlib).
Produces round2_* tables + feeds congruence_verdict_r2.md numbers.
"""
from __future__ import annotations

import json
import math
from pathlib import Path

# Import maps and forward models from congruence_score_r2
import congruence_score_r2 as R

OUT = Path(__file__).resolve().parent / "results"
OUT.mkdir(exist_ok=True)


def main():
    report = R.run()
    R.save(report)

    # Also write a compact model comparison for verdict
    rows = []
    for ev in report["evaluations"]:
        for s in ev["scores"]:
            rows.append((ev["map_tag"], ev["channel"], s["name"], s["L_fit"], s["S"]))

    with open(OUT / "round2_model_rank.tsv", "w") as f:
        f.write("map\tchannel\tmodel\tL_fit\tS\trank_S\n")
        # group by map
        from collections import defaultdict

        by = defaultdict(list)
        for r in rows:
            by[(r[0], r[1])].append(r)
        for key, lst in by.items():
            lst_sorted = sorted(lst, key=lambda x: x[4])
            for i, r in enumerate(lst_sorted):
                f.write(
                    f"{r[0]}\t{r[1]}\t{r[2]}\t{r[3]:.8e}\t{r[4]:.8e}\t{i+1}\n"
                )

    print("offline_r2_eval done")
    print(json.dumps(report["summary"], indent=2))
    return report


if __name__ == "__main__":
    main()
