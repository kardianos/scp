#!/usr/bin/env python3
# PRESTRESS scorer — one log -> one verdict row (appended to SCORES.tsv).
#
# Definitions (frozen for the campaign; same for every candidate):
#   m0        seeded free dense mass (first # LUMP Emfree)
#   t_death   first t where the largest lump's mass < 0.3*m0
#             (deep decay; consistent with the skirt-death phenomenology);
#             blank = censored (alive at end)
#   late_dmdt linear slope of Emfree over the trailing 500 t.u. of the run
#   rough_share  rough/(rough+evap+back_s) from the final # CONV ledger
#   gates/parasites echoed from the # net: seed line
#
#   score_net.py RUN.log [RUN2.log ...]

import os
import re
import sys

SCORES = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                      "runs", "SCORES.tsv")
HDR = ("name\tnv\tne\tgate_min\tgate_mean\tparasites\tgpar_max\tm0\t"
       "t_end\tt_death\tlate_dmdt\trough_share\talive\n")


def score(path):
    name = os.path.basename(path).rsplit(".", 1)[0]
    nv = ne = par = None
    gmin = gmean = gpar = None
    lumps = []          # (t, emfree, biggest_m)
    conv = None
    with open(path) as f:
        for line in f:
            if line.startswith("# net:"):
                m = re.search(r"nv=(\d+) ne=(\d+) worst_pick=\S+ gates "
                              r"min=([\d.eE+-]+) mean=([\d.eE+-]+) "
                              r"parasites=(\d+) gpar_max=([\d.eE+-]+)", line)
                if m:
                    nv, ne = int(m.group(1)), int(m.group(2))
                    gmin, gmean = float(m.group(3)), float(m.group(4))
                    par, gpar = int(m.group(5)), float(m.group(6))
            elif line.startswith("# LUMP"):
                m = re.search(r"t=([\d.]+) n=(\d+) Emfree=([\d.eE+-]+)", line)
                if not m:
                    continue
                t, emf = float(m.group(1)), float(m.group(3))
                big = 0.0
                for mm in re.finditer(r"\| m=([\d.eE+-]+)", line):
                    big = max(big, float(mm.group(1)))
                lumps.append((t, emf, big))
            elif line.startswith("# CONV"):
                m = re.search(r"cond=(\S+) evap=(\S+) rough=(\S+) back_s=(\S+)",
                              line)
                if m:
                    conv = tuple(float(x) for x in m.groups())
    if not lumps:
        return f"{name}\tNO-LUMP-DATA"
    m0 = lumps[0][2] if lumps[0][2] > 0 else lumps[0][1]
    t_end = lumps[-1][0]
    t_death = ""
    for t, emf, big in lumps:
        if big < 0.3 * m0:
            t_death = f"{t:.0f}"
            break
    tail = [(t, emf) for t, emf, _ in lumps if t >= t_end - 500]
    if len(tail) >= 3:
        n = len(tail)
        sx = sum(t for t, _ in tail); sy = sum(e for _, e in tail)
        sxx = sum(t * t for t, _ in tail); sxy = sum(t * e for t, e in tail)
        den = n * sxx - sx * sx
        slope = (n * sxy - sx * sy) / den if den else 0.0
    else:
        slope = 0.0
    rough_share = ""
    if conv:
        loss = conv[1] + conv[2] + conv[3]
        if loss > 0:
            rough_share = f"{conv[2] / loss:.3f}"
    alive = "1" if not t_death else "0"
    fmt = lambda v: "" if v is None else (f"{v:.4f}" if isinstance(v, float) else str(v))
    return (f"{name}\t{fmt(nv)}\t{fmt(ne)}\t{fmt(gmin)}\t{fmt(gmean)}\t"
            f"{fmt(par)}\t{fmt(gpar)}\t{m0:.4f}\t{t_end:.0f}\t{t_death}\t"
            f"{slope:.3e}\t{rough_share}\t{alive}")


def main():
    rows = [score(p) for p in sys.argv[1:]]
    new = not os.path.exists(SCORES)
    with open(SCORES, "a") as f:
        if new:
            f.write(HDR)
        for r in rows:
            f.write(r + "\n")
    for r in rows:
        print(r)


if __name__ == "__main__":
    main()
