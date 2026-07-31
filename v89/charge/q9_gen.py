#!/usr/bin/env python3
"""Q9 — weak-link I-V (CHARGE.md §7.5): net + cfg generator.

Finds ONE linked pair on the standard foam (seed 20260727) with link
length d nearest the m=1 pair rung d* = pi/omega at the tuning-curve
load, then emits a family of 2-vertex nets scanning the BARE pitch
difference dw at fixed rung SUM: omega_a + omega_b = 2*pi*C/d exactly
(the lock's sum mode stays on-rung; dw is a pure difference-mode
drive — the Adler/Josephson variable).

DICHOTOMY GUARD: nothing here is a charge. Loads and phases only —
the pitch law x = (w2/omega - 1)/q maps pitches to seeded loads; the
kernel's pin_net holds them; slip_diag counts the existing phase's
holonomy. No new fields, no carriers.
"""
import math, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
FOAM = os.path.join(HERE, "..", "prestress", "foam", "foam_s20260727.tsv")
LAWS = os.path.join(HERE, "..", "battery", "laws_V2g.cfg")
NETS = os.path.join(HERE, "nets")
RUNS = os.path.join(HERE, "runs")
W2, Q, C = 2.9, 1.2, 1.0
LINK_FAC = 1.15

DW_SCAN = [0.00, 0.03, 0.06, 0.09, 0.12, 0.15, 0.18, 0.21, 0.25, 0.30, 0.40, 0.55]
T_RUN = 600.0
T_SMOKE = 20.0


def load_foam():
    cells = []
    with open(FOAM) as f:
        hdr = f.readline().split()
        ix, iy, iz, ir = hdr.index("x"), hdr.index("y"), hdr.index("z"), hdr.index("r")
        for line in f:
            p = line.split()
            if len(p) < 5:
                continue
            cells.append((float(p[ix]), float(p[iy]), float(p[iz]), float(p[ir])))
    return cells


def find_bridge(cells):
    """Linked pair, centered, d nearest the rung window."""
    best = None
    n = len(cells)
    for i in range(n):
        xi, yi, zi, ri = cells[i]
        if not (7.0 < xi < 17.0 and 7.0 < yi < 17.0 and 7.0 < zi < 17.0):
            continue
        for j in range(i + 1, n):
            xj, yj, zj, rj = cells[j]
            if abs(xj - xi) > 2.0 or abs(yj - yi) > 2.0 or abs(zj - zi) > 2.0:
                continue
            d = math.dist((xi, yi, zi), (xj, yj, zj))
            if d >= LINK_FAC * (ri + rj):
                continue
            if not (1.36 <= d <= 1.47):
                continue
            # rung self-consistency: at d, the m=1 tuning-curve load
            om = math.pi * C / d
            x = (W2 / om - 1.0) / Q
            if x < 0.10:
                continue
            score = abs(d - 1.408)
            if best is None or score < best[0]:
                best = (score, i, j, d)
    return best


def emit(i, j, d, cells, dw, name):
    om_bar = math.pi * C / d           # rung sum: (om_a+om_b)*d = 2*pi*C
    om_a, om_b = om_bar + dw / 2.0, om_bar - dw / 2.0
    xa = (W2 / om_a - 1.0) / Q
    xb = (W2 / om_b - 1.0) / Q
    assert xa > 0.07 and xb > 0.07, (dw, xa, xb)
    # phases: forward gate open at seed: th_u - om_a*d/C - th_w = 0
    th_w = 0.0
    th_u = math.fmod(om_a * d / C, 2.0 * math.pi)
    xi, yi, zi, _ = cells[i]
    xj, yj, zj, _ = cells[j]
    net = os.path.join(NETS, f"{name}.net")
    with open(net, "w") as f:
        f.write(f"# Q9 weak link: bridge d={d:.4f} om_bar={om_bar:.4f} dw={dw:.3f}\n")
        f.write(f"# om_a={om_a:.4f} om_b={om_b:.4f} xa={xa:.4f} xb={xb:.4f}\n")
        f.write(f"V {xi:.4f} {yi:.4f} {zi:.4f} {xa:.6f} {th_u:.6f}\n")
        f.write(f"V {xj:.4f} {yj:.4f} {zj:.4f} {xb:.6f} {th_w:.6f}\n")
        f.write("E 0 1\n")
    return net


def emit_cfg(net, name, T):
    cfg = os.path.join(RUNS, f"{name}.cfg")
    with open(cfg, "w") as f:
        f.write(open(LAWS).read())
        f.write(f"""
L=24
T={T}
dt=0.02
seed=20260727
init=net
net_file={os.path.abspath(net)}
pin_net=1
slip_diag=1
edge_sink=0
diag_every=250
snap_every=0
""")
    return cfg


def main():
    os.makedirs(NETS, exist_ok=True)
    os.makedirs(RUNS, exist_ok=True)
    cells = load_foam()
    hit = find_bridge(cells)
    if not hit:
        sys.exit("no bridge pair found")
    _, i, j, d = hit
    om = math.pi * C / d
    print(f"# bridge: cells {i},{j} d={d:.4f} om_bar={om:.4f} x_bar={(W2/om-1)/Q:.4f}")
    jobs = []
    for dw in DW_SCAN:
        name = f"q9_dw{dw:.2f}".replace(".", "p")
        net = emit(i, j, d, cells, dw, name)
        cfg = emit_cfg(net, name, T_RUN)
        jobs.append(cfg)
    for dw, nm in ((0.0, "smoke_dw0"), (0.55, "smoke_dw055")):
        net = emit(i, j, d, cells, dw, nm)
        emit_cfg(net, nm, T_SMOKE)
    print("\n".join(jobs))


if __name__ == "__main__":
    main()
