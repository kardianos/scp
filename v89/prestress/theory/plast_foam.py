#!/usr/bin/env python3
"""PLAST diagnosis — is the frozen foam the obstruction? (P15, prestress/PLASTICITY.md)

Quantifies, from the real foam geometry (foam_s20260727.tsv, the standard
campaign foam), what strut-length uniformity the FROZEN foam can supply for
small consonant shells, the gate cost of the residual spread, and therefore
what a plasticity law would have to deliver (delta-d magnitudes, timescale).

Pure geometry + the kernel's own formulas (gate_of, rung, lens overlap).
No simulator run. Kernel formulas transcribed from v89/cellfab.c:
  candidate link:  d < 1.15*(r_i + r_j)          (build_field)
  live channel:    d < r_i + r_j  (lens overlap A > 0)
  gate:            G(psi) = ((1+cos psi)/2)^p_gate,  p_gate = 8
  pair rung (1:1): (w_i + w_j) d / C = 2 pi m ;  equal pair: w d / C = pi m
  pitch:           w = w2 / (1 + q x),  w2 = 2.9, q = 1.2, C = 1
"""
import numpy as np
import sys

FOAM = "/home/d/code/scp/v89/prestress/foam/foam_s20260727.tsv"
W2, Q, C, PGATE = 2.9, 1.2, 1.0, 8.0
L = 24.0

def gate(psi):
    return (0.5 * (1.0 + np.cos(psi)))**PGATE

def wrap(a):
    return (a + np.pi) % (2*np.pi) - np.pi

d = np.loadtxt(FOAM, skiprows=1)
x, y, z, r = d[:,1], d[:,2], d[:,3], d[:,4]
N = len(x)
print(f"# foam: {N} cells, box L={L}, r mean={r.mean():.4f} spread=[{r.min():.4f},{r.max():.4f}]")

# ---- link table by the kernel rule (O(N^2) is fine at N=9741) ----
# candidate: d < 1.15 (ri + rj); live: d < ri + rj
P = np.stack([x, y, z], axis=1)
links_i, links_j, links_d = [], [], []
# cell binning for speed
from collections import defaultdict
bs = 2.3   # > max cut 1.15*2*rmax ~ 2.09
bins = defaultdict(list)
for i in range(N):
    bins[(int(x[i]/bs), int(y[i]/bs), int(z[i]/bs))].append(i)
for i in range(N):
    bi = (int(x[i]/bs), int(y[i]/bs), int(z[i]/bs))
    for ox in (-1,0,1):
        for oy in (-1,0,1):
            for oz in (-1,0,1):
                for j in bins.get((bi[0]+ox, bi[1]+oy, bi[2]+oz), ()):
                    if j <= i: continue
                    dd = np.sqrt((x[i]-x[j])**2 + (y[i]-y[j])**2 + (z[i]-z[j])**2)
                    if dd < 1.15*(r[i]+r[j]):
                        links_i.append(i); links_j.append(j); links_d.append(dd)
li = np.array(links_i); lj = np.array(links_j); ld = np.array(links_d)
live = ld < (r[li] + r[lj])
print(f"# links: candidates={len(ld)} (mean_deg={2*len(ld)/N:.2f}), live={live.sum()} "
      f"(mean_deg={2*live.sum()/N:.2f})")
print(f"# link lengths: dbar={ld.mean():.4f} sigma={ld.std():.4f} "
      f"[{ld.min():.3f},{ld.max():.3f}]  live: dbar={ld[live].mean():.4f} sigma={ld[live].std():.4f}")

# fractional jitter of live links
print(f"# live-link fractional spread sigma_d/dbar = {ld[live].std()/ld[live].mean():.4f}")

# ---- the rung landscape on the raw foam ----
# a uniform-omega structure at pitch w picks up per-link rung offset
# delta_l = 2 w d_l / C - 2 pi (m=1).  If the builder tunes w to the MEAN
# link length of the chosen strut set (w = pi C / dbar_set), the residual
# delta is set purely by the SPREAD:  delta_l = 2 pi (d_l/dbar - 1).
# Gate cost of a fractional spread eps: G^2(2 pi eps).
for eps in (0.01, 0.02, 0.05, 0.075, 0.10, 0.15):
    dl = 2*np.pi*eps
    print(f"#   spread {eps*100:5.1f}%  ->  delta={dl:.3f} rad  G^2={gate(dl/2)**2:.4f}")

# ---- cube pick experiment: what does the foam actually supply? ----
def pick_cube(cx0, cy0, cz0, a, refine=3, rng=None, exclude_r=0.0):
    """Greedy nearest-vertex pick + edge-uniformity refinement, transcribed
    from cellfab.c init=shell (the H-series seeder)."""
    h = 0.5*a
    tgt = np.array([[cx0+((k&1)*2-1)*h, cy0+((k&2)//2*2-1)*h, cz0+((k&4)//4*2-1)*h]
                    for k in range(8)])
    pick = np.full(8, -1, dtype=int)
    for k in range(8):
        dd = np.sum((P - tgt[k])**2, axis=1)
        if exclude_r > 0:
            core = np.sum((P - np.array([cx0,cy0,cz0]))**2, axis=1) < exclude_r**2
            dd[core] = 1e30
        for q in pick[:k]:
            if q >= 0: dd[q] = 1e30
        pick[k] = np.argmin(dd)
    # refinement passes (keep swaps that reduce edge-length spread)
    for _ in range(refine):
        for k in range(8):
            cand = np.where(np.sum((P - tgt[k])**2, axis=1) < 1.44)[0]
            best_sc, best_c = 1e30, pick[k]
            for i in cand:
                if exclude_r > 0:
                    if np.sum((P[i] - np.array([cx0,cy0,cz0]))**2) < exclude_r**2:
                        continue
                if any(pick[q] == i for q in range(8) if q != k):
                    continue
                sc = 0.0
                for b in range(3):
                    k2 = k ^ (1 << b)
                    e = np.linalg.norm(P[pick[k2]] - P[i])
                    sc += (e - a)**2
                if sc < best_sc:
                    best_sc, best_c = sc, i
            pick[k] = best_c
    return pick

def cube_report(pick, a, label, verbose=True):
    """Edge stats + seed-gate model with BFS phase seeding (as the kernel does):
    tree edges get exact forward gates; co-tree edges carry the cycle defects."""
    edges = []
    for k in range(8):
        for b in range(3):
            k2 = k ^ (1 << b)
            if k2 < k: continue
            edges.append((k, k2, np.linalg.norm(P[pick[k]] - P[pick[k2]])))
    dlen = np.array([e[2] for e in edges])
    abar = dlen.mean()
    om = np.pi*C/abar                       # pi-rung pitch at the mean edge
    xload = (W2/om - 1.0)/Q                 # tuning-curve load
    # per-edge rung offset (both-gate condition): delta = 2 w d / C - 2 pi
    delta = 2*om*dlen/C - 2*np.pi
    # BFS phase seeding from vertex 0 over cube edges (kernel order)
    th = np.full(8, np.nan); th[0] = 0.0
    seen = [0]; queue = [0]
    tree = set()
    while queue:
        k = queue.pop(0)
        for b in range(3):
            k2 = k ^ (1 << b)
            if k2 in seen: continue
            dd = np.linalg.norm(P[pick[k2]] - P[pick[k]])
            th[k2] = (th[k] - om*dd/C) % (2*np.pi)
            tree.add((min(k,k2), max(k,k2)))
            seen.append(k2); queue.append(k2)
    gates_f = []
    for (k, k2, dd) in edges:
        psi_f = wrap(th[k] - om*dd/C - th[k2])
        gates_f.append(gate(psi_f))
    gates_f = np.array(gates_f)
    # parasites: face diagonals inside the candidate ceiling
    npar, par_d = 0, []
    for kA in range(8):
        for kB in range(kA+1, 8):
            hb = bin(kA ^ kB).count("1")
            if hb != 2: continue            # face diagonal
            dd = np.linalg.norm(P[pick[kA]] - P[pick[kB]])
            if dd < 1.15*(r[pick[kA]] + r[pick[kB]]):
                npar += 1; par_d.append(dd)
    if verbose:
        print(f"# cube a={a:<5} {label}: abar={abar:.3f} spread={dlen.std():.3f} "
              f"({100*dlen.std()/abar:.1f}%) range=[{dlen.min():.3f},{dlen.max():.3f}]")
        print(f"#   pitch om={om:.4f} load x={xload:.3f}  "
              f"delta: mean|.|={np.abs(delta).mean():.3f} max|.|={np.abs(delta).max():.3f} rad")
        print(f"#   equal-split gate model G^2(delta/2)^2 -> mean={np.mean(gate(delta/2)**2):.3f} "
              f"min={np.min(gate(delta/2)**2):.3f}")
        print(f"#   BFS-seeded forward gates: mean={gates_f.mean():.3f} min={gates_f.min():.4f} "
              f"(co-tree edges carry cycle defects)")
        print(f"#   parasitic in-ceiling diagonals: {npar}"
              + (f"  d=[{min(par_d):.2f},{max(par_d):.2f}]" if npar else ""))
        print(f"#   plasticity must deliver |dd| per edge: mean={np.abs(dlen-abar).mean():.3f} "
              f"max={np.abs(dlen-abar).max():.3f}  (in d units)")
    return dlen, delta, gates_f, npar

print("\n# ==== the H-series cube, center of box (the actual seeds used) ====")
for a, excl in ((1.25, 0.0), (1.5, 0.0), (1.5, 0.8)):
    lbl = "core-excluded" if excl > 0 else "plain"
    pick = pick_cube(L/2, L/2, L/2, a, exclude_r=excl)
    cube_report(pick, a, lbl)

print("\n# ==== ensemble: 40 cube placements across the box (a=1.5) ====")
rng = np.random.default_rng(20260727)
sp, mg, gmin_all, npar_all = [], [], [], []
for t in range(40):
    c0 = rng.uniform(6, L-6, size=3)
    pick = pick_cube(*c0, 1.5)
    dlen, delta, gates_f, npar = cube_report(pick, 1.5, "", verbose=False)
    sp.append(dlen.std()/dlen.mean()); mg.append(gate(delta/2).mean()**2)
    gmin_all.append(gate(delta/2).min()**2); npar_all.append(npar)
sp = np.array(sp); mg = np.array(mg); gmin_all = np.array(gmin_all)
print(f"# edge spread: median={np.median(sp)*100:.1f}% best={sp.min()*100:.1f}% "
      f"worst={sp.max()*100:.1f}%")
print(f"# equal-split G^2: median mean-gate={np.median(mg):.3f}; median min-gate={np.median(gmin_all):.3f}")
print(f"# parasite count: median={np.median(npar_all):.0f} max={max(npar_all)}")

# ---- best equal-length bipartite subgraph the foam offers at all (P3 view):
# how tight can ANY 12-link matched set get? Sample link-length percentile
# window widths needed to harvest 12 links from one neighbourhood.
print("\n# ==== the ladder cost of the residual, in lifetime terms ====")
# leak attribution (MASS.md M0): rough fraction of an off-comb delivery
# R = 2|det|G_r/(det^2+G_r^2) — but for a SEEDED consonant structure the
# relevant leak is gate churn: the per-link roughness radiated scales with
# traffic x (1 - G^2). Residual spread eps -> (1-G^2(2 pi eps)):
for eps in (0.005, 0.01, 0.02, 0.05, 0.075):
    dl = 2*np.pi*eps
    print(f"#   spread {eps*100:5.1f}% -> residual churn factor (1-G^2) = {1-gate(dl/2)**2:.4f}")
print("# (the structural floor: plasticity that halves the spread quarters the")
print("#  small-delta churn, since 1-G^2 ~ p_gate * delta^2 / 4 at small delta)")
