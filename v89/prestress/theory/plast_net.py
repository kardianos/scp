#!/usr/bin/env python3
"""PLAST network model — does the d-flow converge to exact consonance on
shells, and what does it do to parasites?

Cells carry phases th_k advancing at uniform pitch om0 (loads fixed — the
geometry question isolated); links carry plastic lengths d_l. Dynamics:
  entrain   th_j <- pulled toward retarded tail th_i - om d/C at rate
            kth * G(psi_f); symmetrically th_i.
  plastic   dd/dt = -kp * geo(d) * dV/dd,  V = 1 - G(psi_f)G(psi_b)
            (uniform amplitudes: the gate-free urge reduces to geo(d);
             geo -> 0 at the pinch d = 2 r0 — sealed links feel no force)

Graphs: rings N=3..6 (odd = frustrated), cube (12 edges), cube + candidate
face diagonals (the H1 parasite population). d seeded with the foam's
measured +-15% jitter.
No Euclidean embedding is enforced on d — v89 has no background; d is a
link retardation, constrained only by cycle closure integers.
"""
import numpy as np

W2, Q, C, PG = 2.9, 1.2, 1.0, 8.0
R0 = 0.85
DT = 0.02
KTH = 0.35          # entrainment rate (locks ~10 t.u., e6-like)
OM0 = W2/(1+Q*0.32) # 2.0954, d* = pi/om = 1.4993

def wrap(a): return (a + np.pi) % (2*np.pi) - np.pi
def gate(p): return (0.5*(1+np.cos(p)))**PG
def dgate(p):
    h = 0.5*(1+np.cos(p))
    return -0.5*PG*h**(PG-1)*np.sin(p)
def geo(d):
    if d >= 2*R0: return 0.0
    a2 = R0*R0 - 0.25*d*d
    if a2 <= 0: return 0.0
    return (np.pi*a2/(np.pi*R0*R0))*(2*R0/d)

def evolve(V, E, d0, kp, T=600.0, seed=3, th0=None, dstar=None):
    rng = np.random.default_rng(seed)
    th = rng.uniform(0, 2*np.pi, V) if th0 is None else th0.copy()
    d = d0.copy()
    nE = len(E)
    for s in range(int(T/DT)):
        pf = np.array([wrap(th[i] - OM0*d[l]/C - th[j]) for l,(i,j) in enumerate(E)])
        pb = np.array([wrap(th[j] - OM0*d[l]/C - th[i]) for l,(i,j) in enumerate(E)])
        gf, gb = gate(pf), gate(pb)
        dth = np.zeros(V)
        for l,(i,j) in enumerate(E):
            g0 = geo(d[l])
            if g0 <= 0: continue
            dth[j] += KTH*g0*gf[l]*pf[l]*DT
            dth[i] += KTH*g0*gb[l]*pb[l]*DT
            if kp > 0:
                dV = (OM0*dgate(pf[l])*gb[l] + OM0*gf[l]*dgate(pb[l]))/C
                d[l] += -kp*g0*dV*DT
                if d[l] < 0.5: d[l] = 0.5
        th = (th + OM0*DT + dth) % (2*np.pi)
    pf = np.array([wrap(th[i] - OM0*d[l]/C - th[j]) for l,(i,j) in enumerate(E)])
    pb = np.array([wrap(th[j] - OM0*d[l]/C - th[i]) for l,(i,j) in enumerate(E)])
    return th, d, pf, pb

def report(tag, E, d, pf, pb, cycles=None, split=None):
    gg = gate(pf)*gate(pb)
    live = np.array([geo(x) > 0 for x in d])
    churn = sum(geo(d[l])*(1-gg[l]) for l in range(len(E)) if live[l])
    msg = (f"# {tag}: gg min={gg[live].min() if live.any() else 0:.4f} "
           f"mean={gg[live].mean() if live.any() else 0:.4f} "
           f"max|psi_f|={np.abs(pf[live]).max() if live.any() else 0:.4f} "
           f"churn={churn:.4f}")
    if split is not None:
        a, b = split
        gga = gg[a]; ggb = gg[b]; livb = live[b]
        msg += (f" | edges gg={gga.mean():.4f} (min {gga.min():.4f})"
                f" | parasites: live={int(livb.sum())}/{len(b)}"
                f" sealed={int((~livb).sum())}")
        if livb.any():
            msg += f" live_gg={ggb[livb].mean():.3f} live_d={[round(d[l],2) for l in np.array(b)[livb]]}"
    print(msg)
    if cycles:
        for cyc in cycles:
            s = sum(OM0*d[l]/C for l in cyc)/(2*np.pi)
            print(f"#     cycle closure sum(om d/C)/2pi = {s:.4f}")

DSTAR = np.pi*C/OM0
print(f"# net model: om0={OM0:.4f} d*={DSTAR:.4f} pinch=2r0={2*R0} p_gate={PG:.0f}")
rng = np.random.default_rng(20260728)

print("\n# ==== rings: even = exact consonance reachable; odd = frustrated ====")
for N in (3, 4, 5, 6):
    Vv = N
    E = [(k, (k+1) % N) for k in range(N)]
    d0 = DSTAR*(1 + rng.uniform(-0.15, 0.15, N))
    th, d, pf, pb = evolve(Vv, E, d0, kp=0.5, T=600)
    gg = gate(pf)*gate(pb)
    pred = np.pi/(2*N)   # equal-split prediction for the odd-cycle defect (per gate)
    tagx = f"ring{N} kp=0.5"
    report(tagx, E, d, pf, pb, cycles=[list(range(N))])
    if N % 2 == 1:
        print(f"#     odd-N prediction: per-link |psi| = pi/(2N) = {pred:.4f}, "
              f"gg = {gate(pred)**2:.4f}  (measured above)")
    c0 = evolve(Vv, E, d0, kp=0.0, T=600)
    gg0 = gate(c0[2])*gate(c0[3])
    print(f"#     frozen control (kp=0): gg min={gg0.min():.4f} mean={gg0.mean():.4f}")

print("\n# ==== convergence rate vs kappa_p (ring6, 15% jitter) ====")
for kp in (0.05, 0.2, 0.5, 2.0, 8.0):
    E = [(k, (k+1) % 6) for k in range(6)]
    d0 = DSTAR*(1 + rng.uniform(-0.15, 0.15, 6))
    # time to all-gates > 0.99
    th = np.random.default_rng(5).uniform(0, 2*np.pi, 6)
    d = d0.copy(); t99 = np.nan
    for s in range(int(600/DT)):
        pf = np.array([wrap(th[i]-OM0*d[l]/C-th[j]) for l,(i,j) in enumerate(E)])
        pb = np.array([wrap(th[j]-OM0*d[l]/C-th[i]) for l,(i,j) in enumerate(E)])
        gf, gb = gate(pf), gate(pb)
        if np.isnan(t99) and (gf*gb > 0.99).all(): t99 = s*DT
        dth = np.zeros(6)
        for l,(i,j) in enumerate(E):
            g0 = geo(d[l])
            if g0 <= 0: continue
            dth[j] += KTH*g0*gf[l]*pf[l]*DT
            dth[i] += KTH*g0*gb[l]*pb[l]*DT
            dV = (OM0*dgate(pf[l])*gb[l] + OM0*gf[l]*dgate(pb[l]))/C
            d[l] += -kp*g0*dV*DT
            if d[l] < 0.5: d[l] = 0.5
        th = (th + OM0*DT + dth) % (2*np.pi)
    print(f"#   kp={kp:4}: t(all gg>0.99) = {'%.0f' % t99 if not np.isnan(t99) else 'never (600)'}")

print("\n# ==== the cube, edges only (bipartite benchmark) ====")
CUBEV = 8
CE = [(k, k ^ (1 << b)) for k in range(8) for b in range(3) if (k ^ (1 << b)) > k]
faces = []
for ax in range(3):
    for s0 in (0, 1):
        vs = [k for k in range(8) if (k >> ax) & 1 == s0]
        a, b_, c, dd = vs
        cyc = [a, b_, dd, c]
        fl = []
        for q in range(4):
            u, v = cyc[q], cyc[(q+1) % 4]
            for l, (i, j) in enumerate(CE):
                if {i, j} == {u, v}: fl.append(l)
        if len(fl) == 4: faces.append(fl)
d0 = DSTAR*(1 + rng.uniform(-0.15, 0.15, 12))
th, d, pf, pb = evolve(CUBEV, CE, d0, kp=0.5, T=600)
report("cube kp=0.5", CE, d, pf, pb, cycles=faces[:2])
c0 = evolve(CUBEV, CE, d0, kp=0.0, T=600)
gg0 = gate(c0[2])*gate(c0[3])
print(f"#     frozen control: gg min={gg0.min():.4f} mean={gg0.mean():.4f}")
print(f"#     plastic d spread: {d.std()/d.mean()*100:.2f}%  (seeded {d0.std()/d0.mean()*100:.1f}%)")

print("\n# ==== cube + parasitic face diagonals (the H1 population) ====")
# a=1.25-class cube: edges near d*, diagonals near sqrt(2) d* with jitter;
# candidate ceiling 1.15*2r = 1.955, live if d < 1.7
for trial in range(3):
    d0e = DSTAR*(1 + rng.uniform(-0.15, 0.15, 12))
    DIAG = []
    for k in range(8):
        for k2 in range(k+1, 8):
            if bin(k ^ k2).count("1") == 2: DIAG.append((k, k2))
    d0d = np.sqrt(2)*DSTAR*0.85*(1 + rng.uniform(-0.15, 0.15, len(DIAG)))  # jitter-shortened population
    keep = d0d < 1.955
    DIAG = [DIAG[q] for q in range(len(DIAG)) if keep[q]]
    d0d = d0d[keep]
    E = CE + DIAG
    d0 = np.concatenate([d0e, d0d])
    nliv0 = int((d0d < 1.7).sum())
    th, d, pf, pb = evolve(8, E, d0, kp=0.5, T=1200, seed=30+trial)
    report(f"cube+diag t{trial} (seeded: {len(DIAG)} candidates, {nliv0} live, "
           f"d_par={np.round(d0d,2)})", E, d, pf, pb,
           split=(list(range(12)), list(range(12, len(E)))))
    # where did the live parasites go?
    for q in range(12, len(E)):
        if geo(d0[q]) > 0:
            print(f"#     parasite {E[q]}: d {d0[q]:.3f} -> {d[q]:.3f} "
                  f"({'SEALED' if geo(d[q])<=0 else 'live gg=%.3f' % (gate(pf[q])*gate(pb[q]))})")
