#!/usr/bin/env python3
"""PLAST network model, round 2 — the decisive cases with kernel-faithful
protocol: (a) all-live jitter, (b) seeded-lock phases (kernel practice),
(c) ambient phase noise (the kernel's churn/tumble analog), (d) lock-time
hardening (option iv). Companion to plast_net.py; see PLASTICITY.md.
"""
import numpy as np

W2, Q, C, PG = 2.9, 1.2, 1.0, 8.0
R0 = 0.85
DT = 0.02
KTH = 0.35
OM0 = W2/(1+Q*0.32)
DSTAR = np.pi*C/OM0

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

def bfs_phases(V, E, d0):
    """seedlock along a spanning tree (kernel shell/net practice)"""
    th = np.full(V, np.nan); th[0] = 0.0
    adj = {k: [] for k in range(V)}
    for l, (i, j) in enumerate(E):
        adj[i].append((j, l)); adj[j].append((i, l))
    q = [0]
    while q:
        u = q.pop(0)
        for (v, l) in adj[u]:
            if not np.isnan(th[v]): continue
            th[v] = (th[u] - OM0*d0[l]/C) % (2*np.pi)
            q.append(v)
    th[np.isnan(th)] = 0.0
    return th

def evolve(V, E, d0, kp, T=600.0, seed=3, th0=None, noise=0.0,
           harden_tau=0.0, gg_thr=0.9):
    rng = np.random.default_rng(seed)
    th = rng.uniform(0, 2*np.pi, V) if th0 is None else th0.copy()
    d = d0.copy(); nE = len(E)
    tlock = np.zeros(nE)          # accumulated locked time (hardening)
    t99 = np.nan
    for s in range(int(T/DT)):
        pf = np.array([wrap(th[i] - OM0*d[l]/C - th[j]) for l,(i,j) in enumerate(E)])
        pb = np.array([wrap(th[j] - OM0*d[l]/C - th[i]) for l,(i,j) in enumerate(E)])
        gf, gb = gate(pf), gate(pb)
        liv = np.array([geo(x) > 0 for x in d])
        if np.isnan(t99) and liv.any() and (gf[liv]*gb[liv] > 0.99).all():
            t99 = s*DT
        dth = np.zeros(V)
        for l, (i, j) in enumerate(E):
            g0 = geo(d[l])
            if g0 <= 0: continue
            dth[j] += KTH*g0*gf[l]*pf[l]*DT
            dth[i] += KTH*g0*gb[l]*pb[l]*DT
            if kp > 0:
                kp_l = kp/(1.0 + tlock[l]/harden_tau) if harden_tau > 0 else kp
                dV = (OM0*dgate(pf[l])*gb[l] + OM0*gf[l]*dgate(pb[l]))/C
                d[l] += -kp_l*g0*dV*DT
                if d[l] < 0.5: d[l] = 0.5
            if gf[l]*gb[l] > gg_thr: tlock[l] += DT
            else: tlock[l] = max(0.0, tlock[l] - DT)
        if noise > 0:
            dth += noise*np.sqrt(DT)*rng.standard_normal(V)
        th = (th + OM0*DT + dth) % (2*np.pi)
    pf = np.array([wrap(th[i] - OM0*d[l]/C - th[j]) for l,(i,j) in enumerate(E)])
    pb = np.array([wrap(th[j] - OM0*d[l]/C - th[i]) for l,(i,j) in enumerate(E)])
    return th, d, pf, pb, t99, tlock

def stat(E, d, pf, pb):
    gg = gate(pf)*gate(pb)
    liv = np.array([geo(x) > 0 for x in d])
    churn = sum(geo(d[l])*(1-gg[l]) for l in range(len(E)) if liv[l])
    return gg, liv, churn

print(f"# net2: om0={OM0:.4f} d*={DSTAR:.4f} pinch={2*R0}")

print("\n# ==== (A) ring6 all-live jitter, random phases + ambient noise: rate vs kp ====")
print("# kp     t99(noise=0)  t99(noise=0.02)   [3 seeds each; '-' = >600]")
rng = np.random.default_rng(7)
for kp in (0.05, 0.2, 0.5, 2.0, 8.0):
    r0l, r1l = [], []
    for sd in range(3):
        E = [(k, (k+1) % 6) for k in range(6)]
        d0 = DSTAR*(1 + np.random.default_rng(40+sd).uniform(-0.12, 0.12, 6))
        _,_,_,_, t0, _ = evolve(6, E, d0, kp, T=600, seed=200+sd, noise=0.0)
        _,_,_,_, t1, _ = evolve(6, E, d0, kp, T=600, seed=200+sd, noise=0.02)
        r0l.append('%.0f' % t0 if not np.isnan(t0) else '-')
        r1l.append('%.0f' % t1 if not np.isnan(t1) else '-')
    print(f"#  {kp:4}   {','.join(r0l):>12}   {','.join(r1l):>12}")

print("\n# ==== (B) cube edges, SEEDED-LOCK phases (kernel practice), live jitter ====")
CE = [(k, k ^ (1 << b)) for k in range(8) for b in range(3) if (k ^ (1 << b)) > k]
for sd in range(3):
    d0 = DSTAR*(1 + np.random.default_rng(50+sd).uniform(-0.12, 0.12, 12))
    th0 = bfs_phases(8, CE, d0)
    _, d, pf, pb, t99, _ = evolve(8, CE, d0, 0.5, T=600, th0=th0, noise=0.02)
    gg, liv, churn = stat(CE, d, pf, pb)
    c0 = evolve(8, CE, d0, 0.0, T=600, th0=th0, noise=0.02)
    gg0, _, churn0 = stat(CE, c0[1], c0[2], c0[3])
    print(f"# seed{sd}: plastic gg min={gg.min():.4f} mean={gg.mean():.4f} churn={churn:.4f} "
          f"t99={'%.0f' % t99 if not np.isnan(t99) else '-'} "
          f"dspread {d0.std()/d0.mean()*100:.1f}%->{d.std()/d.mean()*100:.2f}%"
          f" | frozen gg min={gg0.min():.4f} mean={gg0.mean():.4f} churn={churn0:.4f}")

print("\n# ==== (C) cube + live diagonals, SEEDED-LOCK phases: desert protection ====")
DIAG_ALL = [(k, k2) for k in range(8) for k2 in range(k+1, 8)
            if bin(k ^ k2).count('1') == 2]
for sd in range(3):
    rngd = np.random.default_rng(60+sd)
    d0e = DSTAR*(1 + rngd.uniform(-0.12, 0.12, 12))
    d0d = np.sqrt(2)*DSTAR*0.85*(1 + rngd.uniform(-0.15, 0.15, 12))
    keep = d0d < 1.955
    DIAG = [DIAG_ALL[q] for q in range(12) if keep[q]]
    d0d = d0d[keep]
    E = CE + DIAG
    d0 = np.concatenate([d0e, d0d])
    th0 = bfs_phases(8, CE, d0e)          # seedlock over INTENDED edges only
    nliv0 = int((d0d < 1.7).sum())
    for tag, htau in (("plain", 0.0), ("harden tau=20", 20.0)):
        _, d, pf, pb, t99, tlk = evolve(8, E, d0, 0.5, T=1200, th0=th0,
                                        noise=0.02, harden_tau=htau, seed=90+sd)
        gg, liv, churn = stat(E, d, pf, pb)
        eg = gg[:12]; pg_ = gg[12:]; plv = liv[12:]
        pdrift = [f"{d0[12+q]:.2f}->{d[12+q]:.2f}" for q in range(len(DIAG)) if geo(d0[12+q]) > 0]
        print(f"# seed{sd} {tag}: edges gg min={eg.min():.3f} mean={eg.mean():.3f} | "
              f"parasites live {int(plv.sum())}/{len(DIAG)} (seeded live {nliv0}) "
              f"gg_live={pg_[plv].mean() if plv.any() else 0:.3f} churn={churn:.3f}")
        if tag == "plain":
            print(f"#        live-parasite d: {', '.join(pdrift)}")

print("\n# ==== (C2) same but RANDOM phases (natural formation): who wins? ====")
for sd in range(2):
    rngd = np.random.default_rng(60+sd)
    d0e = DSTAR*(1 + rngd.uniform(-0.12, 0.12, 12))
    d0d = np.sqrt(2)*DSTAR*0.85*(1 + rngd.uniform(-0.15, 0.15, 12))
    keep = d0d < 1.955
    DIAG = [DIAG_ALL[q] for q in range(12) if keep[q]]
    E = CE + DIAG
    d0 = np.concatenate([d0e, d0d[keep]])
    for tag, htau in (("plain", 0.0), ("harden tau=20", 20.0)):
        _, d, pf, pb, t99, tlk = evolve(8, E, d0, 0.5, T=1200, th0=None,
                                        noise=0.02, harden_tau=htau, seed=90+sd)
        gg, liv, churn = stat(E, d, pf, pb)
        eg = gg[:12]; pg_ = gg[12:]; plv = liv[12:]
        print(f"# seed{sd} {tag}: edges gg min={eg.min():.3f} mean={eg.mean():.3f} | "
              f"parasites live {int(plv.sum())}/{len(DIAG)} gg_live={pg_[plv].mean() if plv.any() else 0:.3f} "
              f"churn={churn:.3f}")

print("\n# ==== (D) odd rings with noise: equal split or seam concentration? ====")
for N in (3, 5):
    E = [(k, (k+1) % N) for k in range(N)]
    d0 = DSTAR*(1 + np.random.default_rng(70+N).uniform(-0.12, 0.12, N))
    _, d, pf, pb, t99, _ = evolve(N, E, d0, 0.5, T=1200, noise=0.02, seed=80+N)
    gg, liv, churn = stat(E, d, pf, pb)
    print(f"# ring{N}: gg={np.round(np.sort(gg),3)} churn={churn:.4f} "
          f"om*d/pi={np.round(OM0*d/np.pi,3)}")
    print(f"#     equal-split prediction gg={gate(np.pi/(2*N))**2:.3f} each; "
          f"concentration = one dark link, rest locked")
