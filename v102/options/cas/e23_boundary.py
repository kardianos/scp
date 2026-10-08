#!/usr/bin/env python3
"""E2 / E3 — the lock's boundary (../PREREG.md §3 E2, E3; combination K-a).

Geometry: ring cells 0..3 (n = 4); outside cell o = 4 linked to ring cell 0
(the boundary link); a lead 5..4+LEAD hangs from o.  Links: 0..3 ring,
4 = (0, o), then the lead chain.  Lead links carry no conversion (at most two
outside quanta; documented).  Transport amplitude 0.9 on every link, except
T0, which shuts the ring and boundary links.

usage: e23_boundary.py E2 RULE CEIL            (empty outside: leak)
       e23_boundary.py E3 RULE CEIL X ENERGY   (one photon arrives; X = X0|X1,
                                                ENERGY = cold|res)
"""
import sys, time
import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import expm_multiply, eigsh
sys.path.insert(0, __file__.rsplit('/', 1)[0])
import graph as G

t0 = time.time()
mode, rule, ceil = sys.argv[1], sys.argv[2], sys.argv[3]
X = sys.argv[4] if mode == 'E3' else 'X0'
energy = sys.argv[5] if mode == 'E3' else None
n = 4; g = 0.1*G.DELTA; TOUT = 0.9
LEAD = 40 if mode == 'E2' else 70
o = n; first = n + 1; ncell = n + 1 + LEAD
links = G.ring_links(n) + [(0, o), (o, first)] + [(first + k, first + k + 1) for k in range(LEAD - 1)]
L = len(links)
conv = [True]*n + [X == 'X1'] + [False]*(L - n - 1)
if rule == 'T0':
    t = [0.0]*(n + 1) + [TOUT]*(L - n - 1); trule = 'T2'
else:
    t = [TOUT]*L; trule = rule
if ceil == 'lock':
    gkw = dict(ceil='fixed', caps=[2]*n + [None]*(ncell - n))
elif ceil == 'geo':
    gkw = dict(ceil='geo')
else:
    gkw = dict(ceil='held')
# absorbing potential on the last 15 lead sites (E2 only)
onsite = None
if mode == 'E2':
    k0 = ncell - 15
    eta = {k: 1.0*((k - k0 + 1)/15)**2 for k in range(k0, ncell)}
    onsite = lambda N, k: -1j*eta[k]*N if k in eta else 0.0
outmax = 1 if mode == 'E2' else 2
gr = G.Graph(ncell, links, g, conv=conv, t=t, rule=trule, onsite=onsite, cmax=None,
             core=set(range(n)) | {o}, outmax=outmax, **gkw)

# ---- the R7 lock: isolated ring ground state (T0, local), S^z = 0
ring = G.Graph(n, G.ring_links(n), g, ceil='fixed', caps=[2]*n)
Hr = ring.build(G.ring_seeds(n, n, sz=0), keep=lambda s: G.sz_of(s) == 0).toarray() - n*G.EF[2]*np.eye(len(ring.states))
wr, vr = np.linalg.eigh(Hr); gs = vr[:, 0]; E0 = wr[0]

def embed(rs, photon=None):
    cells = list(rs[0]) + [(0, 0)]*(ncell - n)
    if photon is not None:
        x, pol = photon; cells[x] = (1, 0) if pol == 0 else (0, 1)
    return (tuple(cells), tuple(rs[1]) + (0,)*(L - n))

seeds = [embed(s) for s in ring.states]
if mode == 'E3':
    seeds = [embed(s, (x, 0)) for s in ring.states for x in range(first, ncell)]
H = gr.build(seeds)
H = H - n*G.EF[2]*sp.identity(H.shape[0])
D = H.shape[0]
st = gr.states
N, C = G.obs_arrays(st, ncell, L)
Nring = N[:, :n].sum(1); Nout = N[:, first:].sum(1); No = N[:, o]
crg = C[:, :n].sum(1); cb = C[:, n]
herm = abs(H - H.conj().T).max() if onsite is None else float('nan')
reach_out = int((N[:, n:].sum(1) > (1 if mode == 'E3' else 0)).sum())
tag = f"{mode} rule={rule} ceil={ceil}" + (f" {X} {energy}" if mode == 'E3' else '')
print(f"# {tag}: dim={D}  states with ring field outside the ring: {reach_out}  hermitian(no CAP)={herm}  [{time.time()-t0:.0f}s]", flush=True)

v0 = np.zeros(D, complex)
if mode == 'E2':
    for k, s in enumerate(ring.states):
        kk = gr.idx.get(embed(s))
        if kk is not None: v0[kk] = gs[k]
    print(f"  R7 ground-state weight on states this ceiling forbids: {1 - np.linalg.norm(v0)**2:.3e} (dropped, renormalized)")
    T = 300.0; NT = 61
else:
    # Gaussian packet on the lead moving toward o (decreasing site index).
    # The engine's hop is +t, so eps(k) = +2t cos k and v = -2t sin k: a packet
    # exp(+i k x) with k in (0, pi) moves toward o.
    if energy == 'cold':
        k0 = np.pi - 0.30                            # eps = -2t cos 0.3, near the band bottom
    else:
        eres = G.OMEGA - G.EF[2]                     # resonant kinetic energy
        k0 = np.arccos(eres/(2*TOUT))
    xc = first + 30; sig = 5.0
    xs = np.arange(first, ncell)
    phi = np.exp(-(xs - xc)**2/(4*sig**2))*np.exp(1j*k0*(xs - xc))
    phi /= np.linalg.norm(phi)
    for k, s in enumerate(ring.states):
        for x, a in zip(xs, phi):
            kk = gr.idx.get(embed(s, (x, 0)))
            if kk is not None: v0[kk] += gs[k]*a
    print(f"  R7 ground-state weight on states this ceiling forbids: {1 - np.linalg.norm(v0)**2:.3e} (dropped, renormalized)")
    vel = 2*TOUT*np.sin(k0)
    T = (2*30 + 10)/vel; NT = 41
    print(f"  photon: k0={k0:.4f}  kinetic eps={2*TOUT*np.cos(k0):+.4f}  speed toward o={vel:.3f}  T={T:.1f}")
v0 /= np.linalg.norm(v0)
traj = expm_multiply(-1j*H, v0, start=0, stop=T, num=NT, endpoint=True)
ts = np.linspace(0, T, NT)
P = np.abs(traj)**2
norm = P.sum(1)
if mode == 'E2':
    ring_field = (P@(Nring + 2*crg + cb))/np.maximum(norm, 1e-300)
    for k in range(0, NT, 10):
        print(f"  t={ts[k]:6.1f}: survival(no quantum absorbed)={norm[k]:.6f}  P(quantum outside ring+o)={P[k]@(Nout>0):.6f}  "
              f"P(o occupied)={P[k]@(No>0):.6f}  ring field (normalized, /8)={ring_field[k]/8:.6f}")
    lost = 1 - norm[-1] + P[-1]@(Nout > 0)
    print(f"  RESULT {tag}: escaped by T={T:.0f} (absorbed + in lead) = {lost:.6f}  retained = {1-lost:.6f}  "
          f"first-escape rate ~ {-np.log(max(norm[-1],1e-300))/T:.3e}/time  [{time.time()-t0:.0f}s]")
else:
    pf = P[-1]
    tot_core = N[:, :o + 1].sum(1) + 2*C[:, :n].sum(1) + 2*cb      # quanta in ring + o (field + converted)
    refl = pf[(tot_core == 8) & (Nout == 1)].sum()
    capt = pf[(tot_core == 9)].sum()
    loss = pf[(tot_core <= 7)].sum()
    # ring fidelity in the reflected part: ring cells + ring links, o empty, boundary c = 0
    rid = {s: k for k, s in enumerate(ring.states)}
    amp = {}
    vf = traj[-1]
    for k, s in enumerate(st):
        if abs(vf[k]) < 1e-14: continue
        cells, cs = s
        if cells[o] != (0, 0) or cs[n] != 0 or Nout[k] != 1: continue
        rk = rid.get((cells[:n], cs[:n]))
        if rk is None: continue
        key = cells[n + 1:]
        amp.setdefault(key, np.zeros(len(ring.states), complex))[rk] += vf[k]
    fid = sum(abs(gs@a)**2 for a in amp.values())
    eng = sum(np.real(a.conj()@(Hr@a)) for a in amp.values())
    wref = sum(np.real(a.conj()@a) for a in amp.values())
    print(f"  RESULT {tag}: norm={norm[-1]:.6f}  reflected(ring intact)={refl:.6f}  captured(9 quanta in ring+o)={capt:.6f}  "
          f"ring field lost={loss:.6f}  ring ground fidelity in reflected part={fid/max(wref,1e-300):.6f}  "
          f"ring energy above E0 in reflected part={eng/max(wref,1e-300)-E0:+.3e}  [{time.time()-t0:.0f}s]")
    for k in range(0, NT, 10):
        print(f"    t={ts[k]:6.1f}: P(o occ)={P[k]@(No>0):.5f}  <c_boundary>={P[k]@cb:.5f}  P(9 in core)={P[k][(tot_core==9)].sum():.5f}  "
              f"<N ring>={P[k]@Nring:.5f}")
