#!/usr/bin/env python3
"""K2a A/B reference (../K2a_PREREG.md G2, G3): the kernel module lockq
against lo.py on the same open chain (R6 end cells at both ends, mobile
bulk, T-floor t = tJ0 * J0).

  G2: the zero-charge, no-hole ground energy (in J0).
  G3: a seeded jump (photon mode 0) at bulk site s on that ground state,
      then exact evolution to time tL; the hole-density, S^z (and I^z)
      profiles.
  K2b B5: the same with photon mode nu (default 0).

usage: k2a_ab.py M n tJ0 [s tL [nu]]"""
import sys, os
import numpy as np
from scipy.sparse.linalg import eigsh, expm_multiply
here = __file__.rsplit('/', 1)[0]
sys.path.insert(0, here + '/../../perm/cas')
import graphm as GM, lo
lo.P11 = os.environ.get('P11', 'zero')     # K2d: the (1,1) link pricing (default: the original)

M, n, tJ0 = int(sys.argv[1]), int(sys.argv[2]), float(sys.argv[3])
P = GM.P_W if M == 4 else GM.P_2
g = 0.1*GM.DELTA; J0 = 2*g*g/GM.DELTA; t = tJ0*J0
locs = [lo.Local(M, 'end')] + [lo.Local(M, 'mob') for _ in range(n - 2)] + [lo.Local(M, 'end')]
lf = lo.labfun_factory(M, True)
T = lo.chain_terms(locs, P, g, t=t, ring=False)
k = len(locs[0].lab[0])
z0 = tuple([0]*k) + (0,)
H0, b0, dims, strides = lo.sector(locs, T, z0, lf)
w, v = eigsh(H0, k=2, which='SA', tol=1e-13)
o = np.argsort(w); w = w[o]; v = v[:, o]
print(f"# K2a A/B M={M} n={n} t={tJ0} J0: sector {z0} dim {len(b0)}")
print(f"G2 E_gs/J0 = {w[0]/J0:.12f}   (next level {w[1]/J0:.6f})")
if len(sys.argv) > 5:
    s, tL = int(sys.argv[4]), float(sys.argv[5])
    nu = int(sys.argv[6]) if len(sys.argv) > 6 else 0
    g0 = v[:, 0]
    L = locs[s]
    mu0 = [mu for mu in range(M) if P[mu, nu] != 0][0]
    dig = (b0//strides[s]) % dims[s]
    out = {}
    for kk, c in enumerate(b0):
        occ = L.occ[dig[kk]]
        if L.N[dig[kk]] != 2 or occ[mu0] == 0: continue
        o2 = list(occ); o2[mu0] -= 1; kd = L.idx[tuple(o2)]
        c2 = c + (kd - dig[kk])*strides[s]
        out[c2] = out.get(c2, 0.0) + P[mu0, nu]*np.sqrt(occ[mu0])*g0[kk]
    lab_mu = np.array(GM.MODES_W[mu0]) if M == 4 else np.array([1 if mu0 == 0 else -1])
    tgt = tuple(-lab_mu) + (1,)
    H1, b1, _, _ = lo.sector(locs, T, tgt, lf)
    v1 = np.zeros(len(b1), complex); pos = {c: i for i, c in enumerate(b1)}
    for c, x in out.items(): v1[pos[c]] += x
    v1 /= np.linalg.norm(v1)
    vt = expm_multiply(-1j*tL*H1, v1)
    p = np.abs(vt)**2
    nh = []; sz = []; iz = []
    for site in range(n):
        Ls = locs[site]; d = (b1//strides[site]) % dims[site]
        nh.append(float(p@(Ls.N[d] == 1)) if Ls.kind == 'mob' else 0.0)
        sz.append(float(p@(Ls.lab[d, 0]/2.0)))
        if M == 4: iz.append(float(p@(Ls.lab[d, 1]/2.0)))
    E1 = float(np.real(np.vdot(vt, H1@vt)))
    print(f"G3 jump at site {s}, photon mode {nu} -> sector {tuple(int(x) for x in tgt)} dim {len(b1)}; <H>/J0 = {E1/J0:.12f} at tL = {tL}")
    print("G3 nh: " + " ".join(f"{x:.10f}" for x in nh))
    print("G3 sz: " + " ".join(f"{x:+.10f}" for x in sz))
    if M == 4: print("G3 iz: " + " ".join(f"{x:+.10f}" for x in iz))
