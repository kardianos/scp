#!/usr/bin/env python3
"""K1-2 — two holes made by two jumps in an isolated lock (../K1_PREREG.md).
Leading-order lock, mobile holes (T-floor, t = 0.1 J0).  The lock's ground
state (no hole), then two jumps J_nu (photon modes of opposite charges) at
cells a and b; evolve in the closed lock to T = 200/J0.  Report the
weight in the bound pair (the two-hole ground state), at t = 0 and averaged
over time, and the energy above the pair ground state.
usage: k1_closed.py M n a b"""
import sys, time
import numpy as np
from scipy.sparse.linalg import eigsh, expm_multiply
here = __file__.rsplit('/', 1)[0]
sys.path.insert(0, here + '/../../perm/cas')
import graphm as GM, lo

t0 = time.time()
M, n, a, b = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])
P = GM.P_W if M == 4 else GM.P_2
g = 0.1*GM.DELTA; J0 = 2*g*g/GM.DELTA; t = 0.1*J0
mob = [lo.Local(M, 'mob') for _ in range(n)]
lf = lo.labfun_factory(M, True)
T = lo.chain_terms(mob, P, g, t=t, ring=True)
k = len(mob[0].lab[0])
z0 = tuple([0]*k) + (0,)
H0, b0, dims, strides = lo.sector(mob, T, z0, lf)
w0, v0 = eigsh(H0, k=1, which='SA', tol=1e-12); g0 = v0[:, 0]
modes = list(range(M))
nu_a, nu_b = 0, M - 1                      # opposite charges: (+,+) and (-,-) for W; + and - for 3+2
qa = lo.Local(M, 'hole').lab[nu_a]; qb = lo.Local(M, 'hole').lab[nu_b]
def jump(vec, basis, site, nu):
    """apply sum_mu P[mu,nu] b_{site,mu} to vec (on basis); returns {code: amp}"""
    L = mob[site]; out = {}
    dig = (basis//strides[site]) % dims[site]
    for kk, c in enumerate(basis):
        if abs(vec[kk]) < 1e-15: continue
        occ = L.occ[dig[kk]]
        if L.N[dig[kk]] != 2: continue
        for mu in modes:
            p = P[mu, nu]
            if p == 0 or occ[mu] == 0: continue
            o2 = list(occ); o2[mu] -= 1; kd = L.idx[tuple(o2)]
            c2 = c + (kd - dig[kk])*strides[site]
            out[c2] = out.get(c2, 0.0) + p*np.sqrt(occ[mu])*vec[kk]
    return out
st1 = jump(g0, b0, a, nu_a)
tgt1 = tuple(np.array(qa)) + (1,)
H1, b1, _, _ = lo.sector(mob, T, tgt1, lf)
v1 = np.zeros(len(b1)); pos = {c: i for i, c in enumerate(b1)}
for c, x in st1.items(): v1[pos[c]] += x
v1 /= np.linalg.norm(v1)
st2 = jump(v1, b1, b, nu_b)
tgt2 = tuple(np.array(qa) + np.array(qb)) + (2,)
H2, b2, _, _ = lo.sector(mob, T, tgt2, lf)
v2 = np.zeros(len(b2), complex); pos = {c: i for i, c in enumerate(b2)}
for c, x in st2.items(): v2[pos[c]] += x
v2 /= np.linalg.norm(v2)
w2, U2 = eigsh(H2, k=8, which='SA', tol=1e-12); o = np.argsort(w2); w2 = w2[o]; U2 = U2[:, o]
ov = np.abs(U2.T.conj()@v2)**2
E = float(np.real(np.vdot(v2, H2@v2)))
print(f"# K1-2 M={M} n={n} jumps at {a},{b} (modes {nu_a},{nu_b}): two-hole sector {tgt2} dim {len(b2)}")
print(f"  two-hole levels (E - E0 in J0): " + " ".join(f"{(x-w2[0])/J0:.3f}" for x in w2))
print(f"  initial weights on those levels: " + " ".join(f"{x:.4f}" for x in ov))
print(f"  energy of the jump state above the bound pair: {(E - w2[0])/J0:.3f} J0")
# closed evolution: the bound weight in time
ts = np.linspace(0, 200/J0, 21)
traj = expm_multiply(-1j*H2, v2, start=0, stop=ts[-1], num=len(ts), endpoint=True)
bw = np.abs(traj@U2[:, 0].conj())**2
print(f"  bound-pair weight in time: t=0 {bw[0]:.4f}, mean {bw.mean():.4f}, max {bw.max():.4f}, norm drift {abs(np.linalg.norm(traj[-1])-1):.1e}")
print(f"  K1-2 gate (time-mean <= t0 + 0.05): {'PASS' if bw.mean() <= bw[0] + 0.05 else 'FAIL'}  [{time.time()-t0:.0f}s]")
