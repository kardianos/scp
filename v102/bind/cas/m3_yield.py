#!/usr/bin/env python3
"""M3a (../PREREG.md) — P1: a pair born bound.  The binding is parked in the
space quantum of the formation event that completes the pair.

Leading-order ring (lo.py), mobile holes (T-floor, t = 0.1 J0), (1,1) pricing
p11.  First jump (photon mode 0) at cell 0 on the no-hole ground state; second
jump (the complementary mode: 1 for 3+2, 3 for W, i.e. opposite charges) at
cell d.  Reported per d: the jump state's weight on the bound two-hole levels
(B_k = E_k + E0 - 2 E1 < -0.05 J0; their mean pair distance), the mean B of
the jump state, and the exact energy change of the second event as a
detuning delta of its space quantum (window |delta| < 2 s_k = 0.12).
usage: m3_yield.py M n p11"""
import sys, time, functools
print = functools.partial(print, flush=True)
import numpy as np
from scipy.sparse.linalg import eigsh
here = __file__.rsplit('/', 1)[0]
sys.path.insert(0, here + '/../../perm/cas')
import graphm as GM, lo

t0 = time.time()
M, n, p11 = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3]
lo.P11 = p11
P = GM.P_W if M == 4 else GM.P_2
g = 0.1*GM.DELTA; J0 = 2*g*g/GM.DELTA; t = 0.1*J0
HB = GM.EF[2] - GM.EF[1]            # E_f step of a hole
W1E, OM, SK = 1.65, GM.OMEGA, 0.06
mob = [lo.Local(M, 'mob') for _ in range(n)]
lf = lo.labfun_factory(M, True)
T = lo.chain_terms(mob, P, g, t=t, ring=True)
k = len(mob[0].lab[0])
def sect(tgt):
    H, b, dims, strides = lo.sector(mob, T, tgt, lf)
    return H, b, dims, strides
H0, b0, dims, strides = sect(tuple([0]*k) + (0,))
w0, v0 = eigsh(H0, k=1, which='SA', tol=1e-12); E0 = w0[0]; g0 = v0[:, 0]
nu_a, nu_b = 0, M - 1
qa = lo.Local(M, 'hole').lab[nu_a]; qb = lo.Local(M, 'hole').lab[nu_b]
def jump(vec, basis, site, nu, tb):
    L = mob[site]; out = np.zeros(len(tb), complex); pos = {c: i for i, c in enumerate(tb)}
    dig = (basis//strides[site]) % dims[site]
    for kk, c in enumerate(basis):
        if abs(vec[kk]) < 1e-15 or L.N[dig[kk]] != 2: continue
        occ = L.occ[dig[kk]]
        for mu in range(M):
            p = P[mu, nu]
            if p == 0 or occ[mu] == 0: continue
            o2 = list(occ); o2[mu] -= 1; kd = L.idx[tuple(o2)]
            out[pos[c + (kd - dig[kk])*strides[site]]] += p*np.sqrt(occ[mu])*vec[kk]
    return out
H1, b1, _, _ = sect(tuple(qa) + (1,))
w1, _ = eigsh(H1, k=1, which='SA', tol=1e-12); E1 = w1[0]
v1 = jump(g0, b0, 0, nu_a, b1); v1 /= np.linalg.norm(v1)
H2, b2, _, _ = sect(tuple(np.array(qa) + np.array(qb)) + (2,))
K = min(80, H2.shape[0] - 2)
w2, U2 = eigsh(H2, k=K, which='SA', tol=1e-11); o = np.argsort(w2); w2 = w2[o]; U2 = U2[:, o]
B2 = (w2 + E0 - 2*E1)/J0
hp = np.array([[mob[i].N[(c//strides[i]) % dims[i]] == 1 for i in range(n)] for c in b2])
def pd(h):
    a, b = np.where(h)[0]; d = b - a; return min(d, n - d)
dist = np.array([pd(h) for h in hp])
bound = np.where(B2 < -0.05)[0]
mdist = [float(np.abs(U2[:, j])**2 @ dist) for j in bound]
e1 = float(np.real(np.vdot(v1, H1@v1)))
print(f"# M3a P1 yield: {'W' if M == 4 else '3+2'} ring n={n}, p11={p11}, t=0.1 J0; E0={E0/J0:.6f} E1={E1/J0:.6f} J0 (two-hole dim {len(b2)})")
print(f"  bound levels (B < -0.05 J0): {len(bound)}; B = " + " ".join(f"{B2[j]:+.3f}" for j in bound[:12]) + "; mean distance " + " ".join(f"{x:.2f}" for x in mdist[:12]))
print(f"  lowest unbound B: {B2[len(bound)] if len(bound) < K else float('nan'):+.4f}")
for d in range(1, n//2 + 1):
    v2 = jump(v1, b1, d, nu_b, b2); nrm = np.linalg.norm(v2)
    if nrm < 1e-12: print(f"  d={d}: no jump"); continue
    v2 /= nrm
    ov = np.abs(U2.T.conj()@v2)**2
    e2 = float(np.real(np.vdot(v2, H2@v2)))
    delta = W1E - OM - ((-HB) + (e2 - e1))
    print(f"  d={d}: bound weight {ov[bound].sum():.4f}  (kept levels {ov.sum():.3f});  mean B of the jump state {(e2 + E0 - 2*E1)/J0:+.3f} J0;  delta of the event {delta:+.4f} ({'in' if abs(delta) < 2*SK else 'OUT of'} window)")
print(f"done [{time.time()-t0:.0f}s]")
