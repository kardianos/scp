#!/usr/bin/env python3
"""K2e E2 reference: seeded jumps under singlet core (lo.P11 = 'score'), each
projected onto the allowed subspace (no conversion singlet on any (1,1) link)
by the filter exp(-tau H_c), H_c = sum_l Pi_l, tau = 1000 — independent of the
kernel's alternating projections.  Then exact evolution to tL.
usage: k2e_ab.py M n tJ0 tL site1 mode1 [site2 mode2]"""
import sys
import numpy as np
from scipy.sparse.linalg import eigsh, expm_multiply
here = __file__.rsplit('/', 1)[0]
sys.path.insert(0, here + '/../../perm/cas')
import graphm as GM, lo
lo.P11 = 'score'
M, n, tJ0, tL = int(sys.argv[1]), int(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4])
jumps = [(int(sys.argv[k]), int(sys.argv[k + 1])) for k in range(5, len(sys.argv), 2)]
P = GM.P_W if M == 4 else GM.P_2
g = 0.1*GM.DELTA; J0 = 2*g*g/GM.DELTA; t = tJ0*J0
locs = [lo.Local(M, 'end')] + [lo.Local(M, 'mob') for _ in range(n - 2)] + [lo.Local(M, 'end')]
lf = lo.labfun_factory(M, True)
T = lo.chain_terms(locs, P, g, t=t, ring=False)
c2 = float((P*P).sum())
def pi_link(Li, Lj):
    X = 0
    for r in range(M):
        for m in range(M):
            for s_ in range(M):
                for nn in range(M):
                    w = P[r, s_]*P[m, nn]
                    if w: X = X + w*np.kron(Li.P[1]@Li.E[(r, m)], Lj.P[1]@Lj.E[(s_, nn)])
    return X/c2
Tc = [((l, l + 1), pi_link(locs[l], locs[l + 1])) for l in range(n - 1)]
k = len(locs[0].lab[0])
z0 = tuple([0]*k) + (0,)
H0, b0, dims, strides = lo.sector(locs, T, z0, lf)
w, v = eigsh(H0, k=1, which='SA', tol=1e-13)
psi, tgt, basis = v[:, 0].astype(complex), z0, b0
for (site, nu) in jumps:
    L = locs[site]; dig = (basis//strides[site]) % dims[site]; out = {}
    for kk, c in enumerate(basis):
        if L.N[dig[kk]] != 2: continue
        occ = L.occ[dig[kk]]
        for mu in range(M):
            p = P[mu, nu]
            if p == 0 or occ[mu] == 0: continue
            o2 = list(occ); o2[mu] -= 1; kd = L.idx[tuple(o2)]
            c2_ = c + (kd - dig[kk])*strides[site]
            out[c2_] = out.get(c2_, 0.0) + p*np.sqrt(occ[mu])*psi[kk]
    mu = [mu for mu in range(M) if P[mu, nu] != 0][0]
    q = np.array(GM.MODES_W[mu]) if M == 4 else np.array([1 if mu == 0 else -1])
    tgt = tuple(int(x) for x in np.array(tgt[:-1]) - q) + (int(tgt[-1]) + 1,)
    H, basis, _, _ = lo.sector(locs, T, tgt, lf)
    x = np.zeros(len(basis), complex); pos = {c: i for i, c in enumerate(basis)}
    for c, a in out.items(): x[pos[c]] += a
    Hc, _, _, _ = lo.sector(locs, Tc, tgt, lf)
    xp = expm_multiply(-1000.0*Hc, x)
    jj = np.vdot(xp, xp).real/np.vdot(psi, psi).real
    resid = np.linalg.norm(Hc@xp)/np.linalg.norm(xp)
    print(f"jump at site {site}, mode {nu}: sector {tgt} dim {len(basis)}; projected weight {jj:.15f}; |H_c x|/|x| after the filter {resid:.1e}")
    psi = xp/np.linalg.norm(xp)
H = lo.sector(locs, T, tgt, lf)[0]
vt = expm_multiply(-1j*tL*H, psi)
p = np.abs(vt)**2
nh = [float(p@(locs[s].N[(basis//strides[s]) % dims[s]] == 1)) if locs[s].kind == 'mob' else 0.0 for s in range(n)]
sz = [float(p@(locs[s].lab[(basis//strides[s]) % dims[s], 0]/2.0)) for s in range(n)]
iz = [float(p@(locs[s].lab[(basis//strides[s]) % dims[s], 1]/2.0)) for s in range(n)] if M == 4 else []
Hc = lo.sector(locs, Tc, tgt, lf)[0]
print(f"<H>/J0 = {np.vdot(vt, H@vt).real/J0:.12f} at tL = {tL}; singlet weight on (1,1) links {np.vdot(vt, Hc@vt).real:.3e}")
print("nh: " + " ".join(f"{x:.10f}" for x in nh))
print("sz: " + " ".join(f"{x:+.10f}" for x in sz))
if M == 4: print("iz: " + " ".join(f"{x:+.10f}" for x in iz))
