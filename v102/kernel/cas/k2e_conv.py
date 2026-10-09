#!/usr/bin/env python3
"""K2e: how well do 200 sweeps of alternating projections reach the allowed subspace?
Replays a kernel lock history with lockq_p11=2 lockq_proj=1 on the kernel's own path (each
fired jump projected with the 200-sweep rule, as in the kernel and k2d_replay.py), and before
each event projects EVERY candidate jump (each bulk site without a hole, each photon mode) twice:
with the 200-sweep rule; where that rule stops at its cap (201 = not converged), the projection is
continued to deep convergence (removed weight per sweep <= 1e-34 of the norm, at most 100000 sweeps).
Every result of alternating projections is x minus a sum of vectors in the ranges of the Pi_l, so it
is the orthogonal projection onto the allowed subspace (ker H_c, H_c = sum_l Pi_l) once its residual
singlet weight sum_l <Pi_l> is zero; the residual and <x200|x - x200> are reported as the checks.
(A MINRES solve of H_c z = H_c x was tried first as the reference and rejected: in floating point its
z picked up a kernel component, <xd|x - xd> = 5e-4, while x200 was orthogonal to 2e-17.) Reports per event: the sweeps the 200 rule used (201 = not converged), the
largest relative change of a candidate's projected weight |w200 - wdeep|/wdeep (the rate), the
largest residual singlet weight sum_l <Pi_l> / |x|^2 left by the 200 rule, and for the fired jump
the distance between the two normalized projected states.
usage: k2e_conv.py kernel.out [dt_sample]  (dt_sample > 0: also every candidate at t = dt_sample, 2 dt_sample, ... tL)"""
import sys, re
import numpy as np
from scipy.sparse.linalg import eigsh, expm_multiply
here = __file__.rsplit('/', 1)[0]
sys.path.insert(0, here + '/../../perm/cas')
import graphm as GM, lo
txt = open(sys.argv[1]).read()
lo.P11 = 'score'
m = re.search(r'# LOCKQ v102 K2a: M=(\d+) n=(\d+) g/Delta=(\S+) t/J0=(\S+)', txt)
M, n, tJ0 = int(m.group(1)), int(m.group(2)), float(m.group(4))
dtL = float(re.search(r' dt=(\S+)', txt).group(1))*int(re.search(r'# LOCKQ v102 K2a: .* every=(\d+)', txt).group(1))
ev = [(float(a) + dtL, int(b), int(c)) for a, b, c in re.findall(r'# LQ event t=(\S+) site=(\d+) cell=\d+ o=\d+ nu=(\d+)', txt)]
P = GM.P_W if M == 4 else GM.P_2
g = 0.1*GM.DELTA; J0 = 2*g*g/GM.DELTA; t = tJ0*J0
locs = [lo.Local(M, 'end')] + [lo.Local(M, 'mob') for _ in range(n - 2)] + [lo.Local(M, 'end')]
lf = lo.labfun_factory(M, True)
T = lo.chain_terms(locs, P, g, t=t, ring=False)
cache = {}
def sect(tgt):
    if tgt not in cache: cache[tgt] = lo.sector(locs, T, tgt, lf)
    return cache[tgt]
k = len(locs[0].lab[0])
z0 = tuple([0]*k) + (0,)
H0, b0, dims, strides = sect(z0)
w, v = eigsh(H0, k=1, which='SA', tol=1e-13)
psi, tgt, basis, tnow = v[:, 0].astype(complex), z0, b0, 0.0
def rawjump(psi, basis, tgt, site, nu):
    L = locs[site]; dig = (basis//strides[site]) % dims[site]; out = {}
    for kk, c in enumerate(basis):
        if abs(psi[kk]) < 1e-300 or L.N[dig[kk]] != 2: continue
        occ = L.occ[dig[kk]]
        for mu in range(M):
            p = P[mu, nu]
            if p == 0 or occ[mu] == 0: continue
            o2 = list(occ); o2[mu] -= 1; kd = L.idx[tuple(o2)]
            c2 = c + (kd - dig[kk])*strides[site]
            out[c2] = out.get(c2, 0.0) + p*np.sqrt(occ[mu])*psi[kk]
    mu = [mu for mu in range(M) if P[mu, nu] != 0][0]
    q = np.array(GM.MODES_W[mu]) if M == 4 else np.array([1 if mu == 0 else -1])
    nt = tuple(int(x) for x in np.array(tgt[:-1]) - q) + (int(tgt[-1]) + 1,)
    _, b, _, _ = sect(nt); pos = {c: i for i, c in enumerate(b)}
    x = np.zeros(len(b), complex)
    for c, a in out.items(): x[pos[c]] += a
    return x, nt, b
def pi_link(Li, Lj):
    c2 = float((P*P).sum()); X = 0
    for r in range(M):
        for m_ in range(M):
            for s_ in range(M):
                for nn in range(M):
                    wgt = P[r, s_]*P[m_, nn]
                    if wgt: X = X + wgt*np.kron(Li.P[1]@Li.E[(r, m_)], Lj.P[1]@Lj.E[(s_, nn)])
    return X/c2
pis = {}
def project(x, tg, cap, tol):
    if tg not in pis: pis[tg] = [lo.sector(locs, [((l, l + 1), pi_link(locs[l], locs[l + 1]))], tg, lf)[0] for l in range(n - 1)]
    n2 = np.vdot(x, x).real
    for sw in range(1, cap + 1):
        rem = 0.0
        for Pl in pis[tg]:
            y = Pl@x; rem += np.vdot(y, y).real; x = x - y
        if rem <= tol*n2: return x, sw
    return x, cap + 1                      # as the kernel counts: cap + 1 = stopped at the cap
def resid(x, tg): return sum(np.vdot(x, Pl@x).real for Pl in pis[tg])/np.vdot(x, x).real
DTS = float(sys.argv[2]) if len(sys.argv) > 2 else 0.0      # optional: also sample candidates every DTS
diag = re.findall(r'# LQ t=(\S+) tL=(\S+) ', txt)
tL = float(diag[-1][1]) if diag else (ev[-1][0] if ev else 0.0)
def cands(psi, basis, tgt, fire=None):
    st = dict(sw=0, dw=0.0, res=0.0, orth=0.0, unc=0); fired = None
    for s in range(1, n - 1):
        for mu in range(M):
            x, nt, b = rawjump(psi, basis, tgt, s, mu)
            if np.vdot(x, x).real < 1e-300: continue
            x200, sw = project(x.copy(), nt, 200, 1e-30)
            xd, swd = project(x200.copy(), nt, 100000, 1e-34) if sw > 200 else (x200, 0)
            st['orth'] = max(st['orth'], abs(np.vdot(x200, x - x200))/np.vdot(x, x).real)
            w200, wd = np.vdot(x200, x200).real, np.vdot(xd, xd).real
            st['unc'] += sw > 200
            st['sw'] = max(st['sw'], sw); st['res'] = max(st['res'], resid(x200, nt))
            if wd > 1e-300: st['dw'] = max(st['dw'], abs(w200 - wd)/wd)
            if (s, mu) == fire: fired = (x200, xd, nt, b, swd)
    return st, fired
def fmt(st): return (f"candidates' sweeps (200 rule) max {st['sw']}, max|w200-wdeep|/wdeep {st['dw']:.1e}, max residual singlet weight {st['res']:.1e}, "
                     f"max|<x200|x-x200>|/|x|^2 {st['orth']:.1e}, unconverged candidates {st['unc']}")
print(f"# K2e convergence of {sys.argv[1].rsplit('/',1)[-1]}: M={M} n={n} events={len(ev)} tL={tL} sampling every {DTS}", flush=True)
timeline = [(te, 1, site, nu) for (te, site, nu) in ev]
if DTS > 0: timeline += [(k*DTS, 0, None, None) for k in range(1, int(tL/DTS) + 1)]
timeline.sort(key=lambda r: (r[0], r[1]))
worst = dict(sw=0, dw=0.0, res=0.0, dx=0.0, orth=0.0, unc=0)
for (te, isev, site, nu) in timeline:
    if te > tnow: psi = expm_multiply(-1j*(te - tnow)*sect(tgt)[0], psi); tnow = te
    if tgt[-1] >= n - 2: continue
    st, fired = cands(psi, basis, tgt, (site, nu) if isev else None)
    for k in ('sw', 'dw', 'res', 'orth'): worst[k] = max(worst[k], st[k])
    worst['unc'] += st['unc']
    if isev:
        x200, xd, nt, b, swd = fired
        dx = np.linalg.norm(x200/np.linalg.norm(x200) - xd/np.linalg.norm(xd)); worst['dx'] = max(worst['dx'], dx)
        print(f"  event site={site} nu={nu} t={te:.4f}: {fmt(st)}; fired jump: further sweeps {swd}, |x200 - xdeep| {dx:.1e}", flush=True)
        psi, tgt, basis = x200/np.linalg.norm(x200), nt, b
    else:
        print(f"  sample t={te:.4f} (N_h={tgt[-1]}): {fmt(st)}", flush=True)
print(f"  worst: sweeps {worst['sw']}, rate weight change {worst['dw']:.1e}, residual singlet weight {worst['res']:.1e}, max|<x200|x-x200>|/|x|^2 {worst['orth']:.1e}, unconverged candidates {worst['unc']}, fired-state distance {worst['dx']:.1e}")
