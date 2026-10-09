#!/usr/bin/env python3
"""K2d replay: rebuild a kernel lock's history from its event lines.

Reads a freecell output with lockq=1: the LOCKQ header (M, n, t/J0), the
'# LQ event' lines (time, bulk site, photon mode nu) and the last '# LQ t='
diag line.  Replays it exactly in lo.py: the no-hole ground state, closed
evolution between events, the jump sum_mu P[mu,nu] b_(site,mu) at each event.
Checks the final hole, S^z and I^z profiles and <H> against the kernel's
diag line (an A/B of the whole history), then reports the final state:
N_h, adjacency sum_i <n_i n_i+1>, its mean B = <H> - E0 - N_h (E1 - E0), and
its weight on the sector's levels with B_k < -0.05 J0 (bound) — E1 is the
one-hole ground (minimal sector).
usage: k2d_replay.py kernel.out p11"""
import sys, re
import numpy as np
from scipy.sparse.linalg import eigsh, expm_multiply
here = __file__.rsplit('/', 1)[0]
sys.path.insert(0, here + '/../../perm/cas')
import graphm as GM, lo

txt = open(sys.argv[1]).read()
lo.P11 = sys.argv[2] if len(sys.argv) > 2 else 'zero'
PROJ = len(sys.argv) > 3 and sys.argv[3] == 'proj'      # K2e: project each jump (singlet core / core)
m = re.search(r'# LOCKQ v102 K2a: M=(\d+) n=(\d+) g/Delta=(\S+) t/J0=(\S+)', txt)
M, n, tJ0 = int(m.group(1)), int(m.group(2)), float(m.group(4))
# an event fires after its lock step's evolution, at the end of the kernel step in which the lock
# step runs: t + dt, where the line prints that kernel step's start t (lq_step runs after step()).
# For lockq_every = 1 this is t + dt_L, as K2d wrote it; for lockq_every > 1 it is still t + dt.
dtL = float(re.search(r' dt=(\S+)', txt).group(1))
ev = [(float(a) + dtL, int(b), int(c)) for a, b, c in re.findall(r'# LQ event t=(\S+) site=(\d+) cell=\d+ o=\d+ nu=(\d+)', txt)]
# seeded jumps (lockq_force, lockq_force2) are applied at install, t = 0, before any evolution
ev = [(0.0, int(b), int(c)) for b, c in re.findall(r'# LOCKQ seeded jump at site (\d+), photon mode (\d+):', txt)] + ev
diag = re.findall(r'# LQ t=(\S+) tL=(\S+) E_lock=\S+ Ex/J0=(\S+) Nh=(\d+).*?\| nh: ([^|]+)\| sz: ([^|\n]+)(?:\| iz: ([^\n]+))?', txt)
P = GM.P_W if M == 4 else GM.P_2
g = 0.1*GM.DELTA; J0 = 2*g*g/GM.DELTA; t = tJ0*J0
locs = [lo.Local(M, 'end')] + [lo.Local(M, 'mob') for _ in range(n - 2)] + [lo.Local(M, 'end')]
lf = lo.labfun_factory(M, True)
T = lo.chain_terms(locs, P, g, t=t, ring=False)
k = len(locs[0].lab[0])
cache = {}
def sect(tgt):
    if tgt not in cache: cache[tgt] = lo.sector(locs, T, tgt, lf)
    return cache[tgt]
z0 = tuple([0]*k) + (0,)
H0, b0, dims, strides = sect(z0)
w, v = eigsh(H0, k=1, which='SA', tol=1e-13); E0 = w[0]
psi, tgt, basis, tnow = v[:, 0].astype(complex), z0, b0, 0.0
def jump(psi, basis, tgt, site, nu):
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
    H, b, _, _ = sect(nt); pos = {c: i for i, c in enumerate(b)}
    x = np.zeros(len(b), complex)
    for c, a in out.items(): x[pos[c]] += a
    if PROJ: x = project(x, nt)
    return x/np.linalg.norm(x), nt, b
def pi_link(Li, Lj):
    c2 = float((P*P).sum()); X = 0
    if lo.P11 == 'core': return np.kron(Li.P[1], Lj.P[1])
    for r in range(M):
        for m in range(M):
            for s_ in range(M):
                for nn in range(M):
                    w = P[r, s_]*P[m, nn]
                    if w: X = X + w*np.kron(Li.P[1]@Li.E[(r, m)], Lj.P[1]@Lj.E[(s_, nn)])
    return X/c2
pis = {}
def project(x, tg):
    if tg not in pis: pis[tg] = [lo.sector(locs, [((l, l + 1), pi_link(locs[l], locs[l + 1]))], tg, lf)[0] for l in range(n - 1)]
    n2 = np.vdot(x, x).real
    for sw in range(200):
        rem = 0.0
        for Pl in pis[tg]:
            y = Pl@x; rem += np.vdot(y, y).real; x = x - y
        if rem <= 1e-30*n2: break
    return x
nbr = []                                   # hole density on the new hole's neighbours, just before each event
for (te, site, nu) in ev:
    if te > tnow: psi = expm_multiply(-1j*(te - tnow)*sect(tgt)[0], psi); tnow = te
    pp = np.abs(psi)**2
    nbr.append(sum(float(pp@(locs[s].N[(basis//strides[s]) % dims[s]] == 1)) for s in (site - 1, site + 1) if locs[s].kind == 'mob'))
    psi, tgt, basis = jump(psi, basis, tgt, site, nu)
psi_e = psi.copy()
pe = np.abs(psi)**2                        # K3: the profile just after the last event (the walk's start)
nh_e = np.array([float(pe@(locs[s].N[(basis//strides[s]) % dims[s]] == 1)) if locs[s].kind == 'mob' else 0.0 for s in range(n)])
te_last = tnow
tL = float(diag[-1][1]) if diag else tnow
if tL > tnow: psi = expm_multiply(-1j*(tL - tnow)*sect(tgt)[0], psi)
H = sect(tgt)[0]
Ex = float(np.real(np.vdot(psi, H@psi)))
p = np.abs(psi)**2
nh = np.array([float(p@(locs[s].N[(basis//strides[s]) % dims[s]] == 1)) if locs[s].kind == 'mob' else 0.0 for s in range(n)])
sz = np.array([float(p@(locs[s].lab[(basis//strides[s]) % dims[s], 0]/2.0)) for s in range(n)])
Nh = tgt[-1]
print(f"# K2d replay of {sys.argv[1].rsplit('/',1)[-1]}: M={M} n={n} p11={lo.P11}{' proj' if PROJ else ''} events={len(ev)} final sector {tgt} dim {len(basis)} tL={tL}")
if diag:
    d = diag[-1]
    knh = np.array([float(x) for x in d[4].split()]); ksz = np.array([float(x) for x in d[5].split()])
    print(f"  A/B against the kernel's last diag: |<H>/J0 diff| {abs(Ex/J0 - float(d[2])):.1e}, max|nh diff| {np.abs(nh - knh).max():.1e}, max|sz diff| {np.abs(sz - ksz).max():.1e}")
if Nh >= 1:
    q1 = (1, 1, 1) if M == 4 else (1, 1)
    E1 = eigsh(sect(q1)[0], k=1, which='SA', tol=1e-12)[0][0]
    adjp = 0.0
    for s in range(1, n - 2):
        ds = (basis//strides[s]) % dims[s]; dt = (basis//strides[s + 1]) % dims[s + 1]
        adjp += float(p@((locs[s].N[ds] == 1) & (locs[s + 1].N[dt] == 1)))
    Bm = (Ex - E0 - Nh*(E1 - E0))/J0
    print(f"  events (site, nu, hole density on its neighbours): " + " ".join(f"({si},{nu},{x:.2f})" for (te, si, nu), x in zip(ev, nbr))
          + f"; next to an end: {sum(1 for (te, si, nu) in ev if si in (1, n - 2))}")
    line = f"  final: N_h={Nh}, nh = " + " ".join(f"{x:.3f}" for x in nh) + f"; adjacency {adjp:.3f}; mean B {Bm:+.3f} J0"
    if Nh >= 2:
        if H.shape[0] <= 4000:
            wk, Uk = np.linalg.eigh(H.toarray())
        else:
            wk, Uk = eigsh(H, k=60, which='SA', tol=1e-11)
        Bk = (wk - E0 - Nh*(E1 - E0))/J0
        ov = np.abs(Uk.T.conj()@psi)**2
        line += f"; bound weight (B_k < -0.05) {ov[Bk < -0.05].sum():.4f} (levels kept {ov.sum():.3f}, lowest B_k {Bk.min():+.3f})"
        if Nh == 2:                        # K3: hole-hole distance of the bound and the unbound parts
            occ = np.array([[locs[s].kind == 'mob' and locs[s].N[(c//strides[s]) % dims[s]] == 1 for s in range(n)] for c in basis])
            dd = np.array([np.ptp(np.nonzero(r)[0]) if r.sum() == 2 else 0 for r in occ], float)
            Ub = Uk[:, Bk < -0.05]
            def dparts(x):
                xb = Ub@(Ub.T.conj()@x); xu = x - xb
                wb, wu = np.vdot(xb, xb).real, np.vdot(xu, xu).real
                return (float(np.abs(xb)**2@dd/wb) if wb > 1e-12 else float('nan'), wb, float(np.abs(xu)**2@dd/wu) if wu > 1e-12 else float('nan'), wu)
            a, c = dparts(psi_e), dparts(psi)
            pair_line = (f"  pair: <d> of the bound part {a[0]:.3f} -> {c[0]:.3f} (weight {a[1]:.4f} -> {c[1]:.4f}); "
                         f"<d> of the unbound part {a[2]:.3f} -> {c[2]:.3f} (after the last event -> tL)")
    if PROJ:
        sw = sum(np.vdot(psi, Pl@psi).real for Pl in pis.get(tgt, [])) if tgt in pis else float('nan')
        line += f"; (1,1) singlet weight {sw:.2e}"
    print(line)
    if 'pair_line' in globals(): print(pair_line)
    print(f"  after the last event (t={te_last:.4f}): nh = " + " ".join(f"{x:.3f}" for x in nh_e) + f"; profile change to tL, half sum |nh(tL) - nh|: {0.5*np.abs(nh - nh_e).sum():.4f}")
