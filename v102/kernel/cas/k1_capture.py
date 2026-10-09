#!/usr/bin/env python3
"""K1-3 — capture of two holes into their bound pair by cold light (../K1_PREREG.md).

3+2, n = 6 ring (the smallest ring that holds the separation-3 pair),
leading-order lock with mobile holes (T-floor, t = 0.1 J0, floor exact),
two holes.  One boundary link from cell 0 to the end site of a photon lead
(t_out = 0.9).  Exact multichannel S-matrix as ../../options/cas/e5_cool.py:
P1 = photon at o x two-hole lock eigenstates; P2 = the photon converted
with a cell-0 quantum (cell 0 full -> a third hole there; a cell-0 hole
cannot convert its last quantum: the floor), c_b = 1.  Cold photons (the
lowest delta of the band).  Master equation over photons from the least
bound two-hole state below the band of unbound states into the bound pair.
usage: k1_capture.py delta"""
import sys, time, itertools, functools
print = functools.partial(print, flush=True)
import numpy as np
here = __file__.rsplit('/', 1)[0]
sys.path.insert(0, here + '/../../perm/cas')
import graphm as GM, lo

t0 = time.time()
delta = float(sys.argv[1]) if len(sys.argv) > 1 else 0.01
M, n = 2, 6; P = GM.P_2
g = 0.1*GM.DELTA; J0 = 2*g*g/GM.DELTA; TOUT = 0.9; t = 0.1*J0
mob = [lo.Local(M, 'mob') for _ in range(n)]
lf = lo.labfun_factory(M, True)
T = lo.chain_terms(mob, P, g, t=t, ring=True)
# two-hole lock, all S^z sectors (2S^z in {-2, 0, 2})
chE, chV, chB, chQ = [], [], [], []
blocks = {}
for sz in (-2, 0, 2):
    H, basis, dims, strides = lo.sector(mob, T, (sz, 2), lf)
    w, V = np.linalg.eigh(H.toarray())
    blocks[sz] = (w, V, basis, dims, strides)
print(f"# K1-3 capture by cold light: 3+2 n={n} ring, two mobile holes, t=0.1 J0, delta={delta}  "
      f"(two-hole sector dims: {[len(blocks[s][0]) for s in (-2, 0, 2)]})  [{time.time()-t0:.0f}s]", flush=True)
# P2: three holes, one at cell 0 (cell 0 converted), sectors by 2S^z
def p2_sector(sz):
    H, basis, dims, strides = lo.sector(mob, T, (sz, 3), lf)
    # keep only states with a hole at cell 0
    dig0 = (basis//strides[0]) % dims[0]
    keep = mob[0].N[dig0] == 1
    idx = np.where(keep)[0]
    Hk = H[idx][:, idx]
    return Hk, basis[idx], dims, strides
# ring levels and their spin (by recurrence across sectors)
allE = np.concatenate([blocks[s][0] for s in (-2, 0, 2)])
E0 = blocks[0][0][0]
def spin_of(e):
    return 1 if np.any(np.abs(blocks[2][0] - e) < 1e-9) else 0
w0 = blocks[0][0]
S0 = np.array([spin_of(e) for e in w0])
print(f"  two-hole levels (S^z = 0 sector, E - E0 in J0): " + " ".join(f"{(e-E0)/J0:.3f}(S={s})" for e, s in zip(w0[:12], S0[:12])))
# pair distance of each S^z=0 eigenstate (to identify bound vs unbound)
w, V, basis, dims, strides = blocks[0]
hpos = np.array([[mob[i].N[(c//strides[i]) % dims[i]] == 1 for i in range(n)] for c in basis])
dist = np.array([min(abs(np.where(h)[0][1] - np.where(h)[0][0]), n - abs(np.where(h)[0][1] - np.where(h)[0][0])) for h in hpos])
meand = np.array([float(np.abs(V[:, k])**2 @ dist) for k in range(len(w))])
print(f"  mean pair distance of the lowest levels: " + " ".join(f"{meand[k]:.2f}" for k in range(10)))
# --- S-matrix machinery (blocked by total 2S^z of lock + photon) ---
def sigma(x, tt):
    x = np.asarray(x, float); ins = np.abs(x) < 2*tt
    o = np.empty(x.shape, complex); o[ins] = (x[ins] - 1j*np.sqrt(4*tt*tt - x[ins]**2))/2
    xo = x[~ins]; o[~ins] = (xo - np.sign(xo)*np.sqrt(xo**2 - 4*tt*tt))/2
    return o
off = GM.OMEGA - GM.EF[2]
pq = [1, -1]                                    # photon mode 0: +1/2, mode 1: -1/2 (2S^z)
def channel_set(qt):
    chs = []
    for s in (0, 1):
        sz = qt - pq[s]
        if sz in blocks:
            w, V, basis, dims, strides = blocks[sz]
            for k in range(len(w)): chs.append((sz, s, k))
    return chs
cache = {}
def block(qt):
    if qt in cache: return cache[qt]
    chs = channel_set(qt)
    e1 = np.array([blocks[sz][0][k] for (sz, s, k) in chs])
    Hd, bd, dd, sd = p2_sector(qt)
    lam, U = np.linalg.eigh(Hd.toarray()); lam = lam + off
    didx = {c: k for k, c in enumerate(bd)}
    Vm = np.zeros((len(bd), len(chs)))
    for col, (sz, s, k) in enumerate(chs):
        w, Vl, basis, dims, strides = blocks[sz]
        dig0 = (basis//strides[0]) % dims[0]
        rest = basis - dig0*strides[0]
        vec = np.zeros(len(bd))
        for mu in range(M):
            p = P[mu, s]
            if p == 0: continue
            for k0 in range(dims[0]):
                occ = mob[0].occ[k0]
                if mob[0].N[k0] != 2 or occ[mu] == 0: continue
                o2 = list(occ); o2[mu] -= 1; kd = mob[0].idx[tuple(o2)]
                sel = np.where(dig0 == k0)[0]
                for si in sel:
                    c = kd*strides[0] + rest[si]
                    if c in didx: vec[didx[c]] += g*p*np.sqrt(occ[mu])*Vl[si, k]
        Vm[:, col] = vec
    Vt = U.T@Vm
    cache[qt] = (chs, e1, lam, Vt)
    return cache[qt]
KLOW = 40   # cold photons only lower the energy: the lowest KLOW levels per sector carry the flow
levels = [(sz, k) for sz in (-2, 0, 2) for k in range(min(KLOW, len(blocks[sz][0])))]
lid = {lv: i for i, lv in enumerate(levels)}
NL = len(levels)
def trans(eps_list):
    Mx = np.zeros((NL, NL)); lost = np.zeros(NL)
    for eps in eps_list:
        for (sz, k) in levels:
            for s in (0, 1):
                qt = sz + pq[s]
                chs, e1, lam, Vt = block(qt)
                a = chs.index((sz, s, k)); E = blocks[sz][0][k] + eps
                w_ = E - e1
                if not abs(w_[a]) < 2*TOUT: continue
                S1 = sigma(w_, TOUT); D2 = 1.0/(E - lam + 1e-13j)
                A = np.diag(E - e1 - S1) - (Vt.T*D2)@Vt
                r = np.zeros(len(chs), complex); r[a] = 1
                Gc = np.linalg.solve(A, r); Gam = -2*S1.imag
                op = np.where(np.abs(w_) < 2*TOUT)[0]
                col = -1j*np.sqrt(Gam[op])*Gc[op]*np.sqrt(Gam[a]); col[np.where(op == a)[0][0]] += 1
                for kk, c in enumerate(op):
                    szf, sf, kf = chs[c]
                    if (szf, kf) in lid:
                        Mx[lid[(szf, kf)], lid[(sz, k)]] += np.abs(col[kk])**2/2/len(eps_list)
                    else:
                        lost[lid[(sz, k)]] += np.abs(col[kk])**2/2/len(eps_list)
    print(f"  flow out of the kept levels (heating above level {KLOW}): max {lost.max():.2e}")
    Mx = Mx/Mx.sum(0)[None, :]
    return Mx
eps = np.linspace(-2*TOUT + delta/20, -2*TOUT + delta, 5)
Mx = trans(eps)
# start: the least bound two-hole state among the S=0 levels below the first level with mean distance... use the
# highest S^z=0 level whose energy is below E0 + 1 gap of the no-hole ring (the "unbound" band), report its distance
gap_ring = None
Hn, *_ = lo.sector(mob, T, (0, 0), lf); Hn1, *_ = lo.sector(mob, T, (2, 0), lf)
gap_ring = np.linalg.eigvalsh(Hn1.toarray())[0] - np.linalg.eigvalsh(Hn.toarray())[0]
cand = [k for k in range(len(w0)) if w0[k] - E0 < gap_ring]
start = cand[-1]
print(f"  no-hole ring gap {gap_ring/J0:.3f} J0; start = S^z=0 level {start} at {(w0[start]-E0)/J0:.3f} J0 (S={S0[start]}, mean distance {meand[start]:.2f})")
p = np.zeros(NL); p[lid[(0, start)]] = 1.0
band = [k for k in range(min(KLOW, len(w0))) if (w0[k] - E0) < 0.05*J0 and meand[k] > 2.9]
gsel = np.array([lid[(0, k)] for k in band])
print(f"  the bound-pair band: S^z=0 levels {band} (E - E0 < 0.05 J0, mean distance > 2.9)")
hit = None; Pk = Mx.copy(); q = p.copy()
for kexp in range(0, 44):
    q = Pk@p
    if hit is None and q[gsel].sum() >= 0.9: hit = 2**kexp
    if kexp in (0, 5, 10, 15, 20, 25, 30, 35, 40): print(f"    N = 2^{kexp}: bound-pair population {q[gsel].sum():.4f}")
    Pk = Pk@Pk
up = 1 - Mx[np.ix_(gsel, gsel)].sum(0).min()
print(f"  photons to bound-pair population >= 0.9: {'<= %.3g' % hit if hit else 'never within 2^43'};  per-photon probability that the bound pair is excited: {up:.2e}")
print(f"done [{time.time()-t0:.0f}s]")
