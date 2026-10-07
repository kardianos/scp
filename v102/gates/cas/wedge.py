#!/usr/bin/env python3
"""v102 G-wedge instrument.  See ../../PREREG.md §2 and ../WEDGE.md.

W1  invariant conversion operators.  Field: two bosonic modes per cell
    (R1, R2).  Space: a scalar mode.  Which operators that remove field
    quanta are invariant (a) under SO(2) rotations about the link axis,
    (b) under SU(2) (R1)?  Does the on-site singlet survive for bosons?
W2  two-site Schwinger spaces, N quanta per site, N = 1..4: is A^+A
    (conversion, A = b_i1 b_j2 - b_i2 b_j1) affine in S_i.S_j?  Is B^+B
    (transport, B = sum_a b_ja^+ b_ia) affine in S_i.S_j?
M   microscopic ring: cells with N0 = 2, one scalar space mode per link,
    H = Delta sum c^+c + g sum (c^+ A + h.c.), locked-link transport 0.
    Low spectrum vs the derived Heisenberg ring.  Delta < 0 control.
"""
import itertools
import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import eigsh

rng = np.random.default_rng(102)
np.set_printoptions(precision=6, suppress=True)

# ------------------------------------------------------------------ W1
print("# W1  invariant conversion operators")
def rot(th):  # rotation about the link axis acting on the two plane amplitudes
    return np.array([[np.cos(th), -np.sin(th)], [np.sin(th), np.cos(th)]])
def rand_su2():
    a = rng.normal(size=4); a /= np.linalg.norm(a)
    return np.array([[a[0] + 1j*a[1], -a[2] + 1j*a[3]], [a[2] + 1j*a[3], a[0] - 1j*a[1]]])

def invariant_bilinears(group_elems):
    # M with U^T M U = M for all U  (pair operator M_ab b_ia b_jb)
    rows = []
    for U in group_elems:
        T = np.kron(U.T, U.T) - np.eye(4)          # vec(U^T M U) = (U^T kron U^T) vec(M)  (row-major)
        rows.append(T)
    A = np.vstack(rows)
    _, s, vh = np.linalg.svd(A)
    null = vh[np.sum(s > 1e-9):]
    return [v.reshape(2, 2) for v in null]

def invariant_linears(group_elems):
    rows = [U.T - np.eye(2) for U in group_elems]
    _, s, vh = np.linalg.svd(np.vstack(rows))
    return vh[np.sum(s > 1e-9):]

so2 = [rot(t) for t in rng.uniform(0, 2*np.pi, 6)]
su2 = [rand_su2() for _ in range(6)]
P = np.diag([1.0, -1.0])                           # reflection of one plane
for name, G in (("SO(2) about the link", so2), ("SU(2) (R1)", su2)):
    B = invariant_bilinears(G)
    L = invariant_linears(G)
    print(f" {name:22s}: invariant linear forms {len(L)}, invariant bilinear pair forms {len(B)}")
    for M in B:
        M = M/np.abs(M).max()
        sym = np.allclose(M, M.T); anti = np.allclose(M, -M.T)
        par = "even" if np.allclose(P.T@M@P, M) else ("odd" if np.allclose(P.T@M@P, -M) else "mixed")
        kind = "delta-like (symmetric)" if sym else ("epsilon (antisymmetric)" if anti else "mixed")
        print(f"      M = {np.round(M.real, 6).tolist()}  {kind}, reflection-{par}")
# on-site singlet for bosons
def fock_ops(nmax):
    a = np.diag(np.sqrt(np.arange(1, nmax + 1)), 1)
    return a
nm = 4
a = fock_ops(nm); I = np.eye(nm + 1)
b1 = np.kron(a, I); b2 = np.kron(I, a)
print(f" on-site singlet b1 b2 - b2 b1 for bosons: max |entry| = {np.abs(b1@b2 - b2@b1).max():.1e}")
print(f" on-site delta pair b1 b1 + b2 b2:          max |entry| = {np.abs(b1@b1 + b2@b2).max():.2f}")

# ------------------------------------------------------------------ W2
print()
print("# W2  two-site Schwinger spaces: A^+A and B^+B against S_i.S_j")
def two_site(N):
    # basis: (n_i1, n_j1) with n_i2 = N - n_i1 ... but A, B leave the (N,N) space; work in the
    # full two-site Fock space with per-mode cutoff N+1, then restrict to (N,N).
    nmax = N + 1
    a = fock_ops(nmax); I = np.eye(nmax + 1)
    ops = []
    for k in range(4):
        mats = [I]*4; mats[k] = a
        M = mats[0]
        for m in mats[1:]: M = np.kron(M, m)
        ops.append(M)
    bi1, bi2, bj1, bj2 = ops
    ni = bi1.T@bi1 + bi2.T@bi2; nj = bj1.T@bj1 + bj2.T@bj2
    Si = [0.5*(bi1.T@bi2 + bi2.T@bi1), -0.5j*(bi1.T@bi2 - bi2.T@bi1), 0.5*(bi1.T@bi1 - bi2.T@bi2)]
    Sj = [0.5*(bj1.T@bj2 + bj2.T@bj1), -0.5j*(bj1.T@bj2 - bj2.T@bj1), 0.5*(bj1.T@bj1 - bj2.T@bj2)]
    SS = sum(Si[k]@Sj[k] for k in range(3))
    A = bi1@bj2 - bi2@bj1
    B = bj1.T@bi1 + bj2.T@bi2
    keep = np.where((np.abs(np.diag(ni) - N) < 1e-9) & (np.abs(np.diag(nj) - N) < 1e-9))[0]
    r = lambda X: X[np.ix_(keep, keep)]
    return r(SS), r(A.T@A), r(B.T.conj()@B)

print("  N  S     A+A = a + b S.S          resid      B+B = a + b S.S          resid")
for N in range(1, 5):
    SS, AA, BB = two_site(N)
    out = []
    for X in (AA, BB):
        # least squares over matrix entries: X ~ a*1 + b*SS
        D = SS.shape[0]
        Mx = np.stack([np.eye(D).ravel(), SS.ravel()], 1)
        coef, *_ = np.linalg.lstsq(Mx, X.ravel(), rcond=None)
        res = np.abs(X - (coef[0]*np.eye(D) + coef[1]*SS)).max()
        out.append((coef.real, res))
    (ca, ra), (cb, rb) = out
    print(f"  {N}  {N/2:3.1f}   {ca[0]:+8.4f} {ca[1]:+8.4f} S.S    {ra:.1e}     {cb[0]:+8.4f} {cb[1]:+8.4f} S.S    {rb:.1e}")
print("  -(g^2/Delta) A+A  = const + (g^2/Delta) * |b_A| S.S   -> antiferromagnetic when Delta > 0")
print("  -(t^2/U)    B+B  = const - (t^2/U)     * |b_B| S.S   -> ferromagnetic (the E1 competitor)")

# ------------------------------------------------------------------ M
print()
print("# M   microscopic ring: N0 = 2 per cell, scalar space mode per link, transport 0 at the ceiling")
MICRO_S2 = [None]
def micro_ring(n, N0, g, Delta):
    # link l joins cell l and cell (l+1)%n.  c_l = space quanta on link l.
    basis = []
    for cs in itertools.product(range(N0 + 1), repeat=n):
        Ns = [N0 - cs[(i - 1) % n] - cs[i] for i in range(n)]
        if min(Ns) < 0: continue
        for n1 in itertools.product(*[range(Ni + 1) for Ni in Ns]):
            basis.append((cs, n1, tuple(Ns)))
    idx = {(b[0], b[1]): k for k, b in enumerate(basis)}
    rows, cols, vals = [], [], []
    for k, (cs, n1, Ns) in enumerate(basis):
        rows.append(k); cols.append(k); vals.append(Delta*sum(cs))
        for l in range(n):
            i, j = l, (l + 1) % n
            # A_l = b_i1 b_j2 - b_i2 b_j1 ; apply c_l^+ A_l
            for (ai, aj, sgn) in ((1, 2, +1), (2, 1, -1)):
                ni_a = n1[i] if ai == 1 else Ns[i] - n1[i]
                nj_a = n1[j] if aj == 1 else Ns[j] - n1[j]
                if ni_a == 0 or nj_a == 0: continue
                amp = sgn*np.sqrt(ni_a)*np.sqrt(nj_a)*np.sqrt(cs[l] + 1)*g
                cs2 = list(cs); cs2[l] += 1
                n12 = list(n1)
                if ai == 1: n12[i] -= 1
                if aj == 1: n12[j] -= 1
                key = (tuple(cs2), tuple(n12))
                if key not in idx: continue
                k2 = idx[key]
                rows += [k2, k]; cols += [k, k2]; vals += [amp, amp]
    D = len(basis)
    H = sp.csr_matrix((vals, (rows, cols)), shape=(D, D))
    # total S^2 is not needed: use the S^z-resolved spectrum via S^z
    Sz = np.array([0.5*sum(2*n1[i] - Ns[i] for i in range(n)) for (cs, n1, Ns) in basis])
    space = np.array([sum(cs) for (cs, n1, Ns) in basis])
    # total S^+ over field quanta: sum_i b_i1^+ b_i2 (space is scalar)
    r2, c2, v2 = [], [], []
    for k, (cs, n1, Ns) in enumerate(basis):
        for i in range(n):
            n2 = Ns[i] - n1[i]
            if n2 == 0: continue
            n12 = list(n1); n12[i] += 1
            r2.append(idx[(cs, tuple(n12))]); c2.append(k); v2.append(np.sqrt(n2)*np.sqrt(n1[i] + 1))
    Spl = sp.csr_matrix((v2, (r2, c2)), shape=(D, D))
    Szm = sp.diags(Sz)
    S2 = Spl@Spl.T + Szm@Szm - Szm
    MICRO_S2[0] = S2
    return H, Sz, space, D

def heis_ring_levels(n, s, k=8):
    from scipy.sparse import kron, identity, diags, lil_matrix, csr_matrix
    m = np.arange(s, -s - 1, -1); d = len(m)
    Sz = diags(m); Sp = lil_matrix((d, d))
    for q in range(1, d): Sp[q - 1, q] = np.sqrt(s*(s + 1) - m[q]*(m[q] + 1))
    Sp = csr_matrix(Sp)
    at = lambda op, i: kron(kron(identity(d**i), op), identity(d**(n - i - 1)), format='csr')
    Z = [at(Sz, i) for i in range(n)]; Pp = [at(Sp, i) for i in range(n)]
    H = csr_matrix((d**n, d**n))
    for i in range(n):
        j = (i + 1) % n
        H = H + Z[i]@Z[j] + 0.5*(Pp[i]@Pp[j].T + Pp[i].T@Pp[j])
    w = np.linalg.eigvalsh(H.toarray()) if d**n <= 1000 else np.sort(eigsh(H, k=k, which='SA')[0])
    return np.sort(w)

def distinct(w, tol=1e-7):
    out = [w[0]]
    for x in w[1:]:
        if x - out[-1] > tol: out.append(x)
    return np.array(out)

for n in (4, 6):
    hl = distinct(heis_ring_levels(n, 1.0, k=30))
    hr = (hl[2] - hl[0])/(hl[1] - hl[0])
    print(f" n={n}: spin-1 Heisenberg ring, (E2-E0)/(E1-E0) over distinct levels = {hr:.5f}")
    print("    g/Delta   dim    E1-E0      J_eff=(E1-E0)/gapHeis   ratio(E2-E0)/(E1-E0)   rel.miss   <space> in ground   S_tot(ground)")
    for gD, Delta in ((0.05, 1.0), (0.1, 1.0), (0.2, 1.0), (0.4, 1.0), (0.1, -1.0)):
        H, Sz, space, D = micro_ring(n, 2, gD*abs(Delta), Delta)
        w, v = eigsh(H, k=min(40, D - 2), which='SA') if D > 1500 else np.linalg.eigh(H.toarray())
        o = np.argsort(w); w = w[o]; v = v[:, o]
        lv = distinct(w)
        r = (lv[2] - lv[0])/(lv[1] - lv[0])
        Jeff = (lv[1] - lv[0])/(hl[1] - hl[0])
        sp0 = float(np.sum(np.abs(v[:, 0])**2*space))
        g0 = np.where(w < w[0] + 1e-8)[0]
        blk = v[:, g0].T@(MICRO_S2[0]@v[:, g0]); s2 = np.linalg.eigvalsh((blk + blk.T)/2)
        Sg = -0.5 + np.sqrt(0.25 + np.clip(s2, 0, None))
        Sg = f"{Sg.min():.3f}" if np.ptp(Sg) < 1e-6 else f"{Sg.min():.2f}..{Sg.max():.2f}"
        tag = "  (Delta<0 control)" if Delta < 0 else ""
        print(f"    {gD:5.2f}   {D:6d}   {lv[1]-lv[0]:9.6f}   {Jeff:12.6f}            {r:9.5f}         {abs(r-hr)/hr:8.2e}   {sp0:7.4f}          {Sg}{tag}")
print(" J_eff predicted at second order: |b_A| g^2/Delta with |b_A| from W2 (N=2).")
