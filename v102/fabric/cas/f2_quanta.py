#!/usr/bin/env python3
"""F2 — counting lock quanta (reading Q).

(a) One photon on the n = 4 lock in its ground state, by the exact E5
    machinery (../../options/cas/e5_cool.py, executed up to its source
    loop, so the S-matrix code is the same).  Final ring states are sorted:
    ground; exactly one lock quantum (S = 1 and below the two-quantum
    threshold 2 Delta_4); more (other spin-manifold levels); space-carrying
    levels (index >= 81, V's mismatch); destroyed (absorbed, F1 only).
(b) Sudden starts on effective rings (reading V, h2 + r3): the energy above
    the ground state per cell, in units of the infinite-chain gap
    Delta_inf = 0.4357 J (../../lock/ L4): the density of lock quanta made.
usage: f2_quanta.py a GR F | b
"""
import sys, time
import numpy as np
import scipy.sparse as sp
here = __file__.rsplit('/', 1)[0]
sys.path.insert(0, here)
t0 = time.time()
part = sys.argv[1]

if part == 'a':
    gr_, F = sys.argv[2], sys.argv[3]
    e5 = here + '/../../options/cas/e5_cool.py'
    src = open(e5).read()
    head = src[:src.index('sources = []')]
    sys.argv = [e5, gr_, F]
    ns = {'__file__': e5, '__name__': 'e5_head'}
    exec(compile(head, e5, 'exec'), ns)
    er, szf, NR, trans_matrix, G = ns['er'], ns['szf'], ns['NR'], ns['trans_matrix'], ns['G']
    # total spin of each ring level: the largest |m| at which its energy recurs
    S = np.zeros(NR)
    for i in range(NR):
        same = np.abs(er - er[i]) < 1e-9
        S[i] = np.abs(szf[same]).max()/2
    gap4 = er[1] - er[0]
    cls = np.full(NR, 'more', dtype=object)
    cls[0] = 'ground'
    for i in range(1, NR):
        if i >= 81: cls[i] = 'space'
        elif abs(S[i] - 1) < 1e-9 and er[i] - er[0] < 2*gap4 - 1e-12: cls[i] = 'one'
    eres = G.OMEGA - G.EF[2]
    TOUT = 0.9
    srcs = [("cold delta=0.005", np.linspace(-2*TOUT + 0.005/20, -2*TOUT + 0.005, 10)),
            ("resonant eps_res+-0.05", np.linspace(eres - 0.05, eres + 0.05, 11)),
            ("broadband", np.linspace(-2*TOUT + 0.01, 2*TOUT - 0.01, 37))]
    print(f"# F2(a) one photon on the n = 4 lock in its ground state: g/Delta = {gr_}, {F}; one-quantum gap {gap4:.5f}; "
          f"levels classed 'one' = {int((cls=='one').sum())}, 'more' = {int((cls=='more').sum())}, 'space' = {int((cls=='space').sum())}")
    for name, eps in srcs:
        M, absb, worst = trans_matrix(eps)
        cs = M.sum(0) + absb
        if F != 'F1':
            absb = np.zeros(NR); M = M/M.sum(0)[None, :]
        col = M[:, 0]
        out = {c: float(col[cls == c].sum()) for c in ('ground', 'one', 'more', 'space')}
        out['destroyed'] = float(absb[0])
        print(f"  {name:24s}: " + "  ".join(f"{k}={v:.4e}" for k, v in out.items()) +
              f"  (one / not-ground = {out['one']/max(1e-300, 1-out['ground']):.3f})  [{time.time()-t0:.0f}s]", flush=True)

if part == 'b':
    import eff
    h2, r3, co = eff.terms_for(0.1); J = co['J']
    Dinf = 0.4357*J
    print(f"# F2(b) sudden starts on effective rings (reading V), g/Delta = 0.1; Delta_inf = 0.4357 J")
    for n in (8, 10):
        caps = [2]*n
        T = eff.chain_terms(caps, h2, r3, ring=True)
        # full-space H as a sum of S^z sectors
        blocks = []; bases = []
        for tw in range(-2*n, 2*n + 1, 2):
            H, basis, dims, strides = eff.sector_hamiltonian(caps, T, tw)
            blocks.append(H); bases.append(basis)
        from scipy.sparse.linalg import eigsh
        H0 = blocks[n]                          # tw = 0
        E0 = eigsh(H0, k=1, which='SA', tol=1e-11)[0][0]
        allb = np.concatenate(bases); order = np.argsort(allb)
        Hf = sp.block_diag(blocks, format='csr')[order][:, order]
        D = 3**n
        rng = np.random.default_rng(1001)
        starts = [("neel", [np.array([1, 0, 0]) if i % 2 == 0 else np.array([0, 0, 1]) for i in range(n)])]
        for s in (1, 2, 3):
            cv = []
            for i in range(n):
                z = rng.normal(size=2) + 1j*rng.normal(size=2); z /= np.linalg.norm(z)
                cv.append(np.array([z[0]**2, np.sqrt(2)*z[0]*z[1], z[1]**2]))   # n1 = 2, 1, 0
            starts.append((f"jones{s}", cv))
        for name, cv in starts:
            v = cv[0]
            for c in cv[1:]: v = np.kron(v, c)
            e = float(np.real(np.vdot(v, Hf@v)))
            dens = (e - E0)/(n*Dinf)
            print(f"  n={n:2d} {name:6s}: <H> - E0 = {(e-E0)/J:7.3f} J  = {dens:.3f} lock quanta per cell  "
                  f"(mean spacing between quanta {1/max(dens,1e-12):.2f} cells)  [{time.time()-t0:.0f}s]", flush=True)
print(f"done [{time.time()-t0:.0f}s]")
