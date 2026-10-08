#!/usr/bin/env python3
"""E4 — the static consequence of the geometric ceiling C-geo (../PREREG.md §3 E4).

Under C-geo a cell at its ceiling holds z quanta (z = its degree), so its
spin is z/2.
(a) S-law two-cell clusters for caps (z, z), z = 2, 3, 4, 6 (des Cloizeaux):
    the projected two-spin spectrum, the bilinear fit J, and its residual.
(b) Gaps of open spin-S Heisenberg chains (DMRG, compact MPO), S = 1, 3/2, 2,
    with end spins that screen the edge states for integer S
    (S = 1: spin-1/2 ends; S = 2: spin-1 ends).
(c) Mixed degrees: an alternating (z = 2, 3) ring, n = 6, exact: the
    ground-state spin (Lieb–Mattis predicts a ferrimagnet).
"""
import sys, time, itertools
import numpy as np
here = __file__.rsplit('/', 1)[0]
sys.path.insert(0, here); sys.path.insert(0, here + '/../../phase1/cas'); sys.path.insert(0, here + '/../../lock/cas')
import graph as G
from dmrg import DMRG
from mpo3 import mpo_chain
import scipy.sparse as sp

t0 = time.time()
part = sys.argv[1] if len(sys.argv) > 1 else 'all'
g = 0.1*G.DELTA

def SdotS(s1, s2):
    Z1, P1 = G.spin_ops(s1); Z2, P2 = G.spin_ops(s2)
    return np.kron(Z1, Z2) + 0.5*(np.kron(P1, P2.T) + np.kron(P1.T, P2))

if part in ('all', 'a'):
    print(f"# E4(a) S-law two-cell clusters at caps (z, z), g/Delta = 0.1 (g = {g:.5f})")
    for z in (2, 3, 4, 6):
        gr = G.Graph(2, [(0, 1)], g, ceil='fixed', caps=[z, z])
        seeds = [(((a, z - a), (b, z - b)), (0,)) for a in range(z + 1) for b in range(z + 1)]
        H = gr.build(seeds).toarray() - 2*G.EF[z]*np.eye(len(gr.states))
        w, V = np.linalg.eigh(H)
        model = [gr.idx[(((z - a, a), (z - b, b)), (0,))] for a in range(z + 1) for b in range(z + 1)]
        m = len(model)
        X = V[model, :m]
        ev, U = np.linalg.eigh(X@X.T)
        Sm = U@np.diag(ev**-0.5)@U.T
        Heff = Sm@X@np.diag(w[:m])@X.T@Sm; Heff = 0.5*(Heff + Heff.T)
        s = z/2
        SS = SdotS(s, s)
        A = np.stack([np.eye(m).ravel(), SS.ravel()], 1)
        coef, *_ = np.linalg.lstsq(A, Heff.ravel(), rcond=None)
        resid = np.abs(Heff - coef[0]*np.eye(m) - coef[1]*SS).max()
        # projected spectrum by total spin
        ws = np.linalg.eigvalsh(Heff)
        Dz = G.OMEGA - 2*(G.EF[z] - G.EF[z - 1])
        print(f"  z={z} (spin {s}): Delta_z={Dz:.5f}  sep to converted={w[m]-w[m-1]:.5f}  J_bilinear={coef[1]:.6f} "
              f"(= {coef[1]/(2*g*g/Dz):.3f} x 2g^2/Delta_z)  max non-bilinear residual={resid:.2e} ({resid/abs(coef[1]):.3f} J)  "
              f"levels: " + " ".join(f"{x:+.5f}" for x in np.unique(np.round(ws, 9))), flush=True)

if part in ('all', 'b'):
    print(f"\n# E4(b) open spin-S Heisenberg chains, J = 1, gap E(Sz=1) - E(Sz=0)  [{time.time()-t0:.0f}s]")
    for (S, ends, Ls, chi) in ((1.5, None, (8, 16, 32), 48), (1.0, 0.5, (8, 32), 48)):
        for L in Ls:
            spins = ([ends] if ends else []) + [S]*L + ([ends] if ends else [])
            dims = [int(round(2*x + 1)) for x in spins]
            h2 = [SdotS(spins[i], spins[i + 1]) for i in range(len(spins) - 1)] + [None]
            r3 = [None]*len(spins)
            qops = [G.spin_ops(x)[0] for x in spins]
            E = {}
            for sz in (0, 1):
                W, const = mpo_chain(dims, h2, r3, qops, 1.0, sz)
                d = DMRG(W, dims, chi=chi, seed=11 + sz)
                E[sz] = d.run(sweeps=10, etol=1e-9) + const
            print(f"  S={S} ends={ends} L_bulk={L:3d} chi={chi}: gap={E[1]-E[0]:.6f}  gap*L={((E[1]-E[0])*L):.4f}  [{time.time()-t0:.0f}s]", flush=True)

if part in ('all', 'c'):
    print(f"\n# E4(c) mixed degrees: alternating spins 1 and 3/2 (z = 2, 3) on an n = 6 ring, J = 1, exact  [{time.time()-t0:.0f}s]")
    spins = [1.0, 1.5]*3
    n = len(spins); dims = [int(2*s + 1) for s in spins]
    ops = [G.spin_ops(s) for s in spins]
    def at(op, i):
        return sp.kron(sp.kron(sp.identity(int(np.prod(dims[:i]))), sp.csr_matrix(op)), sp.identity(int(np.prod(dims[i + 1:]))), format='csr')
    Z = [at(o[0], i) for i, o in enumerate(ops)]; P = [at(o[1], i) for i, o in enumerate(ops)]
    H = sp.csr_matrix(Z[0].shape)
    for i in range(n):
        j = (i + 1) % n
        H = H + Z[i]@Z[j] + 0.5*(P[i]@P[j].T + P[i].T@P[j])
    w, v = np.linalg.eigh(H.toarray())
    Ztot = sum(Z); Ptot = sum(P)
    S2 = (Ztot@Ztot + 0.5*(Ptot@Ptot.T + Ptot.T@Ptot)).toarray()
    s2 = [float(v[:, k]@S2@v[:, k]) for k in range(8)]
    Sval = [-0.5 + np.sqrt(0.25 + x) for x in s2]
    print(f"  lowest levels: " + " ".join(f"{a:+.5f}(S={b:.2f})" for a, b in zip(w[:8], Sval)))
    print(f"  ground-state spin {Sval[0]:.3f}; Lieb–Mattis: n/2 x |1 - 3/2| = {n/2*0.5:.2f}  (a ferrimagnet: gapless spin waves)")
print(f"done [{time.time()-t0:.0f}s]")
