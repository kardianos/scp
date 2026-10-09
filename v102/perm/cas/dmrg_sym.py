#!/usr/bin/env python3
"""DMRG of ../../phase1/cas/dmrg.py with one change: truncation never splits
a block of (near-)degenerate singular values.  A cut inside a symmetry
multiplet breaks the symmetry (here SU(2) x SU(2) of W), so the kept number
is lowered to the last gap in the singular values below chi (relative gap
> 1e-6).  Everything else is the carried code."""
import sys
import numpy as np
here = __file__.rsplit('/', 1)[0]
sys.path.insert(0, here + '/../../phase1/cas')
import dmrg as _d

def _keep(S, chi, cutoff):
    keep = max(1, min(chi, int(np.sum(S > cutoff*S[0]))))
    while 1 < keep < len(S) and (S[keep - 1] - S[keep]) < 1e-6*S[0]:
        keep -= 1
    return keep

class DMRGSym(_d.DMRG):
    def sweep(self, tol=1e-10, cutoff=1e-12):
        E = None; trunc = 0.0
        for i in range(self.N - 1):
            E, th = self._solve(i, tol)
            Dl, d1, d2, Dr = th.shape
            U, S, V = np.linalg.svd(th.reshape(Dl*d1, d2*Dr), full_matrices=False)
            keep = _keep(S, self.chi, cutoff)
            trunc = max(trunc, float(np.sum(S[keep:]**2)))
            U, S, V = U[:, :keep], S[:keep]/np.linalg.norm(S[:keep]), V[:keep]
            self.A[i] = U.reshape(Dl, d1, keep)
            self.A[i + 1] = (np.diag(S)@V).reshape(keep, d2, Dr)
            self.L[i + 1] = _d._lenv(self.L[i], self.A[i], self.W[i])
        for i in range(self.N - 2, -1, -1):
            E, th = self._solve(i, tol)
            Dl, d1, d2, Dr = th.shape
            U, S, V = np.linalg.svd(th.reshape(Dl*d1, d2*Dr), full_matrices=False)
            keep = _keep(S, self.chi, cutoff)
            trunc = max(trunc, float(np.sum(S[keep:]**2)))
            U, S, V = U[:, :keep], S[:keep]/np.linalg.norm(S[:keep]), V[:keep]
            self.A[i] = (U@np.diag(S)).reshape(Dl, d1, keep)
            self.A[i + 1] = V.reshape(keep, d2, Dr)
            self.R[i + 1] = _d._renv(self.R[i + 2], self.A[i + 1], self.W[i + 1])
        self.last_trunc = trunc
        return E, trunc
