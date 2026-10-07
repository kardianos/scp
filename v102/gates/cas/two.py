#!/usr/bin/env python3
"""v102 G-two instrument.  See ../../PREREG.md §3 and ../TWO.md.

With transport refused at the ceiling, the derived link Hamiltonian is
H = -(g^2/Delta) sum_l A_l^+ A_l = (g^2/Delta) sum_l (2 S_i.S_j - N^2/2)
(W2).  For uniform N per cell (S = N/2) on a ring of n cells, the energy
per quantum is e(N) = (2 E0_Heis(n, S)/n - N^2/2) / N in units g^2/Delta.
If e(N) falls with N, quanta hoard: the occupation runs up to whatever
ceiling exists.  The premises do not fix N by themselves.
"""
import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import eigsh

def heis_E0(n, s):
    m = np.arange(s, -s - 1, -1); d = len(m)
    Sz = sp.diags(m); Sp = sp.lil_matrix((d, d))
    for q in range(1, d): Sp[q - 1, q] = np.sqrt(s*(s + 1) - m[q]*(m[q] + 1))
    Sp = sp.csr_matrix(Sp)
    at = lambda op, i: sp.kron(sp.kron(sp.identity(d**i), op), sp.identity(d**(n - i - 1)), format='csr')
    Z = [at(Sz, i) for i in range(n)]; P = [at(Sp, i) for i in range(n)]
    H = sp.csr_matrix((d**n, d**n))
    for i in range(n):
        j = (i + 1) % n
        H = H + Z[i]@Z[j] + 0.5*(P[i]@P[j].T + P[i].T@P[j])
    keep = np.where(np.abs(sum(Z).diagonal() - (0 if (s*n) % 1 == 0 else 0.5)) < 1e-9)[0]
    Hs = H[keep][:, keep]
    if Hs.shape[0] < 500:
        return np.linalg.eigvalsh(Hs.toarray())[0]
    return eigsh(Hs, k=1, which='SA', tol=1e-10)[0][0]

print("# G-two instrument: energy per quantum vs occupation N (units g^2/Delta)")
print("  N   S     n   E0_Heis/n     e(N) per quantum    gapped? (integer S, Haldane)")
for n in (6, 8):
    for N in (1, 2, 3, 4):
        s = N/2
        e0 = heis_E0(n, s)/n
        eN = (2*e0 - N**2/2)/N
        print(f"  {N}  {s:3.1f}   {n}   {e0:10.5f}    {eN:12.5f}          {'yes' if N % 2 == 0 else 'no (LSM)'}")
print("  e(N) falls monotonically with N: the conversion bond hoards quanta.")
print("  With no ceiling the occupation runs away; with a ceiling N_c, every cell sits at N_c.")
print("  A gapped lock needs N_c even (S integer): N_c = 2 -> xi ~ 6.03, N_c = 4 -> xi ~ 50 (literature).")
