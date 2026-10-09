#!/usr/bin/env python3
"""X4 — S: hole formation per resonant photon, for 3+2 (M = 2, calibration
against ../../options/ E5's exact numbers) and W (M = 4).

n = 4 ring at leading order (reading V for the ring), one boundary link
(cell 0 to the end site o of a photon lead, t_out = 0.9).  P1 = photon at o
(mode nu) x ring eigenstates; P2 = boundary space quantum c_b = 1 x
depleted ring (cell 0 at N = 1) eigenstates.  The photon lead and (decision
S) a 1-D space lead at t_s = s_k are integrated out as self-energies.
Channels are restricted by the conserved charges of ring + photon.
Absorption (= a hole made, the space quantum leaving) per resonant photon
from the ring ground state, polarization averaged, eps in eps_res +- 0.05.
usage: x4_form.py M gr"""
import sys, time
import numpy as np
here = __file__.rsplit('/', 1)[0]
sys.path.insert(0, here)
import graphm as GM, lo

t0 = time.time()
M, gr = int(sys.argv[1]), float(sys.argv[2])
P = GM.P_W if M == 4 else GM.P_2
g = gr*GM.DELTA; n = 4; TOUT, TS = 0.9, 0.06
full = [lo.Local(M, 'full') for _ in range(n)]
dep = [lo.Local(M, 'hole')] + [lo.Local(M, 'full') for _ in range(n - 1)]
lf = lo.labfun_factory(M, False)
Tr = lo.chain_terms(full, P, g, ring=True)
Td = lo.chain_terms(dep, P, g, ring=True)
pmodes = [lo.Local(M, 'hole').lab[k] for k in range(M)]          # charges of a one-quantum state (mode order)
pocc = lo.Local(M, 'hole').occ
eres = GM.OMEGA - GM.EF[2]
off = GM.OMEGA - GM.EF[2]       # P2 offset: Omega + (n-1)E_f(2) + 0 - [n E_f(2) + E_f(1)] + E_f(1)... see below

def sigma(w, t):
    w = np.asarray(w, dtype=float); ins = np.abs(w) < 2*t
    o = np.empty(w.shape, complex)
    o[ins] = (w[ins] - 1j*np.sqrt(4*t*t - w[ins]**2))/2
    wo = w[~ins]; o[~ins] = (wo - np.sign(wo)*np.sqrt(wo**2 - 4*t*t))/2
    return o

# ring ground state
z = tuple([0]*(2 if M == 4 else 1))
H0, b0, dims0, st0 = lo.sector(full, Tr, z, lf)
w0, v0 = np.linalg.eigh(H0.toarray()) if H0.shape[0] < 6000 else None
E0 = w0[0]; g0 = v0[:, 0]
print(f"# X4 formation under S: M={M}, g/Delta={gr}, n={n}; ring ground {E0:+.6f} (sector dim {H0.shape[0]})", flush=True)

def digits(code, dims, strides):
    return [(code//s) % d for s, d in zip(strides, dims)]

absorb = []
for nu in range(M):                   # incoming photon polarization
    qtot = np.array(pmodes[nu])
    # P1 channels: ring sector q_r = qtot - q_photon(nu'), photon mode nu'
    ch_E, ch_vec, ch_info = [], [], []
    for nup in range(M):
        qr = tuple(qtot - np.array(pmodes[nup]))
        H, b, dims, st = lo.sector(full, Tr, qr, lf)
        w, V = np.linalg.eigh(H.toarray())
        for k in range(len(w)):
            ch_E.append(w[k]); ch_info.append((nup, qr, k))
        ch_vec.append((nup, qr, b, dims, st, V))
    ch_E = np.array(ch_E)
    # P2: depleted ring in sector qtot (the space pair is a double singlet)
    Hd, bd, dimsd, std = lo.sector(dep, Td, tuple(qtot), lf)
    lam, U = np.linalg.eigh(Hd.toarray())
    lam = lam + off                      # energy of (c_b, depleted ring) relative to (ring, photon at rest at o)
    dindex = {c: k for k, c in enumerate(bd)}
    # V[p2 basis, channel]: photon nu' at o + ring basis state -> depleted basis state
    nch = len(ch_E); Vm = np.zeros((len(bd), nch))
    col0 = 0
    for (nup, qr, b, dims, st, Vr) in ch_vec:
        # ring basis -> depleted basis: cell 0 loses mode mu, photon nu' absorbed
        Lf = full[0]; Ld = dep[0]
        A = np.zeros((len(bd), len(b)))
        dig0 = (b//st[0]) % dims[0]
        rest = b - dig0*st[0]
        for mu in range(M):
            p = P[mu, nup]
            if p == 0: continue
            for k0 in range(dims[0]):
                occ = Lf.occ[k0]
                if occ[mu] == 0: continue
                o2 = list(occ); o2[mu] -= 1; o2 = tuple(o2)
                kd = Ld.idx[o2]
                sel = np.where(dig0 == k0)[0]
                newc = rest[sel]//st[0]*0 + kd*std[0] + (rest[sel] % st[0])   # same lower digits, new cell-0 digit
                for s, c in zip(sel, newc):
                    if c in dindex: A[dindex[c], s] += g*p*np.sqrt(occ[mu])
        Vm[:, col0:col0 + Vr.shape[1]] = A@Vr
        col0 += Vr.shape[1]
    Vt = U.T@Vm
    # incoming channel: ring ground (sector z) with photon nu
    a_in = [i for i, (nup, qr, k) in enumerate(ch_info) if nup == nu and qr == z and k == 0][0]
    probs = []
    for eps in np.linspace(eres - 0.05, eres + 0.05, 11):
        E = E0 + eps
        S1 = sigma(E - ch_E, TOUT)
        S2 = sigma(E - lam, TS)
        D2 = 1.0/(E - lam - S2 + 1e-13j)
        A = np.diag(E - ch_E - S1) - (Vt.T*D2)@Vt
        Gam = -2*S1.imag
        r = np.zeros(nch, complex); r[a_in] = 1
        Gc = np.linalg.solve(A, r)
        op = np.where(np.abs(E - ch_E) < 2*TOUT)[0]
        col = -1j*np.sqrt(Gam[op])*Gc[op]*np.sqrt(Gam[a_in]); col[np.where(op == a_in)[0][0]] += 1
        probs.append(1 - (np.abs(col)**2).sum())
    absorb.append(np.mean(probs))
    print(f"  photon mode {nu} (charges {tuple(pmodes[nu])}): channels {nch}, P2 dim {len(lam)}; hole made per resonant photon = {np.mean(probs):.4f}  [{time.time()-t0:.0f}s]", flush=True)
print(f"  polarization average: {np.mean(absorb):.4f}")
print(f"done [{time.time()-t0:.0f}s]")
