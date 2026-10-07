#!/usr/bin/env python3
"""v99 reopen — feasibility numbers for a generated-gap selector.

Three checks, no kernel, no law change:

  H  the v97 field hop and self-rotation, transcribed from
     v97/kernel/freecell.c (pass F), act on one complex amplitude per
     cell, are linear, and commute with a global phase.  The field
     sector the law supplies is U(1) and free.
  X  sensitivity of a transmuted size to its coupling, using the exact
     O(4) lattice mass asymptotic of openai/math family 215:
     m_lat(b) ~ 32 exp(pi/4 - 1/2) sqrt(b) exp(-pi b).
  W  winding-sector free energy on a U(1) ring of n sites with bond
     weight exp(K cos d) (exact n-fold convolution).  A generated
     correlation length bounds the wrap from above; it does not pick n.
"""
import numpy as np

rng = np.random.default_rng(20261007)
print("# H  field hop symmetry (v97 freecell.c pass F)")

def self_rot(a1, a2, ang):              # freecell.c:2409-2416
    cc, ss = np.cos(ang), np.sin(ang)
    return cc*a1 + ss*a2, -ss*a1 + cc*a2

def hop(a1i, a2i, a1j, a2j, tau):       # freecell.c:2438-2443
    cc, ss = np.cos(tau), np.sin(tau)
    return (cc*a1i + ss*a2j, cc*a2i - ss*a1j,
            cc*a1j + ss*a2i, cc*a2j - ss*a1i)

def z(a1, a2): return a1 + 1j*a2

worst_form = worst_phase = worst_lin = worst_norm = 0.0
for _ in range(2000):
    a = rng.normal(size=4); tau, ang, al = rng.uniform(-3, 3, 3)
    zi, zj = z(a[0], a[1]), z(a[2], a[3])
    o = hop(*a, tau)
    zi2, zj2 = z(o[0], o[1]), z(o[2], o[3])
    # closed form: z_i' = cos z_i - i sin z_j  (a beamsplitter)
    worst_form = max(worst_form, abs(zi2 - (np.cos(tau)*zi - 1j*np.sin(tau)*zj)),
                     abs(zj2 - (np.cos(tau)*zj - 1j*np.sin(tau)*zi)))
    # global phase commutes
    p = np.exp(1j*al)
    b = hop(*(lambda u, v: (u.real, u.imag, v.real, v.imag))(p*zi, p*zj), tau)
    worst_phase = max(worst_phase, abs(z(b[0], b[1]) - p*zi2), abs(z(b[2], b[3]) - p*zj2))
    # linear (superposition)
    c = rng.normal(size=4)
    s = hop(*(a + c), tau); s1 = hop(*a, tau); s2 = hop(*c, tau)
    worst_lin = max(worst_lin, max(abs(s[k] - s1[k] - s2[k]) for k in range(4)))
    worst_norm = max(worst_norm, abs(abs(zi2)**2 + abs(zj2)**2 - abs(zi)**2 - abs(zj)**2))
    r1, r2 = self_rot(a[0], a[1], ang)
    worst_form = max(worst_form, abs(z(r1, r2) - np.exp(-1j*ang)*zi))
print(f"H1 closed form z_i'=cos z_i - i sin z_j, self-rot = e^(-i ang) z : max err {worst_form:.1e}")
print(f"H2 commutes with global phase                                 : max err {worst_phase:.1e}")
print(f"H3 linear (superposition)                                     : max err {worst_lin:.1e}")
print(f"H4 pairwise |z|^2 conserved                                    : max err {worst_norm:.1e}")
print("H  verdict: one complex amplitude per cell, linear U(1)-covariant hop.")
print("   Real dof per cell = 2. No SU(2) can act on C^1. Field sector is abelian and free.")

print()
print("# X  transmuted size vs coupling, O(4) lattice (family 215, exact asymptotic)")
pref = 32*np.exp(np.pi/4 - 0.5)
print(" beta    xi=1/m_lat   dln(xi)/dln(beta)   note")
for b, note in [(1.0, "strong coupling, formula not valid"), (2.0, ""), (3.0, ""),
                (np.pi, "beta = pi"), (4.0, ""), (2*np.pi, "beta = 2 pi")]:
    m = pref*np.sqrt(b)*np.exp(-np.pi*b)
    print(f" {b:5.3f}  {1/m:12.4g}   {np.pi*b - 0.5:8.2f}            {note}")
xi = lambda b: 1/(pref*np.sqrt(b)*np.exp(-np.pi*b))
print(f"X  xi(2 pi)/xi(2) = {xi(2*np.pi)/xi(2):.3g}   (a factor pi in normalization)")
print(f"X  at beta=2 pi a 1% change in beta moves xi by {100*(np.exp((np.pi*2*np.pi-0.5)*0.01)-1):.0f}%")

print()
print("# W  winding sector on a U(1) ring, bond weight exp(K cos d), d in (-pi, pi]")
print("#    dF(n) = -ln[Z_{W=1}/Z_{W=0}]; Villain/Gaussian estimate 2 pi^2 K_eff / n")
M = 1024                                  # grid points per 2 pi
d = (np.arange(M) - M//2) * (2*np.pi/M)   # (-pi, pi)
from scipy.special import i0, i1
for K in (2.0, 8.0):
    p = np.exp(K*(np.cos(d) - 1)); p /= p.sum()
    keff = 1/(-2*np.log(i1(K)/i0(K)))     # matches 1D XY spin correlation length / 2
    xi1d = 2*keff
    print(f" K={K:4.1f}  xi_1D = {xi1d:6.2f}   (pi^2 xi = {np.pi**2*xi1d:6.1f})")
    print("    n     dF(n)    2pi^2 Keff/n")
    for n in (4, 8, 16, 32, 64, 128, 256):
        L = n*M + M                       # support of the n-fold sum, padded
        f = np.fft.rfft(p, L)
        dens = np.fft.irfft(f**n, L)      # density of the sum, origin at index n*(M//2)
        o = n*(M//2)
        z0, z1 = dens[o], dens[o + M]     # sum = 0 and sum = 2 pi
        dF = -np.log(z1/z0) if z1 > 0 else float('inf')
        print(f"  {n:4d}  {dF:8.4f}   {2*np.pi**2*keff/n:8.4f}")
print("W  verdict: dF falls monotonically with n; the wrap is distinguishable (dF > 1)")
print("   only below ~ pi^2 xi. A generated xi is a ceiling on the lock, not a line.")
