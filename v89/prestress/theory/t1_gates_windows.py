#!/usr/bin/env python3
"""
t1_gates_windows.py — v89 prestress theory, part 1.

Gate algebra, pitch/length windows, reachable rung integers, jitter
tolerance curves. Every constant from battery/laws_V2g.cfg; every law
verified against cellfab.c (gate_of ~line 482; pitch pass 1 ~line 1560;
clock advance pass 6 ~line 2216 current numbering; both-direction dense
gates pass 2 ~lines 1795-1800 old / grep th2 anchors).

All angles in radians unless noted.
"""
import numpy as np
from math import pi, cos, sin, acos, sqrt, log

# ---- laws_V2g constants (verified against the cfg byte-for-byte) ----
C      = 1.0
w2     = 2.9
q      = 1.2          # q_detune
p_gate = 8.0
cap    = 2.5
e_s0   = 1.0
es_floor = 0.05
s_pull = 0.15
dmin   = 1.0
r0     = 0.85
rjit   = 0.06
gm     = 0.10         # gamma_res_m
G_r    = 0.5          # gamma_rough
rough_k = 0.35
kd, km = 1.2, 2.0     # k_dep, k_dep_m
mob_floor = 0.004
lock_floor = 0.005

def gate(psi):
    """cellfab.c gate_of: ((1+cos psi)/2)^p_gate = cos^{2p}(psi/2)."""
    return (0.5 * (1.0 + np.cos(psi))) ** p_gate

def omega(x):        # pitch law, pass 1: w2e = w2/(1 + q*(Em+flload)/cap)
    return w2 / (1.0 + q * x)

def x_of(w):
    return (w2 / w - 1.0) / q

def Es_seeded(x):
    """Seeder bookkeeping (pair/ring/shell/net seeders, verified):
    add = x*cap/(1+s_pull); Em += add(1+s_pull) = x*cap;
    Es -= s_pull*add  =>  Es = e_s0 - s_pull*x*cap/(1+s_pull)."""
    return e_s0 - s_pull * x * cap / (1.0 + s_pull)

def r_live(x):
    """pass 1: cr = cr0 * cbrt(Es/e_s0), nominal cr0 = r0."""
    return r0 * Es_seeded(x) ** (1.0 / 3.0)

print("=" * 72)
print("A. GATE ALGEBRA — anchors vs measured/claimed numbers")
print("=" * 72)
rows = [
    ("phi=pi (m=6/N=12 unwound ring): back gate gate(2phi)=gate(2pi)", gate(2 * pi)),
    ("phi=5pi/6 (comp12 m=5): back gate gate(5pi/3)     [measured 0.100]", gate(5 * pi / 3)),
    ("phi=2pi/3 (comp6 m=2): back gate gate(4pi/3)      [B1: 1.53e-5]", gate(4 * pi / 3)),
    ("phi=pi/2 (diode point): back gate gate(pi) EXACT ZERO", gate(pi)),
    ("desert: gate(2.6184)                              [B1: ~4e-10 @ ~2.6]", gate(2.6184)),
    ("drift duty: <gate> over a full slip = C(16,8)/2^16", 12870 / 65536),
]
for label, v in rows:
    print(f"  {label:66s} = {v:.6g}")
print(f"  <gate> numeric check: {np.mean(gate(np.linspace(0, 2*pi, 20001)[:-1])):.6g}")
# B1 note claimed 'gate >= 0.9 for |psi| <= 0.16' — check
print(f"  gate(0.16) = {gate(0.16):.4f}  (B1 said 0.9 there; actual 0.9 point below)")
for gth in (0.99, 0.95, 0.9, 0.5):
    psi = 2 * acos(gth ** (1.0 / 16.0))
    print(f"  flat top: gate >= {gth}  <=>  |psi| <= {psi:.6f} rad")
print("  quadratic top: gate ~= exp(-p*psi^2/4) = exp(-2 psi^2)   "
      f"(check exp(-2*0.2292^2)={np.exp(-2*0.22915**2):.4f})")
print("  inflection (stiffness validity): psi_infl = 2*atan(1/sqrt(2p-1))"
      f" = {2*np.arctan(1/sqrt(2*p_gate-1)):.4f} rad")

print()
print("=" * 72)
print("B. WINDOWS — pitch, skirt, contact, reachable rung integers")
print("=" * 72)
# skirt: voice dissolves when within the room's 1:1 acceptance 2*Gamma_m
x_skirt = 2 * gm / (q * (w2 - 2 * gm))
print(f"  vacuum skirt x_skirt = 2G/(q(w2-2G)) = {x_skirt:.6f}   [README: 0.0617]")
print(f"  pitch window (soft edges):")
for x in (0.0617, 0.1, 0.12, 0.2, 0.32054, 0.41, 0.43593, 0.5, 0.705, 0.8, 0.9, 1.0):
    print(f"    x={x:7.4f}  omega={omega(x):.4f}  Es_seed={Es_seeded(x):.4f}"
          f"  r_live={r_live(x):.4f}  2r={2*r_live(x):.4f}")
wmin, wmax = omega(0.9), omega(x_skirt)
print(f"  practical omega window (x in [x_skirt, 0.9]): [{wmin:.4f}, {wmax:.4f}]")
print(f"  (brief said ~(1.36, 2.75): exact values above; x=1.0 hard floor gives {omega(1.0):.4f})")

# candidate-link window (build uses NOMINAL radii cr0 with 1.15 margin)
print(f"  candidate-link build cut: d < 1.15*(cr0_i+cr0_j);"
      f" nominal 2*1.15*r0 = {2*1.15*r0:.4f},"
      f" jitter range [{2*1.15*r0*(1-rjit):.4f}, {2*1.15*r0*(1+rjit):.4f}]")
print(f"  LIVE channel needs lens overlap of loaded radii: d < r_i(x)+r_j(x)"
      f"  (A<=0 links carry nothing; flight on pinched links holds)")

# strut ladder: both gates open <=> 2 w d / C = 2 pi m  <=>  d = pi m / omega
print("\n  strut ladder d = pi*m/omega over the pitch window:")
for m in (1, 2, 3):
    dlo, dhi = pi * m / wmax, pi * m / wmin
    print(f"    m={m}: d in [{dlo:.4f}, {dhi:.4f}]"
          + ("   <-- inside link window" if dlo < 2.07 else "   UNREACHABLE (beyond candidate cut 2.07)"))

# self-consistent strut window: d = pi/omega(x) must be < 2*r_live(x)
from scipy.optimize import brentq
f = lambda x: pi / omega(x) - 2 * r_live(x)
x_strut_max = brentq(f, 0.05, 0.95)
print(f"\n  strut existence (nominal radii): pi/omega(x) < 2 r(x)  =>  x <= {x_strut_max:.4f}"
      f"  (omega >= {omega(x_strut_max):.4f}, d <= {pi/omega(x_strut_max):.4f})")
fj = lambda x: pi / omega(x) - 2 * r_live(x) * (1 - rjit)
x_strut_max_j = brentq(fj, 0.05, 0.95)
print(f"  worst-case jitter (-6% both radii):            x <= {x_strut_max_j:.4f}"
      f"  (omega >= {omega(x_strut_max_j):.4f}, d <= {pi/omega(x_strut_max_j):.4f})")
print(f"  STRUT WINDOW: x in [{x_skirt:.4f}, {x_strut_max:.2f}] (design floor 2*x_skirt = {2*x_skirt:.3f}),"
      f" d = pi/omega in [{pi/wmax:.4f}, {pi/omega(x_strut_max):.3f}]")

# diode cables: back gate exactly zero at phi = pi/2 (and 3pi/2)
print("\n  cable phase drop phi = omega*d/C; back gate = gate(2 phi) exactly")
print("  diode phi=pi/2: requires d = pi/(2 omega) >= dmin  =>  omega <= pi/2"
      f" = {pi/2:.4f}  =>  x >= {x_of(pi/2):.4f}")
print(f"    diode window: x in [{x_of(pi/2):.4f}, 0.9], d = pi/(2 omega) in"
      f" [{pi/2/ (pi/2):.4f}, {pi/(2*wmin):.4f}]")
d32lo = 3 * pi / (2 * wmax)
print(f"  second diode phi=3pi/2: d = 3pi/(2 omega) >= {d32lo:.4f} vs live ceiling"
      f" {2*r_live(x_skirt):.4f} at that load  =>  UNREACHABLE")

# per-omega cable phi window
print("\n  reachable phi window vs omega (cables, both ends at x(omega)):")
for x in (0.0617, 0.12, 0.2, 0.32, 0.41, 0.5, 0.705, 0.9):
    w = omega(x); dmax = 2 * r_live(x)
    print(f"    x={x:6.4f} omega={w:.4f}: phi in [{w*dmin:.4f}, {w*dmax:.4f}]"
          f"  (d in [1, {dmax:.3f}])")

print()
print("=" * 72)
print("C. JITTER TOLERANCE — strut length error dd vs gates")
print("=" * 72)
print("""  Tree/seed lock (forward psi=0 exact): the whole per-link defect sits on
  the back gate, psi_b = 2*omega*dd  =>  g_b = cos^16(omega*dd).
  Entrained equilibrium (E6 P2: comma tempered evenly): psi_f = psi_b =
  omega*dd  =>  both gates cos^16(omega*dd/2). Tempering doubles tolerance.""")
for gth in (0.9, 0.99):
    y = acos(gth ** (1.0 / 16.0))
    print(f"  gate >= {gth}:")
    print(f"    tree lock:  |dd| <= {y:.6f}/omega     = {y/pi*100:.2f}% of strut length pi/omega")
    print(f"    tempered:   |dd| <= {2*y:.6f}/omega     = {2*y/pi*100:.2f}% of strut length")
    for w in (1.9424, 2.0944, 2.5139, 2.7002):
        print(f"      omega={w:.4f}: tree {y/w:.4f}, tempered {2*y/w:.4f}")
print("""  NOTE the tolerance as a FRACTION of the strut length is omega-independent:
  3.65%/7.29% (gate 0.9), 1.13%/2.26% (gate 0.99). Blue-noise foam placement
  error is ~10-20% of d  =>  raw seeded shells sit at gates ~0.5-0.7 (measured
  cube: mean 0.6) — annealing/entrainment-grown skins are REQUIRED, or picks
  must be mined for equal lengths (P3).""")

print("=" * 72)
print("D. PAIR TONGUE — G^2(delta), against E6 measured points")
print("=" * 72)
p2 = 2 * p_gate
G2 = lambda d: (0.5 * (1 + np.cos(d / 2))) ** p2
dhalf = 2 * acos(2 * 0.5 ** (1.0 / p2) - 1)
# CONSONANCE quotes delta_half = 2*acos(2^(1-1/(2p)) - 1)
print(f"  half-max of G^2: delta_half = {4*acos(0.5**(1.0/(2*p2))):.4f} rad"
      f" (ladder printed 0.83)")
meas = [(0.17, 0.965), (0.31, 0.899), (0.56, 0.721), (0.70, 0.607), (0.98, 0.277), (1.26, 0.036)]
print("  delta   G2_computed   G2_measured(E6)")
for d, gmeas in meas:
    print(f"  {d:5.2f}   {G2(d):11.3f}   {gmeas:.3f}")
print("  (off-tongue fails faster than the even-split formula — nonlinear unlock, as logged)")

print()
print("=" * 72)
print("E. RING GEOMETRY FACTOR — seeded plane alignment gpl(N)")
print("=" * 72)
print("""  Ring seeder sets both planes transverse to the OUTGOING chord; for a
  planar ring, n1 = z for every voice and n2 = u_out x z. On link k->k+1:
  source axial factor 1 exactly; receiver picks up the polygon turn 2pi/N:
  axj = cos^2(2pi/N), (n2_i.n2_j)^2 = cos^2(2pi/N)  =>  gpl = cos^4(2pi/N).
  If kappa_align settles the axes onto the bisector: gpl -> cos^4(pi/N).""")
for N in (4, 5, 6, 8, 12, 16, 24):
    print(f"    N={N:3d}: seeded gpl = {cos(2*pi/N)**4:.4f}   bisector gpl = {cos(pi/N)**4:.4f}")
print("  N=6 pays a 9x conductance penalty vs N=12 at seed; N=4 is dead (turn=90deg).")
