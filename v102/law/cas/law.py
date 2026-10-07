#!/usr/bin/env python3
"""Part L — what the table implies for t, Delta (and g) on the fiber.

All pitches are in the kernel's rate units (the field self-rotation is
z -> exp(-i w1e dt) z, freecell.c pass F; the hop angle is
field_J*w_hat*dt, freecell.c:2436).  Energies are hbar_eff * rate with
hbar_eff = A0/2pi; every ratio below is independent of A0.

Table (laws_V2g verbatim): w1 = 1.65 (field rung), w2 = 2.9 (dense rung),
q_detune = 1.2, cap = 2.5, quant_A0 = 1.15, field_J = 1.8.
Pitch under load: w/(1 + q x), x = load/cap (freecell.c:1579).
"""
import numpy as np
from scipy.optimize import brentq

w1, w2, q, cap, A0, fJ = 1.65, 2.9, 1.2, 2.5, 1.15, 1.8
hbar = A0/(2*np.pi)

# ---- t: the R1 lift of the kernel hop, hop rate field_J * w_hat.
#      w_hat = w/sqrt(s_i s_j); on a degree-2 ring with equal links w_hat = 1/2.
t_ring = fJ*0.5
print("# t (reading T-a: the fiber inherits the 3-D field transport)")
print(f"  t = field_J * w_hat = {fJ} * 1/2 = {t_ring:.3f}   (rate units)")

# ---- Delta: E(space quantum) - [E_field(2,2) - E_field(1,1)]
def cell_energy(N, w=w1):
    """Rate-unit energy of N quanta sharing one cell's pitch w/(1+q x),
    with x = hbar*e/cap self-consistent (load flattens pitch)."""
    if N == 0: return 0.0
    f = lambda e: e - N*w/(1 + q*hbar*e/cap)
    return brentq(f, 0, N*w + 1)

e1, e2 = cell_energy(1), cell_energy(2)
release = 2*(e2 - e1)            # field energy freed by converting one quantum from each of two full cells
x2 = hbar*e2/cap
print("\n# Delta readings (rate units); the lock needs Delta > 0 and t/Delta < t_c ~ 0.16-0.26")
print(f"  full cell (2 quanta): e = {e2:.4f}, load x = {x2:.4f}, pitch = {w1/(1+q*x2):.4f};  pair release = {release:.4f}")
rows = []
# R-a: space quantum at the dense rung, no load anywhere (bare rungs)
rows.append(("R-a  space quantum at bare dense rung w2; field at bare w1", w2 - 2*w1))
# R-b: space quantum at bare w2; field quanta at their loaded pitch (self-consistent)
rows.append(("R-b  space at bare w2; field loaded at the R6 ceiling", w2 - release))
# R-b': space quantum self-loaded on its link (x_s = hbar*Omega/cap)
Om = brentq(lambda o: o - w2/(1 + q*hbar*o/cap), 0, w2)
rows.append(("R-b' space at self-loaded w2; field loaded", Om - release))
for name, D in rows:
    tag = "no lock (Delta <= 0: converted ground state)" if D <= 0 else (
          f"t/Delta = {t_ring/D:.2f}  ->  {'inside' if t_ring/D < 0.16 else 'OUTSIDE'} the lock phase")
    print(f"  {name:58s} Delta = {D:+.4f}   {tag}")
print("  R-c  space has no pitch (v89 sec 2; kernel Es is a continuum): a pair always finds a")
print("       resonant space state -> Delta = 0 and real decay at the golden-rule rate; no stationary lock")

# ---- what would t have to be?
print("\n# For the lock to exist under R-b, the fiber transport must satisfy t < t_c * Delta:")
D = w2 - release
for tc in (0.16, 0.19, 0.26):
    print(f"  t_c = {tc:.2f}: t < {tc*D:.4f}  i.e. field_J*w_hat must be {t_ring/(tc*D):.0f}x smaller than the 3-D value")

# ---- g
print("\n# g: the table's conversion is a rate law (f_conv = 0.25 of the excess per step, k_dep = 1.2,")
print("#    atoms at the source's pitch, credit lapsing at two atoms; freecell.c:2740).  It moves")
print("#    energy, not amplitude: there is no coherent pair-to-space coupling constant to read.")

# ---- reading T-b (added after the first computation; source v98 PRINCIPLES sec 1,
#      "locked field cannot travel as itself", already discussed in gates/WEDGE.md W4):
#      no transport inside a locked fiber, t = 0.
print("\n# T-b (added after the first computation): v98 sec 1 read as t = 0 inside a lock")
for name, D in rows:
    print(f"  {name:58s} Delta = {D:+.4f}   {'inside the lock phase (t/Delta = 0)' if D > 0 else 'no lock (Delta <= 0)'}")
print("  R-c  no pitch: resonant decay, no stationary lock, whatever t is")
print("\n# Which Delta reading matches the current representation?  v99 sec 1 retired the dense")
print("#  fluid; under pattern_only (v99 KERNEL.md) conversion is field <-> space, and space (Es)")
print("#  carries no pitch: the kernel's load x uses Ee/Em + flight, never Es.  That is R-c.")
print("#  R-a/R-b/R-b' borrow the retired dense rung w2 for the space quantum.")
