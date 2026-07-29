#!/usr/bin/env python3
"""
t5_flux.py — v89 prestress theory, part 5: FLUX AS PRESTRESS.

Wound cable cycles carry one-way circulation sigma (energy/time through
each link). This script derives, with laws_V2g numbers:

  1. the steady circulation + the FLIGHT-LOAD FIXED POINT: flight is load
     (pass 0/1), so circulation re-pitches the ring — quantified seeding
     correction, and the census-invisible flight inventory;
  2. vertex divergence-freedom (self-stress identification);
  3. the linearized phase dynamics: entrainment damping + GYROSCOPIC
     precession — does circulation stiffen? (answer: no, it precesses);
  4. the unwind barrier: restoring-force profile ~ gate(psi)*psi, the
     desert, capture range, slip kinematics;
  5. the vacuum-bleed leak floor formula, per structure class, vs measured.

Kernel facts used (verified): pass 0 flload = 0.5*sum incident lem;
pass 2 want chain; pass 4 deposit + deposit-time entrainment
(mix = f/(f+mob+lock_floor*cap)); pass 5 delivery each completed d/C
cycle + delivery-time entrainment (mix = take/(take+mob_prev+lockf));
sawtooth slot inventory; pass 6 clock advance at w2e.
"""
import numpy as np
from math import pi, cos, sin, sqrt

# laws_V2g
w2, qd, cap = 2.9, 1.2, 2.5
e_s0, s_pull, r0, dmin = 1.0, 0.15, 0.85, 1.0
kd, km = 1.2, 2.0
gm, G_r, rough_k = 0.10, 0.5, 0.35
mob_floor, lock_floor, klock = 0.004, 0.005, 1.0
p_gate = 8.0
KDM = kd * km
lockf = lock_floor * cap

def gate(psi): return (0.5 * (1 + np.cos(psi))) ** p_gate
def omega(x): return w2 / (1 + qd * x)
def x_of(w): return (w2 / w - 1) / qd
def r_live(x): return r0 * (e_s0 - s_pull * x * cap / (1 + s_pull)) ** (1 / 3)
def geo(d, ri, rj):
    if d >= ri + rj: return 0.0
    t = d * d - rj * rj + ri * ri
    a2 = (4 * d * d * ri * ri - t * t) / (4 * d * d)
    A = pi * a2 if a2 > 0 else 0.0
    return (A / (pi * r0 ** 2)) * (2 * r0 / d)

print("=" * 72)
print("1. STEADY CIRCULATION + THE FLIGHT-LOAD FIXED POINT")
print("=" * 72)
print("""  On-lock ring: per-direction throughput (pass-2 chain, gates 1/g_b):
     sigma_f = kd*km*geo*gpl*res(0)*head*Em,   sigma_b = g_b*sigma_f
  Slot inventory is a sawtooth (fills for d/C, delivers whole): mean
  lem = sigma*d/(2C) per slot; flload = 0.5*(2 links)*(sigma_f+sigma_b)
  *d/(2C) = (sigma_f+sigma_b)*d/(2C). PITCH LOAD x = (Em+flload)/cap:
  the seeder puts x*cap ALL in Em => a fresh ring wakes up FLAT of its
  rung by flload/cap and sheds to the rung (F1's one-sided ratchet,
  measured as the early fast 'leak'). Fixed point per class:""")

classes = [
    # name, N, x_rung, d, phi, g_b
    ("ring6  (m=3, pi-rung)", 6, 0.128, 1.25, pi, 1.0),
    ("comp12 (m=5, wound)",   12, 0.32054, 1.25, 5*pi/6, gate(5*pi/3)),
    ("comp6  (m=2, one-way)", 6, 0.43593, 1.10, 2*pi/3, gate(4*pi/3)),
    ("diode8 (m=2, phi=pi/2)", 8, 0.75, 1.047, pi/2, 0.0),
    ("prism3 tube belts",     3, 0.363, 1.037, 2*pi/3, gate(4*pi/3)),
]
print(f"  {'class':24s} {'gpl':>6s} {'sig_f':>7s} {'flload':>7s} {'Em_eq':>6s}"
      f" {'flight%':>7s} {'seed_x_corr':>11s}")
for name, N, xr, d, phi, gb in classes:
    gpl = cos(2 * pi / N) ** 4
    r = r_live(xr)
    g = geo(d, r, r)
    # solve Em + (1+gb)*KDM*g*gpl*Em*(1-Em/cap)*d/2 = xr*cap
    coef = (1 + gb) * KDM * g * gpl * d / 2
    Em = xr * cap
    for _ in range(200):
        Em = xr * cap / (1 + coef * (1 - Em / cap))
    fl = xr * cap - Em
    sf = KDM * g * gpl * Em * (1 - Em / cap)
    print(f"  {name:24s} {gpl:6.4f} {sf:7.4f} {fl:7.4f} {Em:6.4f}"
          f" {100*fl/(xr*cap):6.1f}% {Em/cap:11.4f}")
print("""  => CENSUS CAVEAT QUANTIFIED: 3-27% of a wound ring's rung mass lives
  IN FLIGHT (lem), invisible to the Em-census — part of the measured early
  M drop is inventory transfer, not loss (M0's -0.039 residual, tET check).
  SEEDING RULE: seed Em = x*cap - flload_eq (last column * cap), or accept
  the F1 shed. NOTE ring6's tiny gpl (0.0625) makes its internal exchange
  ~9x weaker than comp12's — slower re-temper, weaker entrainment grip.""")

print()
print("=" * 72)
print("2. DIVERGENCE-FREEDOM — the self-stress identification")
print("=" * 72)
print("""  Every kernel move is exactly paired (deposit debits source at resolve;
  delivery credits receiver; roughness/back_s split at the receiver), so
  per cell:  dEm/dt = in - out + cond - rough_local - evap.
  STEADY STATE (dEm/dt = 0) with no conversions (uniform-omega interior:
  det = 0 => rough = 0; Em < cap => evap = 0; Ee ~ 0 => cond = 0) gives
  in = out at every voice: the dense flux sigma is DIVERGENCE-FREE, i.e.
  sigma lies in ker(B) — the cycle space. That is precisely the FDM
  self-stress: member forces in equilibrium with zero external load.
  It is AUTOMATIC at steady state (conservation), NOT instantaneously.
  A leaking structure's divergence = its conversion/bleed rate — the
  self-stress decays with the leak (G4's moral for space flux, here for
  dense flux). C5's 'internal cycling at constant mass' = a nonzero
  self-stress with zero divergence everywhere INCLUDING the rim.""")

print()
print("=" * 72)
print("3. LINEARIZED PHASE DYNAMICS — damping vs gyroscopic precession")
print("=" * 72)
print("""  Perturb clocks u_k = dtheta_k around lock. Each arriving parcel kicks
  the receiver by klock*mix*(u_send - u_recv) (passes 4c and 5, receiver
  only). Per unit time (deposit-time + delivery-time):
     R = klock*sigma*[1/(Em+lockf) + 1/(Em+lockf+sigma*d/C)] ~ 2*klock*sigma/Em
  Ring:  du_k/dt = R_f (u_{k-1}-u_k) + R_b (u_{k+1}-u_k),  R_b/R_f = g_b.
  Modes u ~ e^{i 2pi j k/N}:
     lambda_j = -(R_f+R_b)(1 - cos(2pi j/N))  +/- i (R_f-R_b) sin(2pi j/N)
  DECAY is set by the SUM (total exchange); NET circulation (R_f - R_b)
  contributes ONLY the imaginary part: perturbations precess around the
  loop while decaying — gyroscopic coupling, derived. Circulation does
  NOT add restoring stiffness. Numbers (slowest mode j=1):""")
for name, N, xr, d, phi, gb in classes:
    gpl = cos(2 * pi / N) ** 4
    r = r_live(xr); g = geo(d, r, r)
    coef = (1 + gb) * KDM * g * gpl * d / 2
    Em = xr * cap
    for _ in range(200): Em = xr * cap / (1 + coef * (1 - Em / cap))
    sf = KDM * g * gpl * Em * (1 - Em / cap); sb = gb * sf
    Rf = klock * sf * (1 / (Em + lockf) + 1 / (Em + lockf + sf * d))
    Rb = klock * sb * (1 / (Em + lockf) + 1 / (Em + lockf + sb * d))
    th1 = 2 * pi / N
    lam = (Rf + Rb) * (1 - cos(th1))
    prec = abs(Rf - Rb) * sin(th1)
    print(f"  {name:24s} R_f={Rf:6.4f} R_b={Rb:6.4f}/t.u."
          f"  decay_1={lam:7.5f}  precession_1={prec:7.5f} rad/t.u.")
print("""  => relock time after a kick ~ 1/decay_1: tens to thousands of t.u. for
  the slowest mode — consistent with rings surviving perturbation but
  re-tempering slowly. ring6's grip is an order weaker than comp12's.""")

print()
print("=" * 72)
print("4. THE UNWIND BARRIER — desert, capture, slip kinematics")
print("=" * 72)
psis = np.linspace(0, pi, 2001)
pull = gate(psis) * psis
istar = np.argmax(pull)
print(f"""  One link's phase error against detune: dpsi/dt = dw - R0*gate(psi)*psi,
  R0 = 2*klock*sigma0/Em (restoring rides the FLUX, which rides the gate).
  max restoring at psi* = {psis[istar]:.3f} rad, gate(psi*) = {gate(psis[istar]):.3f},
  capture: |dw| <= R0 * {pull[istar]:.4f}.""")
for name, N, xr, d, phi, gb in classes[:3]:
    gpl = cos(2 * pi / N) ** 4
    r = r_live(xr); g = geo(d, r, r)
    coef = (1 + gb) * KDM * g * gpl * d / 2
    Em = xr * cap
    for _ in range(200): Em = xr * cap / (1 + coef * (1 - Em / cap))
    sf = KDM * g * gpl * Em * (1 - Em / cap)
    R0 = 2 * klock * sf / Em
    print(f"  {name:24s} capture |dw| <= {R0*pull[istar]:.4f} rad/t.u."
          f"  (clipped by comb unlock at 2*Gamma_m = {2*gm}: eff "
          f"{min(R0*pull[istar], 2*gm):.4f})")
print(f"""  THE DESERT: restoring ~ 0 wherever gate ~ 0: gate(1.5)={gate(1.5):.2e},
  gate(2.618)={gate(2.618):.2e}, gate(pi)=0 exactly. Between psi ~ 1.2 and
  2pi-1.2 a slipping link is neither pushed back NOR pushed on — the slip
  has no motor. To unwind (m -> m-1) a link must CROSS the desert: it
  needs a sustained one-sided detune dw for a time T ~ (2pi-2.4)/dw; at
  the post-unlock drift scale dw ~ 2*Gamma_m = 0.2 that is T ~ {(2*pi-2.4)/0.2:.0f} t.u.
  of coherent detune — but a UNIFORM ring leaks uniformly: pairwise det
  stays 0 and no motor ever appears. VERDICT: winding protection is the
  gate desert + detune-symmetry (B1's argument, now quantitative), NOT
  gyroscopic stiffness. Slip risk concentrates where symmetry breaks:
  unequal voices (census mass steps), rim contacts, strike transients.""")

print()
print("=" * 72)
print("5. THE LEAK FLOOR — vacuum bleed per voice, vs measured")
print("=" * 72)
print("""  A rung voice at pitch w(x) bleeds through EVERY candidate link to the
  vacuum choir (w2 = 2.9). Per rim link (pass-2 chain):
     lam_1 = kd*km*geo*gpl_rim*res(Dw)*<g>*head_vac*sqrt(Em*mob_floor*cap)
  Dw = w2*q*x/(1+q*x) (detune to vacuum), res = gm^2/(gm^2+Dw^2),
  head_vac ~ 1, mob-floor sqrt because the vacuum partner is silent.
  Fraction R(Dw)*rough_k of delivered bleed RADIATES (D->F); the rest is
  dense crumbs (partially recondensable). Bands: <g> in [0.196 random,
  1.0 recruited], gpl_rim in [0.07 random, 0.3 partial], n_rim ~ 8-12.""")
print(f"  {'class':24s} {'Dw':>6s} {'res':>7s} {'lam/voice (band)':>22s} {'rough%':>6s}")
meas = {"ring6  (m=3, pi-rung)": 1.9e-4, "comp12 (m=5, wound)": None}
for name, N, xr, d, phi, gb in classes:
    w = omega(xr)
    Dw = w2 - w
    res = gm ** 2 / (gm ** 2 + Dw ** 2)
    Em0 = xr * cap
    r = r_live(xr)
    gg = geo(1.25, r, r0)         # typical rim link d~1.25, vacuum radius r0
    lam_lo = 10 * KDM * gg * 0.07 * res * 0.196 * sqrt(Em0 * mob_floor * cap)
    lam_hi = 10 * KDM * gg * 0.30 * res * 1.0 * sqrt(Em0 * mob_floor * cap)
    Rr = 2 * Dw * G_r / (Dw ** 2 + G_r ** 2) * rough_k
    print(f"  {name:24s} {Dw:6.3f} {res:7.4f}   [{lam_lo:8.2e}, {lam_hi:8.2e}]"
          f" {100*Rr:5.1f}%")
print("""  Measured anchors: ring6 early leak 1.1e-3/t.u. total (M 1.92, 6 voices)
  = 1.9e-4/voice — a factor ~3.7 BELOW the band's low edge: the formula is
  an UPPER-BOUND family (it overestimates rim recruitment <g>, ignores
  crumb-neighbor headroom loss and the measured recondensation reclaim).
  comp12 late (~3e-4/t.u. absolute at M=0.59) matches the low edge within
  ~2x once Em_eq is taken at its late value (sqrt(Em) scaling). The
  load-bearing claim is the SCALING LAW, not the prefactor:
  lam ~ sqrt(x)*res(Dw(x)) falls monotonically with load —""")
for name, N, xr, d, phi, gb in classes:
    w = omega(xr); Dw = w2 - w
    print(f"    {name:24s} sqrt(x)*res = {sqrt(xr)*gm**2/(gm**2+Dw**2):.5f}")
print("""  ring6 : comp12 : diode8 = 6.6 : 1 : 0.5 per unit mass — heavier rungs
  are QUIETER against the vacuum in the resonance factor alone, on top of
  being farther from the skirt. The floor never reaches zero while any
  candidate link to the choir exists: mob_floor*cap = 0.01 guarantees a
  sqrt(Em*0.01) trickle — the design goal is to make res(Dw)*<g> small
  (heavy, and phase-decoupled rims) and let recondensation recycle crumbs.
  FLUX-MOMENT INSTRUMENT (per-ring circulation signature, MASS open item):
  loop-tangential sum of lem imbalance = N*(sigma_f-sigma_b)*d/(2C):""")
for name, N, xr, d, phi, gb in classes:
    gpl = cos(2 * pi / N) ** 4
    r = r_live(xr); g = geo(d, r, r)
    coef = (1 + gb) * KDM * g * gpl * d / 2
    Em = xr * cap
    for _ in range(200): Em = xr * cap / (1 + coef * (1 - Em / cap))
    sf = KDM * g * gpl * Em * (1 - Em / cap)
    print(f"    {name:24s} predicted loop flux moment = {N*(sf-gb*sf)*d/2:.4f}")
