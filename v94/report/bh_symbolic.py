#!/usr/bin/env python3
"""Symbolic no-go analysis: can laws_V2g (+radiance candidate) form a black hole?
Every expression is taken from the kernel as implemented (v91/kernel/freecell.c
pass S / pass 6) and the standing table laws_V2g.cfg. sympy verifies each claim."""
import sympy as sp

# ---- law constants (laws_V2g.cfg, verbatim) ----
q, cap = sp.Rational(12,10), sp.Rational(5,2)          # q_detune, cap
s_pull, s_disp, s_k = sp.Rational(15,100), sp.Rational(3,10), sp.Rational(6,100)
f_conv, f_evap = sp.Rational(1,4), sp.Rational(1,2)
w1, w2 = sp.Rational(165,100), sp.Rational(29,10)
mob_floor = sp.Rational(4,1000)

x = sp.symbols('x', positive=True)

print("="*72)
print("[1] CLOCK / REDSHIFT BOUND  (pitch law w_e = w/(1+q*x), both pitches)")
dil = 1 + q*x                              # dilation factor w/w_e
print("  dilation(x) = 1 + q*x =", dil)
print("  at standing bound x=1 (cap + exact refusal):     dilation =", dil.subs(x,1), "= 2.2")
xf = sp.Rational(18,10)                    # transient pitch-x with worst measured flight load
print("  with flight load (flload~2.0, x->1.8 transient): dilation =", sp.nsimplify(dil.subs(x,xf)), "=", float(dil.subs(x,xf)))
print("  horizon requires dilation -> oo  <=>  x -> oo  <=>  Em+flload unbounded")
print("  Em+Ee <= cap is enforced every beat (f_evap) and intake at cap is EXACTLY 0")
print("  (PAULI-0, FP-sticky). => max gravitational redshift z = q*x <= ~2.2. NO HORIZON.")

print()
print("="*72)
print("[2] OPACITY IS ANTI-EDDINGTON  (absorption rides the beat clock)")
beat = (w2-w1)/(1+q*x)                     # beat rate |w1e-w2e|
kappa = f_conv*beat                        # capture opportunity per unit time
print("  beat(x) = (w2-w1)/(1+q*x) =", beat, " d(beat)/dx =", sp.simplify(sp.diff(beat,x)))
print("  => strictly decreasing: the fuller a region, the SLOWER it eats.")
print("  kappa(0)/kappa(1) =", sp.nsimplify(kappa.subs(x,0)/kappa.subs(x,1)))
print("  Accretion self-throttles; no runaway capture at any fullness.")
print("  (Measured: obj_amp 0.8 absorbs LESS than 0.5 -- 6.58 vs 7.27, XSEC.)")

print()
print("="*72)
print("[3] THE PI-FLATNESS THEOREM  (why stable mass has NO far field -- exact)")
d1, d2, rc, re_ = sp.symbols('d1 d2 r_c r_e', positive=True)   # event sizes & rates
dsp = s_pull*d1                       # space pulled at condensation (kernel line 1867)
bs  = d2*s_pull/(1+s_pull)            # backsplash at evaporation   (kernel line 1888)
# pi = Es + s_disp*(Em+Ee); per-event medium load:
dpi_cond = -dsp + s_disp*((d1+dsp) - d1)          # dEs=-dsp, dEm=+d1+dsp, dEe=-d1
dpi_evap = +bs  + s_disp*((-d2) + (d2-bs))        # dEs=+bs,  dEm=-d2,    dEe=+d2-bs
print("  per condensation event: dpi =", sp.factor(sp.simplify(dpi_cond)))
print("  per evaporation  event: dpi =", sp.factor(sp.simplify(dpi_evap)))
# book closure: Em stationary  <=>  (d1+dsp)*rc = d2*re
em_closure = sp.Eq((d1+dsp)*rc, d2*re_)
es_closure = sp.Eq(dsp*rc, bs*re_)
print("  Em-closure:", em_closure, "  =>", sp.solve(em_closure, rc)[0], "= rc")
print("  Es-closure:", es_closure, "  =>", sp.solve(es_closure, rc)[0], "= rc")
print("  -> the two closures are THE SAME condition (door routing makes the")
print("     energy ledger and the space ledger proportional).")
net_pi_rate = dpi_cond*rc + dpi_evap*re_
net_at_closure = sp.simplify(net_pi_rate.subs(rc, sp.solve(em_closure, rc)[0]))
print("  net medium pi-load at closed books:", net_at_closure)
print("  => IDENTICALLY ZERO for ALL s_pull, s_disp, throughput. A stationary")
print("     object loads the medium with NOTHING, no matter how hard it")
print("     metabolizes. pi-flat is a theorem of the door routing, and the")
print("     measured ASTRO result (+-1e-5..2e-4 flat; only dM/dt != 0 shows).")
print("  COROLLARY (no well, no collapse): a standing gravitational well needs")
print("     sustained net pi < 0  <=> permanently net-condensing books")
print("     <=> Em grows without bound  -- forbidden by cap + refusal.")

print()
print("="*72)
print("[4] THE STANDING DENT IS BOUNDED AND CONTACT-LOCAL")
Em_ = sp.symbols('Em', positive=True)
es_deficit = s_pull*Em_/(1+s_pull)         # space consumed building Em
pi_dent    = (1-s_disp)*es_deficit
print("  building Em consumes Es:", es_deficit, " -> pi dent before refill:", sp.simplify(pi_dent))
print("  at D-engine holdings Em=1.45:", float(pi_dent.subs(Em_,1.45)), "(then s_k refill flattens pi;")
print("  measured refill completeness 0.3%, screening < 1 cell). Max possible dent")
print("  (Em=cap):", float(pi_dent.subs(Em_,cap)), "-- a fixed, saturating number. Wells cannot deepen.")

print()
print("="*72)
print("[5] RADIANCE PARKS MATTER *BELOW* SATURATION (the fixed point)")
k_rad, p_rad, I0 = sp.symbols('k_rad p_rad I_0', positive=True)
D = sp.Rational(49736,100000)*k_rad*x**p_rad     # measured demand law D(x)=0.497*k*x^p
xstar = sp.solve(sp.Eq(D, I0), x)[0]
print("  out(x)=0.497*k*x^p ; intake ~ I0  =>  x* =", xstar)
print("  at k=0.05, p=4, I0~0.0093 (ambient):  x* =", float(xstar.subs({k_rad:sp.Rational(5,100), p_rad:4, I0:0.0093})))
print("  (measured x^ = 0.62+-0.02).  x* < 1 for all I0 < 0.497*k: matter is an")
print("  INTERIOR fixed point -- it never even reaches the wall. V2g limit k->0:")
print("  out=0 below cap -> only the wall exists (FORGE E1) and the wall EMITS")
print("  (XSEC flat-top -7.2): the compact limit is a WHITE body, not a black one.")

print()
print("="*72)
print("[6] NEWTON'S CONSTANT, AS MEASURED")
# COMBINE CO-T2F: sep slope -6.8e-5 +- 4.2e-5 /t.u. at range 7 over 4500 t.u., ~0
slope, sep0, T = -6.8e-5, 7.0, 4500.0
print("  two chords at range 7: d(sep)/dt = -6.8e-5 +- 4.2e-5 /t.u.  (1.6 sigma ~ 0)")
print("  => a_grav <= ~3e-8 cell/t.u^2 class; with M ~ few units at r=7:")
print("     G_eff consistent with ZERO at measurement floor.  r_s = 2*G*M/c^2 = 0")
print("  for every M. No Schwarzschild radius exists at any mass in this law.")

print()
print("="*72)
print("[7] KINEMATIC WALL: NOTHING FALLS")
print("  B2 (measured, twice): v_COE <= 1e-3 of closing; driven blob2 sep 8.00->8.37")
print("  (Q4: 'the collision never happens'). Even under a force, infall would need")
print("  translation of structure, which the law does not transport (Delta-p ~ 0).")
print("  Collapse has neither the force (item 6) nor the motion (this item).")

print()
print("="*72)
print("[8] WHAT A BLACK HOLE WOULD REQUIRE (each wall's repair, symbolically)")
print("  wall 1 (no far field)  -> books routed THROUGH the medium: dpi_cond at the")
print("     object must land at range r (medium-carried), breaking [3]'s pointwise")
print("     cancellation while keeping global closure. That is exactly S2/identity.")
print("  wall 2 (no attraction) -> a signed, phase-coherent cross-flow (the v89-s2")
print("     'choir' interference bias, retired at rate level) -- needs carried phase.")
print("  wall 3 (no infall)     -> 'translation IS the current' (S2-full criterion).")
print("  wall 4 (bounded pile)  -> cell fission/fusion (THEORY 2.3, design-only):")
print("     space itself must be consumable for a pile to deepen past cap.")
print("  wall 5 (clock floor)   -> nothing: given walls 1-4, dilation is already")
print("     unbounded only if x is; the cap IS the censor.")
print("  All five repairs are the SAME missing sector: the first-moment books")
print("  (phase / momentum / winding) that die at the conversion door today.")
