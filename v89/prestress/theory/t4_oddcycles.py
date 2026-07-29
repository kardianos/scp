#!/usr/bin/env python3
"""
t4_oddcycles.py — v89 prestress theory, part 4: BIPARTITENESS, GENERALIZED.

(a) Can interval rungs (the comb, p*q <= comb_limit=6) close an ODD cycle
    of BOTH-GATE (strut-grade) links where 1:1 antiphase cannot?
    Key identity: an interval link on its rung carries a HALF-PERIOD drop.
    Lock q w_i = p w_j  +  rung (q w_i + p w_j) d = 2 pi m
    =>  phi = q w_i d = pi m   for EVERY interval — the pi is universal.
    Phase closure around a cycle then reads sum_k y_k m_k == 0 (mod 2),
    y = the integer left-null vector of the cycle's M_pq. Exhaustive scan
    of ratio assignments and rung integers against the kernel windows.

(b) Non-uniform omega without a rational lock: psi drifts at det; can the
    periodic gate opening synchronize with the action atoms eps(w) =
    A0*w/2pi to give a stroboscopic lock? Verdict from kernel structure +
    numbers.
"""
import numpy as np
from math import pi, gcd
from fractions import Fraction
from itertools import product

# laws_V2g
w2, q_det, cap = 2.9, 1.2, 2.5
e_s0, s_pull, r0, dmin = 1.0, 0.15, 0.85, 1.0
gm, G_r, rough_k = 0.10, 0.5, 0.35
x_skirt = 2 * gm / (q_det * (w2 - 2 * gm))

def omega(x): return w2 / (1 + q_det * x)
def x_of(w): return (w2 / w - 1) / q_det
def r_live(x): return r0 * (e_s0 - s_pull * x * cap / (1 + s_pull)) ** (1 / 3)

W_LO_HARD, W_HI = omega(1.0), omega(x_skirt)      # x=1.0 absolute cap
W_LO_PRACT = omega(0.9)                            # practical x<=0.9

combs = [(p, qq) for p in range(1, 7) for qq in range(1, 7)
         if p * qq <= 6 and gcd(p, qq) == 1]
ratios = sorted({Fraction(qq, p) for p, qq in combs})
print(f"comb ratios (w_j/w_i = q/p, pq<=6): {[str(r) for r in ratios]}")
print(f"pitch window: hard [{W_LO_HARD:.4f}, {W_HI:.4f}] (x<=1.0),"
      f" practical [{W_LO_PRACT:.4f}, {W_HI:.4f}] (x<=0.9)")
print(f"max in-window pitch ratio: {W_HI/W_LO_HARD:.4f} (hard), {W_HI/W_LO_PRACT:.4f} (practical)")
print()

def leftnull(pq_list):
    """integer left-null vector of the cycle rows q_k th_k - p_k th_(k+1)"""
    N = len(pq_list)
    y = [Fraction(1)]
    for k in range(N - 1):
        p_k, q_k = pq_list[k]
        p_n, q_n = pq_list[k + 1]
        # -y_k p_k + y_(k+1) q_(k+1) = 0
        y.append(y[k] * p_k / q_n)
    # consistency on the closing edge is guaranteed by prod(ratios)=1
    from math import lcm
    L = 1
    for f in y: L = L * f.denominator // gcd(L, f.denominator)
    yi = [int(f * L) for f in y]
    g = 0
    for v in yi: g = gcd(g, v)
    return [v // g for v in yi]

def scan_cycle(N, practical=True, mmax=4, verbose_top=3):
    W_LO = W_LO_PRACT if practical else W_LO_HARD
    found, near = [], []
    for assign in product(range(len(ratios)), repeat=N):
        rs = [ratios[a] for a in assign]
        prod_r = Fraction(1)
        for r in rs: prod_r *= r
        if prod_r != 1: continue
        if all(r == 1 for r in rs): continue          # pure 1:1 handled classically
        # frequencies relative to w0
        rel = [Fraction(1)]
        for r in rs[:-1]: rel.append(rel[-1] * r)
        mx, mn = max(rel), min(rel)
        # w0 window from pitch bounds
        w0lo, w0hi = W_LO / float(mn), W_HI / float(mx)
        if w0lo >= w0hi: continue
        pqs = []
        ok = True
        for k in range(N):
            # lock q w_k = p w_(k+1): ratio r = q/p reduced
            r = rs[k]
            qq, p = r.numerator, r.denominator
            if p * qq > 6: ok = False; break
            pqs.append((p, qq))
        if not ok: continue
        y = leftnull(pqs)
        # scan w0 grid; for each link pick smallest feasible m; then parity
        for w0 in np.linspace(w0lo, w0hi, 60):
            ws = [w0 * float(f) for f in rel]
            ds, ms, feas = [], [], True
            for k in range(N):
                p, qq = pqs[k]
                wsum = qq * ws[k] + p * ws[(k + 1) % N]
                dmax = r_live(x_of(ws[k])) + r_live(x_of(ws[(k + 1) % N]))
                mfeas = [m for m in range(1, mmax + 1)
                         if dmin <= 2 * pi * m / wsum <= dmax]
                if not mfeas: feas = False; break
                ds.append([2 * pi * m / wsum for m in mfeas]); ms.append(mfeas)
            if not feas: continue
            # parity: exists choice m_k with sum y_k m_k even?
            for choice in product(*ms):
                if sum(yk * mk for yk, mk in zip(y, choice)) % 2 == 0:
                    found.append((rs, w0, choice,
                                  [2 * pi * m / (pqs[k][1] * ws[k] + pqs[k][0] * ws[(k + 1) % N])
                                   for k, m in enumerate(choice)],
                                  [x_of(w) for w in ws]))
                    break
            else:
                near.append((rs, w0, [m[0] for m in ms], y))
                continue
            break
    return found, near

for N in (3, 5):
    for practical, label in ((True, "practical x<=0.9"), (False, "hard x<=1.0")):
        found, near = scan_cycle(N, practical)
        print(f"ODD {N}-CYCLE with interval rungs [{label}]:"
              f"  feasible assignments: {len(found)}")
        if found:
            seen = set()
            for rs, w0, ms, ds, xs in found:
                key = tuple(str(r) for r in rs)
                if key in seen: continue
                seen.add(key)
                print(f"   ratios {[str(r) for r in rs]} w0={w0:.3f} m={ms}"
                      f" d={[round(d,3) for d in ds]} x={[round(x,3) for x in xs]}")
                if len(seen) >= 6: break
        if near and not found:
            rs, w0, ms, y = near[0]
            print(f"   (window-feasible but PARITY-blocked example: ratios"
                  f" {[str(r) for r in rs]}, y={y}, m={ms})")
print("""
RESULT — the brief's expectation is REFUTED in one direction: interval
rungs CAN beat 1:1 bipartiteness. The loophole is even-m rungs: an
interval link's rung sum (q w_i + p w_j) = 2 q w_i is TWICE the 1:1 sum
at matched pitch, so m=2 lands at in-window length d = pi*m/(2 w_i) while
its drop phi = pi*m == 0 (mod 2pi) — a 'ferromagnetic' bond in the
antiphase Ising picture. Parity sum y.m mod 2 then closes odd cycles.
1:1-ONLY strut networks remain strictly bipartite (m=1 forced).""")

print("=" * 72)
print("THE FIFTH-TRIANGLE — smallest odd strut-grade cycle (design card)")
print("=" * 72)
print("  v0 --(3:2, m=2)-- v1 --(2:3, m=2)-- v2 --(1:1, m=1)-- v0")
print("  all three links BOTH-GATE on-rung; all lengths equal d = pi/w0")
print("  w0 window   x0=x2      w1=2w0/3  x1        d=pi/w0  contact margin (nom)")
for w0 in (2.091, 2.2, 2.35, 2.5, 2.7):
    w1 = 2 * w0 / 3
    x0, x1 = x_of(w0), x_of(w1)
    if x1 > 1.0 or x0 < x_skirt: continue
    d = pi / w0
    marg = r_live(x0) + r_live(x1) - d
    tag = "" if x1 <= 0.9 else "  (x1 beyond practical)"
    print(f"  w0={w0:.3f}   x0={x0:.3f}   w1={w1:.3f}  x1={x1:.3f}   {d:.4f}   {marg:+.4f}{tag}")
print(f"""  COST (fragility, quantified): the interval links run at comb amplitude
  1/(pq) = 1/6 and acceptance Gamma_m/(pq) = {gm/6:.4f} rad/t.u. (vs {gm} for
  1:1) — 6x weaker exchange, 6x narrower tongue; the measured interval
  lifetime hierarchy (E9: fifth ~20 t.u. by detuning drift) applies to the
  3:2 legs. The heavy voice (x1 ~ 0.75-0.9) has headroom 0.1-0.25 and one
  fluctuation of margin to the cap. PREDICTION: the fifth-triangle forms
  and holds while its 3:2 legs stay in the Gamma/6 tongue, then dies by
  interval drift — lifetime between the naked fifth (~20 t.u.) and a
  closed unison ring (>10^3), closer to the fifth. Worth one fleet probe;
  not a mass-particle candidate.""")

print("=" * 72)
print("(b) BEAT-CONSONANCE / STROBOSCOPIC LOCK — verdict")
print("=" * 72)
print("""KERNEL FACTS (verified):
  1. Dense-to-dense transport is CONTINUOUS — 'transport within a mode is
     continuous; the action atom lives at mode boundaries' (pass 4, the
     comment block before the deposit loop; atoms_fire is called ONLY in
     pass 5 roughness and pass 6 beat conversion/evaporation).
  2. A drifting 1:1 link (uniform d, det = dw) has psi(t) = dw*t: the gate
     opens periodically with duty <gate> = C(16,8)/2^16 = 0.19638.
  3. Deliveries complete every d/C regardless of gate phase; the delivered
     amount is the integral of gated deposits — no event for the beat to
     synchronize onto. The delivery-time entrainment kick IS beat-sampled:
     stroboscopic freezing of the sampled phase needs dw*d/C = 2pi*k.""")
d = 1.25
dw_strobe = 2 * pi / d
print(f"  stroboscopic commensurability needs dw = 2pi*C/d = {dw_strobe:.3f} rad/t.u.")
print(f"  comb acceptance (1:1): res(dw) = gm^2/(gm^2+dw^2) = "
      f"{gm**2/(gm**2+dw_strobe**2):.2e}  (vs 1 on-tune)")
print(f"  measured unlock boundary (S2): |det| ~ 2*Gamma_m = {2*gm}")
print(f"  => the first strobe point sits {dw_strobe/(2*gm):.0f}x beyond the unlock")
print(f"     boundary with 4 orders of magnitude less exchange — nothing to lock.")
print("""  4. The action atoms quantize AMOUNTS of conversion, not timing; the
     credit pool lapses at 2 atoms (atoms_fire) — no phase is stored, so
     there is no atom-clock for the gate window to entrain.
  5. What a detuned link ACTUALLY does at delivery: radiate. Rough fraction
     R(det) = 2|det|*G_r/(det^2+G_r^2) * rough_k of every delivered parcel:""")
for dw in (0.05, 0.1, 0.2, 0.4, 0.8, 1.6):
    R = 2 * dw * G_r / (dw ** 2 + G_r ** 2) * rough_k
    print(f"     det={dw:4.2f}: rough fraction {R:.3f}"
          + ("   <- M0: this channel is 52% of blob loss" if dw == 0.8 else ""))
print("""VERDICT: refuted. Beat-consonance with quantization CANNOT rescue
non-uniform omega: (i) transport carries no atoms to synchronize; (ii) the
strobe commensurability point lies ~20x outside the comb acceptance, where
exchange is already dead; (iii) sustained det is precisely the roughness
channel — the measured particle-killer, not a lock. The ONLY discrete
beat-locks in the kernel are the comb's own p:q rungs (rate 1/(pq),
acceptance Gamma_m/(pq)) — already in the law table, already measured to
die young (E9: unison >> fifth > octave).""")
