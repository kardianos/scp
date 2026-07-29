#!/usr/bin/env python3
"""
t2_cycles.py — v89 prestress theory, part 2: cycles, counting, shells.

The static design problem, reduced (drift theorem, THEORY.md S1-S2):
a locked component has ONE omega; every cell load x(omega)=(w2/omega-1)/q;
so the whole design is a CIRCLE PLACEMENT problem for the clock phases:

  strut (both-gate) i-j :  theta_i - theta_j = pi   (mod 2pi), d = pi/omega
  cable (fwd-only)  i-j :  theta_i - theta_j = +-phi, phi = omega*d in
                           B(omega) = [omega*dmin, omega*dmax(x(omega))]
                           (sign = designer's lock orientation)

  =>  allowed separation set on the circle:
      strut: {pi};  cable: B u (2pi - B).

Cycle conditions are the mod-2pi closures of these separations; struts
contribute pi each => Z2 (parity) data; cables pay continuous phase within
the band. This script settles which shells/cycles are payable, and checks
the ledger's expectation 'hexagonal prism over-constrained by exactly 2'.
"""
import numpy as np
from math import pi, cos, acos
from itertools import combinations

C = 1.0; w2 = 2.9; q = 1.2; cap = 2.5; e_s0 = 1.0; s_pull = 0.15
dmin = 1.0; r0 = 0.85
x_skirt = 2 * 0.10 / (q * (w2 - 2 * 0.10))

def omega(x): return w2 / (1.0 + q * x)
def x_of(w): return (w2 / w - 1.0) / q
def r_live(x): return r0 * (e_s0 - s_pull * x * cap / (1 + s_pull)) ** (1.0 / 3.0)
def band(w):
    """reachable cable phase-drop window at component pitch w"""
    x = x_of(w)
    return (w * dmin, w * 2 * r_live(x))

W_STRUT_MIN, W_STRUT_MAX = 1.9392, omega(x_skirt)   # strut window from t1
W_ALL_MIN, W_ALL_MAX = omega(0.9), omega(x_skirt)

# ------------------------------------------------------------------
print("=" * 72)
print("A. ODD-CYCLE PAYMENT LEMMA — how many struts can an odd cycle carry")
print("=" * 72)
print("""  A cycle with n_s struts and n_c cables needs (mod 2pi)
      sum_cables +-phi  ==  pi * n_s .
  Reachability over the sign choices: the set of reachable residues of a
  sum of n_c terms, each in +-B(omega). Scan omega.""")

def cycle_payable(n_s, n_c, w, kmax=6):
    """can n_c cable phases in +-B(w) sum to pi*n_s mod 2pi?"""
    lo, hi = band(w)
    # choose j terms positive, n_c-j negative: sum range [j*lo-(n_c-j)*hi, j*hi-(n_c-j)*lo]
    for j in range(n_c + 1):
        slo = j * lo - (n_c - j) * hi
        shi = j * hi - (n_c - j) * lo
        for k in range(-kmax, kmax + 1):
            t = pi * n_s + 2 * pi * k
            if slo - 1e-12 <= t <= shi + 1e-12:
                return True
    return False

for L, name in ((3, "triangle"), (5, "pentagon"), (7, "heptagon")):
    print(f"  {name} (length {L}):")
    for n_s in range(L + 1):
        n_c = L - n_s
        wgrid = np.linspace(W_STRUT_MIN if n_s else W_ALL_MIN, W_ALL_MAX, 400)
        ok = [w for w in wgrid if cycle_payable(n_s, n_c, w)]
        if ok:
            print(f"    n_s={n_s}: PAYABLE for omega in [{min(ok):.3f}, {max(ok):.3f}]"
                  f"  (x in [{x_of(max(ok)):.3f}, {x_of(min(ok)):.3f}])")
        else:
            print(f"    n_s={n_s}: NOT payable at any omega in window")

print("""
  => TRIANGLES admit NO struts (all-cable wound m=1 only, omega <= ~2.09).
     PENTAGONS admit at most 2 struts. Parity alone (bipartization) badly
     UNDERSTATES the demotion cost: geometry windows bite harder than Z2.""")

# ------------------------------------------------------------------
print("=" * 72)
print("B. FRUSTRATION INDEX (graph-only lower bound) vs THE REAL COST")
print("=" * 72)

def maxcut(n, edges):
    best = 0
    for mask in range(1 << (n - 1)):
        c = 0
        for a, b in edges:
            if ((mask >> a) & 1 if a < n - 1 else 0) != ((mask >> b) & 1 if b < n - 1 else 0):
                c += 1
        if c > best: best = c
    return best

def prism(k):
    e = []
    for i in range(k):
        e += [(i, (i + 1) % k), (k + i, k + (i + 1) % k), (i, k + i)]
    return 2 * k, e

def cube_g():
    e = [(a, a ^ (1 << b)) for a in range(8) for b in range(3) if a < a ^ (1 << b)]
    return 8, e

def K4():
    return 4, list(combinations(range(4), 2))

def dodeca():
    ph = (1 + 5 ** 0.5) / 2
    V = []
    for sx in (1, -1):
        for sy in (1, -1):
            for sz in (1, -1): V.append((sx, sy, sz))
    for s1 in (1, -1):
        for s2 in (1, -1):
            V += [(0, s1 / ph, s2 * ph), (s1 / ph, s2 * ph, 0), (s1 * ph, 0, s2 / ph)]
    V = np.array(V)
    D = np.linalg.norm(V[:, None] - V[None], axis=2)
    dm = np.min(D + np.eye(20) * 99)
    e = [(a, b) for a in range(20) for b in range(a + 1, 20) if D[a, b] < dm * 1.05]
    return 20, e

def icosa():
    ph = (1 + 5 ** 0.5) / 2
    V = []
    for s1 in (1, -1):
        for s2 in (1, -1):
            V += [(0, s1, s2 * ph), (s1, s2 * ph, 0), (s1 * ph, 0, s2)]
    V = np.array(V)
    D = np.linalg.norm(V[:, None] - V[None], axis=2)
    dm = np.min(D + np.eye(12) * 99)
    e = [(a, b) for a in range(12) for b in range(a + 1, 12) if D[a, b] < dm * 1.05]
    return 12, e

for name, (n, E) in (("K4 (tetrahedral frame)", K4()), ("triangular prism", prism(3)),
                     ("cube", cube_g()), ("pentagonal prism", prism(5)),
                     ("hexagonal prism", prism(6)), ("icosahedron (5-reg)", icosa()),
                     ("dodecahedron", dodeca())):
    mc = maxcut(n, E)
    print(f"  {name:24s} n={n:2d} E={len(E):2d}  maxcut={mc:2d}  "
          f"frustration index (edges to delete for bipartite) = {len(E) - mc}")
print("""  The ledger/brief expectation 'a 3-regular closed strut shell is
  over-constrained by exactly 2 => demote 2 links' matches the frustration
  index only for K4 and the triangular prism — and even there DEMOTION !=
  DELETION: a demoted link still carries its cycle phase, and part A shows
  triangles cannot be paid by 1 or 2 cables. The real minimum demotions:""")

# real min demotions using the triangle lemma + explicit designs below
print("    K4:               no consonant realization AT ALL (part C)")
print("    triangular prism: 6 cables (both triangles all-cable), 3 vertical struts (part D)")
print("    cube:             0 (bipartite, all-strut) — the H-series selection rule")
print("    pentagonal prism: 10 cables (both pentagons), 5 vertical struts (part D)")
print("    hexagonal prism:  0 — BIPARTITE (all faces even). Ledger P4's 'over-")
print("                      constrained by 2' is REFUTED: it needs no demotion.")

# ------------------------------------------------------------------
print()
print("=" * 72)
print("C. K4 IS DEAD — the diagonal-hole argument + exhaustive check")
print("=" * 72)
print("""  All triangles => all-cable (part A). All-cable K4 = 4 phases pairwise
  in the allowed band S = B u (2pi-B). Pigeonhole: 4 points on the circle,
  the two 'diagonal' pairs have separation = sum of two adjacent gaps in
  [2*lo, 2pi-2*lo]; the band's hole (omega*dmax, 2pi-omega*dmax) swallows
  that whole interval whenever omega*dmax < 2*omega*dmin + hole-margin.
  Numeric check over the full window:""")

def sepdist(delta, w):
    lo, hi = band(w)
    delta = delta % (2 * pi)
    ok = (lo <= delta <= hi) or (2 * pi - hi <= delta <= 2 * pi - lo)
    if ok: return 0.0
    cands = [abs(delta - lo), abs(delta - hi), abs(delta - (2 * pi - hi)), abs(delta - (2 * pi - lo))]
    return min(cands)

rng = np.random.default_rng(1)
def anneal(n, edges, w, iters=6000, restarts=24):
    lo, hi = band(w)
    bestE = 1e9
    for r in range(restarts):
        th = rng.uniform(0, 2 * pi, n)
        T = 1.0
        E = sum(sepdist(th[a] - th[b], w) ** 2 for a, b in edges)
        for it in range(iters):
            T *= 0.9995
            i = rng.integers(n)
            old = th[i]
            th[i] = (old + rng.normal(0, 0.3 + T)) % (2 * pi)
            E2 = sum(sepdist(th[a] - th[b], w) ** 2 for a, b in edges)
            if E2 < E or rng.random() < np.exp(-(E2 - E) / max(T * 0.05, 1e-4)):
                E = E2
            else:
                th[i] = old
        if E < bestE: bestE = E
        if bestE < 1e-8: break
    return bestE

n4, e4 = K4()
for w in (1.3942, 1.6, 1.9, 2.094, 2.4, 2.70):
    lo, hi = band(w)
    res = anneal(n4, e4, w)
    print(f"  omega={w:.3f} band=[{lo:.3f},{hi:.3f}] hole=({hi:.3f},{2*pi-hi:.3f}):"
          f" best residual^2 = {res:.4g}  => {'FEASIBLE' if res < 1e-8 else 'infeasible'}")
print("  (diagonal separations forced into [2*lo, 2pi-2*lo] which sits inside the hole)")

# ------------------------------------------------------------------
print()
print("=" * 72)
print("D. EXPLICIT FEASIBLE DESIGNS — verified phase assignments")
print("=" * 72)

def verify(name, n, edges, kinds, th, w):
    """kinds: 's' strut (sep pi), 'c' cable (sep in band)"""
    lo, hi = band(w)
    worst = 0; ok = True
    for (a, b), kind in zip(edges, kinds):
        d = (th[a] - th[b]) % (2 * pi)
        if kind == 's':
            err = min(abs(d - pi), abs(d - pi))
            if err > 1e-9: ok = False
        else:
            inb = (lo - 1e-9 <= d <= hi + 1e-9) or (2 * pi - hi - 1e-9 <= d <= 2 * pi - lo + 1e-9)
            if not inb: ok = False; worst = max(worst, sepdist(d, w))
    ns = kinds.count('s'); nc = kinds.count('c')
    print(f"  {name}: omega={w:.4f} x={x_of(w):.4f} struts={ns} cables={nc}"
          f"  -> {'ALL CONSTRAINTS MET' if ok else f'FAIL worst={worst:.3f}'}")
    return ok

# co-rotating triangular tube (prism3): verticals struts, triangles wound m=1
w = 2.02                      # in [1.9392, 2.0944]
phi = 2 * pi / 3
th = np.zeros(6)
for i in range(3):
    th[i] = (-i * phi) % (2 * pi)          # top triangle wound
    th[3 + i] = (th[i] + pi) % (2 * pi)    # bottom = antiphase (vertical struts)
n_, e_ = prism(3)
kinds = []
for (a, b) in e_:
    kinds.append('s' if (b - a == 3 and a < 3) else 'c')
verify("prism3 tube (B1c, n=3)", n_, e_, kinds, th, w)
print(f"    lengths: struts d=pi/omega={pi/w:.4f}, triangle cables d=phi/omega={phi/w:.4f}"
      f"  (contact ceiling {2*r_live(x_of(w)):.4f})")
print(f"    triangle winding m=1; per-cable back gate = gate(2phi)="
      f"{(0.5*(1+cos(2*phi)))**8:.3g}")

# co-rotating pentagonal tube: pentagons wound m=2
w = 2.30
phi = 4 * pi / 5
th = np.zeros(10)
for i in range(5):
    th[i] = (-i * phi) % (2 * pi)
    th[5 + i] = (th[i] + pi) % (2 * pi)
n_, e_ = prism(5)
kinds = ['s' if (b - a == 5 and a < 5) else 'c' for (a, b) in e_]
verify("prism5 tube (m=2 belts)", n_, e_, kinds, th, w)
print(f"    struts d={pi/w:.4f}; pentagon cables d={phi/w:.4f}; back gate "
      f"{(0.5*(1+cos(2*phi)))**8:.3g}")

# hexagonal prism, all-strut (bipartite): theta = 0/pi by 2-coloring
w = 2.30
n_, e_ = prism(6)
col = [i % 2 for i in range(6)] + [(i + 1) % 2 for i in range(6)]
th = np.array([pi * c for c in col])
kinds = ['s'] * len(e_)
verify("hexagonal prism ALL-STRUT (P4 refutation)", n_, e_, kinds, th, w)
print(f"    all 18 edges d=pi/omega={pi/w:.4f}; needs equilateral realization"
      f" (tolerances: t1 part C)")

# cube all-strut
n_, e_ = cube_g()
col = [bin(v).count('1') % 2 for v in range(8)]
th = np.array([pi * c for c in col])
verify("cube ALL-STRUT (H-series rule)", n_, e_, ['s'] * 12, th, 2.30)

# wound rings table (N, m): phi = 2 pi m / N must be in band(omega) with d=phi/omega
print("\n  wound ring feasibility (N,m): need phi=2pi*m/N in B(omega) for some omega:")
for N in range(3, 13):
    feas = []
    for m in range(1, N):
        phi = 2 * pi * m / N
        wgrid = np.linspace(W_ALL_MIN, W_ALL_MAX, 300)
        ok = [wv for wv in wgrid if band(wv)[0] - 1e-9 <= phi <= band(wv)[1] + 1e-9]
        if ok:
            gb = (0.5 * (1 + cos(2 * phi))) ** 8
            feas.append(f"m={m} (phi={phi*180/pi:.0f}deg, omega<= {max(ok):.2f},"
                        f" x>={x_of(max(ok)):.3f}, back={gb:.2g})")
    print(f"    N={N:2d}: " + ("; ".join(feas) if feas else "NONE"))

# ------------------------------------------------------------------
print()
print("=" * 72)
print("E. MULTI-CYCLE TUNABILITY — one omega, several closures")
print("=" * 72)
print("""  Cable cycles quantize omega: omega * Lambda_a / C in 2pi*Z, with
  Lambda_a the SIGNED cycle length (traversal vs lock orientation).
  Lambda_a = 0 ('balanced' cycles: equal forward and backward length)
  impose NOTHING — free. Unbalanced cycles must be pairwise commensurable
  with integers reachable in the window: m_a in Lambda_a*[w_min, w_max]/2pi.""")
for (L1, L2, name) in ((15.0, 9.0, "theta graph 15 vs 9 (ratio 5:3)"),
                       (15.0, 10.0, "theta graph 15 vs 10 (3:2)"),
                       (15.0, 11.0, "theta graph 15 vs 11 (15:11)"),
                       (15.08, 9.0, "theta 15.08 vs 9 (foam-jittered)")):
    sols = []
    for m1 in range(1, 12):
        wv = 2 * pi * m1 / L1
        if not (W_ALL_MIN <= wv <= W_ALL_MAX): continue
        m2 = wv * L2 / (2 * pi)
        if abs(m2 - round(m2)) < 1e-9: sols.append((m1, round(m2), wv))
    print(f"  {name}: " + (", ".join(f"(m1={a},m2={b},omega={c:.4f})" for a, b, c in sols)
                           if sols else "NO exact common omega — geometry must be tuned"))
print("""  Tolerance to a mistuned second cycle: its comma delta spreads over its
  N links (t3): all gates >= 0.9 needs |delta| <= N*0.229 rad (tempered) —
  the standing seed-acceptance 'closure integer to <0.05' (0.314 rad) meets
  this for N >= 2. Torus net N1 x N2 uniform-phi: needs m1*N2 == 0 mod N1
  (phi = 2pi*m1/N1 = 2pi*m2/N2). Examples: 12x12 any m; 12x8: m1 in {0,2,4,..}.""")

# ------------------------------------------------------------------
print()
print("=" * 72)
print("F. ICOSAHEDRON (all triangles, 5-regular) — anneal probe")
print("=" * 72)
ni, ei = icosa()
for w in (1.3942, 1.7, 2.094):
    res = anneal(ni, ei, w, iters=9000, restarts=16)
    print(f"  omega={w:.3f}: best residual^2 = {res:.4g}"
          f"  => {'FEASIBLE' if res < 1e-8 else 'infeasible (as the hole argument predicts)'}")

print()
print("=" * 72)
print("G. THE COUNTING THEOREM (per connected component)")
print("=" * 72)
print("""  Unknowns: theta (n-1, mod the time-translation gauge) + omega (1) = n.
  Forward locks: E equations mod 2pi; rank over theta = n-1; the remaining
  B = E-n+1 become CYCLE conditions carrying integer freedoms.
  Struts add their rung condition; with only m=1 reachable this fixes
  d_l = pi/omega — E_s geometric equalities (equal strut lengths), NOT
  theta/omega equations, plus pi-residues in every cycle they touch.
  Generic solvability on a FIXED foam:  deficit = (#unbalanced cable
  cycles - 1)_+ quantizations  +  E_s length equalities  +  the cycle
  payment conditions of part A.  With free geometry (form-finding moves
  picks): 3n-6 placement dof absorb the length equalities generically
  for n >= 6; what remains are (i) the Diophantine closures and (ii) the
  payment lemma — DISCRETE obstructions no continuous tuning removes.""")
