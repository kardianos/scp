#!/usr/bin/env python3
"""
t3_laplacian.py — v89 prestress theory, part 3: THE FORCE-DENSITY MATRIX.

Least-squares phase assignment over an infeasible consonant network
solves Q theta = B W phi with Q = B W B^T — a weighted graph Laplacian.

Construction (kernel-faithful):
  * every DIRECTED gate constraint is one row:  psi = theta_i - theta_j - phi
    (mod 2pi), phi = omega_i d/C for lock direction i->j (pass-2 gates).
  * a CABLE contributes its forward row; a STRUT both directions — the row
    pair sums to the pair-rung defect: struts ARE 2-cycles of the directed
    constraint system.
  * W = diag(w_row), w_row = K_row * p_gate/4: gate(psi) ~ 1 - (p/4) psi^2
    on the flat top (gate''(0) = -p/2), so maximizing total gated flow
    sum K*gate(psi) is, to 2nd order, minimizing sum (K p/4) psi^2:
    W IS DERIVED FROM gate''. K_row = k_dep*k_dep_m*geo*gpl*res*head*mob —
    the pass-2 multiplier chain with the gate factored out.
  * phases live on the torus: the LS is solved in a WINDING SECTOR —
    wrapped Gauss-Newton from the tree-lock (seed) initialization, which is
    also what kernel entrainment does dynamically (E6 P5, ~10 t.u.).

Checks: gauge/kernel structure, residual = weighted cycle projection, the
1/w comma law (temperament generalized), seam concentration, jittered-cube
gate budget vs the measured 0.6, the Z2 of interval networks.
"""
import numpy as np
from math import pi, cos, sqrt, acos
import sympy

p_gate = 8.0
rng = np.random.default_rng(7)

def gate(psi): return (0.5 * (1 + np.cos(psi))) ** p_gate
def wrap(a): return (np.asarray(a) + pi) % (2 * pi) - pi

def build(nc, edges, kinds, phis, weights):
    rows, tg, wr = [], [], []
    for (i, j), kind, ph, w in zip(edges, kinds, phis, weights):
        r = np.zeros(nc); r[i] = 1; r[j] = -1
        rows.append(r); tg.append(ph); wr.append(w)
        if kind == 's':
            rows.append(-r); tg.append(ph); wr.append(w)
    return np.array(rows), np.array(tg), np.array(wr)

def tree_lock(nc, edges, phis):
    """kernel seeder: BFS lock recursion theta_j = theta_i - phi (tree edges)"""
    th = np.zeros(nc); seen = {0}
    adj = {v: [] for v in range(nc)}
    for idx, (a, b) in enumerate(edges):
        adj[a].append((b, idx, +1)); adj[b].append((a, idx, -1))
    qq = [0]
    while qq:
        u = qq.pop(0)
        for v, idx, sgn in adj[u]:
            if v not in seen:
                seen.add(v)
                th[v] = th[u] - sgn * phis[idx]
                qq.append(v)
    return th

def solve_wrapped(M, phi, w, th0):
    """wrapped weighted LS (Gauss-Newton on the torus, sector of th0)"""
    th = th0.copy()
    Q = M.T @ (w[:, None] * M)
    Qp = np.linalg.pinv(Q)
    for it in range(300):
        psi = wrap(M @ th - phi)
        g = M.T @ (w * psi)
        if np.linalg.norm(g) < 1e-12: break
        th = th - Qp @ g
    return th, wrap(M @ th - phi), Q

print("=" * 72)
print("A. STRUCTURE CHECKS on a random consonant network")
print("=" * 72)
nc = 9
edges = [(0,1),(1,2),(2,3),(3,4),(4,0),(0,5),(5,6),(6,1),(2,7),(7,8),(8,3)]
kinds = ['c','c','s','c','c','s','c','c','c','s','c']
phis  = rng.uniform(1.5, 3.0, len(edges))
wts   = rng.uniform(0.2, 2.0, len(edges))
M, tg, wr = build(nc, edges, kinds, phis, wts)
th0 = tree_lock(nc, edges, phis)
th, psi, Q = solve_wrapped(M, tg, wr, th0)
ev = np.linalg.eigvalsh(Q)
print(f"  rows={M.shape[0]} (cables {kinds.count('c')}, struts {kinds.count('s')}x2), cells={nc}")
print(f"  ker Q: smallest eigenvalues {ev[:2].round(9)} -> exactly 1 zero mode ="
      f" global clock shift (time-translation gauge); rank(Q) = {np.linalg.matrix_rank(Q)} = n-1")
u, s, vt = np.linalg.svd(M.T)
null_dim = M.shape[0] - int(np.sum(s > 1e-10))
print(f"  cycle space dim (ker M^T) = {null_dim} = rows - (n-1)"
      f"   (graph cycles {len(edges)-nc+1} + strut 2-cycles {kinds.count('s')})")
print(f"  stationarity |M^T W psi| = {np.linalg.norm(M.T @ (wr*psi)):.2e}"
      f"  => residual psi* lies in W^-1 ker(M^T): DEFECTS LIVE IN THE CYCLE SPACE")

print()
print("=" * 72)
print("B. THE 1/w LAW — commas distribute in inverse proportion to stiffness")
print("=" * 72)
print("""  Single cycle with defect delta (mod 2pi): the wrapped-LS residual is
      psi_l = -delta * (1/w_l) / sum_k (1/w_k)     [series compliance]
  Uniform w => even tempering: E6's measured 'comma tempered evenly' is the
  uniform-stiffness special case. Wound ring12 (m=5), comma delta=0.6:""")
N = 12
edges = [(i, (i + 1) % N) for i in range(N)]
kinds = ['c'] * N
delta = 0.60
phis = np.full(N, 2 * pi * 5 / N + delta / N)
wts = np.ones(N)
M, tg, wr = build(N, edges, kinds, phis, wts)
# ring seeder init=7: SEQUENTIAL lock recursion around the loop; the seam
# (last link) carries the whole closure defect — the fleet-v2 measured dump
th0 = np.zeros(N)
for k in range(N - 1): th0[k + 1] = th0[k] - phis[k]
psi0 = wrap(M @ th0 - tg)
th, psi, Q = solve_wrapped(M, tg, wr, th0)
print(f"  seed (seq lock):   11 links psi=0, seam psi={psi0[-1]:+.4f}"
      f" (gate {gate(psi0[-1]):.4f})  <- fixed-ring_x fleet's measured seam dump")
print(f"  wrapped LS:        psi_l = {psi[0]:+.4f} each (predicted {-delta/N:+.4f}),"
      f" gates {gate(psi[0]):.4f}")
wts2 = np.ones(N); wts2[0] = 0.1
M2, tg2, wr2 = build(N, edges, kinds, phis, wts2)
th2, psi2, _ = solve_wrapped(M2, tg2, wr2, th0)
pred = -delta * (1 / wts2) / np.sum(1 / wts2)
print(f"  one soft link (w=0.1): psi_soft={psi2[0]:+.4f} (pred {pred[0]:+.4f}),"
      f" others {psi2[1]:+.4f} (pred {pred[1]:+.4f})")
print(f"  => commas CONCENTRATE on soft links (low flow*stiffness); the ring_m")
print(f"     seeder fix = set delta -> 0 so there is nothing to distribute.")
print(f"  comma cost min sum(w psi^2) = delta^2/sum(1/w):"
      f" uniform {delta**2/12:.5f}, soft-link {delta**2/np.sum(1/wts2):.5f}")

print()
print("=" * 72)
print("C. STRUT 2-CYCLES — the pair rung inside the same formalism")
print("=" * 72)
Delta = 0.3    # length error: phi = pi + Delta, Delta = omega*dd
M, tg, wr = build(2, [(0, 1)], ['s'], [pi + Delta], [1.0])
th0 = np.array([0.0, -(pi + Delta)])
th, psi, Q = solve_wrapped(M, tg, wr, th0)
print(f"  strut with phi = pi + {Delta}: round-trip defect (mod 2pi) = {wrap(-2*(pi+Delta)):+.4f}")
print(f"  wrapped LS splits it evenly: psi = ({psi[0]:+.4f}, {psi[1]:+.4f})"
      f" — E6 P2's tempered comma as the strut-2-cycle residual")
gt, gtr = gate(psi[0]), gate(2 * Delta)
print(f"  tempered gates {gt:.4f} x {gt:.4f} (product {gt*gt:.4f})  vs"
      f"  tree lock 1.0 x {gtr:.4f} (product {gtr:.4f})")
print(f"  => entrainment BEATS seeding on the two-gate product"
      f" (exp(-4x^2) vs exp(-8x^2), x = omega*dd)")

print()
print("=" * 72)
print("D. JITTERED CUBE — gate budget vs foam placement noise")
print("=" * 72)
print("""  Cube, all struts, target phi=pi; jitter: phi_l = pi + omega*dd_l,
  dd_l ~ N(0, sigma_d), omega=2.20 (abar 1.43 class). BFS = the kernel
  shell seeder (tree gates 1, co-tree edges eat whole cycle defects);
  LS = what entrainment anneals toward. 100 foams, ALL 24 directed gates:""")
cube_edges = [(a, a ^ (1 << b)) for a in range(8) for b in range(3) if a < a ^ (1 << b)]
om, abar = 2.20, 1.43
print("  sigma_d/a    BFS mean/min        LS mean/min")
for sig_frac in (0.02, 0.05, 0.08, 0.10, 0.12, 0.15):
    mb, nb, ml, nl = [], [], [], []
    for trial in range(100):
        dd = rng.normal(0, sig_frac * abar, 12)
        phis = pi + om * dd
        M, tg, wr = build(8, cube_edges, ['s'] * 12, phis, np.ones(12))
        th0 = tree_lock(8, cube_edges, phis)
        g0 = gate(wrap(M @ th0 - tg))
        th, psi, _ = solve_wrapped(M, tg, wr, th0)
        g1 = gate(psi)
        mb.append(g0.mean()); nb.append(g0.min())
        ml.append(g1.mean()); nl.append(g1.min())
    print(f"    {sig_frac:.2f}      {np.mean(mb):.3f} / {np.mean(nb):.3f}"
          f"        {np.mean(ml):.3f} / {np.mean(nl):.3f}")
print("""  Measured H1 cube: gates mean 0.60-0.64, min ~0 at seed => the foam's
  effective sigma_d/a ~ 0.10-0.15 for nearest-pick cubes. The LS/annealed
  assignment lifts the dead minimum gate by an order of magnitude and adds
  +0.05-0.1 mean — but CANNOT beat the sum rule: the cycle defect is
  invariant, only its distribution is designable. Gate budget per cycle:
  sum psi = -delta_c exactly. The only real cures are delta_c -> 0
  (equal-length mining P3 / plasticity F) or softer w on sacrificial links.""")

print()
print("=" * 72)
print("E. INTERVAL NETWORKS — integer incidence, omega-gauge, the Z2")
print("=" * 72)
print("""  p:q rows (pass 2): psi_ij = q theta_i - p theta_j - q omega_i d/C.
  Integer coefficient matrix M_pq (+q at i, -p at j). Gauge: theta ->
  theta + omega tau shifts psi_ij by (q w_i - p w_j) tau = 0 at lock:
  ker Q_pq = the OMEGA-VECTOR (uniform-omega nets: the constant vector).
  Solvability of M_pq theta == phi (mod 2pi): for every integer left-null
  vector y of M_pq:  y . phi == 0 (mod 2pi).
  Octave triangle (A-B 1:1, B-C 2 th_B - th_C, C-A th_C - 2 th_A):""")
Mpq = sympy.Matrix([[1, -1, 0], [0, 2, -1], [-2, 0, 1]])
y = sympy.Matrix([[2, 1, 1]])
print(f"  y = (2,1,1):  y.M_pq = {list(y*Mpq)}   (the left-null lattice generator)")
print("""  Both-gate (strut) rungs put HALF-PERIOD targets on the links:
  phi = (pi m1, pi m2, pi m3)  =>  y.phi = pi(2 m1 + m2 + m3) == 0 mod 2pi
  <=>  m2 + m3 == 0  (mod 2)   — a Z2 PARITY on rung integers.
  The Z2 is the ratio (clock period 2pi)/(pair-rung spacing pi) = 2; it is
  the same Z2 that makes 1:1 strut networks bipartite (all m_l = 1: every
  cycle needs an even strut count). Cable targets are full-period free
  (continuous phi), so cables carry no parity — they pay it (t2 part A).""")

print()
print("=" * 72)
print("F. W FROM THE KERNEL — the stiffness of a real ring link")
print("=" * 72)
kd, km, capv = 1.2, 2.0, 2.5
w2v, qv, s_pull, e_s0, r0v = 2.9, 1.2, 0.15, 1.0, 0.85
def r_live(x): return r0v * (e_s0 - s_pull * x * capv / (1 + s_pull)) ** (1 / 3)
def geo(d, ri, rj):
    if d >= ri + rj: return 0.0
    t = d * d - rj * rj + ri * ri
    a2 = (4 * d * d * ri * ri - t * t) / (4 * d * d)
    A = pi * a2 if a2 > 0 else 0.0
    return (A / (pi * r0v ** 2)) * (2 * r0v / d)
for name, x, d, N in (("comp12 ring link", 0.32054, 1.25, 12),
                      ("ring6 link", 0.128, 1.25, 6)):
    r = r_live(x); Em = x * capv
    gpl = cos(2 * pi / N) ** 4
    K = kd * km * geo(d, r, r) * gpl * 1.0 * (1 - Em / capv) * Em
    print(f"  {name}: geo={geo(d,r,r):.4f} gpl={gpl:.4f} head={1-Em/capv:.3f} Em={Em:.3f}"
          f"  K={K:.4f}/t.u.  w = K*p/4 = {K*p_gate/4:.4f} per rad^2")
print("""  (res=1 on-comb, mob_eff=Em for equal voices.) These w set which links
  absorb comma under jitter: ring6's 9x smaller gpl makes its links 9x
  SOFTER — same comma, 9x slower internal exchange to re-temper it.""")

print()
print("=" * 72)
print("G. THE DICTIONARY — force density method <-> consonant network")
print("=" * 72)
print("""  FDM tensegrity                      v89 consonant network
  ----------------------------------  ------------------------------------
  node coordinates x                  clock phases theta (+ pitch omega)
  member force density q = t/l        phase stiffness w_l = K_l p_gate/4
                                        (K_l = kd km geo gpl res head mob)
  equilibrium  D x = f                stationarity  Q theta = B W phi
  external load f                     retardation targets phi = omega d/C
  self-stress (ker of equilibrium     divergence-free circulating flux
    on edge space)                      sigma in ker B — winding (t5)
  mechanism (rigidity kernel)         ker Q beyond the omega-gauge: free
                                        retunings (balanced cycles' free
                                        omega; disconnected components)
  prestress stabilizes mechanisms     flux adds entrainment damping +
                                        gyroscopic precession, NOT static
                                        stiffness; protection = gate desert
  form-finding (choose q, solve x)    choose (topology, m, intervals),
                                        solve (picks, omega, theta, x*)""")
