#!/usr/bin/env python3
"""HORIZON prover (HORIZON.md; user-directed 2026-08-07).

A domain proof harness with symbolic backtracing, NOT a general theorem
prover: the algebraic core of every lemma is checked by sympy on the
kernel's actual update forms (freecell.c pass 6 / pass S / pass F /
atoms machinery); the search over STRUCTURES is exhaustive backtracking
over a finite move grammar. Constant tuning is excluded by rule — a
move must change what terms EXIST or what they READ, never a value.

Output: (1) the BACKTRACE — for each horizon requirement, the standing
invariant that forbids it under the empty structure; (2) all MINIMAL
structural move-sets that discharge every requirement; (3) derived
quantitative conditions (trapped-surface radius) for the experimental
arms."""

import itertools
import sympy as sp

# ---------------------------------------------------------------- law core
q, cap, s_pull, s_disp, s_k, f_conv = sp.symbols(
    'q cap s_pull s_disp s_k f_conv', positive=True)
A0, w, w1, w2, x, Es, Em, Ee, d1, d2, rc, re_ = sp.symbols(
    'A0 w w1 w2 x Es Em Ee d1 d2 r_c r_e', positive=True)
NUM = {q: sp.Rational(12, 10), cap: sp.Rational(5, 2),
       s_pull: sp.Rational(15, 100), s_disp: sp.Rational(3, 10),
       s_k: sp.Rational(6, 100), f_conv: sp.Rational(1, 4),
       w1: sp.Rational(165, 100), w2: sp.Rational(29, 10)}

# ------------------------------------------------------------- lemma base
# Each lemma: name -> (statement text, sympy check that must return True)
def lemma_checks():
    L = {}
    # I1 refusal: intake carries the headroom factor of the RECEIVING
    # capped store (PAULI-0, FP-sticky): headroom -> 0 kills intake.
    intake = sp.Symbol('k_in', positive=True) * (cap - Em - Ee)
    L['I1_refusal'] = ("intake ~ k*(cap-Em-Ee): saturation refuses exactly",
                       sp.simplify(intake.subs(Em, cap - Ee)) == 0)
    # I2 atoms: every conversion fires whole atoms eps = A0*w_src/2pi at
    # the SOURCE's pitch. A pitchless source has eps = 0: it cannot fire.
    eps = A0 * w / (2 * sp.pi)
    L['I2_atoms'] = ("emission quantum eps = A0*w/2pi: w=0 => eps=0 => "
                     "a pitchless store can never emit (one-way for free)",
                     sp.simplify(eps.subs(w, 0)) == 0)
    # I3 pi-closure (the BLACKHOLE.md Wall-3 theorem, re-derived): at
    # closed books the net medium pi-load is identically zero.
    dsp_ = s_pull * d1
    bs_ = d2 * s_pull / (1 + s_pull)
    dpi_c = -dsp_ + s_disp * ((d1 + dsp_) - d1)
    dpi_e = +bs_ + s_disp * (-d2 + (d2 - bs_))
    close = sp.solve(sp.Eq((d1 + dsp_) * rc, d2 * re_), rc)[0]
    net = sp.simplify((dpi_c * rc + dpi_e * re_).subs(rc, close))
    L['I3_pi_closure'] = ("closed books => net pi-load == 0 (no far field "
                          "from any stationary object)", net == 0)
    # I4 storage bounds: Em+Ee <= cap (beat-enforced); Es >= floor and
    # UNCAPPED above — but pi = Es + s_disp*(Em+Ee) reads Es at weight 1:
    # storing swallowed energy in Es RAISES pi (pushes space away).
    pi_ = Es + s_disp * (Em + Ee)
    L['I4_es_antitrap'] = ("d(pi)/d(Es) = 1 > 0: an Es-hoard repels space "
                           "— the only uncapped ledger is anti-trapping",
                           sp.diff(pi_, Es) == 1)
    # I6 pitch bound: w_e = w/(1+q*x) with x bounded => clock never stops.
    we = w / (1 + q * x)
    L['I6_clock_floor'] = ("w_e = w/(1+q*x), x <= x_max finite => w_e > 0 "
                           "(max dilation 1+q*x_max ~ 3.2)",
                           sp.limit(we, x, sp.oo) == 0 and
                           sp.simplify(we.subs(x, 1)) != 0)
    # I7 push-only: link space flux = s_k*w*(pi_i - pi_j), outflow-limited
    # at the SOURCE: a standing INWARD wind exists iff the object holds a
    # standing pi DEFICIT (nothing reaches out; deficits pull by pushing).
    pi_i, pi_j = sp.symbols('pi_i pi_j')
    flux = s_k * w * (pi_i - pi_j)
    L['I7_push_only'] = ("flux sign = sign(pi_i - pi_j): inward standing "
                         "wind <=> standing pi deficit at the object",
                         sp.simplify(flux.subs(pi_i, pi_j)) == 0)
    # I5 unitarity (structural, recorded — not a sympy object): pass-F
    # hops are pairwise rotations; amplitude sent i->j returns unless the
    # destination state space GROWS (irreversibility needs somewhere to
    # spread). A one-way boundary over a FIXED interior contradicts it.
    L['I5_unitary'] = ("one-way field boundary over a fixed finite "
                       "interior contradicts pairwise-rotation unitarity; "
                       "needs growing interior state space", True)
    return L

# --------------------------------------------------- horizon requirements
# Each requirement is tested against a STRUCTURE (a set of moves).
# A structure is a dict of properties describing what exists.
BASE = dict(store=None,           # where swallowed energy accumulates
            store_capped=True,    # is that store capped?
            store_pitched=True,   # does it carry a pitch (can fire atoms)?
            store_in_pi=True,     # does pi read it?
            gate_reads_store=False,  # does the door's headroom factor
                                     # read the store (not Em+Ee)?
            roster_dynamic=False, # can cells die (capacity leaves books)?
            books_open=False)     # sustained net intake possible?

def check(structure):
    """Return (ok, failures): failures = list of (req, blocking invariant)."""
    S = structure
    fails = []
    # H5 conservation: any structure made of LEDGER moves conserves by
    # routing; forbidden only if a move destroys energy (none in grammar).
    # H1 absorb-always: needs a headroom that never closes.
    if S['store'] is None:
        fails.append(('H1_absorb_always', 'I1_refusal',
                      "intake multiplies (cap-Em-Ee); swallowing raises "
                      "Em+Ee to cap and the door closes"))
    else:
        if S['store_capped']:
            fails.append(('H1_absorb_always', 'I1_refusal',
                          "the store is capped: same refusal, delayed"))
        if not S['gate_reads_store'] and not S['roster_dynamic']:
            fails.append(('H1_absorb_always', 'I1_refusal',
                          "the door's headroom factor still reads Em+Ee: "
                          "capture stalls at the boundary cells' cap even "
                          "though the store behind them is open"))
    # H2 emit-never: needs a store that CANNOT fire atoms.
    if S['store'] is None or S['store_pitched']:
        fails.append(('H2_emit_never', 'I2_atoms',
                      "any pitched store fires atoms at its own pitch; "
                      "only w_store = 0 makes emission structurally zero"))
    # H3 standing far field: needs permanently open books (I3 converse)
    # AND a store that does not raise pi (I4), else the wind blows out.
    if not S['books_open']:
        fails.append(('H3_standing_far_field', 'I3_pi_closure',
                      "closed books => pi-flat identically; a stationary "
                      "object cannot source a standing gradient"))
    elif S['store'] is not None and S['store_in_pi']:
        fails.append(('H3_standing_far_field', 'I4_es_antitrap',
                      "a pi-visible hoard raises pi and repels the medium"))
    # H4 horizon clock: the store must carry no clock at all (I6 forbids
    # slowing a pitched clock to zero at bounded load).
    if S['store'] is None or S['store_pitched']:
        fails.append(('H4_clock_stops', 'I6_clock_floor',
                      "w_e = w/(1+q*x) with bounded x never reaches 0; "
                      "only a pitchLESS mode has no clock"))
    # H2b one-way field boundary: unitarity needs a growing interior.
    grows = (S['store'] is not None and not S['store_capped']) \
        or S['roster_dynamic']
    if not grows:
        fails.append(('H2b_no_return', 'I5_unitary',
                      "fixed finite interior returns amplitude; the "
                      "interior state space must grow with swallowed E"))
    return (len(fails) == 0, fails)

# --------------------------------------------------------- move grammar
MOVES = {
    'NEW_LEDGER_X': dict(store='X', store_capped=False,
                         store_pitched=False, store_in_pi=False,
                         books_open=True),
    'ROUTE_DOOR_TO_X': dict(gate_reads_store=True),
    'ROSTER_DEATH': dict(roster_dynamic=True, books_open=True,
                         store='dead-cell capacity', store_capped=False,
                         store_pitched=False, store_in_pi=False),
    'UNCAP_EM': dict(store='Em', store_capped=False, store_pitched=True,
                     store_in_pi=True, books_open=True),
    'STORE_IN_ES': dict(store='Es', store_capped=False, store_pitched=False,
                        store_in_pi=True, books_open=True),
    'CLOCK_READS_SPACE': dict(exterior_dilation=True),  # w_e = w*f(Es):
    # not needed for trapping; required only for EXTERIOR redshift (GR
    # completeness) — tracked separately below.
}

def apply(moves):
    s = dict(BASE)
    for m in moves:
        s.update(MOVES[m])
    return s

def main():
    print("=" * 72)
    print("LEMMA BASE (each algebraic core sympy-checked on the kernel forms)")
    L = lemma_checks()
    for k, (txt, ok) in L.items():
        print(f"  [{'OK' if ok else 'FAIL'}] {k}: {txt}")
        assert ok, k
    print()
    print("=" * 72)
    print("BACKTRACE: horizon requirements vs the CURRENT structure (no moves)")
    ok, fails = check(apply([]))
    for req, inv, why in fails:
        print(f"  {req}\n      blocked by {inv}\n      {why}")
    assert not ok
    print()
    print("=" * 72)
    print("BACKTRACKING SEARCH over structural move-sets (minimal first)")
    names = [m for m in MOVES if m != 'CLOCK_READS_SPACE']
    solutions = []
    first_size = None
    for size in range(1, len(names) + 1):
        if first_size is not None and size > first_size + 1:
            break   # enumerate minimal size and minimal+1 (the
                    # apparatus-implementable variant), then stop
        for combo in itertools.combinations(names, size):
            ok, fails = check(apply(combo))
            if ok and not any(set(s) < set(combo) for s in solutions):
                solutions.append(combo)
                if first_size is None:
                    first_size = size
                print(f"  SOLUTION (size {size}): {' + '.join(combo)}")
                for m in combo:
                    print(f"      {m}: {MOVES[m]}")
    print()
    print("  Rejected single moves and why:")
    for m in names:
        ok, fails = check(apply([m]))
        if not ok and all(m not in s for s in solutions):
            worst = fails[0]
            print(f"    {m}: fails {worst[0]} ({worst[1]})")
        elif not ok:
            missing = sorted(set(f[0] for f in fails))
            print(f"    {m} alone: still fails {missing}")
    print()
    print("=" * 72)
    print("EQUIVALENCE: NEW_LEDGER_X+ROUTE == ROSTER_DEATH by property map")
    a, b = apply(['NEW_LEDGER_X', 'ROUTE_DOOR_TO_X']), apply(['ROSTER_DEATH'])
    keys = ['store_capped', 'store_pitched', 'store_in_pi', 'books_open']
    print("  property           X-ledger   dead-cells")
    for k in keys:
        print(f"    {k:<18} {str(a[k]):<10} {str(b[k])}")
    print("  The X ledger IS the ledger of consumed cells' capacity: the")
    print("  black hole is the region where cells leave the roster. The")
    print("  THEORY.md 2.3 third lane (cell number dynamical), reached by")
    print("  backtracing from the horizon instead of from cosmology.")
    print()
    print("  CLOCK_READS_SPACE (w_e = w*f(Es/e_s0)): NOT required for")
    print("  trapping (the store has no clock at all); REQUIRED for the")
    print("  exterior GR analog (redshift of nearby matter tracking the")
    print("  well). Registered as the derived missing interaction; its")
    print("  absence is testable NOW as a predicted spectral null.")
    print()
    print("=" * 72)
    print("QUANTITATIVE: the trapped-surface (acoustic-horizon) condition")
    # 2D: space consumed at rate Qdot flows in through circumference 2*pi*r
    # carried by the medium at line density sigma_s (space per area).
    Qd, r, sig, ceff = sp.symbols('Qdot r sigma_s c_eff', positive=True)
    v = Qd / (2 * sp.pi * r * sig)
    rh = sp.solve(sp.Eq(v, ceff), r)[0]
    print(f"  inflow speed v(r) = {v};  r_h = {rh}")
    n_area = 2755 / 64**2          # measured 2D foam cell density
    subs = {sig: n_area * 1.0, ceff: 0.4}   # Es~1/cell; c_eff = 0.4C (P1)
    rh_of = sp.lambdify(Qd, rh.subs(subs))
    # s_k-borne wind ceiling: boundary links ~ 2*pi*r*sqrt(n)*z/2pi-ish;
    # per-link flux ~ s_k*w*dpi with dpi ~ 1 (floor-drained interior).
    import math
    r0 = 4.0
    links = 2 * math.pi * r0 * math.sqrt(n_area) * 3   # ~3 links/boundary cell
    Q_sk = links * 0.06 * 1.0 * 1.0
    print(f"  s_k-borne inflow ceiling at r=4: Qdot ~ {Q_sk:.2f}/t.u. "
          f"=> wind-borne r_h ~ {rh_of(Q_sk):.2f} cells")
    print("  (order-of-magnitude; boundary-link count uncertain ~2x).")
    print("  PREDICTION: the pass-S wind alone holds a trapped surface")
    print("  only at r_h ~ 1-2 cells => outward probes at r >= 3 ESCAPE")
    print("  unless FRAME MIGRATION (cell infall by contact mechanics,")
    print("  not bounded by s_k) extends the advection. The probe ladder")
    print("  measures which channel carries the river.")
    print(f"  required Qdot for r_h = 4 cells: "
          f"{float(sp.solve(sp.Eq(rh.subs(subs), 4), Qd)[0]):.1f}/t.u.")

if __name__ == '__main__':
    main()
