#!/usr/bin/env python3
"""
v63/predictivity_audit.py  --  the decision-relevant audit:
rigorous INDEPENDENT-INPUTS  vs  INDEPENDENT-RELATIONS count for v59-v63.

Criterion (set with the user): if the structure reproduces MANY independent,
look-elsewhere-robust SM relations from FEW inputs (ratio >> 1) it is
over-determined = predictive (worth continuing); if relations ≈ inputs it is a
sophisticated fit (stop).

This codifies the project's OWN audits (v61/PROVEN_LEDGER.md,
v59/furey_construction/RIGOR_AUDIT.md) into one number, applying three honest
deflations a naive match-count ignores:
  (D1) SM-tree redundancy: m_W, m_Z follow from (v, g_W, θ_W) by definition --
       not independent reproductions.
  (D2) look-elsewhere: α (5/8 arbitrary templates hit to 0.03%) and quark Koide
       (12 simple rationals per RG band) carry ~zero evidential weight.
  (D3) load-bearing unproven ansätze (v_Higgs=784a², g_W²=5√α, the Z₂×Z₂
       selection rule) are ASSUMPTIONS, not free wins -- they belong on the input
       side or void the relation.

Each entry's evidential weight (0-1) is justified from the cited audits.
"""

SEP = "=" * 78

# (name, precision, independent_of_SM_tree, LE_robust, conditional_on_undriven,
#  weight, justification)  -- weight = effective independent look-elsewhere-robust
#  evidential content as a "relation" (0 = none, 1 = one clean independent relation)
RELATIONS = [
    ("lepton Koide Q=2/3 = dimG₂/dimSpin7", "1e-5", True, True, True, 0.80,
     "tight 5-digit match of a known relation to a FORCED ratio 14/21; but rests on "
     "the undriven D_lepton=28 selection (RIGOR_AUDIT): reproduces Koide GIVEN t²=1/2"),
    ("lepton phase φ=2/9 = Q/3", "7e-6", False, True, True, 0.25,
     "tight, but = Q/3 (the /3 structural via S₃, magnitude = Q) -> mostly redundant "
     "with the Koide entry; mechanism for Q absent"),
    ("sin²θ_W = 2/9 (Pati-Salam 5,2)", "~4%", True, False, True, 0.30,
     "2/9=0.2222 vs ~0.231 (scheme-dep, ~4%); independent of mass 2/9 (E5) but a soft, "
     "simple-rational match"),
    ("m_W", "0.04%", False, True, True, 0.00,
     "D1: SM-tree m_W = g_W·v/2 from (v, g_W) -- not an independent reproduction"),
    ("m_Z", "0.02%", False, True, True, 0.00,
     "D1: SM-tree m_Z = m_W/cosθ_W -- not independent"),
    ("v_Higgs = 28²·a_l²", "0.07%", True, False, True, 0.20,
     "D3: 'looks fitted, no mechanism' (RIGOR_AUDIT). If coincidence it is a 2nd free "
     "scale (an INPUT), not a relation"),
    ("quark Koide (Q_d,Q_u)", "~0.3%", True, False, True, 0.05,
     "D2: RG/scheme-dependent; '12 simple rationals per RG band' (PROVEN_LEDGER) -> "
     "look-elsewhere ~kills it"),
    ("α = 25/(324π²)", "0.03%", True, False, True, 0.05,
     "D2: '5 of 8 arbitrary templates hit α to 0.03%'; g_W²=5√α is α in disguise"),
    ("quark phases φ_d≈14α, φ_u≈−10α", "~1%", True, False, True, 0.15,
     "α-downstream (v63 LOOP_BOUND); structural N but ~1% and look-elsewhere concerns"),
]

# (name, kind, justification)
INPUTS = [
    ("a_lepton", "dimensionful scale",
     "unavoidable units choice; 1 in any theory (RIGOR_AUDIT)"),
    ("α", "dimensionless coupling",
     "genuine value-conjecture, NOT derived (RIGOR_AUDIT, PROVEN_LEDGER)"),
    ("Z₂×Z₂ selection rule D∈{28,35,63}", "undriven assignment",
     "N→(Bit-L,Bit-F) 'not derived'; sets t² per sector, load-bearing for Koide"),
    ("v_Higgs = 28²·a_l² ansatz", "load-bearing, unproven",
     "if coincidence, v_Higgs is a 2nd free scale (RIGOR_AUDIT #1 priority)"),
    ("g_W² = 5√α form", "load-bearing, unproven",
     "no Lagrangian; 'NOT a derivable law'; load-bearing for α(M_Z), m_W, m_Z"),
]


def banner(s):
    print(SEP); print(s); print(SEP)


def main():
    banner("PREDICTIVITY AUDIT  v59-v63  (independent inputs vs independent relations)")

    print("RELATIONS (claimed reproductions), with honest evidential weight:")
    print(f"  {'relation':40s} {'prec':>6s} {'indep':>6s} {'LE-rob':>7s} {'wt':>5s}")
    naive_count = 0
    honest_weight = 0.0
    for name, prec, indep, le, cond, wt, why in RELATIONS:
        naive_count += 1
        honest_weight += wt
        print(f"  {name:40s} {prec:>6s} {str(indep):>6s} {str(le):>7s} {wt:>5.2f}")
        print(f"      └ {why}")

    print()
    print("INPUTS / ASSUMPTIONS:")
    for name, kind, why in INPUTS:
        print(f"  • {name:38s} [{kind}]")
        print(f"      └ {why}")
    n_inputs = len(INPUTS)
    # the two unavoidable/genuine inputs vs the full assumption stack
    genuine_inputs = 2          # a_lepton + α (the project's headline)
    assumption_stack = n_inputs  # incl. the load-bearing unproven ansätze

    banner("THE RATIOS")
    print(f"  naive match-count           : {naive_count} reproductions")
    print(f"  naive 'inputs'              : {genuine_inputs} (a_lepton + α)")
    print(f"  NAIVE ratio                 : {naive_count/genuine_inputs:.1f}×  (looks over-determined)")
    print()
    print(f"  honest independent + LE-robust relation content (Σ weights) : {honest_weight:.2f}")
    print(f"  honest assumption count (incl. load-bearing ansätze)        : {assumption_stack}")
    print(f"  HONEST predictivity ratio   : {honest_weight/assumption_stack:.2f}×")

    banner("VERDICT")
    ratio = honest_weight / assumption_stack
    print(f"""  NAIVE: ~{naive_count} matches from ~{genuine_inputs} inputs reads as ~{naive_count/genuine_inputs:.0f}× over-determined.

  HONEST (after D1 SM-tree redundancy, D2 look-elsewhere, D3 load-bearing
  ansätze counted as assumptions):
    independent, look-elsewhere-robust relation content  ≈ {honest_weight:.1f}
    genuine assumptions (inputs + unproven load-bearing)  = {assumption_stack}
    => predictivity ratio ≈ {ratio:.2f}  ( < 1 )

  So the structure is NOT robustly over-determined.  By the agreed criterion
  (ratio >> 1 = predictive; ≈ 1 or below = fit), it sits on the FIT side: there
  are at least as many undriven assumptions as look-elsewhere-robust relations.
  This matches the project's own bottom line (PROVEN_LEDGER: 'a set of
  striking-but-unexplained numerical coincidences anchored on two genuine inputs').

  THE ONE GENUINE SURVIVOR (does not deflate away):
    lepton Koide Q = 2/3 = dimG₂/dimSpin7 at 1e-5.  A 5-digit match of a 40-year
    unexplained relation to a FORCED group-dimension ratio.  It is conditional on
    the undriven selection rule (D=28 ⇒ t²=1/2 ⇒ Q=2/3), so even this is
    'reproduces Koide given one undriven assignment' -- but it is the strongest,
    hardest-to-dismiss piece, and it is the real residual value of the program.
""")
    return {"naive_count": naive_count, "naive_ratio": naive_count / genuine_inputs,
            "honest_weight": honest_weight, "assumptions": assumption_stack,
            "honest_ratio": ratio}


def run_checks(R):
    banner("EMBEDDED SELF-CHECKS")
    ok = True

    def check(name, cond):
        nonlocal ok
        ok = ok and bool(cond)
        print(f"  [{'PASS' if cond else 'FAIL'}] {name}")

    check("naive ratio looks over-determined (> 3×)", R["naive_ratio"] > 3)
    check("honest predictivity ratio is NOT over-determined (< 1)", R["honest_ratio"] < 1)
    check("m_W, m_Z carry zero independent weight (SM-tree, D1)",
          all(wt == 0.0 for n, *_, wt, _ in
              [(r[0],) + r[1:] for r in RELATIONS] if "m_W" in n or "m_Z" in n))
    check("the load-bearing ansätze are on the input side (D3)",
          any("v_Higgs" in i[0] for i in INPUTS) and any("g_W" in i[0] for i in INPUTS))
    check("one survivor carries the highest weight (lepton Koide)",
          max(r[5] for r in RELATIONS) == RELATIONS[0][5] and RELATIONS[0][5] >= 0.8)

    print(SEP)
    print("ALL CHECKS PASSED" if ok else "SOME CHECKS FAILED")
    print(SEP)
    return ok


if __name__ == "__main__":
    R = main()
    ok = run_checks(R)
    import sys
    sys.exit(0 if ok else 1)
