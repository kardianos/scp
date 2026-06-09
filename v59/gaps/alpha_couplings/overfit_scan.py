#!/usr/bin/env python3
"""
overfit_scan.py  —  Gap G3/G2: HONEST overfitting quantification.

Central question for the rigor ethos: how SPECIAL are the v59 structural forms?
If thousands of "structural-looking" formulas built from small integers and pi
hit alpha to the same precision, then a 0.03% match is numerology, not evidence.

We do three scans and report look-elsewhere / multiplicity numbers:

  A. EW form family:  alpha = p^2 / (q * pi^2)  with small-integer p,q
                      (the v59 form is p=5, q=324). How many beat 0.03%?
  B. IR form family:  -ln(alpha) + c*alpha = a*pi^2/b  for small int a,b and
                      c in a small candidate set. How many hit alpha0?
  C. Density of "simple constants" near alpha:  count how many short
     closed-forms (built from {pi,e,2,3,5,7, integers, sqrt, ln}) land within
     a target window of alpha0^-1 and alpha(M_Z)^-1.  This estimates the
     baseline coincidence rate: P(a random simple form hits within eps).

A match is "evidence" only to the extent (1/N_competitors) is small AND the
form is independently motivated. We report N at each precision band.
"""
import math
import itertools

PI = math.pi
ALPHA0_INV   = 137.035999084
ALPHA0       = 1.0/ALPHA0_INV
ALPHA_MZ_INV = 127.951
ALPHA_MZ     = 1.0/ALPHA_MZ_INV

def reldev(a, b): return abs(a-b)/abs(b)

print("="*78)
print("OVERFITTING SCAN — how unique are the v59 alpha forms?")
print("="*78)

# ---------------------------------------------------------------------------
# A. EW family: alpha = p^2 / (q pi^2).  v59 = 5^2/(324 pi^2).
#    324 = 4*81 = (18)^2 = (2*9)^2.  So q itself is structured (a perfect square).
#    Restrict to the v59-style family alpha = (p/(m pi))^2 = p^2/(m^2 pi^2),
#    i.e. sqrt(alpha) is a RATIONAL multiple of 1/pi.  Count p,m hitting target.
# ---------------------------------------------------------------------------
print("\n[A] EW family: sqrt(alpha) = p/(m*pi),  alpha = p^2/(m^2 pi^2)")
print("    v59 form: p=5, m=18  ->  alpha = 25/(324 pi^2)")
target = ALPHA_MZ
bands = [1e-2, 1e-3, 3e-4, 1e-4]   # 1%, 0.1%, 0.03%, 0.01%
hits = {b: [] for b in bands}
PMAX, MMAX = 30, 60
for p in range(1, PMAX+1):
    for m in range(1, MMAX+1):
        a = (p/(m*PI))**2
        rd = reldev(a, target)
        for b in bands:
            if rd <= b:
                hits[b].append((p, m, a, rd))
for b in bands:
    # keep only reduced (gcd) forms to avoid trivial scalings
    reduced = [(p,m,a,rd) for (p,m,a,rd) in hits[b] if math.gcd(p,m)==1]
    print(f"    within {b*100:6.3f}%: {len(reduced):3d} reduced (p,m) forms "
          f"(p<= {PMAX}, m<= {MMAX})")
    if b <= 3e-4:
        for (p,m,a,rd) in sorted(reduced, key=lambda x:x[3])[:6]:
            flag = "  <-- v59" if (p==5 and m==18) else ""
            print(f"        p={p:2d} m={m:2d}  alpha^-1={1/a:8.3f}  reldev={rd*100:.4f}%{flag}")

# ---------------------------------------------------------------------------
# B. IR family: -ln(alpha) + c*alpha = a*pi^2/b.  v59: c=2, a=1, b=2.
#    The mechanism fixes RHS = pi^2/2 (a=1,b=2). The freedom is the +c*alpha
#    correction (c reverse-engineered). Quantify: how many (a,b,c) hit alpha0
#    to 0.01%, and is c=2,a=1,b=2 special?
# ---------------------------------------------------------------------------
print("\n[B] IR family: -ln(alpha) + c*alpha = a*pi^2/b")
print("    v59 form: a=1, b=2, c=2.  (RHS pi^2/2 is mechanism-fixed; c is the fit.)")
lhs0 = -math.log(ALPHA0)   # = 4.92024...
print(f"    measured -ln(alpha0) = {lhs0:.6f};  needed RHS-LHS correction = {PI**2/2 - lhs0:.6f}")
# For fixed RHS = a pi^2/b, the required c is c = (RHS - (-ln a0))/alpha0.
candidate_c = list(range(-5,9))   # small integers
ir_hits = []
for a in range(1,5):
    for b in range(1,13):
        rhs = a*PI**2/b
        # solve -ln(x)+c*x = rhs for x near alpha0, check reldev to alpha0
        for c in candidate_c:
            # Newton solve
            x = ALPHA0
            ok = True
            for _ in range(100):
                fx = -math.log(x)+c*x-rhs
                fpx = -1/x + c
                if abs(fpx) < 1e-12: ok=False; break
                xn = x - fx/fpx
                if xn <= 0: ok=False; break
                if abs(xn-x) < 1e-15: x=xn; break
                x = xn
            if not ok: continue
            rd = reldev(1/x, ALPHA0_INV)
            if rd <= 1e-3 and math.gcd(a,b)==1:
                ir_hits.append((a,b,c,1/x,rd))
ir_hits.sort(key=lambda t:t[4])
print(f"    (a,b,c) forms (a<=4,b<=12,c in [-5,8]) hitting alpha0 within 0.1%: {len(ir_hits)}")
for (a,b,c,ainv,rd) in ir_hits[:12]:
    flag = "  <-- v59" if (a==1 and b==2 and c==2) else ""
    print(f"        a={a} b={b} c={c:+d}: RHS={a*PI**2/b:.4f}  alpha^-1={ainv:.4f}  reldev={rd*100:.4f}%{flag}")

# ---------------------------------------------------------------------------
# C. Baseline coincidence rate: density of short closed forms near alpha^-1.
#    Generate a large bag of "simple" constants and count how many land within
#    a relative window of each target. This estimates P(random simple form hits).
# ---------------------------------------------------------------------------
print("\n[C] Baseline density: simple closed forms near alpha^-1 targets")
consts = {'1':1.0,'2':2.0,'3':3.0,'5':5.0,'7':7.0,'pi':PI,'e':math.e,
          'pi^2':PI**2,'sqrt2':math.sqrt(2),'sqrt7':math.sqrt(7),
          'ln2':math.log(2),'phi':(1+math.sqrt(5))/2}
vals = list(consts.values())
names = list(consts.keys())
# build pool: a*b, a/b, a*b+c, a*b^2 for small structure
pool = []
ints = list(range(1,40))
# integer * const, integer * const^2, ratios of products
for i in ints:
    for v,nm in zip(vals,names):
        for op,lbl in [(i*v, f"{i}*{nm}"), (i*v*v, f"{i}*{nm}^2"),
                       (i/v if v!=0 else None, f"{i}/{nm}"), (v**i if v<3 and i<8 else None, f"{nm}^{i}")]:
            if op is not None and 50 < op < 300:
                pool.append((op,lbl))
# also products of two consts times integer
for i in range(1,20):
    for (v1,n1),(v2,n2) in itertools.combinations_with_replacement(zip(vals,names),2):
        op = i*v1*v2
        if 50 < op < 300:
            pool.append((op,f"{i}*{n1}*{n2}"))
print(f"    pool size of 'simple' constants in (50,300): {len(pool)}")
# expected count under uniform null: pool covers (50,300), width 250.
# a window of half-width eps*tgt around tgt has width 2*eps*tgt; expected
# hits = (pool_size) * (2*eps*tgt)/250.
for tgt_inv, tname in [(ALPHA0_INV,'alpha0^-1'), (ALPHA_MZ_INV,'alpha(M_Z)^-1')]:
    for eps in (1e-2, 1e-3, 3e-4, 1e-4):
        n = sum(1 for (v,_) in pool if reldev(v,tgt_inv) <= eps)
        exp_null = len(pool) * (2*eps*tgt_inv)/250.0
        print(f"    {tname:<14s} within {eps*100:6.3f}%: {n:4d} hits  "
              f"(expected under uniform-null ~ {exp_null:.2f})")
    # show the closest few to alpha0^-1
    if tname=='alpha0^-1':
        closest = sorted(pool, key=lambda x:reldev(x[0],tgt_inv))[:8]
        print(f"      closest simple forms to {tname}={tgt_inv:.3f}:")
        for v,lbl in closest:
            print(f"        {lbl:<14s} = {v:.4f}  reldev={reldev(v,tgt_inv)*100:.4f}%")

# ---------------------------------------------------------------------------
# D. The REAL freedom: number of *families* (look-elsewhere across templates).
#    Within a fixed 2-param family the v59 forms are nearly unique at 0.03%
#    (scans A,B). The honest overfitting risk is the number of DISTINCT simple
#    templates one could have chosen. Enumerate a representative template set
#    and count how many admit *some* small-integer instance hitting alpha0^-1
#    to <=0.03%. If many templates can, then "alpha has a structural form" is
#    cheap and the specific template is the (post-hoc) choice.
# ---------------------------------------------------------------------------
print("\n[D] Look-elsewhere across TEMPLATES (the real overfitting axis)")
def best_int_fit(template, target_inv, irange=range(1,400), jrange=range(1,60)):
    """template(i,j)->value; return best reldev to target_inv over small i,j."""
    best = (1e9, None)
    for i in irange:
        for j in jrange:
            try:
                v = template(i,j)
            except (ValueError, ZeroDivisionError, OverflowError):
                continue
            if v <= 0: continue
            rd = reldev(v, target_inv)
            if rd < best[0]:
                best = (rd, (i,j))
    return best
templates = {
    "i/(j alpha-style) p^2/(q pi^2)":   lambda i,j: (j*PI)**2 / i**2,   # inverse of family A
    "i exp(pi^2/j)":                    lambda i,j: i*math.exp(PI**2/j),
    "i pi^2 + j":                       lambda i,j: i*PI**2 + j,
    "i pi^2 - j":                       lambda i,j: i*PI**2 - j,
    "(i+ j/pi) pi^2 /-style i*j":       lambda i,j: i*j,
    "i/(j) * e^pi":                     lambda i,j: (i/j)*math.exp(PI),
    "i pi - j":                         lambda i,j: i*PI - j,
    "i ln(j)":                          lambda i,j: i*math.log(j+1),
}
print("    template best small-int fit to alpha0^-1 = 137.036:")
n_good = 0
for nm, tpl in templates.items():
    rd, ij = best_int_fit(tpl, ALPHA0_INV)
    good = rd <= 3e-4
    n_good += good
    print(f"      {nm:<34s} best reldev={rd*100:7.4f}%  at {ij}{'  <=0.03%' if good else ''}")
print(f"    => {n_good}/{len(templates)} arbitrary templates already hit 0.03% with small ints.")

print("\n" + "="*78)
print("INTERPRETATION")
print("="*78)
print("""  - WITHIN a fixed 2-parameter family the v59 forms are surprisingly unique:
      [A] alpha=(p/(m pi))^2: ZERO other reduced (p<=30,m<=60) forms beat 0.03%;
          the v59 (5,18) is the sole sub-0.1% hit.
      [B] -ln a + c alpha = a' pi^2/b': the v59 (a'=1,b'=2,c=2) is the UNIQUE
          small-integer hit to 0.1%.
    This means the *parametrization* is tight, not loose -- a genuine (if modest)
    point in v59's favor that earlier prose understated.
  - The real overfitting risk is the LOOK-ELSEWHERE over which TEMPLATE/family
    you pick [D]: several unrelated simple templates already reach 0.03% with
    small integers. So 'alpha has a clean structural form' is cheap; the content
    is whether the chosen template's integers are FORCED by the same algebra that
    fixes everything else (5=h^v, 2/9, dim Cl(3,1)=16).
  - Bottom line: precision alone is NOT the evidence. The evidence is the SHARED
    PROVENANCE of the integers. v59's EW form reuses 5 and 2/9 (independently
    motivated) -> moderate weight. The IR form's pi^2/2 reuses 16=dim Cl(3,1)
    but the +2 alpha is a free fitted term -> weak. g_W^2=5 sqrt(a) adds no new
    information beyond the EW value (shown in verify_forms.py).""")
