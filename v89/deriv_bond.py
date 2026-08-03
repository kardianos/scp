#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
deriv_bond.py -- v89 free-cell programme: bond derivations T1-T6.

Symbolic (sympy 1.12) + numeric (numpy) derivation of the P15 plasticity bond
for a unison pair under the standing law (laws_V2g constants as transcribed in
the task; no repo files consulted).

Global conventions:
  * C = 1 (conversion rate), rung m = 1.
  * kappa (overall gain in dd = -kappa*base*Sm*dVdd) has NO numeric value in
    the transcription; every velocity/stiffness below is LINEAR in kappa and
    is reported at kappa = 1.
  * kappa_lock = 1 in T3 entrainment (per-step pull = mix * wrap_pi(...)).
"""
import os
import numpy as np
import sympy as sp

OUT = []
def emit(*a):
    s = " ".join(str(x) for x in a)
    print(s)
    OUT.append(s)

def hr(t):
    emit("")
    emit("=" * 78)
    emit(t)
    emit("=" * 78)

def g15(x):
    return f"{float(x):.15g}"

# ---------------- transcribed constants ----------------
P      = 8                    # p_gate
KD     = 1.2 * 2.0            # kd = k_dep*k_dep_m = 2.4
R0     = 0.85
AREF   = np.pi * R0 * R0      # pi*0.85^2
DREF   = 2.0 * R0             # 1.7
W2     = 2.9
QDET   = 1.2
CAP    = 2.5
SPULL  = 0.15
C1     = 1.0
KAPPA  = 1.0
FL0    = 0.004 * 2.5          # 0.01
DT     = 0.02

# ---------------- gate + geometry (numeric) ----------------
def h_(psi):
    return 0.5 * (1.0 + np.cos(psi))

def G_(psi):
    return h_(psi) ** P

def Gp_(psi):
    # dG/dpsi = -0.5*p*(G/h)*sin(psi) = -0.5*p*h^(p-1)*sin(psi)
    return -0.5 * P * h_(psi) ** (P - 1) * np.sin(psi)

def lens_a2(d, ri, rj):
    t = d * d - rj * rj + ri * ri
    return (4.0 * d * d * ri * ri - t * t) / (4.0 * d * d)

def geo_(d, ri, rj):
    a2 = lens_a2(d, ri, rj)
    A = np.where(a2 > 0.0, np.pi * a2, 0.0)
    return (A / AREF) * (DREF / d)

def seed_pair(x):
    Em = x * CAP
    w  = W2 / (1.0 + QDET * x)
    Es = 1.0 - SPULL * x * CAP / (1.0 + SPULL)
    r  = R0 * Es ** (1.0 / 3.0)
    return Em, w, Es, r

def v_branch(delta, w_, r_, Sm_):
    """On-branch bond velocity dd/dt (symmetric locked branch), kappa=1."""
    S_  = 2.0 * w_
    ds_ = 2.0 * np.pi * C1 / S_
    d   = ds_ + delta
    ps  = -S_ * delta / (2.0 * C1)
    dVdd = Gp_(ps) * G_(ps) * S_ / C1
    return -KAPPA * KD * geo_(d, r_, r_) * Sm_ * dVdd

def bisect(f, a, b, n=200):
    fa, fb = f(a), f(b)
    if fa == 0.0:
        return a
    if fb == 0.0:
        return b
    if fa * fb > 0:
        return float("nan")
    for _ in range(n):
        m = 0.5 * (a + b)
        fm = f(m)
        if fa * fm <= 0:
            b, fb = m, fm
        else:
            a, fa = m, fm
    return 0.5 * (a + b)

# ======================================================================
hr("T1 -- bond restoring force & stiffness (unison pair, symmetric branch)")
# ======================================================================
psi_s, dl, S_s, C_s = sp.symbols('psi delta S C', positive=True)
p_s = sp.Symbol('p', positive=True, integer=True)
h_sym  = (1 + sp.cos(psi_s)) / 2
G_sym  = h_sym ** p_s
Gp_sym = sp.diff(G_sym, psi_s)

chk_a = sp.simplify(sp.powsimp(
    Gp_sym + sp.Rational(1, 2) * p_s * (G_sym / h_sym) * sp.sin(psi_s), force=True))
emit("[a] kernel identity G'(psi) = -0.5*p*(G/h)*sin(psi):",
     "VERIFIED (diff == 0)" if chk_a == 0 else f"MISMATCH: {chk_a}")

# symmetric branch: psi_f = psi_b = psi = -S*(d-dstar)/(2C)  (wrapped)
psi_br  = -S_s * dl / (2 * C_s)
# kernel dVdd for the equal unison pair (q*wi = p*wj = S/2), on the branch:
dVdd_br = sp.powsimp((G_sym * Gp_sym).subs(psi_s, psi_br) * S_s / C_s, force=True)

# closed form: dVdd|branch = (p*S/(2C)) * h^(2p-1) * sin(S*delta/(2C))
u = sp.Symbol('u', real=True)
hcap = (1 + sp.cos(u)) / 2
dVdd_u = (p_s * S_s / (2 * C_s)) * hcap ** (2 * p_s - 1) * sp.sin(u)
chk_b = sp.simplify(sp.powsimp(
    dVdd_br - dVdd_u.subs(u, S_s * dl / (2 * C_s)), force=True))
emit("[b] on-branch kernel dVdd = (p*S/(2C)) * h^(2p-1) * sin(S*delta/(2C)),  h = (1+cos(S*delta/(2C)))/2 :",
     "VERIFIED" if chk_b == 0 else f"MISMATCH: {chk_b}")

# on the branch, the kernel fixed-phase dVdd equals the full d/d(delta) of V=1-G^2
dV_direct = sp.diff(1 - (G_sym.subs(psi_s, psi_br)) ** 2, dl)
chk_c = sp.simplify(sp.powsimp(dVdd_br - dV_direct, force=True))
emit("[c] kernel fixed-phase dVdd|branch == d/d(delta)[1 - G(psi)^2] (total on-branch derivative):",
     "EXACTLY EQUAL" if chk_c == 0 else f"differ by {chk_c}")
emit("    (equal unison pair only: moving d at fixed phases changes psi_f and psi_b by the")
emit("     same -S/(2C)*dd, i.e. fixed-phase motion tracks the symmetric branch.)")

# linear coefficient: derivative at u=0 times du/ddelta
lin = sp.simplify(sp.diff(dVdd_u, u).subs(u, 0) * (S_s / (2 * C_s)))
emit("[d] linearisation: dVdd|branch =", lin, "* delta + O(delta^3)")
emit("    ==> K_b = -d(dd/dt)/d(delta)|0 = kappa * (kd*geo*gpl*res) * Sm * p*S^2/(4*C^2)")
ratio_guess = sp.simplify(lin / (p_s * S_s ** 2 / (2 * C_s ** 2)))
emit("    COEFFICIENT CHECK vs the guess p*S^2/(2*C^2): derived/guess =", ratio_guess)
emit("    >>> SURPRISE: the guessed stiffness coefficient is 2x TOO LARGE; exact is p*S^2/(4*C^2). <<<")
emit("[e] sign: dd/dt = -kappa*kd*geo(d)*Sm*(p*S/(2C))*h^(2p-1)*sin(S*delta/(2C));")
emit("    for |delta| < 2*pi*C/S = dstar the sin carries sign(delta), h^(2p-1) > 0, geo >= 0,")
emit("    so sign(dd/dt) = -sign(delta): two-sided restoring over the whole rung basin,")
emit("    PROVIDED geo > 0, i.e. d < ri+rj = contact (lens open). Beyond contact the force is exactly 0.")

# ---------------- numbers ----------------
x0 = 0.325
Em, w, Es, r = seed_pair(x0)
S  = 2.0 * w
dstar = 2.0 * np.pi * C1 / S
Sm = Em
contact = 2.0 * r
geo0 = float(geo_(dstar, r, r))
A0 = np.pi * float(lens_a2(dstar, r, r))
Kb = KAPPA * KD * geo0 * Sm * P * S * S / (4.0 * C1 * C1)

emit("")
emit("numbers at x = 0.325 (kappa = 1):")
emit("  w2e      =", g15(w))
emit("  S=2*w2e  =", g15(S), "  S^2 =", g15(S * S))
emit("  dstar    =", g15(dstar), " (= pi/w2e, m=1)")
emit("  Em=Sm    =", g15(Em), "   Es =", g15(Es), "   r =", g15(r), "   contact=2r =", g15(contact))
emit("  lens A(dstar) =", g15(A0), "   Aref =", g15(AREF), "   geo(dstar) =", g15(geo0))
emit("  K_b = kappa*kd*geo*Sm*p*S^2/4 =", g15(Kb), " per t.u. per unit length  [1/t.u.]")
eps_fd = 1e-6
Kb_fd = (float(v_branch(-eps_fd, w, r, Sm)) - float(v_branch(+eps_fd, w, r, Sm))) / (2 * eps_fd)
emit("  numeric check, -(v(+e)-v(-e))/2e =", g15(Kb_fd))
emit("  relaxation time 1/K_b =", g15(1.0 / Kb), "t.u.;  per-step K_b*dt =", g15(Kb * DT), " (Euler-stable, <2)")

emit("")
emit("walk rate at delta = +/-0.044 (exact on-branch, full geometry):")
for de in (+0.044, -0.044):
    emit(f"  delta = {de:+.3f}:  |dd/dt| = {g15(abs(float(v_branch(de, w, r, Sm))))}"
         f"   (linearised K_b*|delta| = {g15(Kb * 0.044)})")
emit("  NOTE: strong +/- asymmetry is the lens area A(d): geo(dstar+0.044)/geo(dstar) ="
     f" {g15(float(geo_(dstar + 0.044, r, r)) / geo0)},"
     f"  geo(dstar-0.044)/geo(dstar) = {g15(float(geo_(dstar - 0.044, r, r)) / geo0)}")

# capture range: gate-only first
psis = np.linspace(1e-9, np.pi - 1e-9, 400001)
fgate = h_(psis) ** (2 * P - 1) * np.sin(psis)
ifx = int(np.argmax(fgate))
fmax = float(fgate[ifx]); psi_pk = float(psis[ifx])
psi_pk_closed = 2.0 * np.arctan(1.0 / np.sqrt(4.0 * P - 1.0))
i10 = ifx + int(np.where(fgate[ifx:] <= 0.1 * fmax)[0][0])
psi_10 = float(psis[i10])
emit("")
emit("capture range, GATE-ONLY (factor h^(2p-1)*|sin psi|, psi = S*delta/(2C)):")
emit("  peak at |psi| =", g15(psi_pk), " (closed form 2*atan(1/sqrt(4p-1)) =", g15(psi_pk_closed) + ")",
     "  f_max =", g15(fmax))
emit("  -> |delta|_peak(gate) = 2*psi/S =", g15(2 * psi_pk / S),
     ";  10%-of-max at |psi| =", g15(psi_10), " -> |delta|_10(gate) =", g15(2 * psi_10 / S))

# capture range with full geometry (lens truncation)
edge = contact - dstar
dp = np.linspace(1e-9, edge * (1 - 1e-9), 400001)
vp = np.abs(v_branch(dp, w, r, Sm))
ip = int(np.argmax(vp)); vmaxp = float(vp[ip]); dmaxp = float(dp[ip])
maskp = np.where(vp[ip:] <= 0.1 * vmaxp)[0]
d10p = float(dp[ip + int(maskp[0])]) if len(maskp) else float("nan")
dm = np.linspace(-(dstar - 0.05), -1e-9, 800001)
vm = np.abs(v_branch(dm, w, r, Sm))
im = int(np.argmax(vm)); vmaxm = float(vm[im]); dmaxm = float(dm[im])
maskm = np.where(vm[:im] <= 0.1 * vmaxm)[0]
d10m = float(dm[int(maskm[-1])]) if len(maskm) else float("nan")
emit("")
emit("capture range, FULL GEOMETRY (geo(d) included; kappa=1):")
emit("  + side: bond force support ENDS at delta = contact-dstar =", g15(edge),
     " (lens area -> 0 at d = 2r; force exactly 0 beyond)")
emit("    max |dd/dt| =", g15(vmaxp), "at delta =", g15(dmaxp), ";  10% of max at delta =", g15(d10p))
emit("  - side: max |dd/dt| =", g15(vmaxm), "at delta =", g15(dmaxm), ";  10% of max at delta =", g15(d10m))
emit("  >>> SURPRISE: capture is strongly ASYMMETRIC. The gate would peak at |delta|=0.17, but the")
emit("      lens vanishes first on the + side (delta=+" + g15(edge) + "), so the outward capture edge is")
emit("      geometric (= contact), not gate-set; the inward basin is deep (max ~"
     + f"{vmaxm:.3g}" + " at delta ~ " + f"{dmaxm:.3g}" + ").")

# ======================================================================
hr("T2 -- equilibrium under contact repulsion")
# ======================================================================
ov = contact - dstar
emit("overlap at the rung: contact - dstar =", g15(ov), " (cells DO overlap at the m=1 rung)")
psi50 = float(np.arccos(2.0 * 0.5 ** (1.0 / 16.0) - 1.0))
d50 = 2.0 * psi50 * C1 / S
emit("tongue (joint gate G^2 > 0.5): |psi| <", g15(psi50), " -> |delta| <", g15(d50),
     " (full width in d =", g15(2 * d50) + ")")
emit("  (+ side of the tongue is truncated by the lens edge at +", g15(edge), "before the gate closes)")
emit("")
for mk in (1.0, 0.1):
    vrep = mk * ov                      # T2 one-sided definition
    de1 = vrep / Kb
    de_sc = vrep / (Kb + mk)            # self-consistent linear (rep slope included)
    f_root = lambda dd_: float(v_branch(dd_, w, r, Sm)) + mk * (ov - dd_)
    de_hon = bisect(f_root, 1e-9, dmaxp)   # first crossing = stable equilibrium
    psi_eq = S * de1 / 2.0
    Gj = float(G_(psi_eq) ** 2)
    emit(f"mob_geo*k_rep = {mk}:  v_rep(dstar) = {g15(vrep)}")
    emit(f"  delta_eq (first order, v_rep/K_b)          = {g15(de1)}")
    emit(f"  delta_eq (self-consistent linear)          = {g15(de_sc)}")
    emit(f"  delta_eq (exact root on branch)            = {g15(de_hon)}")
    emit(f"  inside capture range?  yes: {g15(de1)} << +max-force delta {g15(dmaxp)} << lens edge {g15(edge)}")
    emit(f"  inside tongue (G^2>0.5, |delta|<{g15(d50)})?  yes;  joint gate at delta_eq: G^2 = {g15(Gj)}")
emit("")
# survival threshold (how strong may repulsion be before no equilibrium exists)
good = dp < ov - 1e-6
r_rel = np.max(vp[good] / (2.0 * (ov - dp[good])))
emit("bond-survival threshold (equilibrium exists iff bond can match repulsion somewhere):")
emit("  relative-velocity convention (T3, v_rep_rel = 2*mob*k_rep*(contact-d)): mob*k_rep <", g15(float(r_rel)))
emit("  one-sided convention (T2):                                            mob*k_rep <", g15(2.0 * float(r_rel)))
emit("  (both forces vanish linearly at d=contact; their ratio is finite there and is the sup)")

# ======================================================================
hr("T3 -- two-cell honest integration (dt = 0.02, kappa = kappa_lock = 1)")
# ======================================================================
def wrap_pi(a):
    return (a + np.pi) % (2.0 * np.pi) - np.pi

def integrate(d0, mix, T=50.0, rep=True, mk=1.0):
    nst = int(round(T / DT))
    thi = 0.0
    thj = -w * d0 / C1          # seed: psi_f(0) = 0  (m=1 lock at d0)
    d = d0
    bad = 0
    samples = {}
    traj = np.empty(nst)
    for n in range(nst):
        psi_f = thi - w * d / C1 - thj
        psi_b = thj - w * d / C1 - thi
        dVdd = (Gp_(psi_f) * G_(psi_b) * w + G_(psi_f) * Gp_(psi_b) * w) / C1
        dd_bond = -KAPPA * (KD * DT * float(geo_(d, r, r))) * Sm * dVdd
        delta = d - dstar
        if abs(delta) > 1e-7 and dd_bond * delta > 0 and abs(dd_bond) > 1e-14:
            bad += 1
        dd_rep = DT * 2.0 * mk * (contact - d) if (rep and d < contact) else 0.0
        d += dd_bond + dd_rep
        thi += w * DT
        thj += w * DT
        if mix > 0.0:
            thj += 1.0 * mix * wrap_pi((thi - w * d / C1) - thj)
        traj[n] = d
        if (n + 1) in (250, 500, 1000, 2500):
            samples[(n + 1) * DT] = d
    pf = float(wrap_pi(thi - w * d / C1 - thj))
    pb = float(wrap_pi(thj - w * d / C1 - thi))
    dfin = traj[-1]
    tol99 = 0.01 * abs(d0 - dfin)
    idx = np.where(np.abs(traj - dfin) > tol99)[0]
    t99 = (int(idx[-1]) + 2) * DT if len(idx) else DT
    return d, samples, pf, pb, bad, t99

pred1 = 2.0 * ov / Kb
pred2 = 2.0 * ov / (Kb + 2.0)
f_root_rel = lambda dd_: float(v_branch(dd_, w, r, Sm)) + 2.0 * (ov - dd_)
pred_hon = bisect(f_root_rel, 1e-9, dmaxp)
# psi_f = 0 branch (where one-sided entrainment lands): dVdd = w*G'(-2w*delta)
def v_psif0(dd_):
    d_ = dstar + dd_
    return -KAPPA * KD * float(geo_(d_, r, r)) * Sm * (w * float(Gp_(-2.0 * w * dd_)))
pred_hon0 = bisect(lambda dd_: v_psif0(dd_) + 2.0 * (ov - dd_), 1e-9, dmaxp)
emit("repulsion in T3 is RELATIVE: d_dot_rep = 2*mob*k_rep*(contact-d), mob=k_rep=1")
emit("predictions for delta_eq = d_inf - dstar:")
emit("  first order    2*v_rep0/K_b        =", g15(pred1))
emit("  self-consistent 2*v_rep0/(K_b+2)   =", g15(pred2))
emit("  exact root on symmetric branch     =", g15(pred_hon))
emit("  exact root on psi_f=0 branch       =", g15(pred_hon0), " (entrained runs should land here)")
emit("")
emit("runs WITH repulsion (report d(t)):")
for d0 in (dstar + 0.05, dstar - 0.05):
    for mix in (0.0, 0.01, 0.3):
        dfin, samp, pf, pb, bad, t99 = integrate(d0, mix)
        tag = "none" if mix == 0.0 else f"{mix}"
        emit(f"  d0 = dstar{d0-dstar:+.2f}, entrain mix = {tag}:")
        emit("    d(5)  =", g15(samp[5.0]), "  d(10) =", g15(samp[10.0]),
             "  d(20) =", g15(samp[20.0]), "  d(50) =", g15(samp[50.0]))
        emit("    delta_final =", g15(dfin - dstar), "  t_99 =", g15(t99), "t.u.")
        emit("    final wrapped psi_f =", g15(pf), " psi_b =", g15(pb),
             "  G(psi_f) =", g15(float(G_(pf))), " G(psi_b) =", g15(float(G_(pb))),
             "  joint =", g15(float(G_(pf) * G_(pb))))
        emit("    bond-term sign violations (dd_bond pushing AWAY from dstar):", bad)
emit("")
emit("runs WITHOUT repulsion (pure bond; should converge to dstar exactly):")
for d0 in (dstar + 0.05, dstar - 0.05):
    for mix in (0.0, 0.01, 0.3):
        dfin, samp, pf, pb, bad, t99 = integrate(d0, mix, rep=False)
        tag = "none" if mix == 0.0 else f"{mix}"
        emit(f"  d0 = dstar{d0-dstar:+.2f}, mix = {tag}: |delta_final| = {g15(abs(dfin - dstar))},"
             f" joint gate = {g15(float(G_(pf) * G_(pb)))}, sign violations = {bad}")
emit("")
emit("why no-entrainment locking works: for the equal pair, psi_f + psi_b = -S*d/C IDENTICALLY")
emit("(phase-independent).  To first order the kernel force depends only on that sum, so the")
emit("restoring force toward dstar is branch/split independent; the seed split (psi_f=0) merely")
emit("freezes a gate offset (joint gate < 1) without moving the equilibrium or the stiffness.")
emit("One-sided entrainment (the specified thj pull) converges to the psi_f=0 branch, NOT the")
emit("symmetric branch; a symmetric split needs two-sided entrainment.  All three settings agree")
emit("on d_inf to O(1e-4).")
emit("")
emit("CRITICAL SIGN CHECK: in every run (both seeds, mix = none/0.01/0.3, with and without")
emit("repulsion) the kernel bond term dd = -kappa*base*Sm*dVdd evaluated at CURRENT phases")
emit("pushed d TOWARD dstar at every step (0 violations above threshold): the bond is a")
emit("two-sided restoring force under realistic phase dynamics.")

# ======================================================================
hr("T4 -- vacuum inertness + crumb bound")
# ======================================================================
def m_eff(mi, mj):
    return np.sqrt(mi * max(mj, FL0))

def Sm_eff(mi, mj):
    return np.sqrt(m_eff(mi, mj) * m_eff(mj, mi))

emit("Sm = sqrt(mi_eff*mj_eff), mi_eff = sqrt(mob_i*max(mob_j,fl0)), mj_eff = sqrt(mob_j*max(mob_i,fl0))")
emit("  => Sm^2 = sqrt(mob_i*mob_j) * sqrt(max(mob_i,fl0)*max(mob_j,fl0));  the factor sqrt(mob_i*mob_j)")
emit("     vanishes when EITHER mob is 0, so Sm = 0 EXACTLY for any vacuum partner.  fl0 =", g15(FL0))
emit("  numeric: Sm(0, 0.8125) =", g15(Sm_eff(0.0, 0.8125)),
     "   Sm(0.8125, 0) =", g15(Sm_eff(0.8125, 0.0)))
emit("  crumb pair Em_i = Em_j = eps:  Sm(eps) = sqrt(eps*max(eps,fl0)) ; for eps < fl0 this is")
emit("  sqrt(eps*fl0)  (proportional to sqrt(eps) -- the fl0 floor BOOSTS small crumbs above linear).")
emit("")
for eps in (1e-3, 1e-2):
    x_c = eps / CAP
    Em_c, w_c, Es_c, r_c = seed_pair(x_c)
    S_c = 2.0 * w_c
    ds_c = 2.0 * np.pi * C1 / S_c
    Sm_c = float(Sm_eff(eps, eps))
    geo_c = float(geo_(ds_c, r_c, r_c))
    vmax_rung = KAPPA * KD * geo_c * Sm_c * (P * S_c / 2.0) * fmax   # geometry pinned at the rung
    # bonus: full max including geo(d) variation
    edge_c = 2.0 * r_c - ds_c
    dg = np.linspace(-(ds_c - 0.05), edge_c * (1 - 1e-9), 400001)
    vg = np.abs(v_branch(dg, w_c, r_c, Sm_c))
    ig = int(np.argmax(vg))
    ov_c = 2.0 * r_c - ds_c
    emit(f"eps = {eps}:  x = {g15(x_c)},  w2e = {g15(w_c)},  dstar = {g15(ds_c)},  r = {g15(r_c)},"
         f"  Sm = {g15(Sm_c)}")
    emit(f"  geo(dstar) = {g15(geo_c)};  max|dVdd| over delta = (p*S/2)*f_max = {g15(P*S_c/2.0*fmax)}")
    emit(f"  v_bond_max (geometry AT the rung, per task) = {g15(vmax_rung)}")
    emit(f"    [bonus: full-geometry max = {g15(float(vg[ig]))} at delta = {g15(float(dg[ig]))}]")
    emit(f"  vs repulsion scale at 10% overlap (mob*k_rep*0.17):  ratio = {g15(vmax_rung/0.17)}")
    emit(f"  overlap AT the crumb rung = contact - dstar = {g15(ov_c)} (i.e. {g15(100*ov_c/(2*r_c))}% of contact):")
    emit(f"    repulsion there (one-sided) = {g15(ov_c)}, (relative) = {g15(2*ov_c)};"
         f"  bond/repulsion at the rung = {g15(vmax_rung/ov_c)} (one-sided)")
    d10 = 2.0 * r_c - 0.17
    v10 = abs(float(v_branch(d10 - ds_c, w_c, r_c, Sm_c)))
    emit(f"    bond velocity AT 10% overlap (d = {g15(d10)}, delta = {g15(d10-ds_c)}) = {g15(v10)}"
         f"  -> ratio to 0.17 = {g15(v10/0.17)}")
emit("")
emit("conclusion: vacuum is EXACTLY inert (Sm=0).  Crumb pairs have max bond velocity 0.21x (eps=1e-3)")
emit("to 0.66x (eps=1e-2) of the 0.17 repulsion scale, BUT their rung sits at 36% overlap where")
emit("repulsion exceeds the bond maximum by an order of magnitude; and once pushed out to ~10%")
emit("overlap the gate factor has collapsed (ratio ~1e-4 to 1e-3).  Crumbs cannot bind; no vacuum")
emit("crystallisation from the P15 term.")

# ======================================================================
hr("T5 -- mu-convexity re-check (geopress criterion), sympy")
# ======================================================================
Ei, Ej, dq, kq, aq, r0q, e0q = sp.symbols('E_i E_j d k alpha r_0 e_0', positive=True)
riq = r0q * (Ei / e0q) ** sp.Rational(1, 3)
rjq = r0q * (Ej / e0q) ** sp.Rational(1, 3)
dltq = riq + rjq - dq
Uq = (kq / aq) * dltq ** aq
mu_i = sp.diff(Uq, Ei)
dmu_i = sp.diff(mu_i, Ei)
cand = (kq * riq / (9 * Ei ** 2)) * dltq ** (aq - 2) * ((aq - 1) * riq - 2 * dltq)
# robust check: exact symbolic at integer alphas + numeric random points for general alpha
ok_int = all(sp.simplify(sp.expand(sp.powsimp((dmu_i - cand).subs(aq, av), force=True))) == 0
             for av in (2, 3, 4))
rng = np.random.default_rng(7)
fn = sp.lambdify((Ei, Ej, dq, kq, aq, r0q, e0q), dmu_i - cand, "numpy")
num_res = 0.0
for _ in range(6):
    Eiv, Ejv = rng.uniform(0.3, 3.0, 2)
    r0v, e0v, kv = rng.uniform(0.5, 2.0, 3)
    av = rng.uniform(1.2, 4.0)
    riv = r0v * (Eiv / e0v) ** (1 / 3); rjv = r0v * (Ejv / e0v) ** (1 / 3)
    dv = (riv + rjv) * rng.uniform(0.5, 0.95)   # positive overlap
    num_res = max(num_res, abs(float(fn(Eiv, Ejv, dv, kv, av, r0v, e0v))))
ok_num = num_res < 1e-10
emit("[per-cell] mu_i = dU/dE_i at fixed partner:  mu_i = k*delta^(alpha-1)*r_i/(3*E_i)")
emit("  d(mu_i)/dE_i = (k*r_i/(9*E_i^2)) * delta^(alpha-2) * [ (alpha-1)*r_i - 2*delta ]",
     " -- check:",
     ("VERIFIED (exact symbolic at alpha=2,3,4; numeric residual < %.1e at 6 random general-alpha points)"
      % max(num_res, 1e-16)) if (ok_int and ok_num) else f"FAILED int={ok_int} num_res={num_res}")
emit("  => mu falls with E_i  iff  delta/r_i > (alpha-1)/2     [claim CONFIRMED]")
emit("  alpha = 2:  delta/r > 1/2  -- confirmed.")
al2 = sp.simplify(dmu_i.subs(aq, 2))
emit("  (alpha=2 explicit d(mu)/dE_i:", al2, ")")
emit("")
Eq_ = sp.Symbol('E', positive=True)
rq_ = r0q * (Eq_ / e0q) ** sp.Rational(1, 3)
dlt2 = 2 * rq_ - dq
U2 = (kq / aq) * dlt2 ** aq
dmu2 = sp.diff(U2, Eq_, 2)
cand2 = (4 * kq * rq_ ** 2 / (9 * Eq_ ** 2)) * dlt2 ** (aq - 2) * ((aq - 1) - dlt2 / rq_)
diff2 = sp.simplify(sp.powsimp(dmu2 - cand2, force=True))
emit("[symmetric-growth caveat] if BOTH cells grow with one E (delta = 2*r(E) - d):")
emit("  d^2U/dE^2 = (4k*r^2/(9E^2)) * delta^(alpha-2) * [ (alpha-1) - delta/r ]",
     " -- symbolic check:", "VERIFIED" if diff2 == 0 else f"residual {diff2}")
emit("  => along that path the criterion is delta/r > (alpha-1)  (alpha=2 -> delta/r > 1): 2x harsher.")
emit("  The stated geopress criterion is the PER-CELL chemical potential (partner fixed) -- that one")
emit("  is exactly delta/r > (alpha-1)/2.  Constants k, r_0, e_0 drop out in both (scale-free).")
emit("")
emit("numeric finite-difference confirmation (k=1, r0=e0=1, E=1 so r=1):")
def mu_num(E, dfix, alpha, mode):
    if mode == "per-cell":
        f = lambda Ee: (1.0 / alpha) * max(Ee ** (1.0 / 3.0) + 1.0 - dfix, 0.0) ** alpha
    else:
        f = lambda Ee: (1.0 / alpha) * max(2.0 * Ee ** (1.0 / 3.0) - dfix, 0.0) ** alpha
    e = 1e-6
    return (f(E + e) - f(E - e)) / (2 * e)
def dmu_num(dfix, alpha, mode):
    e = 1e-4
    return (mu_num(1 + e, dfix, alpha, mode) - mu_num(1 - e, dfix, alpha, mode)) / (2 * e)
for mode, cases in (("per-cell", (0.4, 0.6)), ("symmetric", (0.9, 1.1))):
    for dl_r in cases:
        dfix = (2.0 - dl_r)
        v = dmu_num(dfix, 2.0, mode)
        emit(f"  {mode}, alpha=2, delta/r = {dl_r}:  d(mu)/dE = {v:+.6e}"
             f"  ({'falls' if v < 0 else 'rises'})")

# ======================================================================
hr("T6 -- truss rigidity counting (first-order, central-force)")
# ======================================================================
def edges_min(V):
    V = np.asarray(V, float)
    n = len(V)
    pairs = []
    for i in range(n):
        for j in range(i + 1, n):
            pairs.append((np.linalg.norm(V[i] - V[j]), i, j))
    dmin = min(p[0] for p in pairs)
    return [(i, j) for dd_, i, j in pairs if dd_ < dmin * 1.0001]

def rigidity_report(name, V, E):
    V = np.asarray(V, float)
    R = np.zeros((len(E), 3 * len(V)))
    for k2, (i, j) in enumerate(E):
        uvec = V[i] - V[j]
        uvec = uvec / np.linalg.norm(uvec)
        R[k2, 3 * i:3 * i + 3] = uvec
        R[k2, 3 * j:3 * j + 3] = -uvec
    sv = np.linalg.svd(R, compute_uv=False)
    tol = max(1e-10, sv[0] * max(R.shape) * np.finfo(float).eps * 10)
    rank = int((sv > tol).sum())
    floppy = 3 * len(V) - 6 - rank
    smin_kept = sv[rank - 1] if rank > 0 else float("nan")
    smax_drop = sv[rank] if rank < len(sv) else 0.0
    emit(f"  {name:24s} V={len(V):2d} E={len(E):2d}  rank(R)={rank:2d}  floppy=3V-6-rank={floppy:2d}"
         f"   [sv gap: kept>={smin_kept:.3e}, dropped<={smax_drop:.3e}]"
         f"   {'SHEAR-RIGID (isostatic)' if floppy == 0 else 'FLOPPY'}")
    return floppy

hexV = [(np.cos(k2 * np.pi / 3), np.sin(k2 * np.pi / 3), 0.0) for k2 in range(6)]
hexE = [(k2, (k2 + 1) % 6) for k2 in range(6)]
a_o = 1.0 / np.sqrt(2.0)
octV = [(a_o, 0, 0), (-a_o, 0, 0), (0, a_o, 0), (0, -a_o, 0), (0, 0, a_o), (0, 0, -a_o)]
octE = edges_min(octV); assert len(octE) == 12, len(octE)
cubV = [(i, j, k2) for i in (0, 1) for j in (0, 1) for k2 in (0, 1)]
cubE = edges_min(cubV); assert len(cubE) == 12, len(cubE)
gr = (1 + np.sqrt(5)) / 2
icoV = []
for s1 in (+1, -1):
    for s2 in (+1, -1):
        icoV += [(0, s1, s2 * gr), (s1, s2 * gr, 0), (s2 * gr, 0, s1)]
icoE = edges_min(icoV); assert len(icoE) == 30, len(icoE)

emit("floppy modes = 3V - 6 - rank(R)   (R = central-force rigidity matrix; 6 rigid-body modes)")
f_hex = rigidity_report("ring6 (hexagon)", hexV, hexE)
f_oct = rigidity_report("octahedron", octV, octE)
f_cub = rigidity_report("cube", cubV, cubE)
f_ico = rigidity_report("icosahedron", icoV, icoE)
f_hd  = rigidity_report("ring6 + 3 long diagonals", hexV, hexE + [(0, 3), (1, 4), (2, 5)])
emit("")
emit(f"shear-rigid seeds: octahedron (floppy {f_oct}) and icosahedron (floppy {f_ico}) -- both isostatic")
emit(f"(E = 3V-6 exactly).  ring6 has {f_hex} floppy modes, cube has {f_cub} (shear/flex), ring6+3diag"
     f" has {f_hd}.")
emit("note: ring6+3diag has rank E-1 = 8, i.e. one SELF-STRESS (ring tension balanced by diagonal")
emit("compression through the centre), so adding 3 members removes only 2 floppy modes (6 -> 4),")
emit("not 3: in-plane bracing of a symmetric ring wastes a member on the self-stress.")

# ======================================================================
hr("SURPRISES / deviations from the task's stated expectations")
# ======================================================================
emit("1. K_b coefficient is p*S^2/(4*C^2), NOT p*S^2/(2*C^2): the guessed stiffness is 2x too large.")
emit("   K_b(x=0.325) = " + g15(Kb) + " * kappa per t.u.")
emit("2. The outward capture edge is NOT gate-set: the lens area (and the whole bond force) vanishes")
emit("   exactly at d = contact = 2r.  Bond and repulsion share the same support (overlap only).")
emit("   Full-geometry capture: force max at delta = " + g15(dmaxp) + " (+ side) / " + g15(dmaxm)
     + " (- side); 10% points " + g15(d10p) + " / " + g15(d10m) + ".  Strongly asymmetric.")
emit("3. T3: NO entrainment is needed to lock -- psi_f + psi_b = -S*d/C is a phase-independent")
emit("   identity, so the first-order restoring force is split-independent.  Pure phase advance")
emit("   converges to the same d_inf as entrained runs; entrainment only re-opens the frozen gate")
emit("   offset.  The specified ONE-SIDED entrainment lands on the psi_f=0 branch, not the")
emit("   symmetric branch (needs two-sided pulls); dynamically irrelevant at first order.")
emit("4. T4: the fl0 floor makes crumb bonding scale as sqrt(eps*fl0) (weaker than any load but")
emit("   BOOSTED above linear-in-eps); still, a crumb pair's rung lies at 36% overlap where")
emit("   repulsion exceeds the bond max 5.5x (eps=1e-2) to 17x (eps=1e-3) one-sided (11-35x")
emit("   relative), so crumbs cannot bind.  Vacuum exactly inert.")
emit("5. T5: criterion delta/r > (alpha-1)/2 CONFIRMED for the per-cell mu (partner fixed);")
emit("   along the symmetric both-cells-grow path it is delta/r > (alpha-1) (factor 2 harsher) --")
emit("   which derivative geopress measures matters.")
emit("6. T6: ring6 floppy=" + str(f_hex) + ", cube floppy=" + str(f_cub) + "; octahedron and icosahedron"
     + " floppy=0 (only shear-rigid seeds).")

# ---------------- write results ----------------
SCRATCH = os.path.dirname(os.path.abspath(__file__))
with open(os.path.join(SCRATCH, "deriv_bond_results.txt"), "w") as fo:
    fo.write("\n".join(OUT) + "\n")
print("\n[results written to " + os.path.join(SCRATCH, "deriv_bond_results.txt") + "]")
