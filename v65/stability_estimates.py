#!/usr/bin/env python3
"""
v65 — Back-calculating a STABLE structure from theory (oscillon / phased-frequency).

Premise (from 8 failed self-tuning runs): no STATIC soliton is stable here. Derrick's
theorem says so. The escape is TIME-PERIODICITY: a phased frequency omega provides an
effective-mass barrier (m^2 - omega^2) that holds the lump open against collapse/dispersal.
This script estimates the stabilizing frequency, size, and amplitude so we can SEED a tuned
oscillon instead of searching blindly.

Constants (CLAUDE.md standard): m^2=2.25 (m=1.5), mu=-41.345, kappa=50, delta={0,3.0005,4.4325}.
The potential is V(P)=(mu/2)P^2/(1+kappa P^2), P=phi0*phi1*phi2 -> a QUINTIC restoring force.
The oscillon seed builds phi_a = A*cos(delta_a), so at the core P = A^3 * prod(cos delta_a).
"""
import math

M2, MU, KAPPA = 2.25, -41.345, 50.0
M = math.sqrt(M2)
DELTA = [0.0, 3.0005, 4.4325]
CPROD = math.cos(DELTA[0])*math.cos(DELTA[1])*math.cos(DELTA[2])  # P = A^3 * CPROD

print("="*70)
print("v65 — phased-frequency (oscillon) stability estimates")
print("="*70)
print(f"m = {M:.3f}  (mass gap; oscillon ceiling omega < m)")
print(f"phase product cos(d0)cos(d1)cos(d2) = {CPROD:+.4f}  => P_core = {abs(CPROD):.4f}*A^3")
print()

# ---- 1. Derrick: static is a saddle (quantify the instability) -------------
print("1. DERRICK SCALING  E(lambda) = lambda*E2 + lambda^3*E_V")
print("-"*70)
print("   Extremum: E2 + 3 E_V = 0 ;  d2E/dlambda2 = 6 E_V < 0  => MAXIMUM (unstable).")
print("   => no stable static soliton. Confirmed by every run (collapse OR dispersal).")
print("   The collapse mode is the lambda (breathing) direction; it must be lifted by")
print("   a frequency (oscillon) or a conserved charge (Q-ball).")
print()

# ---- 2. Oscillon size vs frequency: R(omega) = 1/sqrt(m^2 - omega^2) --------
print("2. OSCILLON SIZE vs FREQUENCY   R(omega) = 1/sqrt(m^2 - omega^2)")
print("   (the localized mode's exponential tail e^{-r/R}; sets the lump width)")
print("-"*70)
def R_of(omega): return 1.0/math.sqrt(M2 - omega*omega)
# ---- 3. Amplitude from core nonlinear balance ------------------------------
# Reduced oscillon balance: the nonlinear frequency pull must bridge the gap m^2-omega^2.
# Force on phi0 ~ |mu| * P^2 / [phi0 (1+kappa P^2)^2], P=CPROD*A^3 -> ~ |mu| g5 A^4 / sat,
# with g5 = CPROD^2 (the quintic coefficient) and sat=(1+kappa P^2)^2 the saturation.
# Frequency pull ~ that /A. Set m^2-omega^2 = |mu|*g5*A^4 / sat  (self-consistent in A).
G5 = CPROD*CPROD
def amplitude_for(omega):
    # Unsaturated core balance: m^2-omega^2 = |mu|*g5*A^4  =>  A=(gap/(|mu| g5))^(1/4).
    # (The low-amplitude oscillon root; saturation kappa*P^2 reported as a caveat — when
    #  it is O(1) the estimate is only order-of-magnitude.)
    gap = M2 - omega*omega
    A = (gap / (abs(MU)*G5))**0.25
    P = abs(CPROD)*A**3
    return A, P, KAPPA*P*P

print(f"   {'omega':>6} {'period':>7} {'R(code)':>8} {'R(fm)':>7} {'A_amp':>7} "
      f"{'P_core':>8} {'kP^2':>7}")
CODE_L_FM = 0.5624
cands = []
for omega in [0.9, 1.1, 1.25, 1.35, 1.42, 1.46]:
    R = R_of(omega); A, P, kP2 = amplitude_for(omega)
    per = 2*math.pi/omega
    print(f"   {omega:>6.2f} {per:>7.2f} {R:>8.2f} {R*CODE_L_FM:>7.2f} {A:>7.3f} "
          f"{P:>8.4f} {kP2:>7.3f}")
    cands.append((omega, R, A))
print()

# ---- 4. The theta-sector radiation constraint (phase-locking) --------------
print("4. THETA SECTOR (massless) — the radiation constraint")
print("-"*70)
print("   phi disperson: omega^2 = k^2 + m^2 (massive, bound for omega<m).")
print("   theta dispersion: omega^2 = k^2 (massless). A phi core oscillating at omega")
print("   drives theta at the SAME omega -> k_theta = omega -> theta RADIATES (propagates).")
print(f"   That radiation is a loss channel (the eta*curl coupling, eta=0.5). To minimize it:")
print("   - lower omega (slower drive) reduces theta production, OR")
print("   - exploit the phase offsets so the theta source's net (curl) cancels at low multipole.")
print("   => prefer the LOWER end of the oscillon band (omega ~ 0.9-1.25) for longevity,")
print("      trading width (larger R) for less theta radiation.")
print()

# ---- 5. Recommendation -----------------------------------------------------
print("5. SEED RECOMMENDATIONS (back-calculated tuned oscillons to test)")
print("-"*70)
print("   Our previous seed (A=1.0, sigma=1.5) was OVER-DRIVEN: at sigma=1.5 the matched")
print(f"   frequency is omega~{math.sqrt(M2-1/1.5**2):.2f} (R=1.5), whose tuned amplitude is")
A15,_,_ = amplitude_for(math.sqrt(M2-1/1.5**2))
print(f"   A~{A15:.2f}, NOT 1.0 -> the excess radiated away (the dispersal we saw).")
print("   Tuned candidates to seed (gen_oscillon -A <A> -sigma <R>) and test for COHERENCE")
print("   (P_max stays ~steady & oscillates at the predicted period, instead of decaying):")
for omega, R, A in cands:
    if 1.0 <= omega <= 1.4:
        print(f"     omega={omega:.2f}: -sigma {R:.2f} -A {A:.2f}   (period {2*math.pi/omega:.1f} t.u.)")
print()
print("   Then: launch with a TIME phase (nonzero phi velocity ~ -omega*A*sin) so it starts")
print("   on the oscillon limit cycle, not from rest. A coherent persistent oscillon is the")
print("   stable substrate the self-tuning needs.")
print("="*70)

# self-checks
def ck(n,c):
    print(f"  [{'PASS' if c else 'FAIL'}] {n}")
    return c
print("\nSELF-CHECKS\n"+"-"*11)
ok=True
ok&=ck("R increases as omega->m (wider, shallower)", R_of(1.4) > R_of(0.9))
ok&=ck("tuned amplitude decreases toward the ceiling omega->m", amplitude_for(1.42)[0] < amplitude_for(0.9)[0])
ok&=ck("matched-A for sigma=1.5 is below our seed's A=1.0", A15 < 1.0)
ok&=ck("oscillon band below mass gap", R_of(1.46) > 0 and M2-1.46**2 > 0)
import sys
print("\n"+("All self-checks passed." if ok else "SELF-CHECKS FAILED"))
sys.exit(0 if ok else 1)
