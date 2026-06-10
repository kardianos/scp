#!/usr/bin/env python3
"""v67 theta_dynamics WP-A: bulk theta as medium -> relative mass. Numbers.

Symbolic basis: theta_frame.mac (28/28 PASS). Conventions v66/THEORY.md.
Standard parameters m2=2.25, eta=0.5, mu=-41.345, kappa=50, m_theta=0.

Pieces evaluated here:
 (1) vacuum polariton (V50 limit check): v_ph = sqrt(1-eta^2/m2) = 0.9428
 (2) bath-induced ambient mass shift [thm for coherent Beltrami, F7/F8]:
       dm2(K) = mu eta^4 K^4 B^4 / M2t^4,  M2t = m2 + K^2 - W^2,
       e_bath = (W^2+K^2) B^2,  W = photon-branch polariton frequency W(K).
     Charge-neutral incoherent (v67 bath): Gaussian mean-field
       dm2_neutral ~ mu * n_comp^2  (per-component driven variance n_comp),
     isotropic in components [estimate]; single-mode neutral = 1/2 charged (F12).
 (3) driven Phi variance: n_tot(K) = 2 eta^2 K^2 e / ((W^2+K^2) M2t^2^2) ... see below
 (4) Q-ball core mean-field shift [estimate]:
       d(omega^2) = g'(f0^2) * n_comp,  g(X) = mu X^2/(1+kappa X^3)^2  (F10),
       f0 = f(0) = 0.6405 at omega = 1.39 (v66 table).
 (5) effective-metric / photon-speed corrections, boost anisotropy.
"""
import numpy as np

m2, eta, mu, kap = 2.25, 0.5, -41.345, 50.0
m = np.sqrt(m2)
omega_ball = 1.39
f0 = 0.6405          # core amplitude at omega=1.39 (v66 THEORY table)

def W2_photon(K):
    """photon-branch polariton frequency^2 at wavevector K (F2)."""
    return (2*K**2 + m2 - np.sqrt(m2**2 + 4*eta**2*K**2))/2

def coeffs(K, branch="photon"):
    W2 = W2_photon(K) if branch == "photon" else K**2
    M2t = m2 + K**2 - W2                    # detuning (Lorentz-invariant, F11 note)
    # per unit e_bath:  B^2 = e/(W^2+K^2)
    c0sq_per_e = eta**2*K**2/((W2+K**2)*M2t**2)   # |Phi_drv| per x,y component, per e
    n_tot_per_e = 2*c0sq_per_e                    # two driven components (F7 point)
    dm2_per_e2 = mu*eta**4*K**4/((W2+K**2)**2*M2t**4)  # coherent charged Beltrami (F8b)
    return W2, M2t, c0sq_per_e, n_tot_per_e, dm2_per_e2

print("=== (1) V50 limit check [verified F3] ===")
print(f"v_photon(k->0) = sqrt(1-eta^2/m2) = {np.sqrt(1-eta**2/m2):.4f} c  (V50: 0.9428c)")
print(f"matter-branch c_eff^2 = 1+eta^2/m2 = {1+eta**2/m2:.4f}  (transverse pol. only)")

print("\n=== (2) band scan K in [1.1,1.7] (photon-branch bath) ===")
Ks = np.linspace(1.1, 1.7, 13)
rows = []
for K in Ks:
    W2, M2t, c0e, nte, dme2 = coeffs(K)
    rows.append((K, np.sqrt(W2), M2t, c0e, nte, dme2))
print("  K     W(K)    M2t     c0^2/e   n_tot/e   dm2/e^2(coh)")
for r in rows[::3]:
    print(f"  {r[0]:.2f}  {r[1]:.4f}  {r[2]:.4f}  {r[3]:.5f}  {r[4]:.5f}  {r[5]:+.5f}")
band = np.array(rows)
C_coh   = band[:,5].mean()         # coherent charged dm2 coefficient (K^ component)
n_tot_e = band[:,4].mean()         # total driven |Phi|^2 per unit e
n_comp_e = n_tot_e/3               # isotropic bath: spread over 3 components
C_gauss = mu*n_comp_e**2           # neutral incoherent Gaussian mean-field [estimate]
print(f"band-avg: C_coh = dm2/e^2 = {C_coh:+.5f} (coherent, K^-component only)")
print(f"band-avg: n_tot/e = {n_tot_e:.5f}, n_comp/e = {n_comp_e:.5f}")
print(f"Gaussian neutral-bath: dm2 = mu*n_comp^2 = {C_gauss:+.5f} e^2 (isotropic) [estimate]")
print(f"single-mode neutral (F12c): {C_coh/2:+.5f} e^2 + Floquet at (4K,4W)")

print("\n=== (4) core mean-field slope ===")
X = f0**2
gp = 2*mu*X*(1-2*kap*X**3)/(1+kap*X**3)**3    # g'(X), F10b
print(f"f0={f0}, X=f0^2={X:.4f}, kappa*X^3={kap*X**3:.3f} (>1/2: saturated => g'>0)")
print(f"g'(X) = {gp:+.4f}   => d(omega^2)_core = g'(X)*n_comp = {gp*n_comp_e:+.5f} * e_bath")

print("\n=== (5) PREDICTION TABLE: v67 bath ladder, ball omega=1.39 ===")
print("dm2_amb = C_gauss*e^2 [estimate]; dm/m = dm2/(2 m2)")
print("d_omega_core = g'*n_comp*e/(2 omega) [estimate, factor ~2 syst]")
print("d_omega_amb  = dm2_amb/(2 omega)")
print(" e_bath | n_comp   | dm2_amb    | dm_eff/m  | dw_core   | dw_amb     | dw_tot   | omega_pred")
for e in (0.05, 0.2, 0.8):
    nc = n_comp_e*e
    dm2a = C_gauss*e**2
    dwc = gp*nc/(2*omega_ball)
    dwa = dm2a/(2*omega_ball)
    dwt = dwc + dwa
    print(f"  {e:4.2f}  | {nc:.6f} | {dm2a:+.6f}  | {dm2a/(2*m2):+.2e} | {dwc:+.6f} | {dwa:+.6f}  | {dwt:+.6f} | {omega_ball+dwt:.4f}")
print("(coherent/charged-bath variant: dm2_amb = C_coh*e^2 =",
      f"{C_coh*0.8**2:+.5f} at e=0.8, K^-component only)")

print("\n=== (6) effective metric / photon-speed corrections at e_bath=0.8 ===")
dm2 = C_gauss*0.8**2
vph2_0 = 1-eta**2/m2
vph2_b = 1-eta**2/(m2+dm2)
print(f"photon: v^2 = {vph2_0:.6f} -> {vph2_b:.6f}; d(v^2) = {vph2_b-vph2_0:+.2e}")
print(f"  (eta^2*dm2/m2^2 = {eta**2*dm2/m2**2:+.2e}, F9 agrees)")
print(f"matter g^ij_perp = (1+eta^2/(m2+dm2)) = {1+eta**2/(m2+dm2):.6f} (vac {1+eta**2/m2:.6f})")
print("coherent-bath birefringence (k perp K^): Delta(m_eff^2) between polarizations "
      f"= |C_coh| e^2 = {abs(C_coh)*0.8**2:.5f} at e=0.8")

print("\n=== (7) motion through bath: relative-mass anisotropy [estimate] ===")
print("isotropic on-shell bath, ball velocity v: e_bath -> gamma^2(1+v^2/3) e")
print("dm2(v) = dm2(0) * gamma^4 (1+v^2/3)^2  ~  dm2(0) * (1 + 8 v^2/3)  (O(v^2))")
for v in (0.1, 0.3, 0.5):
    g2 = 1/(1-v**2)
    fac = g2**2*(1+v**2/3)**2
    print(f"  v={v}: dm2(v)/dm2(0) = {fac:.4f}  (O(v^2) approx {1+8*v**2/3:.4f})")

print("\n=== consistency cross-checks ===")
K = 1.4; W2, M2t, c0e, nte, dme2 = coeffs(K)
# numeric polariton roots vs analytic branches
A = np.array([[K**2+m2, 1j*K*eta], [-1j*K*eta, K**2]])
w2 = np.sort(np.linalg.eigvalsh(A))
assert abs(w2[0]-W2_photon(K)) < 1e-12 and abs(w2[1]-((2*K**2+m2+np.sqrt(m2**2+4*eta**2*K**2))/2)) < 1e-12
print(f"2x2 eigencheck at K=1.4: photon W^2={w2[0]:.6f}, matter {w2[1]:.6f} : PASS")
# on-shell bare-massless closed form F8c
K=1.4; e=1.0; B2 = e/(2*K**2); dm2_direct = mu*eta**4*K**4*B2**2/m2**4
assert abs(dm2_direct - mu*eta**4*e**2/(4*m2**4)) < 1e-15
print(f"F8c on-shell closed form: dm2 = mu eta^4 e^2/(4 m2^4) = {mu*eta**4/(4*m2**4):+.5f} e^2 : PASS")
