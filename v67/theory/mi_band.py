#!/usr/bin/env python3
"""v67 modulational-instability band evaluation + full 12-field (48x48) verification.

Phi-only analytic dispersion (verified symbolically in mi_dispersion.mac):
  per eigenmode M of the R-sector matrix C_ab = 4A^4[A^6 Vt'' + Vt'(1-delta_ab)]:
    sigma^4 + sigma^2 (2k^2 + M + 4 w_c^2) + k^2 (k^2 + M) = 0
  growth iff M<0, band 0 < k^2 < -M.
  M_sym  = 4 A^4 mu (1 - 2 kappa A^6) / (1 + kappa A^6)^3   (mode (1,1,1))
  M_anti = -2 mu A^4 / (1 + kappa A^6)^2                    (x2, always >0 for mu<0)

Full verification: linearize ALL 12 fields (u,v,tu,tv) about the condensate in the
corotating frame (U(1) rotation applied to both (u,v) and (tu,tv) doublets makes the
system autonomous; eta-curl coupling commutes with the rotation). Real cos/sin spatial
components -> 24 second-order = 48x48 first-order matrix; growth = max Re eig.
"""
import numpy as np

m2, mu, kappa, eta = 2.25, -41.345, 50.0, 0.5

Vt   = lambda s: 0.5*mu*s/(1.0+kappa*s)
Vtp  = lambda s: 0.5*mu/(1.0+kappa*s)**2
Vtpp = lambda s: -mu*kappa/(1.0+kappa*s)**3

def background(A):
    s = A**6
    wc2 = m2 + 2.0*Vtp(s)*A**4
    return s, wc2, np.sqrt(wc2)

def Msym(A):
    s = A**6
    return 4.0*A**4*mu*(1.0-2.0*kappa*s)/(1.0+kappa*s)**3

def Manti(A):
    return -2.0*mu*A**4/(1.0+kappa*A**6)**2

def sigma2(k2, M, w2):
    B = 2.0*k2 + M + 4.0*w2
    C = k2*(k2+M)
    disc = B*B - 4.0*C
    if disc < 0: return -1.0
    return 0.5*(-B + np.sqrt(disc))

def band_max(A):
    """returns (k_edge, k_max, sigma_max) for the symmetric mode, or None if stable"""
    M = Msym(A); _, w2, _ = background(A)
    if M >= 0: return None
    k2e = -M
    k2g = np.linspace(1e-9, k2e*(1-1e-9), 20001)
    sg2 = np.array([sigma2(k2, M, w2) for k2 in k2g])
    i = int(np.argmax(sg2))
    # golden refine
    lo, hi = k2g[max(i-1,0)], k2g[min(i+1,len(k2g)-1)]
    gr = (np.sqrt(5)-1)/2
    a, b = lo, hi
    c, d = b-gr*(b-a), a+gr*(b-a)
    for _ in range(80):
        if sigma2(c,M,w2) > sigma2(d,M,w2): b, d = d, c; c = b-gr*(b-a)
        else: a, c = c, d; d = a+gr*(b-a)
    k2m = 0.5*(a+b)
    return np.sqrt(k2e), np.sqrt(k2m), np.sqrt(sigma2(k2m, M, w2))

# ---------------- full 48x48 (12 fields, cos/sin spatial pair) ----------------
def cross_mat(khat):
    kx,ky,kz = khat
    return np.array([[0,-kz,ky],[kz,0,-kx],[-ky,kx,0]], float)

def full_matrix(A, k, khat, eta_):
    s, w2, w = background(A)
    p = 4.0*A**10*Vtpp(s)               # C diagonal
    q = 4.0*A**4*Vtp(s)                 # C off-diagonal extra
    C3 = p*np.ones((3,3)) + q*(np.ones((3,3))-np.eye(3))
    I3, Z3 = np.eye(3), np.zeros((3,3))
    K = cross_mat(khat)
    # field blocks (each 6 = [c(3), s(3)]): R, I, RT, IT  -> 24 positions
    # curl on (c,s) pair: (curl v)_c = k K v_s ; (curl v)_s = -k K v_c
    def curl6():
        Cm = np.zeros((6,6)); Cm[0:3,3:6] = k*K; Cm[3:6,0:3] = -k*K; return Cm
    C6 = np.kron(np.eye(2), C3)         # potential acts same on c and s
    L6 = curl6()
    k2 = k*k
    # accelerations: X'' = S X + D X'
    S = np.zeros((24,24)); D = np.zeros((24,24))
    iR, iI, iRT, iIT = slice(0,6), slice(6,12), slice(12,18), slice(18,24)
    # R'' = 2w I' - (k2 + C) R + eta curl(RT)
    S[iR,iR] = -(k2*np.eye(6) + C6);  S[iR,iRT] =  eta_*L6;  D[iR,iI]  =  2*w*np.eye(6)
    # I'' = -2w R' - k2 I + eta curl(IT)
    S[iI,iI] = -k2*np.eye(6);         S[iI,iIT] =  eta_*L6;  D[iI,iR]  = -2*w*np.eye(6)
    # RT'' = 2w IT' + (w2 - k2) RT + eta curl(R)
    S[iRT,iRT] = (w2-k2)*np.eye(6);   S[iRT,iR] =  eta_*L6;  D[iRT,iIT] =  2*w*np.eye(6)
    # IT'' = -2w RT' + (w2 - k2) IT + eta curl(I)
    S[iIT,iIT] = (w2-k2)*np.eye(6);   S[iIT,iI] =  eta_*L6;  D[iIT,iRT] = -2*w*np.eye(6)
    Mfull = np.zeros((48,48))
    Mfull[0:24,24:48] = np.eye(24)
    Mfull[24:48,0:24] = S
    Mfull[24:48,24:48] = D
    return Mfull

def full_growth(A, k, khat, eta_):
    return float(np.max(np.linalg.eigvals(full_matrix(A, k, khat, eta_)).real))

KHATS = {"z":(0,0,1.0), "111":(1/np.sqrt(3),)*3, "1m10":(1/np.sqrt(2),-1/np.sqrt(2),0)}

def full_scan(A, eta_, kmaxscan=2.0, nk=241):
    out = {}
    for name,kh in KHATS.items():
        ks = np.linspace(1e-4, kmaxscan, nk)
        gs = np.array([full_growth(A,k,kh,eta_) for k in ks])
        i = int(np.argmax(gs))
        # refine
        lo,hi = ks[max(i-1,0)], ks[min(i+1,nk-1)]
        kk = np.linspace(lo,hi,41)
        gg = np.array([full_growth(A,k,kh,eta_) for k in kk])
        j = int(np.argmax(gg))
        out[name] = (kk[j], gg[j])
    return out

if __name__ == "__main__":
    Astar = (2.0*kappa)**(-1.0/6.0)
    print(f"threshold A* = (2 kappa)^(-1/6) = {Astar:.6f}   (kappa A*^6 = 1/2)")
    print(f"{'A':>6} {'kap*A^6':>8} {'w_c':>7} {'rho_Q':>7} {'e_0':>7} {'M_sym':>8} {'M_anti':>7}")
    for A in (0.30,0.35,0.40,0.4642,0.50,0.585):
        s,w2,w = background(A)
        e0 = 1.5*(w2+m2)*A*A + Vt(s)
        print(f"{A:6.3f} {kappa*s:8.4f} {w:7.4f} {3*w*A*A:7.4f} {e0:7.4f} {Msym(A):8.4f} {Manti(A):7.4f}")
    print()
    L = 25.0; side = 2*L   # N=192, L=25 box, side 50
    print("Unstable cases (phi-only analytic):")
    print(f"{'A':>6} {'k_edge':>7} {'k_max':>7} {'sig_max':>8} {'lambda':>7} {'(2L/lam)^3':>10} {'1/sig':>6} {'t_frag~3/sig':>12} {'Q_box':>9} {'Q/frag':>7}")
    for A in (0.30,0.35,0.40,0.44,0.45,0.50,0.585):
        r = band_max(A)
        s,w2,w = background(A)
        if r is None:
            print(f"{A:6.3f}   linearly STABLE (M_sym = {Msym(A):+.4f} > 0)")
            continue
        ke,km,sm = r
        lam = 2*np.pi/km; nfrag = (side/lam)**3
        Qbox = 3*w*A*A*side**3
        print(f"{A:6.3f} {ke:7.4f} {km:7.4f} {sm:8.4f} {lam:7.3f} {nfrag:10.1f} {1/sm:6.2f} {3/sm:12.1f} {Qbox:9.0f} {Qbox/nfrag:7.0f}")
    print()
    print("Full 48x48 (all 12 fields, corotating frame), max growth over k per direction:")
    print(f"{'A':>6} {'eta':>5} | " + " | ".join(f"{n}: (k*, sigma*)" for n in KHATS))
    for A in (0.40,0.50,0.585):
        for e_ in (0.0, 0.5):
            sc = full_scan(A, e_)
            row = " | ".join(f"{n}: ({sc[n][0]:.4f}, {sc[n][1]:.5f})" for n in KHATS)
            print(f"{A:6.3f} {e_:5.2f} | {row}")
    # cross-check analytic vs full at eta=0, A=0.40
    r = band_max(0.40); ke,km,sm = r
    g = full_growth(0.40, km, KHATS["z"], 0.0)
    print(f"\ncross-check A=0.40 eta=0: analytic sigma(k_max)={sm:.6f}, 48x48={g:.6f}, diff={abs(sm-g):.2e}")
