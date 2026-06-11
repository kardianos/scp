#!/usr/bin/env python3
"""flavored_qball.py — multi-frequency ("flavored") 3-component Q-ball BVP.

The symmetric Q-ball has Phi_a = f(r) e^{i w t} for all three components
(equal charge partition Q_a = Q/3). This solver asks whether the theory also
contains ASYMMETRIC stationary baryons

    Phi_a = f_a(r) e^{i w_a t},   w_0 != w_1 = w_2

(s = prod f_a^2 stays time-independent, so the ansatz is consistent).
Charge partition Q_a = w_a int f_a^2 dV then differs per component — the
in-model analog of different flavor content (uud vs udd).

System (ungauged, g=0):
    f_a'' + (2/r) f_a' = (m^2 - w_a^2) f_a + 2 Vt'(s) f_a prod_{b!=a} f_b^2
Newton relaxation on a uniform grid, block-tridiagonal Jacobian via
scipy.sparse; continuation in w_0 from the symmetric solution.

Usage: flavored_qball.py [dw_max=0.03] [w_base=1.42]
"""
import sys
import numpy as np
from scipy.sparse import lil_matrix, csr_matrix
from scipy.sparse.linalg import spsolve

M2, MU, KAP = 2.25, -41.345, 50.0
H, R = 0.05, 40.0
N = int(R / H)
r = (np.arange(N) + 1) * H          # r=h..R (avoid r=0; ghost handles f'(0)=0)


def Vp(s):
    return 0.5 * MU / (1.0 + KAP * s) ** 2


def Vpp(s):
    return -KAP * MU / (1.0 + KAP * s) ** 3


def residual(F, w2):
    """F: (3,N). Returns (3,N) residual with f'(0)=0 and f(R)=0 BCs."""
    res = np.zeros_like(F)
    s = F[0] ** 2 * F[1] ** 2 * F[2] ** 2
    vp = Vp(s)
    for a in range(3):
        f = F[a]
        fp = np.empty(N); fm = np.empty(N)
        fp[:-1] = f[1:]; fp[-1] = 0.0            # f(R)=0
        fm[1:] = f[:-1]; fm[0] = f[1]            # ghost: f(-h)=f(h) -> f'(0)=0
        lap = (fp - 2 * f + fm) / H ** 2 + (2.0 / r) * (fp - fm) / (2 * H)
        prest = np.ones(N)
        for b in range(3):
            if b != a:
                prest *= F[b] ** 2
        res[a] = lap - (M2 - w2[a]) * f - 2 * vp * f * prest
    return res


def jacobian(F, w2):
    """Block-tridiagonal Jacobian as sparse matrix (3N x 3N)."""
    s = F[0] ** 2 * F[1] ** 2 * F[2] ** 2
    vp, vpp = Vp(s), Vpp(s)
    J = lil_matrix((3 * N, 3 * N))
    diag_lap = -2.0 / H ** 2
    up = 1.0 / H ** 2 + (2.0 / r) / (2 * H)
    dn = 1.0 / H ** 2 - (2.0 / r) / (2 * H)
    prest = np.empty((3, N))
    for a in range(3):
        p = np.ones(N)
        for b in range(3):
            if b != a:
                p *= F[b] ** 2
        prest[a] = p
    # ds/df_b = 2 f_b * prod_{c!=b} f_c^2 = 2 f_b prest[b]
    for a in range(3):
        base = a * N
        for i in range(N):
            # diagonal block: d res_a[i] / d f_a[i]
            d = diag_lap - (M2 - w2[a]) - 2 * vp[i] * prest[a, i] \
                - 2 * vpp[i] * (2 * F[a, i] * prest[a, i]) * F[a, i] * prest[a, i]
            J[base + i, base + i] = d
            if i + 1 < N:
                J[base + i, base + i + 1] = up[i]
            if i - 1 >= 0:
                J[base + i, base + i - 1] = dn[i]
            else:
                J[base + i, base + 1] += dn[i]   # ghost f(-h)=f(h)
            # off-blocks: d res_a / d f_b, b != a
            for b in range(3):
                if b == a:
                    continue
                pab = np.ones(1)
                # prod_{c != a, c != b} f_c^2
                c = 3 - a - b
                pab = F[c, i] ** 2
                dfb = -2 * vpp[i] * (2 * F[b, i] * prest[b, i]) * F[a, i] * prest[a, i] \
                      - 2 * vp[i] * F[a, i] * 2 * F[b, i] * pab
                J[base + i, b * N + i] = dfb
    return csr_matrix(J)


def solve(F, w2, tol=1e-10, itmax=40):
    for it in range(itmax):
        res = residual(F, w2)
        rn = np.max(np.abs(res))
        if rn < tol:
            return F, True, it, rn
        dF = spsolve(jacobian(F, w2), -res.reshape(-1)).reshape(3, N)
        lam = 1.0
        for _ in range(8):
            Ft = F + lam * dF
            if np.max(np.abs(residual(Ft, w2))) < rn:
                F = Ft
                break
            lam *= 0.5
        else:
            return F, False, it, rn
    return F, False, itmax, rn


def observables(F, w):
    Q = [4 * np.pi * w[a] * np.sum(F[a] ** 2 * r ** 2) * H for a in range(3)]
    s = F[0] ** 2 * F[1] ** 2 * F[2] ** 2
    grad = sum(np.sum(np.gradient(F[a], H) ** 2 * r ** 2) for a in range(3))
    E = 4 * np.pi * H * (
        sum(np.sum((w[a] ** 2 + M2) * F[a] ** 2 * r ** 2) for a in range(3)) * 0.5 * 2
        / 2  # (w^2 f^2 + m^2 f^2) kinetic+mass, real-field convention x3 comps
        + np.sum(0.5 * MU * s / (1 + KAP * s) * r ** 2)
        + 0.5 * grad)
    return Q, E


def main():
    dw_max = float(sys.argv[1]) if len(sys.argv) > 1 else 0.03
    wb = float(sys.argv[2]) if len(sys.argv) > 2 else 1.42

    # symmetric start from the v66 profile resampled
    src = np.loadtxt("/home/d/code/scp/v69/theory/gprofile_w142_g005.txt")
    f0 = np.interp(r, src[:, 0], src[:, 1], right=0.0)
    F = np.vstack([f0, f0, f0])
    w = np.array([wb, wb, wb])
    F, ok, it, rn = solve(F, w ** 2)
    print(f"symmetric re-solve at w={wb}: ok={ok} iters={it} resid={rn:.2e}")
    if not ok:
        sys.exit(1)
    Q, E = observables(F, w)
    print(f"  Q_a = {Q[0]:.2f} {Q[1]:.2f} {Q[2]:.2f}  (sum {sum(Q):.2f})")

    print(f"\ncontinuation: w_0 = wb - dw, w_1 = w_2 = wb (asymmetric flavor)")
    print(f"{'dw':>7s} {'w0':>8s} {'f0(0)':>8s} {'f1(0)':>8s} {'Q_0':>9s} "
          f"{'Q_1':>9s} {'Q_0/Q_1':>8s} {'E':>10s}")
    for dw in np.arange(0.005, dw_max + 1e-9, 0.005):
        w = np.array([wb - dw, wb, wb])
        F, ok, it, rn = solve(F, w ** 2)
        if not ok:
            print(f"{dw:7.3f}  NO CONVERGENCE (resid {rn:.1e}) — branch end")
            break
        Q, E = observables(F, w)
        print(f"{dw:7.3f} {w[0]:8.4f} {F[0,0]:8.4f} {F[1,0]:8.4f} {Q[0]:9.2f} "
              f"{Q[1]:9.2f} {Q[0]/Q[1]:8.4f} {E:10.2f}")
    np.savetxt("/home/d/code/scp/v71/results/flavored_profile_last.txt",
               np.column_stack([r] + [F[a] for a in range(3)]),
               header=f"r f0 f1 f2   (w = {w[0]:.4f} {w[1]:.4f} {w[2]:.4f})")
    print("\nwrote v71/results/flavored_profile_last.txt")


if __name__ == "__main__":
    main()
