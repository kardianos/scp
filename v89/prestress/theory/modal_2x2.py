#!/usr/bin/env python3
"""
modal_2x2.py -- v89 prestress: closed-form rate-level linearization of a
consonant N-ring about the locked uniform state, per Fourier mode n.
Companion (and cross-check) to the exact Floquet computation in
ring_map.py; the derivation chain is written out in ../STABILITY.md S2-S3.

Per-voice state (dth_k, dE_k), mode ansatz ~ exp(i q k), q = 2 pi n / N,
z = e^{iq}. Uniform locked ring: forward gates at the flat top (psi_f=0),
back gates at psi_b = wrap(-2 phi), phi = w d / C the per-link retard.

Rate-level channels kept (all first order):
  gamma   = w2 q_d /(1+q_d x)^2 / cap : pitch stiffness per unit store
  c_h     = 1/(h cap)                 : headroom log-slope (receiver)
  c_m     = 1/(2 E)                   : mobility log-slope (symmetric)
  c_tau   = tau*gamma                 : retard-angle shift per source store
  Gb'     = -p_gate*tan(psi_b/2)      : back-gate log-slope (0 on pi-rung)
  g_e     = klock*F*(1/(E'+F dt)+1/(E'+F tau)) : entrainment gain / err
Flight delay and slot dynamics are NOT in this matrix (the exact map has
them); expect agreement for |lambda| tau << 1 and departures for the fast
branch.
"""
import numpy as np
import math
import cmath

# law table (dense sector)
C = 1.0; W2 = 2.9; QD = 1.2; GM = 0.10; PG = 8
LOCKF = 0.005 * 2.5; KD = 2.4; CAP = 2.5
ES0 = 1.0; ESFL = 0.05; SPULL = 0.15; KLOCK = 1.0
R0 = 0.85; DT = 0.02


def lens_geo(d, r):
    if d >= 2 * r:
        return 0.0
    a2 = r * r - 0.25 * d * d
    A = math.pi * a2
    return (A / (math.pi * R0 * R0)) * (2 * R0 / d)


def fixed_point(d, m_over_N_times_2pi, gpl, nlinks_per_cell=2):
    """store/flight partition of the locked ring: x*cap = store+flload."""
    phi = m_over_N_times_2pi              # w*d/C per link
    w = phi * C / d
    xbar = (W2 / w - 1.0) / QD
    add = xbar * CAP / (1 + SPULL)
    Es = ES0 - min(SPULL * add, ES0 - ESFL)
    r = R0 * (Es / ES0) ** (1 / 3)
    geo = lens_geo(d, r)
    gb = (0.5 * (1 + math.cos(2 * math.pi - 2 * phi))) ** PG \
        if abs(math.sin(phi)) > 1e-12 else 1.0
    gb = (0.5 * (1 + math.cos(-2 * phi))) ** PG
    tau = d / C
    E = xbar * CAP * 0.7
    for _ in range(400):
        head = max(1 - E / CAP, 0.0)
        base = KD * geo * gpl
        Ff = base * 1.0 * head * E
        Fb = base * gb * head * E
        fll = 0.5 * nlinks_per_cell * (Ff + Fb) * tau
        E = 0.5 * E + 0.5 * (xbar * CAP - fll)
        if E < 0.01:
            E = 0.01
    return dict(w=w, xbar=xbar, Es=Es, geo=geo, gb=gb, tau=tau, E=E,
                head=1 - E / CAP, Ff=Ff, Fb=Fb, fll=fll, gpl=gpl,
                psi_b=(math.pi + (-2 * phi)) % (2 * math.pi) - math.pi)


def mode_matrix(fp, q, use_gbslope=True):
    z = cmath.exp(1j * q)
    E, hd, tau = fp['E'], fp['head'], fp['tau']
    Ff, Fb = fp['Ff'], fp['Fb']
    w = fp['w']
    gamma = W2 * QD / (1 + QD * (fp['xbar'])) ** 2 / CAP   # dw per dE
    c_h = 1 / (hd * CAP)
    c_m = 1 / (2 * E)
    c_t = tau * gamma
    psb = fp['psi_b']
    Gbp = -PG * math.tan(psb / 2) if use_gbslope else 0.0
    Ep = E + LOCKF
    ge_per_F = KLOCK * (1 / (Ep) + 1 / (Ep + (Ff + Fb) * tau))
    gef = ge_per_F * Ff
    geb = ge_per_F * Fb
    # d psi_f_k = (1-z) Th + c_t E ; d psi_b_k = (z-1) Th + c_t z E
    # dF^f_k / Ff = -z c_h E + (1+z) c_m E
    # dF^b_k / Fb = -c_h E + (1+z) c_m E + Gbp ((z-1) Th + c_t z E)
    M = np.zeros((2, 2), dtype=complex)
    # phase row: dth' = -gamma E + gef*dpsi_f(k-1) + geb*dpsi_b(k)
    #            + klock*psi_b*d(rate of back arrivals)  [kept via geb Gbp]
    M[0, 0] = gef * (1 / z) * (1 - z) + geb * (z - 1) \
        + geb * psb * Gbp * (z - 1)
    M[0, 1] = -gamma + gef * (1 / z) * c_t + geb * c_t * z \
        + geb * psb * (Gbp * c_t * z - c_h + (1 + z) * c_m)
    # energy row: dE' = [dF^f_{k-1} - dF^f_k] + [dF^b_k - dF^b_{k-1}]
    M[1, 0] = Fb * Gbp * (z - 1) * (1 - 1 / z)
    M[1, 1] = Ff * (1 / z - 1) * (-z * c_h + (1 + z) * c_m) \
        + Fb * (1 - 1 / z) * (-c_h + (1 + z) * c_m + Gbp * c_t * z)
    return M


def pair_rate(d=1.25, gpl=1.0):
    fp = fixed_point(d, math.pi, gpl, nlinks_per_cell=1)
    lam0 = 2 * KD * fp['geo'] * gpl * fp['E'] / CAP
    # flight-delay correction: lam = -lam0 * exp(-lam*tau)
    lam = -lam0
    for _ in range(60):
        lam = -lam0 * math.exp(-lam * fp['tau'])
    return fp, lam0, lam


if __name__ == '__main__':
    print('== pair (d=1.25, m=1, gpl=1): the sign theorem rate ==')
    fp, lam0, lam = pair_rate()
    print(f'   fixed point: store E={fp["E"]:.4f} flload={fp["fll"]:.4f} '
          f'head={fp["head"]:.4f} F={fp["Ff"]:.4f}')
    print(f'   lambda_h (no delay) = {-lam0:+.4f} /t.u.')
    print(f'   lambda_h (with flight delay e^{{-lam tau}}) = {lam:+.4f} /t.u.')
    print(f'   exact map (Floquet): -0.1914 /t.u.')

    print('== ring12 (m=5, phi=150deg, gpl=0.5625): per-mode 2x2 ==')
    fp = fixed_point(1.25, 2 * math.pi * 5 / 12, 0.5625)
    print(f'   fixed point: E={fp["E"]:.4f} fll={fp["fll"]:.4f} '
          f'Ff={fp["Ff"]:.4f} Fb={fp["Fb"]:.4f} psi_b={fp["psi_b"]:+.4f} '
          f"Gb'={-PG*math.tan(fp['psi_b']/2):+.3f}")
    print('   n     with Gb-slope           without (frozen g_b)')
    for n in range(0, 7):
        q = 2 * math.pi * n / 12
        ev1 = np.linalg.eigvals(mode_matrix(fp, q, True))
        ev0 = np.linalg.eigvals(mode_matrix(fp, q, False))
        s1 = ' '.join(f'{e.real:+.4f}{e.imag:+.3f}i' for e in
                      sorted(ev1, key=lambda e: -e.real))
        s0 = ' '.join(f'{e.real:+.4f}{e.imag:+.3f}i' for e in
                      sorted(ev0, key=lambda e: -e.real))
        print(f'   {n}  {s1}   |   {s0}')

    print('== ring6 (m=3, phi=pi, gpl=0.0625): slow branch vs Floquet ==')
    fp = fixed_point(1.25, math.pi, 0.0625)
    print(f'   fixed point: E={fp["E"]:.4f} fll={fp["fll"]:.4f} '
          f'Ff={fp["Ff"]:.4f} (g_b=1: mutual)')
    print('   n   2x2 eigenvalues                Floquet (exact map)')
    flo = {0: '0 (cons.) / 0 (gauge)', 1: '-0.0077 / -0.146',
           2: '-0.0302 / -0.462', 3: '-0.0416 / -0.635'}
    for n in range(0, 4):
        q = 2 * math.pi * n / 6
        ev = np.linalg.eigvals(mode_matrix(fp, q, True))
        s = ' '.join(f'{e.real:+.4f}{e.imag:+.3f}i' for e in
                     sorted(ev, key=lambda e: -e.real))
        print(f'   {n}  {s}    {flo[n]}')


# ---- delayed-arrival characteristic equation (the instability's home) ----
def delayed_roots(fp, q, use_gbslope=True):
    """dE_k/dt = arrivals(t-tau) - outbound(t); phase row instant.
    Solve det(lam I - M(lam)) = 0 by Newton from undelayed seeds."""
    z = cmath.exp(1j * q)
    E, hd, tau = fp['E'], fp['head'], fp['tau']
    Ff, Fb = fp['Ff'], fp['Fb']
    gamma = W2 * QD / (1 + QD * fp['xbar']) ** 2 / CAP
    c_h = 1 / (hd * CAP)
    c_m = 1 / (2 * E)
    c_t = tau * gamma
    psb = fp['psi_b']
    Gbp = -PG * math.tan(psb / 2) if use_gbslope else 0.0
    Ep = E + LOCKF
    ge_per_F = KLOCK * (1 / Ep + 1 / (Ep + (Ff + Fb) * tau))
    gef, geb = ge_per_F * Ff, ge_per_F * Fb
    # phase row (instant)
    Ptt = gef * (1 / z) * (1 - z) + geb * (z - 1) + geb * psb * Gbp * (z - 1)
    PtE = -gamma + gef * (1 / z) * c_t + geb * c_t * z \
        + geb * psb * (Gbp * c_t * z - c_h + (1 + z) * c_m)
    # energy row split
    Mout_t = Fb * Gbp * (1 - 1 / z)
    Mout_E = Ff * (-z * c_h + (1 + z) * c_m) \
        + Fb * (-c_h / z + (1 / z + 1) * c_m + Gbp * c_t)
    Min_t = Fb * Gbp * (z - 1)
    Min_E = Ff * (-c_h + (1 / z + 1) * c_m) \
        + Fb * (-c_h + (1 + z) * c_m + Gbp * c_t * z)

    def f(lam):
        e = cmath.exp(-lam * tau)
        Ett = Min_t * e - Mout_t
        EEE = Min_E * e - Mout_E
        return (lam - Ptt) * (lam - EEE) - PtE * Ett

    roots = []
    M0 = np.array([[Ptt, PtE], [Min_t - Mout_t, Min_E - Mout_E]])
    for lam in np.linalg.eigvals(M0):
        for _ in range(80):
            h = 1e-8
            df = (f(lam + h) - f(lam - h)) / (2 * h)
            if abs(df) < 1e-300:
                break
            step = f(lam) / df
            lam = lam - step
            if abs(step) < 1e-13:
                break
        roots.append(lam)
    return roots


def _delayed_table():
    print('== ring12 delayed-arrival roots (the chiral pump) ==')
    fp = fixed_point(1.25, 2 * math.pi * 5 / 12, 0.5625)
    print('   n    with Gb-slope (delayed)      without')
    for n in range(0, 7):
        q = 2 * math.pi * n / 12
        r1 = sorted(delayed_roots(fp, q, True), key=lambda e: -e.real)
        r0 = sorted(delayed_roots(fp, q, False), key=lambda e: -e.real)
        s1 = ' '.join(f'{e.real:+.4f}{e.imag:+.3f}i' for e in r1)
        s0 = ' '.join(f'{e.real:+.4f}{e.imag:+.3f}i' for e in r0)
        print(f'   {n}  {s1}   |   {s0}')
    print('   (exact map Floquet, dt=0.02: n=1 +0.035-0.640i, '
          'n=2 +0.053-1.054i, n=3 +0.002-1.203i)')


if __name__ == '__main__' and True:
    _delayed_table()
