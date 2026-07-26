#!/usr/bin/env python3
"""ws5 — FABLE seat, v87 B0.

Monte Carlo of the Door-A model with an EXPLICIT pseudorandom event source.

Event source (the operational story, per O1): each run, the fabric is in a
definite configuration; the uncontrolled degrees of freedom are summarised by
lambda in S^2. The experimenter does not control lambda between runs; the
sampling below IS the epistemic distribution rho_p(lambda|b) — the
measurement-dependent piece is that the source distribution is tilted toward
the B-side analyser axis with weight |lambda.b|^p.

Sampling: with b = z-axis, the density factorises: azimuth v uniform on
[0,2pi); z = lambda.b has density proportional to |z|^p on [-1,1], so
|z| = U^(1/(p+1)), sign(z) = +-1 with prob 1/2. Then rotate z-axis -> b.
Responses: A = sgn(a.lambda), B = -sgn(b.lambda). Deterministic given lambda.

PRNG: numpy PCG64, seed printed. N = 2e7 per correlator.
Expected: S(p=1) = 2*sqrt(2) +- ~4/sqrt(4N); S(p=0) = 2; S(p=8) > 2*sqrt(2).
"""
import numpy as np

SEED = 20260726
N = 20_000_000

def rot_from_z(b):
    b = b / np.linalg.norm(b)
    if abs(b[2]) > 0.999999:
        return np.eye(3) * (1 if b[2] > 0 else -1)
    v = np.cross([0, 0, 1], b); s = np.linalg.norm(v); c = b[2]
    vx = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
    return np.eye(3) + vx + vx @ vx * ((1 - c) / s**2)

def sample_lambda(rng, b, pval, n):
    z = rng.random(n) ** (1.0 / (pval + 1.0))
    z *= np.where(rng.random(n) < 0.5, -1.0, 1.0)
    v = rng.random(n) * 2 * np.pi
    r = np.sqrt(1 - z**2)
    lam = np.stack([r*np.cos(v), r*np.sin(v), z], axis=1)
    return lam @ rot_from_z(b).T

def E_mc(rng, a, b, pval, n):
    lam = sample_lambda(rng, b, pval, n)
    A = np.sign(lam @ a)
    B = -np.sign(lam @ b)
    return float(np.mean(A * B))

def unit(theta):
    return np.array([np.sin(theta), 0.0, np.cos(theta)])

def run(pval, angles, label):
    rng = np.random.default_rng(SEED + int(100*pval))
    al, alp, be, bep = angles
    a, ap, b, bp = unit(al), unit(alp), unit(be), unit(bep)
    E = {}
    E['ab']  = E_mc(rng, a,  b,  pval, N)
    E['abp'] = E_mc(rng, a,  bp, pval, N)
    E['apb'] = E_mc(rng, ap, b,  pval, N)
    E['apbp']= E_mc(rng, ap, bp, pval, N)
    S = E['ab'] + E['abp'] + E['apb'] - E['apbp']
    sig = 2.0 / np.sqrt(N)   # each E has sd <= 1/sqrt(N); S combines 4
    print(f"{label}: E = {E['ab']:+.5f} {E['abp']:+.5f} {E['apb']:+.5f} "
          f"{E['apbp']:+.5f}   S = {S:.6f} +- {sig:.6f}")
    return S, sig

print(f"PRNG: numpy PCG64, base seed {SEED}, N = {N:,} per correlator\n")

# p = 1 at the quantum-optimal angles (ws1): a=0, a'=pi/2, b=pi+pi/4, b'=pi-pi/4
ang_q = (0.0, np.pi/2, np.pi + np.pi/4, np.pi - np.pi/4)
S1, s1 = run(1.0, ang_q, "p=1, quantum angles")
tsir = 2*np.sqrt(2)
print(f"   vs 2sqrt2 = {tsir:.6f}: deviation {abs(S1-tsir):.6f} "
      f"({abs(S1-tsir)/s1:.2f} sigma)")
assert abs(S1 - tsir) < 5*s1

# p = 0 (uniform lambda — the MI-respecting local model) at ITS optimum:
# sawtooth E = -1+2theta/pi; optimum e.g. all four |angles| -> boundary of
# monotone branch; the sawtooth model attains exactly S=2 at generic non-
# degenerate angles: use theta_ab=theta_ab'=theta_a'b=pi (E=-1... ) — simplest
# verified optimum: a=0, a'=0, b=pi, b'=pi gives S = -(-1)-... use LP-style
# angles: a=0,a'=pi, b=pi, b'=0: E(0,pi)=1... just scan numerically instead.
best0 = (-9, None)
rngscan = np.random.default_rng(3)
for _ in range(2000):
    x = rngscan.uniform(0, 2*np.pi, 4)
    th = lambda i, j: np.abs((x[i]-x[j]+np.pi) % (2*np.pi) - np.pi)
    Es = [-1+2*th(0,2)/np.pi, -1+2*th(0,3)/np.pi, -1+2*th(1,2)/np.pi,
          -1+2*th(1,3)/np.pi]
    Sv = Es[0]+Es[1]+Es[2]-Es[3]
    if Sv > best0[0]:
        best0 = (Sv, x)
print(f"\np=0 sawtooth analytic scan optimum: S = {best0[0]:.6f} (must be <= 2)")
S0, s0 = run(0.0, tuple(best0[1]), "p=0, its optimal angles")
assert S0 < 2 + 5*s0

# p = 8: overshoot demonstration (the family has no internal Tsirelson wall)
# angles for the near-step correlator: a=0, a'=100deg, b=50deg, b'=-50deg
ang_s = (0.0, np.radians(100), np.radians(50), np.radians(-50))
S8, s8 = run(8.0, ang_s, "p=8, step-adapted angles")
print(f"   p=8 overshoots Tsirelson: |S| = {abs(S8):.4f} > 2sqrt2 = {tsir:.4f}")
assert abs(S8) > tsir   # sign is a labelling convention; |S| is the object

print("\nws5 OK: MC matches closed form; p=1 hits Tsirelson, p=0 local, "
      "p=8 overshoots.")
