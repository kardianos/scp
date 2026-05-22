#!/usr/bin/env python3
"""
02_silent_direction.py

Test: for ξ ∈ S³ ⊂ ℍ, do there exist directions of motion in S³ that leave the
Brannen eigenvalues invariant?

Motivation: the simple Cosserat coupling on the Brannen phase δ (from
01_cosserat_test.py) gives ∂(ln m_e)/∂δ ≈ -51, ruled out by atomic-clock
bounds at 8 orders of magnitude.  The escape is for the strain to push ξ in a
"silent" direction on S³ that doesn't enter the spectrum.

The v59 paper notes: "20 random quaternionic ξ on the constraint S³ gave
Koide values in [0.6666637, 0.6666694]" — Q is essentially invariant on S³.
But individual eigenvalues might vary, OR there might be a 1-parameter family
of silent rotations.

This script:
  1. Builds the 12×12 real Brannen mass-operator over ℍ (from 04_findings.md).
  2. Tabulates the 3 unique eigenvalues for many ξ on S³.
  3. Identifies the ANGULAR DOF on S³ that does NOT shift individual eigenvalues.
  4. Computes the response of each eigenvalue to motion in each independent
     direction on S³.  The "silent direction" is the one with zero response.
"""

import numpy as np

# ----------------------------------------------------------------------
# Quaternion arithmetic (from 04_quaternionic_constraint.py)
# ----------------------------------------------------------------------
def qmul(p, q):
    a1, b1, c1, d1 = p
    a2, b2, c2, d2 = q
    return np.array([
        a1*a2 - b1*b2 - c1*c2 - d1*d2,
        a1*b2 + b1*a2 + c1*d2 - d1*c2,
        a1*c2 - b1*d2 + c1*a2 + d1*b2,
        a1*d2 + b1*c2 - c1*b2 + d1*a2,
    ])

def qconj(p):
    a, b, c, d = p
    return np.array([a, -b, -c, -d])

def L_xi(xi):
    """4×4 matrix of left-multiplication by quaternion ξ."""
    a, b, c, d = xi
    return np.array([
        [a, -b, -c, -d],
        [b,  a, -d,  c],
        [c,  d,  a, -b],
        [d, -c,  b,  a],
    ])

def build_M(a_scalar, xi):
    """12×12 real matrix of M(ξ) = a_scalar·I + ξ·S + ξ̄·S^T."""
    I4 = np.eye(4)
    Lx = L_xi(xi)
    Lxs = L_xi(qconj(xi))
    M = np.zeros((12, 12))
    for k in range(3):
        M[4*k:4*k+4, 4*k:4*k+4] = a_scalar * I4
    M[0:4, 4:8] = Lx
    M[4:8, 0:4] = Lxs
    M[4:8, 8:12] = Lx
    M[8:12, 4:8] = Lxs
    M[8:12, 0:4] = Lx
    M[0:4, 8:12] = Lxs
    return M

def unique_eigs(M, tol=1e-7):
    """Return the unique eigenvalues of M (grouped at precision `tol`)."""
    raw = np.sort(np.linalg.eigvalsh(M))
    out = []
    for x in raw:
        if not out or abs(x - out[-1]) > tol:
            out.append(x)
    return np.array(out)

# ----------------------------------------------------------------------
# Step 1.  Parameterize S³ ⊂ ℍ in Hopf coordinates.
# ----------------------------------------------------------------------
# Standard parameterization:
#   ξ = (1/√2) · (cos(η/2) e^{iα},  sin(η/2) e^{iβ})    in ℂ²
#     = (1/√2) · (cos(η/2)cos α + i cos(η/2)sin α + j sin(η/2)cos β + k sin(η/2)sin β)
# Components: a = cos(η/2)cos α / √2
#             b = cos(η/2)sin α / √2
#             c = sin(η/2)cos β / √2
#             d = sin(η/2)sin β / √2
# (η ∈ [0, π], α ∈ [0, 2π], β ∈ [0, 2π])
#
# The Hopf fiber is parameterized by ψ = (α + β)/2 (or α − β) — varying ψ
# while keeping η and the orthogonal phase fixed should be a U(1) action.

def hopf_to_quat(eta, alpha, beta):
    """Hopf-coordinate quaternion on S³ of radius 1/√2."""
    r = 1.0 / np.sqrt(2)
    return r * np.array([
        np.cos(eta/2) * np.cos(alpha),
        np.cos(eta/2) * np.sin(alpha),
        np.sin(eta/2) * np.cos(beta),
        np.sin(eta/2) * np.sin(beta),
    ])

# ----------------------------------------------------------------------
# Step 2.  Identify the "physical" Brannen phase direction.
# ----------------------------------------------------------------------
# For ξ ∈ ℂ (pure complex), η = 0, β = arbitrary (but c=d=0 means sin(η/2)=0,
# so β is irrelevant). The phase is α.  This is the standard Brannen phase
# δ = 2/9 corresponds to α = 2/9, η = 0.

delta = 2/9  # the Brannen phase (within m_τ precision)

# Reference point: ξ = (cos δ, sin δ, 0, 0) / √2
xi_ref = hopf_to_quat(eta=0, alpha=delta, beta=0)
M_ref = build_M(1.0, xi_ref)
eig_ref = unique_eigs(M_ref)
print("=" * 70)
print("Reference point: ξ in pure-ℂ subspace at Brannen phase δ = 2/9")
print("=" * 70)
print(f"  ξ = {xi_ref}")
print(f"  |ξ|² = {np.dot(xi_ref, xi_ref):.6f}  (should be 0.5)")
print(f"  3 unique eigenvalues: {eig_ref}")
Q_ref = (eig_ref**2).sum() / eig_ref.sum()**2
print(f"  Koide Q = {Q_ref:.10f}  (expect 2/3 = {2/3:.10f})")

# ----------------------------------------------------------------------
# Step 3.  Scan each Hopf coordinate axis independently.
# ----------------------------------------------------------------------
print()
print("=" * 70)
print("Step 3: Sensitivity of each unique eigenvalue to motion on S³")
print("=" * 70)
print("""
For each Hopf coordinate (η, α, β), apply a small displacement ε from the
reference point.  Compute (eig - eig_ref)/ε.  A direction with zero response
for all 3 eigenvalues is a "silent" direction (Cosserat-strain candidate).
""")

eps = 1e-5
displacements = {
    'η (polar)': (eps, 0.0, 0.0),
    'α (1st phase)': (0.0, eps, 0.0),
    'β (2nd phase)': (0.0, 0.0, eps),
    'α + β (Hopf fiber)': (0.0, eps/2, eps/2),
    'α − β (orthogonal)': (0.0, eps/2, -eps/2),
}

results = {}
for name, (deta, dalpha, dbeta) in displacements.items():
    xi_pert = hopf_to_quat(eta=0+deta, alpha=delta+dalpha, beta=0+dbeta)
    M_pert = build_M(1.0, xi_pert)
    eig_pert = unique_eigs(M_pert)
    if len(eig_pert) != 3:
        print(f"  WARNING: {name} gave {len(eig_pert)} unique eigenvalues, skipping.")
        continue
    dE = (eig_pert - eig_ref) / eps
    Q_pert = (eig_pert**2).sum() / eig_pert.sum()**2
    dQ = (Q_pert - Q_ref) / eps
    results[name] = (dE, dQ)
    print(f"  {name:<28} | dE/dε for 3 eigenvalues: "
          f"({dE[0]:+8.4f}, {dE[1]:+8.4f}, {dE[2]:+8.4f}) | dQ/dε = {dQ:+8.4f}")

# Identify silent directions:
print()
print("-" * 70)
print("Silent-direction summary")
print("-" * 70)
max_response = {}
for name, (dE, dQ) in results.items():
    max_response[name] = np.max(np.abs(dE))
    sym = "✓ SILENT" if max_response[name] < 1e-3 else "  active"
    print(f"  {name:<28} | max |dE/dε| = {max_response[name]:.3e}  {sym}")

# ----------------------------------------------------------------------
# Step 4.  Larger sweep: do silent directions remain silent at finite displacement?
# ----------------------------------------------------------------------
print()
print("=" * 70)
print("Step 4: Finite-displacement test of candidate silent directions")
print("=" * 70)

# Sweep the candidate silent direction over a wider range
silent_candidates = [name for name in results if max_response[name] < 1e-3]
if not silent_candidates:
    print("  No silent direction found at infinitesimal level.")
else:
    print("  Sweeping each silent candidate over [0, 2π] and tracking the spectrum.")
    print("  If the eigenvalues remain INVARIANT, the silent direction is exact.")
    print()
    angles = np.linspace(0, 2*np.pi, 13)
    for name in silent_candidates:
        deta, dalpha, dbeta = displacements[name]
        # Direction unit vector in (η, α, β) space
        norm = np.sqrt(deta**2 + dalpha**2 + dbeta**2)
        u = (deta/norm, dalpha/norm, dbeta/norm)
        print(f"\n  Direction: {name}")
        print(f"    (Δη, Δα, Δβ) per unit angle = {u}")
        eig_max_dev = 0.0
        for theta in angles:
            xi = hopf_to_quat(
                eta=0 + u[0]*theta,
                alpha=delta + u[1]*theta,
                beta=0 + u[2]*theta
            )
            M = build_M(1.0, xi)
            eig = unique_eigs(M)
            if len(eig) == 3:
                dev = np.max(np.abs(eig - eig_ref))
                eig_max_dev = max(eig_max_dev, dev)
        print(f"    Max |eig − eig_ref| over sweep [0, 2π]: {eig_max_dev:.6e}")
        if eig_max_dev < 1e-6:
            print(f"    → EXACT silent direction (eigenvalues invariant on full circle)")
        else:
            print(f"    → not exactly silent at finite displacement")

# ----------------------------------------------------------------------
# Step 5.  Connection to v59 step 6 finding: Q invariant on S³
# ----------------------------------------------------------------------
print()
print("=" * 70)
print("Step 5: Random S³ sweep — Q invariance (independent check)")
print("=" * 70)
np.random.seed(42)
Qs = []
eig_spread = []
for trial in range(200):
    xi = np.random.randn(4)
    xi = xi / np.linalg.norm(xi) / np.sqrt(2)
    M = build_M(1.0, xi)
    eig = unique_eigs(M, tol=1e-7)
    if len(eig) == 3 and eig.sum() > 0:
        Q = (eig**2).sum() / eig.sum()**2
        Qs.append(Q)
        eig_spread.append(eig)
Qs = np.array(Qs)
eig_spread = np.array(eig_spread)
print(f"  {len(Qs)} random points on S³ of radius 1/√2.")
print(f"  Koide Q: mean = {Qs.mean():.10f}, σ = {Qs.std():.3e}")
print(f"  Eigenvalues per axis (min, max, range):")
for k in range(3):
    arr = eig_spread[:, k]
    print(f"    e_{k}: min={arr.min():+8.4f}, max={arr.max():+8.4f}, "
          f"range={arr.max()-arr.min():.3e}")

# ----------------------------------------------------------------------
# Step 6.  Assessment
# ----------------------------------------------------------------------
print()
print("=" * 70)
print("Step 6: Assessment")
print("=" * 70)
print(f"""
The reference Brannen phase point on S³:
    ξ = ({xi_ref[0]:.4f}, {xi_ref[1]:.4f}, {xi_ref[2]:.4f}, {xi_ref[3]:.4f})
with 3 unique eigenvalues {eig_ref}.

Random-S³ sweep shows individual eigenvalues vary BY ORDER UNITY across the
3-sphere: range of e_0 is {eig_spread[:,0].max()-eig_spread[:,0].min():.3f},
range of e_1 is {eig_spread[:,1].max()-eig_spread[:,1].min():.3f},
range of e_2 is {eig_spread[:,2].max()-eig_spread[:,2].min():.3f}.

So eigenvalues are NOT all invariant on S³ — only Q is.  The "silent
directions" (if any) at the Brannen point would need to be a 1- or 2-parameter
subset; the rest of S³ is "active" and would shift lepton masses if traversed.

Conclusion: simple Cosserat coupling to a generic direction on S³ is ruled
out (Step 1 result holds), but coupling restricted to a silent direction
might survive.  Whether that silent direction is *physically* realised by
density coupling is the next question — needs a mechanism that picks out
the silent direction (e.g. the U(1) Hopf-fiber identification with EM gauge).
""")

# Save data
np.savez('/home/d/code/scp/v59/cosserat_experiment/02_silent.npz',
         xi_ref=xi_ref, eig_ref=eig_ref, Q_ref=Q_ref,
         Qs_random=Qs, eig_spread=eig_spread,
         eps=eps,
         responses={name: r[0] for name, r in results.items()})

print("Saved data to 02_silent.npz")
