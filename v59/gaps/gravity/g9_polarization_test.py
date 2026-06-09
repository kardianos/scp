#!/usr/bin/env python3
"""
v59/gaps/gravity/g9_polarization_test.py

G9 -- SCALAR vs TENSOR.  The decisive LIGO test.

In the v59 OBE the long-range force is carried by a connection Omega(x).
Two readings of "what propagates" give very different polarization content:

  (R-scalar)  The established gravity result is  box Omega_grav = f_g rho_grav,
              rho_grav = Tr(M^dag M) = Sum m  (a SCALAR source, FINDINGS_obe_bridge
              section 4 + obe_radial_test.py).  A massless SCALAR field on R^{3,1}
              has ONE on-shell polarization: helicity h = 0.  A massless VECTOR
              (the connection A_mu, EM-like) has h = +-1.  Neither has h = +-2.
              => LIGO (which detects the transverse-traceless h = +-2 strain)
                 is NOT reproduced.  This is the gap.

  (R-tensor)  v59 hopes a genuine spin-2 mode hides in the Lambda^2 = so(8)
              bivector structure of L.  This script tests, by EXPLICIT
              polarization decomposition, whether the bivector connection
              Omega in Lambda^2 -- when promoted to a SPACETIME field and
              coupled to a test excitation -- can carry a transverse-traceless
              (h = +-2) mode.

The method (helicity / little-group SO(2) decomposition for a massless mode
propagating along z):
  * A massless field's physical polarizations are labelled by the helicity
    eigenvalue under rotations about the propagation axis z.
  * We build the SO(2)_z rotation generator on each candidate field space,
    diagonalize it, and read off the helicity content {h}.
  * A spin-2 GRADED mode requires |h| = 2 to be present AND to be
    transverse-traceless (the two LIGO polarizations h_+ , h_x).

We test, in order:
  1. a spacetime SCALAR  S(x)          -> {0}
  2. a spacetime VECTOR  A_mu(x)       -> {0(x2 gauge/long), +-1}  physical {+-1}
  3. a spacetime symmetric tensor h_mu_nu (the graviton)   -> contains +-2
  4. the v59 object: Omega is an INTERNAL bivector (so(8), 28-dim) carried by a
     spacetime SCALAR amplitude (box Omega = rho).  Its spacetime helicity is
     that of the carrier = 0.  The internal index does NOT add spacetime
     helicity.  <-- the crux.
  5. could a spacetime BIVECTOR F_mu_nu (a 2-form, like EM field strength)
     supply +-2?  Decompose a real antisymmetric spacetime 2-form. -> {+-1, 0}
     (NO +-2; antisymmetric rank-2 is spin-1 content, not spin-2).
  6. the ONLY route: a SYMMETRIC spacetime tensor.  Can the so(8) connection
     source one?  Test whether the bilinear  Omega_[mu a] Omega_[nu]^a  (two
     spacetime-vector-valued bivector legs contracted on the internal index)
     yields a symmetric TT spacetime tensor with an h=+-2 part.
"""

import numpy as np

np.set_printoptions(precision=4, suppress=True)

# ---------------------------------------------------------------------------
# Helicity machinery: rotation about z by angle ph acts on a field carrying
# given spacetime indices.  We build the generator J_z = d/dph |_{ph=0} of the
# representation and read its imaginary eigenvalues (= helicities h).
# ---------------------------------------------------------------------------

def Rz(ph):
    """4x4 spatial rotation about z (Minkowski indices t,x,y,z)."""
    c, s = np.cos(ph), np.sin(ph)
    return np.array([[1, 0, 0, 0],
                     [0, c, -s, 0],
                     [0, s,  c, 0],
                     [0, 0,  0, 1]], dtype=float)

def Jz_vector():
    """Generator on a 4-vector V_mu."""
    eps = 1e-6
    return (Rz(eps) - Rz(-eps)) / (2 * eps)

def helicities_from_generator(J):
    """Eigenvalues of -i J are the helicities (J is real antisymmetric in the
    rotation plane => eigenvalues +-i*h)."""
    w = np.linalg.eigvals(J)
    # helicities are imag parts (the generator J has eigenvalues i*h)
    h = np.imag(w)
    return np.round(np.sort(h), 6)


print("=" * 74)
print("G9: POLARIZATION CONTENT of candidate propagating modes (helicity along z)")
print("=" * 74)

# 1. SCALAR
print("\n[1] spacetime SCALAR S(x):  helicities = {0}")
print("    A massless scalar (box S = rho) has ONE polarization, h = 0.")
print("    => This is what  box Omega_grav = f_g rho_grav  gives if Omega's")
print("       spacetime carrier is a scalar amplitude.  NO h=+-2.")

# 2. VECTOR
Jz_v = Jz_vector()
hv = helicities_from_generator(Jz_v)
print(f"\n[2] spacetime VECTOR A_mu:  helicities = {hv}")
print("    Physical (massless, transverse) = {+1, -1}; the h=0 entries are")
print("    the timelike + longitudinal (gauge / non-propagating).  NO h=+-2.")

# 3. SYMMETRIC RANK-2 TENSOR (the graviton h_mu_nu) on the spatial 2x2 transverse block
# Build J_z acting on symmetric 2-tensors over the transverse plane (x,y).
def Jz_on_sym2_transverse():
    """Generator on symmetric traceless 2x2 tensors in the (x,y) plane.
    Basis: h_+ = diag(1,-1)/sqrt2 ,  h_x = [[0,1],[1,0]]/sqrt2 ."""
    eps = 1e-6
    R2 = lambda a: np.array([[np.cos(a), -np.sin(a)], [np.sin(a), np.cos(a)]])
    basis = [np.array([[1, 0], [0, -1]]) / np.sqrt(2),
             np.array([[0, 1], [1, 0]]) / np.sqrt(2)]
    def rot_tensor(M, a):
        R = R2(a)
        return R @ M @ R.T
    J = np.zeros((2, 2))
    for j, B in enumerate(basis):
        dM = (rot_tensor(B, eps) - rot_tensor(B, -eps)) / (2 * eps)
        # express dM in the basis
        for i, A in enumerate(basis):
            J[i, j] = np.sum(A * dM)
    return J

Jz_t = Jz_on_sym2_transverse()
ht = helicities_from_generator(Jz_t)
print(f"\n[3] symmetric traceless transverse TENSOR (h_+, h_x):  helicities = {ht}")
print("    => contains h = +-2.  THIS is the LIGO spin-2 / graviton signature.")
print("    (h_+ and h_x are the LIGO polarizations; they form an h=+-2 doublet.)")

# 4. v59 object: INTERNAL so(8) bivector carried by a SCALAR spacetime amplitude
print("\n[4] v59 OBE object:  Omega(x) = (so(8)-bivector) x (scalar carrier)")
print("    box Omega^{ab} = f_g rho_grav delta^{ab}-projected source.")
print("    Spacetime helicity = helicity of the CARRIER. rho_grav = Tr(M^dag M)")
print("    is a Lorentz SCALAR (a generation-space trace), so the carrier is a")
print("    scalar => helicities = {0}.  The internal a,b indices are SO(8)")
print("    (internal), NOT spacetime Lorentz, so they add NO spacetime helicity.")
print("    => v59 gravity-as-density is h=0 SCALAR.  Confirms the gap.")

# 5. spacetime BIVECTOR / 2-form F_mu_nu (could the connection be a spacetime 2-form?)
def Jz_on_2form():
    """Generator on antisymmetric spacetime 2-forms F_mu_nu (6-dim: tx,ty,tz,xy,xz,yz)."""
    eps = 1e-6
    comps = [(0,1),(0,2),(0,3),(1,2),(1,3),(2,3)]
    def rot_2form_vec(vec, a):
        # reconstruct F, rotate, flatten
        F = np.zeros((4,4))
        for k,(mu,nu) in enumerate(comps):
            F[mu,nu] = vec[k]; F[nu,mu] = -vec[k]
        R = Rz(a)
        Fp = R @ F @ R.T
        return np.array([Fp[mu,nu] for (mu,nu) in comps])
    J = np.zeros((6,6))
    for k in range(6):
        e = np.zeros(6); e[k]=1.0
        dvec = (rot_2form_vec(e, eps) - rot_2form_vec(e, -eps))/(2*eps)
        J[:,k] = dvec
    return J

Jz_F = Jz_on_2form()
hF = helicities_from_generator(Jz_F)
print(f"\n[5] spacetime 2-form / bivector F_mu_nu:  helicities = {hF}")
print("    An antisymmetric rank-2 spacetime tensor decomposes as spin-1")
print("    content: helicities {+-1, 0,...}.  NO h=+-2.  A spacetime BIVECTOR")
print("    connection (EM-like) therefore CANNOT carry the LIGO tensor mode.")

# 6. The only route: a SYMMETRIC spacetime tensor sourced by the bivector.
#    Test: take Omega as a spacetime-VECTOR-valued internal bivector,
#    Omega_mu^{a} (a = internal so(8) index, mu = spacetime). Build the
#    gauge-invariant bilinear  T_mu_nu = Omega_mu^a Omega_nu^a (internal trace),
#    which IS symmetric in (mu,nu).  Does its transverse part contain h=+-2?
print("\n[6] BILINEAR route:  T_mu_nu = Omega_mu^a Omega_nu^a  (sum over internal a)")
print("    This is symmetric in (mu,nu) by construction (a stress-tensor-like")
print("    bilinear). A symmetric spacetime 2-tensor's TRANSVERSE-TRACELESS part")
print("    carries h=+-2.  So IF the connection has a spacetime-VECTOR index mu")
print("    (Omega_mu^a, a Yang-Mills-like connection 1-form valued in so(8)),")
print("    its bilinear can source a spin-2 TT tensor.")
print("    BUT: T_mu_nu is QUADRATIC in Omega -> it is the *energy-momentum* of")
print("    the gauge field, i.e. it sources gravity as a SOURCE, it is not a")
print("    propagating massless graviton itself. Quadratic => O(g^2)=O(alpha)")
print("    suppressed and short-range (no massless pole from a composite).")

# Numerically confirm the symmetric-tensor helicity content one more way:
# random symmetric traceless transverse tensor, project to helicity basis.
print("\n--- numeric confirmation: full symmetric 4x4 tensor h_mu_nu ---")
def Jz_on_sym4():
    eps = 1e-6
    # basis of symmetric 4x4 (10-dim)
    idx = [(i,j) for i in range(4) for j in range(i,4)]
    def vec_to_M(v):
        M = np.zeros((4,4))
        for k,(i,j) in enumerate(idx):
            M[i,j]=v[k]; M[j,i]=v[k]
        return M
    def M_to_vec(M):
        return np.array([M[i,j] for (i,j) in idx])
    J = np.zeros((10,10))
    for k in range(10):
        e=np.zeros(10); e[k]=1.0
        M=vec_to_M(e)
        R=Rz(eps); Rm=Rz(-eps)
        dM=(R@M@R.T - Rm@M@Rm.T)/(2*eps)
        J[:,k]=M_to_vec(dM)
    return J
hsym = helicities_from_generator(Jz_on_sym4())
print(f"  symmetric 4x4 tensor helicities = {hsym}")
print("  => the full symmetric tensor DOES contain h=+-2 (the +-2 entries).")
print("     The transverse-traceless projection isolates exactly {+2,-2} (LIGO).")

print("\n" + "=" * 74)
print("VERDICT on G9 (scalar vs tensor)")
print("=" * 74)
print("- A massless SCALAR carrier (box Omega = rho_grav, rho_grav a Lorentz")
print("  scalar) => h=0 ONLY.  This is the v59 established gravity mode. NO LIGO.")
print("- The INTERNAL so(8)/Lambda^2 index does NOT supply spacetime helicity:")
print("  helicity is a SPACETIME (little-group) label; internal indices are inert")
print("  under spatial rotations. So 'gravity lives in the 28-dim bivector L'")
print("  does NOT by itself give h=+-2.")
print("- A spacetime 2-form/bivector connection (EM-like, F_mu_nu) carries h=+-1,")
print("  not +-2. So even a *spacetime* bivector connection is spin-1, not spin-2.")
print("- The ONLY h=+-2 carrier is a SYMMETRIC spacetime rank-2 tensor h_mu_nu.")
print("  The bivector structure can SOURCE one quadratically (T_mu_nu =")
print("  Omega_mu^a Omega_nu^a), but that is the gauge-field stress tensor")
print("  (a source term), not a fundamental massless graviton => O(alpha) and")
print("  short-range, not the long-range LIGO mode.")
print("- CONCLUSION: as currently formulated, v59 gravity is SCALAR. To pass LIGO")
print("  it needs a genuine symmetric-tensor massless field h_mu_nu (a graviton")
print("  in the spacetime sector), which the present internal-bivector OBE does")
print("  NOT contain. This is a FATAL gap for LIGO consistency unless an")
print("  emergent/induced-metric (spacetime-symmetric) mode is added.")
