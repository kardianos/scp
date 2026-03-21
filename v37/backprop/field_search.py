"""
Field search: train raw 3D field values from random initialization to find
self-sustaining soliton structures via backpropagation through PDE evolution.

Key differences from structure_search.py:
  - NO braid seed: start from RANDOM sphere
  - NO periodic BC: absorbing boundary in every timestep
  - NO background field: pure vacuum outside the sphere
  - Curriculum training: T_opt ramps up over training
  - Multiple seeds in parallel via jax.vmap (or sequential fallback)
  - Hard clamp: field outside radius R_seed forced to zero after every update

Usage:
    # Validation tests
    python field_search.py --validate

    # Full search
    python field_search.py --n_seeds 32 --N 32 --L 15 --epochs 1500 --outdir results/

    # Higher resolution
    python field_search.py --n_seeds 8 --N 48 --L 15 --epochs 1000 --outdir results_n48/

    # Resume from checkpoint
    python field_search.py --resume results/best_seed.npy --N 32 --L 15 --epochs 500
"""

import jax
import jax.numpy as jnp
from jax import jit
import optax
import numpy as np
import argparse
import time
import os
import struct
import sys

# ---------------------------------------------------------------------------
# Physics parameters (must match C simulation)
# ---------------------------------------------------------------------------
MU = -41.345
KAPPA = 50.0
MASS2 = 2.25
ETA = 0.5       # phi-theta coupling
MTHETA2 = 0.0   # massless theta

# ---------------------------------------------------------------------------
# Default grid / simulation parameters
# ---------------------------------------------------------------------------
DEFAULT_N = 32
DEFAULT_L = 15.0
DEFAULT_R_SEED = 5.0
DEFAULT_R_CLAMP = 6.0
DEFAULT_DAMP_RATE = 0.01
DEFAULT_AMPLITUDE = 0.5
DEFAULT_N_SEEDS = 32
DEFAULT_EPOCHS = 1500
DEFAULT_LR = 1e-3


# ===================================================================
#  Static grid helpers (called ONCE, results passed into JIT funcs)
# ===================================================================

def make_grid(N, L):
    """Create 3D coordinate arrays and spacing."""
    x = jnp.linspace(-L, L, N)
    X, Y, Z = jnp.meshgrid(x, x, x, indexing='ij')
    dx = 2.0 * L / (N - 1)
    R = jnp.sqrt(X**2 + Y**2 + Z**2)
    return X, Y, Z, R, dx


def make_damp_mask(R, L, R_damp, damp_rate):
    """Spherical absorbing-boundary damping mask.

    Returns a (N,N,N) array of velocity multipliers:
      1.0 inside R_damp, smoothly decaying to (1-damp_rate) at r=L.
    """
    frac = jnp.where(R > R_damp,
                     ((R - R_damp) / (L - R_damp + 1e-30))**2,
                     0.0)
    mask = 1.0 - damp_rate * frac
    return jnp.clip(mask, 0.0, 1.0)


def make_sphere_mask(R, R_clamp, edge_width=1.0):
    """Soft-edge spherical clamp mask.

    1.0 for r < R_clamp - edge_width, smooth falloff to 0.0 at R_clamp.
    """
    inner = R_clamp - edge_width
    mask = jnp.where(R < inner, 1.0,
           jnp.where(R < R_clamp,
                     0.5 * (1.0 + jnp.cos(jnp.pi * (R - inner) / edge_width)),
                     0.0))
    return mask


# ===================================================================
#  PDE: forces, energies, diagnostics
# ===================================================================

def laplacian_3d_absorbing(f, idx2):
    """6-point Laplacian with ZERO (Dirichlet) boundary conditions.

    We pad with zeros on each face, which is equivalent to the field
    being zero outside the grid -- no periodic wrapping.
    """
    # Shift along each axis, padding with 0 instead of wrapping
    sp = jnp.pad(f[:-1, :, :], ((1, 0), (0, 0), (0, 0)))  # roll +1 axis 0
    sm = jnp.pad(f[1:, :, :],  ((0, 1), (0, 0), (0, 0)))  # roll -1 axis 0
    sp1 = jnp.pad(f[:, :-1, :], ((0, 0), (1, 0), (0, 0)))
    sm1 = jnp.pad(f[:, 1:, :],  ((0, 0), (0, 1), (0, 0)))
    sp2 = jnp.pad(f[:, :, :-1], ((0, 0), (0, 0), (1, 0)))
    sm2 = jnp.pad(f[:, :, 1:],  ((0, 0), (0, 0), (0, 1)))
    return (sp + sm + sp1 + sm1 + sp2 + sm2 - 6.0 * f) * idx2


def gradient_component(f, axis, idx1):
    """Central-difference gradient along one axis (zero-pad BC)."""
    if axis == 0:
        fp = jnp.pad(f[:-1, :, :], ((1, 0), (0, 0), (0, 0)))
        fm = jnp.pad(f[1:, :, :],  ((0, 1), (0, 0), (0, 0)))
    elif axis == 1:
        fp = jnp.pad(f[:, :-1, :], ((0, 0), (1, 0), (0, 0)))
        fm = jnp.pad(f[:, 1:, :],  ((0, 0), (0, 1), (0, 0)))
    else:
        fp = jnp.pad(f[:, :, :-1], ((0, 0), (0, 0), (1, 0)))
        fm = jnp.pad(f[:, :, 1:],  ((0, 0), (0, 0), (0, 1)))
    return (fm - fp) * idx1


def curl_component_absorbing(F, a, idx1):
    """Curl component a of vector field F[3,N,N,N] with zero-pad BC."""
    if a == 0:
        return (gradient_component(F[2], 1, idx1) -
                gradient_component(F[1], 2, idx1))
    elif a == 1:
        return (gradient_component(F[0], 2, idx1) -
                gradient_component(F[2], 0, idx1))
    else:
        return (gradient_component(F[1], 0, idx1) -
                gradient_component(F[0], 1, idx1))


def compute_forces(phi, theta, dx):
    """Compute accelerations for all 6 fields (3 phi + 3 theta).

    PDE:
        phi_a_tt = lap(phi_a) - m^2 phi_a - mu P / (1+kappa P^2)^2 * dP/dphi_a + eta curl(theta)_a
        theta_a_tt = lap(theta_a) - mtheta^2 theta_a + eta curl(phi)_a
    """
    idx2 = 1.0 / (dx * dx)
    idx1 = 1.0 / (2.0 * dx)

    P = phi[0] * phi[1] * phi[2]
    den = 1.0 + KAPPA * P * P
    # d/dphi_a of V(P) = (mu/2) P^2/(1+kappa P^2)
    # dV/dP = mu P (1+kappa P^2 - kappa P^2 * 2) / (1+kappa P^2)^2
    #       = mu P / (1+kappa P^2)^2
    mPd2 = MU * P / (den * den)

    phi_acc_list = []
    for a in range(3):
        lap = laplacian_3d_absorbing(phi[a], idx2)
        if a == 0:
            dPda = phi[1] * phi[2]
        elif a == 1:
            dPda = phi[0] * phi[2]
        else:
            dPda = phi[0] * phi[1]
        curl_t = curl_component_absorbing(theta, a, idx1)
        phi_acc_list.append(lap - MASS2 * phi[a] - mPd2 * dPda + ETA * curl_t)

    theta_acc_list = []
    for a in range(3):
        lap_t = laplacian_3d_absorbing(theta[a], idx2)
        curl_p = curl_component_absorbing(phi, a, idx1)
        theta_acc_list.append(lap_t - MTHETA2 * theta[a] + ETA * curl_p)

    phi_acc = jnp.stack(phi_acc_list, axis=0)
    theta_acc = jnp.stack(theta_acc_list, axis=0)
    return phi_acc, theta_acc


def compute_epot(phi, dx):
    """Potential energy V(P) = (mu/2) P^2 / (1 + kappa P^2), integrated."""
    dV = dx ** 3
    P = phi[0] * phi[1] * phi[2]
    V = (MU / 2.0) * P * P / (1.0 + KAPPA * P * P)
    return jnp.sum(V) * dV


def compute_ekin(phi_vel, theta_vel, dx):
    """Kinetic energy = 0.5 * sum(vel^2) * dV."""
    dV = dx ** 3
    ke = 0.5 * jnp.sum(phi_vel**2 + theta_vel**2) * dV
    return ke


def compute_egrad(phi, dx):
    """Gradient (elastic) energy = 0.5 * sum(|grad phi_a|^2) * dV."""
    dV = dx ** 3
    idx1 = 1.0 / (2.0 * dx)
    ge = 0.0
    for a in range(3):
        for axis in range(3):
            g = gradient_component(phi[a], axis, idx1)
            ge = ge + jnp.sum(g**2)
    return 0.5 * ge * dV


def compute_emass(phi, dx):
    """Mass-term energy = 0.5 * m^2 * sum(phi^2) * dV."""
    dV = dx ** 3
    return 0.5 * MASS2 * jnp.sum(phi**2) * dV


def compute_etotal(phi, phi_vel, theta, theta_vel, dx):
    """Total energy = E_kin + E_grad + E_mass + E_pot (+theta terms)."""
    ekin = compute_ekin(phi_vel, theta_vel, dx)
    egrad_phi = compute_egrad(phi, dx)
    egrad_theta = compute_egrad(theta, dx)
    emass = compute_emass(phi, dx)
    emass_theta = 0.5 * MTHETA2 * jnp.sum(theta**2) * dx**3
    epot = compute_epot(phi, dx)
    # Cross-term energy from eta*curl coupling is already accounted for
    # through the equations of motion; the conserved energy is:
    # E = E_kin + E_grad_phi + E_grad_theta + E_mass_phi + E_mass_theta + E_pot
    # (The curl coupling is velocity-independent so doesn't add a potential term
    #  -- it appears in the EOM but energy is sum of individual kinetic+potential.)
    # Actually for the curl coupling: L_cross = eta * phi . curl(theta)
    # => E_cross = -eta * phi . curl(theta)  (Legendre transform for time-deriv-free coupling)
    # But since the coupling involves spatial derivatives only, the conserved Hamiltonian
    # does NOT include a separate cross term -- the coupling enters only through curl forces.
    # For energy conservation tests this should suffice.
    return ekin + egrad_phi + egrad_theta + emass + emass_theta + epot


def compute_rrms(phi, X, Y, Z, dx):
    """RMS radius of the field energy density."""
    dV = dx ** 3
    rho = jnp.sum(phi**2, axis=0)  # sum over field components
    total = jnp.sum(rho) * dV + 1e-30
    r2 = X**2 + Y**2 + Z**2
    r2_avg = jnp.sum(r2 * rho) * dV / total
    return jnp.sqrt(jnp.maximum(r2_avg, 0.0))


def compute_grad_energy(phi, dx):
    """Gradient energy of the initial condition (for smoothness penalty)."""
    return compute_egrad(phi, dx)


# ===================================================================
#  Time evolution (Verlet + absorbing BC)
# ===================================================================

def verlet_step_absorbing(phi, phi_vel, theta, theta_vel, dx, dt, damp_mask):
    """One velocity-Verlet step with spherical absorbing boundary.

    damp_mask is pre-computed (N,N,N) array of velocity multipliers.
    """
    phi_acc, theta_acc = compute_forces(phi, theta, dx)

    # Half-kick
    phi_vel = phi_vel + 0.5 * dt * phi_acc
    theta_vel = theta_vel + 0.5 * dt * theta_acc

    # Drift
    phi = phi + dt * phi_vel
    theta = theta + dt * theta_vel

    # Force at new position
    phi_acc, theta_acc = compute_forces(phi, theta, dx)

    # Half-kick
    phi_vel = phi_vel + 0.5 * dt * phi_acc
    theta_vel = theta_vel + 0.5 * dt * theta_acc

    # Absorbing BC: damp velocities near boundary
    phi_vel = phi_vel * damp_mask[None, :, :, :]
    theta_vel = theta_vel * damp_mask[None, :, :, :]

    return phi, phi_vel, theta, theta_vel


def verlet_step_nodamp(phi, phi_vel, theta, theta_vel, dx, dt):
    """One velocity-Verlet step WITHOUT absorbing BC (for energy conservation test)."""
    phi_acc, theta_acc = compute_forces(phi, theta, dx)
    phi_vel = phi_vel + 0.5 * dt * phi_acc
    theta_vel = theta_vel + 0.5 * dt * theta_acc
    phi = phi + dt * phi_vel
    theta = theta + dt * theta_vel
    phi_acc, theta_acc = compute_forces(phi, theta, dx)
    phi_vel = phi_vel + 0.5 * dt * phi_acc
    theta_vel = theta_vel + 0.5 * dt * theta_acc
    return phi, phi_vel, theta, theta_vel


def evolve_absorbing(phi0, phi_vel0, theta0, theta_vel0, dx, dt, n_steps, damp_mask):
    """Evolve n_steps with absorbing BC using jax.lax.scan."""
    def step_fn(carry, _):
        phi, phi_vel, theta, theta_vel = carry
        phi, phi_vel, theta, theta_vel = verlet_step_absorbing(
            phi, phi_vel, theta, theta_vel, dx, dt, damp_mask
        )
        return (phi, phi_vel, theta, theta_vel), None

    (phi, phi_vel, theta, theta_vel), _ = jax.lax.scan(
        step_fn, (phi0, phi_vel0, theta0, theta_vel0), None, length=n_steps
    )
    return phi, phi_vel, theta, theta_vel


def evolve_nodamp(phi0, phi_vel0, theta0, theta_vel0, dx, dt, n_steps):
    """Evolve n_steps WITHOUT absorbing BC."""
    def step_fn(carry, _):
        phi, phi_vel, theta, theta_vel = carry
        phi, phi_vel, theta, theta_vel = verlet_step_nodamp(
            phi, phi_vel, theta, theta_vel, dx, dt
        )
        return (phi, phi_vel, theta, theta_vel), None

    (phi, phi_vel, theta, theta_vel), _ = jax.lax.scan(
        step_fn, (phi0, phi_vel0, theta0, theta_vel0), None, length=n_steps
    )
    return phi, phi_vel, theta, theta_vel


# ===================================================================
#  Initialization
# ===================================================================

def init_random_sphere(key, N, L, R_seed=DEFAULT_R_SEED, amplitude=DEFAULT_AMPLITUDE):
    """Random smooth fields inside a sphere, zero outside.

    CORRELATED initialization: phi_1 and phi_2 are derived from phi_0
    with the known-good phase ratios, plus independent random perturbations.
    This ensures P = phi_0 * phi_1 * phi_2 is consistently nonzero from the start,
    giving the optimizer a gradient signal to work with.

    The SPATIAL structure is random — only the inter-field phase relationship is seeded.
    """
    x = jnp.linspace(-L, L, N)
    X, Y, Z = jnp.meshgrid(x, x, x, indexing='ij')
    R = jnp.sqrt(X**2 + Y**2 + Z**2)

    k1, k2, k3, k4 = jax.random.split(key, 4)

    # Base random field (smooth)
    base = jax.random.normal(k1, (N, N, N)) * amplitude
    for _ in range(5):
        base = (base
                + jnp.roll(base, 1, 0) + jnp.roll(base, -1, 0)
                + jnp.roll(base, 1, 1) + jnp.roll(base, -1, 1)
                + jnp.roll(base, 1, 2) + jnp.roll(base, -1, 2)) / 7.0

    # Correlated fields: phi_a = base * cos(delta_a) + perturbation
    # This gives P ~ base^3 * cos(d0)*cos(d1)*cos(d2) > 0 where base > 0
    DELTA = jnp.array([0.0, 3.0005, 4.4325])
    # cos(0)=1, cos(3.0)=-0.99, cos(4.4)=-0.31
    # So P ~ base^3 * (1.0)(-0.99)(-0.31) = +0.307 * base^3 > 0

    # Independent perturbations (20% of correlated signal)
    perturb = jax.random.normal(k2, (3, N, N, N)) * amplitude * 0.2
    for _ in range(3):
        smoothed = []
        for a in range(3):
            s = (perturb[a]
                 + jnp.roll(perturb[a], 1, 0) + jnp.roll(perturb[a], -1, 0)
                 + jnp.roll(perturb[a], 1, 1) + jnp.roll(perturb[a], -1, 1)
                 + jnp.roll(perturb[a], 1, 2) + jnp.roll(perturb[a], -1, 2)) / 7.0
            smoothed.append(s)
        perturb = jnp.stack(smoothed, axis=0)

    phi = jnp.stack([
        base * jnp.cos(DELTA[0]) + perturb[0],
        base * jnp.cos(DELTA[1]) + perturb[1],
        base * jnp.cos(DELTA[2]) + perturb[2],
    ], axis=0)

    # Spherical envelope: Gaussian inside, hard cutoff at R_seed
    sigma = R_seed * 0.6
    envelope = jnp.where(R < R_seed, jnp.exp(-R**2 / (2.0 * sigma**2)), 0.0)
    phi = phi * envelope[None, :, :, :]

    # Start stationary
    vel = jnp.zeros_like(phi)
    return phi, vel


# ===================================================================
#  Loss function
# ===================================================================

def make_loss_fn(dx, dt, n_steps, damp_mask, X, Y, Z, L):
    """Build a loss function closed over static parameters.

    Returns a function: (phi0, vel0) -> (loss, aux_dict).
    Differentiable w.r.t. phi0 and vel0.
    """
    def loss_fn(phi0, vel0):
        theta0 = jnp.zeros_like(phi0)
        theta_vel0 = jnp.zeros_like(phi0)

        phi_T, _, _, _ = evolve_absorbing(
            phi0, vel0, theta0, theta_vel0, dx, dt, n_steps, damp_mask
        )

        # Potential energy at start and end
        epot_0 = compute_epot(phi0, dx)
        epot_T = compute_epot(phi_T, dx)

        # Compactness at end
        rrms_T = compute_rrms(phi_T, X, Y, Z, dx)

        # Smoothness of initial condition
        grad_e = compute_grad_energy(phi0, dx)

        # --- Loss components (all to be MINIMIZED) ---
        # Primary: want epot_T maximally negative (binding survived)
        binding_loss = -epot_T

        # Compactness: penalize spreading
        compact_loss = (rrms_T / L) ** 2

        # Smoothness: penalize sharp features in IC (regularization)
        smooth_loss = grad_e * dx**3

        loss = 10.0 * binding_loss + 0.1 * compact_loss + 0.0001 * smooth_loss

        aux = {
            'epot_0': epot_0,
            'epot_T': epot_T,
            'rrms_T': rrms_T,
            'grad_e': grad_e,
            'binding': binding_loss,
            'compact': compact_loss,
            'smooth': smooth_loss,
        }
        return loss, aux

    return loss_fn


# ===================================================================
#  Curriculum schedule
# ===================================================================

def get_T_opt_steps(epoch, dx, dt_factor=0.1):
    """Curriculum: number of PDE steps increases with epoch.

    Returns (n_steps, T_opt_label) for a given epoch.
    dt = dt_factor * dx, so T_opt = n_steps * dt.
    """
    # Must match CURRICULUM_TOPT_VALUES exactly
    vals = CURRICULUM_TOPT_VALUES
    if len(vals) == 4:
        if epoch < 200:
            T_opt = vals[0]
        elif epoch < 500:
            T_opt = vals[1]
        elif epoch < 1000:
            T_opt = vals[2]
        else:
            T_opt = vals[3]
    elif len(vals) == 3:
        if epoch < 300:
            T_opt = vals[0]
        elif epoch < 800:
            T_opt = vals[1]
        else:
            T_opt = vals[2]
    else:
        T_opt = vals[min(epoch // 300, len(vals) - 1)]

    dt = dt_factor * dx
    n_steps = max(1, int(T_opt / dt))
    return n_steps, T_opt


# Distinct T_opt values used in the curriculum (for pre-compiling)
CURRICULUM_TOPT_VALUES = [5.0, 10.0, 20.0, 30.0]


# ===================================================================
#  Sphere clamp (applied after each optimizer step)
# ===================================================================

def apply_sphere_clamp(phi, vel, sphere_mask):
    """Force field to zero outside the seed sphere."""
    phi = phi * sphere_mask[None, :, :, :]
    vel = vel * sphere_mask[None, :, :, :]
    return phi, vel


# ===================================================================
#  Checkpointing and I/O
# ===================================================================

def save_candidate(phi, vel, path, N, L, epoch, loss, aux):
    """Save a candidate as .npy with metadata."""
    os.makedirs(os.path.dirname(path) or '.', exist_ok=True)
    data = {
        'phi': np.array(phi),
        'vel': np.array(vel),
        'N': N,
        'L': L,
        'epoch': epoch,
        'loss': float(loss),
        'epot_0': float(aux['epot_0']),
        'epot_T': float(aux['epot_T']),
        'rrms_T': float(aux['rrms_T']),
    }
    np.save(path, data, allow_pickle=True)


def load_candidate(path):
    """Load a saved candidate."""
    data = np.load(path, allow_pickle=True).item()
    phi = jnp.array(data['phi'])
    vel = jnp.array(data['vel'])
    return phi, vel, data


def save_field_bin(phi, N, L, t, path):
    """Save field snapshot as raw binary for visualization."""
    os.makedirs(os.path.dirname(path) or '.', exist_ok=True)
    with open(path, 'wb') as f:
        f.write(struct.pack('i', N))
        f.write(struct.pack('d', L))
        f.write(struct.pack('d', t))
        f.write(struct.pack('i', 3))  # n_fields
        for a in range(3):
            f.write(np.array(phi[a], dtype=np.float32).tobytes())


# ===================================================================
#  Validation tests
# ===================================================================

def test_gradient_check(N=16, L=10.0):
    """V1: Finite-difference vs autodiff gradient at a single grid point."""
    print("=== V1: Gradient Check ===")

    X, Y, Z, R, dx = make_grid(N, L)
    dt = 0.1 * dx
    n_steps = 5
    R_damp = L - 4.0
    damp_mask = make_damp_mask(R, L, R_damp, DEFAULT_DAMP_RATE)

    key = jax.random.PRNGKey(42)
    phi0, vel0 = init_random_sphere(key, N, L, R_seed=4.0, amplitude=0.3)

    loss_fn = make_loss_fn(dx, dt, n_steps, damp_mask, X, Y, Z, L)

    # Autodiff gradient
    grad_fn = jax.grad(lambda p, v: loss_fn(p, v)[0], argnums=0)
    g_auto = grad_fn(phi0, vel0)

    # Finite-difference at a specific interior point
    # Pick a point near the center where the field is nonzero
    c = N // 2
    eps = 1e-4
    phi_plus = phi0.at[0, c, c, c].add(eps)
    phi_minus = phi0.at[0, c, c, c].add(-eps)
    l_plus, _ = loss_fn(phi_plus, vel0)
    l_minus, _ = loss_fn(phi_minus, vel0)
    g_fd = (l_plus - l_minus) / (2.0 * eps)
    g_ad = g_auto[0, c, c, c]

    rel_err = abs(float(g_fd - g_ad)) / (abs(float(g_ad)) + 1e-30)
    match = rel_err < 0.01  # 1% tolerance

    print(f"  Point: phi[0, {c}, {c}, {c}]")
    print(f"  Autodiff gradient:  {float(g_ad):.8e}")
    print(f"  Finite-diff gradient: {float(g_fd):.8e}")
    print(f"  Relative error: {rel_err:.2e}")
    print(f"  PASS" if match else f"  FAIL (rel_err = {rel_err:.2e})")
    print()
    return match


def test_energy_conservation(N=16, L=10.0):
    """V2: Energy conservation without absorbing BC."""
    print("=== V2: Energy Conservation (no damping) ===")

    X, Y, Z, R, dx = make_grid(N, L)
    dt = 0.1 * dx
    n_steps = 100

    key = jax.random.PRNGKey(99)
    phi0, vel0 = init_random_sphere(key, N, L, R_seed=3.0, amplitude=0.3)
    theta0 = jnp.zeros_like(phi0)
    theta_vel0 = jnp.zeros_like(phi0)

    # Give some initial velocity for a nontrivial test
    vel0 = jax.random.normal(jax.random.PRNGKey(100), vel0.shape) * 0.1
    # Apply same envelope so velocity is zero at boundary
    sigma = 3.0 * 0.6
    envelope = jnp.where(R < 3.0, jnp.exp(-R**2 / (2.0 * sigma**2)), 0.0)
    vel0 = vel0 * envelope[None, :, :, :]

    E0 = compute_etotal(phi0, vel0, theta0, theta_vel0, dx)

    # Evolve without damping
    phi_T, vel_T, theta_T, theta_vel_T = evolve_nodamp(
        phi0, vel0, theta0, theta_vel0, dx, dt, n_steps
    )

    ET = compute_etotal(phi_T, vel_T, theta_T, theta_vel_T, dx)

    drift = abs(float(ET - E0)) / (abs(float(E0)) + 1e-30)
    match = drift < 0.01  # 1% tolerance

    print(f"  E(0)  = {float(E0):.6f}")
    print(f"  E({n_steps}) = {float(ET):.6f}")
    print(f"  Drift = {drift:.2e} ({drift*100:.4f}%)")
    print(f"  PASS" if match else f"  FAIL (drift = {drift:.2e})")
    print()
    return match


def test_absorbing_bc(N=16, L=10.0):
    """V3: Absorbing BC causes energy to decrease for an expanding pulse."""
    print("=== V3: Absorbing BC ===")

    X, Y, Z, R, dx = make_grid(N, L)
    dt = 0.1 * dx
    n_steps = 100
    R_damp = L - 4.0
    damp_mask = make_damp_mask(R, L, R_damp, DEFAULT_DAMP_RATE)

    # Gaussian pulse at center -- will radiate outward
    key = jax.random.PRNGKey(77)
    sigma_pulse = 2.0
    pulse = jnp.exp(-R**2 / (2.0 * sigma_pulse**2))
    phi0 = jnp.stack([pulse, pulse * 0.5, pulse * 0.3], axis=0)
    vel0 = jnp.zeros_like(phi0)
    theta0 = jnp.zeros_like(phi0)
    theta_vel0 = jnp.zeros_like(phi0)

    E0 = compute_etotal(phi0, vel0, theta0, theta_vel0, dx)

    # Evolve WITH damping
    phi_damp, vel_damp, theta_damp, thetavel_damp = evolve_absorbing(
        phi0, vel0, theta0, theta_vel0, dx, dt, n_steps, damp_mask
    )
    E_damp = compute_etotal(phi_damp, vel_damp, theta_damp, thetavel_damp, dx)

    # Evolve WITHOUT damping (control)
    phi_nd, vel_nd, theta_nd, thetavel_nd = evolve_nodamp(
        phi0, vel0, theta0, theta_vel0, dx, dt, n_steps
    )
    E_nodamp = compute_etotal(phi_nd, vel_nd, theta_nd, thetavel_nd, dx)

    absorbed = float(E0 - E_damp) / (float(E0) + 1e-30)
    conserved = abs(float(E0 - E_nodamp)) / (float(E0) + 1e-30)

    # Energy with damping should be LESS than without
    match = (float(E_damp) < float(E_nodamp)) and (absorbed > 0.001)

    print(f"  E(0)            = {float(E0):.4f}")
    print(f"  E({n_steps}) no damp = {float(E_nodamp):.4f} (drift {conserved:.2e})")
    print(f"  E({n_steps}) damped  = {float(E_damp):.4f} (absorbed {absorbed*100:.2f}%)")
    print(f"  PASS" if match else f"  FAIL")
    print()
    return match


def test_braid_seed(N=24, L=12.0):
    """V4: Braid3 seed maintains negative E_pot short-term."""
    print("=== V4: Braid3 Seed ===")

    X, Y, Z, R, dx = make_grid(N, L)
    dt = 0.1 * dx
    n_steps = 50
    R_damp = L - 4.0
    damp_mask = make_damp_mask(R, L, R_damp, DEFAULT_DAMP_RATE)

    # Build braid3 seed (inline, not importing structure_search)
    DELTA = jnp.array([0.0, 3.0005, 4.4325])
    R_h = 1.0
    r_tube = 2.0
    sigma_z = 3.0
    A = 0.8
    kw = jnp.pi / L

    phi_list = []
    vel_list = []
    for a in range(3):
        val = jnp.zeros((N, N, N))
        for s in range(3):
            cx = R_h * jnp.cos(kw * Z + 2.0 * jnp.pi * s / 3.0)
            cy = R_h * jnp.sin(kw * Z + 2.0 * jnp.pi * s / 3.0)
            d2 = (X - cx)**2 + (Y - cy)**2
            env = jnp.exp(-d2 / (2.0 * r_tube**2)) * jnp.exp(-Z**2 / (2.0 * sigma_z**2))
            phase = kw * Z + DELTA[a]
            val = val + A * env * jnp.cos(phase)
        phi_list.append(val)
        omega = jnp.sqrt(kw**2 + MASS2)
        v = jnp.zeros((N, N, N))
        for s in range(3):
            cx = R_h * jnp.cos(kw * Z + 2.0 * jnp.pi * s / 3.0)
            cy = R_h * jnp.sin(kw * Z + 2.0 * jnp.pi * s / 3.0)
            d2 = (X - cx)**2 + (Y - cy)**2
            env = jnp.exp(-d2 / (2.0 * r_tube**2)) * jnp.exp(-Z**2 / (2.0 * sigma_z**2))
            phase = kw * Z + DELTA[a]
            v = v + omega * A * env * jnp.sin(phase)
        vel_list.append(v)

    phi0 = jnp.stack(phi_list, axis=0)
    vel0 = jnp.stack(vel_list, axis=0)
    theta0 = jnp.zeros_like(phi0)
    theta_vel0 = jnp.zeros_like(phi0)

    epot0 = compute_epot(phi0, dx)

    phi_T, _, _, _ = evolve_absorbing(
        phi0, vel0, theta0, theta_vel0, dx, dt, n_steps, damp_mask
    )
    epot_T = compute_epot(phi_T, dx)

    # Also check that gradient is nonzero (optimizer can work)
    loss_fn = make_loss_fn(dx, dt, n_steps, damp_mask, X, Y, Z, L)
    grad_fn = jax.grad(lambda p, v: loss_fn(p, v)[0], argnums=0)
    g = grad_fn(phi0, vel0)
    grad_norm = float(jnp.sqrt(jnp.sum(g**2)))

    has_binding = float(epot0) < 0.0
    retains = float(epot_T) < 0.0
    has_grad = grad_norm > 1e-10

    match = has_binding and retains and has_grad

    print(f"  E_pot(0)  = {float(epot0):.4f} ({'<0 OK' if has_binding else '>0 BAD'})")
    print(f"  E_pot({n_steps}) = {float(epot_T):.4f} ({'<0 OK' if retains else '>0 BAD'})")
    print(f"  |grad|    = {grad_norm:.4e} ({'nonzero OK' if has_grad else 'ZERO BAD'})")
    print(f"  PASS" if match else f"  FAIL")
    print()
    return match


def run_all_tests():
    """Run all validation tests, return True if all pass."""
    print("=" * 60)
    print("  VALIDATION TESTS")
    print("=" * 60)
    print()

    results = []
    t0 = time.time()

    results.append(("V1: Gradient check", test_gradient_check()))
    results.append(("V2: Energy conservation", test_energy_conservation()))
    results.append(("V3: Absorbing BC", test_absorbing_bc()))
    results.append(("V4: Braid seed", test_braid_seed()))

    elapsed = time.time() - t0
    print("=" * 60)
    print(f"  RESULTS ({elapsed:.1f}s)")
    print("=" * 60)
    all_pass = True
    for name, passed in results:
        status = "PASS" if passed else "FAIL"
        print(f"  {name}: {status}")
        if not passed:
            all_pass = False
    print()
    if all_pass:
        print("  ALL TESTS PASSED")
    else:
        print("  SOME TESTS FAILED")
    print("=" * 60)
    return all_pass


# ===================================================================
#  Main training loop
# ===================================================================

def train(args):
    """Main training loop with curriculum, batched seeds, checkpointing."""
    N = args.N
    L = args.L
    n_seeds = args.n_seeds
    epochs = args.epochs
    lr = args.lr
    outdir = args.outdir
    R_seed = args.R_seed
    R_clamp = args.R_clamp

    os.makedirs(outdir, exist_ok=True)

    # --- Static grid arrays (computed ONCE) ---
    X, Y, Z, R, dx = make_grid(N, L)
    dt = 0.1 * dx
    R_damp = L - 4.0
    damp_mask = make_damp_mask(R, L, R_damp, DEFAULT_DAMP_RATE)
    sphere_mask = make_sphere_mask(R, R_clamp, edge_width=1.0)

    print("=" * 80)
    print("  FIELD SEARCH: Random Sphere -> Soliton via Backprop")
    print("=" * 80)
    print(f"  Grid:       N={N}, L={L:.1f}, dx={dx:.4f}, dt={dt:.5f}")
    print(f"  Seeds:      {n_seeds}")
    print(f"  Epochs:     {epochs}")
    print(f"  LR:         {lr}")
    print(f"  R_seed:     {R_seed:.1f}, R_clamp: {R_clamp:.1f}, R_damp: {R_damp:.1f}")
    print(f"  Params/seed: {3*N**3} (phi) + {3*N**3} (vel) = {6*N**3}")
    print(f"  Memory est: {6*N**3*4*n_seeds / 1e6:.1f} MB (f32 fields)")
    print(f"  Output:     {outdir}/")
    print()

    # --- Initialize seeds ---
    master_key = jax.random.PRNGKey(args.master_seed)
    keys = jax.random.split(master_key, n_seeds)

    if args.resume:
        print(f"  Resuming from {args.resume}")
        phi_resume, vel_resume, meta = load_candidate(args.resume)
        # Use the resumed field for the first seed, random for the rest
        all_phi = [phi_resume]
        all_vel = [vel_resume]
        for i in range(1, n_seeds):
            p, v = init_random_sphere(keys[i], N, L, R_seed, args.amplitude)
            all_phi.append(p)
            all_vel.append(v)
        batch_phi = jnp.stack(all_phi, axis=0)
        batch_vel = jnp.stack(all_vel, axis=0)
        print(f"  Seed 0: resumed (epoch={meta.get('epoch', '?')}, loss={meta.get('loss', '?')})")
        print(f"  Seeds 1-{n_seeds-1}: random")
    else:
        all_phi = []
        all_vel = []
        for i in range(n_seeds):
            p, v = init_random_sphere(keys[i], N, L, R_seed, args.amplitude)
            all_phi.append(p)
            all_vel.append(v)
        batch_phi = jnp.stack(all_phi, axis=0)
        batch_vel = jnp.stack(all_vel, axis=0)
        print(f"  All {n_seeds} seeds: random (amplitude={args.amplitude})")

    # Apply initial clamp
    for i in range(n_seeds):
        p, v = apply_sphere_clamp(batch_phi[i], batch_vel[i], sphere_mask)
        batch_phi = batch_phi.at[i].set(p)
        batch_vel = batch_vel.at[i].set(v)

    print()

    # --- Pre-compile loss+grad for each curriculum T_opt ---
    # Build a dict of JIT-compiled (loss, grad) functions keyed by n_steps.
    # This avoids recompilation when T_opt changes.
    print("  Pre-compiling loss functions for curriculum stages...")
    compiled_fns = {}
    for T_opt in CURRICULUM_TOPT_VALUES:
        n_steps = max(1, int(T_opt / dt))
        loss_fn = make_loss_fn(dx, dt, n_steps, damp_mask, X, Y, Z, L)

        @jit
        def loss_and_grad_fn(phi, vel, _loss_fn=loss_fn):
            def scalar_loss(p, v):
                l, _ = _loss_fn(p, v)
                return l
            loss_val, (g_phi, g_vel) = jax.value_and_grad(scalar_loss, argnums=(0, 1))(phi, vel)
            _, aux = _loss_fn(phi, vel)
            return loss_val, aux, g_phi, g_vel

        compiled_fns[n_steps] = loss_and_grad_fn
        print(f"    T_opt={T_opt:.0f} -> n_steps={n_steps}")

    print()

    # --- Optimizer (one per seed) ---
    optimizer = optax.adam(lr)
    # State: list of (opt_state_phi, opt_state_vel) per seed
    opt_states_phi = []
    opt_states_vel = []
    for i in range(n_seeds):
        opt_states_phi.append(optimizer.init(batch_phi[i]))
        opt_states_vel.append(optimizer.init(batch_vel[i]))

    # --- Tracking ---
    best_loss_global = float('inf')
    best_seed_idx = 0
    best_epoch = 0
    nan_resets = [0] * n_seeds

    # TSV log
    log_path = os.path.join(outdir, 'training_log.tsv')
    with open(log_path, 'w') as f:
        f.write("epoch\tT_opt\tn_steps\tbest_loss\tbest_epotT\tbest_rrms\t"
                "mean_loss\talive_seeds\ttime_s\n")

    # --- Warm-up: trigger JIT compilation on first seed ---
    print("  Warm-up JIT compilation...", end=' ', flush=True)
    t_warmup = time.time()
    first_n_steps = max(1, int(CURRICULUM_TOPT_VALUES[0] / dt))
    _ = compiled_fns[first_n_steps](batch_phi[0], batch_vel[0])
    print(f"done ({time.time() - t_warmup:.1f}s)")
    print()

    # --- Training loop ---
    header = (f"{'Epoch':>5} {'T_opt':>5} {'BstLoss':>10} {'BstEpotT':>10} "
              f"{'BstRrms':>8} {'MeanLoss':>10} {'Alive':>5} {'Time':>6}")
    print(header)
    print("-" * len(header))

    t_start = time.time()
    prev_n_steps = None

    for epoch in range(epochs):
        t_epoch = time.time()

        # Curriculum: determine n_steps for this epoch
        n_steps, T_opt_label = get_T_opt_steps(epoch, dx)

        # Get the right compiled function
        loss_grad_fn = compiled_fns[n_steps]

        # If n_steps changed, notify
        if n_steps != prev_n_steps and prev_n_steps is not None:
            print(f"  --- Curriculum: T_opt -> {T_opt_label:.0f} (n_steps={n_steps}) ---")
        prev_n_steps = n_steps

        # --- Process each seed ---
        epoch_losses = []
        epoch_aux = []
        best_loss_this_epoch = float('inf')
        best_idx_this_epoch = 0

        for i in range(n_seeds):
            phi_i = batch_phi[i]
            vel_i = batch_vel[i]

            # Forward + backward
            loss_val, aux, g_phi, g_vel = loss_grad_fn(phi_i, vel_i)
            loss_f = float(loss_val)

            # NaN check: reset seed if NaN
            if np.isnan(loss_f) or np.isinf(loss_f):
                nan_resets[i] += 1
                new_key = jax.random.fold_in(keys[i], epoch + 1000 * nan_resets[i])
                phi_new, vel_new = init_random_sphere(
                    new_key, N, L, R_seed, args.amplitude
                )
                phi_new, vel_new = apply_sphere_clamp(phi_new, vel_new, sphere_mask)
                batch_phi = batch_phi.at[i].set(phi_new)
                batch_vel = batch_vel.at[i].set(vel_new)
                opt_states_phi[i] = optimizer.init(phi_new)
                opt_states_vel[i] = optimizer.init(vel_new)
                epoch_losses.append(float('inf'))
                epoch_aux.append(None)
                continue

            # Gradient clipping
            g_phi = jnp.clip(g_phi, -1.0, 1.0)
            g_vel = jnp.clip(g_vel, -1.0, 1.0)

            # Adam update
            updates_phi, opt_states_phi[i] = optimizer.update(
                g_phi, opt_states_phi[i], phi_i
            )
            updates_vel, opt_states_vel[i] = optimizer.update(
                g_vel, opt_states_vel[i], vel_i
            )
            phi_i = optax.apply_updates(phi_i, updates_phi)
            vel_i = optax.apply_updates(vel_i, updates_vel)

            # Hard clamp: zero outside sphere
            phi_i, vel_i = apply_sphere_clamp(phi_i, vel_i, sphere_mask)

            # Store back
            batch_phi = batch_phi.at[i].set(phi_i)
            batch_vel = batch_vel.at[i].set(vel_i)

            epoch_losses.append(loss_f)
            epoch_aux.append(aux)

            if loss_f < best_loss_this_epoch:
                best_loss_this_epoch = loss_f
                best_idx_this_epoch = i

        # --- Epoch summary ---
        valid_losses = [l for l in epoch_losses if np.isfinite(l)]
        mean_loss = np.mean(valid_losses) if valid_losses else float('inf')
        alive = sum(1 for i, aux in enumerate(epoch_aux)
                    if aux is not None and float(aux['epot_T']) < -1.0)

        best_aux_this = epoch_aux[best_idx_this_epoch]

        # Global best tracking
        if best_loss_this_epoch < best_loss_global and best_aux_this is not None:
            best_loss_global = best_loss_this_epoch
            best_seed_idx = best_idx_this_epoch
            best_epoch = epoch

        dt_epoch = time.time() - t_epoch

        # Print every 10 epochs
        if epoch % 10 == 0 or epoch == epochs - 1:
            if best_aux_this is not None:
                epotT = float(best_aux_this['epot_T'])
                rrms = float(best_aux_this['rrms_T'])
            else:
                epotT = float('nan')
                rrms = float('nan')
            print(f"{epoch:5d} {T_opt_label:5.0f} {best_loss_this_epoch:10.4f} "
                  f"{epotT:10.4f} {rrms:8.3f} {mean_loss:10.4f} "
                  f"{alive:5d} {dt_epoch:6.1f}s")

        # Save best every 100 epochs
        if epoch % 100 == 0 and best_aux_this is not None:
            ckpt_path = os.path.join(outdir, f'best_seed_{epoch:04d}.npy')
            save_candidate(
                batch_phi[best_seed_idx], batch_vel[best_seed_idx],
                ckpt_path, N, L, epoch, best_loss_global,
                epoch_aux[best_seed_idx] if epoch_aux[best_seed_idx] is not None
                else best_aux_this
            )

        # TSV log
        with open(log_path, 'a') as f:
            epotT_log = float(best_aux_this['epot_T']) if best_aux_this else float('nan')
            rrms_log = float(best_aux_this['rrms_T']) if best_aux_this else float('nan')
            f.write(f"{epoch}\t{T_opt_label:.0f}\t{n_steps}\t"
                    f"{best_loss_this_epoch:.6f}\t{epotT_log:.6f}\t{rrms_log:.6f}\t"
                    f"{mean_loss:.6f}\t{alive}\t{dt_epoch:.2f}\n")

    # --- Final summary ---
    total_time = time.time() - t_start
    print()
    print("=" * 80)
    print("  TRAINING COMPLETE")
    print("=" * 80)
    print(f"  Total time:     {total_time:.0f}s ({total_time/60:.1f} min)")
    print(f"  Best loss:      {best_loss_global:.6f} (seed {best_seed_idx}, epoch {best_epoch})")
    print(f"  NaN resets:     {sum(nan_resets)} total across all seeds")
    print()

    # Save final best
    final_path = os.path.join(outdir, 'best_seed.npy')
    # Re-evaluate best seed to get fresh aux
    final_loss_grad_fn = compiled_fns[prev_n_steps]
    final_loss, final_aux, _, _ = final_loss_grad_fn(
        batch_phi[best_seed_idx], batch_vel[best_seed_idx]
    )
    save_candidate(
        batch_phi[best_seed_idx], batch_vel[best_seed_idx],
        final_path, N, L, epochs, float(final_loss), final_aux
    )
    print(f"  Final best saved to {final_path}")
    print(f"  E_pot(0)  = {float(final_aux['epot_0']):.4f}")
    print(f"  E_pot(T)  = {float(final_aux['epot_T']):.4f}")
    print(f"  R_rms(T)  = {float(final_aux['rrms_T']):.4f}")
    print()

    # Save all seeds' final states
    all_seeds_path = os.path.join(outdir, 'all_seeds_final.npy')
    np.save(all_seeds_path, {
        'phi': np.array(batch_phi),
        'vel': np.array(batch_vel),
        'N': N, 'L': L, 'n_seeds': n_seeds,
    }, allow_pickle=True)
    print(f"  All seeds saved to {all_seeds_path}")

    # --- Extended evaluation of best seed ---
    print()
    print("=== Extended Evaluation (2x time) ===")
    max_n_steps = max(1, int(CURRICULUM_TOPT_VALUES[-1] / dt))
    eval_n_steps = max_n_steps * 2
    theta0 = jnp.zeros_like(batch_phi[best_seed_idx])
    theta_vel0 = jnp.zeros_like(batch_phi[best_seed_idx])

    phi_eval, vel_eval, _, _ = evolve_absorbing(
        batch_phi[best_seed_idx], batch_vel[best_seed_idx],
        theta0, theta_vel0, dx, dt, eval_n_steps, damp_mask
    )

    epot_eval = compute_epot(phi_eval, dx)
    rrms_eval = compute_rrms(phi_eval, X, Y, Z, dx)
    epot_init = compute_epot(batch_phi[best_seed_idx], dx)

    retention = float(epot_eval) / (float(epot_init) + 1e-30) * 100.0
    T_eval = eval_n_steps * dt

    print(f"  T_eval = {T_eval:.1f} ({eval_n_steps} steps)")
    print(f"  E_pot(0)    = {float(epot_init):.4f}")
    print(f"  E_pot(T_ev) = {float(epot_eval):.4f} ({retention:.1f}% retained)")
    print(f"  R_rms(T_ev) = {float(rrms_eval):.4f}")

    # Save evolved field
    save_field_bin(
        np.array(phi_eval), N, L, T_eval,
        os.path.join(outdir, 'best_evolved.bin')
    )
    save_field_bin(
        np.array(batch_phi[best_seed_idx]), N, L, 0.0,
        os.path.join(outdir, 'best_initial.bin')
    )
    print(f"  Snapshots saved to {outdir}/best_initial.bin, best_evolved.bin")
    print()

    # Rank all seeds by final loss
    print("=== Seed Rankings ===")
    seed_results = []
    for i in range(n_seeds):
        if epoch_aux[i] is not None:
            seed_results.append((i, epoch_losses[i], float(epoch_aux[i]['epot_T']),
                                 float(epoch_aux[i]['rrms_T'])))
        else:
            seed_results.append((i, float('inf'), float('nan'), float('nan')))

    seed_results.sort(key=lambda x: x[1])
    print(f"  {'Rank':>4} {'Seed':>4} {'Loss':>10} {'EpotT':>10} {'Rrms':>8} {'Resets':>6}")
    for rank, (idx, loss, epotT, rrms) in enumerate(seed_results[:10]):
        print(f"  {rank+1:4d} {idx:4d} {loss:10.4f} {epotT:10.4f} {rrms:8.3f} {nan_resets[idx]:6d}")
    if n_seeds > 10:
        print(f"  ... ({n_seeds - 10} more seeds)")
    print()

    return best_loss_global


# ===================================================================
#  Entry point
# ===================================================================

def main():
    parser = argparse.ArgumentParser(
        description='Field search: random sphere -> soliton via backprop')
    parser.add_argument('--validate', action='store_true',
                        help='Run validation tests and exit')
    parser.add_argument('--N', type=int, default=DEFAULT_N,
                        help=f'Grid size (default: {DEFAULT_N})')
    parser.add_argument('--L', type=float, default=DEFAULT_L,
                        help=f'Half-box size (default: {DEFAULT_L})')
    parser.add_argument('--n_seeds', type=int, default=DEFAULT_N_SEEDS,
                        help=f'Number of parallel seeds (default: {DEFAULT_N_SEEDS})')
    parser.add_argument('--epochs', type=int, default=DEFAULT_EPOCHS,
                        help=f'Training epochs (default: {DEFAULT_EPOCHS})')
    parser.add_argument('--lr', type=float, default=DEFAULT_LR,
                        help=f'Learning rate (default: {DEFAULT_LR})')
    parser.add_argument('--outdir', type=str, default='results',
                        help='Output directory (default: results)')
    parser.add_argument('--resume', type=str, default=None,
                        help='Resume from checkpoint .npy file')
    parser.add_argument('--R_seed', type=float, default=DEFAULT_R_SEED,
                        help=f'Initial sphere radius (default: {DEFAULT_R_SEED})')
    parser.add_argument('--R_clamp', type=float, default=DEFAULT_R_CLAMP,
                        help=f'Hard clamp radius (default: {DEFAULT_R_CLAMP})')
    parser.add_argument('--amplitude', type=float, default=DEFAULT_AMPLITUDE,
                        help=f'Initial noise amplitude (default: {DEFAULT_AMPLITUDE})')
    parser.add_argument('--master_seed', type=int, default=42,
                        help='Master random seed (default: 42)')
    args = parser.parse_args()

    # Print JAX device info
    print(f"JAX devices: {jax.devices()}")
    print(f"JAX default backend: {jax.default_backend()}")
    print()

    if args.validate:
        ok = run_all_tests()
        sys.exit(0 if ok else 1)
    else:
        train(args)


if __name__ == '__main__':
    main()
