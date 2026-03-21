"""
Differentiable PDE-based structure search via backpropagation.

Treats the initial field configuration as learnable weights.
The PDE evolution is the forward pass. The loss measures binding survival.
JAX autodiff gives exact gradients via the adjoint method.

Usage:
    python structure_search.py [--N 32] [--T_opt 50] [--n_epochs 200] [--seed braid|random|knot]
"""

import jax
import jax.numpy as jnp
from jax import grad, jit, vmap
import optax
import numpy as np
import argparse
import time
import os
import struct

# Physics parameters
MU = -41.345
KAPPA = 50.0
MASS2 = 2.25
ETA = 0.5
A_BG = 0.1

# Phase offsets (from v28 optimization)
DELTA = jnp.array([0.0, 3.0005, 4.4325])


def make_grid(N, L):
    """Create 3D coordinate grid."""
    x = jnp.linspace(-L, L, N)
    X, Y, Z = jnp.meshgrid(x, x, x, indexing='ij')
    dx = 2 * L / (N - 1)
    return X, Y, Z, dx


def laplacian_3d(f, idx2):
    """6-point Laplacian with periodic BC."""
    return (jnp.roll(f, 1, 0) + jnp.roll(f, -1, 0) +
            jnp.roll(f, 1, 1) + jnp.roll(f, -1, 1) +
            jnp.roll(f, 1, 2) + jnp.roll(f, -1, 2) -
            6.0 * f) * idx2


def curl_component(F, a, idx1):
    """Curl component a of vector field F[3, N, N, N]."""
    if a == 0:
        return (jnp.roll(F[2], -1, 1) - jnp.roll(F[2], 1, 1) -
                jnp.roll(F[1], -1, 2) + jnp.roll(F[1], 1, 2)) * idx1
    elif a == 1:
        return (jnp.roll(F[0], -1, 2) - jnp.roll(F[0], 1, 2) -
                jnp.roll(F[2], -1, 0) + jnp.roll(F[2], 1, 0)) * idx1
    else:
        return (jnp.roll(F[1], -1, 0) - jnp.roll(F[1], 1, 0) -
                jnp.roll(F[0], -1, 1) + jnp.roll(F[0], 1, 1)) * idx1


def compute_forces(phi, theta, dx):
    """Compute accelerations for all 6 fields."""
    idx2 = 1.0 / (dx * dx)
    idx1 = 1.0 / (2.0 * dx)

    P = phi[0] * phi[1] * phi[2]
    den = 1.0 + KAPPA * P * P
    mPd2 = MU * P / (den * den)

    phi_acc = jnp.zeros_like(phi)
    theta_acc = jnp.zeros_like(theta)

    for a in range(3):
        lap = laplacian_3d(phi[a], idx2)
        if a == 0:
            dPda = phi[1] * phi[2]
        elif a == 1:
            dPda = phi[0] * phi[2]
        else:
            dPda = phi[0] * phi[1]

        curl_t = curl_component(theta, a, idx1)
        phi_acc = phi_acc.at[a].set(
            lap - MASS2 * phi[a] - mPd2 * dPda + ETA * curl_t
        )

    for a in range(3):
        lap_t = laplacian_3d(theta[a], idx2)
        curl_p = curl_component(phi, a, idx1)
        theta_acc = theta_acc.at[a].set(lap_t + ETA * curl_p)

    return phi_acc, theta_acc


def verlet_step(phi, phi_vel, theta, theta_vel, dx, dt):
    """One Verlet integration step."""
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

    return phi, phi_vel, theta, theta_vel


def evolve(phi0, phi_vel0, theta0, theta_vel0, dx, dt, n_steps):
    """Run n_steps of Verlet evolution. Differentiable via JAX."""
    def step_fn(carry, _):
        phi, phi_vel, theta, theta_vel = carry
        phi, phi_vel, theta, theta_vel = verlet_step(
            phi, phi_vel, theta, theta_vel, dx, dt
        )
        return (phi, phi_vel, theta, theta_vel), None

    (phi, phi_vel, theta, theta_vel), _ = jax.lax.scan(
        step_fn, (phi0, phi_vel0, theta0, theta_vel0), None, length=n_steps
    )
    return phi, phi_vel, theta, theta_vel


def compute_epot(phi, dx):
    """Compute total V(P) potential energy."""
    dV = dx ** 3
    P = phi[0] * phi[1] * phi[2]
    V = (MU / 2.0) * P * P / (1.0 + KAPPA * P * P)
    return jnp.sum(V) * dV


def compute_pint(phi, dx):
    """Integrated |P|."""
    dV = dx ** 3
    P = phi[0] * phi[1] * phi[2]
    return jnp.sum(jnp.abs(P)) * dV


def compute_rrms(phi, X, Y, Z, dx):
    """RMS radius of field energy."""
    dV = dx ** 3
    rho = phi[0]**2 + phi[1]**2 + phi[2]**2
    total = jnp.sum(rho) * dV + 1e-30
    r2 = X**2 + Y**2 + Z**2
    r2_avg = jnp.sum(r2 * rho) * dV / total
    return jnp.sqrt(r2_avg)


def compute_phase_coherence(phi, dx):
    """Phase coherence between field pairs."""
    dV = dx ** 3
    amp01 = jnp.sqrt(phi[0]**2 + phi[1]**2) + 1e-30
    cos01 = jnp.sum(phi[0] / amp01 * amp01) * dV
    sin01 = jnp.sum(phi[1] / amp01 * amp01) * dV
    total = jnp.sum(amp01) * dV + 1e-30
    coh = jnp.sqrt(cos01**2 + sin01**2) / total
    return coh


def loss_fn(phi0, phi_vel0, X, Y, Z, dx, dt, n_steps, L):
    """
    Loss function: maximize binding survival, minimize spreading.

    We evolve for n_steps, then measure:
    - E_pot retention (want high)
    - Compactness (want small R_rms)
    - P_int retention (want high)
    """
    theta0 = jnp.zeros_like(phi0)
    theta_vel0 = jnp.zeros_like(phi0)

    # Initial metrics
    epot0 = compute_epot(phi0, dx)
    pint0 = compute_pint(phi0, dx)

    # Evolve
    phi_T, _, _, _ = evolve(phi0, phi_vel0, theta0, theta_vel0, dx, dt, n_steps)

    # Final metrics
    epot_T = compute_epot(phi_T, dx)
    pint_T = compute_pint(phi_T, dx)
    rrms_T = compute_rrms(phi_T, X, Y, Z, dx)

    # Loss components (all to be MINIMIZED)
    # 1. Binding loss: want epot_T / epot0 close to 1
    #    epot is negative, so epot_T / epot0 > 0 means binding retained
    binding_loss = -epot_T / (jnp.abs(epot0) + 1e-10)  # negative = good retention

    # 2. Compactness: penalize large R_rms relative to box
    compact_loss = (rrms_T / L) ** 2

    # 3. P_int retention
    pint_loss = -pint_T / (pint0 + 1e-10)

    # 4. Smoothness: penalize high gradient energy in initial condition
    grad_energy = 0.0
    idx1 = 1.0 / (2.0 * dx)
    for a in range(3):
        gx = (jnp.roll(phi0[a], -1, 0) - jnp.roll(phi0[a], 1, 0)) * idx1
        gy = (jnp.roll(phi0[a], -1, 1) - jnp.roll(phi0[a], 1, 1)) * idx1
        gz = (jnp.roll(phi0[a], -1, 2) - jnp.roll(phi0[a], 1, 2)) * idx1
        grad_energy += jnp.sum(gx**2 + gy**2 + gz**2)
    smooth_loss = grad_energy * dx**3

    # Weighted combination
    total = (1.0 * binding_loss +
             0.3 * compact_loss +
             0.5 * pint_loss +
             0.001 * smooth_loss)

    return total, {
        'binding': binding_loss,
        'compact': compact_loss,
        'pint': pint_loss,
        'smooth': smooth_loss,
        'epot0': epot0,
        'epot_T': epot_T,
        'pint0': pint0,
        'pint_T': pint_T,
        'rrms_T': rrms_T,
    }


def seed_braid3(N, L, dx):
    """Seed with braid3(z) — the known-good structure."""
    X, Y, Z, _ = make_grid(N, L)
    R_h = 1.0
    r_tube = 2.0
    sigma_z = 3.0
    A = 0.8
    kw = jnp.pi / L
    omega = jnp.sqrt(kw**2 + MASS2)
    k_bg = jnp.pi / L
    omega_bg = jnp.sqrt(k_bg**2 + MASS2)

    phi = jnp.zeros((3, N, N, N))
    phi_vel = jnp.zeros((3, N, N, N))

    for a in range(3):
        val = jnp.zeros((N, N, N))
        vel = jnp.zeros((N, N, N))
        for s in range(3):
            cx = R_h * jnp.cos(kw * Z + 2 * jnp.pi * s / 3)
            cy = R_h * jnp.sin(kw * Z + 2 * jnp.pi * s / 3)
            d2 = (X - cx)**2 + (Y - cy)**2
            env = jnp.exp(-d2 / (2 * r_tube**2)) * jnp.exp(-Z**2 / (2 * sigma_z**2))
            phase = kw * Z + DELTA[a]
            val = val + A * env * jnp.cos(phase)
            vel = vel + omega * A * env * jnp.sin(phase)

        # Background
        ph_bg = k_bg * Z + 2 * jnp.pi * a / 3
        val = val + A_BG * jnp.cos(ph_bg)
        vel = vel + omega_bg * A_BG * jnp.sin(ph_bg)

        phi = phi.at[a].set(val)
        phi_vel = phi_vel.at[a].set(vel)

    return phi, phi_vel


def seed_random(N, L, dx, key):
    """Random seed with smooth spatial structure."""
    # Low-pass filtered random fields
    phi = jax.random.normal(key, (3, N, N, N)) * 0.3

    # Smooth by applying several rounds of averaging (poor man's low-pass)
    for _ in range(5):
        for a in range(3):
            phi = phi.at[a].set(
                (phi[a] +
                 jnp.roll(phi[a], 1, 0) + jnp.roll(phi[a], -1, 0) +
                 jnp.roll(phi[a], 1, 1) + jnp.roll(phi[a], -1, 1) +
                 jnp.roll(phi[a], 1, 2) + jnp.roll(phi[a], -1, 2)) / 7.0
            )

    # Add compact envelope centered at origin
    X, Y, Z, _ = make_grid(N, L)
    R = jnp.sqrt(X**2 + Y**2 + Z**2)
    envelope = jnp.exp(-R**2 / (2 * 4.0**2))  # sigma=4

    phi = phi * envelope[None, :, :, :]

    # Add background
    k_bg = jnp.pi / L
    for a in range(3):
        ph_bg = k_bg * Z + 2 * jnp.pi * a / 3
        phi = phi.at[a].set(phi[a] + A_BG * jnp.cos(ph_bg))

    phi_vel = jnp.zeros_like(phi)
    return phi, phi_vel


def save_field_bin(phi, phi_vel, N, L, t, path):
    """Save field snapshot for visualization."""
    os.makedirs(os.path.dirname(path) or '.', exist_ok=True)
    with open(path, 'wb') as f:
        f.write(struct.pack('i', N))
        f.write(struct.pack('d', L))
        f.write(struct.pack('d', t))
        f.write(struct.pack('i', 3))  # n_fields
        for a in range(3):
            f.write(np.array(phi[a], dtype=np.float32).tobytes())


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--N', type=int, default=32)
    parser.add_argument('--L', type=float, default=12.0)
    parser.add_argument('--T_opt', type=float, default=30.0,
                        help='Evolution time for loss evaluation')
    parser.add_argument('--n_epochs', type=int, default=300)
    parser.add_argument('--lr', type=float, default=0.001)
    parser.add_argument('--seed', type=str, default='braid',
                        choices=['braid', 'random'])
    parser.add_argument('--outdir', type=str, default='backprop_results')
    parser.add_argument('--device', type=str, default='gpu')
    args = parser.parse_args()

    N, L = args.N, args.L
    X, Y, Z, dx = make_grid(N, L)
    dt = 0.1 * dx
    n_steps = int(args.T_opt / dt)

    print(f"=== Structure Search via Backprop ===")
    print(f"N={N}, L={L}, dx={dx:.4f}, dt={dt:.5f}")
    print(f"T_opt={args.T_opt}, n_steps={n_steps}")
    print(f"Seed: {args.seed}, epochs: {args.n_epochs}, lr: {args.lr}")
    print(f"Parameters per field: {N**3}, total: {3*N**3}")
    print()

    # Initialize
    if args.seed == 'braid':
        phi0, phi_vel0 = seed_braid3(N, L, dx)
        print("Seeded from braid3(z)")
    else:
        key = jax.random.PRNGKey(42)
        phi0, phi_vel0 = seed_random(N, L, dx, key)
        print("Seeded from smoothed random")

    # Make phi0 the learnable parameter, phi_vel0 also learnable
    # Package as a single parameter vector
    params = {
        'phi': phi0,
        'vel': phi_vel0,
    }

    # JIT-compile the loss function
    @jit
    def loss_and_grad(params):
        loss_val_fn = lambda p: loss_fn(
            p['phi'], p['vel'], X, Y, Z, dx, dt, n_steps, L
        )
        (loss_val, aux), grads = jax.value_and_grad(loss_val_fn, has_aux=True)(params)
        return loss_val, aux, grads

    # Optimizer
    optimizer = optax.adam(args.lr)
    opt_state = optimizer.init(params)

    os.makedirs(args.outdir, exist_ok=True)

    # Training loop
    print(f"{'Epoch':>5} {'Loss':>10} {'Bind':>10} {'E_pot0':>10} {'E_potT':>10} "
          f"{'P_int0':>8} {'P_intT':>8} {'Rrms':>6} {'Time':>6}")
    print("-" * 90)

    best_loss = float('inf')
    t0 = time.time()

    for epoch in range(args.n_epochs):
        t_epoch = time.time()

        loss_val, aux, grads = loss_and_grad(params)

        # Clip gradients
        grads = jax.tree.map(lambda g: jnp.clip(g, -1.0, 1.0), grads)

        updates, opt_state = optimizer.update(grads, opt_state, params)
        params = optax.apply_updates(params, updates)

        dt_epoch = time.time() - t_epoch

        if epoch % 10 == 0 or epoch == args.n_epochs - 1:
            bind_ret = -float(aux['epot_T']) / (abs(float(aux['epot0'])) + 1e-10)
            print(f"{epoch:5d} {float(loss_val):10.4f} {bind_ret:10.4f} "
                  f"{float(aux['epot0']):10.2f} {float(aux['epot_T']):10.2f} "
                  f"{float(aux['pint0']):8.1f} {float(aux['pint_T']):8.1f} "
                  f"{float(aux['rrms_T']):6.2f} {dt_epoch:6.1f}s")

            # Save best
            if float(loss_val) < best_loss:
                best_loss = float(loss_val)
                np.save(os.path.join(args.outdir, 'best_phi.npy'),
                        np.array(params['phi']))
                np.save(os.path.join(args.outdir, 'best_vel.npy'),
                        np.array(params['vel']))

            # Save snapshot
            save_field_bin(np.array(params['phi']),
                          np.array(params['vel']),
                          N, L, float(epoch),
                          os.path.join(args.outdir, f'epoch_{epoch:04d}.bin'))

    total_time = time.time() - t0
    print(f"\nTotal time: {total_time:.0f}s ({total_time/60:.1f} min)")
    print(f"Best loss: {best_loss:.6f}")

    # Final evaluation at longer time
    print("\n=== Final evaluation (2x time) ===")
    n_eval = n_steps * 2
    theta0 = jnp.zeros_like(params['phi'])
    theta_vel0 = jnp.zeros_like(params['phi'])
    phi_T, _, _, _ = evolve(params['phi'], params['vel'],
                            theta0, theta_vel0, dx, dt, n_eval)
    epot_final = compute_epot(phi_T, dx)
    epot_init = compute_epot(params['phi'], dx)
    pint_final = compute_pint(phi_T, dx)
    rrms_final = compute_rrms(phi_T, X, Y, Z, dx)
    print(f"E_pot: {float(epot_init):.2f} -> {float(epot_final):.2f} "
          f"({100*float(epot_final)/float(epot_init):.1f}% retained)")
    print(f"P_int: {float(compute_pint(params['phi'], dx)):.1f} -> {float(pint_final):.1f}")
    print(f"R_rms: {float(rrms_final):.2f}")

    # Save final state
    save_field_bin(np.array(params['phi']), np.array(params['vel']),
                   N, L, 0.0,
                   os.path.join(args.outdir, 'final_init.bin'))
    save_field_bin(np.array(phi_T), np.array(jnp.zeros_like(phi_T)),
                   N, L, float(args.T_opt * 2),
                   os.path.join(args.outdir, 'final_evolved.bin'))

    print(f"\nResults saved to {args.outdir}/")


if __name__ == '__main__':
    main()
