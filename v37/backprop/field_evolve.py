"""
Evolutionary field search: Game-of-Life style.

Instead of backprop, use selection pressure:
1. Generate many random compact field configs
2. Evolve each for T_eval with absorbing BC (forward only, no grad)
3. Measure fitness = how much binding survived
4. Keep top survivors, mutate, repeat

Like training an LLM with REINFORCE instead of supervised learning.

Usage:
    python field_evolve.py --pop 256 --N 32 --L 15 --generations 500 --outdir results/
"""

import jax
import jax.numpy as jnp
from jax import jit, vmap
import numpy as np
import argparse
import time
import os
import struct

# Physics
MU = -41.345
KAPPA = 50.0
MASS2 = 2.25
ETA = 0.5
MTHETA2 = 0.0
DELTA = jnp.array([0.0, 3.0005, 4.4325])


def make_grid(N, L):
    x = jnp.linspace(-L, L, N)
    X, Y, Z = jnp.meshgrid(x, x, x, indexing='ij')
    R = jnp.sqrt(X**2 + Y**2 + Z**2)
    dx = 2*L/(N-1)
    return X, Y, Z, R, dx


def make_damp_mask(R, L, R_damp, rate=0.01):
    raw = jnp.where(R > R_damp, 1.0 - rate*((R - R_damp)/(L - R_damp))**2, 1.0)
    return jnp.clip(raw, 0.0, 1.0)


def laplacian(f, idx2):
    """Laplacian with zero-pad (Dirichlet) BC."""
    padded = jnp.pad(f, 1, mode='constant', constant_values=0.0)
    return (padded[2:,1:-1,1:-1] + padded[:-2,1:-1,1:-1] +
            padded[1:-1,2:,1:-1] + padded[1:-1,:-2,1:-1] +
            padded[1:-1,1:-1,2:] + padded[1:-1,1:-1,:-2] -
            6.0 * f) * idx2


def curl_a(F, a, idx1):
    """Curl component a with zero-pad BC."""
    def d(f, axis):
        p = jnp.pad(f, 1, mode='constant', constant_values=0.0)
        slc_p = [slice(1,-1)]*3; slc_m = [slice(1,-1)]*3
        slc_p[axis] = slice(2, None); slc_m[axis] = slice(None, -2)
        return (p[tuple(slc_p)] - p[tuple(slc_m)]) * idx1
    if a == 0: return d(F[2], 1) - d(F[1], 2)
    if a == 1: return d(F[0], 2) - d(F[2], 0)
    return d(F[1], 0) - d(F[0], 1)


def compute_forces(phi, theta, dx):
    idx2 = 1.0/(dx*dx)
    idx1 = 1.0/(2.0*dx)
    P = phi[0]*phi[1]*phi[2]
    den = 1.0 + KAPPA*P*P
    mPd2 = MU*P/(den*den)

    phi_acc = jnp.zeros_like(phi)
    theta_acc = jnp.zeros_like(theta)
    for a in range(3):
        lap = laplacian(phi[a], idx2)
        dPda = phi[(a+1)%3]*phi[(a+2)%3]
        ct = curl_a(theta, a, idx1)
        phi_acc = phi_acc.at[a].set(lap - MASS2*phi[a] - mPd2*dPda + ETA*ct)

    for a in range(3):
        lap_t = laplacian(theta[a], idx2)
        cp = curl_a(phi, a, idx1)
        theta_acc = theta_acc.at[a].set(lap_t - MTHETA2*theta[a] + ETA*cp)

    return phi_acc, theta_acc


def step(phi, pv, theta, tv, dx, dt, damp):
    pa, ta = compute_forces(phi, theta, dx)
    pv = pv + 0.5*dt*pa; tv = tv + 0.5*dt*ta
    phi = phi + dt*pv; theta = theta + dt*tv
    pa, ta = compute_forces(phi, theta, dx)
    pv = pv + 0.5*dt*pa; tv = tv + 0.5*dt*ta
    pv = pv * damp[None,:,:,:]; tv = tv * damp[None,:,:,:]
    return phi, pv, theta, tv


def evolve(phi0, pv0, dx, dt, n_steps, damp):
    theta0 = jnp.zeros_like(phi0)
    tv0 = jnp.zeros_like(phi0)
    def body(carry, _):
        p, pv, t, tv = carry
        p, pv, t, tv = step(p, pv, t, tv, dx, dt, damp)
        return (p, pv, t, tv), None
    (phi_T, pv_T, theta_T, tv_T), _ = jax.lax.scan(body, (phi0, pv0, theta0, tv0), None, length=n_steps)
    return phi_T, pv_T


def compute_epot(phi, dx):
    dV = dx**3
    P = phi[0]*phi[1]*phi[2]
    V = (MU/2.0)*P*P/(1.0+KAPPA*P*P)
    return jnp.sum(V)*dV


def compute_pint(phi, dx):
    dV = dx**3
    return jnp.sum(jnp.abs(phi[0]*phi[1]*phi[2]))*dV


def compute_rrms(phi, X, Y, Z, dx):
    dV = dx**3
    rho = jnp.sum(phi**2, axis=0)
    total = jnp.sum(rho)*dV + 1e-30
    r2 = X**2+Y**2+Z**2
    return jnp.sqrt(jnp.sum(r2*rho)*dV/total)


def fitness(phi_T, dx, X, Y, Z, L):
    """Higher is better. Combines binding + compactness."""
    ep = compute_epot(phi_T, dx)  # negative is good
    rrms = compute_rrms(phi_T, X, Y, Z, dx)
    pint = compute_pint(phi_T, dx)
    # Fitness: maximize -ep (binding), penalize spread
    return -ep - 0.01*(rrms/L)**2, ep, pint, rrms


def random_sphere(key, N, L, R_seed, amplitude, R):
    """Correlated random: base field * phase offsets + perturbation."""
    k1, k2 = jax.random.split(key)

    # Smooth random base
    base = jax.random.normal(k1, (N, N, N)) * amplitude
    for _ in range(5):
        base = (base + jnp.roll(base,1,0)+jnp.roll(base,-1,0)
                +jnp.roll(base,1,1)+jnp.roll(base,-1,1)
                +jnp.roll(base,1,2)+jnp.roll(base,-1,2)) / 7.0

    # Phase-correlated fields
    phi = jnp.stack([
        base * jnp.cos(DELTA[0]),
        base * jnp.cos(DELTA[1]),
        base * jnp.cos(DELTA[2]),
    ], axis=0)

    # Small independent perturbation
    pert = jax.random.normal(k2, (3, N, N, N)) * amplitude * 0.15
    for _ in range(3):
        smoothed = []
        for a in range(3):
            s = (pert[a]+jnp.roll(pert[a],1,0)+jnp.roll(pert[a],-1,0)
                 +jnp.roll(pert[a],1,1)+jnp.roll(pert[a],-1,1)
                 +jnp.roll(pert[a],1,2)+jnp.roll(pert[a],-1,2))/7.0
            smoothed.append(s)
        pert = jnp.stack(smoothed, axis=0)
    phi = phi + pert

    # Envelope
    sigma = R_seed * 0.6
    env = jnp.where(R < R_seed, jnp.exp(-R**2/(2*sigma**2)), 0.0)
    phi = phi * env[None,:,:,:]

    vel = jnp.zeros_like(phi)
    return phi, vel


def mutate(key, phi, vel, R, R_clamp, mutation_rate=0.1):
    """Add smooth random perturbation to a candidate."""
    k1, k2 = jax.random.split(key)
    N = phi.shape[1]

    # Smooth perturbation
    pert = jax.random.normal(k1, phi.shape) * mutation_rate
    for _ in range(3):
        smoothed = []
        for a in range(3):
            s = (pert[a]+jnp.roll(pert[a],1,0)+jnp.roll(pert[a],-1,0)
                 +jnp.roll(pert[a],1,1)+jnp.roll(pert[a],-1,1)
                 +jnp.roll(pert[a],1,2)+jnp.roll(pert[a],-1,2))/7.0
            smoothed.append(s)
        pert = jnp.stack(smoothed, axis=0)

    # Clamp to sphere
    mask = jnp.where(R < R_clamp, 1.0, 0.0)
    phi_new = (phi + pert) * mask[None,:,:,:]

    # Small velocity perturbation too
    vpert = jax.random.normal(k2, vel.shape) * mutation_rate * 0.3
    for _ in range(2):
        smoothed = []
        for a in range(3):
            s = (vpert[a]+jnp.roll(vpert[a],1,0)+jnp.roll(vpert[a],-1,0)
                 +jnp.roll(vpert[a],1,1)+jnp.roll(vpert[a],-1,1)
                 +jnp.roll(vpert[a],1,2)+jnp.roll(vpert[a],-1,2))/7.0
            smoothed.append(s)
        vpert = jnp.stack(smoothed, axis=0)
    vel_new = (vel + vpert) * mask[None,:,:,:]

    return phi_new, vel_new


def crossover(key, phi_a, vel_a, phi_b, vel_b):
    """Blend two parents."""
    alpha = jax.random.uniform(key, minval=0.3, maxval=0.7)
    return alpha*phi_a + (1-alpha)*phi_b, alpha*vel_a + (1-alpha)*vel_b


def save_candidate(phi, vel, path, N, L):
    os.makedirs(os.path.dirname(path) or '.', exist_ok=True)
    np.savez_compressed(path, phi=np.array(phi), vel=np.array(vel), N=N, L=L)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--pop', type=int, default=128, help='Population size')
    parser.add_argument('--N', type=int, default=32)
    parser.add_argument('--L', type=float, default=15.0)
    parser.add_argument('--T_eval', type=float, default=10.0, help='Eval time per generation')
    parser.add_argument('--generations', type=int, default=500)
    parser.add_argument('--R_seed', type=float, default=5.0)
    parser.add_argument('--R_clamp', type=float, default=6.0)
    parser.add_argument('--amplitude', type=float, default=0.5)
    parser.add_argument('--mutation_rate', type=float, default=0.1)
    parser.add_argument('--elite_frac', type=float, default=0.1)
    parser.add_argument('--outdir', type=str, default='results_evolve')
    parser.add_argument('--seed', type=int, default=42)
    args = parser.parse_args()

    N, L = args.N, args.L
    X, Y, Z, R, dx = make_grid(N, L)
    dt = 0.1 * dx
    n_steps = max(1, int(args.T_eval / dt))
    R_damp = L - 4.0
    damp = make_damp_mask(R, L, R_damp)
    n_elite = max(2, int(args.pop * args.elite_frac))

    print("=" * 70)
    print("  EVOLUTIONARY FIELD SEARCH")
    print("=" * 70)
    print(f"  N={N} L={L} dx={dx:.4f} dt={dt:.5f}")
    print(f"  Population: {args.pop}, Elite: {n_elite}")
    print(f"  T_eval={args.T_eval} ({n_steps} steps)")
    print(f"  Generations: {args.generations}")
    print(f"  R_seed={args.R_seed} R_clamp={args.R_clamp} R_damp={R_damp:.1f}")
    print(f"  Mutation rate: {args.mutation_rate}")
    print()

    os.makedirs(args.outdir, exist_ok=True)

    # JIT compile evolve + fitness
    print("  JIT compiling...", end=' ', flush=True)
    t0 = time.time()

    @jit
    def eval_one(phi, vel):
        phi_T, _ = evolve(phi, vel, dx, dt, n_steps, damp)
        fit, ep, pint, rrms = fitness(phi_T, dx, X, Y, Z, L)
        return fit, ep, pint, rrms

    # Warmup
    key = jax.random.PRNGKey(0)
    test_phi, test_vel = random_sphere(key, N, L, args.R_seed, args.amplitude, R)
    _ = eval_one(test_phi, test_vel)
    print(f"done ({time.time()-t0:.1f}s)")
    print()

    # Initialize population
    print("  Initializing population...", end=' ', flush=True)
    master_key = jax.random.PRNGKey(args.seed)
    keys = jax.random.split(master_key, args.pop)

    pop_phi = []
    pop_vel = []
    for i in range(args.pop):
        p, v = random_sphere(keys[i], N, L, args.R_seed, args.amplitude, R)
        pop_phi.append(p)
        pop_vel.append(v)
    print("done")
    print()

    # Log
    log_path = os.path.join(args.outdir, 'evolution_log.tsv')
    with open(log_path, 'w') as f:
        f.write("gen\tbest_fit\tbest_epot\tbest_pint\tbest_rrms\tmean_fit\talive\ttime_s\n")

    header = f"{'Gen':>4} {'BstFit':>10} {'BstEpot':>10} {'BstPint':>8} {'BstRrms':>7} {'MeanFit':>10} {'Alive':>5} {'Time':>6}"
    print(header)
    print("-" * len(header))

    best_ever_fit = -1e30
    best_ever_phi = None
    best_ever_vel = None

    for gen in range(args.generations):
        t_gen = time.time()

        # Evaluate all
        fits = []; epots = []; pints = []; rrmss = []
        for i in range(len(pop_phi)):
            f, ep, pi, rr = eval_one(pop_phi[i], pop_vel[i])
            fits.append(float(f))
            epots.append(float(ep))
            pints.append(float(pi))
            rrmss.append(float(rr))

        fits = np.array(fits)
        epots = np.array(epots)
        pints = np.array(pints)
        rrmss = np.array(rrmss)

        # Sort by fitness (descending)
        order = np.argsort(-fits)
        best_idx = order[0]
        best_fit = fits[best_idx]
        alive = np.sum(epots < -0.5)

        dt_gen = time.time() - t_gen

        if gen % 10 == 0 or gen == args.generations - 1:
            print(f"{gen:4d} {best_fit:10.4f} {epots[best_idx]:10.4f} "
                  f"{pints[best_idx]:8.2f} {rrmss[best_idx]:7.3f} "
                  f"{np.mean(fits):10.4f} {alive:5d} {dt_gen:6.1f}s")

        with open(log_path, 'a') as f:
            f.write(f"{gen}\t{best_fit:.6f}\t{epots[best_idx]:.6f}\t"
                    f"{pints[best_idx]:.4f}\t{rrmss[best_idx]:.4f}\t"
                    f"{np.mean(fits):.6f}\t{alive}\t{dt_gen:.1f}\n")

        # Save best
        if best_fit > best_ever_fit:
            best_ever_fit = best_fit
            best_ever_phi = pop_phi[best_idx]
            best_ever_vel = pop_vel[best_idx]
            save_candidate(best_ever_phi, best_ever_vel,
                          os.path.join(args.outdir, f'best_gen{gen:04d}.npz'), N, L)

        # Checkpoint every 50 gens
        if gen % 50 == 0:
            save_candidate(pop_phi[best_idx], pop_vel[best_idx],
                          os.path.join(args.outdir, f'checkpoint_gen{gen:04d}.npz'), N, L)

        # Selection + reproduction
        elite_phi = [pop_phi[order[i]] for i in range(n_elite)]
        elite_vel = [pop_vel[order[i]] for i in range(n_elite)]

        new_phi = list(elite_phi)  # keep elites unchanged
        new_vel = list(elite_vel)

        gen_key = jax.random.PRNGKey(args.seed * 10000 + gen)
        sub_keys = jax.random.split(gen_key, args.pop)

        for i in range(n_elite, args.pop):
            k1, k2, k3 = jax.random.split(sub_keys[i], 3)

            # Pick parent(s) from elite
            p1_idx = int(jax.random.randint(k1, (), 0, n_elite))

            if jax.random.uniform(k2) < 0.3 and n_elite > 1:
                # Crossover
                p2_idx = int(jax.random.randint(k2, (), 0, n_elite))
                child_phi, child_vel = crossover(k3, elite_phi[p1_idx], elite_vel[p1_idx],
                                                  elite_phi[p2_idx], elite_vel[p2_idx])
            else:
                child_phi = elite_phi[p1_idx]
                child_vel = elite_vel[p1_idx]

            # Mutate
            # Adaptive mutation: higher rate for early gens, lower later
            mr = args.mutation_rate * max(0.1, 1.0 - gen / args.generations)
            child_phi, child_vel = mutate(k3, child_phi, child_vel, R, args.R_clamp, mr)
            new_phi.append(child_phi)
            new_vel.append(child_vel)

        pop_phi = new_phi
        pop_vel = new_vel

    # Final summary
    print()
    print("=" * 70)
    print("  EVOLUTION COMPLETE")
    print("=" * 70)
    print(f"  Best fitness: {best_ever_fit:.6f}")
    f, ep, pi, rr = eval_one(best_ever_phi, best_ever_vel)
    print(f"  E_pot(T): {float(ep):.4f}")
    print(f"  P_int(T): {float(pi):.4f}")
    print(f"  R_rms(T): {float(rr):.4f}")

    save_candidate(best_ever_phi, best_ever_vel,
                  os.path.join(args.outdir, 'best_final.npz'), N, L)
    print(f"  Saved to {args.outdir}/best_final.npz")


if __name__ == '__main__':
    main()
