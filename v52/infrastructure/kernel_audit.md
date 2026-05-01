# CUDA Kernel Audit: sfa/sim/scp_sim.cu vs Theory

## Summary

The CUDA kernel force terms are **algebraically correct** -- they match the
theoretical equations derived in Maxima and formalized in Lean. The 99.99% energy
loss is caused by the **absorbing boundary layer**, which is active by default
even when the user expects periodic boundary conditions.

---

## 1. Force-by-Force Verification

### 1.1 Phi Forces (kernel lines 684-711)

The kernel computes:
```
acc_phi[a] = lap - m^2*phi[a] - dVdP*dPda - t2c*phi[a]
           + eta*curl(theta)_a + chiral
           - alpha_cs * curl(M)_a       // Cosserat strain
           - 2.0 * curl(Q)_a            // Curl^2 hardening
```

**Theory (V50/C4 two-pass formulation):**
```
F_phi_a = Lap(phi)_a - m^2*phi_a - V'(P)*dP/dphi_a
        + eta*curl(theta)_a
        - alpha*curl(M)_a
        - 2*curl(Q)_a
```

| Term | Kernel | Theory | Match |
|------|--------|--------|-------|
| Laplacian | `(phi[a][n_ip]+...−6*phi[a])*idx2` (L685-686) | standard 7-point stencil | YES |
| Mass | `-me2*phi[a][idx]` (L707) | `-m^2*phi_a` | YES |
| Potential | `-dVdP*dPda` (L707), dVdP = `mu*P/(den*den)` (L631) | `-V'(P)*dP/dphi_a` where V'(P) = mu*P/(1+kappa*P^2)^2 | YES |
| Curl coupling | `+d_ETA*ct` (L708), ct = curl(theta)_a (L690-693) | `+eta*curl(theta)_a` | YES |
| Cosserat strain | `-d_ALPHA_CS*curl_M[a]` (L709) | `-alpha*curl(M)_a` | YES |
| Hardening | `-2.0*curl_Q[a]` (L710) | `-2*curl(Q)_a` | YES |

The Cosserat strain phi force is proven correct in `lean/V50C3/Cosserat.lean`:
- `strainMomentum` (L61): EL momentum pi_{ja} = alpha * Sigma_b M_b * eps_{bja}
- `levi_civita_sign_identity` (L134-148): proves eps_{bja} = -eps_{ajb} (all 9 cases)
- Result: F_phi = -alpha*curl(M), with the **minus sign** verified by exhaustive case analysis.

The hardening phi force is proven structurally in `lean/V50C4/CurlHardening.lean`:
- `hardenQ` (L25): Q_a = (beta/2)|theta|^2 * curl(phi)_a
- Phi force = -2*curl(Q) (follows from Euler-Lagrange with curl_smul identity)

The algebraic equivalence of two-pass and single-pass formulations is proven in
`lean/SinglePass.lean` (L1-93) and verified numerically in `v50/em_wave/single_pass.mac`
to 10^{-18} precision.

### 1.2 Theta Forces (kernel lines 714-733)

The kernel computes:
```
acc_theta[a] = lapt - mtheta2*theta[a] + eta*curl(phi)_a
             + 2*alpha_cs*M[a]            // Cosserat restoring
             - beta_h*|curl(phi)|^2*theta[a]  // Hardening
```

**Theory:**
```
F_theta_a = Lap(theta)_a - m_theta^2*theta_a + eta*curl(phi)_a
          + 2*alpha*M_a
          - beta*|curl(phi)|^2*theta_a
```

| Term | Kernel | Theory | Match |
|------|--------|--------|-------|
| Laplacian | `(th[a][n_ip]+...−6*th[a])*idx2` (L715-716) | standard 7-point stencil | YES |
| Mass | `-d_MTHETA2*th[a][idx]` (L732) | `-m_theta^2*theta_a` | YES |
| Curl coupling | `+d_ETA*cp` (L732), cp = curl(phi)_a (L717-720) | `+eta*curl(phi)_a` | YES |
| Cosserat | `+2.0*d_ALPHA_CS*M[a][idx]` (L726) | `+2*alpha*M_a` | YES |
| Hardening | `-d_BETA_H*curl_sq*th[a][idx]` (L730) | `-beta*|curl(phi)|^2*theta_a` | YES |

Theta Cosserat force proven in `lean/V50C3/Cosserat.lean`:
- `thetaForceStrain` (L33): F_theta = 2*alpha*M
- `thetaForce_is_neg_deriv` (L44): this equals -dE_strain/dtheta_a
- `strain_force_linear` (L223-231): at theta = curl(phi)/2 + delta_theta, force = -2*alpha*delta_theta (restoring)

Hardening theta force proven in `lean/V50C4/CurlHardening.lean`:
- `hardenThetaForce` (L21): F_theta = -beta*|curl(phi)|^2*theta_a

### 1.3 Intermediate Fields (kernel lines 562-599)

**Pass 1 (compute_intermediates_kernel):**
- Mismatch: `M_a = curl(phi)_a/2 - theta_a` (L586-588) -- matches `lean/V50C3/Cosserat.lean` `mismatch` definition
- Hardening Q: `Q_a = (beta/2)*|theta|^2*curl(phi)_a` (L596-598) -- matches `lean/V50C4/CurlHardening.lean` `hardenQ` definition

**Pass 2 (compute_forces_kernel):** Takes curls of M and Q using the standard 2-point central difference stencil, consistent with how curl(phi) and curl(theta) are computed.

### 1.4 Potential V(P) and its Derivative

Kernel (L629-631):
```c
double P = p0*p1*p2, P2 = P*P;
double den = 1.0 + keff*P2;
double dVdP = d_MU * P / (den*den);
```

Theory: V(P) = (mu/2)*P^2/(1+kappa*P^2)
- V'(P) = mu*P/(1+kappa*P^2)^2  -- matches kernel exactly

---

## 2. Energy Conservation Analysis

### 2.1 Lagrangian Derivability

Every force term in the kernel is derivable from a well-defined Lagrangian/energy:

| Force | Energy density | Lagrangian? |
|-------|---------------|-------------|
| Lap(phi) | (1/2)\|grad(phi)\|^2 | YES (kinetic gradient) |
| -m^2*phi | (1/2)m^2\|phi\|^2 | YES (mass) |
| -V'(P)*dP/dphi | V(P) | YES (potential) |
| +eta*curl(theta) | -eta*phi.curl(theta) | YES (curl coupling) |
| -alpha*curl(M) | alpha\|M\|^2 | YES (Cosserat strain, proven in Lean) |
| -2*curl(Q) | (beta/2)\|theta\|^2\|curl(phi)\|^2 | YES (hardening, proven in Lean) |
| Lap(theta) | (1/2)\|grad(theta)\|^2 | YES |
| -m_theta^2*theta | (1/2)m_theta^2\|theta\|^2 | YES |
| +eta*curl(phi) | -eta*phi.curl(theta) | YES (same coupling term) |
| +2*alpha*M | alpha\|M\|^2 | YES (same strain energy) |
| -beta\|curl(phi)\|^2*theta | (beta/2)\|theta\|^2\|curl(phi)\|^2 | YES (same hardening energy) |

**All forces are Hamiltonian.** The system has a conserved energy functional:
```
E = E_phi_kin + E_theta_kin + E_grad + E_mass + E_pot + E_tgrad + E_tmass + E_coupling
  + E_strain + E_harden [+ E_chiral if kappa_h != 0]
```

The Cosserat force is NOT dissipative: the `+2*alpha*M` on theta is a **restoring spring**
(toward theta = curl(phi)/2), not damping. It oscillates, exchanging energy between phi
and theta sectors. This is proven in `lean/V50C3/Cosserat.lean` theorem `theta_restoring_sign` (L70-76).

### 2.2 Velocity Verlet / Leapfrog

The kernel uses standard velocity Verlet (lines 1128-1168):
```
half-kick:   v += (dt/2) * a
drift:       x += dt * v
forces:      a = F(x)        [with intermediates M, Q recomputed]
half-kick:   v += (dt/2) * a
boundary:    damp v in absorbing shell
```

Velocity Verlet is symplectic for Hamiltonian systems and conserves energy to O(dt^2)
with no secular drift. **This is correct** for the force terms above.

The one non-symplectic element is the **absorbing boundary** (line 1165-1167),
which multiplies velocities by a damping factor < 1 near the boundary.

### 2.3 Diagnostics Energy Computation

The GPU reduction kernel (`reduce_diagnostics_kernel`, L878-1026) computes:

```
E_total = E_phi_kin + E_theta_kin + E_grad + E_mass + E_pot + E_tgrad + E_tmass + E_coupling
```

Where E_pot includes (L945-961):
- V(P) potential
- alpha_cs * |M|^2 (Cosserat strain energy)
- (beta_h/2) * |theta|^2 * |curl(phi)|^2 (hardening energy)
- kappa_h * P^2 * phi.curl(phi) (chiral energy)

And E_coupling = -eta * phi . curl(theta) (L968).

**The energy computation includes ALL terms.** There is no missing energy contribution.

---

## 3. ROOT CAUSE: Absorbing Boundary Layer

### 3.1 The Problem

The default configuration (`bc_type=0`) enables absorbing boundary damping:
```c
c.damp_width = 3.0;  c.damp_rate = 0.01;   // LINE 64
```

The absorbing_boundary_kernel (L766-787) damps velocities in a spherical shell:
```c
if (r > R_damp) {
    double s = (r - R_damp) / d_DAMP_WIDTH;
    double damp = 1.0 - d_DAMP_RATE * s * s;
    vp0[idx] *= damp; vp1[idx] *= damp; ...   // ALL 6 velocity fields
    vt0[idx] *= damp; vt1[idx] *= damp; ...
}
```

where `R_damp = L - damp_width`. This runs EVERY timestep (L1165).

With periodic indexing, the grid wraps around. Waves that leave one face enter from
the opposite face -- they hit the absorbing shell every time they cross the boundary.
Over long simulations, this removes virtually all kinetic energy.

### 3.2 Why 99.99% Loss

For a grid of size 2L with damp_width=3.0, the fraction of the grid volume inside
the absorbing shell is approximately:
```
V_shell/V_total = 1 - (1 - 3/L)^3   (for spherical shell in cube)
```
For L=50: V_shell/V_total ~ 17%. For L=15: V_shell/V_total ~ 47%.

Any wave energy crossing this region loses a fraction `damp_rate * s^2` of its
kinetic energy PER TIMESTEP. Over thousands of timesteps, this compounds to
near-total energy loss.

### 3.3 How to Get True Periodic BC

Set in the config file:
```
damp_width=0
damp_rate=0
```

This causes the absorbing_boundary_kernel to return immediately (L772):
```c
if (d_DAMP_WIDTH <= 0 || d_DAMP_RATE <= 0) return;
```

Alternatively, use `bc_switch_time=T` to switch from absorbing to periodic at time T
(L1839-1844), which sets d_DAMP_RATE=0 on the GPU.

### 3.4 The Grid IS Periodic

The indexing uses modular arithmetic everywhere:
```c
int ip=(i+1)%N, im=(i-1+N)%N, jp=(j+1)%N, jm=(j-1+N)%N, kp=(k+1)%N, km=(k-1+N)%N;
```

The Laplacian, curl, and all spatial derivatives correctly wrap around. The grid
topology is T^3 (3-torus). The absorbing boundary is a LAYER ON TOP of this
periodic topology that damps velocities.

---

## 4. Secondary Findings

### 4.1 No Asymmetries in Force Terms

The phi-theta coupling is symmetric:
- Phi sees `+eta*curl(theta)` from energy `-eta*phi.curl(theta)`
- Theta sees `+eta*curl(phi)` from the same energy term
- Integration by parts: phi.curl(theta) = theta.curl(phi) + boundary terms (zero for periodic BC)

The Cosserat strain is symmetric:
- E_strain = alpha*|M|^2 = alpha*|curl(phi)/2 - theta|^2
- Both phi and theta forces are exact variational derivatives (proven in Lean)

The hardening is symmetric:
- E_harden = (beta/2)|theta|^2|curl(phi)|^2
- Theta force: -beta*|curl(phi)|^2*theta (correct algebraic derivative)
- Phi force: -2*curl(Q) where Q = (beta/2)|theta|^2*curl(phi) (correct EL derivative, proven in Lean via curl_smul identity)

### 4.2 The Two-Pass Algorithm is Exact

The two-pass implementation (compute M and Q first, then take their curls) is
algebraically identical to the single-pass expansion. This is:
- Proven in `lean/SinglePass.lean` using curl_curl and curl_smul axioms
- Verified numerically in `v50/em_wave/single_pass.mac` to 10^{-18} precision

The two-pass approach avoids computing grad(div(phi)) and grad(|theta|^2) x curl(phi)
explicitly, which would require additional derivative computations. Instead, it
computes the same result using only curls of vector fields.

### 4.3 Mode 3 Term II (L634-638)

The `t2c` term handles the density-dependent kappa mode:
```c
if (d_MODE == 3 && d_kappa_gamma > 0) {
    double D = 1.0 + d_kappa_gamma*sig + d_KAPPA*P2;
    t2c = d_MU * d_kappa_gamma * d_KAPPA * P2 * P2 / (D*D);
}
```
This adds `-t2c*phi[a]` to the phi force (L707). This is the derivative of the
modified potential V(P,sigma) = (mu/2)*P^2*(1+kappa_gamma*sigma)/(1+kappa_gamma*sigma+kappa*P^2)
with respect to sigma, times d(sigma)/d(phi_a) = 2*phi_a. The factor of 2 is absorbed
into the expression. This is only active for mode=3 and does not affect the standard
Cosserat equations.

---

## 5. Recommendations

1. **To test energy conservation with periodic BC**: Set `damp_width=0` and `damp_rate=0`
   in the config file. Energy drift should then be O(dt^2) per orbit, without secular loss.

2. **Default behavior warning**: The kernel's default `bc_type=0` with `damp_width=3.0`
   and `damp_rate=0.01` is **intentionally dissipative** -- it's designed for compact
   objects that emit radiation. This is correct for proton formation/relaxation but
   misleading for energy conservation tests.

3. **No kernel modifications needed**: The force terms are correct and energy-conserving.
   The dissipation is entirely from the absorbing boundary layer, which is a feature
   (not a bug) for its intended use case.

---

## Files Referenced

- **CUDA kernel**: `/home/d/code/scp/sfa/sim/scp_sim.cu`
- **Lean proofs**: `lean/V50C3/Cosserat.lean`, `lean/V50C4/CurlHardening.lean`, `lean/SinglePass.lean`, `lean/V44/Equations.lean`, `lean/ScpLib/Energy.lean`
- **Maxima derivations**: `v50/em_wave/single_pass.mac`, `v52/core_interaction.mac`, `v52/equilibrium_seed.mac`, `v52/verify_cancellation.mac`, `v52/theta_equilibrium.mac`
- **Force summary in Lean**: `lean/V50C3/Forces.lean` (complete C3 equations, L111-133)
