# Kinetic inertia sketch: \(\tfrac12(E_\star/c^2)v^2\) — v0

**Approach:** A  
**Round:** 2  
**Status:** partial derivation; coefficient 1 open  
**Parent:** `free_response_kernel_v0.md` §6

---

## 1. Claim (killable)

For a stable lock with unlock ledger \(E_\star\), free-medium resistance to slow acceleration satisfies

\[
E_{\mathrm{kin}}
\;=\;
\frac{1}{2}\frac{E_\star}{c^2}\,v^2
\;+\;O(v^4/c^2)
\]

in the free rest frame of the ambient medium, with \(c\) free locality.

Equivalent: \(m_{\mathrm{inertial}}=E_\star/c^2\).

---

## 2. Free-response energy (F1)

Quasistatic free capacity potential:

\[
-\sigma_0\nabla^2\psi = s\rho_b,
\quad
U[\psi]=\frac{\sigma_0}{2}\int|\nabla\psi|^2\,dV=\frac{s}{2}\int\psi\rho_b\,dV.
\]

Static self-energy \(U(0)\) depends on lock form factor (UV). Subtract it for inertia.

Moving lock \(\rho_b(x-vt)\): free medium cannot rearrange faster than \(c\). In frequency domain the free operator is the free wave/relaxational operator; expansion for small \(v/c\) yields

\[
U(v)=U(0)+\frac{1}{2}\mu v^2+\cdots.
\]

**EM analogy (structural, not dualist photons):** electromagnetic mass \(\mu_{\mathrm{em}}\sim U_{\mathrm{field}}/c^2\). Here \(U_{\mathrm{field}}\to U[\psi]\) of free continuum.

---

## 3. Scaling argument (why \(E_\star/c^2\))

1. Linear response: \(\psi\propto s E_\star/(\sigma_0 R)\) near characteristic lock size \(R\).  
2. \(U(0)\propto s^2 E_\star^2/(\sigma_0 R)\) (self-energy).  
3. Retardation / boost corrections replace \(R\) by kinematic factors involving \(c\); the kinetic piece scales as \(U(0)\cdot(v/c)^2\) times \(O(1)\) geometric factor \(\xi\):

\[
\frac{1}{2}\mu v^2 \sim \xi\,\frac{U(0)}{c^2}\,v^2.
\]

4. **Problem:** \(\mu\sim U(0)/c^2\) tracks **self-energy**, not necessarily unlock ledger \(E_\star\).  
   Classical EM has the same tension (electromagnetic mass vs bare mass).

### 3.1 Two ways to force \(\mu=E_\star/c^2\)

| Route | Idea | Status |
|-------|------|--------|
| **R1 Universality** | Free constitutive normalization chosen so that the only free energy scale available to inertia is \(E_\star\) (bound ledger), not independent \(U(0)\). Self-energy renormalized into \(E_\star\). | Needs precise free action + renormalization rule |
| **R2 Drag definition** | Define \(m\) from low-frequency drag of free medium on accelerating lock; show equality to unlock work \(/c^2\) by virtual work | Best for B numeric |

Round 2 adopts **R2 as operational** for congruence T5; **R1** remains theoretical open.

---

## 4. Virtual-work sketch (R2)

Displace lock center by \(\delta x\) quasistatically:

- Unlock ledger fixed if Stable(L) holds (no unlock).  
- Free state \(\psi\) must shift with the lock: \(\delta U = \int (\delta\psi/\delta x_{\mathrm{L}})\,(\cdots)\).  
- External free-medium force \(F_{\mathrm{ext}}=-\partial U/\partial x_{\mathrm{L}}\).

For steady velocity (no radiation in relaxational model), power \(F_{\mathrm{ext}}v\) balances free dissipation if any; for pure potential free energy (F1 quasistatic elliptic), force derives from comoving \(U\).

Accelerating: effective Lagrangian \(L=\tfrac12\mu\dot x^2 - U_{\mathrm{ext}}\).  
**Missing step:** prove \(\mu=E_\star/c^2\) from free wave speed \(c\) entering the free kinetic term of \(\psi\):

\[
\mathcal{L}_\psi \supset \frac{\sigma_0}{2c^2}(\partial_t\psi)^2 - \frac{\sigma_0}{2}|\nabla\psi|^2
\]

(hyperbolic free capacity). Then standard field-momentum arguments give inertia \(\propto\) free field energy \(/c^2\). Identifying free field energy cost of the lock with \(E_\star\) is the monist step (budget conversion).

---

## 5. What B should measure

```text
1. Form lock; measure E_star = ∫ ρ_b
2. Apply weak body force or prescribed acceleration a of lock center
   while free ψ relaxes
3. m_inertial = F_ext / a  (fit low-a linear regime)
4. Compare m_inertial vs E_star/c^2 vs M_ray from exterior ℓ
```

**Pass:** all three agree within \(\varepsilon_{\mathrm{triad}}\).  
**Kill:** stable disagreement under free_relaxation origin.

---

## 6. Round-3 status

- No coefficient closure this round (priority = 3D free multipole congruence).  
- B R2 correctly **deferred** tautological \(a=F/m_L\) “inertia.”  
- 3D F1 self-energy more singular (\(\propto 1/R_{\mathrm{lock}}\)) than 2D log — renormalization into \(E_\star\) more important before claiming \(\mu=E_\star/c^2\).  
- Next useful step: free-drag force from **\(\psi\)-wake** at prescribed velocity (not prescribed \(m\)).

## 7. Bottom line

- Kinetic inertia is **free-response cost of moving a lock**, not primitive mass.  
- Scaling \(\mu\sim E/c^2\) is natural; **exact** \(E=E_\star\) needs renormalization or virtual-work closure.  
- Congruence treats triad equality as **T5 killable test**, not yet a theorem.
