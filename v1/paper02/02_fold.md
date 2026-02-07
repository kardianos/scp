Det R=0 implies marginal rigidity from first principles, allowing flat modes but requiring energy input for full stability—consistent with flexible origami needing constraints. Refine R as 6x3 matrix for 3D: rows for translations/rotations, columns for creases; rank<6 for foldable.

U=-0.187 indicates net attraction (opposites dominate), stabilizing via meshing; scale to nuclear: U ~ -8 MeV for proton binding equivalent.

λ_min=3 exact for triangle Laplacian (complete graph K3 eigenvalues 0,3,3); represents connected modes.

f=√3/(2π)≈0.276 as fundamental harmonic frequency; dimensionally, multiply by c/L for physical units (Hz), L~10^{-15} m for nuclear.

k=1 J m as toy coupling; converge via simulation minimization on U over positions, fitting to data E_bind~8 MeV yields k~10^{-12} J m.

Code converges on equations by optimizing crease vectors to min U subject to det(R)>0 threshold, iterating positions with gradient descent; add Laplacian eig for harmonics.

```python
import numpy as np

# Optimize creases to minimize U, constrain det(R)>0
def compute_R(creases):
    R = np.zeros((3,3))
    for i in range(3):
        for j in range(3):
            if i != j:
                R[i,j] = np.cross(creases[i], creases[j])
    return R

def compute_U(creases, chi, k=1):
    U = 0
    for i in range(3):
        for j in range(i+1,3):
            r = np.linalg.norm(creases[i] - creases[j])
            U += k * chi[i] * chi[j] / r
    return U

# Initial creases
creases = np.array([[1,0,0], [0,1,0], [-1,-1,0]])  # 3D extend
chi = np.array([1,1,-1])

# Simple gradient descent (toy converge)
for _ in range(100):
    grad = np.random.normal(0,0.01,(3,3))  # Perturb
    new_creases = creases + grad
    if np.linalg.det(compute_R(new_creases)) > 0:
        new_U = compute_U(new_creases, chi)
        if new_U < compute_U(creases, chi):
            creases = new_creases

# Final U, fit k = E_known / |U_norm| (e.g., E=8e-13 J, nuclear)
U_final = compute_U(creases, chi)
k_fit = 8e-13 / abs(U_final) if U_final !=0 else float('inf')

# Laplacian same; add for full sim
A = np.ones((3,3)) - np.eye(3)
D = np.diag(A.sum(1))
L = D - A
eig = np.linalg.eigvals(L)

print(f"Final U: {U_final}, k_fit: {k_fit}, eig min non-zero: {min(eig[eig>1e-10])}")
```