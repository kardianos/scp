import numpy as np

# Optimize creases to minimize U, constrain det(R)>0
def compute_R(creases):
    R = np.zeros((3,3))
    for i in range(3):
        for j in range(3):
            if i != j:
                R[i,j] = np.dot(creases[i], creases[j])  # Fix: use dot for scalar
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