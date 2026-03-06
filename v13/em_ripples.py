import sympy as sp
import os
import json

print("HFKT v13: Avenue 1 - Electromagnetic Ripples from Bulk Topology")
print("===============================================================")
print("Projecting Cl+(3,0,1) rotor dynamics onto the Faraday Bivector F")

os.makedirs("v13/results", exist_ok=True)

# Define space-time coordinates
t, x, y, z = sp.symbols('t x y z', real=True)
r = sp.sqrt(x**2 + y**2 + z**2)

# Define a time-harmonic, spatially decaying localized knot (Rotor)
omega = sp.Symbol('omega', real=True)
k = sp.Symbol('k', real=True) 

# Simple oscillating bump form: theta(r, t) = A * exp(-r) * cos(omega*t - k*r)
A = sp.Symbol('A', real=True)
theta = A * sp.exp(-r) * sp.cos(omega*t - k*r)

# Unit vector field for the knot profile
n_x = -y/r
n_y = x/r
n_z = 0 

q0 = sp.cos(theta/2)
q1 = n_x * sp.sin(theta/2)
q2 = n_y * sp.sin(theta/2)
q3 = n_z * sp.sin(theta/2)

R = [q0, q1, q2, q3]
R_tilde = [q0, -q1, -q2, -q3]

def quat_mult(a, b):
    w1, x1, y1, z1 = a
    w2, x2, y2, z2 = b
    return [
        w1*w2 - x1*x2 - y1*y2 - z1*z2,
        w1*x2 + x1*w2 + y1*z2 - z1*y2,
        w1*y2 - x1*z2 + y1*w2 + z1*x2,
        w1*z2 + x1*y2 - y1*x2 + z1*w2
    ]

def deriv(Q, var):
    return [sp.diff(q, var) for q in Q]

R_dot = deriv(R, t)

# Geometric electric field proxy: E ~ (d_t R) R~
E_quat = quat_mult(R_dot, R_tilde)
E_vec = E_quat[1:] 

dx_R = deriv(R, x)
dy_R = deriv(R, y)
dz_R = deriv(R, z)

# Geometric magnetic field proxy: B ~ curl{(nabla R) R~}
B_vec = [
    quat_mult(dy_R, R_tilde)[3] - quat_mult(dz_R, R_tilde)[2],
    quat_mult(dz_R, R_tilde)[1] - quat_mult(dx_R, R_tilde)[3],
    quat_mult(dx_R, R_tilde)[2] - quat_mult(dy_R, R_tilde)[1]
]

def div(V):
    return sp.diff(V[0], x) + sp.diff(V[1], y) + sp.diff(V[2], z)

div_B = div(B_vec)
div_E = div(E_vec)

# Evaluate specific points over a radial distance grid
eval_params = {A: 0.1, omega: 1, k: 1, t: 0}

out_tsv = "v13/results/em_ripple_radials.tsv"
print(f"\nStreaming intermediate radial evaluations to {out_tsv}...")

with open(out_tsv, "w") as f:
    f.write("Radius\tDiv_E\tDiv_B\tE_x\tB_x\n")
    
    # Evaluate radially from core to far-field
    for r_idx in range(1, 101):
        r_val = r_idx * 0.1  # 0.1 to 10.0
        
        # Test points along x-axis for simplicity (y=0, z=0)
        c_point = {**eval_params, x: r_val, y: 0.001, z: 0.001} # slight off-axis to avoid zero divisions
        
        div_B_val = float(div_B.subs(c_point).evalf())
        div_E_val = float(div_E.subs(c_point).evalf())
        ex_val = float(E_vec[0].subs(c_point).evalf())
        bx_val = float(B_vec[0].subs(c_point).evalf())
        
        f.write(f"{r_val:.1f}\t{div_E_val:.8f}\t{div_B_val:.8f}\t{ex_val:.8f}\t{bx_val:.8f}\n")

print("Streaming complete.")
