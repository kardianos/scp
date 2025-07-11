import numpy as np
from numpy.fft import fftn, ifftn, fftfreq

def get_laplacian_op(N, dx):
    k = 2 * np.pi * fftfreq(N, d=dx)
    kx, ky, kz = np.meshgrid(k, k, k)
    return -(kx**2 + ky**2 + kz**2)

def step(phi, vel, dt, lap_op, k_param, lam):
    accel = np.real(ifftn(lap_op * fftn(phi))) - k_param * phi - lam * phi**3
    new_vel = vel + 0.5 * dt * accel
    new_phi = phi + dt * new_vel
    new_accel = np.real(ifftn(lap_op * fftn(new_phi))) - k_param * new_phi - lam * new_phi**3
    new_vel = new_vel + 0.5 * dt * new_accel
    return new_phi, new_vel

def print_slice(phi, t, N):
    print(f"Time: {t:.2f}")
    slice_ = phi[:, :, N//2]  # XY slice at Z center
    minv = np.min(slice_)
    maxv = np.max(slice_)
    if maxv > minv:
        norm = (slice_ - minv) / (maxv - minv)
    else:
        norm = np.zeros_like(slice_)
    chars = ' .:-=+*#%@'
    for i in range(N):
        line = ''
        for j in range(N):
            idx = int(norm[i, j] * (len(chars) - 1))
            line += chars[idx] * 2  # Wider characters for better visibility
        print(line)
    print("\n")

N = 32
L = 10.0
dx = L / N
dt = 0.01
num_steps = 2000
print_every = 20
k_param = 1.0  # Linear term coefficient (related to mass/harmonic)
lam = 0.1      # Nonlinear coefficient for repulsion
A = 1.0        # Amplitude

x = np.linspace(0, L, N, endpoint=False)
X, Y, Z = np.meshgrid(x, x, x, indexing='ij')

# Initial standing wave for a 'particle'
phi = A * np.sin(2 * np.pi * X / L) * np.sin(2 * np.pi * Y / L) * np.sin(2 * np.pi * Z / L)
vel = np.zeros_like(phi)

lap_op = get_laplacian_op(N, dx)

t = 0.0
print_slice(phi, t, N)

emitted = False

for step_num in range(num_steps):
    if not emitted and t >= 1.0:
        phi *= 0.5  # Simulate drop to lower energy state
        # Add light emission as a propagating massless wave (ignore nonlinear for it)
        # Simple plane wave in x direction
        phi += 0.5 * np.sin(2 * np.pi * (X / L - (t - 1.0))) * np.exp(-((Y - L/2)**2 + (Z - L/2)**2) / (L/4)**2)
        emitted = True

    phi, vel = step(phi, vel, dt, lap_op, k_param, lam)
    t += dt
    if (step_num + 1) % print_every == 0:
        print_slice(phi, t, N)
