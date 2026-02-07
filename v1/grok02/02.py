import numpy as np
from numpy.fft import fftn, ifftn, fftfreq
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from mpl_toolkits.mplot3d import Axes3D

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

# Parameters
N = 32  # Grid size (keep small for real-time performance)
L = 10.0  # Box length
dx = L / N
dt = 0.01  # Time step
k_param = 1.0  # Linear coefficient
lam = 0.1  # Nonlinear coefficient
A = 1.0  # Initial amplitude
threshold = 0.2  # Voxel display threshold (adjust for visibility)
frame_interval = 50  # ms between animation frames

# Grid for field (centers)
x = np.linspace(0, L, N, endpoint=False)
X, Y, Z = np.meshgrid(x, x, x, indexing='ij')

# Grid for voxel edges
x_edges = np.linspace(0, L, N + 1)
X_edges, Y_edges, Z_edges = np.meshgrid(x_edges, x_edges, x_edges, indexing='ij')

# Initial standing wave
phi = A * np.sin(2 * np.pi * X / L) * np.sin(2 * np.pi * Y / L) * np.sin(2 * np.pi * Z / L)
vel = np.zeros_like(phi)

lap_op = get_laplacian_op(N, dx)

# Time tracking
t = 0.0
emitted = False

# Set up figure and 3D axis
fig = plt.figure(figsize=(10, 10))
ax = fig.add_subplot(111, projection='3d')
ax.set_xlim(0, L)
ax.set_ylim(0, L)
ax.set_zlim(0, L)
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
ax.set_title('Volumetric Field Simulation')

# Initial voxel plot (empty placeholder)
initial_mask = np.zeros((N, N, N)) > 0
voxels = ax.voxels(X_edges, Y_edges, Z_edges, initial_mask, facecolors='blue', edgecolors='k', alpha=0.5)

def update(frame):
    global phi, vel, t, emitted, voxels
    
    # Simulate emission at t=1.0
    if not emitted and t >= 1.0:
        phi *= 0.5  # Drop to lower harmonic
        # Add propagating massless wave (plane wave in x, Gaussian in yz)
        phi += 0.5 * np.sin(2 * np.pi * (X / L - (t - 1.0))) * np.exp(-((Y - L/2)**2 + (Z - L/2)**2) / (L/4)**2)
        emitted = True
    
    # Evolve the field
    phi, vel = step(phi, vel, dt, lap_op, k_param, lam)
    t += dt
    
    # Clear previous voxels
    for voxel in list(voxels.values()):
        voxel.remove()
    
    # Compute voxel mask and colors
    mask = np.abs(phi) > threshold
    colors = np.empty(phi.shape + (3,), dtype=float)
    norm_phi = (phi - phi.min()) / (phi.max() - phi.min() + 1e-10)  # Normalize for coloring
    colors[..., 0] = norm_phi  # Red channel based on value
    colors[..., 1] = 1 - norm_phi  # Green inverse
    colors[..., 2] = 0.5  # Blue fixed
    
    # Plot new voxels
    voxels = ax.voxels(X_edges, Y_edges, Z_edges, mask, facecolors=colors, edgecolor='k', alpha=0.5)
    
    ax.set_title(f'Volumetric Field Simulation - Time: {t:.2f}')
    return list(voxels.values())

# Animate continuously (repeat=True for indefinite loop until window closed)
ani = FuncAnimation(fig, update, interval=frame_interval, blit=False, repeat=True)

plt.show()