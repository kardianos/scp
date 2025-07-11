import numpy as np
from scipy.fft import fftn, ifftn, fftfreq
import struct
import time

# Parameters (toy model)
N = 32  # Grid size (XYZ)
T_steps = 100  # Time steps
dx = 10.0 / N
dt = 0.01
k = 1.0  # Restoring (mass-like)
lam = 0.1  # Nonlinear (repulsion)
gamma = 0.01  # Damping for unstable (set >0 for decays)
T = 1.0  # Temperature (noise scale)

# Chiral factors: Array for field, e.g., +1/-1 for chiral, 0 for achiral
chiral = np.zeros((N, N, N))  # Default achiral (neutrino-like)
# Example: Proton-like (locked +1 in region)
chiral[10:15, 10:15, 10:15] = 1.0
# Electron-like (dual sub-modes, say -1 with flexibility)
chiral[20:25, 20:25, 20:25] = -1.0

# Laplacian operator (precompute)
k_vec = 2 * np.pi * fftfreq(N, d=dx)
kx, ky, kz = np.meshgrid(k_vec, k_vec, k_vec)
lap_op = -(kx**2 + ky**2 + kz**2)

# Initial field (standing waves as Gaussians)
phi = np.exp(-((np.indices((N,N,N))[0]-N/4)**2 + (np.indices((N,N,N))[1]-N/4)**2 + (np.indices((N,N,N))[2]-N/4)**2) / (N/8)**2)
vel = np.zeros_like(phi)

# File setup
file_path = 'field.scpv'
with open(file_path, 'wb') as f:
    f.write(b'SCPV\x00\x00\x00\x00')  # Magic with NULL
    # Initial tokens
    def write_token(type_id, name, value_bytes):
        name_bytes = name.encode('utf-8')
        f.write(struct.pack('>q h', type_id, len(name_bytes)) + name_bytes + struct.pack('>q', len(value_bytes)) + value_bytes)
    
    write_token(0, 'height', struct.pack('>i', N))
    write_token(1, 'width', struct.pack('>i', N))
    write_token(2, 'length', struct.pack('>i', N))
    write_token(3, 'encoding', b'float64')  # Field data type

    buffer = bytearray()
    flush_interval = 10
    step = 0

    while step < T_steps:
        # Add temperature noise
        noise = np.random.normal(0, T, phi.shape)
        phi += noise
        
        # "Cool" by emitting null rotors (propagating waves ∝ T^4)
        emission_amp = T**4 / 1000  # Scaled for sim
        emission = emission_amp * np.sin(2 * np.pi * (np.indices((N,N,N))[0]/N - step*dt))  # Simple x-dir wave
        phi -= emission
        
        # Step evolution
        accel = np.real(ifftn(lap_op * fftn(phi))) - k * phi - lam * phi**3 * chiral - gamma * vel  # Chiral in nonlinear, damping
        new_vel = vel + 0.5 * dt * accel
        new_phi = phi + dt * new_vel
        new_accel = np.real(ifftn(lap_op * fftn(new_phi))) - k * new_phi - lam * new_phi**3 * chiral - gamma * new_vel
        new_vel += 0.5 * dt * new_accel
        
        phi, vel = new_phi, new_vel
        
        # Buffer field data token
        field_bytes = phi.astype(np.float64).tobytes()
        name_bytes = f't{step}'.encode('utf-8')
        buffer.extend(struct.pack('>q h', 4, len(name_bytes)) + name_bytes + struct.pack('>q', len(field_bytes)) + field_bytes)
        
        step += 1
        if step % flush_interval == 0:
            f.write(buffer)
            buffer.clear()
            f.flush()
    
    if buffer:
        f.write(buffer)
        f.flush()