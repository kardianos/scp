# V42: Deuterium — Two Phase-Confined Baryons

## Goal

Bind a UUD (proton) and UDD (neutron) baryon at N=512, L=100, T=500.
Each baryon is a phase-confined 3-braid composite from V41.
The inter-baryon binding comes from the residual depletion interaction.

## Setup

- Grid: N=512, L=100 (dx=0.39, similar to V41's effective resolution)
- GPU: Tesla V100-32GB (19.3 GB for 18 f64 arrays at N=512)
- Two baryons separated by ~40 code units (center to center)
- Each baryon: 3 phase-confined braids (Δ = {0, 2π/3, 4π/3})
- Proton (UUD) at (-20, 0, 0), Neutron (UDD) at (+20, 0, 0)
- T=500, snap_dt=100, f32 output (~6 frames × 6.4 GB ≈ 38 GB)
- Absorbing BC: damp_width=10, damp_rate=0.005

## Seed Generation

The seed generator creates both baryons in one field:
1. Background field at A_bg=0.1
2. UUD proton at cx=-20: three braids xyz with chirality UUD, phases {0, 2π/3, 4π/3}
3. UDD neutron at cx=+20: three braids xyz with chirality UDD, phases {0, 2π/3, 4π/3}
4. Pre-loaded θ at 80% equilibrium
5. Contracting velocity profile per baryon

## SFA Size Management

- Simulation output: f32, snap_dt=100 → ~6 frames → ~38 GB
- Download: convert to f16 with every-other-frame → ~3 frames × 3.2 GB = ~10 GB
- Or: convert to f16 all frames → 6 × 3.2 GB = ~19 GB
- Disk requirement: 50+ GB on instance
