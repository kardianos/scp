# V25 Phase 5: Binary Oscillon Gravitational Wave Emission

## Thesis

V25 Phase 2 showed the spherically symmetric oscillon produces 98% monopole
(spin-0) radiation and only 0.5% quadrupole (spin-2). This is physically
correct: a spherical source cannot emit spin-2. Gravitational waves require
an ASPHERICAL, TIME-VARYING mass distribution.

The cleanest aspherical source: TWO oscillons orbiting each other. A binary
system has a time-varying quadrupole moment Q_{ij}(t) that radiates spin-2
gravitational waves at twice the orbital frequency.

In GR: P_GW = (G/5c⁵)⟨d³Q_{ij}/dt³⟩². The GW strain h ~ (2G/c⁴)(Ë_{ij}/r)
where Ë is the second time derivative of the quadrupole moment.

In our model: the TT strain σ^TT_{ij} from the elastic displacement field
plays the role of h. If two oscillons orbit, σ^TT should show quadrupolar
radiation at 2Ω (twice the orbital frequency).

## Setup

Two 3D oscillons placed at (0, 0, +D/2) and (0, 0, -D/2) with TANGENTIAL
velocities to create a circular orbit in the z=0 plane... actually, that
puts them on the z-axis with velocity in x or y. Let me set up properly:

**Configuration**: Two oscillons in the xy-plane:
- Oscillon A at (+D/2, 0, 0) with velocity (0, +v_orb, 0)
- Oscillon B at (-D/2, 0, 0) with velocity (0, -v_orb, 0)

This creates an orbit in the xy-plane. The orbital frequency Ω = v_orb/(D/2).
The quadrupole moment rotates at Ω, and GW are emitted at 2Ω.

**Orbital velocity**: For a bound orbit, v_orb² = F(D)·D/(2M) where F(D)
is the inter-oscillon force and M is the oscillon mass. From V24-S4
(3D merger at D=20): the oscillons attract. The orbital velocity needs to
BALANCE the attraction.

For a first test: don't try to find exact circular orbit. Instead, give
tangential velocity and let the system evolve. If the oscillons spiral
(energy loss from GW emission), that's even MORE physical — it's the
GW inspiral.

## Method

### Phase 5a: Establish the orbit

1. Use the V25 code (v25.c) with elastic + pairwise coupling
2. Grid: N=128, L=30 (need room for two oscillons + far field)
3. Equilibrate a SINGLE oscillon at N=96, L=15, t=500 (from Phase 1)
4. Place two copies at x=±D/2 with D=16 (close to V23-D equilibrium)
5. Give tangential velocities v_y = ±v_orb
6. Try v_orb = {0.01, 0.03, 0.05, 0.1}
7. Evolve t=500. Track positions. Does an orbit form?

### Phase 5b: Measure GW emission

8. At the best v_orb (stable or quasi-stable orbit): measure the TT
   strain on a sphere at R=20 from the center of mass
9. Decompose into angular multipoles (l=0,1,2,3,4)
10. KEY: is the l=2 (quadrupole) component dominant?
11. Measure the GW frequency: should be 2Ω
12. Measure h₊ and h× at the equator (θ=π/2) and poles (θ=0)
   - At equator: expect maximum h₊ and h×
   - At poles: expect ZERO (no GW along the orbital axis)
   This angular pattern IS the spin-2 signature

### Phase 5c: GW power and inspiral

13. Measure the total GW power: P_GW = ∫ (dE_strain/dt) r² dΩ
14. Does the orbital separation decrease? (inspiral from GW energy loss)
15. Does the orbital frequency increase? (chirp signal)
16. Compare P_GW with the quadrupole formula: P ∝ Ω⁶D⁴M²

## Expected GW Signature

For an equal-mass binary in the xy-plane:
- h₊(θ=π/2) oscillates at 2Ω (equatorial, maximum)
- h×(θ=π/2) oscillates at 2Ω (equatorial, maximum, 90° phase shifted)
- h₊(θ=0) = 0 (along orbital axis, no GW)
- h×(θ=0) = 0

This is EXACTLY the GR prediction for binary GW emission. If our model
reproduces this pattern, the TT strain IS spin-2 gravitational radiation.

## Practical Concerns

1. **Grid size**: Two oscillons at D=16 in a box of L=30 with far-field
   measurement at R=20. The oscillons are 10 units from the boundary.
   Absorbing BC should handle outgoing GW.

2. **Orbital stability**: Without exact force balance, the orbit won't
   be perfectly circular. Expect elliptical or spiral trajectory.
   This is fine — GW are emitted from ANY orbital motion.

3. **Resolution**: At N=128, L=30: dx=0.47. Oscillon core ~5 units =
   ~11 points. Marginal but workable (V24-S4 used N=96, L=25 successfully).

4. **Runtime**: Two oscillons at N=128 for t=500: ~15-30 min with OpenMP.

## Reference Code

- v25/src/v25.c (Phase 1-4 code, base for this)
- v24/proca_3d/src/proca_3d.c (3D two-oscillon tracking)
- v21/src/triad3d.c (3D solver)

## Output

- `src/binary.c` — binary oscillon GW simulator
- `data/binary_vorb{V}_ts.tsv` — position/separation vs time
- `data/binary_gw.tsv` — TT strain at 4+ angles vs time
- `data/binary_multipoles.tsv` — angular decomposition of GW
- `RESULTS.md`

## Parameters

μ=-20, κ=20, m=1.0, λ_pw=0.5, η=0.1, λ_L=0.1
α_g=0.001 (self-consistent metric, from Phase 3)
D=16, v_orb scan: {0.01, 0.03, 0.05, 0.1}
N=128, L=30, t=500

Compile: `gcc -O3 -fopenmp -Wall -o binary src/binary.c -lm`
