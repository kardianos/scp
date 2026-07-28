# Fabric-as-particle media pack

Simulation captures + cinematic stills/animations from v72/v73 SFAs,
illustrating the process-form / fabric-as-particle idea (v73 PROCESS.md).

## Pipeline

1. **volview** (xvfb) — volume ray-march snapshots (field / charge / flavor / gauge)
2. **sfa_slice + render_slices.py** — scientific multi-quantity panels
3. **custom extract** — single-quantity cinematic panels (`scientific/*_rho2.png`, `*_s.png`)
4. **Imagine image_edit** — fabric interpretation of simulation geometry
5. **Imagine image_to_video** — 6s / 720p smooth animations

## SFA sources

| Asset | SFA | Physics |
|---|---|---|
| ball | `v72/results/qflow_w145_eta025.sfa` | Stationary η-soliton (internal clock process) |
| ring | `v73/results/ring_R4.sfa` | Spinning Q-ring (real-space fabric circuit) |
| pair | `/space/scp/v72/qfi_p025.sfa` | Interacting pair (inter-particle correlation) |
| boost | `/space/scp/v73/fx_boost.sfa` | Pass-through: uptake / layment of fabric |
| swirl | `/space/scp/v73/swring_k.sfa` | η-driven ring fission (process split) |
| gring | `/space/scp/v73/gring_k.sfa` | Gauged spinning ring (charge + spin + EM) |

## Stills (cinematic)

| File | Idea |
|---|---|
| `stills/ball_fabric_clock.jpg` | Particle at rest = pure internal re-constitution; envelope quiet, clock turns |
| `stills/ring_fabric_circuit.jpg` | Particle as closed real-space circulation; L_z = nQ |
| `stills/pair_fabric_threads.jpg` | Two processes sharing a medium; correlation threads |
| `stills/boost_fabric_pass.jpg` | Motion = conducting the process through fabric; not transporting matter |
| `stills/swirl_fabric_fission.jpg` | Unclosed / unstable circuit splits into two balls |

## Videos (6s, 720p)

| File | Motion |
|---|---|
| `videos/ball_fabric_clock.mp4` | Push-in; internal filament rotation |
| `videos/ring_fabric_circuit.mp4` | Orbit; toroidal current circulation |
| `videos/pair_fabric_threads.mp4` | Push-in; soft phase shimmer between cores |
| `videos/boost_fabric_pass.mp4` | Drift + leading uptake / trailing layment |
| `videos/swirl_fabric_fission.mp4` | Pull-back; lobes separate, threads thin |

## Scientific (faithful to data)

- `scientific/*_rho2.png`, `*_s.png`, `*_rhoQ.png` — single-panel simulation slices
- `scientific/{ball,ring,pair,boost,swirl,gring}.png` — full multi-quantity time panels
- Raw volview webps also under `images/sfa_captures/`
