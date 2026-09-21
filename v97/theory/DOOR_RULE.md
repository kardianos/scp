# The door selection rule is the Fourier transform of the door profile
(2026-08-13, PLAN step 3b — Go numerics `calc/door.go`, tests green)

**Pre-registered prediction (PLAN §3b): an oriented edge door
concentrates transfer at Δm = ±1 (a dipole selection rule). VERDICT:
REFUTED.** The correct closed form is stronger and simpler, and it
reshapes what the P3 slot-door face can and cannot buy.

## The theorem

On C_N with harmonics χ_m(v) = e^{i2πmv/N}/√N, any **diagonal** door
profile d_v has matrix element

    M(m′←m) = (1/N) Σ_v d_v e^{i2π(m−m′)v/N}  —  the DFT of d at Δm.

So the selection rule in Δm **is** the Fourier transform of the door's
spatial profile, and locality bounds selectivity (uncertainty): a door
supported on w adjacent sites has support of width ≳ N/w in Δm.

| door | profile | \|M(Δm)\| | reading |
|---|---|---|---|
| vertex | δ₀ | 1/N, flat | the measured winding scrambler (qgraph L5) |
| edge, unsigned | δ₀+δ₁ | (2/N)·\|cos(πΔm/N)\| | mild low-pass, still broad |
| **edge, oriented** | δ₀−δ₁ | (2/N)·\|sin(πΔm/N)\| | **anti-selective**: elastic channel closed (Δm=0 forbidden), transfer RISES to the zone boundary |
| wound k | e^{i2πkv/N} | δ_{Δm,−k mod N} | exact pump (SHEAR-G) |
| link current J₀₁ | i(\|1⟩⟨0\|−\|0⟩⟨1\|) | (2/N)·\|sin(π(m+m′)/N)\| | selects on the **sum** m+m′, not Δm |

## What this decides

1. **No door on one link can preserve winding by linear selection.**
   The oriented edge door is the *worst* linear choice for an m-carrier:
   it forbids the elastic channel and pumps hardest toward the zone
   boundary (the neutral, current-free m = N/2 singlet — qgraph Thm
   6.2's dead mode). Orientation buys a **first moment** (sep — the
   leads measured exactly that); it does not buy m-protection.
2. **m-selectivity requires extension around the cycle.** Only a
   profile coherent over the whole loop (wound door) is Δm-sharp.
   This is a theory-level argument for H3's requirement that the
   machine be the *cycle* — the door must be the cycle's door, not a
   link's — and for the helix as the polarized version.
3. **The nonlinear escape is open, not closed.** promo=3's
   |cos Δθ|^p closure gate couples the door to the mode's own phase
   field — outside this linear analysis. Whether a phase-gated local
   door evades the uncertainty bound is a computable next question
   (small numerics, same package), and is the honest theoretical
   content of the "cycle-closure" half of P3.

Prediction MISS recorded; the uncertainty statement replaces it as
the pre-registered expectation for any future slot-door face: **an
oriented one-link door must scramble winding at least as fast as the
vertex door unless phase-gated or extended.**
