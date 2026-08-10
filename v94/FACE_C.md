# v93 — FACE C: the condensation lane (mass formation)

2026-08-08. RESUME §7-C / RING_DNLS route B. The unitary channel + the law's
q_detune detuning spontaneously condense a spread packet into long-lived
dense hoards (B.1/B.2 mapped the existence region at T=80). This face
characterizes the condensation as MASS FORMATION: long-T stability, the hoard
spectrum, and whether the V3a radiance tax selects hoard sizes. Strang
(hop_order=1, Face A); RING_DNLS route B used sequential. T=300, exp=blob.

## C.1 Condensation is stable mass formation (persists to t=300)

| arm (qd, amp) | Em_max@0 | Em_max@40 | Em_max@300 | PR@300 | rms@300 |
|---|---|---|---|---|---|
| law (1.2, 0.5) | 0.56 | 1.81 | **2.17** | 203 | 6.0 |
| deep (3.6, 0.5) | 0.56 | 1.72 | **1.94** | 210 | 5.2 |
| deep (3.6, 2) | 2.25 | 7.90 | **7.24** | 250 | 4.5 |
| melt ctl (qd=0) | 0.56 | 0.78 | **0.52** | 1246 | 7.9 |

- The self-trapped hoards HOLD their peak to t=300 (stable equilibria, not
  transient). The deep corner (qd=3.6, amp=2) holds Em_max~7 (3x cap) with a
  frozen envelope (rms 4.30->4.51) — a dense stable object from a spread
  packet with NO lock/gate/door (the additive law cannot do this). This
  echoes QUENCH ("dynamics created what gates could not").
- The melt control (qd=0) disperses (Em_max->0.52, PR->1246, rms->7.9):
  q_detune is the essential binding nonlinearity.

## C.2 The hoard spectrum is broad/hierarchical (not monodisperse)

Final-frame Em distribution of cells >0.5 (the hoards):

| arm | n_hoards | hoard_Em_tot | spectrum shape |
|---|---|---|---|
| law (1.2,0.5) | 38 | 34.9 | broad 0.5->2.17, few large |
| deep (3.6,2) | **189** | 244.7 (69% of Em) | **long-tailed 0.5->7.2** |
| melt ctl | 1 | 0.5 | (dispersed) |

Deep corner: 99 hoards in [0.5,1.06), 46 in [1.06,1.62), 19 in [1.62,2.19),
11 in [2.19,2.75), 8 in [2.75,3.31), then rare giants to 7.24. A clustering
hierarchy — many small hoards, few large — power-law-like, NOT a single
characteristic mass.

## C.3 Radiance SELECTS hoard sizes — truncates at the fixed point

V3a radiance (k_rad=0.05, p_rad=4; fixed point x*=0.62, Em*=1.55):

| arm | Em_max@300 no-rad | Em_max@300 V3a | rad@300 | verdict |
|---|---|---|---|---|
| law (1.2,0.5) | 2.17 | **1.57** | 4.7 | -> fixed pt 1.55 |
| deep (3.6,0.5) | 1.94 | **1.32** | 2.3 | eroded below fp |
| deep (3.6,2) | 7.24 | **6.24** | 89.9 | deep self-trap resists |

- **Radiance drives the hoard peak to the fixed point.** Law corner Em_max
  2.17->1.57 (~1.55 fixed point); the spectrum's >1.75 hoards VANISH under
  V3a (33 hoards all <=1.62). The tax equilibrates the hoard peak to x*.
- The deep self-trap (qd=3.6, amp=2) RESISTS: Em_max 7.24->6.24 (still 2.5x
  cap) despite rad=89.9 — the conservative self-trap holds above cap. The
  mid-size hoards (2.4-6) erode, but the deepest one survives. Spectrum goes
  bimodal: many small (<2) + one giant (6.24).
- Radiance is a mass-dependent erosion: it truncates the hoard spectrum at
  the fixed point for mild/moderate hoards, but deep self-traps survive.

## Synthesis

1. **Spontaneous condensation is a genuine mass-formation mechanism.** Long-
   lived (t>=300), stable, hierarchical dense hoards form from a spread
   packet under the unitary channel + the law's own detuning, with no lock,
   gate, or door. The additive law cannot do this (B.5). Creation-adjacent.
2. **The hoard spectrum is broad/power-law-like**, not a characteristic mass
   — a clustering hierarchy from many small to few large hoards.
3. **Radiance is the mass selector.** The V3a tax truncates the spectrum at
   the radiance fixed point (Em*~1.55), eroding hoards above it — except the
   deepest self-traps, which hold above cap. This is "space flows: pressure
   pushes" implemented as hoard-size selection.
4. Connects to QUENCH (dynamics-created matter) and to the frozen-ring
   internal condensation of Face A §A.4 (same q_detune rich-get-richer).

## Files
- Meters: `report/analyze_hoard.py` (final-frame Em spectrum).
- Script: `runs/faceC.sh` -> `runs/faceC/{c_*,*.melt}`.
