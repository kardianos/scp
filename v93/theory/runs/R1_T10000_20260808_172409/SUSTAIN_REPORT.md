# SUSTAIN REPORT

T=10000 seed=20260802 kernel=`/home/d/code/scp/v93/fabrun` n=8

## Criterion

A condensate **sustains** if over the last half of [0,T]:
- `|slope|` small vs early growth, **or** Em stays above a floor with bounded PR;
- not pure melt (Em→0 / PR→N);
- preferably net book `cond - evap - rad` near zero if doors are open.

| arm | Em_end | Em_max | slope/1k | cond | evap | rad | fd | verdict |
|---|---:|---:|---:|---:|---:|---:|---:|---|
| R1_S00 | 148.3 | 1.71 | +0.00 | 0.0 | 0.47 | 0.00 | +90.8 | **HOLD** |
| R1_S00_amp12 | 206.7 | 1.80 | -0.03 | 4.6 | 21.70 | 0.00 | +0.0 | **HOLD** |
| R1_S00_rad | 105.4 | 0.75 | -0.54 | 3.3 | 0.00 | 47.11 | +0.0 | **HOLD** |
| R1_S00_s111 | 146.8 | 1.68 | +0.00 | 0.0 | 0.00 | 0.00 | +0.0 | **HOLD** |
| R1_S04 | 147.3 | 1.99 | -0.00 | 0.1 | 1.55 | 0.00 | +0.0 | **HOLD** |
| R1_S04_rad | 113.9 | 0.66 | -0.72 | 11.2 | 0.57 | 47.12 | +0.0 | **HOLD** |
| R1_S14 | 147.7 | 1.56 | -0.00 | 0.0 | 1.05 | 0.00 | +741.3 | **HOLD** |
| R1_S16 | 148.1 | 1.66 | +0.00 | 0.0 | 0.70 | 0.00 | +0.0 | **HOLD** |

## Isolation map

| switch | arms |
|---|---|
| transport uni vs add | S00 vs S01 |
| gate | S00 vs S02 |
| q_detune trap | S03 / S00 / S04 |
| radiance outtake | S00 vs S05; S04 vs S06 |
| field fuel (noise) | S00 vs S07; S05 vs S08 |
| workfn door | S00 vs S10; S07 vs S11 |
| metabolic (open door+fuel+tax) | S12 |
| progressive dense phase | S14, S15, S16 |
| genesis long shear | G01 vs G02 |

## Bohmian-method aspects

See SUSTAIN_BOHMIAN.md / report body: progressive meters (J, fd, support), not pilot particles. S14–S16 vary phase progression on dense amplitude.
