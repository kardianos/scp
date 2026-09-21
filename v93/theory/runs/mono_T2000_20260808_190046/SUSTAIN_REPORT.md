# SUSTAIN REPORT

T=2000 seed=20260802 kernel=`/home/d/code/scp/v93/fabrun` n=8

## Criterion

A condensate **sustains** if over the last half of [0,T]:
- `|slope|` small vs early growth, **or** Em stays above a floor with bounded PR;
- not pure melt (Em→0 / PR→N);
- preferably net book `cond - evap - rad` near zero if doors are open.

| arm | Em_end | Em_max | slope/1k | cond | evap | rad | fd | verdict |
|---|---:|---:|---:|---:|---:|---:|---:|---|
| M_qd12_kr00 | 148.3 | 2.08 | +0.00 | 0.0 | 0.47 | 0.00 | +0.0 | **HOLD** |
| M_qd12_kr02 | 124.8 | 1.52 | -6.95 | 0.8 | 0.00 | 24.93 | +0.0 | **HOLD** |
| M_qd12_kr05 | 113.2 | 1.31 | -7.93 | 3.3 | 0.00 | 39.37 | +0.0 | **HOLD** |
| M_qd12_kr10 | 105.4 | 0.84 | -7.11 | 4.9 | 0.00 | 49.03 | +0.0 | **HOLD** |
| M_qd36_kr00 | 147.7 | 2.51 | -0.43 | 0.0 | 1.11 | 0.00 | +0.0 | **HOLD** |
| M_qd36_kr02 | 132.1 | 1.32 | -4.72 | 6.2 | 0.00 | 23.76 | +0.0 | **HOLD** |
| M_qd36_kr05 | 122.5 | 1.27 | -6.58 | 10.3 | 0.57 | 37.52 | +0.0 | **HOLD** |
| M_qd36_kr10 | 116.3 | 1.07 | -5.75 | 16.8 | 0.11 | 51.68 | +0.0 | **HOLD** |

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
