# SUSTAIN REPORT

T=600 seed=20260802 kernel=`/home/d/code/scp/v93/fabrun` n=9

## Criterion

A condensate **sustains** if over the last half of [0,T]:
- `|slope|` small vs early growth, **or** Em stays above a floor with bounded PR;
- not pure melt (Em→0 / PR→N);
- preferably net book `cond - evap - rad` near zero if doors are open.

| arm | Em_end | Em_max | slope/1k | cond | evap | rad | fd | verdict |
|---|---:|---:|---:|---:|---:|---:|---:|---|
| P2_ctrl_kr0 | 6258.2 | 7.76 | +2618.97 | 35726.6 | 34770.03 | 0.00 | +0.0 | **GROW** |
| P2_xb03_d0 | 4214.3 | 5.26 | -17.38 | 30056.6 | 21188.96 | 9285.46 | +0.0 | **HOLD** |
| P2_xb03_d1 | 4568.4 | 4.52 | +5.64 | 29094.8 | 22937.91 | 11090.06 | +0.0 | **HOLD** |
| P2_xb06_d0 | 3883.7 | 4.96 | +89.59 | 46731.0 | 41986.07 | 7707.60 | +0.0 | **GROW** |
| P2_xb06_d1 | 5360.8 | 7.14 | +259.73 | 41388.3 | 48444.38 | 12798.62 | +0.0 | **GROW** |
| P2_xb10_d0 | 3228.6 | 5.61 | +525.06 | 68551.7 | 68349.17 | 6059.00 | +0.0 | **GROW** |
| P2_xb10_d1 | 5521.5 | 6.97 | +403.59 | 43022.6 | 52262.12 | 13380.12 | +0.0 | **GROW** |
| P2_xb15_d0 | 2229.5 | 5.80 | +1073.37 | 99262.6 | 103562.17 | 5312.91 | +0.0 | **GROW** |
| P2_xb15_d1 | 6098.2 | 8.67 | +1108.73 | 46171.0 | 61713.50 | 13948.51 | +0.0 | **GROW** |

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
