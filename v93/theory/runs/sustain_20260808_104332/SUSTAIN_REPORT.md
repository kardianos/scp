# SUSTAIN REPORT

T=80 seed=20260802 kernel=`/home/d/code/scp/v93/fabrun` n=2

## Criterion

A condensate **sustains** if over the last half of [0,T]:
- `|slope|` small vs early growth, **or** Em stays above a floor with bounded PR;
- not pure melt (Em→0 / PR→N);
- preferably net book `cond - evap - rad` near zero if doors are open.

| arm | Em_end | Em_max | slope/1k | cond | evap | rad | fd | verdict |
|---|---:|---:|---:|---:|---:|---:|---:|---|
| S00_uni_base | 148.8 | 2.17 | +0.00 | 0.0 | 0.00 | 0.00 | +1.1 | **incomplete-series** |
| S03_qd0_melt | 148.8 | 0.67 | +0.00 | 0.0 | 0.00 | 0.00 | +0.0 | **incomplete-series** |

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
