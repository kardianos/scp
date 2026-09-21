# SUSTAIN REPORT

T=60 seed=20260802 kernel=`/home/d/code/scp/v93/fabrun` n=0

## Criterion

A condensate **sustains** if over the last half of [0,T]:
- `|slope|` small vs early growth, **or** Em stays above a floor with bounded PR;
- not pure melt (Em→0 / PR→N);
- preferably net book `cond - evap - rad` near zero if doors are open.

| arm | Em_end | Em_max | slope/1k | cond | evap | rad | fd | verdict |
|---|---:|---:|---:|---:|---:|---:|---:|---|

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
