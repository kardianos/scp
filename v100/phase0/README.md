# v100 phase 0

Pre-registered 2026-09-21, then run the same day. No kernel.
No mesh.

| | file | what it is |
|---|---|---|
| gates | `PREREG.md` | locked before the run |
| spin | `SPIN.md` | G2, written without the profile |
| profile | `PROFILE.md` | the hedgehog numbers |
| decision | `DECIDE.md` | G1 pass, G2 fail; hedgehog refused; phase 1 not opened |
| instrument | `cas/profile.py`, `cas/virial.mac` | |
| logs | `cas/out/` | |

```
python3 v100/phase0/cas/profile.py
maxima --very-quiet -b v100/phase0/cas/virial.mac < /dev/null
```
