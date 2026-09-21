# v98 Phase 2 — size spectrum (2026-08-13)

**Authorized by** `phase1/DECIDE.md` §5. **No new law. No
kernel face.** The live tree is still v97. This phase is a
meter on existing FCS streams.

The claim is: *is the histogram a line, a hierarchy, or a
tax ceiling?* Not *is there a dipole?*

---

## 0. What was built

| path | role |
|---|---|
| `size/` | census + `N_eff` + unweighted Laplacian 3+1 (CHARTEXP Jacobi, read-only) |
| `cmd/sizespec` | CLI |
| `phase2/BARS.md` | pre-registered bars |
| `phase2/H0.md` | honesty grades |
| `phase2/runs/` | raw H0 prints |

```
cd v98 && go test ./size/
go run ./cmd/sizespec -h0 broad|line6|ceiling FILE.fcs
go run ./cmd/sizespec FILE.fcs          # per-frame TSV
```

`N_eff = (Σ Em)² / Σ Em²`. Chart write-back into `w2e` is
not implemented and must not be.

---

## 1. H0

All three named classes pass on the streams the charter
named. The three opposite pairings fail. Details and
numbers: `phase2/H0.md`.

| class | named stream | H0 |
|---|---|---|
| broad | Face C `c_deep2` | PASS |
| line6 | vacuum C6 `ringA_strang_live` | PASS |
| ceiling | `rad_tax_S05` / MONO `M_qd12_kr05` | PASS |

The meter is honest enough to use.

---

## 2. Interpretation (what Phase 2 does *not* say)

- A vacuum C6 line is **not** a nucleon. It is the H0
  existence proof that a line is visible. v97: loop
  topology, Φ = 0, no I/O split.
- Face C being broad is **not** new physics. It is the
  Face C / MONO fact restated in `N` instead of `Em_max`.
- Tax being a ceiling is **not** a species. Inventory
  falls; the histogram stays a zoo.
- Phase 2 does **not** close the 1B gap (`ρ` unforced).
- Phase 2 does **not** authorize Phase 3.

The world’s proton, if it appeared in a stream, would have
to grade **line6-shaped** (one compact `N`, stable, a
discrete Laplacian book) *and* survive the unfold filter
(coarse = one size, fine = field). No existing stream does
that except the vacuum ring, which fails the unfold (it
does not resolve into field; it is six cells of parked Em)
and fails H1 (no live conversion).

---

## 3. What is next

Phase 3 remains unauthorized (`DECIDE.md` §5). A face
would have to force `ρ` or state a selector 1A/1B already
refused, and it would need explicit sign-off.

Until then the size meter is the scoreboard: run it on any
new stream before claiming a species.

---

## 4. Reproduce H0

```
cd /home/d/code/scp/v98
go test ./size/
go run ./cmd/sizespec -h0 broad   ../v93/runs/faceC/c_deep2.fcs
go run ./cmd/sizespec -h0 line6   ../v93/runs/faceA/ringA_strang_live.fcs
go run ./cmd/sizespec -h0 ceiling /space/scp/v93/fcs/rad_tax_S05.fcs
# negatives:
go run ./cmd/sizespec -h0 line6   ../v93/runs/faceC/c_deep2.fcs
go run ./cmd/sizespec -h0 broad   ../v93/runs/faceA/ringA_strang_live.fcs
go run ./cmd/sizespec -h0 ceiling ../v93/runs/faceA/ringA_strang_live.fcs
```
