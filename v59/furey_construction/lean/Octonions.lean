/-
  v59/furey_construction/lean/Octonions.lean

  Minimal octonion algebra: the multiplication table via the Fano plane.

  Octonions are not associative.  We record the multiplication structure
  (signs and products) for the 8 basis elements e_0=1, e_1, ..., e_7.
-/

namespace SCPv59.Octonions

/-- Octonion basis indices: 0 (identity) through 7. -/
def OctIndex : Type := Fin 8

/-- Fano-plane triples: each triple (a, b, c) means e_a * e_b = +e_c
    (and cyclic), with reversed products giving the negative.

    Standard Cayley convention. -/
def fano_triples : List (Nat × Nat × Nat) :=
  [ (1, 2, 3),
    (1, 4, 5),
    (1, 7, 6),
    (2, 4, 6),
    (2, 5, 7),
    (3, 4, 7),
    (3, 6, 5) ]

/-- THEOREM: There are 7 Fano lines (Fano-plane axiom). -/
theorem n_fano_triples : fano_triples.length = 7 := by
  unfold fano_triples
  rfl

/-- THEOREM: Total point-line incidences = 7 lines × 3 points = 21. -/
def n_incidences : Nat := fano_triples.length * 3

theorem incidences_is_21 : n_incidences = 21 := by
  unfold n_incidences fano_triples
  rfl

/-- The 21 incidences match dim Spin(7) (= 21).
    This is a STRUCTURAL claim, not derived here -- it is the empirical
    finding of v59 step 9 (octonionic_extension/01_findings.md).
    The number 21 appears in TWO places:
      - As the total incidences of the Fano plane.
      - As the dimension of Spin(7), the natural symmetry of the octonion imaginary sector.
    These are not numerically distinct; they are the same integer 21 arising
    in connected ways. -/
example : n_incidences = 3 * 7 := by
  unfold n_incidences fano_triples
  rfl

end SCPv59.Octonions
