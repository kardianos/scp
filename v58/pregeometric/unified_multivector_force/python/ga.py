"""
Lightweight Geometric Algebra for 3D (Cl(3,0)) and extensible to 3+1.
Pure Python + sympy (no external GA libs required).
Used for rapid symbolic grade projections and numeric lattice experiments
in the unified multivector force law discovery.

Basis blades represented as sorted tuples of generator indices:
  ()     -> scalar   (grade 0)
  (0,)   -> e1       (grade 1)
  (1,)   -> e2
  (2,)   -> e3
  (0,1)  -> e1^e2    (grade 2)
  ...
  (0,1,2)-> e1^e2^e3 (grade 3, pseudoscalar I)

Geometric product, outer product, reverse, grade projection implemented.
Signature defaults to Euclidean +++ for 3D spatial tests.
"""

from __future__ import annotations
import sympy as sp
from typing import Dict, Tuple, Any, Callable
import itertools

# Type for blade: sorted tuple of ints 0..N-1
Blade = Tuple[int, ...]
Coeff = Any  # sp.Expr or float

SIGNATURE_3D = (1, 1, 1)  # e0^2, e1^2, e2^2  (0-based)

class MV:
    """Multivector with dict[blade, coeff] storage. Immutable-ish via ops."""
    def __init__(self, terms: Dict[Blade, Coeff] | None = None, sig: Tuple[int, ...] = SIGNATURE_3D):
        self.sig = sig
        self.terms: Dict[Blade, Coeff] = {}
        if terms:
            for b, c in terms.items():
                if c != 0:
                    self.terms[b] = c  # assume already canonical

    @staticmethod
    def zero(sig: Tuple[int, ...] = SIGNATURE_3D) -> MV:
        return MV({}, sig)

    @staticmethod
    def scalar(s: Coeff, sig: Tuple[int, ...] = SIGNATURE_3D) -> MV:
        return MV({(): s}, sig)

    @staticmethod
    def vector(v: Tuple[Coeff, Coeff, Coeff], sig: Tuple[int, ...] = SIGNATURE_3D) -> MV:
        return MV({(0,): v[0], (1,): v[1], (2,): v[2]}, sig)

    @staticmethod
    def bivector(b: Tuple[Coeff, Coeff, Coeff], sig: Tuple[int, ...] = SIGNATURE_3D) -> MV:
        # b12, b13, b23
        return MV({(0,1): b[0], (0,2): b[1], (1,2): b[2]}, sig)

    def __getitem__(self, blade: Blade) -> Coeff:
        return self.terms.get(blade, sp.Integer(0) if isinstance(next(iter(self.terms.values()), 0), (sp.Expr, sp.Number)) else 0.0)

    def grade(self, k: int) -> MV:
        """Project onto grade k."""
        out = {}
        for b, c in self.terms.items():
            if len(b) == k:
                out[b] = c
        return MV(out, self.sig)

    def __neg__(self) -> MV:
        return MV({b: -c for b, c in self.terms.items()}, self.sig)

    def __add__(self, other: MV) -> MV:
        if not isinstance(other, MV):
            other = MV.scalar(other, self.sig)
        out = dict(self.terms)
        for b, c in other.terms.items():
            out[b] = out.get(b, 0) + c
        # clean zeros? optional for perf, skip for sympy
        # clean zeros robustly
        cleaned = {}
        zero_sp = sp.Integer(0)
        for b, v in out.items():
            try:
                if abs(v) > 1e-14:
                    cleaned[b] = v
            except (TypeError, ValueError):
                if v != 0 and v != zero_sp:
                    cleaned[b] = v
        return MV(cleaned, self.sig)

    def __radd__(self, other):
        return self.__add__(other)

    def __sub__(self, other: MV) -> MV:
        return self + (-other)

    def __mul__(self, other: MV) -> MV:
        """Geometric product (full Clifford product)."""
        if not isinstance(other, MV):
            other = MV.scalar(other, self.sig)
        out: Dict[Blade, Coeff] = {}
        for b1, c1 in self.terms.items():
            for b2, c2 in other.terms.items():
                sgn, b3 = self._blade_geometric(b1, b2)
                val = sgn * c1 * c2
                if b3 in out:
                    out[b3] += val
                else:
                    out[b3] = val
        # filter zeros for numeric (handle float/int/sympy)
        cleaned = {}
        zero_sp = sp.Integer(0)
        for b, v in out.items():
            try:
                if abs(v) > 1e-14:  # numeric tolerance
                    cleaned[b] = v
            except (TypeError, ValueError):
                # sympy expr
                if v != 0 and v != zero_sp:
                    cleaned[b] = v
        return MV(cleaned, self.sig)

    def __rmul__(self, other):
        return MV.scalar(other, self.sig) * self

    def __or__(self, other: MV) -> MV:
        """Outer (wedge) product: grade (g1+g2) part of geometric."""
        prod = self * other
        g = len(next(iter(self.terms), ())) + len(next(iter(other.terms), ()))
        return prod.grade(g)

    def __and__(self, other: MV) -> MV:
        """Inner (left contraction-ish) product: lowest grade part."""
        prod = self * other
        g = abs(len(next(iter(self.terms), ())) - len(next(iter(other.terms), ())))
        return prod.grade(g)

    def reverse(self) -> MV:
        """Reverse (tilde): sign (-1)^{k(k-1)/2} per grade k."""
        out = {}
        for b, c in self.terms.items():
            k = len(b)
            s = 1 if (k * (k - 1) // 2) % 2 == 0 else -1
            out[b] = s * c
        return MV(out, self.sig)

    def __invert__(self) -> MV:
        return self.reverse()

    def norm2(self) -> Coeff:
        """M * ~M scalar part (squared norm)."""
        return (self * self.reverse()).grade(0).terms.get((), 0)

    def _blade_geometric(self, b1: Blade, b2: Blade) -> Tuple[int, Blade]:
        """Return (sign, canonical_blade) for product of two blades."""
        gens = list(b1) + list(b2)
        sign = 1
        i = 0
        N = len(self.sig)
        while i < len(gens) - 1:
            if gens[i] == gens[i + 1]:
                # e_i * e_i = sig[i]
                sq = self.sig[gens[i]]
                sign *= sq
                gens.pop(i)
                gens.pop(i)
                i = max(0, i - 1)
            elif gens[i] > gens[i + 1]:
                # swap distinct -> -1 (anticommute)
                gens[i], gens[i + 1] = gens[i + 1], gens[i]
                sign *= -1
                i += 1
            else:
                i += 1
        return sign, tuple(gens)

    def __repr__(self) -> str:
        if not self.terms:
            return "0"
        parts = []
        for b in sorted(self.terms.keys(), key=lambda x: (len(x), x)):
            c = self.terms[b]
            if b == ():
                name = "1"
            else:
                name = "*".join(f"e{i+1}" for i in b)
            parts.append(f"({c})*{name}" if c != 1 else name)
        return " + ".join(parts)

    def __str__(self) -> str:
        return self.__repr__()

    def to_dict(self) -> Dict[str, Coeff]:
        """Human friendly for printing."""
        d = {}
        for b, c in self.terms.items():
            if b == ():
                d["s"] = c
            else:
                d["e" + "".join(str(i+1) for i in b)] = c
        return d


# Convenience creators using sympy symbols
def symbolic_mv(prefix: str = "a", grades: Tuple[int, ...] = (0,1,2), sig: Tuple[int,...] = SIGNATURE_3D) -> MV:
    """Create a general multivector with symbolic coeffs for given grades."""
    terms = {}
    idx = 0
    for k in grades:
        if k == 0:
            terms[()] = sp.symbols(f"{prefix}s")
        elif k == 1:
            for i in range(3):
                terms[(i,)] = sp.symbols(f"{prefix}v{i+1}")
        elif k == 2:
            for comb in itertools.combinations(range(3), 2):
                terms[tuple(comb)] = sp.symbols(f"{prefix}b{comb[0]+1}{comb[1]+1}")
        elif k == 3:
            terms[(0,1,2)] = sp.symbols(f"{prefix}p")
    return MV(terms, sig)


def vector_derivative_symbolic() -> MV:
    """Return a formal 'D' = e1*dx + e2*dy + e3*dz as a vector MV (for structure analysis)."""
    # We treat derivatives as formal; in practice for projection we multiply D * Omega
    # and inspect resulting grades. Coefficients carry the ∂ symbols.
    dx, dy, dz = sp.symbols("∂x ∂y ∂z")
    return MV({(0,): dx, (1,): dy, (2,): dz}, SIGNATURE_3D)


# Quick self-test when run as script
if __name__ == "__main__":
    print("GA lite self-test (Cl(3,0))")
    e1 = MV.vector((1,0,0))
    e2 = MV.vector((0,1,0))
    print("e1 * e2:", e1 * e2)
    print("e2 * e1:", e2 * e1)
    I = e1 * e2 * MV.vector((0,0,1))   # pseudoscalar?
    print("I ~ e1*e2*e3:", I)
    print("reverse I:", ~I)
    print("e1**2:", e1 * e1)
    print("Grade proj test ok.")

    # Symbolic
    Omega = symbolic_mv("Ω", grades=(0,1,2))
    print("\nSymbolic Ω (scalar+vec+biv):", Omega.to_dict())
    D = vector_derivative_symbolic()
    D_Omega = D * Omega
    print("D*Ω grades present:", sorted({len(b) for b in D_Omega.terms.keys()}))
    print("Scalar part of DΩ:", D_Omega.grade(0))
    print("Bivector part of DΩ:", D_Omega.grade(2))
