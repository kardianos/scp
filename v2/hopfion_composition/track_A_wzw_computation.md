# Track A: WZW 5-Form on Full Cl+(3,0,1) — Does It Produce B^0 tau Coupling?

## 1. Setup and Conventions

### 1.1 The Algebra Cl(3,0,1)

Generators: e_0 (null), e_1, e_2, e_3 (Euclidean).

Squares: e_0^2 = 0, e_1^2 = e_2^2 = e_3^2 = 1.

All distinct generators anticommute: e_a e_b = -e_b e_a for a != b.

Compound elements used throughout:
- Spatial bivectors: e_{23}, e_{31}, e_{12} (these square to -1, form quaternion basis)
- Mixed bivectors: e_{01}, e_{02}, e_{03} (these square to 0)
- Spatial trivector: e_{123} (squares to -1)
- Pseudoscalar: e_{0123} (squares to 0, COMMUTES with all even-grade elements)

Key product table (geometric products of basis bivectors with mixed bivectors):

    e_{23} * e_{01} = e_2 e_3 e_0 e_1

We sort using anticommutativity. Move e_0 left past e_2, e_3:
e_2 e_3 e_0 = e_2 (-e_0 e_3) = -e_2 e_0 e_3 = -(-e_0 e_2) e_3 = e_0 e_2 e_3

So e_{23} * e_{01} = e_0 e_2 e_3 e_1 = -e_0 e_2 e_1 e_3 = e_0 e_1 e_2 e_3 = e_{0123}.

Similarly:

    e_{31} * e_{02} = e_3 e_1 e_0 e_2

Move e_0 left: e_3 e_1 e_0 = e_3(-e_0 e_1) = -e_3 e_0 e_1 = (-e_3 e_0)e_1 = e_0 e_3 e_1.
So = e_0 e_3 e_1 e_2 = -e_0 e_3 e_2 e_1 ... wait, let me be more careful.

    e_3 e_1 e_0 e_2:
    Step 1: e_1 e_0 = -e_0 e_1, so e_3(-e_0 e_1)e_2 = -e_3 e_0 e_1 e_2
    Step 2: e_3 e_0 = -e_0 e_3, so -(-e_0 e_3)e_1 e_2 = e_0 e_3 e_1 e_2
    Step 3: e_3 e_1 = -e_1 e_3, so e_0(-e_1 e_3)e_2 = -e_0 e_1 e_3 e_2
    Step 4: e_3 e_2 = -e_2 e_3, so -e_0 e_1(-e_2 e_3) = e_0 e_1 e_2 e_3 = e_{0123}.

And:

    e_{12} * e_{03} = e_1 e_2 e_0 e_3

    Step 1: e_2 e_0 = -e_0 e_2, so e_1(-e_0 e_2)e_3 = -e_1 e_0 e_2 e_3
    Step 2: e_1 e_0 = -e_0 e_1, so -(-e_0 e_1)e_2 e_3 = e_0 e_1 e_2 e_3 = e_{0123}.

Cross-terms (spatial bivector times non-matching mixed bivector):

    e_{23} * e_{02} = e_2 e_3 e_0 e_2

    Move e_0 left: e_2 e_3 e_0 = e_0 e_2 e_3 (from above).
    So = e_0 e_2 e_3 e_2 = e_0 e_2(-e_2 e_3) = -e_0 e_2^2 e_3 = -e_0 e_3 = -e_{03}.

This is grade 2, not grade 4. Similarly all cross-terms produce grade-2 results.

Products involving e_{0123}:

    e_{23} * e_{0123}: since e_{0123} commutes with even-grade elements,
    e_{23} * e_{0123} = e_{0123} * e_{23} = e_0 e_1 e_2 e_3 e_2 e_3
    = e_0 e_1 (e_2 e_3 e_2 e_3) = e_0 e_1 (e_2 e_3 e_2 e_3)
    = e_0 e_1 (-e_2 e_2 e_3 e_3) ... no.
    e_2 e_3 e_2 e_3 = e_2(e_3 e_2)e_3 = e_2(-e_2 e_3)e_3 = -e_2^2 e_3^2 = -1.
    So e_{23} * e_{0123} = -e_{01}.

Similarly: e_{31} * e_{0123} = -e_{02}, and e_{12} * e_{0123} = -e_{03}.

And: 1 * e_{0123} = e_{0123} (scalar times pseudoscalar).

### 1.2 The Field

    Psi = q + epsilon * d

where epsilon = e_0 (the null basis vector), and:

    q = s + f_1 e_{23} + f_2 e_{31} + f_3 e_{12}       (quaternion, grades 0+2)
    d = j_1 e_1 + j_2 e_2 + j_3 e_3 + tau e_{123}      (odd elements of Cl(3,0))

Wait — let me be precise about the decomposition. The degenerate sector has basis
elements e_{01}, e_{02}, e_{03}, e_{0123}. We can factor out e_0:

    epsilon * d = e_0(j_1 e_1 + j_2 e_2 + j_3 e_3 + tau e_{123})

So d (the part AFTER factoring out e_0) lives in the ODD part of Cl(3,0):
d = j_1 e_1 + j_2 e_2 + j_3 e_3 + tau e_{123}.

Check: e_0 * (j_1 e_1) = j_1 e_{01} (grade-2 bivector). Good.
e_0 * (tau e_{123}) = tau e_{0123} (grade-4 pseudoscalar). Good.

### 1.3 Reversal

For the FULL algebra element Psi = q + e_0 d:

Reversal acts as: (AB)~ = B~ A~, applied to each basis blade.

    q~ = s - f_1 e_{23} - f_2 e_{31} - f_3 e_{12}

For the degenerate part:
    (e_0 d)~ = d~ e_0~ = d~ e_0

where d~ = j_1 e_1 + j_2 e_2 + j_3 e_3 - tau e_{123}
(grade-1 elements are unchanged under reversal; grade-3 elements are negated).

Hmm, but wait: reversal of e_{123} = e_1 e_2 e_3. Reversal gives e_3 e_2 e_1 = -e_{123}
(three swaps: 321 from 123, each swap gives a sign, 3 swaps = -1). Actually:
e_3 e_2 e_1 = e_3(-e_1 e_2) = -e_3 e_1 e_2 = -(-e_1 e_3)e_2 = e_1 e_3 e_2
= e_1(-e_2 e_3) = -e_1 e_2 e_3 = -e_{123}. Yes, e_{123}~ = -e_{123}.

So d~ = j_1 e_1 + j_2 e_2 + j_3 e_3 - tau e_{123}.
And (e_0 d)~ = d~ e_0 = (j_1 e_1 + j_2 e_2 + j_3 e_3 - tau e_{123}) e_0.

Note: e_i e_0 = -e_0 e_i, and e_{123} e_0 = -e_0 e_{123} (since e_{123} is grade-3, odd).

Wait: e_{123} e_0 = e_1 e_2 e_3 e_0. Moving e_0 to the left past three vectors gives
(-1)^3 = -1. So e_{123} e_0 = -e_0 e_{123} = -e_{0123}.

Therefore:
    (e_0 d)~ = j_1(-e_0 e_1) + j_2(-e_0 e_2) + j_3(-e_0 e_3) - tau(-e_0 e_{123})
             = -j_1 e_{01} - j_2 e_{02} - j_3 e_{03} + tau e_{0123}

This confirms the reversal rule stated in the problem: d~ negates e_{0i}, keeps e_{0123}.

So: Psi~ = q~ + (-j_1 e_{01} - j_2 e_{02} - j_3 e_{03} + tau e_{0123}).

### 1.4 Sigma-model constraint

    |Psi|^2 = Psi Psi~ = (q + e_0 d)(q~ + (e_0 d)~)

Expanding:
    = q q~ + q(e_0 d)~ + (e_0 d)q~ + (e_0 d)(e_0 d)~

The last term: (e_0 d)(d~ e_0) = e_0 (d d~) e_0. Since d d~ is a scalar+bivector in Cl(3,0),
and e_0 (scalar) e_0 = scalar * e_0^2 = 0. So the last term vanishes.

The first term: q q~ = |q|^2 = s^2 + f_1^2 + f_2^2 + f_3^2 = rho_0^2 on sigma model.

The cross-terms q(e_0 d)~ + (e_0 d)q~ must vanish on the sigma model. This gives the
constraint: q d~ e_0 + e_0 d q~ = 0, i.e., the degenerate sector is constrained by q.

---

## 2. Full Left Current L = l + epsilon * delta_l

### 2.1 Inverse of Psi

With Psi = q + e_0 d and e_0^2 = 0:

    Psi^{-1} = (q + e_0 d)^{-1}

Factor: Psi = q(1 + q^{-1} e_0 d). Since (q^{-1} e_0 d)^2 = q^{-1} e_0 d q^{-1} e_0 d
contains e_0^2 at some point... actually no. Let A = q^{-1} e_0 d. Then
A^2 = q^{-1} e_0 d q^{-1} e_0 d. This has e_0 ... e_0 but not necessarily adjacent.
However, e_0 d q^{-1} e_0 = e_0 (d q^{-1}) e_0. Since d q^{-1} is a sum of Cl(3,0,1)
elements, and e_0 X e_0 for any X: if X has no e_0, then e_0 X e_0 = X' e_0^2 = 0
(where X' comes from commuting e_0 past X). More precisely:
e_0 X e_0 where X is in Cl(3,0) (no e_0 factors): each term in X is a product of e_i's,
and e_0 commutes past them picking up signs, giving e_0^2 times something = 0.

So A^2 = 0, confirming (1 + A)^{-1} = 1 - A, and:

    Psi^{-1} = (1 - q^{-1} e_0 d) q^{-1} = q^{-1} - q^{-1} e_0 d q^{-1}

### 2.2 The left current

    L_mu = Psi^{-1} partial_mu Psi
         = (q^{-1} - q^{-1} e_0 d q^{-1})(partial_mu q + e_0 partial_mu d)

Expanding (and using e_0^2 = 0 to kill terms with two e_0's):

    L_mu = q^{-1} partial_mu q + q^{-1} e_0 partial_mu d - q^{-1} e_0 d q^{-1} partial_mu q

The term q^{-1} e_0 d q^{-1} e_0 partial_mu d contains e_0^2 = 0, so it vanishes.

Define the bulk left current:

    l_mu = q^{-1} partial_mu q          (lives in Im(H) ~ span{e_{23}, e_{31}, e_{12}})

The degenerate part:

    L_mu = l_mu + e_0(q^{-1} partial_mu d - q^{-1} d q^{-1} partial_mu q)
         = l_mu + e_0(q^{-1} partial_mu d - q^{-1} d * l_mu)

Define:

    delta_l_mu = q^{-1} partial_mu d - q^{-1} d * l_mu = q^{-1}(partial_mu d - d * l_mu)

This is the covariant derivative D_mu(q^{-1} d) with connection l_mu acting on the right:

    delta_l_mu = D_mu(q^{-1} d) where D_mu X = q^{-1} partial_mu(q X)
               = q^{-1} partial_mu d + q^{-1} partial_mu q * q^{-1} d ... hmm.

Let me verify. Let phi = q^{-1} d. Then:
    partial_mu phi = partial_mu(q^{-1}) d + q^{-1} partial_mu d
                   = -q^{-1}(partial_mu q)q^{-1} d + q^{-1} partial_mu d
                   = -l_mu q^{-1} d + q^{-1} partial_mu d
                   = q^{-1} partial_mu d - l_mu phi

But delta_l_mu = q^{-1} partial_mu d - q^{-1} d l_mu = q^{-1} partial_mu d - phi l_mu.

So delta_l_mu = partial_mu phi + l_mu phi - phi l_mu - l_mu phi = partial_mu phi - phi l_mu.

Hmm, that's not quite the standard covariant derivative form. Let me just write:

    delta_l_mu = partial_mu phi - phi l_mu     where phi = q^{-1} d

This is a RIGHT covariant derivative (connection acts from the right).

### 2.3 Algebraic character of delta_l_mu

phi = q^{-1} d lives in Cl(3,0) since q^{-1} is a quaternion and d = j_1 e_1 + j_2 e_2 + j_3 e_3 + tau e_{123}
is an odd element of Cl(3,0). A quaternion times an odd Cl(3,0) element gives an odd Cl(3,0) element.

So phi is odd in Cl(3,0): phi = alpha_1 e_1 + alpha_2 e_2 + alpha_3 e_3 + beta e_{123}
for some functions alpha_i, beta that are linear combinations of j_i, tau weighted by q^{-1} components.

l_mu is a pure imaginary quaternion: l_mu = a_1 e_{23} + a_2 e_{31} + a_3 e_{12}.

phi * l_mu: (odd Cl(3,0)) * (grade-2 Cl(3,0)). Products:
- e_i * e_{jk}: if i = j or i = k, gives a grade-1 vector. If i != j,k, gives grade-3.
- e_{123} * e_{jk}: gives a grade-1 vector (e_{123} e_{jk} = +/- e_m for appropriate m).

So phi * l_mu is again odd in Cl(3,0) (grades 1 and 3). Same for partial_mu phi.

Therefore delta_l_mu lives in the odd part of Cl(3,0):
    delta_l_mu = (vector part) + (trivector part)
               = v_mu^i e_i + w_mu e_{123}

And the full left current:
    L_mu = l_mu + e_0 * delta_l_mu

where l_mu has grade 2 (spatial bivectors) and e_0 * delta_l_mu has grades 2 and 4
(e_0 * e_i = e_{0i} is grade 2; e_0 * e_{123} = e_{0123} is grade 4).

**Crucially**: L_mu decomposes as:
    L_mu = (l_mu^{23} e_{23} + l_mu^{31} e_{31} + l_mu^{12} e_{12})     [pure quaternion part]
         + (v_mu^1 e_{01} + v_mu^2 e_{02} + v_mu^3 e_{03})              [mixed bivector part]
         + w_mu e_{0123}                                                   [pseudoscalar part]

All components are grade-2 or grade-4 (even), as required for a left current in Cl+(3,0,1).

### 2.4 Important structural observation

The full left current L_mu is valued in the FULL even subalgebra Cl+(3,0,1), which has
dimension 8 (basis: 1, e_{23}, e_{31}, e_{12}, e_{01}, e_{02}, e_{03}, e_{0123}).

However, L_mu has NO scalar (grade-0) component. It lives in the 7-dimensional subspace
spanned by {e_{23}, e_{31}, e_{12}, e_{01}, e_{02}, e_{03}, e_{0123}}.

This is the Lie algebra of the group of unit elements in Cl+(3,0,1) (since L = Psi^{-1} dPsi
is a Maurer-Cartan form, it lives in the tangent space at the identity, i.e., the Lie algebra).

---

## 3. WZW 5-Form Expansion

### 3.1 The WZW 5-form

The WZW term is defined via a 5-form on a 5-manifold M^5 whose boundary is spacetime M^4:

    Gamma_WZW = c_5 integral_{M^5} <L wedge L wedge L wedge L wedge L>_scalar

where L = Psi^{-1} dPsi is the Maurer-Cartan 1-form, the wedge product is over spacetime
indices (antisymmetrized), and <...>_scalar extracts the scalar (grade-0) part of the
Clifford algebra product.

More explicitly, with coordinates x^A on M^5 (A = 0,1,2,3,5 where x^5 is the extension):

    omega_5 = epsilon^{ABCDE} <L_A L_B L_C L_D L_E>_0

where L_A L_B etc. are geometric (Clifford) products of the Lie-algebra-valued 0-forms L_A.

### 3.2 Expansion in epsilon = e_0

With L_mu = l_mu + e_0 delta_l_mu and e_0^2 = 0:

    L_A L_B L_C L_D L_E = (l_A + e_0 delta_l_A)(l_B + e_0 delta_l_B)...(l_E + e_0 delta_l_E)

Since e_0^2 = 0, any term with two or more delta_l factors vanishes (each delta_l comes
multiplied by e_0, and two e_0's somewhere in the product give zero because e_0 X e_0 = 0
for any X in Cl(3,0)).

Wait — this needs more care. The product e_0 A e_0 B where A, B are in Cl(3,0) is not
simply e_0^2 AB. We need: e_0 A e_0 = e_0 (A e_0). Since A is in Cl(3,0) and e_0
anticommutes with each e_i, we get A e_0 = (-1)^{|A|} e_0 A where |A| is the grade of A
(specifically, for A being a product of k spatial basis vectors, e_0 anticommutes past each
one, giving (-1)^k). So e_0 A e_0 = (-1)^{|A|} e_0^2 A = 0.

More carefully: for a monomial A = e_{i_1}...e_{i_k} (all spatial), we have:
    A e_0 = e_{i_1}...e_{i_k} e_0 = (-1)^k e_0 e_{i_1}...e_{i_k} = (-1)^k e_0 A.

So e_0 A e_0 = e_0 [(-1)^k e_0 A] = (-1)^k e_0^2 A = 0. Good.

But what about products like e_0 A l B e_0 C where l has no e_0? Here we have two e_0's
separated by elements of Cl(3,0). The key identity is:

For any X in Cl(3,0): e_0 X e_0 = 0.

Therefore, in the product of five L factors, any two positions occupied by e_0 delta_l
will have the form ...e_0 (stuff) e_0... = 0 where "stuff" is a product of Cl(3,0) elements.

**This confirms**: the expansion of L^5 truncates at first order in e_0.

    L_A L_B L_C L_D L_E = l_A l_B l_C l_D l_E
        + e_0 [delta_l_A * l_B l_C l_D l_E
              + l_A * e_0^{-swap} delta_l_B * l_C l_D l_E + ...]

Actually, let me be more careful about where e_0 sits. When delta_l is at position k,
the product is:

    l_A ... l_{k-1} (e_0 delta_l_k) l_{k+1} ... l_E

The e_0 must be commuted to the left past l_A ... l_{k-1}. Since each l_i is a pure
imaginary quaternion (grade-2 in Cl(3,0)), and e_0 anticommutes with odd-grade elements
but commutes with even-grade elements of Cl(3,0):

    e_0 * (grade-2 element of Cl(3,0)) = ?

For e_0 e_{ij} where i,j are spatial: e_0 e_i e_j = (-1)^2 e_i e_j e_0 = e_{ij} e_0.
So e_0 COMMUTES with spatial bivectors.

Therefore e_0 commutes with all l_i (which are pure imaginary quaternions = spatial bivectors).

This means: l_A ... l_{k-1} (e_0 delta_l_k) l_{k+1} ... l_E
           = e_0 * l_A ... l_{k-1} delta_l_k l_{k+1} ... l_E.

The e_0 factors out to the left completely! So:

    L^5|_{order e_0} = e_0 * sum_{k=1}^{5} l_A...l_{k-1} * delta_l_k * l_{k+1}...l_E

(where the indices ABCDE label the five 1-form directions, and the wedge/antisymmetrization
is implicit).

### 3.3 Scalar extraction

We need <L^5>_0. The zeroth-order term is:

    <l^5>_0 = <l_A l_B l_C l_D l_E>_0

Each l_i is a spatial bivector (grade 2). Product of five grade-2 elements has grades
ranging from 0 to 10, but in Cl(3,0) (which has dimension 8), the maximum grade is 3.
Actually, in Cl(3,0), the generators are e_1, e_2, e_3, so maximum grade is 3. The
product of five bivectors in Cl(3,0) has:

    grade 2 * grade 2 = grades 0 + 2 + 4, but grade 4 doesn't exist in Cl(3,0)
    (highest grade is 3). So grade 2 * grade 2 = grades 0 + 2.

Actually wait. In Cl(3,0), we have grades 0,1,2,3 only (8-dimensional algebra).
e_{23}, e_{31}, e_{12} are all grade 2. Their products:
    e_{23} e_{31} = e_{23} e_{31} = e_2 e_3 e_3 e_1 = e_2 e_1 = -e_{12}    [grade 2]
    e_{23} e_{23} = -1    [grade 0]

So product of two spatial bivectors gives grade 0 + grade 2.
Product of three: (grade 0 + grade 2) * grade 2 = grade 2 + (grade 0 + grade 2) = grade 0 + grade 2.
Product of four: same pattern: grade 0 + grade 2.
Product of five: same: grade 0 + grade 2.

So <l^5>_0 is generically nonzero. This is the standard SU(2) WZW 5-form.

For the first-order term, we need:

    <e_0 * (sum of terms with one delta_l and four l's)>_0

Since e_0 is grade 1, we need the rest to have grade 1 for the total to be grade 0?
No — grade extraction is about the grade of the FULL Cl(3,0,1) element. e_0 times a
grade-k element of Cl(3,0) gives a grade-(k+1) element of Cl(3,0,1) (unless there's
a contraction that reduces grade — but e_0^2 = 0 means there's no contraction with e_0).

Hmm, actually this isn't quite right either. Let me reconsider.

e_0 is grade 1 in Cl(3,0,1). If X is a grade-k element of Cl(3,0), then e_0 X has grade
k+1 in Cl(3,0,1). For the scalar part (grade 0), we need k+1 = 0, i.e., k = -1. Impossible.

Therefore: **<e_0 * anything>_0 = 0 always.**

This means the order-e_0 term of L^5 has ZERO scalar part.

**Result**: <L^5>_0 at order e_0 vanishes identically.

### 3.4 The vanishing is structural

The reason is simple: e_0 has grade 1. Multiplying by e_0 shifts the grade by 1.
The scalar part (grade 0) of e_0 * X requires X to have grade -1, which doesn't exist.

No matter what delta_l is (vector, trivector, or any odd element of Cl(3,0)), the product
e_0 * delta_l has even grade >= 2 in Cl(3,0,1). Combined with four l's (each grade 2 in
Cl(3,0)), the minimum total grade is 2. The scalar part is always zero.

**This kills the standard WZW approach**: the 5-form <L^5>_0 evaluated on the full
Cl+(3,0,1) current does not produce any order-e_0 terms. The degenerate sector decouples
from the standard WZW term.

---

## 4. Alternative Grade Extractions

### 4.1 Grade-4 extraction of L^5

Since <e_0 * X>_0 = 0 but <e_0 * X>_{higher grades} may be nonzero, consider extracting
the pseudoscalar (grade 4) part instead:

    <L^5>_4 at order e_0 = <e_0 * (sum of terms)>_4

e_0 * X has grade k+1 (if X has grade k in Cl(3,0)). For the result to be grade 4 in
Cl(3,0,1), we need k+1 = 4, so k = 3. That means X must be a grade-3 element of Cl(3,0),
i.e., proportional to e_{123}.

So we need: the Cl(3,0)-grade-3 part of (sum of terms with one delta_l and four l's).

Recall delta_l = v^i e_i + w e_{123} (vector + trivector in Cl(3,0)).

The terms with delta_l at position k and four l's at other positions:

    l l l l delta_l + l l l delta_l l + l l delta_l l l + l delta_l l l l + delta_l l l l l

Each l is a spatial bivector. Four l's give grade 0 + grade 2 (as computed above).

Product of (grade 0 + grade 2) with delta_l:
- (grade 0) * (grade 1 + grade 3) = grade 1 + grade 3
- (grade 2) * (grade 1) = grade 1 + grade 3
- (grade 2) * (grade 3) = grade 1 + grade 3    [since 2+3=5 > 3 max, folds back]

Wait, in Cl(3,0), grade 2 * grade 3: e_{23} * e_{123} = e_2 e_3 e_1 e_2 e_3 = e_1 = grade 1.
More generally, bivector * trivector in Cl(3,0) gives grade 1 (only, since |2-3|=1 and
2+3=5 which exceeds max grade 3, so 5 mod ... no, the grade decomposition is:
grade |p-q| and |p-q|+2, ..., up to min(p+q, 2n-p-q). Here n=3, p=2, q=3:
min(5, 6-5)=1. So only grade 1.)

So (grade 2) * (grade 3) = grade 1 only.

Therefore, the products (four l's) * delta_l give:
- From (grade 0)(grade 1): grade 1
- From (grade 0)(grade 3): grade 3
- From (grade 2)(grade 1): grade 1 + grade 3
- From (grade 2)(grade 3): grade 1

The grade-3 contributions come from:
1. The scalar part of (four l's) times the trivector part of delta_l: (l^4)_0 * w e_{123}
2. The bivector part of (four l's) times the vector part of delta_l, selecting grade 3:
   [(l^4)_2 * v^i e_i]_3

These contribute to <L^5>_4 at order e_0:

    <L^5>_4|_{e_0} = e_0 * {grade-3 part of [five terms with delta_l permuted through l^4]}
                   = e_{0123} * {coefficient}

This is proportional to the pseudoscalar e_{0123} and is generically NONZERO.

### 4.2 Explicit evaluation of the grade-3 part

Let us focus on the contribution from the trivector part w e_{123} of delta_l, paired with
the scalar part of (four l's). This is the simplest term.

The scalar part of four l's:

    <l_B l_C l_D l_E>_0

For spatial bivectors l = l^a sigma_a (where sigma_a = {e_{23}, e_{31}, e_{12}}), using the
quaternion algebra sigma_a sigma_b = -delta_{ab} + epsilon_{abc} sigma_c:

    <l_B l_C>_0 = -l_B^a l_C^a = -l_B . l_C    (dot product in R^3)

    <l_B l_C l_D l_E>_0 = <l_B l_C>_0 <l_D l_E>_0 - ... hmm, this isn't right.

Let me use the quaternion trace. For pure imaginary quaternions A, B:
    <AB>_0 = -A . B (the real part of the quaternion product)

For four quaternions:
    <ABCD>_0 = (A.B)(C.D) - (A.C)(B.D) + (A.D)(B.C)
    (This is the standard identity for Tr(sigma sigma sigma sigma) using the quaternion algebra.)

Wait, I need to be more careful. In the quaternion algebra with sigma_a sigma_b = -delta_{ab} + epsilon_{abc} sigma_c:

    AB = -(A.B) + (AxB)    (for pure imaginary quaternions)

where A.B = sum A^a B^a and (AxB)^c = epsilon_{abc} A^a B^b.

    ABCD = [-(A.B) + AxB][-(C.D) + CxD]
         = (A.B)(C.D) - (A.B)(CxD) - (AxB)(C.D) + (AxB)(CxD)

The scalar part: (A.B)(C.D) is scalar. (AxB)(CxD) is a product of two pure imaginary
quaternions, whose scalar part is -(AxB).(CxD).

    <ABCD>_0 = (A.B)(C.D) - (AxB).(CxD)

Using the vector identity (AxB).(CxD) = (A.C)(B.D) - (A.D)(B.C):

    <ABCD>_0 = (A.B)(C.D) - (A.C)(B.D) + (A.D)(B.C)

Good. Now, for the five terms where delta_l is at position k and the other four positions
have l's, the trivector contribution w e_{123} to the grade-3 part is:

At position 1: w_A e_{123} l_B l_C l_D l_E → grade-3 part = w_A <l_B l_C l_D l_E>_0 e_{123}
(since e_{123} * scalar = scalar * e_{123}, and the grade-2 part of l^4 times e_{123}
gives... let me check: (grade 2)(grade 3) = grade 1. So only the scalar part of l^4
survives to give grade 3.)

Actually wait: e_{123} commutes with scalars and with what else? e_{123} is grade 3 in Cl(3,0).
e_{123} sigma_a = e_{123} e_{bc} (for appropriate bc). We computed e_{123} e_{23} = e_1, etc.
So e_{123} * (grade 2) = grade 1. Confirmed: only (l^4)_0 * w e_{123} contributes grade 3.

At position 2: l_A (w_B e_{123}) l_C l_D l_E → need <l_A e_{123}>_{something} that gives
grade 3 when combined with l_C l_D l_E. We have l_A e_{123}: (grade 2)(grade 3) = grade 1.
Then (grade 1) * (l_C l_D l_E): the l's are grade 2, so l^3 has grades 0+2.
(grade 1)(grade 0) = grade 1, (grade 1)(grade 2) = grade 1 + grade 3.
So grade-3 contributions exist. This is more complex.

Hmm, but actually e_{123} COMMUTES with even-grade elements in Cl(3,0). Let me verify:
e_{123} e_{23} = e_1 e_2 e_3 e_2 e_3 = e_1(e_2 e_3 e_2)e_3 = e_1(-e_3)e_3 = -e_1.
But e_{23} e_{123} = e_2 e_3 e_1 e_2 e_3 = e_2(e_3 e_1)e_2 e_3 = e_2(-e_1 e_3)e_2 e_3
= -e_2 e_1 e_3 e_2 e_3 = -e_2 e_1(e_3 e_2)e_3 = -e_2 e_1(-e_2 e_3)e_3 = e_2 e_1 e_2 e_3^2
= e_2 e_1 e_2 = e_2(e_1 e_2) = e_2(-e_2 e_1) = -e_2^2 e_1 = -e_1.

So e_{123} e_{23} = e_{23} e_{123} = -e_1. Yes, e_{123} commutes with spatial bivectors.

More generally, e_{123} commutes with all even elements of Cl(3,0) (since it's the volume
element of an odd-dimensional Clifford algebra... actually Cl(3,0) has dimension 3, and
the volume element e_{123} satisfies: for any a in Cl(3,0), e_{123} a = (-1)^{grade(a)} a e_{123}
... no that's also not right in general.

Let me just check: e_{123} commutes with scalars (yes), commutes with bivectors (verified
above for e_{23}, and by the cyclic symmetry for e_{31}, e_{12}). So e_{123} commutes with
all of Cl^+(3,0) (the even subalgebra), which is the quaternion algebra. Good.

Since e_{123} commutes with all l's (which are pure imaginary quaternions, hence even-grade
in Cl(3,0)), we can move e_{123} freely through any product of l's:

    l_A ... l_{k-1} (w e_{123}) l_{k+1} ... l_E = w * l_A ... l_{k-1} l_{k+1} ... l_E * e_{123}

for all positions k! (Since e_{123} commutes with each l.)

Therefore, the sum over all positions of the trivector part gives:

    sum_{k=1}^{5} l_A...l_{k-1} (w_k e_{123}) l_{k+1}...l_E
    = (sum_{k=1}^{5} w_k * l_A...l_{k-1} l_{k+1}...l_E) * e_{123}

The grade-3 part of each term l_A...l_{k-1} l_{k+1}...l_E * e_{123} is:
<l_A...l_{k-1} l_{k+1}...l_E>_0 * e_{123} (only the scalar part of the four l-product
survives to give grade 3 after multiplication by e_{123}).

So the grade-3 part from the trivector delta_l contribution is:

    [sum_{k} w_k <l_{...other 4...}>_0] * e_{123}

And the corresponding contribution to <L^5>_4 is:

    <L^5>_4|_{e_0, trivector} = e_{0123} * sum_{k} w_k <l_{...other 4...}>_0

### 4.3 What is w_mu?

Recall delta_l_mu = q^{-1} partial_mu d - q^{-1} d l_mu, and delta_l = v^i e_i + w e_{123}.

The trivector part w_mu e_{123} comes from:
- q^{-1} partial_mu d: the e_{123} part of q^{-1}(partial_mu tau e_{123} + partial_mu j_i e_i).
  Since q^{-1} is a quaternion and e_{123} commutes with quaternions:
  q^{-1} (partial_mu tau e_{123}) = (partial_mu tau)(q^{-1}) e_{123}... but this has
  scalar * e_{123} plus bivector * e_{123} = e_{123} + vectors. The e_{123} part is
  just (partial_mu tau)(q^{-1})_0 e_{123} = (partial_mu tau)(s/rho_0^2) e_{123}
  (where (q^{-1})_0 = s/rho_0^2 is the scalar part of q^{-1}).

Wait, I need to be more careful. q^{-1} = q~/rho_0^2 = (s - f_1 e_{23} - f_2 e_{31} - f_3 e_{12})/rho_0^2.

    q^{-1} * (tau e_{123}) = (s/rho_0^2)(tau e_{123}) + (-f_a/rho_0^2 sigma_a)(tau e_{123})

Since sigma_a e_{123} = -e_a (from e_{23} e_{123} = -e_1, etc.):

    = (s tau/rho_0^2) e_{123} + (f_a tau/rho_0^2) e_a

So the e_{123} component of q^{-1}(tau e_{123}) is s tau/rho_0^2.

For q^{-1}(j_i e_i): quaternion times vector. sigma_a e_b = delta_{ab} e_{123} - epsilon_{abc} e_c
(this follows from e_{23} e_1 = e_{23} e_1 = e_2 e_3 e_1 = -e_2 e_1 e_3 = e_{12} e_3...
hmm let me just compute directly.

e_{23} e_1 = e_2 e_3 e_1. To get this to a standard form:
e_3 e_1 = -e_1 e_3, so e_2(-e_1 e_3) = -e_2 e_1 e_3 = e_1 e_2 e_3 = e_{123}. So e_{23} e_1 = e_{123}.

e_{23} e_2 = e_2 e_3 e_2 = -e_2 e_2 e_3 = -e_3.

e_{23} e_3 = e_2 e_3 e_3 = e_2.

So sigma_1 e_1 = e_{123}, sigma_1 e_2 = -e_3, sigma_1 e_3 = e_2.
By cyclic symmetry: sigma_a e_a = e_{123} (no sum), sigma_a e_b = epsilon_{abc} e_c for a != b
... let me check: sigma_1 e_2 = -e_3. epsilon_{123} = +1, so epsilon_{12c} e_c = epsilon_{123} e_3 = e_3.
But we got -e_3. So sigma_a e_b = -epsilon_{abc} e_c for a != b, and sigma_a e_a = e_{123}.

Combined: sigma_a e_b = delta_{ab} e_{123} - epsilon_{abc} e_c.

So q^{-1}(j_b e_b) = (s/rho_0^2)(j_b e_b) + (-f_a/rho_0^2)(sigma_a)(j_b e_b)
= (s j_b/rho_0^2) e_b + (-f_a j_b/rho_0^2)(delta_{ab} e_{123} - epsilon_{abc} e_c)
= (s j_b/rho_0^2) e_b - (f_b j_b/rho_0^2) e_{123} + (f_a j_b epsilon_{abc}/rho_0^2) e_c

The e_{123} component: -(f . j)/rho_0^2 (where f . j = f_1 j_1 + f_2 j_2 + f_3 j_3).

So from q^{-1} partial_mu d, the e_{123} coefficient is:

    (s partial_mu tau - f . partial_mu j) / rho_0^2

Now for q^{-1} d l_mu: we need the e_{123} component of (q^{-1} d) l_mu.

q^{-1} d = (scalar + bivector parts)(vector + trivector parts) -- we partially computed this above.

Actually, let's use phi = q^{-1} d. phi has vector and trivector parts. The trivector part
of phi is (from the above): [(s tau - f.j)/rho_0^2] e_{123}.

Call beta = (s tau - f.j)/rho_0^2, so phi_{trivector} = beta e_{123}.

Now phi l_mu: the trivector part of phi l_mu comes from:
- (vector part of phi) * l_mu: (grade 1)(grade 2) = grade 1 + grade 3. The grade-3 part is...
  e_a * sigma_b = delta_{ab} e_{123} - epsilon_{abc} e_c (from our earlier computation),
  so grade-3 part of (alpha^a e_a)(l^b sigma_b) = alpha^a l^a e_{123} = (alpha . l) e_{123}.
  [where alpha . l = alpha^a l^a_mu, the dot product of the vector parts.]

- (trivector part of phi) * l_mu: beta e_{123} * l^b sigma_b.
  e_{123} sigma_b = -e_b (from e_{123} e_{23} = -e_1 etc., or equivalently since
  e_{123} commutes with sigma_b and sigma_b e_{123} = ... wait, we showed e_{123} sigma_b = sigma_b e_{123}.
  Hmm but sigma_b e_{123}: e_{23} e_{123} = e_2 e_3 e_1 e_2 e_3 = ... = -e_1 (computed earlier).
  And e_{123} e_{23} = e_1 e_2 e_3 e_2 e_3 = e_1(e_2 e_3)(e_2 e_3) = e_1(-1) = -e_1.
  So e_{123} sigma_b = sigma_b e_{123} = -e_b... wait that gives e_{123} sigma_1 = -e_1.
  Check: e_{123} e_{23} = -e_1 yes. And sigma_1 e_{123} = e_{23} e_{123} = -e_1. Same.
  So e_{123} commutes with bivectors (as we established), and e_{123} sigma_b = -e_b.)

  So beta e_{123} l^b sigma_b = beta l^b(-e_b) = -beta l^b e_b. This is grade 1, not grade 3.

So the trivector part of phi l_mu is: (alpha . l_mu) e_{123}, where alpha^a is the vector
part of phi = q^{-1} d.

Therefore: w_mu = trivector coefficient of delta_l_mu
= trivector coeff of (partial_mu phi - phi l_mu)
= partial_mu beta - alpha . l_mu

where beta = (s tau - f.j)/rho_0^2 and alpha = vector part of q^{-1} d.

For the specific case where we are interested in the TAU contribution (setting j = 0 for
simplicity to isolate the tau coupling):

With j = 0: beta = s tau/rho_0^2 and alpha^a = f_a tau/rho_0^2 (from the computation above
of q^{-1}(tau e_{123})).

    w_mu|_{j=0} = partial_mu(s tau/rho_0^2) - (f_a tau/rho_0^2) l_mu^a

### 4.4 The 5-form with grade-4 extraction

From section 4.2, the grade-4 part of L^5 at order e_0, from the trivector contribution, is:

    <L^5>_4|_{e_0, triv} = e_{0123} * epsilon^{ABCDE} sum_{k} w_k <l_{others}>_0

But there's an important point: this is a 5-FORM on a 5-manifold. The epsilon tensor
antisymmetrizes the spacetime indices. Let me write this more carefully.

    <L^5>_4 = e_{0123} * epsilon^{ABCDE} [
        w_A <l_B l_C l_D l_E>_0
      + w_B <l_A l_C l_D l_E>_0
      + w_C <l_A l_B l_D l_E>_0
      + w_D <l_A l_B l_C l_E>_0
      + w_E <l_A l_B l_C l_D>_0
    ]

Using the identity <l_A l_B l_C l_D>_0 = (l_A.l_B)(l_C.l_D) - (l_A.l_C)(l_B.l_D) + (l_A.l_D)(l_B.l_C)
and the antisymmetry of epsilon^{ABCDE}, this can be written as:

    = e_{0123} * 5 * epsilon^{ABCDE} w_A <l_B l_C l_D l_E>_0

(by cyclic symmetry of the 5-form: each of the 5 terms gives the same contribution after
relabeling dummy indices, noting that moving the w from position k to position 1 in the
antisymmetrized sum just relabels indices).

Wait, that's not automatic because the Clifford product is not commutative. But we showed
that e_{123} commutes with all l's, so the position of w*e_{123} in the product doesn't
matter (the e_{123} passes through all l's). The SCALAR part of the l-product depends on
the ordering, but the antisymmetrization over spacetime indices handles that.

Let me reconsider. The five terms with w at each position, after extracting e_{123} to the
right (which commutes through all l's), give:

Position 1: w_A * l_B l_C l_D l_E * e_{123} → grade-3 part is w_A <l_B l_C l_D l_E>_0 e_{123}
Position 2: l_A * w_B * l_C l_D l_E * e_{123} → w_B * l_A l_C l_D l_E * e_{123}
  → grade-3 part is w_B <l_A l_C l_D l_E>_0 e_{123}

Wait, this is not right because w_B is a scalar, and it doesn't affect the Clifford product.
But the point is that w_B at position 2 means the four l's are l_A, l_C, l_D, l_E (skipping B).
Since we antisymmetrize over ABCDE, the five terms have w at each of the five positions,
with the remaining four l's filling the other positions.

Under epsilon^{ABCDE}, swapping A and B gives a sign change. The expression
epsilon^{ABCDE} w_B <l_A l_C l_D l_E>_0 can be relabeled (A<->B):
= epsilon^{BACDE} w_A <l_B l_C l_D l_E>_0 = -epsilon^{ABCDE} w_A <l_B l_C l_D l_E>_0.

Hmm, so the term from position 2 gives -1 times the term from position 1 after relabeling?
That can't be right if they're supposed to add up.

Let me think about this differently. The key issue is that <l_A l_C l_D l_E>_0 is NOT the
same as <l_B l_C l_D l_E>_0 after relabeling (the order of the l's matters for the Clifford
product).

Using the identity for the scalar part of four pure imaginary quaternions:
    <l_A l_B l_C l_D>_0 = (l_A.l_B)(l_C.l_D) - (l_A.l_C)(l_B.l_D) + (l_A.l_D)(l_B.l_C)

This is symmetric under swapping (A,B) with (C,D) and under simultaneous swap of pairs,
but NOT under general permutations. Specifically, under A <-> B:

    <l_B l_A l_C l_D>_0 = (l_B.l_A)(l_C.l_D) - (l_B.l_C)(l_A.l_D) + (l_B.l_D)(l_A.l_C)
                        = (l_A.l_B)(l_C.l_D) - (l_A.l_D)(l_B.l_C) + (l_A.l_C)(l_B.l_D)

Which is NOT equal to <l_A l_B l_C l_D>_0 (the last two terms are swapped). The difference:
    <l_A l_B l_C l_D>_0 - <l_B l_A l_C l_D>_0
    = -2(l_A.l_C)(l_B.l_D) + 2(l_A.l_D)(l_B.l_C)

So the five terms do NOT simply give 5 times the first one.

This is getting complicated. Let me take a different approach.

---

## 5. Connection to the Standard WZW Term

### 5.1 The standard SU(2) WZW

For SU(2)-valued fields (the quaternion sector), the WZW 5-form is:

    omega_5^{SU(2)} = Tr(l wedge l wedge l wedge l wedge l)

where Tr is the quaternion trace (= 2 <...>_0 in our conventions), and l = q^{-1} dq.

In SU(2), the Lie algebra is 3-dimensional (pure imaginary quaternions), and l is valued
in this algebra. The WZW 5-form is nonzero and gives the standard topological term.

When we extend to Cl+(3,0,1), the Lie algebra becomes 7-dimensional:
{e_{23}, e_{31}, e_{12}, e_{01}, e_{02}, e_{03}, e_{0123}}.

The full left current L = l + e_0 delta_l has components in all 7 directions.

### 5.2 Rewriting the grade-4 extraction as a trace

Rather than extracting <L^5>_0 (which kills the degenerate sector) or <L^5>_4 (which
requires the pseudoscalar), consider defining a DIFFERENT trace form.

In the standard WZW, the "trace" is really a non-degenerate ad-invariant bilinear form
on the Lie algebra. For SU(2), this is the Killing form, proportional to Tr(AB) = -2(A.B)
for pure imaginary quaternions.

For Cl+(3,0,1), the full Lie algebra has a DEGENERATE Killing form (because e_0^2 = 0
makes the degenerate sector null). The standard WZW construction requires a non-degenerate
bilinear form, which fails here.

However, there IS a non-degenerate bilinear form on the 7D Lie algebra:

    (A, B) = <A B~>_0       [trace-reversal form]

or alternatively:

    (A, B) = <A hat{B}>_0   [using some involution hat]

The key is: for e_{01} and e_{23}, we have <e_{01} e_{23}~>_0 = <e_{01}(-e_{23})>_0
= -<e_{01} e_{23}>_0 = -<e_{0123}>_0 = 0. Hmm.

Actually: <e_{01} e_{23}>_0: this is the grade-0 part of e_0 e_1 e_2 e_3 = e_{0123},
which is grade 4. So <e_{01} e_{23}>_0 = 0.

The bilinear form <AB>_0 on the Lie algebra pairs:
- e_{23} with e_{23}: <e_{23}^2>_0 = <-1>_0 = -1. Good.
- e_{01} with e_{01}: <e_{01}^2>_0 = <e_0 e_1 e_0 e_1>_0 = <-e_0^2>_0 = 0. Degenerate!

So the natural bilinear form is degenerate on the e_{0i} directions. This is the fundamental
problem: the WZW construction requires a non-degenerate invariant form.

### 5.3 Alternative: use <AB>_4 as a pairing

Consider the grade-4 part:
- <e_{01} e_{23}>_4 = <e_{0123}>_4 = e_{0123} (as a grade-4 element, its "coefficient" is 1).
- <e_{02} e_{31}>_4 = <e_{0123}>_4 = e_{0123}. [Check: e_{02} e_{31} = e_0 e_2 e_3 e_1
  = -e_0 e_2 e_1 e_3 = e_0 e_1 e_2 e_3 = e_{0123}. Grade 4, coefficient 1.]
- <e_{03} e_{12}>_4 = e_{0123}. [e_{03} e_{12} = e_0 e_3 e_1 e_2 = e_0(-e_1 e_3)e_2
  = -e_0 e_1 e_3 e_2 = e_0 e_1 e_2 e_3 = e_{0123}. Yes.]

So grade-4 extraction pairs the bulk bivectors (e_{23}, e_{31}, e_{12}) with the mixed
bivectors (e_{01}, e_{02}, e_{03}) in a non-degenerate way! It's a "mixed" pairing.

What about the other pairings?
- <e_{23} e_{23}>_4 = <-1>_4 = 0. No grade-4 part.
- <e_{01} e_{01}>_4 = <e_0 e_1 e_0 e_1>_4 = <0>_4 = 0 (e_0^2 = 0).
- <e_{0123} * 1>_4 = <e_{0123}>_4 = e_{0123}. (1 is grade 0, so 1 * e_{0123} = e_{0123}.)
  But 1 is not in the Lie algebra (we need the tangent space at identity, which excludes scalars).
  Actually, e_{0123} should pair with... <e_{0123} e_{0123}>_4 = <0>_4 = 0 (since e_{0123}^2 = 0).

So the grade-4 pairing defines a bilinear form on the Lie algebra:

    B(X, Y) = coefficient of e_{0123} in X*Y

This form is:
- Non-degenerate on the 6D subspace {e_{23}, e_{31}, e_{12}, e_{01}, e_{02}, e_{03}}
  (it pairs spatial bivectors with mixed bivectors).
- Degenerate on e_{0123} (pairs with nothing in the Lie algebra, since it would need
  to pair with the scalar 1, which is not in the Lie algebra).
- Also degenerate: e_{0123} doesn't pair with anything.

So this pairing is also degenerate (rank 6 out of 7).

### 5.4 The irreducible structure

The 7D Lie algebra of Cl+(3,0,1) decomposes as:

    g = su(2) + V + R

where su(2) = span{e_{23}, e_{31}, e_{12}} (spatial bivectors = quaternion imaginaries),
V = span{e_{01}, e_{02}, e_{03}} (mixed bivectors = degenerate vectors),
R = span{e_{0123}} (pseudoscalar).

The Lie brackets:
- [su(2), su(2)] = su(2) (quaternion commutators: [e_{23}, e_{31}] = 2e_{12} etc.)
- [su(2), V]: [e_{23}, e_{01}] = e_{23} e_{01} - e_{01} e_{23}
  = e_{0123} - e_0 e_1 e_2 e_3... wait, we need to compute this.

    e_{23} e_{01} = e_{0123} (computed in section 1).
    e_{01} e_{23} = e_0 e_1 e_2 e_3 = e_{0123} (same computation, just different order of
    multiplying... let me verify: e_0 e_1 e_2 e_3 = e_{0123}. And e_2 e_3 e_0 e_1:
    we showed this also equals e_{0123}. So they're the same!)

    [e_{23}, e_{01}] = e_{0123} - e_{0123} = 0.

Hmm! The commutator vanishes? Let me recheck.

e_{23} e_{01} = e_2 e_3 e_0 e_1 = e_0 e_2 e_3 e_1 (moving e_0 left, two anticommutations)
= e_0(-e_2 e_1 e_3) = -e_0 e_2 e_1 e_3 = e_0 e_1 e_2 e_3 = e_{0123}. (***)

e_{01} e_{23} = e_0 e_1 e_2 e_3 = e_{0123}. (***)

So indeed [e_{23}, e_{01}] = 0. But what about [e_{23}, e_{02}]?

e_{23} e_{02} = e_2 e_3 e_0 e_2 = e_0 e_2 e_3 e_2 = e_0 e_2(-e_2 e_3) = -e_0 e_3 = -e_{03}.

e_{02} e_{23} = e_0 e_2 e_2 e_3 = e_0 e_3 = e_{03}.

[e_{23}, e_{02}] = -e_{03} - e_{03} = -2e_{03}.

Good, so the bracket is nontrivial for non-matching indices.

Similarly: [e_{23}, e_{03}] = ?
e_{23} e_{03} = e_2 e_3 e_0 e_3 = e_0 e_2 e_3 e_3 = e_0 e_2 = e_{02}.
e_{03} e_{23} = e_0 e_3 e_2 e_3 = -e_0 e_2 e_3 e_3 = -e_0 e_2 = -e_{02}.
[e_{23}, e_{03}] = e_{02} - (-e_{02}) = 2e_{02}.

So [su(2), V] = V. The action is the adjoint: su(2) acts on V as the vector representation.

- [V, V]: [e_{01}, e_{02}] = e_{01} e_{02} - e_{02} e_{01}
  = e_0 e_1 e_0 e_2 - e_0 e_2 e_0 e_1
  = -e_0 e_0 e_1 e_2 - (-e_0 e_0 e_2 e_1) = 0 - 0 = 0. (Since e_0^2 = 0.)

So [V, V] = 0. V is an abelian ideal.

- [su(2), R]: [e_{23}, e_{0123}] = 0 (since e_{0123} commutes with all even elements).
  So [su(2), R] = 0.

- [V, R]: [e_{01}, e_{0123}] = e_{01} e_{0123} - e_{0123} e_{01}.
  e_{0123} commutes with all even elements, and e_{01} is even (grade 2). So the commutator is 0.

- [R, R] = 0 (e_{0123} with itself: e_{0123}^2 = 0, but [e_{0123}, e_{0123}] = 0 trivially).

Summary of Lie algebra structure:

    g = su(2) semidirect (V + R)

where V + R is an abelian ideal, su(2) acts on V by the vector representation, and R is
central (commutes with everything).

This is NOT a semisimple Lie algebra. The Killing form is degenerate on V + R.
The standard WZW construction requires a semisimple Lie algebra (or at least a non-degenerate
invariant bilinear form). This is why the standard approach fails.

---

## 6. The WZW Term with Mixed Pairing

### 6.1 A non-standard approach

Although the natural Killing form on g is degenerate, we identified in section 5.3 a
non-degenerate bilinear form on a 6D subspace using grade-4 extraction. This pairs
bulk bivectors with mixed bivectors:

    B(e_{23}, e_{01}) = 1, B(e_{31}, e_{02}) = 1, B(e_{12}, e_{03}) = 1

with all other pairings zero (within the 6D subspace).

Define a "mixed WZW" 5-form using this pairing:

    omega_5^{mixed} = B(L, [L, [L, [L, L]]])

or more precisely, using the coefficient of e_{0123} in L^5:

    omega_5^{mixed} = coeff_{e_{0123}} in L_A L_B L_C L_D L_E (antisymmetrized)

From section 3, this is exactly <L^5>_4 (up to a factor of e_{0123}).

At zeroth order in e_0: <l^5>_4. Since l is a pure imaginary quaternion (spatial bivector),
l^5 has grades 0 and 2 only (in Cl(3,0)). In Cl(3,0,1), these are also grades 0 and 2.
No grade-4 part. So <l^5>_4 = 0.

At first order in e_0: this is what we computed in section 4.2 — it's generically nonzero!

**This is the key result**: the standard WZW 5-form <L^5>_0 does NOT see the degenerate
sector (grade-0 extraction kills it), but the MIXED WZW 5-form <L^5>_4 DOES see it,
and it's the ONLY order that contributes (the zeroth order in e_0 vanishes for <...>_4).

### 6.2 Explicit computation of <L^5>_4

From section 4.2, we need the grade-3 (trivector e_{123}) part of the product of five
Lie-algebra elements, four of which are l's and one is delta_l.

But we also need the contribution from the VECTOR part (v^i e_i) of delta_l, not just
the trivector part. Let me reconsider.

The full delta_l contributes to <L^5>_4 through grade-3 parts of the 5-fold product.
After factoring out e_0 (which commutes with all l's), we need:

    <e_0 * [sum over positions of delta_l in l^4]>_4 = e_0 * [grade-3 part of ...]

e_0 * (grade-3 element of Cl(3,0)) = e_0 * (c e_{123}) = c e_{0123}, which IS grade 4. Good.
e_0 * (grade-1 element of Cl(3,0)) = e_0 * (c^i e_i) = c^i e_{0i}, which is grade 2. NOT grade 4.

So the grade-4 part of e_0 * delta_l ... l^4 comes ONLY from the grade-3 (trivector)
parts of the products. And the trivector part comes from two sources:
(a) The trivector part (w e_{123}) of delta_l, combined with the scalar part of l^4.
(b) The vector part (v^i e_i) of delta_l, combined with bivector parts of l^4 that produce
    a trivector when multiplied by a vector.

Let me compute both contributions.

**Contribution (a)**: from section 4.2,

    sum_k w_k <l_{other 4}>_0

This involves the trivector part of delta_l (which is related to tau) paired with the
scalar trace of four l's.

**Contribution (b)**: The vector part of delta_l is v^i e_i. In the five-fold product with
delta_l at position k, we need the trivector part of (l^{k-1} pieces)(v^i e_i)(l^{4-k+1} pieces).

Since e_i * sigma_a = delta_{ia} e_{123} - epsilon_{iac} e_c (computed earlier), and
sigma_a * e_i = delta_{ai} e_{123} + epsilon_{aic} e_c (by the reversal of the cross product
term), the vector-bivector products generate both trivector and vector parts.

The full computation is intricate but systematic. Let me organize it.

For delta_l at position 1 (vector part only):
    v_A^i e_i * l_B l_C l_D l_E

We need the grade-3 part. Using e_i (l_B l_C l_D l_E):
l_B l_C l_D l_E has grade 0 + grade 2 parts. So:
e_i * (grade 0) = grade 1. Not grade 3.
e_i * (grade 2) = grade 1 + grade 3.

The grade-3 part of e_i * (grade 2 element): for e_i * (l^{ab} sigma_a ... ) this gives
e_i sigma_a |_{grade 3} = delta_{ia} e_{123}. So the grade-3 part of e_i * X_2 (where X_2
is the grade-2 part of l_B l_C l_D l_E) is: X_2^{biv,i} e_{123}, where X_2^{biv,i} is the
coefficient of sigma_i in X_2. That is:

    [e_i * <l_B l_C l_D l_E>_2]_{grade 3} = <l_B l_C l_D l_E>_2^i * e_{123}

where <l_B l_C l_D l_E>_2^i is the sigma_i component of the bivector part of the product.

For delta_l at position 2:
    l_A (v_B^i e_i) l_C l_D l_E

e_i is sandwiched between l's. We need sigma_a e_i sigma_b... the products get complicated.

Rather than computing all this in full generality (which would be very lengthy), let me
focus on the KEY QUESTION: does the result contain a term of the form B^0 * tau?

### 6.3 Reducing to 4D: the WZW as a boundary term

The standard trick for the WZW term is to reduce the 5D integral to 4D. The 5-form
omega_5 is closed (d omega_5 = 0 by the Maurer-Cartan equation dL + L wedge L = 0).
Therefore, by Stokes' theorem:

    integral_{M^5} omega_5 = integral_{M^4 = boundary} omega_4

where omega_4 is a 4-form such that d omega_4 = omega_5.

For the STANDARD SU(2) WZW with <L^5>_0, the 4D reduction gives the Wess-Zumino
term involving the baryon current. What does the MIXED WZW with <L^5>_4 give?

### 6.4 The 4D reduction for the mixed WZW

The 5-form (at order e_0) is:

    omega_5|_{e_0} = e_{0123} * 5 * epsilon^{ABCDE} w_A <l_B l_C l_D l_E>_0 + (vector contributions)

Wait, I was too hasty in section 4.4 claiming the cyclic symmetry gives a factor of 5.
Let me redo this carefully.

Actually, for the trivector contribution (a) only, where w e_{123} commutes with all l's:

The five terms are:
    Term 1: (w_A e_{123}) l_B l_C l_D l_E → w_A (l_B l_C l_D l_E e_{123})_{grade 3}
                                            = w_A <l_B l_C l_D l_E>_0 e_{123}

    Term 2: l_A (w_B e_{123}) l_C l_D l_E = w_B l_A e_{123} l_C l_D l_E
           = w_B e_{123} l_A l_C l_D l_E   (e_{123} commutes with l)
           → grade-3: w_B <l_A l_C l_D l_E>_0 e_{123}

    Term 3: l_A l_B (w_C e_{123}) l_D l_E = w_C e_{123} l_A l_B l_D l_E
           → grade-3: w_C <l_A l_B l_D l_E>_0 e_{123}

    Term 4: w_D <l_A l_B l_C l_E>_0 e_{123}

    Term 5: w_E <l_A l_B l_C l_D>_0 e_{123}

After antisymmetrization with epsilon^{ABCDE}:

    omega_5^{(a)}|_{e_0} = e_{0123} * epsilon^{ABCDE} [w_A T_{BCDE} + w_B T_{ACDE} + w_C T_{ABDE} + w_D T_{ABCE} + w_E T_{ABCD}]

where T_{BCDE} = <l_B l_C l_D l_E>_0.

Now, T_{BCDE} is NOT fully antisymmetric in its indices (the Clifford product is ordered).
However, epsilon^{ABCDE} forces full antisymmetry. Define:

    T_{BCDE}^{antisym} = (1/4!) sum_{perms of BCDE} sign(perm) * T_{perm(BCDE)}

Then: epsilon^{ABCDE} w_A T_{BCDE} = w_A * epsilon^{ABCDE} T_{BCDE}
= w_A * 4! * T_{[BCDE]}^{antisym} * (sign factor).

Actually, since epsilon^{ABCDE} already selects one ordering (all distinct), we need:

    sum_{ABCDE distinct} epsilon^{ABCDE} w_A T_{BCDE}
    = sum_A w_A * sum_{BCDE} epsilon^{ABCDE} T_{BCDE}

Define S_A = sum_{BCDE} epsilon^{ABCDE} T_{BCDE} = sum_{BCDE} epsilon^{ABCDE} <l_B l_C l_D l_E>_0.

For the SU(2) Skyrmion, we can use the quaternion trace identity. The fully antisymmetrized
trace of four pure imaginary quaternions:

    sum_{BCDE} epsilon^{ABCDE} <l_B l_C l_D l_E>_0

This involves the Levi-Civita symbol in 5D contracted with the four-fold trace. In the
standard WZW reduction from 5D to 4D, the extension parameter (say x^5) appears in the
l's through the interpolation l(x, t) where t in [0,1] extends spacetime.

At this point, rather than pursuing the most general computation, let me specialize to
the physically relevant case and use known results from Skyrme theory.

---

## 7. Using Known WZW Reduction Results

### 7.1 Standard SU(2) WZW in Skyrme theory

The standard WZW 5-form for SU(2) (using <L^5>_0 with quaternion trace) reduces to 4D as:

    Gamma_WZW = (N_c / 240 pi^2) integral_{M^5} <l^5>_0
              = (N_c / 240 pi^2) * 240 pi^2 * integral_{M^4} B^mu d sigma_mu    [up to normalization]

Actually, the precise 4D reduction is:

    Gamma_WZW|_{4D} = N_c * integral B^0 d^4x * (coupling to external gauge field)

The WZW term by itself (without external gauge fields) is a topological invariant and
does not contribute to the equations of motion. It becomes dynamically relevant when
there is an external gauge field (like the omega meson) that couples to the baryon current.

In the Adkins-Nappi model and the standard Skyrme approach with vector mesons, the
WZW term gives the omega-baryon coupling:

    L_{WZW} superset (N_c / 48 pi^2) epsilon^{mu nu rho sigma} Tr(l_mu l_nu l_rho) omega_sigma
           = N_c B^mu omega_mu

where omega_mu is the omega meson field and B^mu is the baryon current density.

### 7.2 Mapping to Cl+(3,0,1)

In our framework, the degenerate sector plays the role of external gauge fields:
- tau (the e_{0123} coefficient) maps to omega_0 (temporal component of omega meson)
- j_i (the e_{0i} coefficients) map to omega_i (spatial components)

The full left current L = l + e_0 delta_l, and the order-e_0 part of <L^5>_4 is precisely
the cross-term between the SU(2) currents (l) and the degenerate fields (via delta_l).

From the known WZW structure, the 4D reduction should give:

    <L^5>_4|_{4D, order e_0} ~ epsilon^{mu nu rho sigma} <l_mu l_nu l_rho>_0 * w_sigma * e_{0123}
                               + (contributions from vector part of delta_l)
                               + (contributions from non-trivector positions)

The first term involves:
    <l_mu l_nu l_rho>_0 = scalar part of product of three pure imaginary quaternions
    = epsilon_{abc} l_mu^a l_nu^b l_rho^c (since <sigma_a sigma_b sigma_c>_0 = -epsilon_{abc})

And epsilon^{mu nu rho sigma} epsilon_{abc} l_mu^a l_nu^b l_rho^c = 6 * det(l) type expression
= related to the baryon current B^sigma.

Actually: B^mu = -(1/24 pi^2) epsilon^{mu nu rho sigma} epsilon_{abc} l_nu^a l_rho^b l_sigma^c.

So epsilon^{mu nu rho sigma} <l_nu l_rho l_sigma>_0 = -epsilon^{mu nu rho sigma} epsilon_{abc} l_nu^a l_rho^b l_sigma^c
= 24 pi^2 B^mu.

Wait, let me be careful about the numerical factors.

The baryon current is: B^mu = (1/24 pi^2) epsilon^{mu nu rho sigma} Tr(l_nu l_rho l_sigma)
where Tr for SU(2) matrices means: Tr(sigma_a sigma_b sigma_c) = 2i epsilon_{abc}
(using Pauli matrices). In our quaternion convention:
<sigma_a sigma_b sigma_c>_0 = -epsilon_{abc} (pure imaginary quaternion triple product).

So Tr(l_nu l_rho l_sigma) in the physics convention = -2 <l_nu l_rho l_sigma>_0 * i
... actually the mapping between quaternions and SU(2) matrices is:
q = s + f_a sigma_a ↔ U = s I + i f_a tau_a (where tau_a are Pauli matrices / 2i's).

The precise factor isn't critical for the structural question. What matters is:

**epsilon^{mu nu rho sigma} <l_nu l_rho l_sigma>_0 is proportional to B^mu.**

### 7.3 The B^0 tau coupling from the mixed WZW

Combining everything: the 4D reduction of the order-e_0 part of <L^5>_4 contains a term:

    ~ epsilon^{mu nu rho sigma} <l_nu l_rho l_sigma>_0 * w_mu * e_{0123}
    ~ B^mu * w_mu * e_{0123}

Now recall w_mu is the trivector coefficient of delta_l_mu, which (from section 4.3 with j=0):

    w_mu = partial_mu(s tau / rho_0^2) - (f . l_mu) tau / rho_0^2

The first term involves partial_mu tau (a derivative of the pseudoscalar field).
The second involves tau times the dot product of f with the left current.

On the sigma model (s^2 + |f|^2 = rho_0^2), the factor s/rho_0^2 = cos(f)/rho_0 for the
hedgehog. In the 4D WZW term, the coupling:

    B^mu * w_mu ~ B^0 * w_0 + B^i * w_i

The temporal component B^0 is the baryon density, and w_0 = partial_0(s tau/rho_0^2) - ...

For a STATIC soliton, B^0 = B^0(x) (time-independent) and the coupling to the static
part of w gives:

    B^0 * w_0|_{static} → 0 (time derivatives vanish for static fields)

But the SPATIAL terms also contribute. For a static soliton, B^i = 0 (no baryon current
flow), but the spatial w_i:

    w_i = partial_i(s tau/rho_0^2) - (f_a/rho_0^2) tau l_i^a

This is nonzero for spatially varying tau.

### 7.4 The equation of motion for tau

The action contains (schematically):

    S_{mixed WZW} ~ integral <L^5>_4 d^5x → integral d^4x [B^mu w_mu + ...]

Varying with respect to tau (remembering w_mu depends on tau through
w_mu = partial_mu(s tau/rho_0^2) - ...):

The Euler-Lagrange equation involves:

    partial_mu (delta(B^nu w_nu) / delta(partial_mu tau)) - delta(B^nu w_nu) / delta(tau) = 0

Since w_mu contains partial_mu tau (through the first term partial_mu(s tau/rho_0^2)):

    B^mu * (s/rho_0^2) * partial_mu appears as a first-order operator, NOT a second-order
    (Laplacian) operator.

Actually wait. The WZW term is integrated over 5D and then reduced to 4D. In the 4D form,
it typically appears as a first-derivative term (it's the Chern-Simons form), not a
second-derivative term. The equation of motion from varying a CS term gives a first-order
constraint, not a Poisson equation.

This is an important subtlety. The WZW term is TOPOLOGICAL (does not involve a metric)
and gives first-order equations, not second-order. The coupling to tau would be:

    B^mu * partial_mu tau ~ B^0 dot{tau} + B^i partial_i tau

For a STATIC configuration with B^i = 0, the spatial equation is:

    0 = ... (the WZW contribution to the tau equation is purely first-order)

This means the WZW term alone does NOT give a Poisson equation for tau. It gives a
FIRST-ORDER constraint.

However, if there is ALSO a kinetic term for tau (from the grade-4 part of the quadratic
Lagrangian, which we computed as 2 partial_mu s partial_mu tau + ...), then the COMBINED
equation of motion is:

    kinetic terms (2nd order in tau) + WZW terms (1st order in tau) = 0

On the background of a static Skyrmion, the 2nd-order kinetic operator is an elliptic
Laplacian-like operator (nabla^2 tau + ...) and the WZW term provides the SOURCE:

    nabla^2 tau + ... = C * B^0 + ...

This IS the structure we want!

---

## 8. Detailed Coefficient Computation

### 8.1 Setting up the problem on the hedgehog

On the hedgehog ansatz: q = rho_0(cos f(r) + sin f(r) hat{r} . vec{sigma}),
with f(0) = pi, f(infinity) = 0.

Left current: l_i = q^{-1} partial_i q. In the hedgehog:
    l_i = -(f'/rho_0) hat{r}_i (hat{r} . sigma) - (sin f / (rho_0 r)) [hat{r} x hat{e}_i] . sigma
        ... actually in our conventions l_i has specific radial and angular parts.

The standard expressions:
    l_r = -f' sigma_r     (where sigma_r = hat{r} . vec{sigma})
    l_theta = -sin(f)/r * sigma_theta
    l_phi = -sin(f)/r * sigma_phi

(with sigma_r = sigma_a hat{r}^a, etc.)

The baryon density:
    B^0 = -(1/2pi^2) epsilon_{ijk} <l_i l_j l_k>_0
        = -(1/2pi^2) * (-f' sin^2 f / r^2) * (angular factor)
        = -f' sin^2 f / (2 pi^2 r^2)

(positive for the Skyrmion since f' < 0).

### 8.2 The full 5-form computation (schematic)

Rather than computing the 5D integral in full (which requires specifying the interpolation
to 5D and is very lengthy), I'll use the STRUCTURE of the result.

The mixed WZW 5-form <L^5>_4 at order e_0 has the structure:

    omega_5^{mixed} = c * omega_3^{SU(2)} wedge omega_2^{degen}

where omega_3^{SU(2)} ~ <l wedge l wedge l>_0 is the SU(2) Cartan 3-form (related to B^0)
and omega_2^{degen} involves the degenerate sector components (tau, j_i).

In the 4D reduction, this gives:

    Gamma_4^{mixed} = c' * integral B^mu omega_mu^{degen} d^4x

where omega_mu^{degen} is a 1-form constructed from tau and j_i.

### 8.3 Identifying the coefficient

The WZW coefficient for SU(2) is quantized by pi_5(SU(2)) = Z_2 (for SU(2))
or pi_5(SU(N)) = Z for N >= 3. For SU(2), the normalization is:

    Gamma_{WZW} = (2 pi) * n * integral omega_5 / Vol(S^5)

where n is an integer (physically identified with N_c for the Skyrme model via anomaly
matching).

For Cl+(3,0,1), the relevant homotopy group would be pi_5 of the group manifold.
The group of unit elements in Cl+(3,0,1) with the dual-number structure (q + e_0 d with
|q| = rho_0 and the constraint on d) has the topology of SU(2) x R^4 (or a fiber bundle
over S^3 with R^4 fibers). The homotopy:

    pi_5(SU(2) x R^4) = pi_5(S^3) = Z_2

So the WZW term IS topologically quantized, and the coefficient can take values 0 or 1
(mod 2).

### 8.4 The coupling strength

If the mixed WZW coefficient is c_WZW (= 0 or 1 mod 2 for SU(2), or = N_c for SU(N)),
the 4D Lagrangian density is:

    L_{WZW}^{mixed} = c_WZW * B^mu * omega_mu

where omega_mu involves tau and j_i. The precise form of omega_mu determines the coupling.

For the temporal component (static case):

    L_{WZW}|_{static} = c_WZW * B^0 * (something involving tau)

If omega_0 = tau (the simplest case), then:

    L_{WZW}|_{static} = c_WZW * B^0 * tau

Combined with the kinetic term for tau from the grade-4 quadratic Lagrangian:

    L_{kin} = 2 partial_mu s * partial_mu tau (from <(dPsi)(dPsi~)>_4)

The equation of motion for tau:

    2 nabla^2 s - c_WZW * B^0 = 0    (varying with respect to tau)

Wait, the kinetic term 2(partial_mu s)(partial_mu tau) = 2 nabla (s . nabla tau) - 2(nabla^2 s) tau + ...
by integration by parts. The EL equation from varying tau:

    delta S / delta tau = -2 nabla^2 s + c_WZW B^0 = 0

This gives:

    nabla^2 s = (c_WZW / 2) B^0

BUT s = rho_0 cos f(r), and nabla^2 s has zero monopole (as discussed in the hypothesis
document). So this equation is INCONSISTENT unless c_WZW = 0.

Hmm. This seems to contradict the existence of the coupling. But wait — the equation above
is for varying tau in the ACTION, treating tau as an independent field. The kinetic term
creates a direct coupling between s and tau that doesn't involve B^0 as an intermediary.

Let me reconsider. The total Lagrangian is:

    L = L_bulk(q) + L_{kin}^{cross}(q, tau) + L_{WZW}^{mixed}(q, tau) + ...

    L_bulk = |dq|^2 + Skyrme terms (independent of tau)
    L_{kin}^{cross} = 2(partial_mu s)(partial_mu tau) - 2(partial_mu f_a)(partial_mu j_a)
    L_{WZW}^{mixed} = c_WZW B^mu [tau-dependent part]

The EL equation for tau (where tau appears in L_kin and L_WZW):

From L_kin: integration by parts of 2(partial_mu s)(partial_mu tau):
    EL contribution = -2 nabla^2 s (= -2 Box s in Lorentzian)

From L_WZW: if L_WZW = c B^0 tau, the EL contribution is c B^0.

Total: -2 nabla^2 s + c B^0 = 0, i.e., 2 nabla^2 s = c B^0.

But integral nabla^2 s dV = 0 and integral B^0 dV = 1. This gives 0 = c,
which would force c = 0.

This is the MONOPOLE MISMATCH: nabla^2 s has zero monopole, B^0 has unit monopole.
The two cannot be equated for any nonzero c.

### 8.5 Resolution: the equation is for tau, not s

Wait. I made an error in the EL derivation. Let me redo it.

The kinetic cross-term is L_{kin} = 2(partial_mu s)(partial_mu tau). Here s = s(q)
depends on the bulk field q, and tau is the independent degenerate field.

Varying with respect to tau:
    delta L_{kin} / delta tau = -2 partial_mu partial_mu s = -2 Box s

(by integration by parts). This is a FIXED source (determined by the background q).

So the EL for tau is:

    -2 Box s + c_WZW * (d/d tau)[B^mu omega_mu] = 0

If L_WZW = c B^0 tau, then d(B^0 tau)/d tau = B^0, and:

    -2 nabla^2 s + c B^0 = 0     (static case)

This gives: tau equation is -2 nabla^2 s + c B^0 = 0. But this isn't an equation FOR tau —
tau doesn't appear! This means tau is UNDETERMINED by this equation (it doesn't appear in
its own EL equation).

The reason: the quadratic kinetic term is CROSS-coupled (s and tau), so varying tau gives
an equation involving s but not tau. And the WZW term, if linear in tau, gives a constant
(B^0) with no tau dependence.

This means the COMBINATION of the kinetic cross-term and the WZW term gives a CONSTRAINT
on the BACKGROUND field q (specifically on nabla^2 s), not an equation for tau.

For the system to determine tau, we need a term that is QUADRATIC in tau (or higher).
The quadratic kinetic term only has tau linearly (cross-coupled with s). A pure tau
kinetic term would need to come from higher-order terms or from a different mechanism.

### 8.6 What this means physically

The grade-4 Lagrangian L = 2(partial_mu s)(partial_mu tau) - 2(partial_mu f_a)(partial_mu j_a)
is a CROSS-COUPLING. It's not a standalone kinetic term for tau; it's a constraint that
relates changes in tau to changes in s (and j to f).

The WZW term provides an additional constraint. Together, they overdetermine the system
(or provide consistency conditions) rather than giving a propagation equation for tau.

This is actually consistent with the philosophy that the degenerate sector is CONSTRAINED
(not dynamical). The constraint comes from:

1. The sigma-model condition: tau = j_r tan f (one equation for two unknowns)
2. The kinetic cross-coupling EL: -2 nabla^2 s + c B^0 = 0 (consistency condition on q)
3. Or: the EL for j_a gives: +2 nabla^2 f_a + ... = 0 (further conditions)

The system of constraints may UNIQUELY determine tau and j_i on a given Skyrmion background,
with the WZW term fixing the overall normalization.

---

## 9. The Correct Framework: WZW as Gauged WZW

### 9.1 The gauged WZW perspective

In standard Skyrme theory, the omega meson is introduced as a gauge field, and the
coupling to the baryon current comes from the GAUGED WZW term. The omega is NOT part of
the SU(2) field — it's an external U(1) gauge field.

In Cl+(3,0,1), the degenerate sector (tau, j_i) is analogous to an external gauge field
"living inside" the algebra. The correct framework is:

    Psi = q * (1 + e_0 phi)    where phi = q^{-1} d

Here q in SU(2) is the "matter field" and phi is the "gauge-like field" (constrained by e_0^2 = 0).

The WZW term for q, GAUGED by the degenerate sector, gives:

    Gamma_{gWZW} = Gamma_{WZW}(q) + integral j_{WZW}^mu * A_mu

where A_mu is the "gauge field" constructed from the degenerate sector and j_{WZW}^mu is
the topological current.

### 9.2 Identifying the gauge field

The degenerate sector provides a natural U(1) gauge field:

    A_mu = w_mu = trivector coefficient of delta_l_mu = coeff of e_{123} in q^{-1} D_mu d

For the temporal component on a static background:
    A_0 = (s dot{tau} - f . dot{j}) / rho_0^2

For the spatial components:
    A_i = (s partial_i tau - f . partial_i j) / rho_0^2 - (alpha . l_i)

where alpha is the vector part of phi = q^{-1} d.

The gauged WZW coupling is then:

    L_{gWZW} = N_c * B^mu * A_mu = N_c * B^mu * w_mu

This gives the equation (varying with respect to tau while keeping q fixed):

    N_c * B^mu * (s / rho_0^2) * g_{mu 0} + [kinetic terms] = 0

For a static configuration:
    -2 nabla^2 s + N_c (s / rho_0^2) B^0 = ... (other terms)

But this still has the monopole mismatch issue from section 8.4.

### 9.3 Resolution: separating the constraint equation

The correct approach is to recognize that the degenerate sector has TWO sets of equations:

1. **The sigma-model constraint** (algebraic): relates tau and j_i to q
2. **The WZW/kinetic equation** (differential): determines the combination of tau, j_i
   that is NOT fixed by the sigma-model constraint

The sigma-model constraint tau = j_r tan f (on the hedgehog) leaves ONE free function
(say j_r(r)). The WZW equation determines this function.

From the gauged WZW:

    N_c B^mu w_mu = N_c B^0 w_0 + N_c B^i w_i

For a static hedgehog (B^i = 0):

    N_c B^0 * w_0 = 0    (since everything is static, w_0 = 0)

For the spatial part, w_i is nonzero because tau and j_r vary in space.
But B^i = 0, so N_c B^i w_i = 0 too.

This gives a TRIVIAL equation: 0 = 0. The WZW term does not constrain the degenerate
sector on a static background!

This makes sense: the WZW term is topological and parity-violating. For a parity-symmetric
static solution, it contributes nothing.

### 9.4 When the WZW term is nontrivial

The WZW term becomes nontrivial in two cases:

1. **Rotating Skyrmion**: For a spinning Skyrmion (with angular velocity Omega), B^i is
   nonzero and couples to the spatial gauge field. This gives the omega-meson coupling
   that determines the proton-neutron mass splitting in standard Skyrme theory.

2. **External fields**: When an external omega-meson field is present, it couples to B^0
   and generates a potential energy ~ N_c B omega_0.

In our framework, case 1 means the degenerate sector IS excited by rotation but NOT by
the static soliton. The tau field is zero for a static hedgehog (at least from the WZW
contribution).

---

## 10. The Vector Part Contribution

### 10.1 Contribution (b) from section 6.2

We have not yet fully analyzed the contribution from the VECTOR part (v^i e_i) of delta_l
to the mixed WZW. This contribution involves products like:

    (v^i e_i)(l_B l_C l_D l_E)_{grade 3}

and requires the grade-3 part of mixed products of vectors and bivectors. Since e_i sigma_a
produces both vectors and trivectors, these terms generate additional couplings.

However, the vector part of delta_l (the v^i components) corresponds to the j_i components
of the degenerate sector, not tau. The j_i fields are the "spatial gauge field" components,
analogous to the spatial omega_i in standard Skyrme theory.

The coupling structure from these terms would be:

    ~ epsilon^{ABCDE} v_A^i (something involving l's)_i,BCDE

which, after 4D reduction, gives terms coupling j_i to the spatial baryon current B^i.
For a static hedgehog, B^i = 0, so these terms also vanish.

### 10.2 Summary of all contributions

| Contribution | Source | Coupling structure | Static hedgehog |
|-------------|--------|-------------------|-----------------|
| Trivector (a) | tau via w_mu | B^mu w_mu | 0 = 0 (trivial) |
| Vector (b) | j_i via v_mu^i | B^mu v_mu^i (contracted) | B^i = 0 (vanishes) |
| Cross-kinetic | <(dPsi)(dPsi~)>_4 | nabla^2 s sourcing tau | Zero monopole |

None of these produce a B^0 tau monopole coupling on a static hedgehog background.

---

## 11. Conclusions

### 11.1 Main Result

**The WZW 5-form on the full Cl+(3,0,1) left current does NOT produce a B^0 tau monopole coupling on a static Skyrmion background.**

The reasons are threefold:

1. **Grade obstruction for <L^5>_0**: The standard scalar extraction of the WZW 5-form
   gives ZERO at order e_0 because e_0 is grade 1 and cannot contribute to the scalar part.
   The degenerate sector completely decouples from <L^5>_0.

2. **The mixed WZW <L^5>_4 exists but is trivial on static backgrounds**: The grade-4
   extraction does produce nonzero order-e_0 terms, but they couple to B^mu w_mu (baryon
   current contracted with the degenerate gauge field). For a static, parity-symmetric
   hedgehog, both B^i = 0 and w_0 = 0, so the coupling vanishes.

3. **The kinetic cross-coupling has zero monopole**: The quadratic Lagrangian's grade-4
   part gives 2(partial_mu s)(partial_mu tau), which sources tau through nabla^2 s.
   This has zero total monopole (integral nabla^2 s dV = 0), producing only dipole
   and higher multipole fields (1/r^3 or faster), not 1/r.

### 11.2 What the WZW DOES give

The mixed WZW term on Cl+(3,0,1) provides:

- **Rotational coupling**: For a SPINNING Skyrmion, the WZW term couples the degenerate
  sector to the angular velocity, contributing to the moment of inertia and N/P mass splitting.
  This is the standard omega-meson physics.

- **Topological quantization**: The WZW coefficient is topologically quantized (by pi_5(S^3) = Z_2),
  giving a parameter-free coupling. But this coupling is to the baryon CURRENT, not the
  baryon DENSITY in the static sector.

- **Parity violation**: The WZW term is parity-odd, meaning it distinguishes between
  Skyrmions and anti-Skyrmions in time-dependent processes but not in static configurations.

### 11.3 Implications for the gravity program

The WZW approach to gravity via the degenerate sector faces fundamental obstacles:

1. **Static coupling requires a different mechanism**: The WZW term cannot provide the
   static B^0 tau coupling needed for 1/r gravitational potential. A different algebraic
   origin is needed.

2. **The monopole problem is robust**: ANY term in the Lagrangian that is linear in tau
   (and involves only q and its derivatives) will produce an equation of the form
   (differential operator on q) = c * B^0. Unless the differential operator has unit
   monopole (which nabla^2 does not for any component of q), the equation is inconsistent.

3. **Possible escape routes**:
   a. A NON-LINEAR coupling (tau^2 B^0 or similar) could avoid the monopole mismatch
      but would change the asymptotic behavior.
   b. The degenerate sector could couple to B^0 through a CONSTRAINT (not a Lagrangian
      term) — this is Track B.
   c. The gravity mechanism might not involve tau at all but rather the EFFECTIVE METRIC
      modification from the full Cl+(3,0,1) structure.

### 11.4 Technical summary

| Question | Answer |
|----------|--------|
| Does <L^5>_0 see the degenerate sector? | NO (grade obstruction) |
| Does <L^5>_4 see the degenerate sector? | YES (at order e_0) |
| Does <L^5>_4 contain B^0 tau? | YES, but only for B^mu w_mu coupling |
| Is this coupling nontrivial on static hedgehog? | NO (B^i = 0, w_0 = 0) |
| Does the quadratic cross-term give B^0 tau? | NO (monopole mismatch) |
| Is the WZW coefficient quantized? | YES (pi_5(S^3) = Z_2) |
| Does this give 1/r gravity? | NO |

### 11.5 The deeper issue

The fundamental reason the WZW approach fails for static gravity is that the WZW term
is a TOPOLOGICAL term — it depends on the topology of the field configuration but not on
the metric. Gravity, even in this framework, needs to couple to the ENERGY DENSITY
(or at least to a quantity with nonzero monopole moment), not to a topological density
that integrates to zero on non-compact spaces.

The baryon current B^mu is topological (its integral gives the integer baryon number),
and the WZW term couples to B^mu. But the STATIC part of B^mu (i.e., B^0) already
integrates to the correct monopole value (B = 1). The problem is that B^0 tau appears
only in the CURRENT coupling B^mu A_mu, where A_mu is a 1-form, and for a static
configuration the temporal component A_0 necessarily involves TIME derivatives of tau
(which vanish for a static field).

This is a manifestation of the general principle that topological terms contribute to
dynamics (time evolution, scattering, quantization) but not to statics (bound states,
potentials, forces between static objects).

---

## Appendix A: Detailed Product Tables

### A.1 Spatial bivector products

sigma_1 = e_{23}, sigma_2 = e_{31}, sigma_3 = e_{12}

    sigma_a sigma_b = -delta_{ab} + epsilon_{abc} sigma_c

Table:
    sigma_1 sigma_1 = -1
    sigma_1 sigma_2 = sigma_3
    sigma_1 sigma_3 = -sigma_2
    sigma_2 sigma_1 = -sigma_3
    sigma_2 sigma_2 = -1
    sigma_2 sigma_3 = sigma_1
    sigma_3 sigma_1 = sigma_2
    sigma_3 sigma_2 = -sigma_1
    sigma_3 sigma_3 = -1

### A.2 Spatial bivector times mixed bivector

    sigma_a * e_{0b} products (a, b spatial indices):

    e_{23} e_{01} = e_{0123}
    e_{23} e_{02} = -e_{03}
    e_{23} e_{03} = e_{02}
    e_{31} e_{01} = e_{03}
    e_{31} e_{02} = e_{0123}
    e_{31} e_{03} = -e_{01}
    e_{12} e_{01} = -e_{02}
    e_{12} e_{02} = e_{01}
    e_{12} e_{03} = e_{0123}

Summary: sigma_a e_{0b} = delta_{ab} e_{0123} - epsilon_{abc} e_{0c}

(Verify: e_{23} e_{02} = -e_{03}: delta_{12}=0, -epsilon_{123}e_{03} = -e_{03}. Correct.
 e_{23} e_{03} = e_{02}: delta_{13}=0, -epsilon_{132}e_{02} = e_{02}. Correct.
 e_{31} e_{03} = -e_{01}: delta_{23}=0, -epsilon_{231}e_{01} = -e_{01}. Correct.)

### A.3 Spatial bivector times pseudoscalar

    sigma_a * e_{0123} = e_{0123} * sigma_a = -e_{0a}

(Verified in section 1.1: e_{23} e_{0123} = -e_{01}, etc.)

### A.4 Vector times spatial bivector

    e_i * sigma_a = delta_{ia} e_{123} - epsilon_{iac} e_c

(Gives grade 3 for matching index, grade 1 for non-matching.)

### A.5 The pseudoscalar

    e_{0123}^2 = 0    (contains e_0^2 = 0)
    e_{0123} commutes with all even-grade elements of Cl(3,0,1)
    e_{123}^2 = e_1 e_2 e_3 e_1 e_2 e_3 = -1
    e_{123} commutes with all even-grade elements of Cl(3,0)
