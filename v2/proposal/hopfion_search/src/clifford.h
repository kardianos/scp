/*
 * clifford.h — Even subalgebra of Cl(3,0,1) (Projective Geometric Algebra)
 *
 * The 8 basis elements of Cl+(3,0,1):
 *   Grade 0: 1         (scalar)
 *   Grade 2: e23, e31, e12  (spatial bivectors)
 *            e01, e02, e03  (degenerate bivectors)
 *   Grade 4: e0123     (pseudoscalar)
 *
 * Field representation: Psi = s + F1*e23 + F2*e31 + F3*e12
 *                            + J1*e01 + J2*e02 + J3*e03 + P*e0123
 *
 * Dual quaternion decomposition: Psi = q + e0*p
 *   q = s + F1*e23 + F2*e31 + F3*e12   (bulk, in Cl+(3,0))
 *   p = J1*e1 + J2*e2 + J3*e3 + P*e123 (weight, in Cl-(3,0))
 *
 * Key property: e0^2 = 0, so <Psi Psi~>_0 = s^2 + F1^2 + F2^2 + F3^2
 *              (blind to degenerate sector)
 */

#ifndef CLIFFORD_H
#define CLIFFORD_H

/* Even multivector in Cl+(3,0,1) */
typedef struct {
    double s;           /* scalar (1) */
    double f1, f2, f3;  /* spatial bivector (e23, e31, e12) */
    double j1, j2, j3;  /* degenerate bivector (e01, e02, e03) */
    double p;           /* pseudoscalar (e0123) */
} Multivector;

/* Geometric product of two even multivectors in Cl+(3,0,1).
 *
 * Uses the dual quaternion product: (q1 + e0*p1)(q2 + e0*p2)
 *   = q1*q2 + e0*(q1*p2 + p1*q2)
 * where q*q is quaternion product (in Cl+(3,0))
 * and q*p, p*q are the mixed products (Cl+(3,0) x Cl-(3,0) -> Cl-(3,0)).
 */
static inline Multivector mv_mul(Multivector a, Multivector b)
{
    Multivector r;

    /* Quaternion product: q1*q2 in Cl+(3,0)
     * Basis: {1, e23, e31, e12} with e23^2 = e31^2 = e12^2 = -1
     * e23*e31 = -e12,  e31*e23 = +e12
     * e23*e12 = +e31,  e12*e23 = -e31
     * e31*e12 = -e23,  e12*e31 = +e23
     */
    r.s  = a.s*b.s  - a.f1*b.f1 - a.f2*b.f2 - a.f3*b.f3;
    r.f1 = a.s*b.f1 + a.f1*b.s  - a.f2*b.f3 + a.f3*b.f2;
    r.f2 = a.s*b.f2 + a.f1*b.f3 + a.f2*b.s  - a.f3*b.f1;
    r.f3 = a.s*b.f3 - a.f1*b.f2 + a.f2*b.f1 + a.f3*b.s;

    /* Dual part: e0*(q1*p2 + p1*q2)
     *
     * q*p product (Cl+(3,0) x Cl-(3,0) -> Cl-(3,0)):
     *   Basis: q={1,e23,e31,e12}, p={e1,e2,e3,e123}
     *   Result mapped to (J1,J2,J3,P) via e0*p -> (J=p_vec, P=p_pseudo)
     *
     *   (qp).e1   = a*x - b*w - c*z + d*y
     *   (qp).e2   = a*y + b*z - c*w - d*x
     *   (qp).e3   = a*z - b*y + c*x - d*w
     *   (qp).e123 = a*w + b*x + c*y + d*z
     *
     * where q=(a,b,c,d)=(s,f1,f2,f3), p=(x,y,z,w)=(j1,j2,j3,P)
     *
     * p*q product (Cl-(3,0) x Cl+(3,0) -> Cl-(3,0)):
     *   (pq).e1   = x*a - w*b + z*c - y*d
     *   (pq).e2   = x*d + y*a - z*b - w*c
     *   (pq).e3   = -x*c + y*b + z*a - w*d
     *   (pq).e123 = x*b + y*c + z*d + w*a
     */

    /* q_a * p_b  (a's quaternion part times b's weight part) */
    double qp_j1 = a.s*b.j1 - a.f1*b.p  - a.f2*b.j3 + a.f3*b.j2;
    double qp_j2 = a.s*b.j2 + a.f1*b.j3 - a.f2*b.p  - a.f3*b.j1;
    double qp_j3 = a.s*b.j3 - a.f1*b.j2 + a.f2*b.j1 - a.f3*b.p;
    double qp_p  = a.s*b.p  + a.f1*b.j1 + a.f2*b.j2 + a.f3*b.j3;

    /* p_a * q_b  (a's weight part times b's quaternion part) */
    double pq_j1 = a.j1*b.s  - a.p*b.f1  + a.j3*b.f2 - a.j2*b.f3;
    double pq_j2 = a.j1*b.f3 + a.j2*b.s  - a.j3*b.f1 - a.p*b.f2;
    double pq_j3 = -a.j1*b.f2 + a.j2*b.f1 + a.j3*b.s  - a.p*b.f3;
    double pq_p  = a.j1*b.f1 + a.j2*b.f2 + a.j3*b.f3 + a.p*b.s;

    /* Dual part = q_a*p_b + p_a*q_b */
    r.j1 = qp_j1 + pq_j1;
    r.j2 = qp_j2 + pq_j2;
    r.j3 = qp_j3 + pq_j3;
    r.p  = qp_p  + pq_p;

    return r;
}

/* Reversion: reverse order of basis vectors
 * Grade 0: +1,  Grade 2: -1,  Grade 4: +1
 * So: s -> s, F -> -F, J -> -J, P -> +P
 */
static inline Multivector mv_rev(Multivector a)
{
    Multivector r;
    r.s  =  a.s;
    r.f1 = -a.f1;  r.f2 = -a.f2;  r.f3 = -a.f3;
    r.j1 = -a.j1;  r.j2 = -a.j2;  r.j3 = -a.j3;
    r.p  =  a.p;
    return r;
}

/* Scalar part extraction: <M>_0 */
static inline double mv_scalar(Multivector a)
{
    return a.s;
}

/* Bulk norm squared: <Psi Psi~>_0 = s^2 + |F|^2
 * This is the norm that appears in the potential V.
 * Note: blind to degenerate sector (J, P) because e0^2 = 0.
 */
static inline double mv_bulk_norm2(Multivector a)
{
    return a.s*a.s + a.f1*a.f1 + a.f2*a.f2 + a.f3*a.f3;
}

/* Weight norm squared: |p|^2 = |J|^2 + tau^2
 * This is the norm that appears in V_D.
 */
static inline double mv_weight_norm2(Multivector a)
{
    return a.j1*a.j1 + a.j2*a.j2 + a.j3*a.j3 + a.p*a.p;
}

/* Commutator: [A, B] = AB - BA */
static inline Multivector mv_commutator(Multivector a, Multivector b)
{
    Multivector ab = mv_mul(a, b);
    Multivector ba = mv_mul(b, a);
    Multivector r;
    r.s  = ab.s  - ba.s;
    r.f1 = ab.f1 - ba.f1;
    r.f2 = ab.f2 - ba.f2;
    r.f3 = ab.f3 - ba.f3;
    r.j1 = ab.j1 - ba.j1;
    r.j2 = ab.j2 - ba.j2;
    r.j3 = ab.j3 - ba.j3;
    r.p  = ab.p  - ba.p;
    return r;
}

/* Linear operations */
static inline Multivector mv_add(Multivector a, Multivector b)
{
    Multivector r;
    r.s  = a.s  + b.s;
    r.f1 = a.f1 + b.f1;  r.f2 = a.f2 + b.f2;  r.f3 = a.f3 + b.f3;
    r.j1 = a.j1 + b.j1;  r.j2 = a.j2 + b.j2;  r.j3 = a.j3 + b.j3;
    r.p  = a.p  + b.p;
    return r;
}

static inline Multivector mv_sub(Multivector a, Multivector b)
{
    Multivector r;
    r.s  = a.s  - b.s;
    r.f1 = a.f1 - b.f1;  r.f2 = a.f2 - b.f2;  r.f3 = a.f3 - b.f3;
    r.j1 = a.j1 - b.j1;  r.j2 = a.j2 - b.j2;  r.j3 = a.j3 - b.j3;
    r.p  = a.p  - b.p;
    return r;
}

static inline Multivector mv_scale(double c, Multivector a)
{
    Multivector r;
    r.s  = c*a.s;
    r.f1 = c*a.f1;  r.f2 = c*a.f2;  r.f3 = c*a.f3;
    r.j1 = c*a.j1;  r.j2 = c*a.j2;  r.j3 = c*a.j3;
    r.p  = c*a.p;
    return r;
}

static inline Multivector mv_zero(void)
{
    Multivector r = {0, 0, 0, 0, 0, 0, 0, 0};
    return r;
}

/* Dot product of two multivectors as 8-vectors (for gradient flow) */
static inline double mv_dot(Multivector a, Multivector b)
{
    return a.s*b.s + a.f1*b.f1 + a.f2*b.f2 + a.f3*b.f3
         + a.j1*b.j1 + a.j2*b.j2 + a.j3*b.j3 + a.p*b.p;
}

#endif /* CLIFFORD_H */
