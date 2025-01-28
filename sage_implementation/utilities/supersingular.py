"""
Helper functions for the supersingular elliptic curve computations from FESTA
"""

# Sage Imports
from sage.all import ZZ, factor, inverse_mod, CRT

# Local imports
from montgomery_isogenies.kummer_line import KummerLine
from utilities.pairing import weil_pairing_power_two
from utilities.fast_roots import sqrt_Fp2, quadratic_roots
from utilities.pairing import weil_pairing_pari

# =========================================== #
#   Extract coefficent from Montgomery curve  #
# =========================================== #

def montgomery_coefficient(E):
    a_inv = E.a_invariants()
    A = a_inv[1]
    if a_inv != (0, A, 0, 1, 0):
        raise ValueError(
            "Parent function assumes the curve E is in the Montgomery model"
        )
    return A


# =========================================== #
# Compute points of order 2^e and Torsion Bases #
# =========================================== #

def generate_kummer_point(E, x_start=0):
    """
    Generate points on a curve E with x-coordinate
    i + x for x in Fp and i is the generator of Fp^2
    such that i^2 = -1.
    """
    F = E.base_ring()
    one = F.one()

    if x_start:
        x = F.gen() + x_start
    else:
        x = F.gen() + one

    K = KummerLine(E)
    A = montgomery_coefficient(E)

    # Try 1000 times then give up, just protection
    # for infinite loops
    for _ in range(1000):
        y2 = x * (x**2 + A * x + 1)
        if y2.is_square():
            yield K(x)
        x += one # F.random_element() ??

    raise ValueError(
        "Generated 1000 points, something is probably going wrong somewhere."
    )


def has_order_2e(P, e):
    Pem1 = P.double_iter(e-1)
    Pe=Pem1.double()
    return (not Pem1.is_zero()) and Pe.is_zero()

def has_mult_order_2e(c, e):
    cem1 = c**(2**(e-1))
    ce = cem1**2
    return (not cem1 == 1) and (ce == 1)

def generate_point_order_2e(E, e, x_start=0):
    """
    Input:  An elliptic curve E / Fp2
            An integer 2^e dividing (p+1)
    Output: A point P of order 2^e.
    """
    p = E.base().characteristic()
    k = (p + 1) // 2**e

    Ps = generate_kummer_point(E, x_start=x_start)
    for G in Ps:
        P = k*G
        # Case when we randomly picked a point in the n-torsion
        if P.is_zero():
            continue
        # Check that P has order exactly 2**e
        if has_order_2e(P, e):
            yield P.curve_point()

    raise ValueError(f"Never found a point P of order 2^e.")


def compute_point_order_2e(E, e, x_start=0):
    """
    Wrapper function around a generator which returns the first
    point of order 2**e
    """
    return generate_point_order_2e(
        E, e, x_start=x_start,
    ).__next__()


def compute_linearly_independent_point_with_pairing_2e(
    E, P, e, x_start=0
):
    """
    Input:  An elliptic curve E / Fp2
            A point P ∈ E[2^e]
    Output: A point Q such that E[2^e] = <P, Q>
            The Weil pairing e(P,Q)
    """
    Qs = generate_point_order_2e(E, e, x_start=x_start)
    for Q in Qs:
        # Make sure the point is linearly independent
        # ToOptimize: in practice we will always call this with 2^e maximal, so we could use Tate pairings here
        pair = weil_pairing_power_two(P, Q, e)
        if has_mult_order_2e(pair, e):
            Q._order = ZZ(2**e)
            return Q, pair
    raise ValueError("Never found a point Q linearly independent to P")


def torsion_basis_with_pairing_2e(E, e):
    """
    Generate basis of E(Fp^2)[2^e] of supersingular curve

    While computing E[2^e] = <P, Q> we naturally compute the
    Weil pairing e(P,Q), which we also return as in some cases
    the Weil pairing is then used when solving the BiDLP
    """

    P = compute_point_order_2e(E, e)
    Q, ePQ = compute_linearly_independent_point_with_pairing_2e(
        E, P, e, x_start=P[0]
    )

    return P, Q, ePQ


def torsion_basis_2e(E, e, montgomery: False):
    """
    Wrapper function around torsion_basis_with_pairing which only
    returns the torsion basis <P,Q> = E[2^e]
    If montgomery=True, returns the basis such that Q is above the
    Montgomery point.
    """
    P, Q, _ = torsion_basis_with_pairing_2e(E, e)
    if montgomery:
        p = E.base().characteristic()
        return fix_torsion_basis_fast(P, Q, p)
    return P, Q


# =============================================== #
#  Ensure Basis <P,Q> of E[2^k] has (0,0) under Q #
# =============================================== #


def fix_torsion_basis(P, Q, e):
    """
    Set the torsion basis P,Q such that
    2^(k-1)Q = (0,0) to ensure that (0,0)
    is never a kernel of a two isogeny
    """
    cofactor = 2 ** (k - 1)

    R = cofactor * P
    if R[0] == 0:
        return Q, P
    R = cofactor * Q
    if R[0] == 0:
        return P, Q
    return P, P + Q

# use the Tate pairing to quickly fix the basis
# warning: this assumes that P, Q are of maximal 2^e order
def fix_torsion_basis_fast(P, Q, p):
    exp=(p**2-1)//2
    if P[0]**exp == 1:
        return Q, P
    elif Q[0]**exp == 1:
        return P, Q
    else:
        R=P+Q
        assert R[0]**exp == 1
        return P, R


# ============================================= #
# Compute points of order 3^f and Torsion Bases #
# ============================================= #

def has_order_3f(P, f):
    Pem1 = (3**(f-1))*P
    Pe=Pem1.__mul__(3)
    return (not Pem1.is_zero()) and Pe.is_zero()

def has_mult_order_3f(c, f):
    cem1 = c**(3**(f-1))
    ce = cem1**3
    return (not cem1 == 1) and (ce == 1)

def generate_point_order_3f(E, f, x_start=0):
    """
    Input:  An elliptic curve E / Fp2
            An integer 3^f dividing (p+1)
    Output: A point P of order 3^f.
    """
    p = E.base().characteristic()
    assert (p+1) % 3**f == 0, 'Did you set up f correctly?'
    k = (p + 1) // 3**f

    Ps = generate_kummer_point(E, x_start=x_start)
    for G in Ps:
        P = k*G
        # Case when we randomly picked a point in the n-torsion
        if P.is_zero():
            continue
        # Check that P has order exactly 2**e
        if has_order_3f(P, f):
            yield P.curve_point()

    raise ValueError(f"Never found a point P of order 3^f.")


def compute_point_order_3f(E, f, x_start=0):
    """
    Wrapper function around a generator which returns the first
    point of order 2**e
    """
    return generate_point_order_3f(
        E, f, x_start=x_start,
    ).__next__()


def compute_linearly_independent_point_with_pairing_3f(
    E, P, f, x_start=0
):
    """
    Input:  An elliptic curve E / Fp2
            A point P ∈ E[3^f]
    Output: A point Q such that E[3^f] = <P, Q>
            The Weil pairing e(P,Q)
    """
    Qs = generate_point_order_3f(E, f, x_start=x_start)
    for Q in Qs:
        # Make sure the point is linearly independent
        # ToOptimize: in practice we will always call this with 3^f maximal, so we could use Tate pairings here
        pair = weil_pairing_pari(P, Q, 3**f)
        if has_mult_order_3f(pair, f):
            Q._order = ZZ(3**f)
            return Q, pair
    raise ValueError("Never found a point Q linearly independent to P")


def torsion_basis_with_pairing_3f(E, f):
    """
    Generate basis of E(Fp^2)[3^f] of supersingular curve

    While computing E[3^f] = <P, Q> we naturally compute the
    Weil pairing e(P,Q), which we also return as in some cases
    the Weil pairing is then used when solving the BiDLP
    """

    P = compute_point_order_3f(E, f)
    Q, ePQ = compute_linearly_independent_point_with_pairing_3f(
        E, P, f, x_start=P[0]
    )

    return P, Q, ePQ


def torsion_basis_3f(E, f, montgomery: False):
    """
    Wrapper function around torsion_basis_with_pairing which only
    returns the torsion basis <P,Q> = E[3^f]
    If montgomery=True, returns the basis such that Q is above the
    Montgomery point.
    """
    P, Q, _ = torsion_basis_with_pairing_3f(E, f)
    if montgomery:
        p = E.base().characteristic()
        return fix_torsion_basis_fast(P, Q, p)
    return P, Q
