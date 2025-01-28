"""
Underneath everything, our isogenies are on the Kummer line

L : x^3 + Ax^2 + x

And we perform our x-only isogenies by working with the x-coordinates
represented projectively as x(P) = (X : Z) which we call KummerPoints

However, for FESTA(+) we always need the full point eventually for either
additions or the (2,2)-isogeny, so we need a way to recover the full point.

The trick we use is that we always evaluate our isogenies on torsion bases,
so we can use the Weil pairing to recover phi(P) up to an overall sign.

This file takes elliptic curves and points on these curves, maps them to the
Kummer line, performs fast x-only isogeny computations and then lifts the
result back to full points on the codomain curves.
"""

# Sage imports
from sage.structure.element import RingElement

# Local Imports
from montgomery_isogenies.kummer_line import KummerLine
from montgomery_isogenies.kummer_isogeny import KummerLineIsogeny
# from utilities.supersingular import torsion_basis
# from utilities.pairing import weil_pairing_pari

# =========================================================== #
#    Compute an isogeny and codomain using x-only algorithms  #
# =========================================================== #


def isogeny_from_scalar_x_only(E, D, m, basis=None):
    """
    Computes a D-degree isogeny from E using
    x-only arithmetic and returns the KummerIsogeny
    together with the codomain curve.

    The isogeny has a kernel K which is computed from
    the canonical basis E[D] = <P,Q> and given scalar(s)
    of the form:
        K = P + [m]Q     or     K = [a]P + [b]Q
    depending on whether m is a scalar, or a length two
    tuple of scalars
    """
    # Allow a precomputed basis
    if not basis:
        raise NotImplementedError('we miss torsion_basis generator')
        P, Q = torsion_basis(E, D)
    else:
        P, Q = basis

    # Allow either an integer or tuple of integers
    if isinstance(m, RingElement) or isinstance(m, int):
        K = P + m * Q
    else:
        assert len(m) == 2
        K = m[0] * P + m[1] * Q

    # Map curve and kernel to Kummer Line
    L = KummerLine(E)
    xK = L(K)

    # Use x-only arithmetic to compute an isogeny
    # and codomain
    phi = KummerLineIsogeny(L, xK, D)

    # Compute the curve from the Kummer Line
    codomain = phi.codomain().curve()

    # Speed up SageMath by setting the order of the curve
    p = E.base_ring().characteristic()
    codomain.set_order((p + 1) ** 2, num_checks=0)

    return phi, codomain


def isogeny_from_kernel_x_only(E, D, kernel):
    """
    Computes a D-degree isogeny from E using
    x-only arithmetic and returns the KummerIsogeny
    together with the codomain curve.
    """
    # Map curve and kernel to Kummer Line
    L = KummerLine(E)
    xK = L(kernel)

    # Use x-only arithmetic to compute an isogeny
    # and codomain
    phi = KummerLineIsogeny(L, xK, D)

    # Compute the curve from the Kummer Line
    codomain = phi.codomain().curve()

    # Speed up SageMath by setting the order of the curve
    p = E.base_ring().characteristic()
    codomain.set_order((p + 1) ** 2, num_checks=0)

    return phi, codomain


# ================================================= #
#    Evaluate an x-only isogeny on a torsion basis  #
# ================================================= #


def lift_image_to_curve(P, Q, PQ):
    """
    Given the image of xP, xQ and x(P + Q) we can
    recover the affine points P, Q up to an overall
    sign.
    """
    # Recover the montgomery coefficent
    A = P.parent().a()

    # Lift the points to the curve
    imP = P.curve_point()
    imQ = Q.curve_point()

    # Get the affine x-coordinates
    xP, xQ, xPQ = P.x(), Q.x(), PQ.x()

    # Compute two pairings
    lhs = A + xP + xQ + xPQ
    rhs = ((imQ[1] - imP[1]) / (xQ - xP)) ** 2

    # Correct the sign
    if lhs == rhs:
        imQ = -imQ

    return imP, imQ


def evaluate_isogeny_x_only(phi, P, Q):
    """
    Given an x-only isogeny phi degree d, and the torsion basis
    <P,Q> = E[n], compute the image of the torsion basis up to
    and overall sign: ±phi(P), ±phi(Q)

    Does this by evaluating KummerPoints with a KummerIsogeny
    and lifts them back to the curve using the Weil pairing
    trick in `lift_image_to_curve`
    """
    # Domain of isogeny
    L0 = phi.domain()

    # Extract x-coordinates from points and convert to KummerPoints
    xP, xQ, xPQ = L0(P[0]), L0(Q[0]), L0((P - Q)[0])

    # Evaluate the isogeny
    ximP, ximQ, ximPQ = phi(xP), phi(xQ), phi(xPQ)

    # Use Weil pairing trick to get y-coordinate back
    imP, imQ = lift_image_to_curve(ximP, ximQ, ximPQ)

    return imP, imQ

def evaluate_isogeny_x_only_one_point(phi, P):

    L0 = phi.domain()

    # Extract x-coordinates from points and convert to KummerPoints
    xP = L0(P[0])

    # Evaluate the isogeny
    ximP = phi(xP)

    imP = ximP.curve_point()

    return imP
