"""
Helper functions for various computations associated to
the quaternion algebra, and ideals and orders of the
quaternion algebra.

Some of these functions could be ported to SageMath. It's a TODO
for when SQISign is being less actively worked on.
"""

# Sage Imports
from sage.all import (
    randint, ZZ, ceil, log, gcd,
    Matrix, vector, prod,
    denominator,
)
from utilities.logging_helper import logging
import time
logger = logging.getLogger(__name__)

# Local imports

# ================================================ #
#  Helpers for elements of the quaternion algebra  #
# ================================================ #


def quadratic_norm(x, y):
    """
    Given two integers x,y, which correspond
    to the element x + ωy ∈ R[ω], returns
    Nrd(x + ωy)

    WARNING: Assumes that ω^2 = -1

    Note: This function implements the norm
    function f(x,y) in the SQISign papers.

    For SQISign, we have ω = i and i^2 = -1
    so f(x,y) = x**2 + y**2 which is more
    efficient
    """
    return ZZ(x) ** 2 + ZZ(y) ** 2


def quaternion_change_basis(γ, O):
    """
    Computes the coefficients of a quaternion γ
    in the basis of a given order O
    """
    O_matrix = Matrix([b.coefficient_tuple() for b in O.basis()])
    γ_vector = vector(γ.coefficient_tuple())
    γO_coeffs = γ_vector * O_matrix.inverse()

    assert γ == sum([a * b for a, b in zip(γO_coeffs, O.basis())])
    return γO_coeffs


def quaternion_basis_gcd(γ, O):
    """
    Computes the gcd of the coefficients of a
    quaternion γ in the basis of a given order O
    """
    γO_coeffs = quaternion_change_basis(γ, O)
    return gcd(γO_coeffs)


# ============================================== #
#  Helpers for ideals of the quaternion algebra  #
# ============================================== #


def multiply_ideals(I, J, beta=None, check=True):
    """
    Computes I*J when O_R(I) ≃ O_L(J)

    If these orders do not match, we must provide
    an isomorphism which takes O_L(J) to O_R(I)
    """

    if beta:
        J = beta ** (-1) * J * beta
    if check:
        assert (
            I.right_order() == J.left_order()
        ), "Orders do not match"
    return I * J


def is_integral(I):
    """
    Checks whether the input ideal is integral.
    """
    return all([b in I.left_order() for b in I.basis()])


def ideal_basis_gcd(I):
    """
    Computes the gcd of the coefficients of
    the ideal written as a linear combination
    of the basis of its left order.
    """
    I_basis = I.basis_matrix()
    O_basis = I.left_order().unit_ideal().basis_matrix()

    # Write I in the basis of its left order
    M = I_basis * O_basis.inverse()
    return gcd((gcd(M_row) for M_row in M))


def is_cyclic(I):
    """
    Computes whether the input ideal is cyclic,
    all the work is done by the helper function
    `ideal_basis_gcd()`.
    """
    return ideal_basis_gcd(I) == 1


def make_cyclic(I):
    """
    Given an ideal I, returns a cyclic ideal by dividing
    out the scalar factor g = ideal_basis_gcd(I)
    """
    g = ideal_basis_gcd(I)
    # Ideal was already cyclic
    if g == 1:
        return I, g

    print(f"DEBUG [make_cyclic]: Ideal is not cyclic, removing scalar factor: {g = }")
    J = I.scale(1 / g)

    return J, g


def reduced_basis(I, check=False):
    """
    Computes the Minkowski reduced basis of the
    input ideal. Note: this produces the same
    result for all ideals in the equivalence class
    so corresponds to the reduced basis of the
    smallest equivalent ideal to I

    Input: an ideal
    Output: A Minkowski-reduced basis

    Optional: when check is True, the basis is
              checked against the Minkowski bounds
    """

    def _matrix_to_gens(M, B):
        """
        Converts from a matrix to generators in the quat.
        algebra
        """
        return [sum(c * g for c, g in zip(row, B)) for row in M]

    B = I.basis()
    G = I.gram_matrix()
    U = G.LLL_gram().transpose()

    reduced_basis_elements = _matrix_to_gens(U, B)

    if check:
        norm_product = 16 * prod([x.reduced_norm() for x in reduced_basis_elements])
        tmp = p**2 * I.norm() ** 4
        assert norm_product <= 4 * tmp, "Minkowski reduced basis is too large"
        assert norm_product >= tmp, "Minkowski reduced basis is too small"

    return reduced_basis_elements


def small_equivalent_ideal(I, reduced_basis_elements=None):
    """
    Computes the Minkowski reduced basis of the
    ideal and returns the smallest equivalent
    ideal J = Iα ~ I.
    """
    nI = I.norm()

    if not reduced_basis_elements:
        reduced_basis_elements = reduced_basis(I)

    b0 = reduced_basis_elements[0]
    if b0.reduced_norm() == nI**2:
        return I
    return I * I.right_order().left_ideal([b0.conjugate() / nI])


def equivalent_left_ideals(I, J):
    """
    SageMath has this impl. for right ideals
    only. To work around this, we can first
    take the conjugate of I and J.

    TODO: write a function which does this without
          conjugates?
    """
    return I.conjugate().is_equivalent(J.conjugate())


def equivalent_right_ideals(I, J):
    """
    Sage can do this for us already only making a
    wrapped so it matches the above.
    """
    return I.is_equivalent(J)


def invert_ideal(I):
    """
    Computes the inverse of the ideal which is
    the conjugate of I divided by its norm
    """
    return I.conjugate().scale(1 / I.norm())


def left_isomorphism(I, J):
    """
    Given two isomorphic left ideals I, J computes
    α such that J = I*α
    """
    B = I.quaternion_algebra()

    if B != J.quaternion_algebra():
        raise ValueError("Arguments must be ideals in the same algebra.")

    if I.left_order() != J.left_order():
        raise ValueError("Arguments must have the same left order.")

    IJ = I.conjugate() * J
    L = reduced_basis(IJ)
    for t in L:
        α = t / I.norm()
        if J == I * α:
            return α

    raise ValueError("Could not find a left isomorphism...")

def right_isomorphism(I, J):
    """
    Given two isomorphic right ideals I, J computes
    α such that J = α*I
    """
    αbar = left_isomorphism(I.conjugate(), J.conjugate())
    return αbar.conjugate()



def chi(a, I, adjust_right_order=False):
    r"""
    From section 3.2 in Antonin's thesis.
    Calculates the equivalent ideal of I, of norm(a)
    Based on the surjection from I \ {0} to the set of equivalent ideals of I
    Obtained by a → I * (a_bar / n(I))
    """
    # by default, J will have same left order, and isomorphic right order.
    # With adjust_right_order = True it will have isomorphic left order and same right order
    beta=(a.conjugate() / I.norm())
    J= I * beta
    if adjust_right_order: #be sure the right order of J match the one on I
        J=beta*J*beta**(-1)
    return J


def chi_inverse(I, J):
    """
    Computes the element α such that

    J = Chi(α, I)
    """
    # Compute the element α
    a = left_isomorphism(I, J)
    assert J == I * a, "Left isomorphism to find the element 'a' failed... why?"
    α = (a * I.norm()).conjugate()
    assert J == chi(α, I), "Something about α is wrong!"
    return α


def scaled_norm(a, I):
    """
    Returns Nrd(a) / n(I), the norm of chi(a, I)
    """
    N = I.norm()
    return a.reduced_norm() / N


def ideal_generator(I, coprime_factor=1):
    """
    Given an ideal I of norm D, finds a generator
    α such that I = O(α,D) = Oα + OD

    Optional: Enure the norm of the generator is coprime
    to the integer coprime_factor
    """
    # This is a stupid hack to get the prime from the ideal.
    p = -I.quaternion_algebra().invariants()[1]
    OI = I.left_order()
    D = ZZ(I.norm())
    bound = ceil(4 * log(p))

    gcd_norm = coprime_factor * D**2

    # Stop infinite loops.
    for _ in range(1000):
        α = sum([b * randint(-bound, bound) for b in I.basis()])
        if gcd(ZZ(α.reduced_norm()), gcd_norm) == D:
            assert I == OI * α + OI * D
            return α
    raise ValueError(f"Cannot find a good α for D = {D}, I = {I}, n(I) = {D}")


def eichler_order_from_ideal(I):
    """
    The Eichler order is the intersection
    of two orders.

    Given an ideal I, we compute the Eichler
    order ℤ + I from the intersection of the
    left and right orders of I

    Proposition 1 (SQISign paper):

    EichlerOrder = O0 ∩ O = OL(I) ∩ OR(I) = ℤ + I.
    """
    return I.left_order().intersection(I.right_order())


# ========================================= #
#  Pushforward and Pullback of ideals to O0 #
# ========================================= #


def pullback_ideal(I, Iτ, O0=None, O1=None, check=True):
    """
    Input: Ideal I with left order O1
           Connecting ideal Iτ with left order O0
           and right order O1
    Output The ideal given by the pullback [Iτ]^* I
    """
    if O0 is None:
        O0=Iτ.left_order()
    if check:
        if O1 is None:
            O1=I.left_order()
        assert I.left_order() == O1
        assert Iτ.left_order() == O0
        assert Iτ.right_order() == O1

    N = ZZ(I.norm())
    Nτ = ZZ(Iτ.norm())

    α = ideal_generator(I)
    return O0 * N + O0 * α * Nτ


def pushforward_ideal(I, Iτ, O0=None, O1=None, check=True):
    """
    Input: Ideal I left order O0
           Connecting ideal Iτ with left order O0
           and right order O1
    Output The ideal given by the pushforward [Iτ]_* I
    """
    if O1 is None:
        O1=Iτ.right_order()
    if check:
        if O0 is None:
            O0=I.left_order()
        assert I.left_order() == O0
        assert Iτ.left_order() == O0
        assert Iτ.right_order() == O1

    N = ZZ(I.norm())
    Nτ = ZZ(Iτ.norm())

    K = I.intersection(O1 * Nτ)
    α = ideal_generator(K)
    return O1 * N + O1 * (α / Nτ)

def any_pushforward_ideal(I1, I2):
    """
    Compute the pushforward ideal of I1 by I2.
    If nrd(I2) is not coprime to nrd(I1), first compute an equivalent ideal
    J2 to it. The result depends on the choice of J2, but we are sure it will be of the same norm as I1.
    """
    J2 = I2
    nI1 = ZZ(I1.norm())
    nI2 = ZZ(I2.norm())
    if gcd(nI2, nI1) != 1:
        logger.info("Pushforward ideals: Computing an equivalent ideal")
        def condition(α):
            n = ZZ(α.reduced_norm())
            return gcd(n//nI2, nI1)==1
        J2, _αcomm, _c=random_equivalent_ideal(I2, condition=condition, adjust_right_order = True)
        assert J2.right_order() == I2.right_order()
    # in multiply_ideals, we use check=False because we already
    # constructed our ideals to have matching orders, so no need to
    # waste a check.
    # Here the check=False is mandatory because (if
    # random_ideal_has_prime_norm is False) we adjusted J2 to
    # have the same right order as I2, so the left order is only
    # isomorphic but not equal to O0, so the check would fail.
    I = pushforward_ideal(I1, J2, check=False)
    assert ZZ(I.norm()) == ZZ(I1.norm())
    return I

def connecting_ideal(O0, O1):
    I = O0 * O1
    I = I * ZZ(denominator(I.norm()))
    return I

# ############################
# Equivalent ideals
# ############################

# enumerate a basis of a quaternion ideal
def enumerate_basis(bound=10):
    import itertools
    for i1 in range(0,bound+1):
        for i2 in range(0, i1+1):
            for i3 in range(0, i2+1):
                for i4 in range(0, i3+1):
                    if i1==0 and i2==0 and i3==0 and i4==0:
                        continue
                    els=set(); base=(i1, i2, i3, i4)
                    all_signs = itertools.product([-1, 1], repeat=4)
                    for signs in all_signs:
                        els=els.union(itertools.permutations((element * sign for element, sign in zip(base, signs))))
                    for el in els:
                        yield el

def find_equivalent_ideal(J, basis=None, generator=None, bound=100, condition=None, max_attempts=10**6, adjust_right_order=False):
    if basis is None:
        basis = reduced_basis(J)
    γ1, γ2, γ3, γ4 = basis
    counter = 0
    nJ = ZZ(J.norm())

    def sample_endo():
        nonlocal counter
        init_time=time.time()
        for i1, i2, i3, i4 in generator():
            counter += 1
            α = i1 * γ1 + i2 * γ2 + i3 * γ3  + i4 * γ4
            if condition:
                if condition(α, ZZ(α.reduced_norm())//nJ):
                    logger.debug(f"Found an equivalent ideal of {ZZ(α.reduced_norm()//J.norm()).nbits()}b in {counter} attempts in {(time.time() - init_time):.3f}s [bound={bound}]")
                    return α
            else:
                return α

    α=sample_endo()
    if α is None:
        I=None
    else:
        I=chi(α, J, adjust_right_order=adjust_right_order)

    return I, α, counter

def smallish_equivalent_ideal(J, bound=100, **kws):
    def generator():
        yield from enumerate_basis(bound=bound)
    return find_equivalent_ideal(J, generator=generator, **kws)

def random_equivalent_ideal(J, bound=100, max_attempts=10**6, **kws):
    def generator():
        for counter in range(0, max_attempts):
            yield randint(-bound, bound), randint(-bound, bound), randint(-bound, bound), randint(-bound, bound)
    return find_equivalent_ideal(J, generator=generator, **kws)
