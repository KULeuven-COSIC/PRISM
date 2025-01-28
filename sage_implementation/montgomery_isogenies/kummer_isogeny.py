"""
Implementation of x-only isogenies between the Kummer Lines of 
Montgomery Curves

Taken from:
https://github.com/jack4818/KummerIsogeny

===========================================================================

USAGE:

phi = KummerLineIsogeny(domain, kernel, degree)

The codomain is accessed using `phi.codomain()` and the  elliptic curve 
can be lifted from the KummerLine `phi.dodomain().curve()`.

Evaluation of the isogeny is done via `phi(xQ)` for some KummerPoint `xQ`.

NOTE:

Where the degree can be composite, but for efficiency needs to be smooth.
For a 2-isogeny, the point P = (0,0) cannot be used as a kernel 

========================================================================

INFO:

Heavily inspired by the SageMath isogeny classes, this file implements
x-only Montgomery isogenies using the KummerLine and KummerPoint classes
from `kummer_line.py` for the (co)domains and kernel points.

The algorithms in this file all come from the following literature:

Vélu-like algorithms:

    Even torsion algorithms from: https://ia.cr/2017/1198
    Computing Isogenies between Montgomery Curves Using the Action of (0, 0)
    Joost Renes

    Odd torsion algorithms: https://ia.cr/2017/504.pdf
    A simple and compact algorithm for SIDH with arbitrary degree isogenies
    Craig Costello and Huseyin Hisil

    Codomain computation for velu formula from: https://ia.cr/2018/782
    A faster way to the CSIDH
    Michael Meyer and Steffen Reith

VéluSqrt for large ell isogenies 

    VéluSqrt: https://velusqrt.isogeny.org/
    Faster computation of isogenies of large prime degree
    Daniel J. Bernstein, Luca De Feo, Antonin Leroux, Benjamin Smith

Future Work: 

- Optimise VéluSqrt, it seems to be underperforming with a threshold of about 1000
  rather than 100 
- Include isomorphisms of Kummer Lines
- allow composition by defining __mul__ on isogenies to create a composite isogeny
"""


# Sage imports
from sage.all import prod, ZZ, PolynomialRing
from sage.rings.generic import ProductTree

# Local imports
from .kummer_line import KummerLine, KummerPoint, pari

# =================================================== #
# Generic class for creating an isogeny between       #
# KummerLines of Montgomery model curves using x-only #
# arithmetic                                          #
# =================================================== #


class KummerLineIsogeny_Generic:
    """
    Generic class for Kummer Line isogenies which we build on top of for
    the Vélu, VéluSqrt and Composite isogeny classes
    """

    def __init__(self):
        self._degree = None
        self._domain = None
        self._codomain = None
        pass

    def __repr__(self):
        return f"Isogeny of degree {(self._degree).factor()} from {self._domain} to {self._codomain}"

    @staticmethod
    def validate_input(domain, kernel, degree, check=True):
        """
        Helper function to check the input to the isogeny class is well-formed
        """
        if not isinstance(domain, KummerLine):
            raise ValueError(f"not a kummer line: {domain}")

        if not isinstance(kernel, KummerPoint):
            raise ValueError(f"not a kummer point: {kernel}")

        if kernel.parent() != domain:
            raise ValueError(f"Kernel {kernel} is not a point on {domain}")

        if check:
            # TODO actually check order with has_order_D function
            assert (
                degree * kernel
            ).is_zero(), "Input point does not have correct order"

    def domain(self):
        """
        Return the domain of the isogeny
        """
        return self._domain

    def codomain(self):
        """
        Return the codomain of the isogeny
        """
        return self._codomain

    def degree(self):
        """
        Return the degree of the isogeny
        """
        return self._degree


# =================================================== #
# Computation of isogenies between Kummer lines using #
# x-only formula by Costello-Hisil-Renes              #
# =================================================== #


class KummerLineIsogeny_Velu(KummerLineIsogeny_Generic):
    """
    Computes prime degree isogenies with Vélu-like formula.

    - When ell is odd, we use Costello-Hisil (https://ia.cr/2017/504)
    - When ell is even, we can use Renes (https://ia.cr/2017/1198) providing
    that the kernel is not (0,0)

    TODO: use isomorphisms to change the model of the curve if (0,0) is a
    kernel point
    """

    def __init__(self, domain, kernel, degree, check=True):
        # Check the input to the isogeny is well-formed
        self.validate_input(domain, kernel, degree, check=check)

        # Set kernel and degree and domain
        self._degree = degree
        self._kernel = kernel
        self._domain = domain

        # Compute the codomain, we need different formula for even and
        # odd degree
        if self._degree == 2:
            # We cannot use the point (0 : 0 : 1) as the kernel for
            # these formula
            assert self._kernel.XZ()[0]

        # Compute the codomain
        self._codomain = self._compute_codomain()

    def __call__(self, P):
        """
        phi(xP) evaluates the Kummer point xP
        """
        if not isinstance(P, KummerPoint):
            raise ValueError
        if self._degree == 2:
            return self._evaluate_isogeny_even(P)
        return self._evaluate_isogeny(P)

    def _precompute_edwards_multiples(self, d):
        """
        These multiples are used in both codomain
        computation and isogeny evaluation. We precompute
        them once during initialisation and we can then
        reuse them for every evaluation
        """
        # Compute the [i]K for i in [1...d]
        K_muls = self._kernel.multiples()
        E_muls = []
        for _ in range(d):
            Ki = next(K_muls)
            KX, KZ = Ki.XZ()
            YE = KX - KZ
            ZE = KX + KZ
            E_muls.append((YE, ZE))
        return E_muls

    def _compute_codomain_constants(self):
        """
        When ell is odd, we compute the codomain using the Meyer and Reith
        Twised Edwards trick (https://ia.cr/2018/782)
        """
        # Extract Montgomery constants
        A, C = self._domain.extract_constants()

        # Compute and store pairs of points for later evaluation
        d = (self._degree - 1) // 2
        self._edwards_multiples = self._precompute_edwards_multiples(d)

        # Convert to twisted Edwards curve parameters (Aed,Ded)
        Ded = C + C
        Aed = A + Ded
        Ded = A - Ded

        # Compute product of Edwards multiples
        prod_Y = 1
        prod_Z = 1
        for EY, EZ in self._edwards_multiples:
            prod_Y *= EY
            prod_Z *= EZ

        # compute prod_Y^8 and prod_Z^8
        prod_Y, prod_Z = prod_Y**2, prod_Z**2
        prod_Y, prod_Z = prod_Y**2, prod_Z**2
        prod_Y, prod_Z = prod_Y**2, prod_Z**2

        # A_new = A_old^ell * prod_Z^8
        # D_new = D_old^ell * prod_Y^8
        Aed = Aed**self._degree * prod_Z
        Ded = Ded**self._degree * prod_Y

        # Change back to Montgomery-parameters
        A = Aed + Ded
        C = Aed - Ded
        A = A + A

        return A, C

    def _compute_codomain_constants_even(self):
        """
        When ell is even, we compute the codomain constants
        using Renes formula
        """
        # Extract kernel point
        XK, ZK = self._kernel.XZ()
        assert XK, "XK Cannot be zero"

        # C = ZK^2
        C = ZK * ZK

        # A = 2*(ZK^2 - 2*XK^2)
        A = XK * XK  # A = XK^2
        A = A + A  # A = 2*XK^2
        A = C - A  # A = ZK^2 - 2XK^2
        A = A + A  # A = 2*(ZK^2 - 2*XK^2)
        return A, C

    def _compute_codomain(self):
        """
        Wrapper function to compute the codomain L = x^3 + x^2A' + x in
        projective coordinates: A' = (A' : C') We use different formula
        depending on whether the isogeny degree ell is even or odd
        """
        # Compute the codomain constants, need different formula for
        # odd and even ell
        if self._degree == 2:
            A_codomain, C_codomain = self._compute_codomain_constants_even()
        else:
            A_codomain, C_codomain = self._compute_codomain_constants()

        # Constuct a new KummerLine
        F = self._domain.base_ring()
        return KummerLine(F, [A_codomain, C_codomain])

    def _evaluate_isogeny(self, P):
        """
        Costello-Hisil (https://ia.cr/2017/504) formula for
        evaluating an odd degree isogeny on the point P
        """
        XP, ZP = P.XZ()
        Psum = XP + ZP
        Pdiff = XP - ZP

        # Loop through the d-multiples, these are
        # precomputed from the codomain computation
        X_new, Z_new = 1, 1
        for EY, EZ in self._edwards_multiples:
            diff_EZ = Pdiff * EZ
            sum_EY = EY * Psum
            X_new *= diff_EZ + sum_EY
            Z_new *= diff_EZ - sum_EY

        # Square and multiple with original
        X_new = X_new**2 * XP
        Z_new = Z_new**2 * ZP

        return self._codomain((X_new, Z_new))

    def _evaluate_isogeny_even(self, P):
        """
        Renes (https://ia.cr/2017/1198) formula for
        evaluating an even degree isogeny on the point P
        """
        XK, ZK = self._kernel.XZ()
        assert XK, "XK cannot be zero"

        XP, ZP = P.XZ()

        T0 = XK + ZK
        T1 = XK - ZK
        T2 = XP + ZP
        T3 = ZP - XP  # Typo in formula: paper says XP - ZP
        T4 = T3 * T0  # (ZP - XP)(XK + ZK)
        T5 = T2 * T1  # (XP + ZP)(XK - ZK)
        T6 = T4 - T5  # (ZP - XP)(XK + ZK) - (XP + ZP)(XK - ZK)
        T7 = T4 + T5  # (ZP - XP)(XK + ZK) + (XP + ZP)(XK - ZK)
        T8 = XP * T6  # XP * ((ZP - XP)(XK + ZK) - (XP + ZP)(XK - ZK))
        T9 = ZP * T7  # ZP * ((ZP - XP)(XK + ZK) + (XP + ZP)(XK - ZK))

        return self._codomain((T8, T9))


# ==================================================== #
# Computation of isogenies between Kummer lines using  #
# VéluSqrt x-only formula by Bernstein, De Feo, Leroux #
# and Smith                                            #
# ==================================================== #


def product_tree_resultant(hI_tree, poly):
    r"""
    Helper function to evaluate a resultant with `h_I` quickly,
    using the product tree, taken from FastEllipticPolynomial
    sage/src/sage/schemes/elliptic_curves/hom_velusqrt.py

    Original author: Lorenz Panny (2022)
    """
    rems = hI_tree.remainders(poly)
    r = prod(rems)
    s = -1 if len(hI_tree) % 2 == 1 == poly.degree() else 1
    assert r.is_constant()
    return s * r[0]


class KummerLineIsogeny_VeluSqrt(KummerLineIsogeny_Generic):
    """
    VéluSqrt for large ell isogenies

    https://velusqrt.isogeny.org/
    Faster computation of isogenies of large prime degree
    Daniel J. Bernstein, Luca De Feo, Antonin Leroux, Benjamin Smith

    TODO: currently seems to be under-performing. I think there are
    further optimisations which can be made. Trying to use more projective
    points to remove inversions just seemed to slow down the polynomial
    computations though.
    """

    def __init__(self, domain, kernel, degree, check=True):
        # Check the input to the isogeny is well-formed
        self.validate_input(domain, kernel, degree, check=check)

        # Set kernel and degree and domain
        self._degree = degree
        self._kernel = kernel
        self._domain = domain

        # We need the domain coefficient for the elliptic
        # resultants.
        self.a = self._domain.a()

        # We need a polynomial ring, so we create it once
        # and store it to self
        k = self._domain.base_ring()
        self.R = PolynomialRing(k, names="Z", implementation="NTL")
        self.Z = self.R.gen()

        # Precompute these values, as of the 30us it takes to compute
        # x + 1, 30us of that is computing R(1) lmao
        self.one = self.R.one()
        self.Ra = self.R(self.a)

        # baby step and giant step params
        b = (self._degree - 1).isqrt() // 2
        c = (self._degree - 1) // (4 * b)
        stop = self._degree - 4 * b * c

        # Pre-compute polynomials which are needed
        # throughout. hI is stored as a product tree
        # for faster resultants
        self.hI_tree = self._hI_precomputation(kernel, b, c)
        self.EJ_parts = self._EJ_precomputation(kernel, b)
        self.hK_data = self._hK_precomputation(kernel, stop)

        # Compute the codomain
        self._codomain = self._compute_codomain()

    def __call__(self, P):
        """
        Evaluate the isogeny phi on the point P
        by using phi(P)
        """
        if not isinstance(P, KummerPoint):
            raise ValueError
        return self._evaluate_isogeny(P)

    def _hI_resultant(self, poly):
        """
        Compute the resultant Res(hI, poly) where
        hI has been computed and stored as a product tree
        """
        return product_tree_resultant(self.hI_tree, poly)

    def _hI_precomputation(self, ker, b, c):
        r"""
        Compute the polynomial

        hI = \Prod (Z - x(Q)) for Q in the set I
        I = {2b(2i + 1) | 0 <= i < c}

        The polynomial is computed using a product tree,
        where the leaves are each factor of the above product
        """
        Q = (b + b) * ker
        step, diff = Q.double(), Q
        leaves = []
        # This uses x-only point addition to generate all points
        # in the set I = {2b(2i + 1) | 0 <= i < c}
        for i in range(c):
            leaves.append(self.Z - Q.x())
            if i < c - 1:
                Q, diff = Q.add(step, diff), Q

        return ProductTree(leaves)

    # def _Fs(self, X1, X2):
    #     """
    #     Elliptic Resultants for Montgomery curves
    #     """
    #     X1X2 = X1 * X2
    #     polys = (
    #         (X1 - X2) ** 2,
    #         -2 * ((X1X2 + 1) * (X1 + X2) + 2 * self.a * X1X2),
    #         (X1X2 - 1) ** 2,
    #     )
    #     return polys

    def _Fs(self, X):
        """
        Elliptic Resultants for Montgomery curves

        Looks janky because sage hates making elements of
        polynomial rings
        """
        # convert to R once, to speed things up
        X = self.R(X)
        Z = self.Z

        # Precompute pieces
        # we want to make new polynomials as little as possible
        # as converting from Fp2 to R is crazy slow
        z1 = Z + X
        z2 = Z - X

        XZ = X * Z
        z3 = XZ + self.one
        z4 = XZ - self.one

        z5 = self.Ra * XZ
        z6 = -(z3 * z1 + z5 + z5)

        polys = (
            z2 * z2,
            z6 + z6,
            z4 * z4,
        )
        return polys

    # def _Fs_projective(self, QX, QZ):
    #     """
    #     Elliptic Resultants for Montgomery curves

    #     Looks janky because sage hates making elements of
    #     polynomial rings
    #     """
    #     # convert to R once, to try and speed things up
    #     QX = self.R(QX)
    #     QZ = self.R(QZ)
    #     Z = self.Z

    #     # Precompute pieces
    #     # we want to make new polynomials as little as possible
    #     # as converting from Fp2 to R is crazy slow
    #     ZQZ = QZ * Z
    #     ZQX = QX * Z

    #     z1p = ZQZ + QX
    #     z2p = ZQX + QZ
    #     z1m = ZQZ - QX
    #     z2m = ZQX - QZ

    #     z4 = self.Ra * QZ
    #     z4 *= ZQX

    #     z6 = -(z1p * z2p + z4 + z4)

    #     polys = (
    #         z1m * z1m,
    #         z6 + z6,
    #         z2m * z2m,
    #     )
    #     return polys

    def _EJ_precomputation(self, ker, b):
        """
        The polynomials for EJ are of the form

        alpha^2 * F0(Z, x(Q)) + alpha * F1(Z, x(Q)) + F2(Z, x(Q))

        For x(Q) in the set J = {1, 3, 5, ..., 2b - 1}

        We cannot precompute the whole polynomial, but we can precompute
        the pieces Fi(Z, x(Q)) and then compute the sum when needed
        """
        Q = ker
        step, diff = Q.double(), Q
        EJ_parts = []
        # This uses x-only point addition to generate all points
        # in the set J = {1, 3, 5, ..., 2b - 1}
        for i in range(b):
            polys = self._Fs(Q.x())
            EJ_parts.append(polys)
            if i < b - 1:
                Q, diff = Q.add(step, diff), Q

        return EJ_parts

    def _hK_precomputation(self, ker, stop):
        r"""
        Compute the polynomial

        hK = \Prod (Z - x(Q)) for Q in the set
        K = {4bc+1, ..., ell-2, ell}
        """
        hK = []
        Q = ker.double()
        step, next_point = Q, Q.double()
        # This uses x-only point addition to generate all points
        # in the set K = {4bc+1, ..., ell-2, ell}
        for i in range(2, stop, 2):
            QX, QZ = Q.XZ()
            hK.append((QX, QZ))

            if i < stop - 1:
                Q, next_point = next_point, next_point.add(step, Q)

        return hK

    def _hK_codomain(self):
        h1, h2 = 1, 1
        for QX, QZ in self.hK_data:
            h1 *= QZ - QX
            h2 *= -(QZ + QX)
        return h1, h2

    def _hK_image(self, alpha):
        h1, h2 = 1, 1
        for QX, QZ in self.hK_data:
            h1 *= QZ - alpha * QX
            h2 *= alpha * QZ - QX
        return h1, h2

    def _compute_codomain_constants(self):
        """
        Compute the codomain constant in projective coordinates
        (A : C) using the VéluSqrt adaptation of the Meyers-Reith
        Twisted Edwards curve trick
        """
        # These are the polynomials for alpha = 1 and alpha = -1
        E0J = prod(F0 + F1 + F2 for F0, F1, F2 in self.EJ_parts)
        E1J = prod(F0 - F1 + F2 for F0, F1, F2 in self.EJ_parts)

        # Compute resultants and evaluate hK at 1 and -1
        R0 = self._hI_resultant(E0J)
        R1 = self._hI_resultant(E1J)
        M0, M1 = self._hK_codomain()

        # We have that
        # d = [(A - 2C)(A + 2C)]^ell * (hS(1) / hS(-1))^8
        # hS = hK R / Delta
        # d = [(A - 2C)(A + 2C)]^ell * (hK R(1) / hK R(-1))^8

        # First compute (hS(1) / hS(-1))^8
        num = R0 * M0
        den = R1 * M1

        # num^8, den^8 with three squares
        num, den = num * num, den * den
        num, den = num * num, den * den
        num, den = num * num, den * den

        # [(A - 2)(A + 2)]^ell
        num *= (self.a - 2) ** self._degree
        den *= (self.a + 2) ** self._degree

        # Compute the new curve y^2 = x^3 + (A:C)x^2 + x
        A_new = num + den
        C_new = den - num  # C =   d - n
        A_new = A_new + A_new  # A = 2(n + d)

        return A_new, C_new

    def _compute_codomain(self):
        """
        Wrapper function to compute the codomain L = x^3 + x^2A' + x in
        projective coordinates: A' = (A' : C')
        """
        A_codomain, C_codomain = self._compute_codomain_constants()
        F = self._domain.base_ring()
        return KummerLine(F, [A_codomain, C_codomain])

    def _evaluate_isogeny(self, P):
        """
        Evaluate the isogeny phi at the point P

        NOTE:

        We're suppose to compute the quotient:

        Res(hI, EJ0(1/alpha)) * hK(1/alpha)
        ----------------------------------- * alpha^ell
          Res(hI, EJ0(alpha)) * hK(alpha)

        But we can use that
            f(1/alpha) = reverse(f(alpha)) * alpha^(-d)
        for a degree d polynomial, where "reverse"
        reverses the coefficients of the polynomial
        to rewrite this quotient as:

        Res(hI, reverse(EJ0(alpha))) * reverse(hK(alpha))
        -------------------------------------------------- * alpha
               Res(hI, EJ0(alpha)) * hK(alpha)
        """

        if P.is_zero():
            return self._codomain((1, 0))

        # x-coordinate of point to evaluate
        alpha = pari(P.x())
        alphaR = self.R(alpha)

        # Compute two polynomials from giant steps
        EJ1 = prod((F0 * alphaR + F1) * alphaR + F2 for F0, F1, F2 in self.EJ_parts)
        EJ0 = EJ1.reverse()

        # Resultants and evaluations
        R0 = self._hI_resultant(EJ0)
        R1 = self._hI_resultant(EJ1)
        M0, M1 = self._hK_image(alpha)

        # Make new point
        R0M0 = R0 * M0
        R1M1 = R1 * M1
        X_new = R0M0 * R0M0 * alpha
        Z_new = R1M1 * R1M1

        return self._codomain((X_new, Z_new))


# =============================================== #
# Compute a composite degree isogeny using x-only #
# isogenies using either Vélu or Vélusqrt         #
#                                                 #
# Note: this code is heavily inspired by the      #
# SageMath classes written by Lorenz Panny        #
# =============================================== #


def evaluate_factored_kummer_isogeny(phi_list, P):
    """
    Given a list of isogenies, evaluates the
    point for each isogeny in the list
    """
    for phi in phi_list:
        P = phi(P)
    return P


def factored_kummer_isogeny(K, P, order, threshold=1000):
    """
    Computes a composite degree isogeny using x-only formula

    - Uses the sparse strategy from the SIDH paper for computing
      prime power degree isogenies
    - Uses VéluSqrt when the prime order isogeny has degree > threshold
    """

    def sparse_isogeny_prime_power(P, l, e, split=0.8, threshold=1000):
        """
        Compute chain of isogenies quotienting
        out a point P of order l**e
        https://trac.sagemath.org/ticket/34239
        """

        if l > threshold:
            KummerLineIsogenyAlgorithm = KummerLineIsogeny_VeluSqrt
        else:
            KummerLineIsogenyAlgorithm = KummerLineIsogeny_Velu

        def recursive_sparse_isogeny(Q, k):
            assert k
            if k == 1:  # base case
                return [KummerLineIsogenyAlgorithm(Q.parent(), Q, l, check=False)]

            k1 = int(k * split + 0.5)
            k1 = max(1, min(k - 1, k1))  # clamp to [1, k-1]

            Q1 = l**k1 * Q
            L = recursive_sparse_isogeny(Q1, k - k1)

            Q2 = evaluate_factored_kummer_isogeny(L, Q)
            R = recursive_sparse_isogeny(Q2, k1)

            return L + R

        return recursive_sparse_isogeny(P, e)

    # Ensure P is a point on E
    if P.parent() != K:
        raise ValueError(f"The supplied kernel must be a point on the line {K}")

    # For computing points
    cofactor = order
    assert (P * order).is_zero()

    # TODO: Deal with isomorphisms
    # Easy option: just use the Sage isomorphisms from K.curve() and map down
    # Better option: just write the isomorphisms of Montgomery curves
    if cofactor == 1:
        raise NotImplementedError(
            "Isomorphisms between Kummer Lines are not yet implemented"
        )

    phi_list = []
    for l, e in cofactor.factor():
        if l < 2 * threshold:
            # Compute point Q of order l^e
            D = ZZ(l**e)
            cofactor //= D

            # Use Q as kernel of degree l^e isogeny
            Q = cofactor * P
            psi_list = sparse_isogeny_prime_power(Q, l, e, threshold=threshold)

            # For the last step, we don't need to put the kernel
            # through the isogeny.
            # TODO: this could be an optional check, to ensure
            # that the image of the kernel is the identity.
            if cofactor != 1:
                P = evaluate_factored_kummer_isogeny(psi_list, P)

            phi_list += psi_list

        # In the above, we take the kernel K with order D and first clear
        # out the cofactor to get a point of order l^e. We then use this
        # as a kernel of prime power order to compute the codomain, pushing
        # this point through the isogeny e-1 times. Then, we take this isogeny
        # chain and push the original kernel K through.
        #
        # Using this method, ee skip multiplying K by the cofactor at each step,
        # at the cost of computing (e-1) additional images. When ell is small this
        # means computing an extra ell xADD but saving log(D / ell) xDBLADD, so this
        # is faster
        #
        # However, when we're in the last few factors of the order, we're using
        # velusqrt, which due to sage slowness with polynomial rings is not very
        # fast, and it's actually faster to halve the number of images we need
        # at the cost of more scalar multiplcations. Ultimately, this is because
        # despite needing only sqrt(ell) additions for the image, the cost of the
        # resultants is so large that log(D/ell) xDBLADD is cheaper.
        else:
            for _ in range(e):
                cofactor //= l
                Q = cofactor * P
                psi = KummerLineIsogeny_VeluSqrt(Q.parent(), Q, l)

                if cofactor != 1:
                    P = psi(P)

                phi_list.append(psi)

    return phi_list


class KummerLineIsogeny(KummerLineIsogeny_Generic):
    """
    Computes composite degree isogenies as a chain of prime
    degree isogenies. Essentially built to emulate
    EllipticCurveHom_composite but using x-only formula
    """

    def __init__(self, domain, kernel, degree, check=True, threshold=1500):
        # Check the input to the isogeny is well-formed
        self.validate_input(domain, kernel, degree, check=check)

        # Compute factored isogeny
        self._phis = factored_kummer_isogeny(
            domain, kernel, degree, threshold=threshold
        )

        # Make immutable
        self._phis = tuple(self._phis)

        # Compute degree, domain and codomain
        self._degree = prod(phi.degree() for phi in self._phis)
        self._domain = self._phis[0].domain()
        self._codomain = self._phis[-1].codomain()

    def __call__(self, P):
        """
        Evaluate the composite isogeny by calling phi(P)
        """
        return evaluate_factored_kummer_isogeny(self._phis, P)

    @classmethod
    def from_factors(cls, maps):
        """
        Sometimes we will have factors of some isogeny from a
        different context and we want to simply collect them
        together to create a single object.

        Built following the same classmethod which appears in
        EllipticCurveHom_composite
        """
        maps = tuple(maps)

        L = maps[0].domain()
        for phi in maps:
            if not isinstance(phi, KummerLineIsogeny_Generic):
                raise TypeError(f"not an kummer-line isogeny: {phi}")
            if phi.domain() != L:
                raise ValueError(f"isogeny has incorrect domain: {phi}")
            L = phi.codomain()

        result = cls.__new__(cls)

        # Make immutable
        result._phis = maps

        # Compute degree, domain and codomain
        result._degree = prod(phi.degree() for phi in result._phis)
        result._domain = result._phis[0].domain()
        result._codomain = result._phis[-1].codomain()

        return result