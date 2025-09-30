from sage.all import (
    QuaternionAlgebra,
    ceil, ZZ,
    cached_method
)

from theta_structures.couple_point import CouplePoint
from sqisign2d_west.dim2_wrapper import Dim2Isogeny
from sqisign2d_west.endomorphisms import (
    check_endomorphism,
    evaluate_endomorphism_type,
    check_diagonal_isogeny,
    diagonal_isogeny,
    FullrepresentInteger
)

class FixedDegreeIsogeny:
    def __init__(self, u, precomps):
        p=precomps.p
        e=precomps.e
        self.precomps = precomps
        self.TrulyRandom = precomps.fixed_degree_isogeny_truly_random
        self.e_target = e - self.precomps.extra_torsion_fixed_degree

        E0 = precomps.E0
        P0, Q0 = precomps.basis

        # TODO: choice?
        s = ceil(p/u + u).nbits() + 25
        s = min(s, self.e_target)
        a, b, c, d = FullrepresentInteger(u * (2 ** s - u), p, TrulyRandom = self.TrulyRandom)

        self.s = s
        self.u = u

        B = self.precomps.B
        i, j, k = B.gens()
        α = a + b * i + c * j + d * k
        Mα = self.precomps.matrix_theta(α)
        P0_, Q0_ = 2 ** (self.e_target -s) * P0, 2 ** (self.e_target -s) * Q0
        (alpha_P0, alpha_Q0) = precomps.evaluate_matrix(Mα, (P0_, Q0_))

        self.α = α
        self.alpha_P0 = alpha_P0
        self.alpha_Q0 = alpha_Q0

        # First compute the pairs of images from P0, Q0
        # using the endomorphism alpha
        P1, P2 = self.u * P0_, self.alpha_P0
        Q1, Q2 = self.u * Q0_, self.alpha_Q0
        self.ker_Phi = ((P1, P2), (Q1, Q2))

        self._Phi = None #compute when needed
        #self.Phi = self._compute_22_chain(P0_, Q0_)

    def Phi(self):
        """
        Compute the dimension 2 isogeny embedding the fixed degree isogeny.
        """
        if self._Phi is None:
            self._Phi = self._compute_22_chain()
        return self._Phi

    def ideal(self):
        O0=self.precomps.O0
        return O0 * self.u + O0 * self.α

    def _compute_22_chain(self):
        """TODO"""
        (P1, P2), (Q1, Q2) = self.ker_Phi

        if not self.TrulyRandom:
            ker_Phi = ((P1, P2), (Q1, Q2))
            Phi = Dim2Isogeny(ker_Phi, self.s, extra = self.precomps.extra_torsion_fixed_degree, precomps=self.precomps)
            return Phi

        else:
            # let's compute the computing the isogeny Φ with kernel
            # <(P1, P2), (Q1, Q2)>
            index = 0
            P1_, P2_ = 2 ** (self.s + 1) * P1, 2 ** (self.s + 1) * P2
            Q1_, Q2_ = 2 ** (self.s + 1) * Q1, 2 ** (self.s + 1) * Q2
            endomorphisms = []
            diagonal_isogenies = []
            type_end = check_endomorphism(P1_, P2_, Q1_, Q2_)

            # We need to remember the endomorphisms for images
            if type_end:
                endomorphisms.append(type_end)

            while type_end:
                index += 1
                P1, P2 = evaluate_endomorphism_type(type_end, P1, P2)
                Q1, Q2 = evaluate_endomorphism_type(type_end, Q1, Q2)
                P1_, P2_ = 2 ** (self.s + 1 - index) * P1, 2 ** (self.s + 1 - index) * P2
                Q1_, Q2_ = 2 ** (self.s + 1 - index) * Q1, 2 ** (self.s + 1 - index) * Q2
                type_end = check_endomorphism(P1_, P2_, Q1_, Q2_)
                # We need to remember the endomorphisms for images
                if type_end:
                    endomorphisms.append(type_end)

            # checking for diagonal isogenies
            flag, R1, R2 = check_diagonal_isogeny(P1_, P2_, Q1_, Q2_)
            while flag:
                index += 1
                φ1, φ2 = diagonal_isogeny(R1, R2)
                F1 = φ1.codomain()
                F2 = φ2.codomain()
                P1, P2 = φ1(P1), φ2(P2)
                Q1, Q2 = φ1(Q1), φ2(Q2)
                temp = [φ1, φ2]
                diagonal_isogenies.append(temp)
                P1_, P2_ = 2 ** (self.s + 1 - index) * P1, 2 ** (self.s + 1 - index) * P2
                Q1_, Q2_ = 2 ** (self.s + 1 - index) * Q1, 2 ** (self.s + 1 - index) * Q2
                flag, R1, R2 = check_diagonal_isogeny(P1_, P2_, Q1_, Q2_)

            # Compute the kernel
            ker_Phi = ((P1, P2), (Q1, Q2))
            Phi = Dim2Isogeny(ker_Phi, self.s - index, extra = self.precomps.extra_torsion_fixed_degree, endomorphisms=endomorphisms, diagonal_isogen=diagonal_isogenies)
            return Phi

    @cached_method
    def images(self):
        return self.Phi().embedded_isogeny(degree=self.u)

    def __call__(self, P):
        _, PI, QI = self.images()
        x, y = self.precomps.canonical_BiDLP(P)
        return x * PI + y * QI
