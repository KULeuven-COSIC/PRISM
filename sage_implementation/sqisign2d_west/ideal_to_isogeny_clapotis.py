import time
import random

from sage.all import (
    ZZ,
    ceil,
    inverse_mod,
    gcd,
    randint,
    cached_method
)

from utilities.sum_of_squares import sum_of_squares
from quaternions.ideals import (
    reduced_basis,
    enumerate_basis
)
from sqisign2d_west.fixed_degree_isogeny import FixedDegreeIsogeny
from sqisign2d_west.endomorphisms import iota
from precomputations import Precomputations
from sqisign2d_west.dim2_wrapper import Dim2Isogeny
from utilities.logging_helper import logging

logger = logging.getLogger(__name__)

class IdealToIsogenyClapotis:
    def __init__(self, I, precomps=None):
        if precomps is None:
            _, j, _ = I.quaternion_algebra().gens()
            p = ZZ(-j**2)
            precomps = Precomputations(p)
        self.precomps = precomps

        self.representation = None
        self.O0 = self.precomps.O0
        self.e_target = precomps.e-self.precomps.extra_torsion_clapotis
        self.I = I
        self.N = ZZ(I.norm())
        self.sums = precomps.sums
        self._count = (None, None)

    def short_equivalent_ideals(self):
        """
        Given an ideal I of norm N, we find two equivalent ideals I1, I2 of coprime odd norms d1, d2
        Also, we output u and v such that u d1 + v d2 = 2 ^ f for some f>0
        """
        I = self.I
        N = self.N
        O0 = self.O0
        sums = self.sums
        init_time = time.time()
        precomps = self.precomps

        def find_u_v(d1, d2):
            #logger.debug(f"u,v search for {d1}, {d2}")
            e_target = self.e_target
            # f1 = d1.nbits()
            # f2 = d2.nbits()
            # f = max(f1, f2)
            d1_inv_mod_d2 = inverse_mod(d1, d2)
            u = (2**e_target * d1_inv_mod_d2) % d2
            v = (2**e_target - u * d1) // d2
            attempts = 0
            while u * d1 < 2 ** e_target and attempts < 100:
                attempts += 1
                val_u=u.valuation(2)
                val_v=v.valuation(2)
                min_val=min(val_u, val_v)
                common=2**min_val
                uu=u//common
                vv=v//common
                logger.debug3(f"u={u}, v={v}; min_val={min_val}")

                # this condition should be automatic: since we divide by the common
                # power two factor, one of u, v is odd; but then since d1,
                # d2 are odd and 2**e even, the other has to be odd too.
                if uu%2 == 1 and vv%2 == 1:
                    logger.debug3(f"u={uu}, v={vv} found in {attempts} attempts")
                    if not sums:
                        assert uu*d1+vv*d2 == 2**(e_target-min_val)
                        return uu, vv, e_target-min_val

                    # look for one sums of two squares
                    x_y_u = sum_of_squares(uu)
                    if x_y_u:
                        logger.debug3(f"u={uu}, v={vv} found in {attempts} attempts")
                        uu = x_y_u
                        return uu, vv, e_target-min_val, 0

                    x_y_v = sum_of_squares(v)
                    if x_y_v:
                        logger.debug3(f"u={uu}, v={vv} found in {attempts} attempts")
                        vv = x_y_v
                        return uu, vv, e_target-min_val, 1

                u += d2
                v -= d1

            return False

        def finding_ideals_of_right_norm(element_list):
            for α1, d1, α2, d2 in element_list:
                if gcd(d1, d2) == 1:
                    output = find_u_v(d1, d2)
                    if output:
                        if sums:
                            u, v, f, flag = output
                            # exhanging the roles of u and v
                            # so u is always a sum of two squares
                            if flag == 1:
                                α1p, d1p = α1, d1
                                α2p, d2p = α2, d2
                                α1, d1 = α2p, d2p
                                α2, d2 = α1p, d1p
                                x, y = v
                                v = u
                                u = x ** 2 + y ** 2
                            else:
                                x, y = u
                                u = x ** 2 + y ** 2
                        else:
                            u, v, f = output

                        #if (u % 2 == 0) or (v % 2 == 0):
                        #    continue

                        θ = α2 * α1.conjugate() / N
                        #print(f"Finding short equivalent ideals required {counter} attempts!")
                        if sums:
                            return α1, d1, d2, [ZZ(x), ZZ(y)], v, f, θ
                        else:
                            return α1, d1, d2, u, v, f, θ
            #print(f"Finding short equivalent ideals failed in {counter} attempts!")
            return False

        γ1, γ2, γ3, γ4 = reduced_basis(I)

        small_elements = []
        small_norms = []
        up_bound = ceil(self.precomps.p ** (0.75))
        all_tuples = []
        count = 0; full_count = 0

        def add_new_element(element):
            nonlocal count
            count += 1
            α, d = element
            small_elements.append(element)
            new_tuples = [(α2, d2, α, d) for (α2,d2) in small_elements]
            all_tuples.extend(new_tuples)
            return finding_ideals_of_right_norm(new_tuples)

        if precomps.enumerate_basis_clapotis_random:
            # temporary bound on the elements we are considering
            # experimenting with skew sets
            if sums:
                bound_m1 = [t for t in range(-4, 5, 1)]
                bound_m2 = [t for t in range(-4, 5, 1)]
                bound_m3 = [t for t in range(-4, 5, 1)]
                bound_m4 = [t for t in range(-4, 5, 1)]
            else:
                bound_m1 = [t for t in range(-2, 3, 1)]
                bound_m2 = [t for t in range(-2, 3, 1)]
                bound_m3 = [t for t in range(-2, 3, 1)]
                bound_m4 = [t for t in range(-2, 3, 1)]

            indices = [(i1, i2, i3, i4)
                        for i1 in bound_m1
                        for i2 in bound_m2
                        for i3 in bound_m3
                        for i4 in bound_m4
                       ]
            for i1, i2, i3, i4 in indices:
                α = i1 * γ1 + i2 * γ2 + i3 * γ3  + i4 * γ4
                d = ZZ(α.reduced_norm()) // N
                full_count += 1
                if not (d in small_norms or d%2 == 0 or d > up_bound):
                    output = add_new_element((α,d))
                    if output:
                        self._count=(count,full_count)
                        logger.debug(f"Clapotis succeeded in {count}/{full_count} samples in {(time.time() - init_time):.3f}s to build an isogeny of degree 2^{output[5]}")
                        return output
            while True:
                d = 0
                while (d in small_norms or d%2 == 0 or d > up_bound):
                    full_count += 1
                    i1, i2, i3, i4 = randint(-10,10), randint(-10,10), randint(-10,10), randint(-10,10)
                    α = i1 * γ1 + i2 * γ2 + i3 * γ3  + i4 * γ4
                    d = ZZ(α.reduced_norm()) // N
                output = add_new_element((α,d))
                if output:
                    self._count=(count,full_count)
                    logger.debug(f"Clapotis succeeded in {count}/{full_count} samples by enlarging the bounds in {(time.time() - init_time):.3f}s to build an isogeny of degree 2^{output[5]}")
                    return output
        else:
            for i1, i2, i3, i4 in enumerate_basis(bound=precomps.enumerate_basis_clapotis_bound):
                α = i1 * γ1 + i2 * γ2 + i3 * γ3  + i4 * γ4
                d = ZZ(α.reduced_norm()) // N
                full_count += 1
                if not (d in small_norms or d%2 == 0 or d > up_bound):
                    output = add_new_element((α,d))
                    if output:
                        self._count=(count,full_count)
                        logger.debug(f"Clapotis succeeded in {count}/{full_count} samples in {(time.time() - init_time):.3f}s to build an isogeny of degree 2^{output[5]}")
                        return output
            raise RuntimeError("Clapotis: failed to find a representation")

    def get_representation(self):
        """
        Given an ideal I of norm N, compute the data needed for the
        isogeny. This function only works on the quaternion side.
        """
        if self.representation is None:
            self.representation = self.short_equivalent_ideals()
        return self.representation

    @cached_method
    def images(self):
        """
        Given an ideal I of norm N, this function computes the isogeny
        φI : E0 -> EI.
        Given the 2^e-torsion canonical basis <P0, Q0>
        we compute φI(P0) and φI(Q0)
        """

        precomps = self.precomps
        e=precomps.e
        P0,Q0=precomps.basis
        e_target = self.e_target
        sums = self.sums
        E0 = precomps.E0

        # let's find two equivalent ideals of short norm d1 and d2
        # also verifying u * d1 + v * d2 == 2 ** f
        α1, d1, d2, u, v, f, θ = self.get_representation()

        if sums:
            Pu = u[0] * P0 + u[1] * iota(P0)
            Qu = u[0] * Q0 + u[1] * iota(Q0)
            u = u[0]**2 + u[1]**2
            Eu = E0
        else:
            φu = FixedDegreeIsogeny(u, self.precomps)
            Eu, Pu, Qu = φu.images()

        # let's compute φv
        φv = FixedDegreeIsogeny(v, self.precomps)
        Ev, Pv, Qv = φv.images()

        # let's evaluate θ
        ## Pu_dual, Qu_dual = u * P0, u * Q0
        ## θP, θQ = θ(Pu_dual, Qu_dual))
        #(P0, Q0) θm = θ(Pu_dual, Qu_dual)
        θm = u * precomps.matrix_theta(θ)

        # let's compute and evaluate φv and evaluate
        # P_, Q_ = imv * θm
        P_, Q_ = self.precomps.evaluate_matrix(θm, (Pv, Qv))

        # let ψ = φv ˚ θ ˚ φu~

        ## ToOptimize: we could precompute the 2**i P0, Q0 and push these through
        # let's prepare the kernel for the 2 dimensional isogeny
        c = e_target - f
        P1 = 2**c * (d1 * u) * Pu
        P2 = 2**c * P_
        Q1 = 2**c * (d1 * u) * Qu
        Q2 = 2**c * Q_

        ker_Phi = ((P1, P2), (Q1, Q2))
        Phi = Dim2Isogeny(ker_Phi, f, extra = self.precomps.extra_torsion_clapotis, precomps=self.precomps)
        EI, LP, LQ = Phi.embedded_isogeny(basis=(Pu, Qu), degree=u*d1, e=e, ePQ=precomps.eP0Q0(e)**u)
        u_inv = inverse_mod(u, 2 ** e)
        PI1 = u_inv * LP
        QI1 = u_inv * LQ

        # let's compute I * I1.conjugate()
        d1_inv = inverse_mod(d1, 2 ** e )
        θI = precomps.matrix_theta(α1)
        PI, QI = self.precomps.evaluate_matrix(θI, (PI1, QI1))
        PI = d1_inv * PI; QI = d1_inv * QI

        return EI, PI, QI

    # this function is not used
    def __call__(self, P):
        _, PI, QI = self.images()
        x, y = self.precomps.canonical_BiDLP(P)
        return x * PI + y * QI
