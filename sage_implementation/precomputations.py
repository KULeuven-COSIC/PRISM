import json

from sage.all import(
    EllipticCurve,
    GF, ZZ, Zmod,
    matrix, vector,
    QuaternionAlgebra,
    randint, is_prime,
    gcd
)

from sqisign2d_west.endomorphisms import FullrepresentInteger, iota
from utilities.supersingular import (
        torsion_basis_with_pairing_2e, torsion_basis_with_pairing_3f)
from utilities.discrete_log import discrete_log_pair_power_two
from utilities.pairing import (
        weil_pairing_pari, tate_pairing_pari, weil_pairing_biextension,
        tate_pairing_biextension)
from utilities.strategy import optimised_strategy
import params
from montgomery_isogenies.kummer_line import KummerLine


class Precomputations:
    def __init__(self, p, role='signer', a = 'max'):
        self.p = ZZ(p)
        self.c = c
        e = (self.p + 1).valuation(2)
        if a == 'max':
            a = e
        assert a in ZZ and a <= e and a > 0
        self.a = a
        D = ZZ(2) ** a
        self.e = e
        self.D = D
        self.tate_exp = (p**2 - 1) // D
        self.aux_integer = D #auxiliary integer for random_ideal_of_given_norm

        ## Isogenies:
        # Do we need extra available torsion for our isogenies?
        # Possible choices: 0, 1 (-> 2 extra sqrts, not implemented), 2 (-> 5 extra sqrts)
        self.extra_torsion_clapotis = 0
        self.extra_torsion_fixed_degree = 0

        ## Pairings and DLPs:
        self.use_tate = True #use Tate for full torsion DLPs?
        self.use_biextension = False #use biextension pairings?
        self.window = [10] #windows for the 2**e finite field DLPs

        if role == 'signer':
            self.precomputations_signer(self.p)
        elif role == 'signer_SIDH':
            self.precomputations_signer_SIDH(self.p)

        self.precompute_strategies()

    def precomputations_signer(self, p):

        self.random_ideal_has_prime_norm = True
        ## Clapotis:
        self.sums = False #try to find a sum of 2 squares in IdealToIsogenyClapotis?

        ## FixedDegreeIsogeny:
        self.fixed_degree_isogeny_truly_random = False # Do we allow our dim 2 isogeny from FixedDegreeIsogeny to start with endomorphisms and diagonal isogenies?

        ## How to enumerate basis of ideals
        self.enumerate_basis_random = False
        self.enumerate_basis_random_bound = 5
        self.enumerate_basis_random_max_attempts = 10000
        self.enumerate_basis_bound = 10 #for enumerate_basis_random = False
        self.enumerate_basis_clapotis_random = self.enumerate_basis_random
        self.enumerate_basis_clapotis_bound = self.enumerate_basis_bound

        a = self.a
        D = self.D

        B = QuaternionAlgebra(-1, -p)
        self.B = B
        i, j, k = B.gens()
        self.O0 = B.quaternion_order([1, i, (i + j) / 2, (1 + k) / 2])

        F = GF(p**2, name="i", modulus=[1, 0, 1])
        E0 = EllipticCurve(F, [1, 0])
        E0.set_order((p + 1) ** 2)
        self.E0 = E0

        # Sage tricks to precompute "fast"
        # Proper implementation won't have them
        K = F.extension(2, name="T")
        E0_K = E0.base_extend(K)
        self._E0_K = E0_K
        i = F.gen()
        assert i ** 2 == -1
        self.F = F
        self.sqrt_minus_one = i

        P0, Q0, eWP0Q0 = torsion_basis_with_pairing_2e(E0, a)
        self.basis = (P0, Q0)

        self.eWP0Q0 = eWP0Q0
        self.eTP0Q0 = self.tate_pairing(P0, Q0, a)

        P0_K, Q0_K = E0_K(P0), E0_K(Q0)

        P0_p = P0_K.division_points(2)[0]
        Q0_p = Q0_K.division_points(2)[0]

        P0_ij2 = E0(E0_K(-P0_p[0], i * P0_p[1]) + E0_K(P0_p[0] ** p, P0_p[1] ** p))
        Q0_ij2 = E0(E0_K(-Q0_p[0], i * Q0_p[1]) + E0_K(Q0_p[0] ** p, Q0_p[1] ** p))

        P0_1k2 = E0(P0_p + E0_K(-P0_p[0] ** p, i * P0_p[1] ** p))
        Q0_1k2 = E0(Q0_p + E0_K(-Q0_p[0] ** p, i * Q0_p[1] ** p))

        self.endomorphism_images = (P0_ij2, Q0_ij2, P0_1k2, Q0_1k2)
        self.matrix_id = matrix(Zmod(D), ((1,0),(0,1)))
        self.matrix_iota = matrix(Zmod(D), self.canonical_BiDLPs((iota(P0), iota(Q0))))
        self.matrix_ij2 = matrix(self.canonical_BiDLPs((P0_ij2, Q0_ij2)))
        self.matrix_1k2 = matrix(self.canonical_BiDLPs((P0_1k2, Q0_1k2)))

    def precomputations_signer_SIDH(self, p):
        assert p % 4 == 3
        assert p.is_prime()
        f = (p + 1).valuation(3)
        self.f = f
        self.odd_torsion = ZZ(3 ** self.f)

        self.precomputations_signer(p)

        P3, Q3, eWP3Q3 = torsion_basis_with_pairing_3f(self.E0, f)
        self.P3 = P3
        self.Q3 = Q3
        self.PQ3 = P3 - Q3
        self.eWP3Q3 = eWP3Q3

        ## Kummer Line stuff
        self._E0 = KummerLine(self.E0)
        self._P3 = self._E0(P3[0])
        self._Q3 = self._E0(Q3[0])
        self._PQ3 = self._E0(self.PQ3[0])

    def precompute_strategies(self):
        # https://stackoverflow.com/questions/1450957/pythons-json-module-converts-int-dictionary-keys-to-strings
        def jsonKeys2int(x):
            if isinstance(x, dict):
                return {int(k):v for k,v in x.items()}
            return x
        # TODO: clean
        import pathlib
        strategy_path = pathlib.Path(__file__).parent.resolve() / 'helpers' / 'strategies_dim2.json'
        self.dim2_strats = json.load( open(strategy_path), object_hook=jsonKeys2int )

    def random_degree(self, prime=None, small=False):
        """
        Return a prime random degree suitable to use for the key
        generation or the commitment (where we will then randomize an ideal
        class of this degree)
        """

        quotient=1
        if prime is None:
            prime = self.random_ideal_has_prime_norm
        if small: # used for fast_commitment
            # In FixedDegreeIsogeny we use: `s = ceil(p/u +u).nbits() + 25`
            # so here we divide by 2**30 to be sure that s<=e
            quotient = 2**30
        # TODO: work out an upper bound
        N = 2 * randint(0, self.p//quotient) + 1
        if prime:
            while not is_prime(N):
                N = 2 * randint(0, self.p//quotient) + 1
        return N

    def random_ideal_of_given_norm(self, N, randomness = True, prime = None):
        """
        Given an odd prime N, it generates a uniformly distributed ideal I of norm N
        - prime: set to True if N is known to be prime to save a gcd
        """
        if prime is None:
            prime = self.random_ideal_has_prime_norm

        B = self.B; O0 = self.O0
        i, j, k = B.gens()

        a, b, c, d = FullrepresentInteger(N * self.aux_integer, self.p, primitive = True)
        γ = a + b * i + c * j + d * k

        if not randomness:
            I = O0 * γ + O0 * N
            assert ZZ(I.norm()) == N
            return I

        α1, α2, α3, α4 = O0.basis()
        while True:
            u1, u2, u3, u4 = randint(0, N), randint(0, N), randint(0, N), randint(0, N)
            α = u1 * α1 + u2 * α2  + u3 * α3 +  u4 * α4
            if prime and ZZ(α.reduced_norm()) % N != 0:
                break
            if (not prime) and gcd(ZZ(α.reduced_norm()), N) == 1:
                break

        I = O0 * (γ * α) + O0 * N
        assert ZZ(I.norm()) == N
        return I

    def random_ideal(self, randomness=True, prime=None, small=False):
        """
        Return a random ideal starting from E0
        """
        N = self.random_degree(prime=prime, small=small)
        I = self.random_ideal_of_given_norm(N, randomness=randomness, prime=prime)
        return I

    def matrix_theta(self, θ):
        """
        Given an endomorphism θ, return the matrix of θ on the canonical
        basis P0, Q0
        """
        x1, x2, x3, x4 = θ.coefficient_tuple()
        c1 = ZZ(x1  - x4)
        c2 = ZZ(x2 - x3)
        c3 = ZZ(2 * x3)
        c4 = ZZ(2 * x4)
        return c1 * self.matrix_id + c2 * self.matrix_iota + c3 * self.matrix_ij2 + c4 * self.matrix_1k2

    def evaluate_matrix(self, M, basis = None):
        if basis is None:
            basis = self.basis
        return (M[0,0] * basis[0] + M[0,1] * basis[1],
                M[1,0] * basis[0] + M[1,1] * basis[1])

    def BiDLP(self, R, P, Q, e=None, ePQ=None):
        '''
            computes a,b such that R = aP + bQ where
            P, Q are a 2^e torision basis.
        '''
        if e is None:
            e = self.e
        if ePQ:
            pair_PQ = ePQ
        else:
            pair_PQ = self.pairing(P, Q, e)

        pair_a = self.pairing(R, Q, e)
        pair_b = self.pairing(R, -P, e)
        a, b = discrete_log_pair_power_two(pair_a, pair_b, pair_PQ, e, window=self.window)
        return a, b

    def BiDLPs(self, Rs, P, Q, e=None, ePQ=None):
        if e is None:
            e = self.e
        if ePQ:
            pair_PQ = ePQ
        else:
            pair_PQ = self.pairing(P, Q, e)
        tate = False
        return (self.BiDLP(R, P, Q, e, ePQ=ePQ) for R in Rs)

    def canonical_BiDLP(self, P):
        """
        Return the vector BiDLP with respect to the canonical basis of 2^a-torsion
        Mostly for testing purposes, it's better to work directly with matrices
        """
        P0, Q0 = self.basis

        # many times we will evaluate endomorphisms on P0 or Q0
        # Trick to skip the DLOGs
        if P == P0:
            return vector((1, 0))

        if P == Q0:
            return vector((0, 1))

        ePQ = self.eP0Q0(self.a)
        x, y = self.BiDLP(P, P0, Q0, self.a, ePQ=ePQ)
        return vector((x, y))

    def canonical_BiDLPs(self, Ps):
        return (self.canonical_BiDLP(P) for P in Ps)

    def strategy_dim2(self, n):
        if n in self.dim2_strats:
            return self.dim2_strats[n]
        else:
            logger.warning(f"Dim 2 strategy not precomputed: {n}")
            return optimised_strategy(n)

    def eP0Q0(self, e):
        if self.use_tate and e == self.e:
            return self.eTP0Q0
        else:
            return self.eWP0Q0

    #P, Q of 2^e torsion
    def weil_pairing(self, P, Q, e):
        if self.use_biextension:
            return weil_pairing_biextension(P, Q, e)
        else:
            return weil_pairing_pari(P, Q, 2**e)

    def tate_pairing(self, P, Q, e, exp=None):
        if exp is None:
            exp = self.tate_exp
        if self.use_biextension:
            return tate_pairing_biextension(P, Q, e, exp=exp)
        else:
            return tate_pairing_pari(P, Q, 2**e)**exp

    def pairing(self, P, Q, e):
        if self.use_tate and e == self.e:
            return self.tate_pairing(P, Q, e)
        else:
            return self.weil_pairing(P, Q, e)

