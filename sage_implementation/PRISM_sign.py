from sage.all import (
    is_prime,
    is_even,
    matrix,
    ZZ,
    Zmod,
    is_pseudoprime,
    GF,
    cached_method,
)

import hashlib
import os

import params
from precomputations import Precomputations
from quaternions.ideals import pushforward_ideal
from sqisign2d_west.ideal_to_isogeny_clapotis import IdealToIsogenyClapotis
from sqisign2d_west.dim2_wrapper import Dim2Isogeny
from utilities.discrete_log import discrete_log_pari
from utilities.pairing import weil_pairing_biextension, weil_pairing_pari
from utilities.supersingular import torsion_basis_2e

def HashToPrime(msg, pk_curve, counter=0):
    """
    Hash (message || pk || counter) into a a-1 bit number
    increasing counter until the output is prime.
    """
    msg += str(pk_curve.j_invariant()).encode()

    while True:
        h = hashlib.sha256(msg + str(counter).encode()).digest()
        q = int.from_bytes(h, 'big') % 2**(params.sign_torsion)
        if is_pseudoprime(q):
            return q, counter
        counter += 1

class PRISM_sign:

    def __init__(self, role, precomps=None):
        """
        Initialize the protocol precomputations.
        """
        p = params.p
        to_patch = [GF(2), GF(3), GF(3**2)]
        for x in to_patch:
            type(x).vector_space = cached_method(type(x).vector_space)
        if not precomps:
            precomps = Precomputations(p, role = role, a = 'max')
        self.precomps = precomps

    def keygen(self):
        """
        Generate public and secret keys
        """
        # Generation of a random ideal Isk of norm N
        Isk = self.precomps.random_ideal()
        self._precompute_secret_key(Isk)

    def sign(self, msg):
        """
        Sign a given message
        """
        # Hash the message to a prime
        q, counter = HashToPrime(msg, self.pk_curve)
        assert q < 2**params.sign_torsion

        Isk, Mpk = self.sk
        precomps = self.precomps

        # Generate and pushforward the ideal of given norm
        Ir_norm = q * (2**params.sign_torsion - q)
        Ir_prime = precomps.random_ideal_of_given_norm(Ir_norm, prime=False)

        Ir = pushforward_ideal(Ir_prime, Isk)
        Iaux = Isk * Ir

        # Compute the isogeny
        phi_aux = IdealToIsogenyClapotis(Iaux, precomps = precomps)
        Eaux, Paux, Qaux = phi_aux.images()
        Pr, Qr = self.precomps.evaluate_matrix(Mpk, (Paux, Qaux))

        # We only need smaller torsion
        torsion_delta = 2**(precomps.a - params.sign_torsion)
        Pr *= torsion_delta
        Qr *= torsion_delta

        # Return the points
        if not params.point_compression_sign:
            return Eaux, (Pr, Qr), counter

        # Point compression
        ind, a_, b_, c_ = self._point_compression(Eaux, Pr, Qr)
        return Eaux, (ind, a_, b_, c_), counter

    def verify(self, msg, sig, pk):
        """
        Given the kernel of an higher dimensional isogeny verify
        that it splits and that the degree is correct.
        """
        precomps = self.precomps

        Er, pts, counter = sig
        q, _ = HashToPrime(msg, pk[0], counter = counter)

        Ppk, Qpk = pk[1]

        P = ZZ(q) * Ppk
        Q = ZZ(q) * Qpk

        if params.point_compression_sign:
            ind, a_, b_, c_ = pts
            Pr, Qr = self._point_decompression(ind, a_, b_, c_, Er, P, Q)

        else:
            Pr, Qr = pts

        # Check response representation
        ker_Phi = ((P, Pr), (Q, Qr))

        Phi = Dim2Isogeny(ker_Phi, params.sign_torsion, precomps = precomps)

        # Check degrees using pairings à la SQIsign2D-East
        idr = Er(0)
        P1, P2 = Phi((Ppk, idr))
        Q1, Q2 = Phi((Qpk, idr))

        D = 2**params.sign_torsion

        # Pari seems to be faster
        W = pk[2]
        W1 = weil_pairing_pari(P1, Q1, D)
        Wq = W ** q
        Wqinv = Wq ** (-1)
        return W1 in [Wq, Wqinv]

    def _precompute_secret_key(self, Isk):
        """
        Various precomputations related to the secret key
        """

        precomps=self.precomps
        a = precomps.e

        # Computing the isogeny associated with the ideal Isk
        phi_sk = IdealToIsogenyClapotis(Isk, precomps = self.precomps)
        Epk, phi_sk_P0, phi_sk_Q0 = phi_sk.images()

        self._precompute_public_key(Epk)
        Ppk, Qpk = self.pk_basis

        Mpk = matrix(
            Zmod(2**a), precomps.BiDLPs((Ppk, Qpk), phi_sk_P0, phi_sk_Q0, a)
        )
        self.sk = Isk, Mpk

        return None

    def _precompute_public_key(self, Epk):
        """
        Various precomputations related to the public key
        """
        precomps = self.precomps

        self.pk_curve = Epk
        Ppk, Qpk = torsion_basis_2e(Epk, precomps.e, montgomery=True)
        self.pk_basis = (Ppk, Qpk)

        torsion_delta = 2**(precomps.a - params.sign_torsion)
        Ppk *= torsion_delta
        Qpk *= torsion_delta
        self.pk_sign_basis = (Ppk, Qpk)
        self.eWPpkQpk = weil_pairing_pari(Ppk, Qpk, 2**params.sign_torsion)

        return None

    def _point_compression(self, E, P, Q):
        """
        Given 2 points P, Q on E represent them with 3 scalars.
        Let E[2^e] = <P0, Q0> and let
            P = a1 * P0 + a2 * Q0
            Q = a3 * P0 + a4 * Q0
        It returns
        - (0, a1, a2, a3) if a1%2 == 1
        - (1, a1, a3, a4) if a1%2 == 0
        """
        e = params.sign_torsion
        P0, Q0 = torsion_basis_2e(E, e, montgomery=True)
        (a1, a2), (a3, a4) = self.precomps.BiDLPs((P, Q), P0, Q0, e)
        if a1%2 == 1:
            return 0, a1 , a2, a3
        else:
            return 1, a1 , a3, a4

    def _point_decompression(self, ind, a_, b_, c_, Er, qPpk, qQpk):
        """
        Decompression: exploit the fact that the kernel of a
        two-dimensional isogeny is isotropic, i.e. if the
        kernel is <(P1, P2), (Q1, Q2)> then
            wp(P1, Q1)*wp(P2, Q2) == 1
        Hence to recover Pr, Qr we need q*Ppk, q*Qpk.
        """
        e = params.sign_torsion
        P0, Q0 = torsion_basis_2e(Er, e, montgomery=True)

        e0 = weil_pairing_biextension(P0, Q0, e)
        ePk = weil_pairing_biextension(qPpk, qQpk, e)

        M = discrete_log_pari(er, ePk, 2**e)
        Z_mod = Zmod(2**e)
        M_inv = Z_mod(M) ** (-1)
        if ind == 0:
            a1, a2, a3 = a_, b_, c_
            a1_inv =  Z_mod(a1) ** (-1)
            a4 = Z_mod(a2 * a3 - M_inv) * a1_inv
            Pout = a1 * Pr + a2 * Qr
            Qout = a3 * Pr + a4 * Qr
        else:
            a1, a3, a4 = a_, b_, c_
            a3_inv =  Z_mod(a3) ** (-1)
            a2 = Z_mod(M_inv + a1 * a4) * a3_inv
            Pout = a1 * Pr + a2 * Qr
            Qout = a3 * Pr + a4 * Qr
        return Pout, Qout

