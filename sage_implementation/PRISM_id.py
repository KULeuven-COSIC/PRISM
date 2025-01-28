import time
import hashlib
import os

from sage.all import (
    is_prime,
    is_even,
    matrix,
    ZZ,
    Zmod,
    GF,
    cached_method,
)

import params
from utilities.supersingular import torsion_basis_2e
from precomputations import Precomputations
from sqisign2d_west.ideal_to_isogeny_clapotis import IdealToIsogenyClapotis
from quaternions.ideals import pushforward_ideal
from sqisign2d_west.dim2_wrapper import Dim2Isogeny
from utilities.pairing import weil_pairing_pari


class PRISM_id:

    def __init__(self, role, precomps=None):
        """
        Initialize the protocol precomputations.
        """
        p = params.p
        # Cache vector spaces to improve sage performace
        to_patch = [GF(2), GF(3), GF(3**2)]
        for x in to_patch:
            type(x).vector_space = cached_method(type(x).vector_space)
        self.id_torsion = params.id_torsion
        if not precomps:
            precomps = Precomputations(p, 122, role = role)
        self.precomps = precomps

    def keygen(self):
        """
        Generate public and secret keys
        """
        # Generation of a random ideal Isk of norm N
        Isk = self.precomps.random_ideal()
        self._precompute_secret_key(Isk)

    def identify(self, q):
        """
        Response to a prime challenge q
        """
        Isk, Mpk = self.sk
        precomps = self.precomps

        # Generate and pushforward the ideal of given norm
        assert q < 2**self.id_torsion

        Ir_norm = q * (2**self.id_torsion - q)
        Ir_prime = precomps.random_ideal_of_given_norm(Ir_norm, prime=False)

        Ir = pushforward_ideal(Ir_prime, Isk)
        Iaux = Isk * Ir

        # Compute the isogeny
        phi_aux = IdealToIsogenyClapotis(Iaux, precomps = precomps)
        Eaux, Paux, Qaux = phi_aux.images()
        Pr, Qr = self.precomps.evaluate_matrix(Mpk, (Paux, Qaux))

        # We only need smaller torsion
        torsion_delta = 2**(precomps.a - self.id_torsion)
        Pr *= torsion_delta
        Qr *= torsion_delta

        # Point compression
        P0, Q0 = torsion_basis_2e(Eaux, self.id_torsion, montgomery=True)
        (a1, a2), (a3, a4) = self.precomps.BiDLPs((Pr, Qr), P0, Q0, self.id_torsion)
        return Eaux, (a1, a2, a3, a4)

    def verify(self, q, res, pk):
        """
        Given the kernel of an higher dimensional isogeny verify
        that it splits and that the degree is correct.
        """
        precomps = self.precomps

        Er, pts = res
        a1, a2, a3, a4 = pts

        Ppk, Qpk = pk[1]

        P = ZZ(q) * Ppk
        Q = ZZ(q) * Qpk

        tic = time.time()
        P0, Q0 = torsion_basis_2e(Er, self.id_torsion, montgomery=True)
        Pr = a1*P0 + a2*Q0
        Qr = a3*P0 + a4*Q0
        print(f'Point decompression (inside verification) : {time.time() - tic:.3f}s')

        # Check response representation
        ker_Phi = ((P, Pr), (Q, Qr))

        Phi = Dim2Isogeny(ker_Phi, self.id_torsion, precomps = precomps)

        # Check degrees using pairings à la SQIsign2D-East
        idr = Er(0)
        P1, P2 = Phi((Ppk, idr))
        Q1, Q2 = Phi((Qpk, idr))

        D = 2**self.id_torsion

        # Pairing check
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
        # c = precomps.c
        e = precomps.e

        # Computing the isogeny associated with the ideal Isk
        phi_sk = IdealToIsogenyClapotis(Isk, precomps = self.precomps)
        Epk, phi_sk_P0, phi_sk_Q0 = phi_sk.images()

        self._precompute_public_key(Epk)
        Ppk, Qpk = self.pk_basis

        Mpk = matrix(
            Zmod(2**e), precomps.BiDLPs((Ppk, Qpk), phi_sk_P0, phi_sk_Q0, e)
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

        torsion_delta = 2**(precomps.a - self.id_torsion)
        Ppk *= torsion_delta
        Qpk *= torsion_delta
        self.pk_id_basis = (Ppk, Qpk)
        self.eWPpkQpk = weil_pairing_pari(Ppk, Qpk, 2**self.id_torsion)

