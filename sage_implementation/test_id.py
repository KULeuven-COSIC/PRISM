import time

from sage.all import random_prime, proof

import params
from PRISM_id import PRISM_id

def run_protocol(signer, verifier):
    init_time = time.time()
    signer.keygen()
    keygen_time = time.time()
    print(f'Key generation took {(keygen_time - init_time):.3f}s')

    id_torsion = params.id_torsion
    q = random_prime(2**id_torsion, lbound=2**(id_torsion-1))

    res = signer.identify(q)

    id_time = time.time()
    print(f'Identification took {(id_time - keygen_time):.3f}s')

    pk = (signer.pk_curve, signer.pk_id_basis, signer.eWPpkQpk)
    assert verifier.verify(q, res, pk)

    verif_time = time.time()
    print(f'Verification took {(verif_time - id_time):.3f}s')
    return keygen_time - init_time, id_time - keygen_time, verif_time - id_time

if __name__ == '__main__':
    proof.all(False) # Speed up sage
    begin_time = time.time()
    signer = PRISM_id('signer')
    verifier = PRISM_id('verifier')
    init_time = time.time()
    print(f'Precomputations took {(init_time - begin_time):.3f}s')
    print('===')

    key_gen = []
    ident = []
    verif =  []
    n_runs = 10

    for _ in range(n_runs):
        a, b, c = run_protocol(signer, verifier)
        key_gen.append(a)
        ident.append(b)
        verif.append(c)
        print('===')

    print(f'Average keygen time = {sum(key_gen)/len(key_gen):.3f}s')
    print(f'Average id time = {sum(ident)/len(ident):.3f}s')
    print(f'Average verification time = {sum(verif)/len(verif):.3f}s')

