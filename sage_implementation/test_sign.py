import time

from sage.all import proof

import params
from precomputations import Precomputations
from PRISM_sign import PRISM_sign

def run_protocol(signer, verifier):
    init_time = time.time()
    signer.keygen()
    keygen_time = time.time()
    print(f'Key generation took {(keygen_time - init_time):.3f}s')

    import random
    msg = b'Important message ' + str(random.randint(1, 100)).encode()
    sig = signer.sign(msg)
    sign_time = time.time()
    print(f'Signing took {(sign_time - keygen_time):.3f}s, counter = {sig[-1]}')

    pk = (signer.pk_curve, signer.pk_sign_basis, signer.eWPpkQpk)
    assert verifier.verify(msg, sig, pk)

    verif_time = time.time()
    print(f'Verification took {(verif_time - sign_time):.3f}s')
    return keygen_time - init_time, sign_time - keygen_time, verif_time - sign_time

if __name__ == '__main__':
    proof.all(False)

    begin_time = time.time()
    signer = PRISM_sign('signer')
    verifier = PRISM_sign('verifier')
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

