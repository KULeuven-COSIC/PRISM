# import pari for fast dlog
import cypari2
from utilities.biextension import fast_twotothem_tate, fast_twotothem_weil

# Make instance of Pari
pari = cypari2.Pari()

# ===================================== #
#   Directly access pairings from Pari  #
# ===================================== #


def weil_pairing_pari(P, Q, D, check=False):
    """
    Wrapper around Pari's implementation of the Weil pairing
    Allows the check of whether P,Q are in E[D] to be optional
    """
    if check:
        nP, nQ = D * P, D * Q
        if nP.is_zero() or nQ.is_zero():
            raise ValueError("points must both be n-torsion")

    return pari.ellweilpairing(P.curve(), P, Q, D)


def tate_pairing_pari(P, Q, D):
    """
    Wrapper around Pari's implementation of the Tate pairing
    NOTE: this is not the reduced Tate pairing, so to make this
    match with SageMath you need

    P.tate_pairing(Q, D, k) == pari.elltatepairing(E, P, Q, D)**((p^k - 1) / D)
    """
    E = P.curve()
    return pari.elltatepairing(E, P, Q, D)

# warning: we assume here that we work on a Montgomery curve and that D=2^e
def tate_pairing_biextension(P, Q, e, exp=1):
    from utilities.supersingular import montgomery_coefficient
    E = P.curve()
    A = montgomery_coefficient(E)
    A24 = (A+2)/4
    PQ=P+Q
    #print(f"P={P}, Q={Q}, P+Q={PQ}")
    # special cases
    F=P.base_ring()
    if PQ.is_zero() or P.is_zero() or Q.is_zero():
        return F(1) #we use the fact that we are on E/F_p^2, E supersingular
    elif Q[0]==0 or PQ[0]==0: #Q=T, or P+Q=T, so e(P,Q)=e(P,T)
        return (P[1]**(2**(e-1)))**exp
    elif P[0]==0:
        return (Q[1]**(2**(e-1)))**(-exp)
    return fast_twotothem_tate(e, A24, P[0], Q[0], PQ[0], exp=exp)

def weil_pairing_biextension(P, Q, e):
    from utilities.supersingular import montgomery_coefficient
    E = P.curve()
    A = montgomery_coefficient(E)
    A24 = (A+2)/4
    PQ=P+Q
    F=P.base_ring()
    # special cases
    if PQ.is_zero() or P.is_zero() or Q.is_zero():
        return F(1)
    elif Q[0]==0:
        return -weil_pairing_biextension(Q, P, e)
    elif PQ[0]==0: #Q=T, or P+Q=T, so e(P,Q)=e(P,T)
        return -weil_pairing_biextension(PQ, P, e)
    elif P[0]==0:
        QQ = 2**(e-1)*Q
        if QQ == P or QQ.is_zero():
            return F(1)
        else:
            return F(-1)
    return fast_twotothem_weil(e, A24, P[0], Q[0], PQ[0])

def tate_pairing_default(P, Q, D, exp=1):
    return tate_pairing_pari(P, Q, D)**exp

def tate_pairing_power_two(P, Q, e, exp=1):
    return tate_pairing_pari(P, Q, 2**e)**exp
    #return tate_pairing_biextension(P, Q, e, exp=exp)

def weil_pairing_default(P, Q, D):
    return weil_pairing_pari(P, Q, D)

def weil_pairing_power_two(P, Q, e):
    return weil_pairing_pari(P, Q, 2**e)
    #return weil_pairing_biextension(P, Q, e)
