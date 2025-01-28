# 2S+2M+1m0+4a
# a=(X+Z)^2, b=(X-Z)^2, c=4XZ, X2=ab, Z2=c(b+A24c)
def double(A24, P):
    X, Z = P
    XpZ=X+Z
    XmZ=X-Z
    a=(XpZ)*(XpZ) #in Sage this is faster than XpZ**2
    b=(XmZ)*(XmZ)
    c=a-b
    X2=a*b
    Z2=c*(b+A24*c)
    return (X2, Z2)

# 2S+3M+6a
# X(P+Q)X(P-Q)=[(XP+ZP)(XQ-ZQ)+(XP-ZP)(XQ+ZQ)]^2
# Z(P+Q)Z(P-Q)=[(XP+ZP)(XQ-ZQ)-(XP-ZP)(XQ+ZQ)]^2
def diff_add(P, Q, ixPmQ):
    XP, ZP = P
    XQ, ZQ = Q
    a=XP+ZP
    b=XP-ZP
    c=XQ+ZQ
    d=XQ-ZQ
    da = d*a
    cb = c*b
    dapcb = da+cb
    damcb = da-cb
    XPQ=dapcb*dapcb * ixPmQ
    ZPQ=damcb*damcb
    return (XPQ, ZPQ)

def ratio(P1, P2, P1bis, P2bis):
    xP1, zP1 = P1
    xP2, zP2 = P2
    xP1b, zP1b = P1bis
    xP2b, zP2b = P2bis
    ## sanity checks
    assert zP1b*xP1 == zP1 * xP1b
    assert zP2b*xP2 == zP2 * xP2b
    #return (xP1 * xP2b / (xP1b * xP2))
    return  (xP1b * xP2)/(xP1 * xP2b)

def fast_twotothem_ladder(m, A24, P, Q, PQ, ixP):
    nQ=Q
    PnQ=PQ
    if m==0:
        return (nQ, PnQ)
    # Doublings in the biextension
    for i in range(0, m):
        PnQ=diff_add(PnQ, nQ, ixP)
        nQ=double(A24, nQ)
    return (nQ, PnQ)

def fast_twotothem_biext_ladder(m, A24, xP, xQ, xPQ):
    ixP = 1/xP
    P = (xP, 1)
    Q = (xQ, 1)
    PQ = (xPQ, 1)
    nQ, PnQ = fast_twotothem_ladder(m, A24, P, Q, PQ, ixP)
    return (P, nQ, PnQ)

def translate(P, T):
    A, B=T
    x, z=P
    if B == 0: # T = (1:0) is the neutral point
        return P
    elif A == 0: # T = (0:1)
        return (z,x)
    else:
        return (A*x-B*z, B*x-A*z)

def fast_twotothem_non_reduced_tate(m, A24, xP, xQ, xPQ):
    neutral = (1, 0)
    P, nQ, PnQ = fast_twotothem_biext_ladder(m-1, A24, xP, xQ, xPQ)
    PnQ = translate(PnQ, nQ)
    nQ = translate(nQ, nQ)
    return ratio(neutral, P, nQ, PnQ)

def fast_twotothem_tate(m, A24, xP, xQ, xPQ, exp=1):
    #p=xP.base_ring().characteristic()
    e = fast_twotothem_non_reduced_tate(m, A24, xP, xQ, xPQ)
    return e**exp

def fast_twotothem_weil(m, A24, xP, xQ, xPQ):
    e1 = fast_twotothem_non_reduced_tate(m, A24, xP, xQ, xPQ)
    e2 = fast_twotothem_non_reduced_tate(m, A24, xQ, xP, xPQ)
    return e1/e2
