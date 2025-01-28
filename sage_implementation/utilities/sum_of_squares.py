from sage.all import gcd, is_pseudoprime, two_squares, Factorization, ZZ


def sum_of_squares_friendly(n):
    """
    We can write any n = x^2 + y^2 providing that there
    are no prime power factors p^k | n such that
    p ≡ 3 mod 4 and k odd.

    We can reject bad n by checking these conditions.
    """


    # We consider the odd part of n and try and determine if there are bad factors
    n_val = n.valuation(2)
    n_odd = n // (2**n_val)
    fact = Factorization([(2, n_val)])

    # If n_odd % 4 == 3, then there is certainly a prime factor p^k | n_odd such that p = 3 mod 4 and k is odd, so this is no good
    if n_odd % 4 == 3:
        return False, fact

    # We don't want any prime factors in n such that p ≡ 3 mod 4
    # Technically, we could allow these if they appear with even
    # exp. but this seems good enough.
    # We take the produce p ≡ 3 mod 4 for p < 500
    bad_prime_prod = 9758835872005600424541432018720696484164386911251445952809873289484557900064682994157437552186699082331048989
    if gcd(bad_prime_prod, n) != 1:
        return False, fact

    cofactor = ZZ(1)
    # Now we remove all good prime factors p ≡ 1 mod 4
    # We take the produce p ≡ 1 mod 4 for p < 500
    good_prime_prod = 6396589037802337253083696670901044588941174775898621074430463772770231888063102454871836215129505
    g = gcd(good_prime_prod, n_odd)
    while g != 1:
        #fact = Factorization([*fact, *g.factor()])
        cofactor *= g
        n_odd = n_odd // g
        g = gcd(n_odd, g)

    # We now have two cases, either n_odd is a prime p ≡ 1 mod 4
    # and we have a solution, or n_odd is a composite which is
    # good, or has an even number of primes p ≡ 3 mod 4. The second
    # case is bad and we will only know by fully factoring it, so
    # to avoid the expensive factoring, we will reject these cases
    if n_odd == 1:
        return True, Factorization([*fact, *cofactor.factor()])

    else:
        return is_pseudoprime(n_odd), Factorization([*fact, *cofactor.factor(), (n_odd, 1)])

def two_squares_factored(factors):
    """
    This is the function `two_squares` from sage, except we give it the
    factorisation of n already.
    """
    # these are the primes in good_prime_prod
    small_squares = {5: (1, 2), 13: (2, 3), 17: (1, 4), 29: (2, 5), 37: (1, 6), 41: (4, 5), 53: (2, 7), 61: (5, 6), 73: (3, 8), 89: (5, 8), 97: (4, 9), 101: (1, 10), 109: (3, 10), 113: (7, 8), 137: (4, 11), 149: (7, 10), 157: (6, 11), 173: (2, 13), 181: (9, 10), 193: (7, 12), 197: (1, 14), 229: (2, 15), 233: (8, 13), 241: (4, 15), 257: (1, 16), 269: (10, 13), 277: (9, 14), 281: (5, 16), 293: (2, 17), 313: (12, 13), 317: (11, 14), 337: (9, 16), 349: (5, 18), 353: (8, 17), 373: (7, 18), 389: (10, 17), 397: (6, 19), 401: (1, 20), 409: (3, 20), 421: (14, 15), 433: (12, 17), 449: (7, 20), 457: (4, 21), 461: (10, 19)}

    F = factors
    for (p,e) in F:
        if e % 2 == 1 and p % 4 == 3:
            raise ValueError("%s is not a sum of 2 squares"%n)

    n = factors.expand()
    from sage.rings.finite_rings.integer_mod import Mod
    a = ZZ.one()
    b = ZZ.zero()
    for (p,e) in F:
        if e >= 2:
            m = p ** (e//2)
            a *= m
            b *= m
        if e % 2 == 1:
            if p == 2:
                # (a + bi) *= (1 + I)
                a,b = a - b, a + b
            else:  # p = 1 mod 4
                if p in small_squares:
                    r,s = small_squares[p]
                else:
                    # Find a square root of -1 mod p.
                    # If y is a non-square, then y^((p-1)/4) is a square root of -1.
                    y = Mod(2,p)
                    while True:
                        s = y**((p-1)/4)
                        if not s*s + 1:
                            s = s.lift()
                            break
                        y += 1
                    # Apply Cornacchia's algorithm to write p as r^2 + s^2.
                    r = p
                    while s*s > p:
                        r,s = s, r % s
                    r %= s

                # Multiply (a + bI) by (r + sI)
                a,b = a*r - b*s, b*r + a*s

    a = a.abs()
    b = b.abs()
    assert a*a + b*b == n
    if a <= b:
        return (a,b)
    else:
        return (b,a)

def sum_of_squares(M):
    """
    Given integers M attempts to compute x,y such that
    M = x^2 + y^2
    """

    M = ZZ(M)
    if M.nbits() <= 32:
        from sage.rings import sum_of_squares
        return sum_of_squares.two_squares_pyx(M)

    b, fact = sum_of_squares_friendly(M)
    if not b:
        #print(f"sum_of_squares: found no solutions for {M}")
        return []

    try:
        sol = two_squares_factored(fact)
        #print(f"sum_of_squares: found solution {M}={sol[0]}^2+{sol[1]}^2")
    except:
        sol = []
    return sol
