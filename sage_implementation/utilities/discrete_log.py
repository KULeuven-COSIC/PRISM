# Sage imports
from sage.all import ZZ

# import pari for fast dlog
import cypari2


# ===================================== #
#  Fast DLP solving using pairings      #
# ===================================== #

# Make instance of Pari
pari = cypari2.Pari()


def discrete_log_pari(a, base, order):
    """
    Wrapper around pari discrete log. Works like a.log(b),
    but allows us to use the optional argument order. This
    is important as we skip Pari attempting to factor the
    full order of Fp^k, which is slow.
    """
    x = pari.fflog(a, base, order)
    return ZZ(x)


def _precompute_baby_steps(base, step, e):
    """
    Helper function to compute the baby steps for
    pohlig_hellman_base and windowed_pohlig_hellman.
    """
    baby_steps = [base]
    for _ in range(e):
        base = base**step
        baby_steps.append(base)
    return baby_steps


def pohlig_hellman_base(a, base, e):
    """
    Solve the discrete log for a = base^x for
    elements base,a of order 2^e using the
    Pohlig-Hellman algorithm.
    """
    baby_steps = _precompute_baby_steps(base, 2, e)

    dlog = 0
    exp = 2 ** (e - 1)

    # Solve the discrete log mod 2, only 2 choices
    # for each digit!
    for i in range(e):
        if a**exp != 1:
            a /= baby_steps[i]
            dlog += 2**i

        if a == 1:
            break

        exp //= 2

    return dlog


def windowed_pohlig_hellman(a, base, e, window):
    """
    Solve the discrete log for a = base^x for
    elements base,a of order 2^e using the
    windowed Pohlig-Hellman algorithm following
    https://ia.cr/2016/963.

    Algorithm runs recursively, computing in windows
    l^wi for window=[w1, w2, w3, ...].

    Runs the base case when window = []
    """
    # Base case when all windows have been used
    if not window:
        return pohlig_hellman_base(a, base, e)

    # Collect the next window
    w, window = window[0], window[1:]
    step = 2**w

    # When the window is not a divisor of e, we compute
    # e mod w and solve for both the lower e - e mod w
    # bits and then the e mod w bits at the end.
    e_div_w, e_rem = divmod(e, w)
    e_prime = e - e_rem

    # First force elements to have order e - e_rem
    a_prime = a ** (2**e_rem)
    base_prime = base ** (2**e_rem)

    # Compute base^(2^w*i) for i in (0, ..., e/w-1)
    baby_steps = _precompute_baby_steps(base_prime, step, e_div_w)

    # Initialise some pieces
    dlog = 0
    if e_prime:
        exp = 2 ** (e_prime - w)

        # Work in subgroup of size 2^w
        s = base_prime ** (exp)

        # Windowed Pohlig-Hellman to recover dlog as
        # alpha = sum l^(i*w) * alpha_i
        for i in range(e_div_w):
            # Solve the dlog in 2^w
            ri = a_prime**exp
            alpha_i = windowed_pohlig_hellman(ri, s, w, window)

            # Update a value and dlog computation
            a_prime /= baby_steps[i] ** (alpha_i)
            dlog += alpha_i * step**i

            if a_prime == 1:
                break

            exp //= step

    # TODO:
    # I don't know if this is a nice way to do
    # this last step... Works well enough but could
    # be improved I imagine...
    exp = 2**e_prime
    if e_rem:
        base_last = base**exp
        a_last = a / base**dlog
        dlog_last = pohlig_hellman_base(a_last, base_last, e_rem)
        dlog += exp * dlog_last

    return dlog

def discrete_log_pair_power_two(pair_a, pair_b, pair_PQ, e, window=None):
    if window is None:
        D = 2**e
        a = discrete_log_pari(pair_a, pair_PQ, D)
        b = discrete_log_pari(pair_b, pair_PQ, D)
    else:
        a = windowed_pohlig_hellman(pair_a, pair_PQ, e, window)
        b = windowed_pohlig_hellman(pair_b, pair_PQ, e, window)
    return a,b

