# Pollard's Rho Algorithm
def rho_seq(x, n):
    return (x^2 + 1) % n
def pollard_rho(n):
    x = y = 2 #x0 is 2
    g = 1
    while g == 1: #using Floyd's cycle finding algorithm
        x = rho_seq(x, n) #f(x) mod n
        y = rho_seq(rho_seq(y, n), n) #f(f(y)) mod n
        g = gcd(abs(x - y), n) #gcd(|x – y|, n)
    return g


# Brent's Rho Algorithm
def brent_rho(n):
    y = 2 #x0 is 2
    r = q = G = 1 #ρ is 2 and u is 0, powers of 2
    m = (n^(1/4) – log(n)) / 2 #estimated length of the non-periodic segment of the sequence
    while G == 1: #using Sedgewick-Szymanski cycle finding algorithm
        x = y
        for i in xrange(r):
            y = rho_seq(y, n) #f(y) mod n
        k = 0
        while G == 1 and k < r:
            ys = y
            for i in xrange(min(m, r-k)):
                y = rho_seq(y, n) #f(y) mod n
                q = (q * abs(x - y)) % n #q *= |x – y| mod n
            G = gcd(q, n)
            k += m
    if G == n: #if we chose an m that is too big, we need to backtrack
        G = 1
        while G == 1:
            ys = rho_seq(ys, n) #f(ys) mod n
            G = gcd(abs(x - ys), n) #gcd(|x – ys|, n)
    if G == n: #if we haven’t found a cycle n must be prime
        return -1
    else:
        return G


# Pollard's p-1 Algorithm
def pollard_p1(n):
    smoothness_bound = round(n^(1/6)) #a reasonable compromise between speed and likelihood of a result
    ring = Integers(n) #ring mod n
    while True:
        k = prod(p^(floor(log(n, p))) for p in prime_range(1, smoothness_bound)) #the product of all primes p < smoothness_bound taken to the largest integer exponent that will yield a result less than n
        a = ring(2) #2 mod n
        q = int(gcd(a^k - 1, n))
        if q > 1 and q < n: #factor found
            return q
        elif q == 1: #no prime factors p for which p-1 is powersmooth for the given smoothness bound
            smoothness_bound += 2
        elif q == n: #all factors were powersmooth for the given smoothness bound
            smoothness_bound -= 2
        else:
            return -1


# Lenstra's Elliptic Curve Algorithm
def lenstra(n):
    q = 1
    while q == 1:
        x, y, a = randint(0, n - 1), randint(0, n - 1), randint(0, n - 1) #random elliptic curve with equation y^2 = x^3 + a*x + b
        b = mod(y^2 - x^3 - x*a, n) #solving for the remaining coefficient
        z = 4*a^3 + 27*b^2 #solve for z
        q = gcd(z, n)
    if q > 1:
        return q
    else:
        return -1
