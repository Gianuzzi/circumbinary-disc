# Possible useful Functions
from numpy import exp

#1. : Trapezoid
def trapezoidal(f, a, b, n=1):
    
    if b < a : raise ValueError("a must be lower than b (received a=%f, b=%f)" %(a,b))
    if not isinstance(n,int): raise TypeError("n must be an integer. (received n={})".format(n))
    if n < 1 : raise ValueError("n must be greater than 1 (received n=%d)" % n)

    if n == 2: return (b - a) * 0.5 * (f(b) + f(a))    
    
    dx  = (0.5 * (b - a)) / float(n)
    x   = a
    val = 0
    while x < b:
        val += f(x) * dx
        x   += dx
    return val


#2. : Simpson
def simpson(f, a, b, n=2):
    
    if b < a: raise ValueError("a must be lower than b (received a=%f, b=%f)" %(a,b))
    if not isinstance(n,int): raise TypeError("n must be an integer. (received n={})".format(n))
    if n < 1: raise ValueError("n must be greater than 1 (received n=%d)" % n)
    if n % 2: raise ValueError("n must be even (received n=%d)" % n)

    dx  = (b - a) / float(n)
    val = f(a) + f(b)

    for i in range(1, n,   2): val += 4 * f(a + i * dx)
    for i in range(2, n-1, 2): val += 2 * f(a + i * dx)

    return val * dx / float(3)

# 3. : Boole:

def boole(f, a, b, n=5):
    
    if b < a: raise ValueError("a must be lower than b (received a=%f, b=%f)" %(a,b))
    if not isinstance(n,int): raise TypeError("n must be an integer. (received n={})".format(n))
    if n < 1: raise ValueError("n must be greater than 1 (received n=%d)" % n)
    if n % 5: raise ValueError("n must be multiple of 5 (received n=%d)" % n)
        
    h   = (b - a) / float(n)
    val = 0
    x   = a
    while  x < b:
        bl = (7  * f(x) + 
              32 * f(x + h) +
              12 * f(x + 2 * h) + 
              32 * f(x + 3 * h) + 
              7  * f(x + 4 * h)) * 2 * h / float(45)
        val += bl
        x   += h * 4 
    return val

# Publication's profile, without normalization
def thun_prof(r, alpha=1.5, Rgap=0.55, DR=0.055):
    return r**(-alpha) / (1 + exp((Rgap - r) / DR))

def surf_prof(r, alpha=1.5, Rgap=0.55, DR=0.055):
    return r * thun_prof(r, alpha, Rgap, DR)
    
