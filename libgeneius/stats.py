import math
import sys

def getSNPpvals( cvg, freq):
    confidence = prob(round(cvg*freq), .5, cvg)
    if confidence > 0.5:
        confidence = 1-confidence
    conf_het = confidence * 200
    p_homo = confidence
    return (conf_het, p_homo)

def __main__():
    cvg = sys.argv[1]
    freq = sys.argvp[2]
    (v1,v2)=getSNPpvals(cvg,freq)
    print cvg, freq, v1, v2
    
def factorial(n): 
    if n < 2: return 1
    return reduce(lambda x, y: x*y, xrange(2, int(n)+1))

def prob(s, p, n):
    x = 1.0 - p
    a = n - s
    b = s + 1.0
    c = a + b - 1.0
    prob = 0
    for j in xrange(a, c + 1):
        prob += factorial(c) / (factorial(j)*factorial(c-j)) \
                * x**j * (1 - x)**(c-j)
    return prob
