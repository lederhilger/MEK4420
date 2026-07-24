from numpy import sqrt
from functools import cached_property

class AGM:
    def __init__(self, a_0, b_0, c_0, ε):
        self.a_0 = a_0
        self.b_0 = b_0
        self.c_0 = c_0
        self.ε = ε

    @cached_property
    def alphabet(self) -> tuple:
        a = [self.a_0]; b = [self.b_0]; c = [self.c_0]
        while c[-1] > self.ε:
            a_N = .5*(a[-1] + b[-1])
            b_N = sqrt(a[-1]*b[-1])
            c_N = .5*(a[-1] - b[-1])
            a.append(a_N); b.append(b_N); c.append(c_N)
        return a, b, c

    @property
    def mean(self) -> float:
        a, b, c = self.alphabet
        return a[-1]

    def archimedes(self) -> tuple:
        a = [self.a_0]; b = [self.b_0]; c = [self.c_0]
        while c[-1] > self.ε:
            a_N = sqrt(a[-1]*b[-1])
            b_N = 2*a_N*b[-1]/(a_N + b[-1])
            c_N = .5*(a[-1] - b[-1])
            a.append(a_N); b.append(b_N); c.append(c_N)
        return a, b

    
def π() -> float:
    a_0 = 1; b_0 = .5*sqrt(2); c_0 = a_0 - b_0
    agm = AGM(a_0, b_0, c_0, 1e-16)
    a, b, c = agm.alphabet
    denominator = 0
    for n in range(len(a)):
        denominator += 2**n * (a[n]**2 - b[n]**2)
    denominator = 1 - denominator
    π = 2*a[-1]**2 / denominator
    return π
