from numpy import array, pi, arcsin, sin, sqrt
from collections import deque
from functools import cached_property
from solve.agm import AGM

class Ellipticity:
    def __init__(self, modulus, ε=1e-16):
        self.k = modulus
        if not 0 <= self.k < 1:
            raise ValueError("Degenerate.")
        self.kprime = sqrt((1+self.k)*(1-self.k))
        self.ε = ε

    @cached_property
    def elliptic_AGM(self) -> tuple:
        agm = AGM(1, self.kprime, self.k, self.ε)
        return agm.alphabet

    def φ(self, x: array):
        # A&S 16.4.3
        a, b, c = self.elliptic_AGM
        N = len(a) - 1
        φ = deque([2**N * a[N] * x])
        for n in range(N):
            φ.appendleft(.5*(arcsin((c[N-n]/a[N-n])*sin(φ[0])) + φ[0]))
        return φ

    def K(self) -> float:
        # A&S 17.6.3
        a, b, c = self.elliptic_AGM
        K = .5*pi/a[-1]
        return K

    def E(self) -> float:
        # A&S 17.6.4
        a, b, c = self.elliptic_AGM
        K = .5*pi/a[-1]
        e = 0
        for n in range(len(c)):
            e += 2**n * c[n]**2
        E = K*(1-.5*e)
        return E

    def incomplete_E(self, x):
        a, b, c = self.elliptic_AGM
        φ = self.φ(x)
        E = x*self.E()/self.K()
        for n in range(1, len(c)):
            E += c[n]*sin(φ[n])
        return E

    def incomplete_E_and_derivative(self, x):
        a, b, c = self.elliptic_AGM
        φ = self.φ(x)
        E = x*self.E()/self.K()
        for n in range(1, len(c)):
            E += c[n]*sin(φ[n])
        derivative = 1 - self.k**2 * sin(φ[0])**2
        return E, derivative
