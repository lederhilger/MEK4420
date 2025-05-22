from numpy import (array, zeros_like,
                   sqrt, cos, sin, arcsin,
                   pi)
from collections import deque
from solve.chebyshov import Chebyshov

class Jacobi(Chebyshov):
    def __init__(self, domain):
        super(Jacobi, self).__init__(domain)
        self.sqrtm = sqrt(.5)
        self.α = arcsin(self.sqrtm)

    def AGM(self, a_0: float, b_0: float, c_0: float, ε = 1e-16) -> tuple:
        a = [a_0]; b = [b_0]; c = [c_0]
        while c[-1] > ε:
            a_N = .5*(a[-1]+b[-1])
            b_N = sqrt(a[-1]*b[-1])
            c_N = .5*(a[-1]-b[-1])
            a.append(a_N); b.append(b_N); c.append(c_N)
        return a,b,c
    
    def archimedes(self, a_0: float, b_0: float, c_0: float, ε = 1e-16) -> tuple:
        a = [a_0]; b = [b_0]; c = [c_0]
        while c[-1] > ε:
            a_N = sqrt(a[-1]*b[-1])
            b_N = 2*a_N*b[-1]/(a_N + b[-1])
            c_N = .5*(a[-1] - b[-1])
            a.append(a_N); b.append(b_N); c.append(c_N)
        return a, b

    def π(self) -> float:
        a_0 = 1; b_0 = .5*sqrt(2); c_0 = a_0 - b_0
        a, b, c = self.AGM(a_0, b_0, c_0)
        denominator = 0
        for n in range(len(a)):
            denominator += 2**n * (a[n]**2 - b[n]**2)
        denominator = 1 - denominator
        π = 2*a[-1]**2 / denominator
        print(f"a: {a}")
        print(f"denominator: {denominator}")
        return π

    def φ(self, a, c, x: array):
        N = len(a) - 1
        φ = deque([2**N * a[N] * x])
        for n in range(N):
            φ.appendleft(.5*(arcsin((c[N-n]/a[N-n])*sin(φ[0])) + φ[0]))
        return φ

    def K(self) -> float:
        # A&S 17.6.3
        a, b, c = self.AGM(1, cos(self.α), sin(self.α))
        K = .5*pi/a[-1]
        return K

    def E(self) -> float:
        a, b, c = self.AGM(1, cos(self.α), sin(self.α))
        K = .5*pi/a[-1]
        e = 0
        for k in range(len(c)):
            e += 2**k * c[k]**2
        E = K*(1-.5*e)
        return E

    def cd(self, x):
        a, b, c = self.AGM(1, sqrt(.5), sqrt(.5))
        φ = self.φ(a, c, x)
        cd = cos(φ[1] - φ[0])
        return cd

    def elliptic_map(self) -> array:
        J = zeros_like(self.θ); N = len(J)
        K = self.K()
        for n in range(N):
            J[n] = -self.cd(2*K*n/(N-1))
        return J

    def inverse_map(self) -> array:
        J = self.elliptic_map()
        Θ = zeros_like(J)
        for n in range(len(Θ)):
            Θ[n] = .5*(self.θ[-1] - self.θ[0])*(J[n] + 1) + self.θ[0]
        return Θ
