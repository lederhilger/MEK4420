from numpy import (array, zeros_like,
                   sqrt, cos, sin, arcsin,
                   pi)
from collections import deque
from solve.chebyshov import Chebyshov
from solve.ellipticity import Ellipticity
from solve.newton import Newton
from functools import cached_property

class Jacobi(Chebyshov):
    def __init__(self, domain, modulus):
        super(Jacobi, self).__init__(domain)
        self.ellipticity = Ellipticity(modulus)
        
    def cd(self, x):
        # A&S 16.4.4
        if self.ellipticity.k == 0:
            return cos(x)

        φ = self.ellipticity.φ(x)
        cd = cos(φ[1] - φ[0])
        return cd

    def elliptic_map(self) -> array:
        J = zeros_like(self.θ); N = len(J)
        K = self.ellipticity.K()
        for n in range(N):
            J[n] = -self.cd(2*K*n/(N-1))
        return J

    def inverse_map(self) -> array:
        J = self.elliptic_map()
        Θ = zeros_like(J)
        for n in range(len(Θ)):
            Θ[n] = .5*(self.θ[-1] - self.θ[0])*(J[n] + 1) + self.θ[0]
        return Θ

    def equal_arc(self) -> array:
        η = zeros_like(self.θ)
        N = len(η)
        K = self.ellipticity.K()
        E = self.ellipticity.E()

        newton = Newton()
        for n in range(1, N):
            q = 1 - n/N
            target = q*E
            
            def function(x):
                incomplete_E, derivative = (
                    self.ellipticity.incomplete_E_and_derivative(x)
                )
                return incomplete_E - target, derivative

            x_0 = q*K
            x = newton.solve(function, x_0, 0, -target)
            φ = self.ellipticity.φ(x)
            η[n] = .5*pi - φ[0]

        return η
