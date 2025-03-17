from numpy import array, zeros_like, cos, pi
class Chebyshov:
    def __init__(self, domain: array):
        self.θ = domain
    
    def unit_map(self) -> array:
        Θ = zeros_like(self.θ)
        for n in range(len(Θ)):
            Θ[n] = 2*(self.θ[n] - self.θ[0])/(self.θ[-1] - self.θ[0]) - 1
        return Θ

    def chebyshov_map(self) -> array:
        Ч = zeros_like(self.θ); N = len(Ч)
        for n in range(N):
            Ч[n] = -cos(2*n*pi/(2*(N-1)))
        return Ч
    
    def inverse_map(self) -> array:
        Ч = self.chebyshov_map()
        Θ = zeros_like(Ч)
        for n in range(len(Θ)):
            Θ[n] = .5*(self.θ[-1] - self.θ[0])*(Ч[n] + 1) + self.θ[0]
        return Θ