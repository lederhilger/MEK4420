from numpy import zeros, arctan2, cos, sin, sqrt, pi
# Plots (0, 2pi) for polar coordinates
class Potentials:
    def __init__(self, a: float, b: float, N: int, domain: tuple):
        self.a = a
        self.b = b
        self.N = N
        self.domain = [domain[n]-pi for n in range(len(domain))]

    def circle_1(self):
        phi = cos(self.domain)
        return phi

    def circle_2(self):
        phi = sin(self.domain)
        return phi
    
    def ellipse_1(self):
        phi = zeros(self.N)
        for n in range(self.N):
            denominator = sqrt(self.b**2 * cos(self.domain[n])**2 + self.a**2 * sin(self.domain[n])**2)
            phi[n] = self.b**2 * cos(self.domain[n]) / denominator
        return phi

    def ellipse_2(self):
        phi = zeros(self.N)
        for n in range(self.N):
            denominator = sqrt(self.b**2 * cos(self.domain[n])**2 + self.a**2 * sin(self.domain[n])**2)
            phi[n] = self.a**2 * sin(self.domain[n]) / denominator
        return phi
