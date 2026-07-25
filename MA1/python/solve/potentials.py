from numpy import zeros, arctan2, cos, sin, sqrt, pi
# Plots (0, 2pi) for polar coordinates
class Potentials:
    def __init__(self, a: float, b: float, N: int, domain: tuple):
        self.a = a
        self.b = b
        self.N = N
        self.domain = domain

    def circle_1(self):
        phi = cos(self.domain)
        return phi

    def circle_2(self):
        phi = sin(self.domain)
        return phi
    
    def ellipse_1(self):
        denominator = sqrt(self.b**2 * cos(self.domain)**2 + self.a**2 * sin(self.domain)**2)
        phi = self.b**2 * cos(self.domain) / denominator
        return phi

    def ellipse_2(self):
        denominator = sqrt(self.b**2 * cos(self.domain)**2 + self.a**2 * sin(self.domain)**2)
        phi = self.a**2 * sin(self.domain) / denominator
        return phi

    def ellipse_6(self):
        denominator = self.b**2 * cos(self.domain)**2 + self.a**2 * sin(self.domain)**2
        phi = self.a*self.b*(self.a**2 - self.b**2)*sin(2*self.domain) / (4*denominator)
        return phi
