from numpy import array, zeros_like, zeros, sqrt
from solve.box import Box

class IntegralEquation(Box):
    def __init__(self, Nx: int, Ny: int, L: float, D: float):
        super(IntegralEquation, self).__init__(Nx, Ny, L, D)
        self.xp = super(IntegralEquation, self).x_p().T
        self.xm = super(IntegralEquation, self).x_m().T
        self.ж = super(IntegralEquation, self).box().T
        self.N = len(self.ж[0])

    def Δx(self) -> array:
        Δx = zeros(self.N); Δy = zeros(self.N)
        for n in range(self.N):
            Δx[n] = self.xp[0][n] - self.xm[0][n]
            Δy[n] = self.xp[1][n] - self.xm[1][n]
        return Δx, Δy

    def normal_vector(self) -> array:
        '''
        Δx = x_m - x_p, Δy = y_m - y_p
        nhat = (Δy, -Δx)
        '''
        nx = zeros(self.N); ny = zeros(self.N)
        δx, δy = self.Δx()
        for n in range(self.N):
            denominator = sqrt(δx[n]**2 + δy[n]**2)
            nx[n] = -δy[n]/denominator
            ny[n] = δx[n]/denominator
        return nx, ny

    def dS(self) -> array:
        δx, δy = self.Δx()
        dS = zeros(self.N)
        for n in range(self.N):
            dS[n] = sqrt(δx[n]**2 + δy[n]**2)
        return dS