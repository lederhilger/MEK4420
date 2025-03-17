from numpy import array, zeros_like, zeros, sqrt, pi, arccos
from scipy import special, exp, log
from solve.box import Box

class IntegralEquation(Box):
    def __init__(self, Nx: int, Ny: int, L: float, D: float):
        super(IntegralEquation, self).__init__(Nx, Ny, L, D)
        self.xp = super(IntegralEquation, self).x_p().T
        self.xm = super(IntegralEquation, self).x_m().T
        self.ж = super(IntegralEquation, self).box().T
        self.N = len(self.ж[0])
        self.κ = 1.2/self.D

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

    def f_1(self, z):
        f_1 = -2*exp(z) * (special.exp1(z) + log(z) - log(-z))
        return f_1
    
    def f_2(self, z):
        f_2 = 2*pi*exp(z)
        return f_2

    def lhs_k(self) -> array:
        '''
            Assembles the left-hand side of the integral equation
            -pi phi_k + int_{S_{mathrm{B}}} phi_k partial_{nhat} green ,dee S
            = int_{S_{mathrm{B}}} green hat{n}_k ,dee S
            This is equation 101
        '''
        x_p, y_p = self.xp
        x_m, y_m = self.xm
        ж, ч = self.ж
        nx, ny = self.normal_vector()
        dS = self.dS()
        dΘ = zeros((self.N,self.N), dtype = 'complex_')
        tol = 1e-8
        for n in range(self.N):
            for m in range(self.N):
                a_x = x_m[m] - ж[n]
                a_y = y_m[m] - ч[n]
                b_x = x_p[m] - ж[n]
                b_y = y_p[m] - ч[n]
                if n == m:
                    argument_1 = -pi
                else:
                    argument_1 = log((a_x + 1j*a_y)/(b_x + 1j*b_y)).imag
                c_1 = y_m[m] - ч[n]
                c_2 = y_p[m] - ч[n]
                argument_2 = log((a_x + 1j*c_1)/(b_x + 1j*c_2))
                з = self.κ*(ч[m] + ч[n] - 1j*(ж[m] - ж[n]))
                argument_3 = (nx[m]*(self.f_1(з).imag + 1j*(self.f_2(з).imag)) + ny[m]*(self.f_1(з).real + 1j*(self.f_2(з).real)))*self.κ*dS[m]
                dΘ[n,m] = argument_1 + argument_2 + argument_3
        return dΘ