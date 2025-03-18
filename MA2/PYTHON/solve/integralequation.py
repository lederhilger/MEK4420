from numpy import array, linspace, zeros, sqrt, pi, arccos, linalg
from scipy import special, exp, log
from solve.box import Box

class IntegralEquation(Box):
    def __init__(self, Nx: int, Ny: int, L: float, D: float, kD: float):
        super(IntegralEquation, self).__init__(Nx, Ny, L, D)
        self.xp = super(IntegralEquation, self).x_p().T
        self.xm = super(IntegralEquation, self).x_m().T
        self.ж = super(IntegralEquation, self).box().T
        self.N = len(self.ж[0])
        self.κ = kD/self.D

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
    
    def phi_0(self) -> array:
        ж, ч = self.ж
        phi_0 = exp(self.κ*(ч - 1j*ж))
        return phi_0
    
    def phi_0_n(self) -> array:
        nx, ny = self.normal_vector()
        phi_0 = self.phi_0()
        phi_0_n = zeros(len(phi_0), dtype = 'complex_')
        for n in range(len(phi_0)):
            phi_0_n[n] = self.κ*(ny[n] - 1j*nx[n])*phi_0[n]
        return phi_0_n

    def f_1(self, z):
        f_1 = -2*exp(z) * (special.exp1(z) + log(z) - log(-z))
        return f_1
    
    def f_2(self, z):
        f_2 = 2*pi*exp(z)
        return f_2

    def lhs_k(self, mode: int) -> array:
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
        dΘ = zeros((self.N, self.N), dtype = 'complex_')
        for n in range(self.N):
            for m in range(self.N):
                a_x = x_m[m] - ж[n]
                a_y = y_m[m] - ч[n]
                b_x = x_p[m] - ж[n]
                b_y = y_p[m] - ч[n]
                if n == m:
                    if mode == 0:
                        argument_1 = pi
                    else:
                        argument_1 = -pi
                else:
                    argument_1 = log((a_x + 1j*a_y)/(b_x + 1j*b_y)).imag
                c_1 = y_m[m] + ч[n]
                c_2 = y_p[m] + ч[n]
                argument_2 = log((a_x + 1j*c_1)/(b_x + 1j*c_2)).imag
                з = self.κ*(ч[m] + ч[n] - 1j*(ж[m] - ж[n]))
                argument_3 = (nx[m]*(self.f_1(з).imag + 1j*(self.f_2(з).imag)) + ny[m]*(self.f_1(з).real + 1j*(self.f_2(з).real)))*self.κ*dS[m]
                dΘ[n,m] = argument_1 + argument_2 + argument_3
        return dΘ

    def rhs_k(self) -> array:
        '''
            Assembles the laft-hand side of the integral equation
            -pi phi_k + int_{S_{mathrm{B}}} phi_k partial_{nhat} green ,dee S
            = int_{S_{mathrm{B}}} green hat{n}_k ,dee S
            
            g_0 is a two point Gauss method
            g_1 and g_2 is a midpoint method

            The full right-hand side, int green hat{n}_k, is found by multiplying dΓ with nk
        '''
        δx, δy = self.Δx()
        dS = self.dS()
        ж, ч = self.ж
        gauss_x_pos = ж + .5*δx/sqrt(3)
        gauss_x_neg = ж - .5*δx/sqrt(3)
        gauss_y_pos = ч + .5*δy/sqrt(3)
        gauss_y_neg = ч - .5*δy/sqrt(3)
        dΓ = zeros((self.N, self.N), dtype = 'complex_')
        for n in range(self.N):
            for m in range(self.N):
                a_x = gauss_x_neg[m] - ж[n]
                a_y = gauss_y_neg[m] - ч[n]
                b_x = gauss_x_pos[m] - ж[n]
                b_y = gauss_y_pos[m] - ч[n]
                g_0 = .5*(log(sqrt(a_x**2 + a_y**2)) + log(sqrt(b_x**2 + b_y**2)))
                g_1 = -log(sqrt((ж[m] - ж[n])**2 + (ч[m] + ч[n])**2))
                з = self.κ*(ч[m] + ч[n] - 1j*(ж[m] - ж[n]))
                g_2 = self.f_1(з).real + 1j*(self.f_2(з).real)
                dΓ[n,m] = (g_0 + g_1 + g_2)*dS[m]
        return dΓ

    def assemble_k(self, mode: int) -> array:
        nx, ny = self.normal_vector()
        lhs = self.lhs_k(mode)
        if mode == 1:
            rhs = self.rhs_k()@nx
        elif mode == 2:
            rhs = self.rhs_k()@ny
        else:
            raise ValueError("Choose modes 1 or 2.")
        phi_k = linalg.solve(lhs,rhs)
        return phi_k
    
    def added_mass_y(self, phi: array) -> array:
        nx, ny = self.normal_vector()
        dS = self.dS()
        a = zeros(len(phi)); b = zeros(len(phi))
        a[0] = (phi[0]*ny[0]*dS[0]).real
        b[0] = (phi[0]*ny[0]*dS[0]).imag
        for n in range(1, len(phi)):
            a[n] = a[n-1] + (phi[n]*ny[n]*dS[n]).real
            b[n] = b[n-1] + (phi[n]*ny[n]*dS[n]).imag
        return a, b
    
    def plot_phi_k(self, phi: array):
        import matplotlib.pyplot as plt
        n = linspace(0, self.N, self.N)
        plt.plot(n, phi.real, label = 'real')
        plt.plot(n, phi.imag, label = 'imaginary')
        plt.legend(); plt.show()

    def plot_phi_0(self):
        import matplotlib.pyplot as plt
        n = linspace(0, self.N, self.N)
        phi_0 = self.phi_0()
        phi_0_n = self.phi_0_n()
        lhs = self.lhs_k(0)
        rhs = self.rhs_k()
        plt.plot(n, (lhs@phi_0).real, '*', label = 'lhs')
        plt.plot(n, (rhs@phi_0_n).real, label = 'rhs')
        plt.legend(); plt.show()
        plt.plot(n, (lhs@phi_0).imag, '*', label = 'lhs')
        plt.plot(n, (rhs@phi_0_n).imag, label = 'rhs')
        plt.legend(); plt.show()

    def plot_added_mass(self, phi: array):
        '''
            This should be plotted with respect to different κ
            Move to another module
        '''
        import matplotlib.pyplot as plt
        a, b = self.added_mass_y(phi)
        a = a/(self.D**2); b = b/(self.D**2 * sqrt(9.8*self.κ))
        return a,b