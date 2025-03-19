from numpy import array, zeros, sqrt, pi, linalg
from scipy import special, exp, log, sin
from solve.box import Box

class IntegralEquation(Box):
    def __init__(self, Nx: int, Ny: int, L: float, D: float, kD: float):
        super(IntegralEquation, self).__init__(Nx, Ny, L, D)
        self.xp = super(IntegralEquation, self).x_p().T
        self.xm = super(IntegralEquation, self).x_m().T
        self.ж = super(IntegralEquation, self).box().T
        self.N = len(self.ж[0])
        self.Δx = super(IntegralEquation, self).Δx(self.N, self.xp, self.xm)
        self.normal_vector = super(IntegralEquation, self).normal_vector(self.N, self.Δx)
        self.dS = super(IntegralEquation, self).dS(self.N, self.Δx)
        self.κ = kD/self.D
    
    def phi_0(self) -> array:
        ж, ч = self.ж
        phi_0 = exp(self.κ*(ч - 1j*ж))
        return phi_0
    
    def phi_0_n(self) -> array:
        nx, ny = self.normal_vector
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
        nx, ny = self.normal_vector
        dS = self.dS
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
        δx, δy = self.Δx
        dS = self.dS
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
        nx, ny = self.normal_vector
        lhs = self.lhs_k(mode)
        if mode == 1:
            rhs = self.rhs_k()@nx
        elif mode == 2:
            rhs = self.rhs_k()@ny
        else:
            raise ValueError("Choose mode 1 or 2.")
        phi_k = linalg.solve(lhs,rhs)
        return phi_k

    def assemble_D(self, mode: int = 1) -> array:
        lhs = self.lhs_k(mode); rhs = -2*pi*self.phi_0()
        phi_D = linalg.solve(lhs, rhs)
        return phi_D
    
    def added_mass(self, mode: int) -> array:
        phi = self.assemble_k(mode)
        nx, ny = self.normal_vector
        dS = self.dS
        a = zeros(len(phi)); b = zeros(len(phi))
        if mode == 1:
            nhat = nx
        elif mode == 2:
            nhat = ny
        else:
            raise ValueError("Choose mode 1 or 2.")
        a[0] = (phi[0]*nhat[0]*dS[0]).real
        b[0] = -(phi[0]*nhat[0]*dS[0]).imag
        for n in range(1, len(phi)):
            a[n] = a[n-1] + (phi[n]*nhat[n]*dS[n]).real
            b[n] = b[n-1] - (phi[n]*nhat[n]*dS[n]).imag
        return a, b
    
    def farfield_amplitudes(self, mode: int) -> tuple:
        phi_0 = self.phi_0(); phi = self.assemble_k(mode)
        nx, ny = self.normal_vector
        dS = self.dS
        A_pos = zeros(len(phi_0), dtype = 'complex_')
        A_neg = zeros(len(phi_0), dtype = 'complex_')
        A_pos[0] = 1j*(phi[0]*self.κ*(ny[0] + 1j*nx[0]) - ny[0])*phi_0[0].conjugate()*dS[0]
        A_neg[0] = 1j*(phi[0]*self.κ*(ny[0] - 1j*nx[0]) - ny[0])*phi_0[0]*dS[0]
        for n in range(1, len(phi_0)):
            A_pos[n] = A_pos[n-1] + 1j*(phi[n]*self.κ*(ny[n] + 1j*nx[n]) - ny[n])*phi_0[n].conjugate()*dS[n]
            A_neg[n] = A_neg[n-1] + 1j*(phi[n]*self.κ*(ny[n] - 1j*nx[n]) - ny[n])*phi_0[n]*dS[n]
        return A_pos, A_neg

    def b_22(self, mode = 2) -> float:
        ω = sqrt(9.8*self.κ)
        A_pos, A_neg = self.farfield_amplitudes(mode)
        A_1 = abs(A_pos[-1])**2; A_2 = abs(A_neg[-1])**2
        b_22 = .5*ω*(A_1 + A_2)/(self.D**2 * sqrt(9.8*self.κ)) # Scaled wrt ω D^2
        return b_22
    
    def X_integral(self, mode) -> float:
        phi_D = self.assemble_D()
        nx, ny = self.normal_vector
        if mode == 1:
            nhat = nx
        elif mode == 2:
            nhat = ny
        else:
            raise ValueError("Choose mode 1 or 2.")
        X = 0
        for n in range(self.N):
            X += phi_D[n]*nhat[n]*self.dS[n] # *(-1j)*sqrt(9.8*self.κ) #Eq. 119
        return X # Not scaled: /(9.8*self.D) (?)
    
    def X_haskind(self, mode) -> float:
        A_neg = self.farfield_amplitudes(mode)[1][-1]
        X = 1j*9.8*A_neg
        return X/(9.8*self.D)
    
    def X_froudekrylov(self) -> float:
        X = 9.8*2*self.L*exp(-self.κ*self.D)*sin(.5*self.κ*self.L)/(self.κ*self.L)
        return X/(9.8*self.D)