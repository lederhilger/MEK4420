from numpy import (ndarray, linalg, sqrt, pi,
                    zeros, angle, linspace, trapz)
from solve.chebyshov import Chebyshov
from solve.quadrature import Quadrature
from solve.potentials import Potentials

class IntegralEquation:
    def __init__(self, N, coordinates, order):
        self.N = N
        self.x = coordinates[0]
        self.z = coordinates[1]
        ж = zeros((2,self.N))
        for n in range(self.N):
            ж[0][n], ж[1][n] = .5*(self.x[n+1]+self.x[n]), .5*(self.z[n+1]+self.z[n])
        self.ж = ж
        self.order = order
    
    def Δx(self) -> ndarray:
        Δx = zeros(self.N); Δz = zeros(self.N)
        for n in range(self.N):
            Δx[n] = self.x[n+1] - self.x[n]
            Δz[n] = self.z[n+1] - self.z[n]
        return Δx, Δz

    def normal_vector(self) -> ndarray:
        '''
        Δx = x_m+1 - x_m, Δz = z_m+1 - z_m
        nhat = (Δz, -Δx)
        '''
        nx = zeros(self.N); nz = zeros(self.N)
        δx, δz = self.Δx()
        for n in range(self.N):
            denominator = sqrt(δx[n]**2 + δz[n]**2)
            nx[n] = -δz[n]/denominator
            nz[n] = δx[n]/denominator
        return nx, nz

    def dS(self) -> ndarray:
        δx, δz = self.Δx()
        dS = zeros(self.N)
        for n in range(self.N):
            dS[n] = sqrt(δx[n]**2 + δz[n]**2)
        return dS
    
    def assemble(self) -> ndarray:
        ж, ч = self.ж
        dΘ = zeros((self.N,self.N))
        for i in range(self.N):
            for j in range(self.N):
                if i == j:
                    dΘ[i,j] = -pi
                else:
                    dΘ[i,j] = -angle(complex(self.x[j+1]-ж[i], self.z[j+1]-ч[i])/complex(self.x[j]-ж[i], self.z[j]-ч[i]))
                    # dΘ[i,j] = arctan2((self.z[j]-ч[i]), (self.x[j]-ж[i])) - arctan2((self.z[j+1]-ч[i]), (self.x[j+1]-ж[i]))
        return dΘ
    
    def assemble_h(self) -> ndarray:
        δx = self.Δx()
        ж = self.ж
        h = Quadrature(self.N, δx, ж, 'Lagrange', self.order).quad()
        return h
    
    def right_hs(self, assemble_h, mode: int) -> ndarray:
        n_x, n_z = self.normal_vector()
        ж, ч = self.ж
        if mode == 1:
            n_i = n_x
        elif mode == 2:
            n_i = n_z
        elif mode == 6:
            n_i = zeros(self.N)
            for n in range(self.N):
                n_i[n] = ж[n]*n_z[n] - ч[n]*n_x[n]
        else:
            raise ValueError("Choose mode: 1, 2, 6")
        right_hs = assemble_h @ n_i
        return right_hs
    
    def solve(self):
        assemble = self.assemble()
        h = self.assemble_h()
        phi_1 = linalg.solve(assemble, self.right_hs(h, 1))
        phi_2 = linalg.solve(assemble, self.right_hs(h, 2))
        phi_6 = linalg.solve(assemble, self.right_hs(h, 6))
        return phi_1, phi_2, phi_6

    def added_mass(self, phi):
        m_11 = 0; m_22 = 0; m_66 = 0
        phi_1, phi_2, phi_6 = phi
        nx, nz = self.normal_vector(); dS = self.dS()
        ж, ч = self.ж
        for j in range(self.N):
            m_11 += phi_1[j]*nx[j]*dS[j]
            m_22 += phi_2[j]*nz[j]*dS[j]
            m_66 += phi_6[j]*(ж[j]*nz[j] - ч[j]*nx[j])*dS[j]
        return m_11, m_22, m_66

    def normal_plot(self):
        import matplotlib.pyplot as plt
        nx, nz = self.normal_vector(); ж, ч = self.ж()
        for n in range(self.N):
            plt.plot((ж[n], nx[n]+ж[n]), (ч[n], nz[n]+ч[n]))
            plt.plot((self.x[n+1], self.x[n]), (self.z[n+1], self.z[n]))
        plt.plot(ж, ч, 'x')
        plt.show()
