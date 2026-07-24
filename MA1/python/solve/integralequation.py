from numpy import (ndarray, linalg, sqrt, pi, angle,
                    zeros, linspace, diff, fill_diagonal, column_stack)
from solve.chebyshov import Chebyshov
from solve.quadrature import Quadrature
from solve.potentials import Potentials
from functools import cached_property

class IntegralEquation:
    def __init__(self, N, coordinates, order):
        self.N = N
        self.x = coordinates[0]
        self.z = coordinates[1]
        self.ж = zeros((2,self.N))
        self.ж[0] = .5*(self.x[1:]+self.x[:-1])
        self.ж[1] = .5*(self.z[1:]+self.z[:-1])
        self.order = order

    @cached_property
    def Δx(self) -> ndarray:
        Δx = diff(self.x)
        Δz = diff(self.z)
        return Δx, Δz

    @cached_property
    def dS(self) -> ndarray:
        δx, δz = self.Δx
        dS = sqrt(δx**2 + δz**2)
        return dS

    @cached_property
    def normal_vector(self) -> ndarray:
        '''
        Δx = x_m+1 - x_m, Δz = z_m+1 - z_m
        nhat = (Δz, -Δx)
        '''
        δx, δz = self.Δx
        denominator = self.dS
        nx = δz/denominator
        nz = -δx/denominator
        return nx, nz
    
    def assemble(self) -> ndarray:
        Ж = self.ж[0] + 1j*self.ж[1]
        X = self.x + 1j*self.z
        start = X[:-1]; end = X[1:]
        numerator = end[None,:] - Ж[:,None]
        denominator = start[None,:] - Ж[:,None]
        dΘ = -angle(numerator/denominator)
        fill_diagonal(dΘ, -pi)
        return dΘ
    
    def assemble_h(self) -> ndarray:
        δx = self.Δx
        ж = self.ж
        h = Quadrature(self.N, δx, ж, 'Lagrange', self.order).quad()
        return h
    
    def right_hs(self, assemble_h) -> ndarray:
        n_x, n_z = self.normal_vector
        ж, ч = self.ж
        modes = column_stack((n_x, n_z, ж*n_z - ч*n_x))
        right_hs = assemble_h @ modes
        return right_hs
    
    def solve(self):
        assemble = self.assemble()
        h = self.assemble_h()
        right_hs = self.right_hs(h)
        phi = linalg.solve(assemble, right_hs)
        return phi[:,0], phi[:,1], phi[:,2]

    def added_mass(self, phi):
        m_11 = 0; m_22 = 0; m_66 = 0
        phi_1, phi_2, phi_6 = phi
        nx, nz = self.normal_vector; dS = self.dS
        ж, ч = self.ж
        m_11 = (phi_1*nx*dS).sum()
        m_22 = (phi_2*nz*dS).sum()
        m_66 = (phi_6*(ж*nz - ч*nx)*dS).sum()
        return m_11, m_22, m_66

    def normal_plot(self):
        import matplotlib.pyplot as plt
        nx, nz = self.normal_vector; ж, ч = self.ж
        fig, ax = plt.subplots()
        ax.plot(self.x, self.z, color='k')
        ax.scatter(ж, ч, marker='x')
        ax.quiver(ж, ч, nx, nz, angles = 'xy', scale_units = 'xy', scale = 1)
        ax.set_aspect('equal')
        plt.show()
