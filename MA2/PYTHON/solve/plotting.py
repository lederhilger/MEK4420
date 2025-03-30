from numpy import linspace
from solve.integralequation import IntegralEquation
import matplotlib.pyplot as plt

class Plotting(IntegralEquation):
    def __init__(self, Nx: int, Ny: int, L: float, D: float, kD: float):
        super(Plotting, self).__init__(Nx, Ny, L, D, kD)

    def plt_box(self):
        box = self.box().T
        plt.axis('equal')
        plt.plot(box[0], box[1], 'x', color = 'k')
        plt.show()

    def plot_phi_k(self, mode: int):
        phi = self.assemble_k(mode)
        n = linspace(0, self.N, self.N)
        plt.plot(n, phi.real, '.', label = 'real', color = 'k')
        plt.plot(n, phi.imag, '*', label = 'imaginary', color = 'k')
        plt.legend(); plt.show()

    def plot_phi_0(self):
        n = linspace(0, self.N, self.N)
        phi_0 = self.phi_0()
        phi_0_n = self.phi_0_n()
        lhs = self.lhs_k(0); rhs = self.rhs_k()
        plt.plot(n, (lhs@phi_0).real, '*', label = 'lhs real', color = 'k')
        plt.plot(n, (rhs@phi_0_n).real, label = 'rhs real', color = 'k')
        plt.plot(n, (lhs@phi_0).imag, 'x', label = 'lhs imag', color = 'k')
        plt.plot(n, (rhs@phi_0_n).imag, '-.', label = 'rhs imag', color = 'k')
        plt.legend(); plt.show()

    def plot_phi_D(self):
        phi_D = self.assemble_D()
        n = linspace(0, self.N, self.N)
        plt.plot(n, phi_D.real, '.', label = 'real', color = 'k')
        plt.plot(n, phi_D.imag, '*', label = 'imaginary', color = 'k')
        plt.legend(); plt.show()