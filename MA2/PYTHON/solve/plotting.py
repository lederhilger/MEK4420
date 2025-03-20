from numpy import array, linspace, sqrt
from solve.integralequation import IntegralEquation

class Plotting(IntegralEquation):
    def __init__(self, Nx: int, Ny: int, L: float, D: float, kD: float):
        super(Plotting, self).__init__(Nx, Ny, L, D, kD)

    def plot_phi_k(self, mode: int):
        import matplotlib.pyplot as plt
        phi = self.assemble_k(mode)
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
        plt.plot(n, (lhs@phi_0).real, '*', label = 'lhs real', color = 'k')
        plt.plot(n, (rhs@phi_0_n).real, label = 'rhs real', color = 'k')
        #plt.legend(); plt.show()
        plt.plot(n, (lhs@phi_0).imag, 'x', label = 'lhs imag', color = 'k')
        plt.plot(n, (rhs@phi_0_n).imag, '-.', label = 'rhs imag', color = 'k')
        plt.legend(); plt.show()

    def plot_added_mass(self, mode: int):
        '''
            This should be plotted with respect to different κ
            Move to another module?
        '''
        import matplotlib.pyplot as plt
        a, b = self.added_mass(mode)
        a = a/(self.D**2); b = b/(self.D**2 * sqrt(9.8*self.κ))
        return a, b