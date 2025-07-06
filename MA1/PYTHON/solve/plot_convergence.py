import matplotlib.pyplot as plt
from numpy import pi
from solve.jacobi import Jacobi

class PlotConvergence:
    def __init__(self, shape: str, a: float, b: float, N: int, number: int, abscissa, **kwargs):
        self.shape = shape
        self.a = a; self.b = b
        self.N = N; self.number = number
        self.abscissa = abscissa
        if kwargs:
            kwargs = {m: kwargs[m] for m in ['m11', 'm22', 'm66'] if m in kwargs}
            self.m = kwargs
        else:
            raise ValueError('Please input added mass')
        self.k = 4*(2*Jacobi(None).K()**2/pi - 1)

    def plot_added_mass(self):
        plt.xscale('log')
        plt.title('Added mass')
        plt.xlabel(r'$N$'); plt.ylabel(r'$\vert m_{\mathrm{num}} - m_{\mathrm{theory}}\vert$')
        real_m11 = {'circle': pi*self.b**2, 'ellipse': pi*self.b**2, 'square': self.k*self.a**2}
        real_m22 = {'circle': pi*self.a**2, 'ellipse': pi*self.a**2, 'square': self.k*self.a**2}
        real_m66 = {'circle': 0, 'ellipse': .125*pi*(self.a**2 - self.b**2)**2, 'square': .725*self.a**4}
        if 'm11' in self.m.keys():
            plt.plot(self.abscissa, abs(self.m['m11'] - real_m11[self.shape]), '*', color = 'k', label = r'${m}_{11}$')
        else: pass
        if 'm22' in self.m.keys():
            plt.plot(self.abscissa, abs(self.m['m22'] - real_m22[self.shape]), 'x', color = 'k', label = r'${m}_{22}$')
        else: pass
        if 'm66' in self.m.keys():
            plt.plot(self.abscissa, abs(self.m['m66'] - real_m66[self.shape]), '.', color = 'k', label = r'${m}_{66}$')
        else: pass
        plt.legend()
        plt.savefig(f"addedmass_{self.shape}_N{self.N*self.number}.pgf", transparent = True, format = "pgf")
        plt.show()