import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FuncFormatter
from numpy import pi, zeros, arctan2
from solve.jacobi import Jacobi
from solve.potentials import Potentials

class PlotConvergence:
    def __init__(self, shape: str, a: float, b: float, N: int, number: int, abscissa, phi, **kwargs):
        self.shape = shape
        self.a = a; self.b = b
        self.N = N; self.number = number
        self.abscissa = abscissa
        self.phi = phi
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

    def pi_axis(self):
        plt.gca().xaxis.set_major_formatter(FuncFormatter(lambda val, pos: r'{:.0g}$\pi$'.format(val/pi) if val != 0 else '0'))
        plt.gca().xaxis.set_major_locator(MultipleLocator(base = pi))
    
    def plot_phi(self, ж):
        N = int(self.abscissa[-1])
        plt.rcParams['text.usetex'] = True
        plt.title(f'$N = {self.N}$')
        # arctan2 \in (-pi, pi). Therefore domain -> domain-pi in potentials.py
        domain = zeros(N)
        for n in range(N):
            domain[n] = arctan2(ж[1][n], ж[0][n]) + pi
        sorted_domain = sorted(domain)
        init = Potentials(self.a, self.b, N, domain)
        count = 1
        if self.shape != 'circle' and self.shape != 'ellipse':
            for phi in self.phi:
                self.pi_axis()
                plt.plot(domain, phi, 'x', color = 'k', markersize = 2, label = r'$\phi_{\mathrm{num}}$')
                plt.legend()
                plt.savefig(f"phi_{count}_N{N}.pgf", transparent = True, format = "pgf")
                plt.show()
                count = count + count*count
        else:
            for phi in self.phi:
                self.pi_axis()
                coordinates = list(zip(domain, phi))
                coordinates = sorted(coordinates, key = lambda coordinate: coordinate[0])
                garbo, sorted_phi = list(zip(*coordinates))
                plt.plot(sorted_domain, sorted_phi, 'x', color = 'k', markersize = 2, label = r'$\phi_{\mathrm{num}}$')
                if count == 1:
                    if self.shape == 'ellipse': call = init.ellipse_1()
                    else: call = init.circle_1()
                    plt.plot(sorted_domain, call, color = 'k', label = r'$\phi_1$')
                elif count == 2:
                    if self.shape == 'ellipse': call = init.ellipse_2()
                    else: call = init.circle_2()
                    plt.plot(sorted_domain, call, color = 'k', label = r'$\phi_2$')
                plt.legend()
                plt.savefig(f"phi_{count}_N{N}.pgf", transparent = True, format = "pgf")
                plt.show()
                count = count + count*count
