from numpy import *
from solve.chebyshov import Chebyshov

class IntegralEquation:
    def __init__(self, N, coordinates):
        self.N = N
        self.x = coordinates

    def ж(self) -> ndarray:
        ж = zeros((2,self.N))
        for n in range(self.N):
            ж[0][n], ж[1][n] = .5*(self.x[0][0][n]+self.x[0][1][n]), .5*(self.x[1][0][n]+self.x[1][1][n])
        return ж
    
    def Δx(self) -> ndarray:
        Δx = zeros(self.N); Δy = zeros(self.N)
        for n in range(self.N):
            Δx[n] = self.x[0][0][n] - self.x[0][1][n]
            Δy[n] = self.x[1][0][n] - self.x[1][1][n]
        return Δx, Δy

    def normal_vector(self) -> ndarray:
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

    def dS(self) -> ndarray:
        δx, δy = self.Δx()
        dS = zeros(self.N)
        for n in range(self.N):
            dS[n] = sqrt(δx[n]**2 + δy[n]**2)
        return dS
    
    def assemble(self) -> ndarray:
        x_p, x_m = self.x[0]
        y_p, y_m = self.x[1]
        ж, ч = self.ж()
        dΘ = zeros((self.N,self.N))
        for i in range(self.N):
            for j in range(self.N):
                if i == j:
                    dΘ[i,j] = -pi
                else:
                    dΘ[i,j] = -angle(complex(x_p[j]-ж[i], y_p[j]-ч[i])/complex(x_m[j]-ж[i], y_m[j]-ч[i]))
                    # dΘ[i,j] = arctan2((y_m[j]-ч[i]), (x_m[j]-ж[i])) - arctan2((y_p[j]-ч[i]), (x_p[j]-ж[i]))
                    # a_x = x_m[j] - ж[i]
                    # a_y = y_m[j] - ч[i]
                    # b_x = x_p[j] - ж[i]
                    # b_y = y_p[j] - ч[i]
                    # argument = (a_x*b_x + a_y*b_y)/sqrt((a_x**2 + a_y**2)*(b_x**2 + b_y**2))
                    # if argument > 1 and abs(argument - 1) < 1e-8:
                    #     dΘ[i,j] = -1*arccos(1)
                    # else:
                    #     dΘ[i,j] = -1*arccos(argument)
        return dΘ
    
    def assemble_h(self) -> ndarray:
        δx, δy = self.Δx()
        ж, ч = self.ж()
        dS = self.dS()
        h = zeros((self.N,self.N))
        for i in range(self.N):
            for j in range(self.N):
                x_m1, y_m1 = .5*δx[j]/sqrt(3) + ж[j], .5*δy[j]/sqrt(3) + ч[j]
                x_m2, y_m2 = -.5*δx[j]/sqrt(3) + ж[j], -.5*δy[j]/sqrt(3) + ч[j]
                h[i,j] = log((x_m1 - ж[i])**2 + (y_m1 - ч[i])**2)
                h[i,j] += log((x_m2 - ж[i])**2 + (y_m2 - ч[i])**2)
                h[i,j] *= .25*dS[j]
        return h
    
    def right_hs(self, mode: int) -> ndarray:
        n_x, n_y = self.normal_vector()
        ж, ч = self.ж()
        if mode == 1:
            n_i = n_x
        elif mode == 2:
            n_i = n_y
        elif mode == 6:
            n_i = zeros(self.N)
            for n in range(self.N):
                n_i[n] = ж[n]*n_y[n] - ч[n]*n_x[n]
        else:
            raise ValueError("Choose mode: 1, 2, 6")
        right_hs = dot(self.assemble_h(), n_i)
        return right_hs
    
    def solve(self):
        assemble = self.assemble()
        phi_1 = linalg.solve(assemble, self.right_hs(1))
        phi_2 = linalg.solve(assemble, self.right_hs(2))
        phi_6 = linalg.solve(assemble, self.right_hs(6))
        return phi_1, phi_2, phi_6
    
    def L2_norm(self):
        phi_1, phi_2, phi_6 = self.solve()
        theta = linspace(2*pi/self.N, 2*pi, self.N)
        mode_1 = sqrt(trapz(((phi_1 + cos(theta))**2), theta))
        mode_2 = sqrt(trapz(((phi_2 + sin(theta))**2), theta))
        return mode_1, mode_2

    def added_mass(self):
        m_11 = 0; m_22 = 0; m_66 = 0
        phi_1, phi_2, phi_6 = self.solve()
        nx, ny = self.normal_vector(); dS = self.dS()
        ж, ч = self.ж()
        for j in range(self.N):
            m_11 += phi_1[j]*nx[j]*dS[j]
            m_22 += phi_2[j]*ny[j]*dS[j]
            m_66 += phi_6[j]*(ж[j]*ny[j] - ч[j]*nx[j])*dS[j]
        return m_11, m_22, m_66
    
    def plot_phi(self):
        phi_1, phi_2, phi_6 = self.solve()
        import matplotlib.pyplot as plt
        plt.rcParams['text.usetex'] = True
        from matplotlib.ticker import MultipleLocator, FuncFormatter
        ax = plt.gca()
        ax.xaxis.set_major_formatter(FuncFormatter(lambda val,pos: r'{:.0g}$\pi$'.format(val/pi) if val !=0 else '0'))
        ax.xaxis.set_major_locator(MultipleLocator(base=pi))
        plt.title(f'$N = {self.N}$')
        dom = zeros(self.N); ж, ч = self.ж()
        for n in range(self.N):
            dom[n] = arctan2(ч[n],ж[n])
        count = 0
        for phi in [phi_1, phi_2, phi_6]:
            count += 1
            plt.plot(dom, phi, 'x', color = 'k', markersize = 2, label = r'$\phi_j$')
            plt.legend()
            plt.savefig(f"phi{count}_N{self.N}.pgf", transparent = True, format = "pgf")
            plt.show()

    def normal_plot(self):
        import matplotlib.pyplot as plt
        nx, ny = self.normal_vector(); ж, ч = self.ж()
        x_p, x_m = self.x[0]; y_p, y_m = self.x[1]
        for n in range(self.N):
            plt.plot((ж[n], nx[n]+ж[n]), (ч[n], ny[n]+ч[n]))
            plt.plot((x_p[n], x_m[n]), (y_p[n], y_m[n]))
        plt.plot(ж, ч, 'x')
        plt.show()