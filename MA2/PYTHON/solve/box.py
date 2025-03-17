from numpy import array, zeros, ones, linspace
from solve.chebyshov import Chebyshov
class Box:
    def __init__(self, Nx: int, Ny: int, L: float, D: float):
        self.Nx = Nx
        self.Ny = Ny
        self.L = L
        self.D = D
    
    def x_p(self):
        x = Chebyshov(linspace(-self.L/2, self.L/2, self.Nx)).inverse_map()
        y_r = Chebyshov(linspace(0, -self.D, self.Ny)).inverse_map()
        y_l = y_r[1:]
        x_p = ones((self.Nx + 2*self.Ny - 3, 2)) #One node missing & removing duplicates in corner: -3
        for n in range(self.Ny - 1):
            x_p[n] = [x[0],y_l[n]]
        for n in range(1, self.Nx):
            x_p[n + self.Ny - 2] = [x[n], y_l[-1]]
        for n in range(self.Ny):
            x_p[self.Nx + self.Ny + n - 3] = [x[-1], y_r[-(n+1)]]
        return x_p
    
    def x_m(self):
        x = Chebyshov(linspace(-self.L/2, self.L/2, self.Nx)).inverse_map()
        y_l = Chebyshov(linspace(0, -self.D, self.Ny)).inverse_map()
        y_r = y_l[1:]
        x_m = zeros((self.Nx + 2*self.Ny - 3, 2))
        for n in range(self.Ny):
            x_m[n] = [x[0],y_l[n]]
        for n in range(self.Nx):
            x_m[n + self.Ny - 1] = [x[n], y_l[-1]]
        for n in range(1, self.Ny - 1):
            x_m[self.Nx + self.Ny + n - 2] = [x[-1], y_r[-(n+1)]]
        return x_m

    def box(self, Nx, Ny) -> array:
        x = linspace(-self.L/2, self.L/2, Nx)
        y = linspace(0, -self.D, Ny)
        box = zeros((Nx + 2*Ny, 2))
        for n in range(Ny):
            box[n] = [x[0],y[n]]
            box[Nx + Ny + n] = [x[-1], y[-(n+1)]]
        for n in range(Nx):
            box[n + Ny] = [x[n], y[-1]]
        return box