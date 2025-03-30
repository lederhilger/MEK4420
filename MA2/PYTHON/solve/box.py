from numpy import array, zeros, ones, linspace, sqrt
from solve.chebyshov import Chebyshov
class Box:
    def __init__(self, Nx: int, Ny: int, L: float, D: float):
        self.Nx = Nx
        self.Ny = Ny
        self.L = L
        self.D = D
    
    def x_p(self) -> array:
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
    
    def x_m(self) -> array:
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

    def box(self) -> array:
        xp = self.x_p(); xm = self.x_m()
        M = self.Nx + 2*self.Ny - 3
        box = zeros((M, 2))
        for n in range(M):
            box[n] = [.5*(xp[n][0] + xm[n][0]), .5*(xp[n][1] + xm[n][1])]
        return box

    def Δx(self, N, xp, xm) -> array:
        Δx = zeros(N); Δy = zeros(N)
        for n in range(N):
            Δx[n] = xp[0][n] - xm[0][n]
            Δy[n] = xp[1][n] - xm[1][n]
        return Δx, Δy

    def normal_vector(self, N, Δx) -> array:
        '''
        Δx = x_m - x_p, Δy = y_m - y_p
        nhat = (Δy, -Δx)
        '''
        nx = zeros(N); ny = zeros(N)
        δx, δy = Δx
        for n in range(N):
            denominator = sqrt(δx[n]**2 + δy[n]**2)
            nx[n] = -δy[n]/denominator
            ny[n] = δx[n]/denominator
        return nx, ny

    def dS(self, N, Δx) -> array:
        δx, δy = Δx
        dS = zeros(N)
        for n in range(N):
            dS[n] = sqrt(δx[n]**2 + δy[n]**2)
        return dS