from numpy import array, zeros, linspace
class Box:
    def __init__(self, Nx: int, Ny: int, L: float, D: float):
        self.Nx = Nx
        self.Ny = Ny
        self.L = L
        self.D = D

    

    def box(self) -> array:
        x = linspace(-self.L/2, self.L/2, self.Nx)
        y = linspace(0, -self.D, self.Ny)
        box = zeros((self.Nx + 2*self.Ny, 2))
        for n in range(self.Ny):
            box[n] = [x[0],y[n]]
            box[n + self.Nx + self.Ny] = [x[-1], y[-n-1]]
        for n in range(self.Nx):
            box[n + self.Ny] = [x[n], y[-1]]
        return box

box = box(10,10,2,1).T
print (box)
import matplotlib.pyplot as plt
plt.plot(box[0], box[1], '*'); plt.show()