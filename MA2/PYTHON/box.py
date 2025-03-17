from numpy import array, zeros, linspace
def box(Nx: int, Ny: int, L: float, D: float) -> array:
    x = linspace(-L/2, L/2, Nx)
    y = linspace(0, -D, Ny)
    box = zeros((Nx + 2*Ny, 2))
    for n in range(Ny):
        box[n] = [x[0],y[n]]
        box[n + Nx + Ny] = [x[-1], y[-n-1]]
    for n in range(Nx):
        box[n + Ny] = [x[n], y[-1]]
    return box

box = box(10,10,2,1).T
print (box)
import matplotlib.pyplot as plt
plt.plot(box[0], box[1], '*'); plt.show()