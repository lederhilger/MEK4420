from solve.integralequation import IntegralEquation
from numpy import *
box = IntegralEquation(100,100, 2, 1)
print(box.ж.T)
print(box.ж.T[200])
zhe = box.plot_phi_0()
#box.plot_phi_k(zhe)
#print(box.plot_added_mass(zhe))
import matplotlib.pyplot as plt
plt.plot(box.ж[0], box.ж[1], '*'); plt.show()