from solve.integralequation import IntegralEquation
from numpy import *
box = IntegralEquation(100,100, 2, 1)
zhe = box.assemble_k('y')
box.plot_phi_k(zhe)
#import matplotlib.pyplot as plt
#plt.plot(zhe[0], zhe[1], '*'); plt.show()