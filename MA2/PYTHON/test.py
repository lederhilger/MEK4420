from solve.integralequation import IntegralEquation
from numpy import *
box = IntegralEquation(10,10, 2, 1)
zhe = box.dS(); print(zhe)
#import matplotlib.pyplot as plt
#plt.plot(zhe[0], zhe[1], '*'); plt.show()