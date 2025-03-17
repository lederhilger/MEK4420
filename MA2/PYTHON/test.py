from solve.integralequation import IntegralEquation
from numpy import *
box = IntegralEquation(5,5, 2, 1)
zhe = box.lhs_k(); print(zhe)
#import matplotlib.pyplot as plt
#plt.plot(zhe[0], zhe[1], '*'); plt.show()