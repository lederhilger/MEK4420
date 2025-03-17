from solve.box import Box

box1 = Box(10,10,2,1).x_p().T
box2 = Box(10,10,2,1).x_m().T
import matplotlib.pyplot as plt
plt.plot(box1[0], box1[1], '*'); plt.show()
plt.plot(box2[0], box2[1], '*'); plt.show()