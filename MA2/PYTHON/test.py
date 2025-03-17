from solve.box import Box

box = Box(10,10, 2, 1)
zhe = box.box().T#; print(zhe)
import matplotlib.pyplot as plt
plt.plot(zhe[0], zhe[1], '*'); plt.show()