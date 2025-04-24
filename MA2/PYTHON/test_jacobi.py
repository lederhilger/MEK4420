from solve.jacobi import Jacobi
from solve.chebyshov import Chebyshov
import matplotlib.pyplot as plt
from scipy.special import ellipj
from numpy import linspace, zeros_like, ones

domain = linspace(0, 1, 101)

#sn, cn, dn, ph = ellipj(domain, .5)
#cd = cn/dn

jacob = Jacobi(domain)
chebyshov = Chebyshov(domain)

x_jac = jacob.inverse_map()
x_cheb = chebyshov.inverse_map()

y = zeros_like(domain)
why = ones(len(y))
plt.plot(x_cheb, y, '.', color = 'k', label = 'Chebyshov')
plt.plot(x_jac, why, '*', color = 'k', label = 'Jacobi')
plt.legend(); plt.show()

jac_diff = []
cheb_diff = []
ens = []
for n in range(len(x_jac)-1):
    jac_diff.append(abs(x_jac[n] - x_jac[n+1]))
    cheb_diff.append(abs(x_cheb[n] - x_cheb[n+1]))
    ens.append(n)

plt.plot(ens, jac_diff, '*', color = 'k', label = "cd(x)")
plt.plot(ens, cheb_diff, '.', color = 'k', label = "cos(x)")
plt.title("Jacobi vs Chebyshov nodes")
plt.legend(); plt.show()