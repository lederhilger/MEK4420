from solve.jacobi import Jacobi
from solve.chebyshov import Chebyshov
import matplotlib.pyplot as plt
from scipy.special import ellipj
from numpy import linspace, zeros_like, ones, pi, sqrt

domain = linspace(0, 1, 101)

#sn, cn, dn, ph = ellipj(domain, .5)
#cd = cn/dn

jacob = Jacobi(domain)
chebyshov = Chebyshov(domain)

# K = pi/2 F(.5;.5;1;m)
# F(a;b;.5(a+b+1);.5) = sqrt(pi)/Gamma^2(3/4)
# K = (pi/2) sqrt(pi) / Gamma^2(3/4)
# Gamma(3/4) = sqrt( sqrt(pi/2) * AGM(sqrt(2); 1; .5*(sqrt(2)-1)) )
Kay = .5*pi*sqrt(2)/jacob.AGM(sqrt(2), 1, .5*(sqrt(2)-1))[0][-1]
K = jacob.K()
print(f"Difference in K: {abs(K - Kay)}")
K_MT = 1.854074677
print(f"Milne-Thomson difference in K: {abs(K - K_MT)}")

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