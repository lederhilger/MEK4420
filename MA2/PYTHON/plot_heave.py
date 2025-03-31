from solve.excitation import Excitation
from solve.arguments import parse_args
from numpy import linspace
import matplotlib.pyplot as plt
if __name__ == "__main__":
    args = parse_args()
    Nx = args.Nx; Ny = args.Ny; D = args.D
    Ls = [.1, 1, 2, 10]
    for L in Ls:
        print(f"L = {L}")
        N = 20; i = 0
        kDs = linspace(0, 2, N+2)[1:-1]
        ξ2_INT = []; ξ2_H1 = []; ξ2_H2 = []; ξ2_FK = []
        ξ2_ROUGH = []; ξ2_ROUGH_CORRECTED = []
        for kD in kDs:
            print(f"{i}/{N}")
            init = Excitation(Nx, Ny, L, D, kD, 2)
            ξ2_INT.append(init.ξ2_full('integral'))
            ξ2_H1.append(init.ξ2_full('haskind1'))
            ξ2_H2.append(init.ξ2_full('haskind2'))
            ξ2_FK.append(init.ξ2_full('froudekrylov'))
            #ξ2_ROUGH.append(init.ξ2_rough())
            ξ2_ROUGH_CORRECTED.append(init.ξ2_rough_correction())
            i += 1
        #plt.plot(kDs, ξ2_ROUGH, '--', color = 'k', label = 'Rough')
        plt.plot(kDs, ξ2_ROUGH_CORRECTED, '-.', color = 'k', label = 'Rough corrected')
        plt.plot(kDs, ξ2_INT, '.', color = 'k', label = 'Integral')
        plt.plot(kDs, ξ2_H1, '*', color = 'k', label = 'Haskind 1')
        plt.plot(kDs, ξ2_H2, 'x', color = 'k', label = 'Haskind 2')
        plt.plot(kDs, ξ2_FK, '^', color = 'k', label = 'Froude--Krylov')
        plt.title(f"L = {L}"); plt.xlabel(r"$\kappa D$"); plt.ylabel(r"$\hat{\xi}$")
        plt.legend(); plt.show()
        print("--------")