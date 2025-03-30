from solve.integralequation import IntegralEquation
from solve.arguments import parse_args
from numpy import linspace
if __name__ == "__main__":
    import matplotlib.pyplot as plt
    args = parse_args()
    Nx = args.Nx; Ny = args.Ny; D = args.D
    Ls = [.1, 1, 2, 10]
    for L in Ls:
        print(f"L: {L}")
        N = 100; i = 0
        kDs = linspace(0, 2, N+2)[1:-1]
        X_integral = []; X_froudekrylov = []
        X_haskind1 = []; X_haskind2 = []
        for kD in kDs:
            init = IntegralEquation(Nx, Ny, L, D, kD)
            X_integral.append(abs(init.X_integral(2)))
            X_froudekrylov.append(abs(init.X_froudekrylov()))
            X_haskind1.append(abs(init.X_haskind1(2)))
            X_haskind2.append(abs(init.X_haskind2(2)))
            i += 1; print(f"{i}/{N}")
        plt.plot(kDs, X_integral, '.', label = 'integral', color = 'k')
        plt.plot(kDs, X_froudekrylov, '*', label = 'froude--krylov', color = 'k')
        plt.plot(kDs, X_haskind1, '--', label = 'haskind 1', color = 'k')
        plt.plot(kDs, X_haskind2, 'x', label = 'haskind 2', color = 'k')
        plt.title('Heave excitation force'); plt.legend(); plt.show()
        print("--------")