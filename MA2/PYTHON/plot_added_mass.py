from solve.integralequation import IntegralEquation
from numpy import linspace, sqrt
from solve.arguments import parse_args
if __name__ == "__main__":
    import matplotlib.pyplot as plt
    args = parse_args()
    Nx = args.Nx; Ny = args.Ny; D = args.D
    Ls = [.1, 1, 2, 10]
    for L in Ls:
        N = 100; i = 0
        kDs = linspace(0, 2, N+2)[1:-1]
        A = []; B = []; B_approx = []
        for kD in kDs:
            print(f"{i}/{N}")
            init = IntegralEquation(Nx, Ny, L, D, kD)
            a, b = init.added_mass(2)
            b_approx = init.b_22()
            A.append(a[-1]); B.append(b[-1]); B_approx.append(b_approx)
            i += 1
        plt.plot(kDs, A, '.', color = 'k', label = 'Added mass')
        plt.plot(kDs, B, '*', color = 'k', label = 'Damping')
        plt.plot(kDs, B_approx, color = 'k', label = 'Approx. damping')
        plt.title(f"L = {L}")
        plt.legend(); plt.show()