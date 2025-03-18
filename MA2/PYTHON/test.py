from solve.integralequation import IntegralEquation
from solve.arguments import parse_args

if __name__ == "__main__":
    #plt.rcParams['text.usetex'] = True
    args = parse_args()
    Nx = args.Nx; Ny = args.Ny; L = args.L; D = args.D; kD = args.kD
    box = IntegralEquation(Nx, Ny, L, D, kD)
    box.plot_phi_0()
    box.plot_phi_k(box.assemble_k(2))