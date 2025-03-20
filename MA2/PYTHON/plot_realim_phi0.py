from solve.plotting import Plotting
from solve.arguments import parse_args
#TESTING kD = 1.2, L/D = 2
if __name__ == "__main__":
    args = parse_args()
    Nx = args.Nx; Ny = args.Ny; L = args.L; D = args.D; kD = args.kD
    Plotting(Nx, Ny, L, D, kD).plot_phi_0()