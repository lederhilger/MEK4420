from solve.plotting import Plotting
from solve.arguments import parse_args
if __name__ == "__main__":
    args = parse_args()
    Nx = args.Nx; Ny = args.Ny; D = args.D; kD = args.kD
    Ls = [.1, 1, 2, 10]
    for L in Ls:
        Plotting(Nx, Ny, L, D, kD).plt_box()