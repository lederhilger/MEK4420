from solve.plotting import Plotting
from solve.arguments import parse_args

if __name__ == "__main__":
    import matplotlib.pyplot as plt
    from numpy import linspace
    #plt.rcParams['text.usetex'] = True
    args = parse_args()
    Nx = args.Nx; Ny = args.Ny; L = args.L; D = args.D; kD = args.kD
    I = 20; interval = linspace(0, 2, I)[1:-1]
    for i in range(I-2):
        print(i)
        box = Plotting(Nx, Ny, L, D, interval[i])
        a, b = box.plot_added_mass(2)
        #plt.plot(interval[i], a[-1], 'x', color = 'k')
        plt.plot(interval[i], b[-1], '*', color = 'k')
        plt.plot(interval[i], box.b_22(), 'o', color = 'k')
    plt.show()