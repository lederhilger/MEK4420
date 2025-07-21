from numpy import linspace, ones, zeros, flip, pi
import matplotlib.pyplot as plt
from solve.integralequation import IntegralEquation
from solve.arguments import parse_args
from solve.jacobi import Jacobi
from solve.chebyshov import Chebyshov
from solve.plot_convergence import PlotConvergence
from progress.bar import Bar

def square(a: float, N: int) -> tuple:
    # line = Jacobi(linspace(-a, a, N + 1)).inverse_map()
    line = Chebyshov(linspace(-a, a, N + 1)).inverse_map()
    # line = linspace(-a, a, N + 1)
    x = ones(4*N + 1); z = ones(4*N + 1)
    x[:N] *= a; z[:N] = line[:N]
    x[N:2*N] = flip(line)[:N]; z[N:2*N] *= a
    x[2*N:3*N] *= -a; z[2*N:3*N] = flip(line)[:N]
    x[3*N:4*N] = line[:N]; z[3*N:4*N] *= -a
    x[-1] = x[0]; z[-1] = z[0]
    return x, z

def test_square(N: int):
    M = 4*N
    abscissa = zeros(number)
    m_11 = zeros(number); m_22 = zeros(number); m_66 = zeros(number)
    k = 4*(2*Jacobi(None).K()**2/pi - 1)
    bar = Bar('Calculating', max = number, fill='#', suffix='%(percent)d%% %(elapsed)ds')
    for i in range(number):
        abscissa[i] = M*(i+1)
        geometry = square(a, (i+1)*N)
        init = IntegralEquation((i+1)*M, geometry)
        phi = init.solve()
        m_11[i], m_22[i], m_66[i] = init.added_mass(phi)
        bar.next()
    bar.finish()
    init_plot = PlotConvergence('square', a, a, M, number, abscissa, phi, m11 = m_11, m22 = m_22, m66 = m_66)
    init_plot.plot_added_mass()
    init_plot.plot_phi(init.ж)

if __name__ == "__main__":
    plt.rcParams['text.usetex'] = True
    args = parse_args()
    a = args.a; N = args.N; number = args.number
    test_square(N)
