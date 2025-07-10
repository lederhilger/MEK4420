from numpy import *
import matplotlib.pyplot as plt
from solve.integralequation import IntegralEquation
from solve.arguments import parse_args
from solve.jacobi import Jacobi
from solve.chebyshov import Chebyshov
from solve.plot_convergence import PlotConvergence
from progress.bar import Bar

def square(a: float, N: int) -> ndarray:
    # x = Jacobi(linspace(-a, a, N + 1)).inverse_map()
    x = linspace(-a, a, N + 1)
    x_p = ones(4*N); x_m = ones(4*N); y_p = ones(4*N); y_m = ones(4*N)
    x_p[:N-1] *= a; y_p[:N-1] = x[1:N]
    x_p[N-1:2*N - 1] = flip(x)[:N]; y_p[N-1:2*N-1] *= a
    x_p[2*N - 1:3*N - 1] *= -a; y_p[2*N - 1:3*N - 1] = flip(x)[:N]
    x_p[3*N-1:4*N] = x; y_p[3*N-1:4*N] *= -a
    x_m[:N] *= a; y_m[:N] = x[:N]
    x_m[N:2*N] = flip(x)[:N]; y_m[N:2*N] *= a
    x_m[2*N:3*N] *= -a; y_m[2*N:3*N] = flip(x)[:N]
    x_m[3*N:4*N] = x[:N]; y_m[3*N:4*N] *= -a
    return (x_p, x_m), (y_p, y_m)

def test_square(N: int):
    M = 4*N
    abscissa = zeros(number)
    m_11 = zeros(number); m_22 = zeros(number); m_66 = zeros(number)
    k = 4*(2*Jacobi(None).K()**2/pi - 1)
    bar = Bar('Calculating', max = number, fill='#', suffix='%(percent)d%% %(elapsed)ds')
    for i in range(number):
        abscissa[i] = M*(i+1)
        init = IntegralEquation((i+1)*M, square(a, (i+1)*N))
        m_11[i], m_22[i], m_66[i] = init.added_mass()
        bar.next()
    bar.finish()
    PlotConvergence('square', a, a, M, number, abscissa, m11 = m_11, m22 = m_22, m66 = m_66).plot_added_mass()

if __name__ == "__main__":
    plt.rcParams['text.usetex'] = True
    args = parse_args()
    a = args.a; N = args.N; number = args.number
    test_square(N)
