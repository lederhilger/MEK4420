from numpy import *
import matplotlib.pyplot as plt
from solve.integralequation import IntegralEquation
from solve.arguments import parse_args
from solve.chebyshov import Chebyshov

def square(a: float, N: int) -> ndarray:
    x = Chebyshov(linspace(-a, a, N//4 + 1)).inverse_map()
    x_p = ones(N); x_m = ones(N); y_p = ones(N); y_m = ones(N)
    x_p[:N//4] *= a; y_p[:N//4] = x[:N//4]
    x_p[N//4:int(N//4 * 2)] *= flip(x)[:N//4]; y_p[N//4:int(N//4 * 2)] *= a
    x_p[int(N//4 * 2):int(N//4 * 3)] *= -a; y_p[int(N//4 * 2):int(N//4 * 3)] = flip(x)[:N//4]
    x_p[int(N//4 * 3):N] = x[:N//4]; y_p[int(N//4 * 3):N] *= -a
    x_m[:N//4] *= a; y_m[:N//4] = x[1:N//4+1]
    x_m[N//4:int(N//4 * 2)] *= flip(x)[1:N//4+1]; y_m[N//4:int(N//4 * 2)] *= a
    x_m[int(N//4 * 2):int(N//4 * 3)] *= -a; y_m[int(N//4 * 2):int(N//4 * 3)] = flip(x)[1:N//4+1]
    x_m[int(N//4 * 3):N] = x[1:N//4+1]; y_m[int(N//4 * 3):N] *= -a
    return (x_p, x_m), (y_p, y_m)

def test_square():
    abscissa = zeros(number)
    m_11 = zeros(number); m_22 = zeros(number); m_66 = zeros(number)
    for i in range(number):
        abscissa[i] = 1/(N*(i+1))
        init = IntegralEquation((i+1)*N, square(a, (i+1)*N))
        m_11[i], m_22[i], m_66[i] = init.added_mass()
        print(f'N = {N*(i+1)}')
        print(f'm_11 = {m_11[i]/(4.754*(a)**2)}'); print(f'm_22 = {m_22[i]/(4.754*(a)**2)}'); print(f'm_66 = {m_66[i]/(.725*(a)**4)}')
        print('------------------------------')
    plt.loglog(abscissa, m_11/(4.754*a**2), '*', color = 'k', label = r'${m}_{11}$')
    plt.loglog(abscissa, m_22/(4.754*a**2), 'x', color = 'k', label = r'${m}_{22}$')
    plt.loglog(abscissa, m_66/(.725*a**4), '.', color = 'k', label = r'${m}_{66}$')
    plt.xlabel(r'$\frac{1}{N}$')
    plt.title('Added mass'); plt.legend()
    plt.savefig(f"addedmass_square_N{N*number}.pdf", transparent = True, format = "pdf")
    plt.show()
    init.plot_phi()
    init.normal_plot()

if __name__ == "__main__":
    plt.rcParams['text.usetex'] = True
    args = parse_args()
    a = args.a; N = args.N; number = args.number
    if N%8 != 0:
        raise ValueError("N must be divisible by 8.")
    test_square()