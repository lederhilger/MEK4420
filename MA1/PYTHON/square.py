from numpy import *
import matplotlib.pyplot as plt
from solve.integralequation import IntegralEquation
from solve.arguments import parse_args

def allocate_square(θ: tuple, x: tuple, y: tuple, a: float, N: int) -> tuple:
    for n in range(N):
        if 0 <= θ[n] < pi/4:
            x[n], y[n] = 2*a, 2*a*tan(θ[n])
        elif pi/4 <= θ[n] < 3*pi/4:
            x[n], y[n] = 2*a/tan(θ[n]), 2*a
        elif 3*pi/4 <= θ[n] < 5*pi/4:
            x[n], y[n] = -2*a, -2*a*tan(θ[n])
        elif 5*pi/4 <= θ[n] < 7*pi/4:
            x[n], y[n] = -2*a/tan(θ[n]), -2*a
        elif 7*pi/4 <= θ[n] <= 2*pi:
            x[n], y[n] = 2*a, 2*a*tan(θ[n])
        else:
            print(θ[n])
            raise ValueError("θ out of bounds")
    return x,y

def square(a: float, N: int) -> tuple:
    θ_p = linspace(2*pi/N, 2*pi, N)#; print(θ_p)
    θ_m = linspace(0, 2*pi*(N - 1)/N, N)#; print(θ_m)
    x_p, y_p = allocate_square(θ_p, zeros_like(θ_p), zeros_like(θ_p), a, N)
    x_m, y_m = allocate_square(θ_m, zeros_like(θ_m), zeros_like(θ_m), a, N)
    return (x_p, x_m), (y_p, y_m)

def test_square():
    abscissa = zeros(number)
    m_11 = zeros(number); m_22 = zeros(number); m_66 = zeros(number)
    for i in range(number):
        abscissa[i] = 1/(N*(i+1))
        init = IntegralEquation((i+1)*N, square(a, (i+1)*N))
        m_11[i], m_22[i], m_66[i] = init.added_mass()
        print(f'N = {N*(i+1)}')
        print(f'm_11 = {m_11[i]/(4.754*(2*a)**2)}'); print(f'm_22 = {m_22[i]/(4.754*(2*a)**2)}'); print(f'm_66 = {m_66[i]/(.725*(2*a)**4)}')
        print('------------------------------')
    plt.loglog(abscissa, m_11/(4.754*(2*a)**2), '*', color = 'k', label = r'${m}_{11}$')
    plt.loglog(abscissa, m_22/(4.754*(2*a)**2), 'x', color = 'k', label = r'${m}_{22}$')
    plt.loglog(abscissa, m_66/(.725*(2*a)**4), '.', color = 'k', label = r'${m}_{66}$')
    plt.xlabel(r'$\frac{1}{N}$')
    plt.title('Added mass'); plt.legend()
    plt.savefig(f"addedmass_square_N{N*number}.pdf", transparent = True, format = "pdf")
    plt.show()
    init.plot_phi()

if __name__ == "__main__":
    plt.rcParams['text.usetex'] = True
    args = parse_args()
    a = args.a; N = args.N; number = args.number
    test_square()