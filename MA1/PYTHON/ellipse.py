from numpy import *
import matplotlib.pyplot as plt
from solve.integralequation import IntegralEquation
from solve.arguments import parse_args

def ellipse(a, b, N):
    θ_p = linspace(2*pi/N, 2*pi, N)
    θ_m = linspace(0, 2*pi*(N - 1)/N, N)
    x_p = zeros_like(θ_p); y_p = zeros_like(θ_p)
    x_m = zeros_like(θ_m); y_m = zeros_like(θ_m)
    for n in range(N):
        x_p[n], y_p[n] = a*cos(θ_p[n]), b*sin(θ_p[n])
        x_m[n], y_m[n] = a*cos(θ_m[n]), b*sin(θ_m[n])
    return (x_p, x_m), (y_p, y_m)

def test_convergence():
    abscissa = zeros(number)
    m_11 = zeros(number); m_22 = zeros(number); m_66 = zeros(number)
    for i in range(number):
        abscissa[i] = 1/(N*(i+1))
        init = IntegralEquation((i+1)*N, ellipse(a,b,(i+1)*N))
        m_11[i], m_22[i], m_66[i] = init.added_mass(); L2 = init.L2_norm()
        print(f'N = {N*(i+1)}')
        print(f'm_11: {m_11[i]/(pi*b**2)}'); print(f'm_22: {m_22[i]/(pi*a**2)}')
        if a == b:
            print(f'm_66: {1 - m_66[i]}')
        else:
            print(f'm_66: {m_66[i]/(.125*pi*(a**2 - b**2)**2)}')
        print(f'|phi_1|_L2 = {L2[0]}'); print(f'|phi_2|_L2 = {L2[1]}')
        print('------------------------------')
    plt.loglog(abscissa, m_11/(pi*b**2), '*', color = 'k', label = r'${m}_{11}$')
    plt.loglog(abscissa, m_22/(pi*a**2), 'x', color = 'k', label = r'${m}_{22}$')
    plt.xlabel(r'$\frac{1}{N}$'); plt.title('Added mass')
    if a == b:
        plt.loglog(abscissa, 1 - m_66, '.', color = 'k', label = r'${m}_{66}$')
        plt.legend()
        plt.savefig(f"addedmass_circle_N{N*number}.pdf", transparent = True, format = "pdf")
    else:
        plt.loglog(abscissa, m_66/(.125*pi*(a**2 - b**2)**2), '.', color = 'k', label = r'${m}_{66}$')
        plt.legend()
        plt.savefig(f"addedmass_a{a}_b{b}_N{N*number}.pdf", transparent = True, format = "pdf")
    plt.show()
    init.plot_phi()
    
if __name__ == "__main__":
    plt.rcParams['text.usetex'] = True
    args = parse_args()
    a = args.a; b = args.b; N = args.N; number = args.number
    test_convergence()