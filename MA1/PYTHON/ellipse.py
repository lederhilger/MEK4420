from numpy import *
import matplotlib.pyplot as plt
from solve.integralequation import IntegralEquation
from solve.arguments import parse_args
from solve.plot_convergence import PlotConvergence
from progress.bar import Bar

def radius(a, b, θ):
    denominator = sqrt(b**2 * (cos(θ)**2) + a**2 * (sin(θ)**2))
    radius = a*b/denominator
    return radius

def ellipse(a, b, N):
    θ_p = linspace(2*pi/N, 2*pi, N)
    θ_m = linspace(0, 2*pi*(N - 1)/N, N)
    x_p = zeros_like(θ_p); y_p = zeros_like(θ_p)
    x_m = zeros_like(θ_m); y_m = zeros_like(θ_m)
    r_p = radius(a,b,θ_p); r_m = radius(a,b,θ_m)
    for n in range(N):
        x_p[n], y_p[n] = r_p[n]*cos(θ_p[n]), r_p[n]*sin(θ_p[n])
        x_m[n], y_m[n] = r_m[n]*cos(θ_m[n]), r_m[n]*sin(θ_m[n])
        # x_p[n], y_p[n] = a*cos(θ_p[n]), b*sin(θ_p[n])
        # x_m[n], y_m[n] = a*cos(θ_m[n]), b*sin(θ_m[n])
    return (x_p, x_m), (y_p, y_m)

def test_convergence():
    abscissa = zeros(number)
    m_11 = zeros(number); m_22 = zeros(number); m_66 = zeros(number)
    bar = Bar('Calculating', max = number, fill='#', suffix='%(percent)d%% %(elapsed)ds')
    for i in range(number):
        abscissa[i] = N*(i+1)
        init = IntegralEquation((i+1)*N, ellipse(a,b,(i+1)*N))
        m_11[i], m_22[i], m_66[i] = init.added_mass()
        bar.next()
    bar.finish()
    if a == b: shape = 'circle'
    else: shape = 'ellipse'
    PlotConvergence(shape, a, b, N, number, abscissa, m11 = m_11, m22 = m_22, m66 = m_66).plot_added_mass()
    # init.plot_phi()
    # ax = plt.gca(); ax.set_aspect('equal')
    # ж, ч = init.ж()
    # for n in range(len(ж)):
    #     print(arctan2(ч[n], ж[n]) - arctan2(ч[n-1], ж[n-1]))
    # plt.plot(init.ж()[0], init.ж()[1], 'x'); plt.show()

if __name__ == "__main__":
    plt.rcParams['text.usetex'] = True
    args = parse_args()
    a = args.a; b = args.b; N = args.N; number = args.number
    test_convergence()