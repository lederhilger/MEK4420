from numpy import linspace, zeros, sqrt, pi, sin, cos
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
    θ = linspace(0, 2*pi, N+1)
    r = radius(a, b, θ)
    x = zeros(N+1); z = zeros(N+1)
    for n in range(N+1):
        x[n], z[n] = r[n]*cos(θ[n]), r[n]*sin(θ[n])
    return x, z

def test_convergence():
    abscissa = zeros(number)
    m_11 = zeros(number); m_22 = zeros(number); m_66 = zeros(number)
    bar = Bar('Calculating', max = number, fill='#', suffix='%(percent)d%% %(elapsed)ds')
    for i in range(number):
        abscissa[i] = N*(i+1)
        geometry = ellipse(a,b,(i+1)*N)
        init = IntegralEquation((i+1)*N, geometry)
        phi = init.solve()
        m_11[i], m_22[i], m_66[i] = init.added_mass(phi)
        bar.next()
    bar.finish()
    if a == b: shape = 'circle'
    else: shape = 'ellipse'
    init_plot = PlotConvergence(shape, a, b, N, number, abscissa, phi, m11 = m_11, m22 = m_22, m66 = m_66)
    init_plot.plot_added_mass()
    init_plot.plot_phi(init.ж)

if __name__ == "__main__":
    plt.rcParams['text.usetex'] = True
    args = parse_args()
    a = args.a; b = args.b; N = args.N; number = args.number
    test_convergence()
