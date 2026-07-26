from numpy import linspace, zeros, sqrt, pi, sin, cos, append, flip
from solve.integralequation import IntegralEquation
from solve.arguments import parse_args
from solve.plot_convergence import PlotConvergence, relative_error
from progress.bar import Bar
from solve.jacobi import Jacobi

def radius(a, b, θ):
    denominator = sqrt(b**2 * (cos(θ)**2) + a**2 * (sin(θ)**2))
    radius = a*b/denominator
    return radius

def ellipse(a, b, N, parametrization="eccentric"):
    θ = linspace(0, 2*pi, N+1)
    if parametrization == "polar":
        r = radius(a, b, θ)
        return r*cos(θ), r*sin(θ)
    if parametrization == "eccentric":
        return a*cos(θ), b*sin(θ)
    raise ValueError("Incorrect parametrization.")

def equal_ellipse(a, b, N):
    η = Jacobi(
        domain = linspace(0, .5*pi, N, endpoint = False),
        modulus = sqrt(1 - (b/a)**2)
    ).equal_arc()
    η = append(η, .5*pi)
    xquad = a*cos(η); zquad = b*sin(η)
    x = zeros(4*N + 1); z = zeros(4*N + 1)
    x[:N] = xquad[:N]; z[:N] = zquad[:N]
    x[N:2*N] = -flip(xquad)[:N]
    z[N:2*N] = flip(zquad)[:N]
    x[2*N:3*N] = -xquad[:N]
    z[2*N:3*N] = -zquad[:N]
    x[3*N:4*N] = flip(xquad)[:N]
    z[3*N:4*N] = -flip(zquad)[:N]
    x[-1] = x[0]
    z[-1] = z[0]
    return x, z

def test_convergence(parametrization, Ns, plot=False):
    results = []
    number = len(Ns)
    abscissa = zeros(number)
    m_11 = zeros(number); m_22 = zeros(number); m_66 = zeros(number)
    M_11 = pi*b**2; M_22 = pi*a**2; M_66 = .125*pi*(a**2 - b**2)**2
    bar = Bar('Calculating', max = number, fill='#', suffix='%(percent)d%% %(elapsed)ds')
    for i, N in enumerate(Ns):
        if parametrization == "equal_arc":
            if N%4!=0: raise ValueError("4 must divide N.")
            geometry = equal_ellipse(a, b, N//4)
        else:
            geometry = ellipse(a, b, N, parametrization)
        abscissa[i] = N
        init = IntegralEquation(N, geometry, order)
        phi = init.solve()
        m_11[i], m_22[i], m_66[i] = init.added_mass(phi)
        results.append({
            "parametrization":parametrization,
            "N":int(N),
            "m_11":float(m_11[i]),
            "m_22":float(m_22[i]),
            "m_66":float(m_66[i]),
            "relm_11":float(relative_error(m_11[i], M_11)),
            "relm_22":float(relative_error(m_22[i], M_22)),
            "relm_66":(
                None if M_66 == 0
                else float(relative_error(m_66[i], M_66))
            )
        })
        bar.next()
    bar.finish()
    if plot:
        shape = "circle" if a == b else "ellipse"
        plotter = PlotConvergence(
            shape = shape,
            a = a,
            b = b,
            N = int(Ns[0]),
            number = len(Ns),
            abscissa = abscissa,
            phi = phi,
            parametrization = parametrization,
            m11 = m_11,
            m22 = m_22,
            m66 = m_66
        )
        plotter.plot_added_mass()
        plotter.plot_phi(init.ж)
    return results

if __name__ == "__main__":
    args = parse_args()
    a = args.a; b = args.b
    n = args.N; number = args.number; order = args.order
    parametrizations = ("polar", "eccentric", "equal_arc")
    Ns = tuple(args.N*(i+1) for i in range(args.number))
    results = []
    plotting = args.plot == "yes"
    for parametrization in parametrizations:
        results.extend(test_convergence(parametrization, Ns, plotting))
    import json
    data = {
        "a":a,
        "b":b,
        "order":order,
        "results":results
    }
    with open("ellipse_convergence.json", "w") as output:
        json.dump(data, output, indent=2)
