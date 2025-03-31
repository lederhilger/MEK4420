from solve.integralequation import IntegralEquation
from scipy import exp, sin, sqrt

class Excitation(IntegralEquation):
    def __init__(self, Nx: int, Ny: int, L: float, D: float, kD: float, mode: int):
        super(Excitation, self).__init__(Nx, Ny, L, D, kD)
        self.mode = mode
        self.a, self.b = self.added_mass(self.mode)
        self.phi_D = self.assemble_D()
        self.phi = self.assemble_k(self.mode)

    def X_integral(self) -> float:
        nx, ny = self.normal_vector
        if self.mode == 1:
            nhat = nx
        elif self.mode == 2:
            nhat = ny
        else:
            raise ValueError("Choose mode 1 or 2.")
        X = 0
        for n in range(self.N):
            X += self.phi_D[n]*nhat[n]*self.dS[n]*(-1j)*self.ω #Eq. 119
        return X/(self.g*self.D)
    
    def X_haskind1(self) -> float:
        nx, ny = self.normal_vector
        if self.mode == 1:
            nhat = nx
        elif self.mode == 2:
            nhat = ny
        else:
            raise ValueError("Choose mode 1 or 2.")
        phi_0 = self.phi_0(); phi_0_n = self.phi_0_n()
        X = 0
        for n in range(self.N):
            X += (self.phi[n]*phi_0_n[n] - phi_0[n]*nhat[n])*self.dS[n]
        X *= 1j*self.ω
        return X/(self.g*self.D)

    def X_haskind2(self) -> float:
        A_neg = self.farfield_amplitudes(self.mode)[1][-1]
        X = 1j*self.g*A_neg
        return X/(self.g*self.D)
    
    def X_froudekrylov(self) -> float:
        X = 2*self.g*exp(-self.κ*self.D)*sin(.5*self.κ*self.L)/self.κ
        return X/(self.g*self.D)
    
    def resonance_frequency(self):
        ω_n = sqrt((self.g/self.D)/(1 + (self.a[-1]/(self.L*self.D))))
        return ω_n

    def ξ2_rough(self) -> float:
        X = self.X_froudekrylov()
        denominator = self.g*self.L - self.κ*self.g*( self.D*self.L ) + 1j*self.ω*self.b_22()
        ξ = abs(X)/abs(denominator)
        return abs(ξ)
    
    def ξ2_rough_correction(self) -> float:
        X = self.X_froudekrylov()
        denominator = self.g*self.L - self.κ*self.g*( self.D*self.L + self.a[-1]) + 1j*self.ω*self.b_22()
        ξ = abs(X)/abs(denominator)
        return abs(ξ)*self.g**2
    
    def ξ2_full(self, option: str) -> float:
        """
            Equation 128 for ξ_2/A
        """
        options = {
            "integral": self.X_integral(),
            "haskind1": self.X_haskind1(),
            "haskind2": self.X_haskind2(),
            "froudekrylov": self.X_froudekrylov()
        }
        print(f"X {option}")
        X = options[option]*self.g*self.D; print(f"X: {X}")
        c_22 = self.L; m = c_22*self.D
        denominator = (c_22 - self.κ*(m + self.a[-1]) + 1j*self.κ*self.b[-1])
        print(f"denominator: {denominator}")
        ξ = abs(X)/abs(denominator)
        return ξ