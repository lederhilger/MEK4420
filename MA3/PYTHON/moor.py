from numpy import arcsin, sqrt, sin, cos, pi
class Moor:
    def __init__(self, h, L, r_1, r_2, EA, m):
        self.L = L
        self.h = h
        self.α = arcsin(self.h/self.L)
        self.m = m
        self.r_1 = r_1
        self.r_2 = r_2
        self.ϱ = 1030
        self.EA = EA

    def c_22(self) -> float:
        g = 9.81; S = pi*self.r_1**2
        c_22 = self.ϱ*g*S
        return c_22

    def k(self) -> tuple:
        k_heave = 3*self.EA*sin(self.α)*cos(self.α)/self.L + self.c_22()
        k_surge = 2*self.EA*cos(self.α)**2/self.L
        return k_heave, k_surge
    
    def T(self) -> tuple:
        m_11 = self.m; m_22 = 2*pi*self.r_2**3/3
        k_heave, k_surge = self.k()
        T_heave = 2*pi*sqrt((self.m+m_22)/k_heave)
        T_surge = 2*pi*sqrt((self.m+m_11)/k_surge)
        return T_heave, T_surge
    
if __name__ == "__main__":
    import matplotlib.pyplot as plt
    h = 320
    r_1 = 5; r_2 = 9
    m = 2.3864e7
    EA = 6e8
    anchor_r = 800; rs = []
    T_H = []; T_S = []
    while anchor_r < 850:
        rs.append(anchor_r)
        L = sqrt(h**2 + anchor_r**2)
        init = Moor(h, L, r_1, r_2, EA, m)
        T_heave, T_surge = init.T()
        T_H.append(T_heave/20); T_S.append(T_surge/40)
        anchor_r += 1
    plt.plot(rs, T_H, '.', color = 'k', label = 'Heave')
    plt.plot(rs, T_S, '*', color = 'k', label = 'Surge')
    plt.ylabel('Normalized Period'); plt.xlabel('Anchor radius')
    plt.legend(); plt.show()
