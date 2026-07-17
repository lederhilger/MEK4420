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

    def Iy(self) -> float:
        # From 3DFloat
        #Iz = 2.0496e8
        Iy = 2.1547e10
        return Iy

    def Izz(self) -> float:
        # From 3DFloat
        Izz = 7.3463e10
        return Izz

    def k(self) -> tuple:
        g = 9.81
        k_heave = 3*self.EA*sin(self.α)**2/self.L + self.c_22()
        k_surge = 2*self.EA*cos(self.α)**2/self.L
        k_yaw = 3*r_fairlead*pretension*cos(self.α)
        k_pitch = self.m*(y_B - y_G) + .25*self.ϱ*g*pi*r_2**4
        print(f"k_heave: {k_heave}"); print(f"k_surge: {k_surge}")
        return k_heave, k_surge
    
    def T(self) -> tuple:
        m_11 = self.m; m_22 = 2*pi*self.r_2**3/3
        k_heave, k_surge = self.k()
        T_heave = 2*pi*sqrt((self.m+m_22)/k_heave)
        T_surge = 2*pi*sqrt((self.m+m_11)/k_surge)
        return T_heave, T_surge

def T_yaw(L, h, pretension, r_fairlead) -> float:
    α = arcsin(h/L)
    # From 3DFloat
    I_y = 2.0496e8
    k_yaw = 3*r_fairlead*pretension*cos(α)
    print(f"k_yaw: {k_yaw}")
    T_yaw = 2*pi*sqrt(I_y/k_yaw)
    return T_yaw

def T_yaw_3DF(L, h, pretension, r_fairlead) -> float:
    α = arcsin(h/L)
    # From 3DFloat
    I_y = 2.0496e8
    k_yaw = 3*r_fairlead*pretension*cos(α)
    M_3DF = 1.352607e7
    k_3DF = M_3DF*18/pi
    I_3DF = 2.0496e8
    T_yaw_3DF = 2*pi*sqrt(I_3DF/(k_3DF+k_yaw))
    return T_yaw_3DF

def T_pitch(m, r) -> float:
    g = 9.81; ϱ = 1030; S_33 = .25*pi*r**4
    # From 3DFloat
    y_G = -65.704; y_B = -53.48
    I_zz = 7.3463e10
    I_zz0 = I_zz + m*y_G**2
    k_pitch = m*g*(y_B - y_G) + ϱ*g*S_33
    print(f"k_pitch: {k_pitch}")
    T_pitch = 2*pi*sqrt(I_zz0/k_pitch)
    return T_pitch
    
if __name__ == "__main__":
    import matplotlib.pyplot as plt
    h = 320
    r_1 = 5; r_2 = 9
    r_fairlead = 16
    m = 2.3864e7
    EA = 6e8
    #EA = 2e8
    pretension = 2e6
    y_B = -53.48; y_G = -65.704
    anchor_r = 50; rs = []
    T_H = []; T_S = []; T_Y = []; T_P = []
    while anchor_r < 1010:
        rs.append(anchor_r)
        L = sqrt(h**2 + anchor_r**2)
        init = Moor(h, L, r_1, r_2, EA, m)
        print(f"Anchor radius: {anchor_r}")
        T_heave, T_surge = init.T()
        #T_y = T_yaw(L, h, pretension, r_fairlead)
        T_y = T_yaw_3DF(L, h, pretension, r_fairlead)
        T_p = T_pitch(m, r_1)
        T_H.append(T_heave/20); T_S.append(T_surge/40)
        T_Y.append(T_y); T_P.append(T_p)
        print(f"T_heave: {T_heave}"); print(f"T_surge: {T_surge}")
        print(f"T_yaw: {T_y}"); print(f"T_pitch: {T_p}")
        print("-------------------------------------")
        anchor_r += 10
    plt.plot(rs, T_H, '.', color = 'k', label = 'Heave')
    plt.plot(rs, T_S, '*', color = 'k', label = 'Surge')
    #plt.plot(rs, T_Y, '^', color = 'k', label = 'Yaw')
    #plt.plot(rs, T_P, '-.', color = 'k', label = 'Pitch')
    plt.ylabel('Normalized period'); plt.xlabel('Anchor radius')
    plt.legend(); plt.show()
