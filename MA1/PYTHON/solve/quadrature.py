from numpy import sqrt, zeros, real, log, abs

class Quadrature:
    def __init__(self, N, δx, ж, method: str, order: int):
        self.N = N
        self.δx = [complex(*z) for z in zip(δx[0],δx[1])]
        self.ж = [complex(*z) for z in zip(ж[0],ж[1])]
        self.dS = [abs(z) for z in self.δx]
        self.method = "Lagrange"
        if order != 2 and order != 4:
            raise ValueError("Order 2 or 4")
        else:
            self.order = order

    def weights(self):
        if self.order == 2: return [1, 1]
        elif self.order == 4:
            a = (18+sqrt(30))/36
            b = (18-sqrt(30))/36
            return [a, b, b, a]

    def nodes(self):
        if self.order == 2: return [-1/sqrt(3), 1/sqrt(3)]
        elif self.order == 4:
            a = sqrt(3/7 - 2/7 * sqrt(6/5))
            b = sqrt(3/7 + 2/7 * sqrt(6/5))
            return [-a, -b, b, a]

    def lagrange(self):
        w = self.weights()
        ξ = self.nodes()
        quad = zeros((self.N, self.N))
        for i in range(self.N):
            for j in range(self.N):
                sum = 0
                for k in range(self.order):
                    sum += w[k]*log(.5*self.δx[j]*ξ[k] + self.ж[j] - self.ж[i])
                quad[i,j] = .5*self.dS[j]*real(sum)
        return quad

    def quad(self):
        if self.method == "Lagrange":
            quad = self.lagrange()
        return quad
