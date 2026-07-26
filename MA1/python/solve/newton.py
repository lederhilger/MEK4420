from math import isfinite

class Newton:
    def __init__(self, ε=1e-14, maxit=50):
        if ε <= 0: raise ValueError("Positive ε required.")
        if maxit <= 0: raise ValueError("Positive max iterations required.")
        self.ε = ε
        self.maxit = maxit

    def solve(self, function, x_0, left, f_left):
        if not left < x_0: raise ValueError("Domain not oriented.")
        if not isfinite(f_left): raise ValueError("f(left) must be finite.")
        if abs(f_left) <= self.ε: return left

        right = x = x_0
        f, derivative = function(x)

        if not isfinite(f): raise ValueError("f(x_0) must be finite.")
        if abs(f) <= self.ε: return x
        if not f_left < 0 < f: raise ValueError("Root is not bracketed.")

        for _ in range(self.maxit):
            if f < 0: left, f_left = x, f
            else: right = x

            scale = max(1, abs(left), abs(right))
            midpoint = .5*(left + right)
            if right - left <= self.ε*scale: return midpoint

            if derivative != 0 and isfinite(derivative):
                candidate = x - f/derivative
            else: candidate = midpoint
            if not isfinite(candidate) or not left < candidate < right:
                candidate = midpoint

            x = candidate
            f, derivative = function(x)
            if not isfinite(f): raise RuntimeError("f(x) is infinite.")
            if abs(f) <= self.ε: return x
        
        raise RuntimeError("Did not converge.")
