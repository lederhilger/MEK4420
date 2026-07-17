import numpy as np
from numpy.polynomial.legendre import leggauss

class Quadrature:
    def __init__(self, N, δx, ж, method: str, order: int):
        self.N = N
        self.δx = np.asarray(δx[0], dtype=float) + 1j*np.asarray(δx[1], dtype=float)
        self.ж = np.asarray(ж[0], dtype=float) + 1j*np.asarray(ж[1], dtype=float)
        self.dS = np.abs(self.δx)
        self.method = method
        if order not in (2, 4):
            raise ValueError("Order 2 or 4")
        self.order = order
        self._refined_order = max(16, 4*self.order)

    def weights(self, order = None):
        _, weights = self.rule(order)
        return weights

    def nodes(self, order = None):
        nodes, _ = self.rule(order)
        return nodes

    def rule(self, order = None):
        quadrature_order = self.order if order is None else order
        nodes, weights = leggauss(quadrature_order)
        return nodes, weights

    def panel_distance(self):
        start = self.ж - .5*self.δx
        segment = self.δx
        segment_norm = np.maximum(self.dS**2, np.finfo(float).eps)
        offset = self.ж[:, None] - start[None, :]
        τ = np.real(offset*np.conjugate(segment)[None, :]) / segment_norm[None, :]
        τ = np.clip(τ, 0., 1.)
        projection = start[None, :] + τ*segment[None, :]
        return np.abs(self.ж[:, None] - projection)

    def integrate(self, order):
        ξ, w = self.rule(order)
        evaluation_points = self.ж[None, :] + .5*self.δx[None, :]*ξ[:, None]
        kernel = np.log(np.abs(evaluation_points[:, None, :] - self.ж[None, :, None]))
        return .5*self.dS[None, :]*np.tensordot(w, kernel, axes = (0, 0))

    def diagonal(self):
        return self.dS*(np.log(.5*self.dS) - 1.)

    def lagrange(self):
        quad = self.integrate(self.order)
        np.fill_diagonal(quad, self.diagonal())

        near_mask = self.panel_distance() <= (.5*self.dS)[None, :]
        np.fill_diagonal(near_mask, False)
        if np.any(near_mask):
            refined = self.integrate(self._refined_order)
            quad[near_mask] = refined[near_mask]
        return quad

    def quad(self):
        if self.method == "Lagrange":
            return self.lagrange()
        raise ValueError(f"Unknown quadrature method: {self.method}")
