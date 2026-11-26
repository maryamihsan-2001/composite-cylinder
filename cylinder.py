"""Layered composite thick cylinder: pressure and uniform temperature change.

Axisymmetric, generalised plane strain with twist, closed ends.
Vector order: [z (axial), theta (hoop), r (radial), theta-z (shear)].
Units: mm, N, MPa, K.
"""
from dataclasses import dataclass, field
from typing import Callable

import numpy as np

TOL = 1e-8  # tolerance for detecting the resonant cases lambda = 1 or 2


@dataclass(frozen=True)
class Material:
    E1: float
    E2: float
    E3: float
    G12: float
    nu12: float
    nu13: float
    nu23: float
    alpha1: float = 0.0  # thermal expansion coefficients [1/K]
    alpha2: float = 0.0
    alpha3: float = 0.0


CARBON_EPOXY = Material(E1=135e3, E2=10e3, E3=10e3, G12=5e3,
                        nu12=0.27, nu13=0.27, nu23=0.45,
                        alpha1=-0.5e-6, alpha2=27e-6, alpha3=27e-6)


@dataclass(frozen=True)
class Layer:
    r_in: float
    r_out: float
    C: np.ndarray  # 4x4 stiffness in cylinder axes
    alpha: np.ndarray = field(default_factory=lambda: np.zeros(4))


@dataclass
class Solution:
    layers: list
    eps0: float    # axial strain
    kappa: float   # twist per unit length (rad/mm)
    kinematics: Callable  # (layer index, r) -> (u, du/dr)
    dT: float = 0.0

    def u(self, i, r):
        return self.kinematics(i, r)[0]

    def stress(self, i, r):
        u, du = self.kinematics(i, r)
        strain = np.array([self.eps0, u / r, du, self.kappa * r])
        layer = self.layers[i]
        return layer.C @ (strain - self.dT * layer.alpha)


def ply_stiffness(mat, angle_deg):
    """Stiffness of a ply wound at angle_deg from the cylinder axis."""
    S = np.array([
        [1 / mat.E1, -mat.nu12 / mat.E1, -mat.nu13 / mat.E1],
        [-mat.nu12 / mat.E1, 1 / mat.E2, -mat.nu23 / mat.E2],
        [-mat.nu13 / mat.E1, -mat.nu23 / mat.E2, 1 / mat.E3],
    ])
    C = np.zeros((4, 4))
    C[:3, :3] = np.linalg.inv(S)
    C[3, 3] = mat.G12

    c, s = np.cos(np.radians(angle_deg)), np.sin(np.radians(angle_deg))
    # strain transformation: cylinder axes -> material axes
    T = np.array([
        [c * c, s * s, 0, c * s],
        [s * s, c * c, 0, -c * s],
        [0, 0, 1, 0],
        [-2 * c * s, 2 * c * s, 0, c * c - s * s],
    ])
    return T.T @ C @ T


def ply_expansion(mat, angle_deg):
    """Thermal expansion coefficients in cylinder axes."""
    c, s = np.cos(np.radians(angle_deg)), np.sin(np.radians(angle_deg))
    return np.array([
        c * c * mat.alpha1 + s * s * mat.alpha2,
        s * s * mat.alpha1 + c * c * mat.alpha2,
        mat.alpha3,
        2 * c * s * (mat.alpha1 - mat.alpha2),
    ])


def build_layers(r_inner, layup, mat):
    """layup: list of (angle_deg, thickness), ordered inside to outside."""
    layers, r = [], r_inner
    for angle, thickness in layup:
        layers.append(Layer(r, r + thickness, ply_stiffness(mat, angle),
                            ply_expansion(mat, angle)))
        r += thickness
    return layers


def axial_force(layers, p_in, p_out):
    """Closed-end axial resultant."""
    return np.pi * (p_in * layers[0].r_in**2 - p_out * layers[-1].r_out**2)


def _gauss(a, b, n):
    x, w = np.polynomial.legendre.leggauss(n)
    return 0.5 * (a + b) + 0.5 * (b - a) * x, 0.5 * (b - a) * w


def _solve_scaled(A, b):
    """Dense solve with column and row equilibration."""
    col = 1 / np.abs(A).max(axis=0)
    Ac = A * col
    row = 1 / np.abs(Ac).max(axis=1)
    return col * np.linalg.solve(Ac * row[:, None], b * row)


# ---------------------------------------------------------------- analytical
def _particular(k, lam, r, r_ref):
    """Particular solution of u'' + u'/r - lam^2 u/r^2 = k/r."""
    if abs(lam - 1) < TOL:
        log = np.log(r / r_ref)
        return k / 2 * r * log, k / 2 * (log + 1)
    a = k / (1 - lam**2)
    return a * r, a


def _basis(layer, r):
    """u and du/dr as linear functions of [c1, c2, eps0, kappa]."""
    C = layer.C
    lam = np.sqrt(C[1, 1] / C[2, 2])
    f1 = (r / layer.r_out) ** lam
    f2 = (r / layer.r_in) ** (-lam)
    log = np.log(r / layer.r_in)

    g, dg = _particular((C[0, 1] - C[0, 2]) / C[2, 2], lam, r, layer.r_in)

    k = (C[1, 3] - 2 * C[2, 3]) / C[2, 2]        # twist forcing
    if abs(lam - 2) < TOL:
        h, dh = k / 4 * r**2 * log, k / 4 * r * (2 * log + 1)
    else:
        b = k / (4 - lam**2)
        h, dh = b * r**2, 2 * b * r

    u = np.array([f1, f2, g, h])
    du = np.array([lam * f1 / r, -lam * f2 / r, dg, dh])
    return u, du


def _thermal_basis(layer, r, dT):
    """Particular u and du/dr caused by a uniform temperature change."""
    C = layer.C
    beta = C @ layer.alpha * dT
    lam = np.sqrt(C[1, 1] / C[2, 2])
    k = (beta[2] - beta[1]) / C[2, 2]
    return _particular(k, lam, r, layer.r_in)


def solve_analytical(layers, p_in, p_out=0.0, dT=0.0):
    n = len(layers)
    A = np.zeros((2 * n + 2, 2 * n + 2))
    b = np.zeros(2 * n + 2)

    def dofs(i):
        return [2 * i, 2 * i + 1, 2 * n, 2 * n + 1]

    def evaluate(i, r):
        layer = layers[i]
        u, du = _basis(layer, r)
        u0, du0 = _thermal_basis(layer, r, dT)
        strain = np.array([[0, 0, 1, 0], u / r, du, [0, 0, 0, r]], dtype=float)
        S = layer.C @ strain
        S0 = layer.C @ (np.array([0, u0 / r, du0, 0]) - dT * layer.alpha)
        return u, S, u0, S0

    # radial stress on the inner and outer surfaces
    surfaces = [(0, layers[0].r_in, p_in), (n - 1, layers[-1].r_out, p_out)]
    for row, (i, r, p) in enumerate(surfaces):
        _, S, _, S0 = evaluate(i, r)
        A[row, dofs(i)] = S[2]
        b[row] = -p - S0[2]

    # continuity of u and radial stress at each interface
    eq = 2
    for i in range(n - 1):
        r = layers[i].r_out
        for j, sign in ((i, 1), (i + 1, -1)):
            u, S, u0, S0 = evaluate(j, r)
            A[eq, dofs(j)] += sign * u
            A[eq + 1, dofs(j)] += sign * S[2]
            b[eq] -= sign * u0
            b[eq + 1] -= sign * S0[2]
        eq += 2

    # axial force and zero torque
    b[eq] = axial_force(layers, p_in, p_out)
    for j, layer in enumerate(layers):
        for r, w in zip(*_gauss(layer.r_in, layer.r_out, 20)):
            _, S, _, S0 = evaluate(j, r)
            A[eq, dofs(j)] += 2 * np.pi * r * w * S[0]
            A[eq + 1, dofs(j)] += 2 * np.pi * r**2 * w * S[3]
            b[eq] -= 2 * np.pi * r * w * S0[0]
            b[eq + 1] -= 2 * np.pi * r**2 * w * S0[3]

    x = _solve_scaled(A, b)

    def kinematics(i, r):
        u, du = _basis(layers[i], r)
        u0, du0 = _thermal_basis(layers[i], r, dT)
        return u @ x[dofs(i)] + u0, du @ x[dofs(i)] + du0

    return Solution(layers, x[-2], x[-1], kinematics, dT)


# ---------------------------------------------------------------- FEM
def _shape(r1, r2, r):
    """Quadratic 1D shape functions and their r-derivatives."""
    s = (2 * r - r1 - r2) / (r2 - r1)
    N = np.array([s * (s - 1) / 2, 1 - s * s, s * (s + 1) / 2])
    dN = np.array([s - 0.5, -2 * s, s + 0.5]) * 2 / (r2 - r1)
    return N, dN


def solve_fem(layers, p_in, p_out=0.0, elems_per_layer=8, dT=0.0):
    m = elems_per_layer
    n_u = 2 * len(layers) * m + 1
    K = np.zeros((n_u + 2, n_u + 2))   # nodal u, then eps0 and kappa
    F = np.zeros(n_u + 2)

    for i, layer in enumerate(layers):
        edges = np.linspace(layer.r_in, layer.r_out, m + 1)
        for k in range(m):
            e = i * m + k
            dofs = [2 * e, 2 * e + 1, 2 * e + 2, n_u, n_u + 1]
            for r, w in zip(*_gauss(edges[k], edges[k + 1], 4)):
                N, dN = _shape(edges[k], edges[k + 1], r)
                B = np.zeros((4, 5))
                B[0, 3] = 1
                B[1, :3] = N / r
                B[2, :3] = dN
                B[3, 4] = r
                dV = 2 * np.pi * r * w
                K[np.ix_(dofs, dofs)] += dV * B.T @ layer.C @ B
                F[dofs] += dV * B.T @ (layer.C @ layer.alpha * dT)

    F[0] += 2 * np.pi * layers[0].r_in * p_in
    F[n_u - 1] -= 2 * np.pi * layers[-1].r_out * p_out
    F[n_u] += axial_force(layers, p_in, p_out)

    x = _solve_scaled(K, F)

    def kinematics(i, r):
        edges = np.linspace(layers[i].r_in, layers[i].r_out, m + 1)
        k = int(np.clip(np.searchsorted(edges, r, side="right") - 1, 0, m - 1))
        e = i * m + k
        N, dN = _shape(edges[k], edges[k + 1], r)
        ue = x[2 * e:2 * e + 3]
        return N @ ue, dN @ ue

    return Solution(layers, x[-2], x[-1], kinematics, dT)


# ---------------------------------------------------------------- comparison
def gauss_points(layer, m):
    """Two Gauss points per element, where FEM stresses are most accurate."""
    edges = np.linspace(layer.r_in, layer.r_out, m + 1)
    return np.concatenate([_gauss(a, b, 2)[0]
                           for a, b in zip(edges[:-1], edges[1:])])


def relative_errors(exact, approx, m):
    """Maximum error divided by the maximum exact magnitude."""
    points = [(i, r) for i, layer in enumerate(exact.layers)
              for r in gauss_points(layer, m)]
    errors = {}
    for name in ("u", "stress"):
        e = np.array([getattr(exact, name)(i, r) for i, r in points])
        a = np.array([getattr(approx, name)(i, r) for i, r in points])
        errors[name] = np.abs(a - e).max() / np.abs(e).max()
    return errors
