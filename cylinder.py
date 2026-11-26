"""Layered composite thick cylinder under internal/external pressure.

Axisymmetric, generalised plane strain with twist, closed ends.
Vector order: [z (axial), theta (hoop), r (radial), theta-z (shear)].
Units: mm, N, MPa.
"""
from dataclasses import dataclass
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


CARBON_EPOXY = Material(E1=135e3, E2=10e3, E3=10e3, G12=5e3,
                        nu12=0.27, nu13=0.27, nu23=0.45)


@dataclass(frozen=True)
class Layer:
    r_in: float
    r_out: float
    C: np.ndarray  # 4x4 stiffness in cylinder axes


@dataclass
class Solution:
    layers: list
    eps0: float    # axial strain
    kappa: float   # twist per unit length (rad/mm)
    kinematics: Callable  # (layer index, r) -> (u, du/dr)

    def u(self, i, r):
        return self.kinematics(i, r)[0]

    def stress(self, i, r):
        u, du = self.kinematics(i, r)
        strain = np.array([self.eps0, u / r, du, self.kappa * r])
        return self.layers[i].C @ strain


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


def build_layers(r_inner, layup, mat):
    """layup: list of (angle_deg, thickness), ordered inside to outside."""
    layers, r = [], r_inner
    for angle, thickness in layup:
        layers.append(Layer(r, r + thickness, ply_stiffness(mat, angle)))
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
def _basis(layer, r):
    """u and du/dr as linear functions of [c1, c2, eps0, kappa]."""
    C = layer.C
    lam = np.sqrt(C[1, 1] / C[2, 2])
    f1 = (r / layer.r_out) ** lam
    f2 = (r / layer.r_in) ** (-lam)
    log = np.log(r / layer.r_in)

    k = (C[0, 1] - C[0, 2]) / C[2, 2]            # axial strain forcing
    if abs(lam - 1) < TOL:
        g, dg = k / 2 * r * log, k / 2 * (log + 1)
    else:
        a = k / (1 - lam**2)
        g, dg = a * r, a

    k = (C[1, 3] - 2 * C[2, 3]) / C[2, 2]        # twist forcing
    if abs(lam - 2) < TOL:
        h, dh = k / 4 * r**2 * log, k / 4 * r * (2 * log + 1)
    else:
        b = k / (4 - lam**2)
        h, dh = b * r**2, 2 * b * r

    u = np.array([f1, f2, g, h])
    du = np.array([lam * f1 / r, -lam * f2 / r, dg, dh])
    return u, du


def solve_analytical(layers, p_in, p_out=0.0):
    n = len(layers)
    A = np.zeros((2 * n + 2, 2 * n + 2))
    b = np.zeros(2 * n + 2)

    def dofs(i):
        return [2 * i, 2 * i + 1, 2 * n, 2 * n + 1]

    def evaluate(i, r):
        u, du = _basis(layers[i], r)
        strain = np.array([[0, 0, 1, 0], u / r, du, [0, 0, 0, r]], dtype=float)
        return u, layers[i].C @ strain

    # radial stress on the inner and outer surfaces
    A[0, dofs(0)] = evaluate(0, layers[0].r_in)[1][2]
    A[1, dofs(n - 1)] = evaluate(n - 1, layers[-1].r_out)[1][2]
    b[0], b[1] = -p_in, -p_out

    # continuity of u and radial stress at each interface
    eq = 2
    for i in range(n - 1):
        r = layers[i].r_out
        for j, sign in ((i, 1), (i + 1, -1)):
            u, S = evaluate(j, r)
            A[eq, dofs(j)] += sign * u
            A[eq + 1, dofs(j)] += sign * S[2]
        eq += 2

    # axial force and zero torque
    for j, layer in enumerate(layers):
        for r, w in zip(*_gauss(layer.r_in, layer.r_out, 20)):
            S = evaluate(j, r)[1]
            A[eq, dofs(j)] += 2 * np.pi * r * w * S[0]
            A[eq + 1, dofs(j)] += 2 * np.pi * r**2 * w * S[3]
    b[eq] = axial_force(layers, p_in, p_out)

    x = _solve_scaled(A, b)

    def kinematics(i, r):
        u, du = _basis(layers[i], r)
        return u @ x[dofs(i)], du @ x[dofs(i)]

    return Solution(layers, x[-2], x[-1], kinematics)


# ---------------------------------------------------------------- FEM
def _shape(r1, r2, r):
    """Quadratic 1D shape functions and their r-derivatives."""
    s = (2 * r - r1 - r2) / (r2 - r1)
    N = np.array([s * (s - 1) / 2, 1 - s * s, s * (s + 1) / 2])
    dN = np.array([s - 0.5, -2 * s, s + 0.5]) * 2 / (r2 - r1)
    return N, dN


def solve_fem(layers, p_in, p_out=0.0, elems_per_layer=8):
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
                K[np.ix_(dofs, dofs)] += 2 * np.pi * r * w * B.T @ layer.C @ B

    F[0] = 2 * np.pi * layers[0].r_in * p_in
    F[n_u - 1] = -2 * np.pi * layers[-1].r_out * p_out
    F[n_u] = axial_force(layers, p_in, p_out)

    x = _solve_scaled(K, F)

    def kinematics(i, r):
        edges = np.linspace(layers[i].r_in, layers[i].r_out, m + 1)
        k = int(np.clip(np.searchsorted(edges, r, side="right") - 1, 0, m - 1))
        e = i * m + k
        N, dN = _shape(edges[k], edges[k + 1], r)
        ue = x[2 * e:2 * e + 3]
        return N @ ue, dN @ ue

    return Solution(layers, x[-2], x[-1], kinematics)


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
