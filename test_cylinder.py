import numpy as np
import pytest

from cylinder import (CARBON_EPOXY, Material, build_layers, relative_errors,
                      solve_analytical, solve_fem)

LAYUP = [(55, 5.0), (-55, 5.0), (90, 5.0), (0, 5.0), (30, 5.0)]


@pytest.mark.parametrize("solver, tol", [(solve_analytical, 1e-6),
                                         (solve_fem, 1e-2)])
def test_isotropic_matches_lame(solver, tol):
    E, nu, a, b, p = 200e3, 0.3, 50.0, 100.0, 100.0
    iso = Material(E, E, E, E / (2 * (1 + nu)), nu, nu, nu)
    layers = build_layers(a, [(ang, 12.5) for ang in (0, 30, 60, 90)], iso)
    sol = solver(layers, p)

    A = p * a**2 / (b**2 - a**2)
    B = p * a**2 * b**2 / (b**2 - a**2)
    for i, layer in enumerate(layers):
        r = 0.5 * (layer.r_in + layer.r_out) + 0.1
        sr, st = A - B / r**2, A + B / r**2
        u = r * (st - nu * (sr + A)) / E
        assert sol.u(i, r) == pytest.approx(u, rel=1e-4)
        # element-interior stresses; the default FEM mesh has 8 per layer
        np.testing.assert_allclose(sol.stress(i, r), [A, st, sr, 0],
                                   atol=tol if solver is solve_analytical else 0.1)


def test_boundary_and_interface_conditions():
    layers = build_layers(50.0, LAYUP, CARBON_EPOXY)
    sol = solve_analytical(layers, 50.0, 5.0)
    assert sol.stress(0, layers[0].r_in)[2] == pytest.approx(-50.0)
    assert sol.stress(4, layers[-1].r_out)[2] == pytest.approx(-5.0)
    for i in range(len(layers) - 1):
        r = layers[i].r_out
        assert sol.u(i, r) == pytest.approx(sol.u(i + 1, r))
        assert sol.stress(i, r)[2] == pytest.approx(sol.stress(i + 1, r)[2])


def test_fem_converges_to_analytical():
    layers = build_layers(50.0, LAYUP, CARBON_EPOXY)
    exact = solve_analytical(layers, 50.0)
    errs = [relative_errors(exact, solve_fem(layers, 50.0, elems_per_layer=m), m)
            for m in (2, 4, 8)]
    for key in ("u", "stress"):
        assert errs[0][key] > errs[1][key] > errs[2][key]
        assert errs[2][key] < 1e-4
