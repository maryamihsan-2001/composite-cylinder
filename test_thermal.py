import numpy as np
import pytest

from cylinder import (CARBON_EPOXY, Material, build_layers, ply_expansion,
                      relative_errors, solve_analytical, solve_fem)
from failure_study import tsai_wu_ratio
from thermal_study import pressure_factor

LAYUP = [(55, 5.0), (-55, 5.0), (90, 5.0), (0, 5.0), (30, 5.0)]
DT = -160.0


@pytest.mark.parametrize("solver", [solve_analytical, solve_fem])
def test_isotropic_free_expansion(solver):
    E, nu, a = 200e3, 0.3, 1.2e-5
    iso = Material(E, E, E, E / (2 * (1 + nu)), nu, nu, nu, a, a, a)
    layers = build_layers(50.0, [(ang, 10.0) for ang in (0, 30, 90)], iso)
    sol = solver(layers, 0.0, dT=DT)
    assert sol.eps0 == pytest.approx(a * DT)
    for i, layer in enumerate(layers):
        r = layer.r_in + 3.1
        assert sol.u(i, r) == pytest.approx(a * DT * r)
        np.testing.assert_allclose(sol.stress(i, r), 0, atol=1e-5)


def test_fem_matches_analytical_with_temperature():
    layers = build_layers(50.0, LAYUP, CARBON_EPOXY)
    exact = solve_analytical(layers, 50.0, dT=DT)
    fem = solve_fem(layers, 50.0, elems_per_layer=8, dT=DT)
    errors = relative_errors(exact, fem, 8)
    assert errors["u"] < 1e-4
    assert errors["stress"] < 1e-4


def test_surfaces_and_superposition():
    layers = build_layers(50.0, LAYUP, CARBON_EPOXY)
    thermal = solve_analytical(layers, 0.0, dT=DT)
    pressure = solve_analytical(layers, 50.0)
    both = solve_analytical(layers, 50.0, dT=DT)
    assert thermal.stress(0, layers[0].r_in)[2] == pytest.approx(0, abs=1e-6)
    assert thermal.stress(4, layers[-1].r_out)[2] == pytest.approx(0, abs=1e-6)
    for i, layer in enumerate(layers):
        r = layer.r_in + 2.0
        np.testing.assert_allclose(
            both.stress(i, r), thermal.stress(i, r) + pressure.stress(i, r),
            atol=1e-6)


def test_expansion_rotation():
    m = CARBON_EPOXY
    c, s = np.cos(np.radians(30)), np.sin(np.radians(30))
    R = np.array([[c, s], [-s, c]])
    e = R.T @ np.diag([m.alpha1, m.alpha2]) @ R
    np.testing.assert_allclose(
        ply_expansion(m, 30), [e[0, 0], e[1, 1], m.alpha3, 2 * e[0, 1]])


def test_zero_residual_matches_original_criterion():
    state = (300.0, 20.0, -10.0, 15.0)
    assert pressure_factor((0, 0, 0, 0), state) == pytest.approx(
        tsai_wu_ratio(*state))
