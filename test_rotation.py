"""Independent checks of the ply stiffness and stress rotations."""
import numpy as np
import pytest

from angle_study import ply_stress
from cylinder import CARBON_EPOXY as M
from cylinder import ply_stiffness

ANGLES = [-75, -55, -30, 15, 45, 60]

# material compliances
S11, S22, S33 = 1 / M.E1, 1 / M.E2, 1 / M.E3
S12, S13, S23 = -M.nu12 / M.E1, -M.nu13 / M.E1, -M.nu23 / M.E2
S66 = 1 / M.G12


def compliance(angle):
    """Compliance in cylinder axes [z, theta, r, z-theta] from the code."""
    return np.linalg.inv(ply_stiffness(M, angle))


def check(actual, expected):
    np.testing.assert_allclose(actual, expected, rtol=1e-10, atol=1e-14)


def test_zero_degrees():
    check(compliance(0), [[S11, S12, S13, 0],
                          [S12, S22, S23, 0],
                          [S13, S23, S33, 0],
                          [0, 0, 0, S66]])


def test_ninety_degrees():
    # axial and hoop properties swap
    check(compliance(90), [[S22, S12, S23, 0],
                           [S12, S11, S13, 0],
                           [S23, S13, S33, 0],
                           [0, 0, 0, S66]])


@pytest.mark.parametrize("angle", ANGLES)
def test_matches_textbook_transformation(angle):
    c, s = np.cos(np.radians(angle)), np.sin(np.radians(angle))
    mix = (2 * S12 + S66) * s**2 * c**2
    p = 2 * S11 - 2 * S12 - S66
    q = 2 * S22 - 2 * S12 - S66

    b11 = S11 * c**4 + mix + S22 * s**4
    b22 = S11 * s**4 + mix + S22 * c**4
    b12 = (S11 + S22 - S66) * s**2 * c**2 + S12 * (s**4 + c**4)
    b16 = p * s * c**3 - q * s**3 * c
    b26 = p * s**3 * c - q * s * c**3
    b66 = (2 * (2 * S11 + 2 * S22 - 4 * S12 - S66) * s**2 * c**2
           + S66 * (s**4 + c**4))
    b13 = S13 * c**2 + S23 * s**2
    b23 = S13 * s**2 + S23 * c**2
    b36 = 2 * (S13 - S23) * s * c

    check(compliance(angle), [[b11, b12, b13, b16],
                              [b12, b22, b23, b26],
                              [b13, b23, S33, b36],
                              [b16, b26, b36, b66]])


@pytest.mark.parametrize("angle", ANGLES)
def test_ply_stress_matches_tensor_rotation(angle):
    sz, st, sr, tau = 40.0, 130.0, -20.0, 25.0
    c, s = np.cos(np.radians(angle)), np.sin(np.radians(angle))
    R = np.array([[c, s], [-s, c]])
    local = R @ np.array([[sz, tau], [tau, st]]) @ R.T
    check(ply_stress([sz, st, sr, tau], angle),
          [local[0, 0], local[1, 1], local[0, 1]])
