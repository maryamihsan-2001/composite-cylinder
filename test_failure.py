import pytest

from failure_study import S12, XC, XT, YC, YT, tsai_wu_ratio


@pytest.mark.parametrize("state", [
    (XT, 0, 0, 0), (-XC, 0, 0, 0),
    (0, YT, 0, 0), (0, -YC, 0, 0),
    (0, 0, YT, 0), (0, 0, -YC, 0),
    (0, 0, 0, S12),
])
def test_uniaxial_strength_gives_unit_factor(state):
    assert tsai_wu_ratio(*state) == pytest.approx(1.0)
