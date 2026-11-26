"""Effect of cure cool-down on Tsai-Wu first-ply failure pressure."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from angle_study import PAIRS, PRESSURE, R_INNER, THICKNESS, ply_stress
from cylinder import CARBON_EPOXY, build_layers, solve_analytical
from failure_study import F66, LINEAR, QUAD

DT = -160.0  # cooling from 180 C to 20 C


def failure_index(state):
    """Tsai-Wu index of (s1, s2, s3, t12); failure at 1."""
    s, t = np.array(state[:3]), state[3]
    return LINEAR @ s + s @ QUAD @ s + F66 * t**2


def pressure_factor(residual, mech):
    """f at which residual + f * mech fails; 0 if residual alone fails."""
    sT, tT = np.array(residual[:3]), residual[3]
    sM, tM = np.array(mech[:3]), mech[3]
    a = sM @ QUAD @ sM + F66 * tM**2
    b = LINEAR @ sM + 2 * sT @ QUAD @ sM + 2 * F66 * tT * tM
    c = failure_index(residual) - 1
    if c >= 0:
        return 0.0
    if a == 0:
        return np.inf
    return (-b + np.sqrt(b * b - 4 * a * c)) / (2 * a)


def solve_states(angle, pressure, dT):
    """Radii and ply-axis stress states through the wall."""
    layup = [(angle, THICKNESS), (-angle, THICKNESS)] * PAIRS
    layers = build_layers(R_INNER, layup, CARBON_EPOXY)
    sol = solve_analytical(layers, pressure, dT=dT)
    radii, states = [], []
    for i, (layer, (ply, _)) in enumerate(zip(layers, layup)):
        for r in np.linspace(layer.r_in, layer.r_out, 21):
            stress = sol.stress(i, r)
            s1, s2, t12 = ply_stress(stress, ply)
            radii.append(r)
            states.append((s1, s2, stress[2], t12))
    return radii, states


def study(angle):
    _, mech = solve_states(angle, PRESSURE, 0.0)
    _, res = solve_states(angle, 0.0, DT)
    zero = (0, 0, 0, 0)
    cold = min(pressure_factor(zero, m) for m in mech) * PRESSURE
    cured = min(pressure_factor(t, m) for t, m in zip(res, mech)) * PRESSURE
    return cold, cured, max(failure_index(t) for t in res)


def main():
    angles = np.arange(0, 91)
    without, with_t, index = np.array([study(a) for a in angles]).T

    print("angle  no thermal[MPa]  with thermal[MPa]  cooling index")
    for a in (0, 15, 30, 45, 55, 60, 75, 90):
        print(f"{a:>5}{without[a]:>17.2f}{with_t[a]:>19.2f}{index[a]:>15.3f}")
    for name, p in (("No thermal", without), ("With thermal", with_t)):
        print(f"{name}: highest {p.max():.2f} MPa at {angles[p.argmax()]} deg")
    print("Angles failing on cooling alone:", angles[index >= 1].tolist())

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4.5))
    ax1.plot(angles, without, label="no thermal stress")
    ax1.plot(angles, with_t, label=f"cooled, dT = {DT:.0f} K")
    ax1.set(xlabel="winding angle [deg]",
            ylabel="first-ply failure pressure [MPa]")

    radii, res = solve_states(55, 0.0, DT)
    res = np.array(res)
    for col, label in ((1, "transverse"), (2, "radial"), (3, "shear")):
        ax2.plot(radii, res[:, col], label=label)
    ax2.set(xlabel="radius [mm]", ylabel="residual stress [MPa]",
            title="Cooling only, 55 deg")
    for ax in (ax1, ax2):
        ax.grid(alpha=0.3)
        ax.legend()
    fig.tight_layout()
    fig.savefig("thermal_study.png", dpi=150)


if __name__ == "__main__":
    main()
