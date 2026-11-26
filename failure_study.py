"""First-ply failure pressure: Tsai-Wu and maximum stress criteria."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from angle_study import PRESSURE, ply_stress, solve_angle
from cylinder import CARBON_EPOXY, build_layers, solve_analytical
from main import LAYUP, P_INNER, R_INNER

# Typical T300/5208 carbon/epoxy strengths [MPa], illustrative
XT, XC = 1500.0, 1500.0
YT, YC = 40.0, 246.0
S12 = 68.0

F11, F22 = 1 / (XT * XC), 1 / (YT * YC)
F12 = -0.5 * np.sqrt(F11 * F22)
F23 = -0.5 * F22
LINEAR = np.array([1/XT - 1/XC, 1/YT - 1/YC, 1/YT - 1/YC])
QUAD = np.array([[F11, F12, F12],
                 [F12, F22, F23],
                 [F12, F23, F22]])
F66 = 1 / S12**2


def tsai_wu_ratio(s1, s2, s3, t12):
    """Load factor that brings this stress state to Tsai-Wu failure."""
    s = np.array([s1, s2, s3])
    a = s @ QUAD @ s + F66 * t12**2
    b = LINEAR @ s
    d = b + np.sqrt(b * b + 4 * max(a, 0.0))
    return 2 / d if d > 0 else np.inf


def max_stress_ratio(s1, s2, s3, t12):
    """Load factor and failure mode for the maximum stress criterion."""
    used = {
        "fibre tension": s1 / XT, "fibre compression": -s1 / XC,
        "matrix tension": s2 / YT, "matrix compression": -s2 / YC,
        "radial tension": s3 / YT, "radial compression": -s3 / YC,
        "shear": abs(t12) / S12,
    }
    mode = max(used, key=used.get)
    return 1 / used[mode], mode


def first_ply_failure(layup, layers, sol, pressure):
    """Lowest failure pressure over all layers (linear elastic scaling)."""
    tw, ms = [], []
    for i, (layer, (angle, _)) in enumerate(zip(layers, layup)):
        for r in np.linspace(layer.r_in, layer.r_out, 21):
            stress = sol.stress(i, r)
            s1, s2, t12 = ply_stress(stress, angle)
            state = (s1, s2, stress[2], t12)
            tw.append((tsai_wu_ratio(*state) * pressure, i + 1, r))
            factor, mode = max_stress_ratio(*state)
            ms.append((factor * pressure, i + 1, r, mode))
    return min(tw), min(ms)


def main():
    angles = np.arange(0, 91)
    results = [first_ply_failure(*solve_angle(a), PRESSURE) for a in angles]
    tw = np.array([res[0][0] for res in results])
    ms = np.array([res[1][0] for res in results])

    print("angle  Tsai-Wu[MPa]  layer  MaxStress[MPa]  layer  mode")
    for a in (0, 15, 30, 45, 55, 60, 75, 90):
        t, m = results[a]
        print(f"{a:>5}{t[0]:>14.2f}{t[1]:>7}{m[0]:>16.2f}{m[1]:>7}  {m[3]}")

    for name, p in (("Tsai-Wu", tw), ("Max stress", ms)):
        best = np.argmax(p)
        print(f"{name}: highest {p[best]:.2f} MPa at {angles[best]} deg")

    layers = build_layers(R_INNER, LAYUP, CARBON_EPOXY)
    t, m = first_ply_failure(LAYUP, layers,
                             solve_analytical(layers, P_INNER), P_INNER)
    print(f"\nBaseline layup: Tsai-Wu {t[0]:.2f} MPa (layer {t[1]}), "
          f"max stress {m[0]:.2f} MPa (layer {m[1]}, {m[3]})")

    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.plot(angles, tw, label="Tsai-Wu")
    ax.plot(angles, ms, "--", label="maximum stress")
    ax.set(xlabel="winding angle [deg]",
           ylabel="first-ply failure pressure [MPa]")
    ax.grid(alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig("failure_study.png", dpi=150)


if __name__ == "__main__":
    main()
