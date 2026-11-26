"""Winding angle study for a [+a/-a] x 5 composite cylinder."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from cylinder import (CARBON_EPOXY, build_layers, relative_errors,
                      solve_analytical, solve_fem)

R_INNER = 50.0    # mm
THICKNESS = 2.5   # mm per layer
PAIRS = 5         # 10 layers, 25 mm wall
PRESSURE = 50.0   # MPa


def solve_angle(angle):
    layup = [(angle, THICKNESS), (-angle, THICKNESS)] * PAIRS
    layers = build_layers(R_INNER, layup, CARBON_EPOXY)
    return layup, layers, solve_analytical(layers, PRESSURE)


def ply_stress(stress, angle):
    """Cylinder-axis stresses -> fibre, transverse and in-plane shear."""
    sz, st, _, tau = stress
    c, s = np.cos(np.radians(angle)), np.sin(np.radians(angle))
    return (c*c*sz + s*s*st + 2*c*s*tau,
            s*s*sz + c*c*st - 2*c*s*tau,
            c*s*(st - sz) + (c*c - s*s)*tau)


def summarise(angle):
    layup, layers, sol = solve_angle(angle)
    stress, local = [], []
    for i, (layer, (ply, _)) in enumerate(zip(layers, layup)):
        for r in np.linspace(layer.r_in, layer.r_out, 21):
            s = sol.stress(i, r)
            stress.append(s)
            local.append(ply_stress(s, ply))
    stress, local = np.abs(stress), np.abs(local)
    return {
        "angle": angle,
        "u_in": sol.u(0, layers[0].r_in),
        "u_out": sol.u(len(layers) - 1, layers[-1].r_out),
        "eps0": sol.eps0,
        "twist": sol.kappa,
        "hoop": stress[:, 1].max(),
        "axial": stress[:, 0].max(),
        "fibre": local[:, 0].max(),
        "transverse": local[:, 1].max(),
        "shear": local[:, 2].max(),
    }


def main():
    angles = np.arange(0, 91)
    data = [summarise(a) for a in angles]
    col = {key: np.array([d[key] for d in data]) for key in data[0]}

    print("".join(f"{key:>12}" for key in col))
    for a in (0, 15, 30, 45, 55, 60, 75, 90):
        print("".join(f"{col[key][a]:>12.4g}" for key in col))

    print("\nAngle giving the smallest magnitude:")
    for key in ("u_in", "u_out", "eps0", "hoop", "fibre",
                "transverse", "shear"):
        print(f"  {key:<12}{angles[np.argmin(np.abs(col[key]))]} deg")

    print("\nFEM check (4 elements per layer):")
    for a in (0, 30, 55, 90):
        _, layers, exact = solve_angle(a)
        fem = solve_fem(layers, PRESSURE, elems_per_layer=4)
        err = relative_errors(exact, fem, 4)
        print(f"  {a:>2} deg: u {err['u']:.2e}, stress {err['stress']:.2e}")

    np.savetxt("angle_study.csv", np.column_stack(list(col.values())),
               delimiter=",", header=",".join(col), comments="")

    fig, ax = plt.subplots(2, 3, figsize=(15, 8))

    def panel(a, series, title):
        for key, label in series:
            a.plot(angles, col[key], label=label)
        a.set(title=title, xlabel="winding angle [deg]")
        a.grid(alpha=0.3)
        a.legend()

    panel(ax[0, 0], [("u_in", "inner"), ("u_out", "outer")],
          "Radial displacement [mm]")
    panel(ax[0, 1], [("eps0", "axial strain")], "Axial strain [-]")
    panel(ax[0, 2], [("twist", "twist")], "Twist [rad/mm]")
    panel(ax[1, 0], [("hoop", "hoop"), ("axial", "axial")],
          "Peak |stress|, cylinder axes [MPa]")
    panel(ax[1, 1], [("fibre", "fibre"), ("transverse", "transverse"),
                     ("shear", "shear")], "Peak |stress|, ply axes [MPa]")

    for angle in (0, 30, 55, 75, 90):
        _, layers, sol = solve_angle(angle)
        points = [(i, r) for i, layer in enumerate(layers)
                  for r in np.linspace(layer.r_in, layer.r_out, 21)]
        ax[1, 2].plot([r for _, r in points],
                      [sol.stress(i, r)[1] for i, r in points],
                      label=f"{angle} deg")
    ax[1, 2].set(title="Hoop stress through wall [MPa]", xlabel="radius [mm]")
    ax[1, 2].grid(alpha=0.3)
    ax[1, 2].legend()

    fig.tight_layout()
    fig.savefig("angle_study.png", dpi=150)


if __name__ == "__main__":
    main()
