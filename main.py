from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from cylinder import (CARBON_EPOXY, build_layers, gauss_points,
                      relative_errors, solve_analytical, solve_fem)

R_INNER = 50.0                      # mm
LAYUP = [(55, 5.0), (-55, 5.0), (90, 5.0), (55, 5.0), (-55, 5.0)]
P_INNER = 50.0                      # MPa
ELEMS = 4                           # elements per layer
OUT = Path("results")


def main():
    OUT.mkdir(exist_ok=True)
    layers = build_layers(R_INNER, LAYUP, CARBON_EPOXY)
    exact = solve_analytical(layers, P_INNER)
    fem = solve_fem(layers, P_INNER, elems_per_layer=ELEMS)

    rows = [
        ("u inner [mm]", lambda s: s.u(0, layers[0].r_in)),
        ("u outer [mm]", lambda s: s.u(len(layers) - 1, layers[-1].r_out)),
        ("axial strain [-]", lambda s: s.eps0),
        ("twist [rad/mm]", lambda s: s.kappa),
    ]
    print(f"{'quantity':<20}{'analytical':>16}{'FEM':>16}{'abs diff':>12}")
    for name, get in rows:
        a, f = get(exact), get(fem)
        print(f"{name:<20}{a:>16.8e}{f:>16.8e}{abs(a - f):>12.2e}")

    # through-thickness stresses
    titles = ["axial", "hoop", "radial", "shear (theta-z)"]
    fig, axes = plt.subplots(2, 2, figsize=(10, 7), sharex=True)
    for c, (ax, title) in enumerate(zip(axes.flat, titles)):
        for i, layer in enumerate(layers):
            r = np.linspace(layer.r_in, layer.r_out, 60)
            rg = gauss_points(layer, ELEMS)
            ax.plot(r, [exact.stress(i, x)[c] for x in r], "k-",
                    label="analytical" if i == 0 else None)
            ax.plot(rg, [fem.stress(i, x)[c] for x in rg], "ro", ms=3,
                    label="FEM" if i == 0 else None)
            ax.axvline(layer.r_out, color="0.8", lw=0.8, zorder=0)
        ax.set_title(f"{title} stress [MPa]")
        ax.grid(alpha=0.3)
    axes[0, 0].legend()
    for ax in axes[1]:
        ax.set_xlabel("radius [mm]")
    fig.tight_layout()
    fig.savefig(OUT / "stresses.png", dpi=150)

    # mesh convergence
    meshes = [1, 2, 4, 8, 16]
    errs = [relative_errors(exact, solve_fem(layers, P_INNER, elems_per_layer=m), m)
            for m in meshes]
    print("\nelements/layer   error u      error stress")
    for m, e in zip(meshes, errs):
        print(f"{m:>14}   {e['u']:.3e}    {e['stress']:.3e}")

    fig, ax = plt.subplots(figsize=(6, 4))
    for key in ("u", "stress"):
        ax.loglog(meshes, [e[key] for e in errs], "o-", label=key)
    ax.set_xlabel("elements per layer")
    ax.set_ylabel("max relative error")
    ax.grid(which="both", alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(OUT / "convergence.png", dpi=150)
    print(f"\nPlots saved in {OUT}/")


if __name__ == "__main__":
    main()
