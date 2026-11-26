# Composite thick-walled cylinder

Layered, closed-end composite cylinder under internal/external pressure.
An analytical solution is compared with a 1D finite element model.

## Model
- Axisymmetric, generalised plane strain: u(r), axial strain e0, twist k
- Perfectly bonded, linear elastic, cylindrically orthotropic layers
- Ply angle measured from the cylinder axis (90 deg = hoop)
- Away from end effects; no thermal loads or failure prediction

## Analytical solution
In each layer, radial equilibrium gives

    u'' + u'/r - L^2 u/r^2 = (C12 - C13) e0 / (C33 r) + (C26 - 2 C36) k / C33
    L^2 = C22 / C33        (1 = axial, 2 = hoop, 3 = radial, 6 = shear)

The closed-form solution is `c1 r^L + c2 r^-L` plus particular terms
(logarithmic when L = 1 or 2). The constants, e0 and k follow from surface
pressures, interface continuity, axial force and zero torque. Resultant
integrals are evaluated by Gauss quadrature.

## FEM
Three-node quadratic elements in r, with e0 and k as global unknowns.

## Run
    pip install -r requirements.txt
    pytest
    python main.py

## Results
![Stresses](stresses.png)

![Convergence](convergence.png)
## Verification
All 4 tests pass. With 4 elements per layer:

| Quantity | Analytical | FEM |
|---|---|---|
| Inner radial displacement (mm) | 0.100329 | 0.100329 |
| Outer radial displacement (mm) | 0.040036 | 0.040036 |
| Axial strain | 1.09487e-3 | 1.09487e-3 |
| Twist (rad/mm) | -1.7800e-7 | -1.7800e-7 |

Maximum relative error falls from about 2.5e-4 (1 element per layer)
to about 7e-8 (16 elements per layer), approaching third-order convergence.

Case: inner radius 50 mm, layup [55/-55/90/55/-55] from inside out,
5 mm per layer, 50 MPa internal pressure, closed ends.
