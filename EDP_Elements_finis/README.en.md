[🇫🇷 Version française](README.md)

# Finite Elements P1/P2/P3 — 2D Convection-Reaction-Diffusion

> Implementation and validation of Lagrange (Pk) finite elements for an anisotropic 2D convection-reaction-diffusion problem.

**Context**: Take-home assignment, Finite Elements course (Polytech Nice-Sophia, M2).

## What it does
Numerical solution of the elliptic problem
`-div(A∇u) + b·∇u + c·u = f` on a unit square, with the known exact solution
`u(x,y) = sin(πx)·cos(2πy)` used to **validate convergence** of P1, P2 and P3 elements.

- Assembly of stiffness/mass matrices on a 2D mesh
- Computation of basis functions and their derivatives (P2, P3)
- **Convergence order** study in the L² norm depending on element degree
- Test cases: isotropic/anisotropic diffusion, convection, reaction

## Stack
Python (NumPy, SciPy `sparse`), Scilab (`.sce`). LaTeX report in `projet/`.

## Main files
- `DM_FABRE_EM.py` — main code
- `Laplacian_Pk_2D_2025_V_AF.py` — 2D Pk solver
- `compute_P2_derivatives.py`, `compute_P3_derivatives.py` — high-order bases
- `test_*.py` — validation and convergence tests
- `png/` — result figures · `DM_FABRE_EM.pdf` — report

## Run
```bash
python DM_FABRE_EM.py
```
