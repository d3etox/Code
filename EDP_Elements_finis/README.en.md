[🇫🇷 Version française](README.md)

# Finite Elements P1/P2/P3 — 2D Convection-Reaction-Diffusion & Navier-Stokes

> Implementation and validation of Lagrange (Pk) finite elements for 2D elliptic problems, up to the incompressible Navier-Stokes equations.

**Context**: Finite Elements course and take-home assignment (Polytech Nice-Sophia, M2).

## What it does
- **Take-home — 2D convection-reaction-diffusion**: solving `-div(A∇u) + b·∇u + c·u = f` on a square, with the known exact solution `u(x,y) = sin(πx)·cos(2πy)` to **validate convergence** of P1, P2 and P3 elements (isotropic/anisotropic, convection, reaction).
- **Pk assembly**: P2/P3 basis functions and derivatives, **convergence order** study in the L² norm.
- **1D elements** (`TDTP_01_Pk_1D.py`, `TP_Cours1.ipynb`): P1/P2 Lagrange local matrices.
- **2D Navier-Stokes** (`TP_NavierStokes.ipynb`): solving the incompressible equations, validated on the **Taylor-Green** vortex (unsteady analytical solution).

## Stack
Python (NumPy, SciPy `sparse`), Scilab (`.sce`), Jupyter. LaTeX report in `projet/`.

## Main files
- `DM_FABRE_EM.py` — take-home code (convection-reaction-diffusion)
- `Laplacian_Pk_2D_2025_V_AF.py` — 2D Pk solver
- `compute_P2_derivatives.py`, `compute_P3_derivatives.py`, `P3_analysis.py` — high-order bases
- `TP_NavierStokes.ipynb` — 2D Navier-Stokes (Taylor-Green)
- `TDTP_01_Pk_1D.py`, `TP_Cours1.ipynb`, `TP1.ipynb` — 1D-element labs
- `test_*.py` — validation and convergence tests
- `png/` — result figures · `DM_FABRE_EM.pdf` — report

## Run
```bash
python DM_FABRE_EM.py
```
