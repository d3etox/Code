[🇫🇷 Version française](README.md)

# Finite Volume Method

> Finite-volume discretization: from a simple 1D elliptic PDE up to a stratigraphic model of sedimentary basin formation.

**Context**: Finite Volumes course labs (M2), based on code by R. Masson.

## What it does
- **TP1 / TP2**: solving the 1D elliptic equation `-u'' = f` with finite volumes, Dirichlet and Neumann conditions; convergence study.
- **Stratigraphic model** (`*_strat_litho.py`): modelling diffusive sediment transport over large space and time scales,
  `∂_t h + div(grad(ψ(b))) = 0`, with **TPFA** fluxes on unstructured meshes. Simulation of basin formation (single lithology).

## Stack
Python (NumPy, Matplotlib).

## Files
- `TP1.py`, `TP2.py`, `TP2_offi.py` — 1D elliptic equation
- `1_strat_litho.py`, `2_strat_litho.py` — stratigraphic model

## Run
```bash
python TP1.py
python 1_strat_litho.py
```
