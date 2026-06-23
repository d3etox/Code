[🇫🇷 Version française](README.md)

# Numerical schemes for the advection equation

> Comparison of finite-difference schemes for the linear transport equation `∂_t u + a ∂_x u = 0`.

**Context**: numerical methods project (mathematics for physics).

## What it does
Implementation and comparison of 6 schemes on a discontinuous initial condition:
- Upwind (backward / forward differencing)
- Centered scheme
- **Lax-Friedrichs**
- **Lax-Wendroff**
- Second-order upwind

Analysis of **numerical diffusion**, **oscillations** (Gibbs phenomenon) and **convergence order** (L² error for several `dx` resolutions).

## Stack
Python (NumPy, Matplotlib).

## Files
- `maths_physique.py` — schemes and visualization
- `maths_physique2.py` — convergence study (L² error vs resolution)

## Run
```bash
python maths_physique.py
```
