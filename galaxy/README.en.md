[🇫🇷 Version française](README.md)

# N-body galaxy simulation

> Gravitational simulation of the interaction between hundreds of bodies (stars/planets) forming a galaxy.

**Context**: personal project (extension of my preparatory-class galaxy-simulation project).

## What it does
N-body simulation of a rotating galactic disk:
- Gravitational force computation with **softening** to avoid singularities
- Two methods: naive **O(N²)** computation and a **quadrant approximation** (Barnes-Hut spirit) to speed up large N
- **Symplectic** (Euler) integration conserving energy
- Total-energy tracking, trajectory export (CSV) and galaxy rendering (PNG)

## Stack
Python (NumPy, Matplotlib). Package architecture (`core.py` + `run.py`).

## Files
- `core.py` — physics engine (`Body`, `Simulation` classes)
- `run.py` — galaxy generation and simulation loop
- `out/` — figures and exported positions

## Run
```bash
python -m galaxy.run
```
> Note: the output path in `run.py` is hard-coded — adapt it to your environment (or make it relative).
