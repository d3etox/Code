[🇫🇷 Version française](README.md)

# Korteweg-de Vries equation — Solitons & shoaling

> Numerical study of solitary-wave run-up (solitons) over a seabed of varying depth, using the Korteweg-de Vries (KdV) equation.

**Context**: M1 nano-project (Mathematics and modelling of living systems), presented in an oral defense.

## What it does
Solving the KdV equation `∂_t u + α u ∂_x u + β ∂_xxx u = 0` with an **FTCS** (Forward-Time Centered-Space) scheme:
- Propagation of a soliton `u = A·sech²(...)` over a flat bottom (analytical vs numerical validation)
- Effect of a **sloping bottom** on wave amplitude → **shoaling** phenomenon
- Numerical **shoaling coefficient** compared to the analytical prediction

## Stack
Python (NumPy, Matplotlib).

## Files
- `uneven_bottom_KdV.py` — propagation over a sloping bottom
- `shoaling_coefficient_KdV.py` — shoaling coefficient (analytical vs numerical)
- `ordre_KdV.py` — scheme order study
- `SLIDES_*.pdf` — defense slides · `*.png` — results

## Run
```bash
python uneven_bottom_KdV.py
```
