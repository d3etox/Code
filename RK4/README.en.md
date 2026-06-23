[🇫🇷 Version française](README.md)

# Nodal thermal model — Linearized Runge-Kutta 4

> Solving a thermal network (nodal model) with a linearized Runge-Kutta 4 scheme, applied to space thermal engineering.

**Context**: from my internship / apprenticeship at Dorea (thermal modelling for the European Space Agency).

## What it does
Simulation of the time evolution of temperature for a system of **30 coupled thermal nodes**:
- **Conductive** exchanges (node-to-node conductance matrix) and **capacitive** effects
- Non-linear radiation (Stefan-Boltzmann T⁴ law)
- Localized power input, boundary node at fixed temperature
- **RK4** integration with linearization of the radiative term

This extends the linearization methods developed during my internship (presented at ESA's ESTEW 2024 workshop).

## Stack
Python (NumPy, Matplotlib). Original version also in Fortran 90.

## Files
- `RK4_linearized.py` — main solver (parameters at the top of the file)

## Run
```bash
python RK4_linearized.py
```
