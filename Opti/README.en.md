[🇫🇷 Version française](README.md)

# Optimization

> Optimization, root-finding and metaheuristic algorithms — from course labs to a full constrained black-box optimization project.

**Context**: Optimization course (Master in Mathematical Engineering) and an M2 Numerical Engineering project.

## Flagship project — Optimal wind-farm layout
The [`Parc-eolien`](Parc-eolien) subfolder contains a **constrained black-box optimization** project: placing `N` wind turbines to **maximize annual energy production**, subject to spacing constraints (wake effect) and an allowed-zone constraint, on an objective function with no analytical gradient.

- **Genetic algorithm** built *from scratch* (tournament selection, barycentric crossover, Gaussian mutation, elitism, penalization, fitness cache)
- **Gradient methods** (gradient ascent, adaptive step) with step-sensitivity study
- **Comparison** against the reference derivative-free solver **NOMAD** (MADS algorithm)

→ See the [detailed project README](Parc-eolien/README.md).

## Course labs
- **Newton's method** for root finding (`f(x)=0`) and quadratic convergence study
- Visualization of functions and their roots
- Successive labs exploring different optimization methods

## Stack
Python (NumPy, Matplotlib), Jupyter.

## Contents
- [`Parc-eolien/`](Parc-eolien) — wind-farm optimization project (genetic algorithm + gradient + NOMAD)
- `TP1.ipynb`, `TP2.ipynb`, `TP3.ipynb` — course labs
- `test.py` — Newton implementation

## Run
Open the notebooks in Jupyter, or see the `Parc-eolien` folder for the optimization project.
