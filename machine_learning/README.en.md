[🇫🇷 Version française](README.md)

# Scientific ML — Physics-Informed Neural Networks (PINNs)

> Solving partial differential equations with neural networks (Physics-Informed Neural Networks) using DeepXDE, and studying the depth/width trade-off of networks.

**Context**: M2 Numerical Engineering projects (deep learning for simulation).

## What it does
- **1D Maxwell equations**: solving the (E, H) system with a PINN, including a heterogeneous medium (variable permittivity `εr`) and custom boundary conditions.
- **Burgers equation**: non-linear solving, comparison against a reference dataset.
- **Deep vs Shallow** (`tradeoff_layers_neurons.py`): study of the influence of architecture (deep vs shallow network) on approximation quality, at comparable parameter count — illustrates the universal approximation theorem.

## Stack
Python, **DeepXDE**, TensorFlow / PyTorch (backend), NumPy, Matplotlib.

## Files
- `problem1.py`, `problem2.py` — 1D Maxwell (PINN)
- `maxwell_custom_bc.py` — heterogeneous medium + custom BC
- `tradeoff_layers_neurons.py` — depth vs width study
- `testburger*.py`, `testdde*.py` — Burgers & DDE · `dataset/` · `*.mp4` — animations

## Run
```bash
pip install deepxde
python problem1.py
```
