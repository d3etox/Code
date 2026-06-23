[🇫🇷 Version française](README.md)

# Dimensionality reduction & clustering of temperature fields

> Analysis and reduction of three-dimensional temperature-field data: singular value decomposition and clustering.

**Context**: M2 Numerical Engineering project (interpolation / model reduction), extending my MAM4 project on 3D-field interpolation via Fourier transform.

## What it does
- **SVD / power iteration** (`SVD.py`): computing the dominant eigenvalues of a matrix — a building block of dimensionality reduction (POD/PCA).
- **Clustering** (`clustering.py`): grouping points of a 3D temperature field with **K-means**, **spectral clustering**, **MeanShift** and **Birch**, with 3D visualization.

Skills: data analysis, numerical linear algebra, scikit-learn.

## Stack
Python (NumPy, scikit-learn, Matplotlib).

## Files
- `SVD.py` — power iteration method
- `clustering.py` — comparison of clustering algorithms
- `*.csv` — temperature datasets

## Run
```bash
python clustering.py
```
