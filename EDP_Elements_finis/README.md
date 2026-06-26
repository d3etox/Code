[🇬🇧 English version](README.en.md)

# Éléments finis P1/P2/P3 — Convection-Réaction-Diffusion 2D & Navier-Stokes

> Implémentation et validation d'éléments finis de Lagrange (Pk) pour des problèmes elliptiques 2D, jusqu'aux équations de Navier-Stokes incompressibles.

**Contexte** : cours et devoir maison d'Éléments Finis (Polytech Nice-Sophia, M2).

## Ce que fait le projet
- **Devoir maison — convection-réaction-diffusion 2D** : résolution de `-div(A∇u) + b·∇u + c·u = f` sur un carré, avec solution exacte connue `u(x,y) = sin(πx)·cos(2πy)` pour **valider la convergence** des éléments P1, P2 et P3 (cas isotrope/anisotrope, convection, réaction).
- **Assemblage Pk** : fonctions de base et dérivées P2/P3, étude de l'**ordre de convergence** en norme L².
- **Éléments 1D** (`TDTP_01_Pk_1D.py`, `TP_Cours1.ipynb`) : matrices locales P1/P2 Lagrange.
- **Navier-Stokes 2D** (`TP_NavierStokes.ipynb`) : résolution des équations incompressibles, validée sur l'écoulement de **Taylor-Green** (solution analytique instationnaire).

## Stack
Python (NumPy, SciPy `sparse`), Scilab (`.sce`), Jupyter. Rapport LaTeX dans `projet/`.

## Fichiers principaux
- `DM_FABRE_EM.py` — code du devoir maison (convection-réaction-diffusion)
- `Laplacian_Pk_2D_2025_V_AF.py` — solveur Pk 2D
- `compute_P2_derivatives.py`, `compute_P3_derivatives.py`, `P3_analysis.py` — bases d'ordre élevé
- `TP_NavierStokes.ipynb` — Navier-Stokes 2D (Taylor-Green)
- `TDTP_01_Pk_1D.py`, `TP_Cours1.ipynb`, `TP1.ipynb` — TP éléments 1D
- `test_*.py` — tests de validation et de convergence
- `png/` — figures de résultats · `DM_FABRE_EM.pdf` — rapport

## Lancer
```bash
python DM_FABRE_EM.py
```
