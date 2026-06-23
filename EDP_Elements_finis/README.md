[🇬🇧 English version](README.en.md)

# Éléments finis P1/P2/P3 — Convection-Réaction-Diffusion 2D

> Implémentation et validation d'éléments finis de Lagrange (Pk) pour un problème de convection-réaction-diffusion anisotrope en 2D.

**Contexte** : Devoir maison du cours d'Éléments Finis (Polytech Nice-Sophia, M2).

## Ce que fait le projet
Résolution numérique du problème elliptique
`-div(A∇u) + b·∇u + c·u = f` sur un domaine carré, avec solution exacte connue
`u(x,y) = sin(πx)·cos(2πy)` permettant de **valider la convergence** des éléments P1, P2 et P3.

- Assemblage des matrices de rigidité/masse sur maillage 2D
- Calcul des fonctions de base et de leurs dérivées (P2, P3)
- Étude de l'**ordre de convergence** en norme L² selon le degré des éléments
- Cas tests : diffusion isotrope/anisotrope, convection, réaction

## Stack
Python (NumPy, SciPy `sparse`), Scilab (`.sce`). Rapport LaTeX dans `projet/`.

## Fichiers principaux
- `DM_FABRE_EM.py` — code principal du devoir
- `Laplacian_Pk_2D_2025_V_AF.py` — solveur Pk 2D
- `compute_P2_derivatives.py`, `compute_P3_derivatives.py` — bases d'ordre élevé
- `test_*.py` — tests de validation et de convergence
- `png/` — figures de résultats · `DM_FABRE_EM.pdf` — rapport

## Lancer
```bash
python DM_FABRE_EM.py
```
