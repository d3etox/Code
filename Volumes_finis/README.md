[🇬🇧 English version](README.en.md)

# Méthode des Volumes Finis

> Discrétisation par volumes finis : d'une EDP elliptique 1D simple jusqu'à un modèle stratigraphique de formation de bassins sédimentaires.

**Contexte** : TP du cours de Volumes Finis (M2), encadré à partir de codes de R. Masson.

## Ce que fait le projet
- **TP1 / TP2** : résolution de l'équation elliptique 1D `-u'' = f` par volumes finis, avec conditions de Dirichlet et Neumann ; étude de convergence.
- **Modèle stratigraphique** (`*_strat_litho.py`) : modélisation du transport sédimentaire diffusif à grande échelle d'espace et de temps,
  `∂_t h + div(grad(ψ(b))) = 0`, avec flux **TPFA** sur maillages non structurés. Simulation de la formation de bassins (lithologie unique).

## Stack
Python (NumPy, Matplotlib).

## Fichiers
- `TP1.py`, `TP2.py`, `TP2_offi.py` — équation elliptique 1D
- `1_strat_litho.py`, `2_strat_litho.py` — modèle stratigraphique

## Lancer
```bash
python TP1.py
python 1_strat_litho.py
```
