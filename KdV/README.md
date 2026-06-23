[🇬🇧 English version](README.en.md)

# Équation de Korteweg-de Vries — Solitons & shoaling

> Étude numérique de la montée des ondes solitaires (solitons) sur un fond marin de profondeur variable, via l'équation de Korteweg-de Vries (KdV).

**Contexte** : nano-projet de M1 (Mathématiques et modélisation du vivant), présenté en soutenance.

## Ce que fait le projet
Résolution de l'équation KdV `∂_t u + α u ∂_x u + β ∂_xxx u = 0` par schéma **FTCS** (Forward-Time Centered-Space) :
- Propagation d'un soliton `u = A·sech²(...)` sur fond plat (validation analytique vs numérique)
- Effet d'un **fond incliné** (pente sous-marine) sur l'amplitude de la vague → phénomène de **shoaling**
- Calcul du **coefficient de shoaling** numérique comparé à la prédiction analytique

## Stack
Python (NumPy, Matplotlib).

## Fichiers
- `uneven_bottom_KdV.py` — propagation sur fond incliné
- `shoaling_coefficient_KdV.py` — coefficient de shoaling (analytique vs numérique)
- `ordre_KdV.py` — étude de l'ordre du schéma
- `SLIDES_*.pdf` — slides de soutenance · `*.png` — résultats

## Lancer
```bash
python uneven_bottom_KdV.py
```
