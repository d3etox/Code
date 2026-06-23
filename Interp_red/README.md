[🇬🇧 English version](README.en.md)

# Réduction de dimension & clustering de champs de température

> Analyse et réduction de données de champs de température tridimensionnels : décomposition en valeurs singulières et regroupement (clustering).

**Contexte** : projet de M2 Ingénierie Numérique (interpolation / réduction de modèle), prolongement de mon projet de MAM4 sur l'interpolation de champs 3D par transformée de Fourier.

## Ce que fait le projet
- **SVD / méthode de la puissance itérée** (`SVD.py`) : calcul des valeurs propres dominantes d'une matrice — brique de base de la réduction de dimension (POD/PCA).
- **Clustering** (`clustering.py`) : regroupement de points d'un champ de température 3D par **K-means**, **clustering spectral**, **MeanShift** et **Birch**, avec visualisation 3D.

Compétences : analyse de données, algèbre linéaire numérique, scikit-learn.

## Stack
Python (NumPy, scikit-learn, Matplotlib).

## Fichiers
- `SVD.py` — méthode de la puissance itérée
- `clustering.py` — comparaison d'algorithmes de clustering
- `*.csv` — jeux de données de température

## Lancer
```bash
python clustering.py
```
