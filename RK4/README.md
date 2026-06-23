[🇬🇧 English version](README.en.md)

# Modèle thermique nodal — Runge-Kutta 4 linéarisé

> Résolution d'un réseau thermique (modèle nodal) par schéma de Runge-Kutta 4 linéarisé, appliqué à la thermique spatiale.

**Contexte** : issu de mon stage / alternance chez Dorea (modélisation thermique pour l'Agence Spatiale Européenne).

## Ce que fait le projet
Simulation de l'évolution temporelle de la température d'un système de **30 nœuds** thermiques couplés :
- Échanges **conductifs** (matrice de conductance entre nœuds) et **capacitifs**
- Rayonnement non linéaire (loi de Stefan-Boltzmann en T⁴)
- Apport de puissance localisé, nœud frontière à température imposée
- Intégration par **RK4** avec linéarisation du terme radiatif

Ce travail prolonge les méthodes de linéarisation développées pendant mon stage (présenté au workshop ESTEW 2024 de l'ESA).

## Stack
Python (NumPy, Matplotlib). Version originale également en Fortran 90.

## Fichiers
- `RK4_linearized.py` — solveur principal (paramètres en tête de fichier)

## Lancer
```bash
python RK4_linearized.py
```
