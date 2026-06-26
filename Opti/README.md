[🇬🇧 English version](README.en.md)

# Optimisation

> Algorithmes d'optimisation, de recherche de racines et métaheuristiques — des TP de cours à un projet complet d'optimisation sous contraintes en boîte noire.

**Contexte** : cours d'Optimisation (Master Ingénierie Mathématique) et projet de M2 Ingénierie Numérique.

## Projet phare — Implantation optimale d'un parc éolien
Le sous-dossier [`Parc-eolien`](Parc-eolien) contient un projet d'**optimisation sous contraintes en boîte noire** : placer `N` éoliennes pour **maximiser la production annuelle d'énergie**, sous contraintes d'espacement (effet de sillage) et de zone autorisée, sur une fonction objectif sans gradient analytique.

- **Algorithme génétique** codé *from scratch* (sélection par tournoi, croisement barycentrique, mutation gaussienne, élitisme, pénalisation, cache de fitness)
- **Méthodes de gradient** (montée de gradient, pas adaptatif) avec étude de sensibilité au pas
- **Comparaison** au solveur sans dérivées de référence **NOMAD** (algorithme MADS)

→ Voir le [README détaillé du projet](Parc-eolien/README.md).

## TP de cours
- **Méthode de Newton** pour la recherche de racines (`f(x)=0`) et étude de la convergence quadratique
- Visualisation de fonctions et de leurs racines
- TP successifs explorant différentes méthodes d'optimisation

## Stack
Python (NumPy, Matplotlib), Jupyter.

## Contenu
- [`Parc-eolien/`](Parc-eolien) — projet d'optimisation d'un parc éolien (algorithme génétique + gradient + NOMAD)
- `TP1.ipynb`, `TP2.ipynb`, `TP3.ipynb` — TP du cours
- `test.py` — implémentation de Newton

## Lancer
Ouvrir les notebooks dans Jupyter, ou consulter le dossier `Parc-eolien` pour le projet d'optimisation.
