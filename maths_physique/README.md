[🇬🇧 English version](README.en.md)

# Schémas numériques pour l'équation d'advection

> Comparaison de schémas aux différences finies pour l'équation de transport linéaire `∂_t u + a ∂_x u = 0`.

**Contexte** : projet de méthodes numériques (mathématiques pour la physique).

## Ce que fait le projet
Implémentation et comparaison de 6 schémas sur une condition initiale discontinue :
- Upwind (décentré amont / aval)
- Schéma centré
- **Lax-Friedrichs**
- **Lax-Wendroff**
- Upwind d'ordre 2

Analyse de la **diffusion numérique**, des **oscillations** (Gibbs) et de l'**ordre de convergence** (erreur L² pour différentes résolutions `dx`).

## Stack
Python (NumPy, Matplotlib).

## Fichiers
- `maths_physique.py` — schémas et visualisation
- `maths_physique2.py` — étude de convergence (erreur L² vs résolution)

## Lancer
```bash
python maths_physique.py
```
