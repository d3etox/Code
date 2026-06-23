[🇬🇧 English version](README.en.md)

# Simulation N-corps d'une galaxie

> Simulation gravitationnelle de l'interaction entre des centaines de corps (étoiles/planètes) formant une galaxie.

**Contexte** : projet personnel (prolongement de mon projet de prépa sur la simulation de galaxie).

## Ce que fait le projet
Simulation N-corps d'un disque galactique en rotation :
- Calcul des forces gravitationnelles avec **adoucissement** (softening) pour éviter les singularités
- Deux méthodes : calcul naïf **O(N²)** et **approximation par quadrants** (esprit Barnes-Hut) pour accélérer sur grand N
- Intégration **symplectique** (Euler) conservant l'énergie
- Suivi de l'énergie totale, export des trajectoires (CSV) et rendu de la galaxie (PNG)

## Stack
Python (NumPy, Matplotlib). Architecture en package (`core.py` + `run.py`).

## Fichiers
- `core.py` — moteur physique (classes `Body`, `Simulation`)
- `run.py` — génération de la galaxie et boucle de simulation
- `out/` — figures et positions exportées

## Lancer
```bash
python -m galaxy.run
```
> Note : le chemin de sortie dans `run.py` est en dur — à adapter à ton environnement (ou rendre relatif).
