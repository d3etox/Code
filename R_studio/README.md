[🇬🇧 English version](README.en.md)

# Quantification d'incertitude & analyse de sensibilité (R)

> Analyse de sensibilité globale et propagation d'incertitudes sur le modèle Borehole, par méthodes de Monte-Carlo et indices de Sobol.

**Contexte** : cours d'analyse de sensibilité de B. Iooss (Polytech Nice-Sophia, 2026).

## Ce que fait le projet
Sur le modèle de référence **Borehole** (débit d'eau à travers un forage, 8 variables d'entrée) :
- Propagation d'incertitude par **Monte-Carlo** (étude de convergence, histogrammes)
- **Indices de Sobol** : décomposition de la variance pour hiérarchiser l'influence des paramètres d'entrée
- Scatterplots et diagnostics graphiques

Compétences : statistiques, quantification d'incertitude (UQ), R.

## Stack
**R** (RStudio).

## Fichiers
- `borehole.R` — analyse complète (auteur : Charles Fabre)
- `TP4-*.R`, `TP5sol-*.R` — TP du cours
- `*.png`, `Rplots.pdf` — résultats (convergence MC, indices de Sobol…)

## Lancer
```r
source("borehole.R")
```
