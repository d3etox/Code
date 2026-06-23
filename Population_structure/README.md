[🇬🇧 English version](README.en.md)

# Modélisation de la dynamique des prions

> Modélisation numérique des mécanismes de prolifération des prions, responsables de l'encéphalopathie spongiforme (maladie de Creutzfeldt-Jakob).

**Contexte** : projet de M1 (Mathématiques et modélisation du vivant).

## Ce que fait le projet
Étude d'un système d'équations différentielles décrivant l'interaction entre la protéine saine **PrPC** et sa forme pathologique **PrPSc**, capable de se répliquer en « contaminant » la forme saine.
- Démonstration d'existence/positivité des solutions sur un domaine admissible
- Modèle de croissance logistique (solution analytique vs numérique)
- Intégration numérique du système et étude de la dynamique de population

## Stack
Python (NumPy, Matplotlib), Jupyter.

## Fichiers
- `Modélisation_Prion_CharlesFabre.ipynb` — notebook principal (théorie + simulations)
- `test.py`, `test.ipynb` — prototypes

## Lancer
Ouvrir `Modélisation_Prion_CharlesFabre.ipynb` dans Jupyter.
