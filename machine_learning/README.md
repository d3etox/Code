[🇬🇧 English version](README.en.md)

# IA scientifique — Réseaux de neurones informés par la physique (PINNs)

> Résolution d'équations aux dérivées partielles par réseaux de neurones (Physics-Informed Neural Networks) avec DeepXDE, et étude du compromis profondeur/largeur des réseaux.

**Contexte** : projets de M2 Ingénierie Numérique (deep learning pour la simulation).

## Ce que fait le projet
- **Équations de Maxwell 1D** : résolution du système (E, H) par PINN, y compris avec milieu hétérogène (permittivité `εr` variable) et conditions aux limites personnalisées.
- **Équation de Burgers** : résolution non linéaire, comparaison à un jeu de données de référence.
- **Deep vs Shallow** (`tradeoff_layers_neurons.py`) : étude de l'influence de l'architecture (réseau profond vs peu profond) sur la qualité de l'approximation, à nombre de paramètres comparable — illustre le théorème d'approximation universelle.

## Stack
Python, **DeepXDE**, TensorFlow / PyTorch (backend), NumPy, Matplotlib.

## Fichiers
- `problem1.py`, `problem2.py` — Maxwell 1D (PINN)
- `maxwell_custom_bc.py` — milieu hétérogène + CL custom
- `tradeoff_layers_neurons.py` — étude profondeur vs largeur
- `testburger*.py`, `testdde*.py` — Burgers & EDD · `dataset/` · `*.mp4` — animations

## Lancer
```bash
pip install deepxde
python problem1.py
```
