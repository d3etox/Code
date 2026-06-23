[🇬🇧 English version](README.en.md)

# Code — Portfolio de calcul scientifique & simulation numérique

Recueil de mes projets en **modélisation mathématique**, **méthodes numériques**, **calcul haute performance** et **IA scientifique**, réalisés pendant mon Master Ingénierie Mathématique (Université Côte d'Azur / Polytech Nice-Sophia) et en projets personnels.

> Profil : étudiant en M2 Ingénierie Numérique, en alternance chez Dorea (thermique spatiale, ESA). Intérêt marqué pour le HPC et la simulation numérique.
> Stack principale : **Python** (NumPy, SciPy, Matplotlib, PyTorch, DeepXDE), **Fortran 90 / MPI**, **C++/Qt**, **Julia**, **R**.

---

## Aperçu des projets par thème

### Résolution d'EDP & méthodes numériques
| Projet | Description | Méthode |
|---|---|---|
| [`EDP_Elements_finis`](EDP_Elements_finis) | Convection-réaction-diffusion 2D anisotrope | Éléments finis P1/P2/P3 |
| [`Volumes_finis`](Volumes_finis) | Équation elliptique 1D & modèle stratigraphique de bassins sédimentaires | Volumes finis (TPFA) |
| [`maths_physique`](maths_physique) | Équation d'advection : comparaison de schémas | Upwind, Lax-Friedrichs, Lax-Wendroff |
| [`KdV`](KdV) | Propagation de solitons sur fond marin variable (shoaling) | Korteweg-de Vries, FTCS |
| [`Modèle Spatio-temporel`](Modèle%20Spatio-temporel) | Diffusion & réaction-diffusion (Fisher-KPP, Kolmogorov) | Différences finies, modélisation du vivant |
| [`Population_structure`](Population_structure) | Dynamique de prolifération des prions (Creutzfeldt-Jakob) | Système d'EDO, schémas d'intégration |
| [`RK4`](RK4) | Modèle thermique nodal (thermique spatiale, stage Dorea) | Runge-Kutta 4 linéarisé |

### Calcul haute performance (HPC)
| Projet | Description | Méthode |
|---|---|---|
| [`HPC`](HPC) | Décomposition de domaine parallèle (problème de Poisson) | Méthode de Schwarz, **Fortran 90 + MPI** |

### IA scientifique & réduction de modèle
| Projet | Description | Méthode |
|---|---|---|
| [`machine_learning`](machine_learning) | Résolution d'EDP par réseaux de neurones (Maxwell 1D, Burgers) | **PINNs** (DeepXDE), deep vs shallow |
| [`Neural_network`](Neural_network) | Bases de PyTorch & deep learning | Tenseurs, réseaux profonds |
| [`Interp_red`](Interp_red) | Réduction de dimension & clustering sur champs de température 3D | SVD, K-means, clustering spectral |

### Optimisation & contrôle
| Projet | Description | Méthode |
|---|---|---|
| [`Opti`](Opti) | Recherche de racines et optimisation | Méthode de Newton |
| [`Julia`](Julia) | Contrôle optimal & commande prédictive | **MPC** (JuMP/Ipopt), reinforcement learning |

### Quantification d'incertitude & statistiques
| Projet | Description | Méthode |
|---|---|---|
| [`R_studio`](R_studio) | Analyse de sensibilité (modèle Borehole) | Monte-Carlo, **indices de Sobol** |

### Projets personnels (graphisme & jeux)
| Projet | Description | Méthode |
|---|---|---|
| [`galaxy`](galaxy) | Simulation N-corps d'une galaxie | Gravitation, approximation par quadrants |
| [`raytracer`](raytracer) | Moteur de rendu 3D par lancer de rayons | Ray tracing, accélération **Numba** |
| [`Perso`](Perso) | Puissance 4 (IA minimax) & prototype de jeu de cartes | Minimax, interfaces graphiques |

---

## Organisation

Chaque dossier contient son propre `README.md` détaillant le contexte, les méthodes employées et la manière de lancer le code.

La plupart des projets Python reposent sur `numpy`, `scipy` et `matplotlib`. Les dépendances spécifiques (PyTorch, DeepXDE, Numba, MPI…) sont indiquées dans le README du projet concerné.
