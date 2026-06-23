[🇬🇧 English version](README.en.md)

# Calcul Haute Performance — Méthode de Schwarz (MPI)

> Résolution parallèle d'un problème de Poisson 1D par décomposition de domaine (méthode de Schwarz), implémentée en Fortran 90 avec MPI.

**Contexte** : cours de Calcul Haute Performance (M2 Ingénierie Numérique).

## Ce que fait le projet
Décomposition du domaine en sous-domaines répartis sur plusieurs processus MPI :
- **Résolution globale** (référence séquentielle)
- **Méthode de Schwarz** parallèle : chaque processus résout son sous-domaine et échange ses conditions aux interfaces avec ses voisins jusqu'à convergence
- **Benchmark** de scalabilité sur N processus (mesure du speed-up)

Démontre la maîtrise de la **programmation parallèle distribuée** (échanges `MPI_Send`/`MPI_Recv`, `MPI_Bcast`) — compétence clé pour la simulation à grande échelle.

## Stack
**Fortran 90 + MPI**, Python (post-traitement / visualisation).

## Fichiers
- `Schwartz_advanced.f90` — programme principal (3 modes : global / Schwarz MPI / benchmark)
- `plot_schwartz*.py`, `tracer.py` — visualisation des solutions
- `solution*.dat` — sorties

## Lancer
```bash
mpif90 Schwartz_advanced.f90 -o schwartz
mpirun -np 4 ./schwartz
python plot_schwartz.py
```
