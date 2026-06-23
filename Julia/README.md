[🇬🇧 English version](README.en.md)

# Contrôle optimal & commande prédictive (Julia)

> Problèmes de contrôle optimal et de commande prédictive (MPC) résolus en Julia.

**Contexte** : cours de contrôle / optimisation (Master Ingénierie Mathématique).

## Ce que fait le projet
- **Commande prédictive (MPC)** (`MPC.jl`) : pilotage optimal d'un système (ex. navigation dans un courant) formulé comme un problème d'optimisation sous contraintes, résolu avec **JuMP** et **Ipopt**.
- **Reinforcement learning** (`reinforcement.jl`) et autres scripts d'optimisation/contrôle.
- Découverte de l'écosystème scientifique Julia (JuMP, OrdinaryDiffEq, Plots).

## Stack
**Julia** (JuMP, Ipopt, OrdinaryDiffEq, LinearAlgebra, Plots).

## Fichiers
- `MPC.jl` — commande prédictive
- `reinforcement.jl`, `metro.jl`, `compare.jl`, `script*.jl`, `tp2.jl`

## Lancer
```julia
julia MPC.jl
```
