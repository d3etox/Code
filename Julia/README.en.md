[🇫🇷 Version française](README.md)

# Optimal control & model predictive control (Julia)

> Optimal-control and model-predictive-control (MPC) problems solved in Julia.

**Context**: control / optimization course (Master in Mathematical Engineering).

## What it does
- **Model Predictive Control (MPC)** (`MPC.jl`): optimal steering of a system (e.g. navigation in a current) formulated as a constrained optimization problem, solved with **JuMP** and **Ipopt**.
- **Reinforcement learning** (`reinforcement.jl`) and other optimization/control scripts.
- Exploration of Julia's scientific ecosystem (JuMP, OrdinaryDiffEq, Plots).

## Stack
**Julia** (JuMP, Ipopt, OrdinaryDiffEq, LinearAlgebra, Plots).

## Files
- `MPC.jl` — model predictive control
- `reinforcement.jl`, `metro.jl`, `compare.jl`, `script*.jl`, `tp2.jl`

## Run
```julia
julia MPC.jl
```
