[🇫🇷 Version française](README.md)

# High-Performance Computing — Schwarz method (MPI)

> Parallel solving of a 1D Poisson problem via domain decomposition (Schwarz method), implemented in Fortran 90 with MPI.

**Context**: High-Performance Computing course (M2 Numerical Engineering).

## What it does
Domain split into subdomains distributed across several MPI processes:
- **Global solve** (sequential reference)
- Parallel **Schwarz method**: each process solves its subdomain and exchanges interface conditions with its neighbours until convergence
- **Benchmark** of scalability over N processes (speed-up measurement)

Demonstrates mastery of **distributed parallel programming** (`MPI_Send`/`MPI_Recv`, `MPI_Bcast`) — a key skill for large-scale simulation.

## Stack
**Fortran 90 + MPI**, Python (post-processing / visualization).

## Files
- `Schwartz_advanced.f90` — main program (3 modes: global / Schwarz MPI / benchmark)
- `plot_schwartz*.py`, `tracer.py` — solution plotting
- `solution*.dat` — outputs

## Run
```bash
mpif90 Schwartz_advanced.f90 -o schwartz
mpirun -np 4 ./schwartz
python plot_schwartz.py
```
