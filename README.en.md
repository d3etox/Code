[🇫🇷 Version française](README.md)

# Code — Scientific computing & numerical simulation portfolio

A collection of my projects in **mathematical modelling**, **numerical methods**, **high-performance computing** and **scientific machine learning**, carried out during my Master's in Mathematical Engineering (Université Côte d'Azur / Polytech Nice-Sophia) and as personal projects.

> Profile: M2 student in Numerical Engineering, apprentice at Dorea (space thermal modelling, ESA). Strong interest in HPC and numerical simulation.
> Main stack: **Python** (NumPy, SciPy, Matplotlib, PyTorch, DeepXDE), **Fortran 90 / MPI**, **C++/Qt**, **Julia**, **R**.

---

## Project overview by theme

### PDE solving & numerical methods
| Project | Description | Method |
|---|---|---|
| [`EDP_Elements_finis`](EDP_Elements_finis) | 2D anisotropic convection-reaction-diffusion | Finite elements P1/P2/P3 |
| [`Volumes_finis`](Volumes_finis) | 1D elliptic equation & stratigraphic basin model | Finite volumes (TPFA) |
| [`maths_physique`](maths_physique) | Advection equation: scheme comparison | Upwind, Lax-Friedrichs, Lax-Wendroff |
| [`KdV`](KdV) | Soliton propagation over a variable seabed (shoaling) | Korteweg-de Vries, FTCS |
| [`Modèle Spatio-temporel`](Modèle%20Spatio-temporel) | Diffusion & reaction-diffusion (Fisher-KPP, Kolmogorov) | Finite differences, mathematical biology |
| [`Population_structure`](Population_structure) | Prion proliferation dynamics (Creutzfeldt-Jakob) | ODE system, time integration |
| [`RK4`](RK4) | Nodal thermal model (space thermal, Dorea internship) | Linearized Runge-Kutta 4 |

### High-Performance Computing (HPC)
| Project | Description | Method |
|---|---|---|
| [`HPC`](HPC) | Parallel domain decomposition (Poisson problem) | Schwarz method, **Fortran 90 + MPI** |

### Scientific ML & model reduction
| Project | Description | Method |
|---|---|---|
| [`machine_learning`](machine_learning) | Solving PDEs with neural networks (Maxwell 1D, Burgers) | **PINNs** (DeepXDE), deep vs shallow |
| [`Neural_network`](Neural_network) | PyTorch & deep learning basics | Tensors, deep networks |
| [`Interp_red`](Interp_red) | Dimensionality reduction & clustering on 3D temperature fields | SVD, K-means, spectral clustering |

### Optimization & control
| Project | Description | Method |
|---|---|---|
| [`Opti`](Opti) | Wind-farm optimization & root finding | Genetic algorithm, gradient, NOMAD, Newton |
| [`Julia`](Julia) | Optimal control & model predictive control | **MPC** (JuMP/Ipopt), reinforcement learning |

### Uncertainty quantification & statistics
| Project | Description | Method |
|---|---|---|
| [`R_studio`](R_studio) | Global sensitivity analysis (Borehole model) | Monte Carlo, **Sobol indices** |

### Personal projects (graphics & games)
| Project | Description | Method |
|---|---|---|
| [`galaxy`](galaxy) | N-body galaxy simulation | Gravitation, quadrant approximation |
| [`raytracer`](raytracer) | 3D ray-tracing rendering engine | Ray tracing, **Numba** acceleration |
| [`Perso`](Perso) | Connect 4 (minimax AI) & card-game prototype | Minimax, GUIs |

---

## Organization

Each folder contains its own `README.md` detailing the context, the methods used and how to run the code.

Most Python projects rely on `numpy`, `scipy` and `matplotlib`. Project-specific dependencies (PyTorch, DeepXDE, Numba, MPI…) are listed in each project's README.
