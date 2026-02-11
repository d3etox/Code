import numpy as np
import matplotlib.pyplot as plt
import glob

# Lire tous les fichiers de rang
files = sorted(glob.glob('solution_rank*.dat'))

x_all = []
u_all = []

for f in files:
    data = np.loadtxt(f)
    x_all.append(data[:,0])
    u_all.append(data[:,1])

x = np.concatenate(x_all)
u = np.concatenate(u_all)

plt.figure(figsize=(12,5))
plt.plot(x, u, label='Solution u(x)', marker='o', markersize=3)

# Paramètres Schwarz
N = len(x) - 2
nprocs = len(files)
overlap_ratio = 0.2
Nwork = N // nprocs
overlap = max(2, int(overlap_ratio * Nwork))

for rank in range(nprocs):
    start = rank * Nwork
    end = start + Nwork
    plt.axvspan(x[start], x[start+overlap], color='orange', alpha=0.3, label='Overlap Left' if rank==0 else "")
    plt.axvspan(x[end-overlap], x[end], color='green', alpha=0.3, label='Overlap Right' if rank==0 else "")
    plt.axvspan(x[start+overlap], x[end-overlap], color='blue', alpha=0.1, label='Zone interne' if rank==0 else "")

plt.xlabel('x')
plt.ylabel('u(x)')
plt.title('Méthode de Schwarz MPI 1D')
plt.legend()
plt.grid(True)
plt.show()
