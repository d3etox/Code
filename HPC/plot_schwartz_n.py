import numpy as np
import glob
import matplotlib.pyplot as plt

# Liste de tous les fichiers solution_rank*.dat
files = sorted(glob.glob("solution_rank*.dat"))

x_all = []
u_all = []

for file in files:
    data = np.loadtxt(file)
    x_all.extend(data[:, 0])
    u_all.extend(data[:, 1])

x_all = np.array(x_all)
u_all = np.array(u_all)

# Tri par x au cas où les fichiers ne sont pas parfaitement contigus
sort_idx = np.argsort(x_all)
x_all = x_all[sort_idx]
u_all = u_all[sort_idx]

# Plot
plt.figure(figsize=(8,5))
plt.plot(x_all, u_all, '-o', markersize=3, label='Solution Schwarz MPI')
plt.xlabel('x')
plt.ylabel('u(x)')
plt.title('Solution globale reconstruite à partir des fichiers locaux')
plt.grid(True)
plt.legend()
plt.show()
