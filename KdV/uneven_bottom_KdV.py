import numpy as np
import matplotlib.pyplot as plt
import time
import sys 

# Paramètres du domaine
L = 100     # Longueur du domaine spatial
T = 100       # Durée de la simulation
dt = 0.0001  # Pas temporel
dx = 0.1     # Pas spatial
Nx = int(L/dx)   # Nombre de points dans l'espace
Nt = int(T/dt)   # Nombre de pas de temp
x = np.linspace(0, L, Nx)  # Domaine spatial

# Paramètres KdV
alpha = 0.01     # Paramètre non-linéaire
beta = 0.01      # Paramètre dispersif
delta = 0.7      # Paramètre de fond incliné
h0 = -2          # Hauteur du fond avant la pente
h1 = -1.5        # Hauteur du fond après la pente
x0 = 20          # Debut de la pente
x1 = 70          # Fin de la pente
pente = (h1 - h0)/(x1 - x0)
k = np.zeros(Nx) # Pente du fond marin
k[x0:x1] = pente

# Initialisation : une vague solitaire (soliton)
A = 1  # Amplitude
u_4 = np.zeros((Nt,Nx))
u_4[0, :] = A * (1 / np.cosh(np.sqrt(3 / 4 * (alpha / beta) * A) * (x - 10)))**2

# Hauteur du fond marin
def depth(x):
    if x<x0 :
        return h0
    elif x>x1 :
        return h1
    else :
        return (h0 + pente * (x - x0))

start_time = time.time()
# Vectorisation des calculs pour chaque n
for n in range(Nt-1):
    # Termes calculés avec des fenêtres de voisinage sur u (FTCS 4eme ordre)
    u_4jm1 = np.roll(u_4[n], 1)  # Décalage de 1 pour u[j-1]
    u_4jm2 = np.roll(u_4[n], 2)  # Décalage de 2 pour u[j-2]
    u_4jm3 = np.roll(u_4[n], 3)  # Décalage de 3 pour u[j-3]
    u_4jp1 = np.roll(u_4[n], -1)  # Décalage de -1 pour u[j+1]
    u_4jp2 = np.roll(u_4[n], -2)  # Décalage de -2 pour u[j+2]
    u_4jp3 = np.roll(u_4[n], -3)  # Décalage de -3 pour u[j+3]

    # Calcul de u[n+1, :] (FTCS 4eme ordre)
    u_4[n + 1, :] = u_4[n, :] \
        - dt / dx * (1/12 * u_4jm2 - 2/3 * u_4jm1 + 2/3 * u_4jp1 - 1/12 * u_4jp2) \
        - dt / dx * 3/4 * alpha * (1/12 * u_4jm2**2 - 2/3 * u_4jm1**2 + 2/3 * u_4jp1**2 - 1/12 * u_4jp2**2) \
        - dt / (dx**3) * 1/6 * beta * (1/8 * u_4jm3 - u_4jm2 + 13/8 * u_4jm1 - 13/8 * u_4jp1 + u_4jp2 - 1/8 * u_4jp3) \
        + dt/4 * delta * (k * u_4[n, :] + 2 * k * x * (1/12 * u_4jm2 - 2/3 * u_4jm1 + 2/3 * u_4jp1 - 1/12 * u_4jp2))
    
    # Progression et temps estimé restant
    elapsed_time = time.time() - start_time
    progress = ((n+1) / (Nt-1)) * 100
    time_per_iteration = elapsed_time / (n+1)
    remaining_iterations = ((Nt-1) - (n+1))
    estimated_time_remaining = time_per_iteration * remaining_iterations
    minutes, seconds = divmod(int(estimated_time_remaining), 60)
    sys.stdout.write('\rProgression: [{:.0f}%] Temps restant : {} minutes {} secondes'.format(progress, minutes, seconds))
    sys.stdout.flush()

# Calcul du temps pour que la vague atteigne x1
t_x1 = (x1 - 10) / (1 + (1 / 2) * alpha * A)
n_x1 = int(t_x1 / dt)  # Indice temporel correspondant

# Affichage des résultats
fig, axs = plt.subplots(2, 1, figsize=(12, 6), gridspec_kw={'height_ratios': [1.5, 1]}, sharex=True)

# Premier graphique : élévation de la surface (solutions analytiques et numériques superposées)
for n in range(0, Nt, int(5 / dt)):  # Affichage périodique
    axs[0].plot(x, u_4[n, :])

axs[0].axhline(1, color='k', linestyle='--', linewidth=1, label="Amplitude de référence à x=10")  # Ligne horizontale à y=1
axs[0].axhline(np.max(u_4[n_x1,:]), color='r', linestyle='--', linewidth=1, label=f"Amplitude à x₁=70")  # Ligne horizontale pour l'amplitude
axs[0].set_title("Profil du soliton pris périodiquement")
axs[0].set_ylabel("Hauteur du soliton (η)")
axs[0].grid(True)

# Placer la légende en dessous du premier graphique
axs[0].legend(ncol=2, fontsize=10, loc='upper center', bbox_to_anchor=(0.5, 0))  # Position entre les graphes

# Deuxième graphique : profondeur du fond marin
h = np.array([depth(xi) for xi in x])  # Calcul de la profondeur
axs[1].plot(x, h, 'k-', linewidth=2)
axs[1].set_title("Profil du fond marin")
axs[1].set_ylabel("Profondeur")
axs[1].set_xlabel("x")
axs[1].grid(True)

# Ajustement des espaces entre les graphiques
plt.tight_layout()
plt.show()
