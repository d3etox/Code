import numpy as np
import matplotlib.pyplot as plt
import time
import sys 

# Paramètres du domaine
L = 100     # Longueur du domaine spatial
T = 100       # Durée de la simulation
dt = 0.005  # Pas temporel
dx = 0.1     # Pas spatial
Nx = int(L/dx)   # Nombre de points dans l'espace
Nt = int(T/dt)   # Nombre de pas de temps
x = np.linspace(0, L, Nx)  # Domaine spatial

# Paramètres KdV
alpha = 0.01     # Paramètre non-linéaire
beta = 0.01      # Paramètre dispersif
delta = 0.7      # Paramètre de fond incliné
h0 = -2          # Hauteur du fond avant la pente
x0 = 20          # Début de la pente
x1 = 70          # Fin de la pente
A = 1            # Amplitude de la vague initiale

# Vecteur de h1
h1_values = np.linspace(-2, -1.4, 13)  # 13 points entre -2 et -1.4
Ks_analytical = []
Ks_numerical = []

# Boucle sur les différentes valeurs de h1
for i, h1 in enumerate(h1_values):
    pente = (h1 - h0) / (x1 - x0)
    k = np.zeros(Nx)  # Pente du fond marin
    k[x0:x1] = pente

    # Initialisation : une vague solitaire (soliton)
    u_4 = np.zeros((Nt, Nx))
    u_4[0, :] = A * (1 / np.cosh(np.sqrt(3 / 4 * (alpha / beta) * A) * (x - 10)))**2

    # Simulation
    start_time = time.time()  # Temps de départ pour chaque simulation
    for n in range(Nt-1):
        u_4jm1 = np.roll(u_4[n], 1)  # Décalage de 1 pour u[j-1]
        u_4jm2 = np.roll(u_4[n], 2)  # Décalage de 2 pour u[j-2]
        u_4jm3 = np.roll(u_4[n], 3)  # Décalage de 3 pour u[j-3]
        u_4jp1 = np.roll(u_4[n], -1)  # Décalage de -1 pour u[j+1]
        u_4jp2 = np.roll(u_4[n], -2)  # Décalage de -2 pour u[j+2]
        u_4jp3 = np.roll(u_4[n], -3)  # Décalage de -3 pour u[j+3]

        u_4[n + 1, :] = u_4[n, :] \
            - dt / dx * (1/12 * u_4jm2 - 2/3 * u_4jm1 + 2/3 * u_4jp1 - 1/12 * u_4jp2) \
            - dt / dx * 3/4 * alpha * (1/12 * u_4jm2**2 - 2/3 * u_4jm1**2 + 2/3 * u_4jp1**2 - 1/12 * u_4jp2**2) \
            - dt / (dx**3) * 1/6 * beta * (1/8 * u_4jm3 - u_4jm2 + 13/8 * u_4jm1 - 13/8 * u_4jp1 + u_4jp2 - 1/8 * u_4jp3) \
            + dt / 4 * delta * (k * u_4[n, :] + 2 * k * x * (1/12 * u_4jm2 - 2/3 * u_4jm1 + 2/3 * u_4jp1 - 1/12 * u_4jp2))
        
        # Progression de la simulation interne
        elapsed_time = time.time() - start_time
        progress = ((n + 1) / (Nt - 1)) * 100
        time_per_iteration = elapsed_time / (n + 1)
        remaining_iterations = ((Nt - 1) - (n + 1))
        estimated_time_remaining = time_per_iteration * remaining_iterations
        minutes, seconds = divmod(int(estimated_time_remaining), 60)

        sys.stdout.write(
            f"\rSimulation {i+1}/{len(h1_values)} | Progression: [{progress:.0f}%] Temps restant : {minutes} minutes {seconds} secondes"
        )
        sys.stdout.flush()

    # Calcul des coefficients Ks
    Ks_analytical.append((h0 / h1)**(1/4))
    Ks_numerical.append(np.max(u_4[Nt-1]) / A)

# Tracé des résultats
plt.figure(figsize=(10, 6))
plt.plot(h1_values - h0, Ks_analytical, 'k--', label="Analytical $K_s$")
plt.plot(h1_values - h0, Ks_numerical, 'b*', label="Numerical $K_s$")
plt.xlabel("Slope's Height $(h_1 - h_0)$")
plt.ylabel("Shoaling Coefficient ($K_s$)")
plt.legend()
plt.grid(True)
plt.title("Shoaling Coefficient vs Slope's Height")
plt.show()
