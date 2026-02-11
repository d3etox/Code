import numpy as np
import matplotlib.pyplot as plt
import sys
import time

# Paramètres du domaine spatial et temporel
L = 50       # Longueur du domaine spatial
T = 1        # Durée de la simulation

# Paramètres KdV
alpha = 0.01   # Paramètre non-linéaire
beta = 0.00625 # Paramètre dispersif
A = 1          # Amplitude du soliton

# Fonction pour la solution analytique
def soliton(x, t, A, alpha, beta):
    c = 1 + (1 / 2) * alpha * A
    kappa = np.sqrt(3 / 4 * (alpha / beta) * A)
    return A * (1 / np.cosh(kappa * (x - 10 - c * t)))**2

# Test de convergence en espace
dx_values = np.linspace(0.1,0.01,10)  # Différents pas spatiaux
errors = []
i = 0
for dx in dx_values:
    i += 1
    Nx = int(L / dx)
    x = np.linspace(0, L, Nx)  # Domaine spatial
    dt = dx**4
    Nt = int(T / dt)
    # Initialisation
    u = np.zeros((Nt,Nx))
    u_err = np.zeros((Nt,Nx))
    u[0, :] = soliton(x, 0, A, alpha, beta)

    start_time = time.time()  # Temps de départ pour chaque simulation
    # Boucle temporelle
    for n in range(Nt - 1):
        # Décalages pour le schéma FTCS 4ème ordre
        u_jm1 = np.roll(u[n], 1)  # u[j-1]
        u_jm2 = np.roll(u[n], 2)  # u[j-2]
        u_jm3 = np.roll(u[n], 3)  # u[j-3]
        u_jp1 = np.roll(u[n], -1) # u[j+1]
        u_jp2 = np.roll(u[n], -2) # u[j+2]
        u_jp3 = np.roll(u[n], -3) # u[j+3]

        # Schéma FTCS 4ème ordre
        u[n+1, :] = u[0, :] \
            - dt / dx * (1/12 * u_jm2 - 2/3 * u_jm1 + 2/3 * u_jp1 - 1/12 * u_jp2) \
            - dt / dx * 3/4 * alpha * (1/12 * u_jm2**2 - 2/3 * u_jm1**2 + 2/3 * u_jp1**2 - 1/12 * u_jp2**2) \
            - dt / (dx**3) * 1/6 * beta * (1/8 * u_jm3 - u_jm2 + 13/8 * u_jm1 - 13/8 * u_jp1 + u_jp2 - 1/8 * u_jp3)
        
        u_err = u[n+1, :] - soliton(x, n*dt, A, alpha, beta)

        # Progression de la simulation interne
        elapsed_time = time.time() - start_time
        progress = ((n + 1) / (Nt - 1)) * 100
        time_per_iteration = elapsed_time / (n + 1)
        remaining_iterations = ((Nt - 1) - (n + 1))
        estimated_time_remaining = time_per_iteration * remaining_iterations
        minutes, seconds = divmod(int(estimated_time_remaining), 60)

        sys.stdout.write(
            f"\rSimulation {i}/{len(dx_values)} | Progression: [{progress:.0f}%] Temps restant : {minutes} minutes {seconds} secondes"
        )
        sys.stdout.flush()

    # Calcul de l'erreur L2
    error = np.sqrt(np.sum((u - u_err)**2))
    errors.append(error)

# Analyse des erreurs
plt.figure(figsize=(8, 6))
dx_values_log = np.log(dx_values)
errors_log = np.log(errors)

# Ajustement linéaire pour estimer l'ordre
coeffs = np.polyfit(dx_values_log, errors_log, 1)
order = -coeffs[0]

plt.loglog(dx_values, errors, 'o-', label=f"Pente = {order:.2f}")
plt.xlabel("dx (pas spatial)")
plt.ylabel("Erreur L2")
plt.title("Test de convergence : erreur en fonction de dx")
plt.grid(True, which="both", linestyle="--", linewidth=0.5)
plt.legend()
plt.show()

print(f"Ordre en espace estimé : {order:.2f}")