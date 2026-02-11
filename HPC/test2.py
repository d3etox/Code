import numpy as np
import time
import matplotlib.pyplot as plt

def run_sim(N):

    # Paramètres physiques
    c = 1.0
    L = 1.0
    T = 0.1     # réduit pour éviter des temps énormes

    # Discrétisation
    dx = L / (N - 1)
    beta = 0.9
    dt = beta * dx / c
    nt = int(T / dt)

    x = np.linspace(0.0, L, N)

    # Allocation
    p_prev = np.zeros(N)
    p = np.zeros(N)
    p_next = np.zeros(N)

    # Condition initiale
    x0 = 0.5 * L
    sigma = 0.05
    p_prev[:] = np.exp(-((x - x0)**2) / (2.0 * sigma**2))

    beta2 = (c * dt / dx)**2

    # Construction de p^1
    p[:] = p_prev
    p[0] = p_prev[0] + 0.5 * beta2 * (2*p_prev[1] - 2*p_prev[0])
    p[-1] = p_prev[-1] + 0.5 * beta2 * (2*p_prev[-2] - 2*p_prev[-1])

    for i in range(1, N-1):
        lap = p_prev[i-1] - 2*p_prev[i] + p_prev[i+1]
        p[i] = p_prev[i] + 0.5 * beta2 * lap

    # Boucle en temps (non vectorisée)
    start = time.time()
    for _ in range(1, nt):

        for i in range(1, N-1):
            p_next[i] = (
                -p_prev[i]
                + 2*(1-beta2)*p[i]
                + beta2*(p[i-1] + p[i+1])
            )

        p_next[0] = -p_prev[0] + 2*(1-beta2)*p[0] + 2*beta2*p[1]
        p_next[-1] = -p_prev[-1] + 2*(1-beta2)*p[-1] + 2*beta2*p[-2]

        p_prev, p, p_next = p, p_next, p_prev

    end = time.time()
    return end - start


# --- tailles testées ---
Ns = np.linspace(201, 3001, 10)
Ns = Ns.astype(int)  
times = [run_sim(N) for N in Ns]


# --- tracé ---
plt.figure(figsize=(7,4))
plt.plot(Ns, times, marker='o', color="orange")
plt.xlabel("N : taille du maillage")
plt.ylabel("Temps d'exécution (s)")
plt.title("Complexité du schéma d'onde non vectorisé (en espace)")
plt.grid(True)
plt.show()
