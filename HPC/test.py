import numpy as np
import time
import matplotlib.pyplot as plt

def run_sim(N):
    # --- paramètres physiques ---
    c = 1.0
    L = 1.0
    T = 0.05  # petit temps pour des timings rapides

    # --- discrétisation ---
    dx = L / (N - 1)
    beta = 0.9
    dt = beta * dx / c
    nt = int(T / dt)

    # --- allocation ---
    x = np.linspace(0.0, L, N)
    p_prev = np.exp(-((x - 0.5 * L) ** 2) / (2.0 * 0.05**2))
    p = p_prev.copy()
    p_next = np.zeros(N)
    beta2 = (c * dt / dx) ** 2

    # p^1
    p[1:-1] = p_prev[1:-1] + 0.5*beta2*(p_prev[:-2] - 2*p_prev[1:-1] + p_prev[2:])
    p[0] = p_prev[0] + 0.5*beta2*(2*p_prev[1] - 2*p_prev[0])
    p[-1] = p_prev[-1] + 0.5*beta2*(2*p_prev[-2] - 2*p_prev[-1])

    # --- boucle temporelle ---
    start = time.time()
    for _ in range(1, nt):
        p_next[1:-1] = (
            -p_prev[1:-1]
            + 2*(1-beta2)*p[1:-1]
            + beta2*(p[:-2] + p[2:])
        )
        # bords
        p_next[0] = -p_prev[0] + 2*(1-beta2)*p[0] + 2*beta2*p[1]
        p_next[-1] = -p_prev[-1] + 2*(1-beta2)*p[-1] + 2*beta2*p[-2]

        # permutation
        p_prev, p, p_next = p, p_next, p_prev

    end = time.time()
    return end - start

# --- valeurs testées ---
Ns = np.linspace(201, 3001, 10)
Ns = Ns.astype(int)  

# mesures
times = []

for N in Ns:
    t = run_sim(N)
    times.append(t)

times = np.array(times)

# --- tracé uniquement en fonction de l'espace ---
plt.figure(figsize=(8,5))
plt.plot(Ns, times, marker='o', label="Temps de calcul vs taille du maillage")
plt.xlabel("N (taille du maillage)")
plt.ylabel("Temps de calcul (s)")
plt.title("Complexité expérimentale en fonction de l'espace")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()
