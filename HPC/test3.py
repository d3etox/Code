import numpy as np
import time

def run_sim(N):
    # --- paramètres physiques ---
    c = 1.0
    L = 1.0
    Tfinal = 1 # petit temps pour des timings rapides

    N = 1001                   # nombre de points (N > 1)
    dx = L / (N - 1)
    beta = 0.9                  # nombre de Courant: beta = c dt / dx < 1
    dt = beta * dx / c
    nt = int(Tfinal / dt + 0.5)  # arrondi à l'entier le plus proche
    beta2 = (c * dt / dx)**2

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
    return nt, end - start

# --- simulation unique ---
N = 1001
nt, t = run_sim(N)
print(f"N = {N:5d}, nt = {nt:6d}, temps = {t:10.6f}")
