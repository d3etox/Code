import numpy as np
import matplotlib.pyplot as plt

# Paramètres du problème
K = 1
alpha = 2
t0 = 0
tf = 100
dt = 0.001
y0 = 0.5
C0 = (y0 - K) / y0
beta = 2
nu = 3
gamma = 1.5
F0 = 6
P0 = 10

# Fonction logistique
def logistic(y, t):
    return alpha * y * (K - y)

# Solution analytique pour la logistique
def f(t):
    t = np.asarray(t) # Assurer la compatibilité avec des vecteurs
    return K / (1 - C0 * np.exp(-alpha * K * t))

def periodique(X,P):
    return P,-X
def f1(X,P):
    return(X^2+P^2)
# Modèle de Volterra
def Volterra(F, P):
    dF = alpha * F - beta * F * P
    dP = -nu * P + gamma * F * P
    return dF, dP

# Méthode d'Euler pour une dimension
def Euler(f, y0, t0, tf, dt):
    n = int((tf - t0) / dt)
    Y = np.zeros(n + 1)
    Y[0] = y0
    t = t0
    for i in range(n):
        Y[i + 1] = Y[i] + dt * f(Y[i], t)
        t += dt
    return Y

# Méthode d'Euler pour deux dimensions
def Euler2D(f, y0, t0, tf, dt):
    n = int((tf - t0) / dt)
    Y = np.zeros((n + 1, 2)) # Matrice pour stocker F et P
    Y[0, :] = y0
    t = t0
    for i in range(n):
        dF, dP = f(Y[i, 0], Y[i, 1]) # Appliquer la fonction Volterra
        Y[i + 1, 0] = Y[i, 0] + dt * dF # Mise à jour de F
        Y[i + 1, 1] = Y[i, 1] + dt * dP # Mise à jour de P
        t += dt
    return Y

# Calcul des temps
time = [t0 + i * dt for i in range(int((tf - t0) / dt) + 1)]

# Résolution du système de Volterra avec Euler2D
result = Euler2D(periodique, [F0, P0], t0, tf, dt)

# Création des graphiques
fig, ax = plt.subplots(2, 1, figsize=(8, 12))

# Graphique 1 : Évolution des populations
ax[0].plot(time, result[:, 0], label="Proies (F)")
ax[0].plot(time, result[:, 1], label="Prédateurs (P)")
ax[0].set_xlabel('Temps')
ax[0].set_ylabel('Population')
ax[0].legend()
ax[0].set_title('Évolution des populations')

# Graphique 2 : Portrait de phase
ax[1].plot(result[:, 0], result[:, 1], color='orange')
ax[1].set_xlabel('Proies (F)')
ax[1].set_ylabel('Prédateurs (P)')
ax[1].set_title('Portrait de phase (Proies vs Prédateurs)')

plt.tight_layout()
plt.show()