import deepxde as dde
import numpy as np
import matplotlib.pyplot as plt

# =========================
# Paramètres du problème
# =========================
pi = np.pi
m = 2
omega = m * pi
T_final = 2*pi/omega


# =========================
# Domaine espace-temps
# =========================
geom = dde.geometry.Interval(0.0, 1.0)
timedomain = dde.geometry.TimeDomain(0.0, T_final)
geomtime = dde.geometry.GeometryXTime(geom, timedomain)

# =========================
# Équations de Maxwell 1D
# =========================
def maxwell_pde(x, y):
    """
    y[:,0] = E(x,t)
    y[:,1] = H(x,t)
    """
    E = y[:, 0:1]
    H = y[:, 1:2]

    dE_t = dde.grad.jacobian(y, x, i=0, j=1)
    dE_x = dde.grad.jacobian(y, x, i=0, j=0)

    dH_t = dde.grad.jacobian(y, x, i=1, j=1)
    dH_x = dde.grad.jacobian(y, x, i=1, j=0)

    eq1 = dE_t - dH_x
    eq2 = dH_t - dE_x

    return [eq1, eq2]

# =========================
# Conditions aux limites
# =========================
def boundary_left(x, on_boundary):
    return on_boundary and np.isclose(x[0], 0.0)

def boundary_right(x, on_boundary):
    return on_boundary and np.isclose(x[0], 1.0)

bc_left = dde.icbc.DirichletBC(
    geomtime,
    lambda x: 0.0,
    boundary_left,
    component=0,   # E
)

bc_right = dde.icbc.DirichletBC(
    geomtime,
    lambda x: 0.0,
    boundary_right,
    component=0,   # E
)

# =========================
# Conditions initiales
# =========================
def E_init(x):
    return np.sin(omega * x[:, 0:1])

def H_init(x):
    return np.zeros((len(x), 1))

ic_E = dde.icbc.IC(
    geomtime,
    E_init,
    lambda x, on_initial: on_initial,
    component=0,
)

ic_H = dde.icbc.IC(
    geomtime,
    H_init,
    lambda x, on_initial: on_initial,
    component=1,
)

# =========================
# Données PINN
# =========================
data = dde.data.TimePDE(
    geomtime,
    maxwell_pde,
    [bc_left, bc_right, ic_E, ic_H],
    num_domain=4000,
    num_boundary=200,
    num_initial=200,
)

# =========================
# Réseau de neurones
# =========================
net = dde.nn.FNN(
    [2] + [64] * 4 + [2],
    "tanh",
    "Glorot normal",
)

model = dde.Model(data, net)

# =========================
# Entraînement
# =========================
model.compile(
    "adam",
    lr=1e-3,
)

losshistory, train_state = model.train(epochs=4000)

# Affinage L-BFGS (recommandé)
model.compile("L-BFGS")
model.train()

# =========================
# Solution analytique (validation)
# =========================
def exact_solution(x):
    E = np.sin(omega * x[:, 0:1]) * np.cos(omega * x[:, 1:2])
    H = np.cos(omega * x[:, 0:1]) * np.sin(omega * x[:, 1:2])
    return np.hstack((E, H))

# =========================
# Erreur L2
# =========================
X_test = geomtime.random_points(2000)
y_pred = model.predict(X_test)
y_exact = exact_solution(X_test)

error_E = dde.metrics.l2_relative_error(y_exact[:, 0], y_pred[:, 0])
error_H = dde.metrics.l2_relative_error(y_exact[:, 1], y_pred[:, 1])


print(f"Erreur relative L2 sur E : {error_E:.3e}")
print(f"Erreur relative L2 sur H : {error_H:.3e}")


# =========================
# Paramètres physiques
# =========================
L = 1.0      # longueur de la cavité en mètres
c0 = 3e8     # vitesse de la lumière en m/s

# =========================
# Conversion adim -> physique
# =========================
x_plot = np.linspace(0, 1, 200)       # coordonnées adimensionnées
x_phys = x_plot * L                    # coordonnées physiques

t_plot = [0.0, 0.25, 0.5, 0.75, 1.0] # temps adimensionnés
t_phys_plot = np.array(t_plot) / c0   # temps en secondes

# =========================
# Plot 1 : Champ E(x,t) à différents instants
# =========================

# Plot E(x,t) à différents instants
plt.figure(figsize=(10, 5))
for t_adim, t_phys in zip(t_plot, t_phys_plot):
    X = np.hstack((x_plot.reshape(-1,1), np.full((len(x_plot),1), t_adim)))
    y = model.predict(X)
    plt.plot(x_phys, y[:,0], label=f"E, t={t_phys:.2e} s")
plt.xlabel("x (m)")
plt.ylabel("E(x,t)")
plt.title("Champ E(x,t) à différents instants physiques")
plt.legend()
plt.show()

# Plot H(x,t) à différents instants
plt.figure(figsize=(10, 5))
for t_adim, t_phys in zip(t_plot, t_phys_plot):
    X = np.hstack((x_plot.reshape(-1,1), np.full((len(x_plot),1), t_adim)))
    y = model.predict(X)
    plt.plot(x_phys, y[:,1], label=f"H, t={t_phys:.2e} s")
plt.xlabel("x (m)")
plt.ylabel("H(x,t)")
plt.title("Champ H(x,t) à différents instants physiques")
plt.legend()
plt.show()

# =========================
# Plot 2 : Convergence du PINN
# =========================

# Plot de l'erreur L2 relative (E et H) selon l'epoch


# Correction : utiliser les steps de losshistory et calculer l'erreur L2 sur les checkpoints d'entraînement
epochs = np.array(losshistory.steps)
error_E_epoch = []
error_H_epoch = []
error_sum_epoch = []
if hasattr(losshistory, 'steps') and hasattr(losshistory, 'loss_train'):
    for _ in epochs:
        y_pred = net.predict(data.train_x)
        y_exact = exact_solution(data.train_x)
        err_E = dde.metrics.l2_relative_error(y_exact[:, 0], y_pred[:, 0])
        err_H = dde.metrics.l2_relative_error(y_exact[:, 1], y_pred[:, 1])
        error_E_epoch.append(err_E)
        error_H_epoch.append(err_H)
        error_sum_epoch.append(err_E + err_H)

    plt.figure(figsize=(10, 5))
    plt.plot(epochs, error_E_epoch, label="Erreur L2 relative E")
    plt.plot(epochs, error_H_epoch, label="Erreur L2 relative H")
    plt.plot(epochs, error_sum_epoch, label="Erreur L2 additionnée (E+H)")
    plt.xlabel("Epoch")
    plt.ylabel("Erreur L2 relative")
    plt.title("Erreur L2 relative par epoch (E, H, E+H)")
    plt.legend()
    plt.show()
