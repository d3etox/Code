import deepxde as dde
import numpy as np
import matplotlib.pyplot as plt
import tensorflow as tf

# =========================
# Paramètres du problème
# =========================
pi = np.pi
T_final = 5

# =========================
# Domaine espace-temps
# =========================
geom = dde.geometry.Interval(0.0, 6.0)
timedomain = dde.geometry.TimeDomain(0.0, T_final)
geomtime = dde.geometry.GeometryXTime(geom, timedomain)

# =========================
# Équations de Maxwell 1D
# =========================
def epsr(x):
    return tf.where(x[:, 0:1] <= 4.0, 1.0, 2.0)

def maxwell_pde(x, y):
    E = y[:, 0:1]
    H = y[:, 1:2]
    dE_t = dde.grad.jacobian(y, x, i=0, j=1)
    dE_x = dde.grad.jacobian(y, x, i=0, j=0)
    dH_t = dde.grad.jacobian(y, x, i=1, j=1)
    dH_x = dde.grad.jacobian(y, x, i=1, j=0)
    epsr_val = epsr(x)
    eq1 = dE_t - dH_x
    eq2 = dH_t - dE_x
    return [eq1, eq2]

# =========================
# Conditions aux limites personnalisées
# =========================
def boundary_left(x, on_boundary):
    return on_boundary and np.isclose(x[0], 0.0)

def boundary_right(x, on_boundary):
    return on_boundary and np.isclose(x[0], 6.0)

bc_left_H = dde.icbc.OperatorBC(
    geomtime,
    lambda x, y, X: y[0] - y[1],
    boundary_left,
)

bc_right_H = dde.icbc.OperatorBC(
    geomtime,
    lambda x, y, X: y[0] + y[1],
    boundary_right,
)

# =========================
# Conditions initiales
# =========================
alpha = 10.0
xg = 3.0

def E_init(x):
    return np.exp(-alpha * (x[:, 0:1] - xg)**2)

def H_init(x):
    return -np.exp(-alpha * (x[:, 0:1] - xg)**2)

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
    [bc_left_H, bc_right_H, ic_E, ic_H],
    num_domain=4000,
    num_boundary=200,
    num_initial=200,
)

# =========================
# Réseau de neurones
# =========================
net = dde.nn.FNN(
    [2] + [64]*4 + [2],
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
# Plots et diagnostics
# =========================
c0 = 3e8  # m/s

if losshistory is not None:
    train_loss = np.array(losshistory.loss_train)
    plt.figure()
    if train_loss.ndim > 1 and train_loss.shape[1] >= 2:
        plt.semilogy(losshistory.steps, train_loss[:, 0], label="Train PDE E")
        plt.semilogy(losshistory.steps, train_loss[:, 1], label="Train PDE H")
    if losshistory.loss_test:
        test_loss = np.array(losshistory.loss_test)
        if test_loss.ndim > 1 and test_loss.shape[1] >= 2:
            plt.semilogy(losshistory.steps, test_loss[:, 0], '--', label="Test PDE E")
            plt.semilogy(losshistory.steps, test_loss[:, 1], '--', label="Test PDE H")
    plt.xlabel("Step")
    plt.ylabel("Loss")
    plt.title("Convergence de l'entraînement : PDE E et H")
    plt.legend()
    plt.show()

x_plot = np.linspace(0, 6, 300)
t_snap = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0]
print("\nAffichage de l'onde progressive E(x, t) à différents temps physiques:")

# Animation de l'onde progressive E(x, t) sur t ∈ [0, T_final]
import matplotlib.animation as animation

T_anim = T_final
frames = 500
t_vals = np.linspace(0, T_anim, frames)

fig, ax = plt.subplots()
line, = ax.plot([], [], lw=2)
ax.set_xlim(0, 6)
ax.set_ylim(-1.2, 1.2)
ax.set_xlabel("x (m)")
ax.set_ylabel("E(x, t)")
ax.set_title("Animation de l'onde progressive E(x, t)")

def init():
    line.set_data([], [])
    return line,

def animate(i):
    t = t_vals[i]
    X = np.hstack((x_plot.reshape(-1, 1), np.full((len(x_plot), 1), t)))
    y = model.predict(X)
    line.set_data(x_plot, y[:, 0])
    ax.set_title(f"Onde progressive : E(x, t), t = {t:.2f}")
    return line,

ani = animation.FuncAnimation(fig, animate, frames=frames, init_func=init, blit=True, interval=50)
plt.show()

# Pour sauvegarder la vidéo (optionnel, nécessite ffmpeg ou imagemagick)
ani.save('machine_learning/onde_progressive_custom_bc.mp4', writer='ffmpeg', fps=60)
