import deepxde as dde
import numpy as np
import matplotlib.pyplot as plt

# =========================
# Paramètres du problème
# =========================
pi = np.pi
m = 2
omega = m * pi
T_final = 2 * pi / omega

geom = dde.geometry.Interval(0.0, 1.0)
timedomain = dde.geometry.TimeDomain(0.0, T_final)
geomtime = dde.geometry.GeometryXTime(geom, timedomain)

def maxwell_pde(x, y):
    E = y[:, 0:1]
    H = y[:, 1:2]
    dE_t = dde.grad.jacobian(y, x, i=0, j=1)
    dE_x = dde.grad.jacobian(y, x, i=0, j=0)
    dH_t = dde.grad.jacobian(y, x, i=1, j=1)
    dH_x = dde.grad.jacobian(y, x, i=1, j=0)
    eq1 = dE_t - dH_x
    eq2 = dH_t - dE_x
    return [eq1, eq2]

def boundary_left(x, on_boundary):
    return on_boundary and np.isclose(x[0], 0.0)
def boundary_right(x, on_boundary):
    return on_boundary and np.isclose(x[0], 1.0)

def E_init(x):
    return np.sin(omega * x[:, 0:1])
def H_init(x):
    return np.zeros((len(x), 1))

def exact_solution(x):
    E = np.sin(omega * x[:, 0:1]) * np.cos(omega * x[:, 1:2])
    H = np.cos(omega * x[:, 0:1]) * np.sin(omega * x[:, 1:2])
    return np.hstack((E, H))

layer_options = [1, 2, 3, 4, 5, 6]
neurons_options = [8, 16, 32, 64, 128]
results = []

for num_layers in layer_options:
    for num_neurons in neurons_options:
        print(f"Test: {num_layers} couches, {num_neurons} neurones/couche")
        bc_left = dde.icbc.DirichletBC(geomtime, lambda x: 0.0, boundary_left, component=0)
        bc_right = dde.icbc.DirichletBC(geomtime, lambda x: 0.0, boundary_right, component=0)
        ic_E = dde.icbc.IC(geomtime, E_init, lambda x, on_initial: on_initial, component=0)
        ic_H = dde.icbc.IC(geomtime, H_init, lambda x, on_initial: on_initial, component=1)
        data = dde.data.TimePDE(
            geomtime, maxwell_pde, [bc_left, bc_right, ic_E, ic_H],
            num_domain=1000, num_boundary=100, num_initial=100
        )
        net = dde.nn.FNN([2] + [num_neurons] * num_layers + [2], "tanh", "Glorot normal")
        model = dde.Model(data, net)
        model.compile("adam", lr=1e-3)
        losshistory, train_state = model.train(epochs=1000, display_every=100)
        model.compile("L-BFGS")
        model.train()
        X_test = geomtime.random_points(500)
        y_pred = model.predict(X_test)
        y_exact = exact_solution(X_test)
        error_E = dde.metrics.l2_relative_error(y_exact[:, 0], y_pred[:, 0])
        error_H = dde.metrics.l2_relative_error(y_exact[:, 1], y_pred[:, 1])
        results.append({
            "layers": num_layers,
            "neurons": num_neurons,
            "error_E": error_E,
            "error_H": error_H
        })
        print(f"Erreur E: {error_E:.3e}, Erreur H: {error_H:.3e}")

# Affichage des résultats
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

df = pd.DataFrame(results)
plt.figure(figsize=(10, 6))
sns.heatmap(df.pivot("layers", "neurons", "error_E"), annot=True, fmt=".1e", cmap="viridis")
plt.title("Erreur L2 relative E selon le nombre de couches et de neurones")
plt.xlabel("Neurones par couche")
plt.ylabel("Nombre de couches")
plt.show()

plt.figure(figsize=(10, 6))
sns.heatmap(df.pivot("layers", "neurons", "error_H"), annot=True, fmt=".1e", cmap="viridis")
plt.title("Erreur L2 relative H selon le nombre de couches et de neurones")
plt.xlabel("Neurones par couche")
plt.ylabel("Nombre de couches")
plt.show()
