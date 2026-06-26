import random
import math
import windfarm_eval
import matplotlib.pyplot as plt
import numpy as np

# -----------------------------
# Paramètres
# -----------------------------
n_turbines = 10
x_min, x_max = 104500, 107500
y_min, y_max = 1043500, 1045600
diameter = 80
min_dist = 2 * diameter
eps = 10             # pas pour gradient
alpha = 200          # learning rate
max_iter = 50        # nombre max d'itérations

param_file = "instances/1/param.txt"
x0_file = "instances/1/x0_gradient.txt"
x0p_file = "instances/1/x0p_generated.txt"
x0m_file = "instances/1/x0m_generated.txt"

# -----------------------------
# Fonction pour générer une turbine valide
# -----------------------------
def add_turbine(turbines):
    max_attempts = 5000
    for _ in range(max_attempts):
        x = random.uniform(x_min, x_max)
        y = random.uniform(y_min, y_max)
        if all(math.hypot(tx - x, ty - y) >= min_dist for tx, ty in turbines):
            return turbines + [(x, y)]
    return None

# -----------------------------
# Génération d'une configuration initiale
# -----------------------------
def generate_initial_turbines():
    turbines = []
    for n in range(n_turbines):
        turbines = add_turbine(turbines)
        if turbines is None:
            return None
    return turbines

# -----------------------------
# Gradient par différences finies
# -----------------------------
def compute_gradient(turbines, eps=eps):
    grad = np.zeros(2 * n_turbines)
    for i in range(n_turbines):
        for j in range(2):
            turb_plus = turbines.copy()
            turb_minus = turbines.copy()

            if j == 0:  # x
                turb_plus[i] = (turb_plus[i][0] + eps, turb_plus[i][1])
                turb_minus[i] = (turb_minus[i][0] - eps, turb_minus[i][1])
            else:       # y
                turb_plus[i] = (turb_plus[i][0], turb_plus[i][1] + eps)
                turb_minus[i] = (turb_minus[i][0], turb_minus[i][1] - eps)

            coords_plus = [c for t in turb_plus for c in t]
            coords_minus = [c for t in turb_minus for c in t]

            # Écrire fichiers temporaires pour windfarm_eval
            with open(x0p_file, "w") as f:
                f.write(str(coords_plus))
            with open(x0m_file, "w") as f:
                f.write(str(coords_minus))

            eap_plus, a_plus, b_plus = windfarm_eval.windfarm_eval(param_file, x0p_file)
            eap_minus, a_minus, b_minus = windfarm_eval.windfarm_eval(param_file, x0m_file)

            grad[2*i + j] = (eap_plus+a_plus+b_plus - (eap_minus+a_minus+b_minus)) / (2*eps)
    return grad

# -----------------------------
# Gradient ascent
# -----------------------------
turbines = None
while turbines is None:
    turbines = generate_initial_turbines()

# Sauvegarde de la configuration initiale
initial_turbines = np.array(turbines, dtype=float)

# EAP initial
coords_init = [c for t in initial_turbines for c in t]
with open(x0_file, "w") as f:
    f.write(str(coords_init))
eap_init, spacing_init, placing_init = windfarm_eval.windfarm_eval(param_file, x0_file)
print(f"Initial: EAP = {eap_init}, spacing={spacing_init}, placing={placing_init}")

# Optimisation
turbines = np.copy(initial_turbines)
for iteration in range(max_iter):
    grad = compute_gradient(turbines.tolist())
    turbines += alpha * grad.reshape(-1,2)

    # Clip pour rester dans le terrain
    turbines[:,0] = np.clip(turbines[:,0], x_min, x_max)
    turbines[:,1] = np.clip(turbines[:,1], y_min, y_max)

    coords = [c for t in turbines for c in t]
    with open(x0_file, "w") as f:
        f.write(str(coords))
    eap, spacing, placing = windfarm_eval.windfarm_eval(param_file, x0_file)
    print(f"Iteration {iteration+1}: EAP = {eap}, spacing={spacing}, placing={placing}")

# -----------------------------
# Plot comparatif avant/après
# -----------------------------
x_init, y_init = initial_turbines[:,0], initial_turbines[:,1]
x_final, y_final = turbines[:,0], turbines[:,1]

fig, ax = plt.subplots(figsize=(8,6))

# Initial (en rouge)
ax.scatter(x_init, y_init, color="red", label="Initial")
for (x, y) in zip(x_init, y_init):
    circle = plt.Circle((x, y), diameter/2, fill=False, linestyle="--", color="red", alpha=0.4)
    ax.add_patch(circle)

# Final (en bleu)
ax.scatter(x_final, y_final, color="blue", label="Optimisé")
for (x, y) in zip(x_final, y_final):
    circle = plt.Circle((x, y), diameter/2, fill=False, linestyle="--", color="blue", alpha=0.4)
    ax.add_patch(circle)

ax.set_xlim(x_min, x_max)
ax.set_ylim(y_min, y_max)
ax.set_aspect("equal")
ax.legend()
ax.set_title(f"Disposition des turbines\nAvant (EAP={eap_init:.2f}) vs Après (EAP={eap:.2f})")
plt.show()
