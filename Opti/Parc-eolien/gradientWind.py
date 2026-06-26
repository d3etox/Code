import numpy as np
import matplotlib.pyplot as plt
import windfarm_eval
import windValid

# -----------------------------
# Paramètres
# -----------------------------
n_max = 4
eps_values = [1, 5, 10]
param_file = "instances/1/param.txt"
x0p_file = "instances/1/x0p_generated.txt"
x0m_file = "instances/1/x0m_generated.txt"

# -----------------------------
# Fonction gradient coordonnée
# -----------------------------
def finite_diff_grad(turbSimu, eps):
    grad = np.zeros(2*n_max)
    for i in range(n_max):
        for j in range(2):  # x ou y
            turbxPlus = turbSimu.copy()
            turbxMoins = turbSimu.copy()

            turbxPlus[i] = (turbxPlus[i][0] + (1-j)*eps,
                            turbxPlus[i][1] + j*eps)
            turbxMoins[i] = (turbxMoins[i][0] - (1-j)*eps,
                             turbxMoins[i][1] - j*eps)

            coordsPlus = [c for t in turbxPlus for c in t]
            coordsMoins = [c for t in turbxMoins for c in t]

            with open(x0p_file, "w") as fp:
                fp.write(str(coordsPlus))
            with open(x0m_file, "w") as fm:
                fm.write(str(coordsMoins))

            f_plus, _, _ = windfarm_eval.windfarm_eval(param_file, x0p_file)
            f_moins, _, _ = windfarm_eval.windfarm_eval(param_file, x0m_file)

            grad[2*i+j] = (f_plus - f_moins) / (2*eps)
    return grad

# -----------------------------
# Calcul pour eps = 1, 5, 10
# -----------------------------
turbines = windValid.valid(n_max)
grad_eps = {eps: finite_diff_grad(turbines, eps) for eps in eps_values}

# -----------------------------
# Affichage tableau en console
# -----------------------------
print("\n📊 Gradient estimé pour chaque coordonnée (selon eps):\n")
print(f"{'Coord':<8} {'eps=1':>15} {'eps=5':>15} {'eps=10':>15}")
print("-"*60)
for idx in range(2*n_max):
    values = [grad_eps[eps][idx] for eps in eps_values]
    print(f"{idx+1:<8} {values[0]:>15.6f} {values[1]:>15.6f} {values[2]:>15.6f}")

# -----------------------------
# Plot des 8 gradients
# -----------------------------
fig, axes = plt.subplots(2, 4, figsize=(18, 8), sharex=True)
axes = axes.flatten()

for idx in range(2*n_max):
    grads = [grad_eps[eps][idx] for eps in eps_values]
    axes[idx].plot(eps_values, grads, marker="o", linestyle="-")
    axes[idx].set_title(f"Coord {idx+1}")
    axes[idx].set_xlabel("eps")
    axes[idx].set_ylabel("∂EAP/∂coord")

plt.tight_layout()
plt.show()
