import random
import math
import windfarm_eval
import matplotlib.pyplot as plt

# -----------------------------
# Fonction pour ajouter 1 turbine
# -----------------------------
def add_turbine(turbines, x_min, x_max, y_min, y_max, min_dist, max_attempts=5000):
    attempts = 0
    while attempts < max_attempts:
        x = random.uniform(x_min, x_max)
        y = random.uniform(y_min, y_max)

        too_close = any(math.hypot(tx - x, ty - y) < min_dist for tx, ty in turbines)

        if not too_close:
            return turbines + [(x, y)]  # nouvelle liste avec la turbine en plus

        attempts += 1

    return None  # échec

# -----------------------------
# Paramètres
# -----------------------------
#n_max = 20
x_min, x_max = 104500, 107500
y_min, y_max = 1043500, 1045600
diameter = 80
min_dist = 2 * diameter

x0_file = 'instances/1/x0_generated.txt'
param_file = 'instances/1/param.txt'

# -----------------------------
# Process incrémental avec validation
# -----------------------------

def valid(n_max):
	turbines = []
	n = 1
	while n <= n_max:
		new_turbines = add_turbine(turbines, x_min, x_max, y_min, y_max, min_dist)

		if new_turbines is None:
			print(f"⚠️ Impossible de placer une {n}ᵉ turbine après trop d’essais. On arrête.")
			break

		coords = [coord for turbine in new_turbines for coord in turbine]

		with open(x0_file, 'w') as f:
			f.write(str(coords))

		eap, spacing_constraint, placing_constraint = windfarm_eval.windfarm_eval(param_file, x0_file)

		if spacing_constraint == 0 and placing_constraint == 0:
			turbines = new_turbines
			print(f"✅ Turbine {n} placée. EAP = {eap}")
			n += 1
		else:
			print(f"❌ Turbine {n} non valide → spacing_constraint = {spacing_constraint}, placing_constraint = {placing_constraint}")
	return turbines