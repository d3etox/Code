import random
import math
import windfarm_eval
import matplotlib.pyplot as plt

# -----------------------------
# PARAMÈTRES
# -----------------------------
POP_SIZE = 40
N_TURBINES = 10
X_MIN, X_MAX = 104500, 107500
Y_MIN, Y_MAX = 1043500, 1045600
DIAMETER = 80
MUT_PROB = 0.2
CROSS_PROB = 0.5
GENERATIONS = 100
min_dist = 2 * DIAMETER
param_file = "instances/1/param.txt"
x0_file = "instances/1/x0_genetic.txt"

# Cache pour éviter les réévaluations coûteuses
fitness_cache = {}

# -----------------------------
# ÉVALUATION
# -----------------------------
def evaluate(ind):
    key = tuple(ind)
    if key in fitness_cache:
        return fitness_cache[key]

    with open(x0_file, "w") as f:
        f.write(str(ind))
    eap, spacing, placing = windfarm_eval.windfarm_eval(param_file, x0_file)

    # Fitness : pénalisation si contraintes violées
    if spacing != 0 or placing != 0:
        fitness = -(spacing + placing)
    else:
        fitness = eap

    fitness_cache[key] = fitness
    return fitness

# -----------------------------
# Fonction pour générer une turbine
# -----------------------------
def add_turbine(turbines):
    max_attempts = 5000
    attempts = 0
    while attempts < max_attempts:
        x = random.uniform(X_MIN, X_MAX)
        y = random.uniform(Y_MIN, Y_MAX)
        too_close = any(math.hypot(tx - x, ty - y) < min_dist for tx, ty in turbines)
        if not too_close:
            return turbines + [(x, y)]
        attempts += 1
    return None

# -----------------------------
# GÉNÉRATION D’UN INDIVIDU VALIDE
# -----------------------------
def generate_individual():
    turbines = []
    n = 1
    while n <= N_TURBINES:
        new_turbines = add_turbine(turbines)
        if new_turbines is None:
            print("⚠️ Impossible de placer la turbine, on recommence tout l'individu")
            turbines = []
            n = 1
            continue

        coords = [coord for t in new_turbines for coord in t]

        with open(x0_file, "w") as f:
            f.write(str(coords))

        try:
            eap, spacing, placing = windfarm_eval.windfarm_eval(param_file, x0_file)
        except Exception as ex:
            print("Erreur windfarm_eval:", ex)
            turbines = []
            n = 1
            continue

        if spacing == 0 and placing == 0:
            turbines = new_turbines
            print(f"✅ Turbine {n} placée. EAP intermédiaire = {eap:.2f}")
            n += 1
        else:
            print(f"❌ Turbine {n} invalide (spacing={spacing}, placing={placing}). Regénération...")

    print("🎯 Individu valide généré avec EAP final =", eap)
    return [float(c) for t in turbines for c in t]

# -----------------------------
# SÉLECTION : tournoi
# -----------------------------
def tournament_selection(pop, k=3):
    best = random.choice(pop)
    for _ in range(k - 1):
        challenger = random.choice(pop)
        if evaluate(challenger) > evaluate(best):
            best = challenger
    return best

# -----------------------------
# CROISEMENT : moyenne pondérée
# -----------------------------
def crossover(p1, p2):
    if random.random() < CROSS_PROB:
        alpha = random.random()
        return [alpha * a + (1 - alpha) * b for a, b in zip(p1, p2)]
    return p1.copy()

# -----------------------------
# MUTATION : perturbation gaussienne
# -----------------------------
def mutate(ind, sigma=500):
    for i in range(len(ind)):
        if random.random() < MUT_PROB:
            ind[i] += random.gauss(0, sigma)
            if i % 2 == 0:
                ind[i] = min(max(ind[i], X_MIN), X_MAX)
            else:
                ind[i] = min(max(ind[i], Y_MIN), Y_MAX)
    return ind

# -----------------------------
# INITIALISATION
# -----------------------------
population = []
for i in range(POP_SIZE):
    print("\n" + "-"*40)
    print(f"🧬 Génération de l’individu {i+1}/{POP_SIZE}")
    print("-"*40)
    ind = generate_individual()
    population.append(ind)

# -----------------------------
# BOUCLE PRINCIPALE
# -----------------------------
fitness_history = []

for gen in range(GENERATIONS):
    new_pop = []
    for _ in range(POP_SIZE):
        # Sélection
        p1 = tournament_selection(population)
        p2 = tournament_selection(population)

        # Croisement + mutation
        child = mutate(crossover(p1, p2))

        new_pop.append(child)
    
    # Sélection + croisement + mutation déjà faits, new_pop rempli
    evaluated_children = [(ind, evaluate(ind)) for ind in new_pop]

    # On garde 1 élite de la population précédente
    old_evaluated = [(ind, evaluate(ind)) for ind in population]
    elite = max(old_evaluated, key=lambda x: x[1])

    # Fusionner enfants + élite et garder les meilleurs POP_SIZE
    combined = evaluated_children + [elite]
    combined_sorted = sorted(combined, key=lambda x: x[1], reverse=True)
    population = [ind for ind, f in combined_sorted[:POP_SIZE]]

    # Meilleur individu et fitness
    best, best_fit = combined_sorted[0]
    fitness_history.append(best_fit)
    print(f"Génération {gen+1}: fitness = {best_fit:.2f}")


# -----------------------------
# RÉSULTAT FINAL
# -----------------------------
print("✅ Meilleur individu :", best)
print("✅ Fitness finale :", evaluate(best))

# -----------------------------
# PLOT FINAL
# -----------------------------
best_coords = [(best[i], best[i + 1]) for i in range(0, len(best), 2)]
x_coords, y_coords = zip(*best_coords)

fig, ax = plt.subplots(figsize=(8, 6))
ax.scatter(x_coords, y_coords, color="blue", label="Turbines")

for (x, y) in best_coords:
    circle = plt.Circle((x, y), DIAMETER / 2, fill=False, linestyle="--", color="gray")
    ax.add_patch(circle)

ax.set_xlim(X_MIN, X_MAX)
ax.set_ylim(Y_MIN, Y_MAX)
ax.set_aspect("equal", adjustable="box")
ax.set_xlabel("X (m)")
ax.set_ylabel("Y (m)")
ax.set_title("Disposition optimale des turbines (selfmade GA)")
ax.legend()
plt.grid(True)
plt.show()

# -----------------------------
# PLOT DE CONVERGENCE
# -----------------------------
plt.figure(figsize=(7, 4))
plt.plot(fitness_history, color="darkred")
plt.title("Évolution de la meilleure fitness")
plt.xlabel("Génération")
plt.ylabel("Fitness (EAP)")
plt.grid(True)
plt.show()
