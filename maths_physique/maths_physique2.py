import numpy as np
import matplotlib.pyplot as plt

def u0(x):
    return np.where(x <= 0.25, 1.0, 0.0)

def s1(u, a, dx, dt):
    return u - a*dt/dx*(u - np.roll(u, 1))

def s2(u, a, dx, dt):
    return u - a*dt/(2*dx)*(np.roll(u, 1) - np.roll(u, -1))

def s3(u, a, dx, dt):
    return u - a*dt/dx*(u - np.roll(u, -1))

def s4(u, a, dx, dt):
    #schéma 4 (Lax-Friedrichs) :
    u = 0.5*(np.roll(u, -1) + np.roll(u, 1)) - a*dt/(2*dx)*(np.roll(u, 1) - np.roll(u, -1))
    return u

def s5(u, a, dx, dt):
    #schéma 5 (Lax-Wendroff) :
    u = u - a*dt/(2*dx)*(np.roll(u, 1) - np.roll(u, -1)) + (a*dt/dx)**2/2*(np.roll(u, 1) - 2*u + np.roll(u, -1))
    return u

def solveur(scheme, u0, a, dx, dt, T):
    x = np.linspace(0, 1, int(1/dx), endpoint=False)
    u = u0(x)
    nsteps = int(T/dt)
    for n in range(nsteps):
        u = scheme(u, a, dx, dt)
    return x, u

# Paramètres
T = 2.0
a = 0.25
dx_list = [1/50, 1/100, 1/200, 1/400]  # différentes résolutions
dt = 1e-4  # très petit pour éviter instabilité

# Calcul des erreurs L2
def L2_error(u_num, u_exact):
    return np.sqrt(np.sum((u_num - u_exact)**2)/len(u_num))

schemes = {'s1': s1, 's2': s2, 's3': s3, 's4': s4, 's5': s5}
errors = {name: [] for name in schemes}

for dx in dx_list:
    x_exact = np.linspace(0, 1, int(1/dx), endpoint=False)
    u_ex = u0((x_exact - a*T) % 1)
    for name, scheme in schemes.items():
        x_num, u_num = solveur(scheme, u0, a, dx, dt, T)
        err = L2_error(u_num, u_ex)
        print(f"Scheme {name}, dx={dx:.5f}, L2 error={err:.5e}")
        errors[name].append(err)

# Tracé de l'ordre de convergence
plt.figure(figsize=(8,6))
for name in schemes:
    plt.loglog(dx_list, errors[name], marker='o', label=name)

# Référence pente 1 et 2
dx_ref = np.array(dx_list)
plt.loglog(dx_ref, dx_ref, 'k--', label='pente 1')
plt.loglog(dx_ref, dx_ref**2, 'k-.', label='pente 2')

plt.gca().invert_xaxis()
plt.xlabel('dx')
plt.ylabel('L2 error')
plt.title("Ordre de convergence des schémas")
plt.legend()
plt.show()
