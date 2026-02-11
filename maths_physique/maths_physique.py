import numpy as np
import matplotlib.pyplot as plt
import time as time

def u0(x):
    return np.where(x <= 0.25, 1.0, 0.0)

def s1(u, a, dx, dt):
    #schéma 1 :
    u = u - a*dt/dx*(np.roll(u, 1) - u)
    return u

def s2(u, a, dx, dt):
    #schéma 2 :
    u = u - a*dt/(2*dx)*(np.roll(u, 1) - np.roll(u, -1))
    return u

def s3(u, a, dx, dt):
    #schéma 3:
    u = u - a*dt/dx*(u- np.roll(u, -1))
    return u

def s4(u, a, dx, dt):
    #schéma 4 (Lax-Friedrichs) :
    u = 0.5*(np.roll(u, -1) + np.roll(u, 1)) - a*dt/(2*dx)*(np.roll(u, 1) - np.roll(u, -1))
    return u

def s5(u, a, dx, dt):
    #schéma 5 (Lax-Wendroff) :
    u = u - a*dt/(2*dx)*(np.roll(u, 1) - np.roll(u, -1)) + (a*dt/dx)**2/2*(np.roll(u, 1) - 2*u + np.roll(u, -1))
    return u

def s6(u, a, dx, dt):
    #schéma 6 (Upwind amélioré) :
    u = u - a*dt/dx*(u - np.roll(u, 1)) + (a*dt/dx)**2/2*(u - 2*np.roll(u, 1) + np.roll(u, 2))
    return u

N = 4000
L = 1.0 
dx = L/N
T = 2.0
a = 0.25
x = np.linspace(0, L, N, endpoint=False)
dt = 1e-3
print("Paramètres : N={}, dx={}, dt={}, T={}, a={}".format(N, dx, dt, T, a))
print("Condition CFL s1 et s3 : a*dt/dx={} <= 1 pour stabilité".format(a*dt/dx))
print("Condition CFL s4 : (a*dt/dx)**2)={} <= 3 pour stabilité".format((a*dt/dx)**2))

def solveur(scheme, u0, a, dx, dt, T):
    u = u0(x)
    nsteps = int(T/dt)
    for _ in range(nsteps):
        u = scheme(u, a, dx, dt)
    return u
t0 = time.time()
u1 = solveur(s1, u0, a, dx, dt, T)
t1 = time.time()
print("Temps de calcul schéma 1 : {:.4f} secondes".format(t1 - t0))
u2 = solveur(s2, u0, a, dx, dt, T)
t2 = time.time()
print("Temps de calcul schéma 2 : {:.4f} secondes".format(t2 - t1))
u3 = solveur(s3, u0, a, dx, dt, T)
t3 = time.time()
print("Temps de calcul schéma 3 : {:.4f} secondes".format(t3 - t2))
u4 = solveur(s4, u0, a, dx, dt, T)
t4 = time.time()
print("Temps de calcul schéma 4 : {:.4f} secondes".format(t4 - t3))
u5 = solveur(s5, u0, a, dx, dt, T)
t5 = time.time()
print("Temps de calcul schéma 5 : {:.4f} secondes".format(t5 - t4))
u6 = solveur(s6, u0, a, dx, dt, T)
t6 = time.time()
print("Temps de calcul schéma 6 : {:.4f} secondes".format(t6 - t5)) 

u_exacte = u0((x - a*T) % L)

fig, axes = plt.subplots(6, 1, figsize=(8, 10))  # 3 lignes, 1 colonne

axes[0].plot(x, u1, label='Schéma 1', linestyle='-')
axes[0].plot(x, u_exacte, label='Solution exacte', color='black')
axes[0].set_title("Schéma 1")
axes[0].legend()

axes[1].plot(x, u2, label='Schéma 2', linestyle='-')
axes[1].plot(x, u_exacte, label='Solution exacte', color='black')
axes[1].set_title("Schéma 2")
axes[1].legend()

axes[2].plot(x, u3, label='Schéma 3', linestyle='-')
axes[2].plot(x, u_exacte, label='Solution exacte', color='black')
axes[2].set_title("Schéma 3")
axes[2].legend()

axes[3].plot(x, u4, label='Schéma 4 (Lax-Friedrichs)', linestyle='-')
axes[3].plot(x, u_exacte, label='Solution exacte', color='black')
axes[3].set_title("Schéma 4 (Lax-Friedrichs)")
axes[3].legend()

axes[4].plot(x, u5, label='Schéma 5 (Lax-Wendroff)', linestyle='-')
axes[4].plot(x, u_exacte, label='Solution exacte', color='black')
axes[4].set_title("Schéma 5 (Lax-Wendroff)")
axes[4].legend()

axes[5].plot(x, u6, label='Schéma 6 (Upwind amélioré)', linestyle='-')
axes[5].plot(x, u_exacte, label='Solution exacte', color='black')
axes[5].set_title("Schéma 6 (Upwind amélioré)")
axes[5].legend()    

for ax in axes:
    ax.set_xlabel("x")
    ax.set_ylabel("u(x,T)")

plt.tight_layout()
plt.show()