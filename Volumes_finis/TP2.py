import numpy as np
import matplotlib.pyplot as plt

def f(x, VT, mu):
    num = VT * x**2
    den = ((1 - x)**2)/mu + x**2
    return num / den

def F(S, VT, mu):
    N = len(S)
    Flux = np.zeros(N)
    Flux[0] = f(1, VT, mu) - f(S[0], VT, mu)
    for i in range(1, N):
        Flux[i] = f(S[i-1], VT, mu) - f(S[i], VT, mu)
    return Flux

def VF(VT, mu, L, N, tf):
    h = L/tf
    X = np.linspace(h/2, L-h/2, N)
    Ah = np.zeros((N,N))
    Sh = F(X, VT, mu)*h

    for i in range(N):
        if i == 0:
            Ah[i,i] = 1/h
        else :
            Ah[i,i-1] = -1/h
            Ah[i,i] = 1/h

    Uh = np.linalg.solve(Ah,Sh)

    return X,Uh

tf = 3600*24*360*10 
VT = 10e-6 
mu = 10 
L = 1000
N = 1000

def solve_BL(VT, mu, L, N, tf, CFL=1.0):
    h = L / N

    # vitesse max = max |f'(s)| sur [0,1]
    s = np.linspace(0,1,1000)
    df = np.gradient(f(s,VT,mu), s)
    vmax = np.max(np.abs(df))

    dt = CFL * h / vmax

    nt = int(tf / dt) + 1

    S = np.zeros(N)   # condition initiale

    for n in range(nt):
        S = S - dt/h * F(S, VT, mu)

    x = np.linspace(h/2, L-h/2, N)
    return x, S, dt, nt

for N in [10, 100, 400]:
    x, S, dt, nt = solve_BL(VT, mu, L, N, tf, CFL=1)
    plt.plot(x, S, label=f"N={N}")

x, S, dt, nt = solve_BL(VT, mu, L, N, tf, CFL=1.5)

plt.show()