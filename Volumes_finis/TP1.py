#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug  9 09:20:51 2024

@author: masssonr

resolution par la methode VF de l'equation elliptique 1D

- u''(x) = f(x) sur (0,L)
u(0) = uD
-u'(L) = g

"""

import numpy as np
import matplotlib.pyplot as plt

# longueur du domaine 
L=1.0 

# A COMPLETER 

def u(x):
    s = np.exp(np.sin(np.pi * x))   
    return s

def up(x): 
    s = np.pi * np.exp(np.sin(np.pi * x)) * np.cos(np.pi * x)
    return s 

def f(x):
    s = -np.pi**2 * np.exp(np.sin(np.pi * x)) * (np.cos(np.pi * x)**2 - np.sin(np.pi * x))
    return s


uD = u(0)
g = -up(L)

# calcul de X, Ah, Sh et Uh 

def VF(f,uD,g,L,N):

    h = L/N
    X = np.linspace(h/2, L-h/2, N)
    Ah = np.zeros((N,N))
    Sh = f(X)*h
    
    for i in range(N):
        if i == 0:
            Ah[i,i] = 3/h
            Ah[i,i+1] = -1/h
            Sh[i] += 2*uD/h
        elif i == N-1:
            Ah[i,i-1] = -1/h
            Ah[i,i] = 1/h
            Sh[i] += -g
        else:
            Ah[i,i-1] = -1/h
            Ah[i,i] = 2/h
            Ah[i,i+1] = -1/h
            
    Uh = np.linalg.solve(Ah,Sh)

    return X,Uh

########################




#nombre de mailles
N= 5
N= 10
N= 20


X,Uh = VF(f,uD,g,L,N)


#plot des solutions exactes et VF 
# plt.figure(1)
# plt.clf()
# Xfine = np.linspace(0,L,200)
# plt.plot(Xfine,u(Xfine), label="solution_exacte") 
# plt.plot(X,Uh, label="solution_approchee")
# plt.legend(loc="upper left")
# plt.ylim(0.8, 3.5)
# plt.show()


##########################


# etude de la convergence du schema fct de h = L/N 
Nmesh = [10, 20, 40, 80, 160, 320, 640]

sizeh = np.zeros(len(Nmesh))
erreurl2 = np.zeros(len(Nmesh))
erreurh10 = np.zeros(len(Nmesh))

for k, N in enumerate(Nmesh):

    h = L / N
    sizeh[k] = h

    X, Uh = VF(f, uD, g, L, N)

    # erreur L2 discrete
    e = u(X) - Uh
    erreurl2[k] = np.sqrt(h * np.sum(e**2))

    # erreur H1,0 discrete
    erreurh10[k] = np.sqrt((0 - (Uh[-1] - u(X)[-1]))**2/h + np.sum((e[:-1] - e[1:])**2) / h)


# ---------------- plot erreurs ----------------

plt.figure(2)
plt.loglog(sizeh, erreurl2, "-xb", label="erreur L2")
plt.loglog(sizeh, erreurh10, "-xr", label="erreur H1,0")
plt.xlabel("pas du maillage h")
plt.ylabel("erreur")
plt.legend()
plt.show()


# ---------------- ordres de convergence ----------------

pente_l2 = np.polyfit(np.log(sizeh), np.log(erreurl2), 1)[0]
print("ordre de convergence L2 =", pente_l2)

pente_h10 = np.polyfit(np.log(sizeh), np.log(erreurh10), 1)[0]
print("ordre de convergence H1,0 =", pente_h10)