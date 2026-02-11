#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

@author: masssonr

    Stratigraphic model with diffusive sediment transport 
    Modelize the formation of sedimentary basins at large space and time scales 

    Single lithology Model 


  d_t h(x,t) + div(grad(psi(b(x,t)))) = 0 on (0,Lx)x(0,tf)

  h(x,0) = hinit(x)   on (0,Lx)

  grad(psi).n = g0 for x=0,  grad(psi).n = g1 for x = Lx  

  b(x,t) = hsea(t) - h(x,t)


  Finite Volume discretization using TPFA fluxes on unstructured meshes 

  Ouputs: * h(x,t), b(x,t)

   
          * hs(x,t) = min_(t<=s<=tf) h(x,s) (sediment layers at each time t, taking into account erosions)



"""

import numpy as np
import matplotlib.pyplot as plt
import sys



#data
Lx=2
#space discretization
N=100
Nint = N-1
Nbound = 2
dx=Lx/N

#time discretization
tf = 1.5 # final simulation time
ndt = 50
dt = tf/ndt # initial time step

#Newton convergence 
Newtmax = 50 # maximum number of Newton iterations 
eps=1.0e-6 #  stopping criteria 


Km = 1
Kc = 10

g0 = -20
g1 = 0


def f_hsea(t):
    s = 25 + 5*np.cos(12*t) 
    return s

def f_psi(u):
    if (u<0):
        s = Kc*u
    else: 
        s = Km*u
    return s

def f_psip(u):
    if (u<0):
        s = Kc
    else: 
        s = Km
    return s


def f_hinit(x):
    s = 25*np.exp(-8*x/Lx) + 10
    return s


def residual(h,b,h0,dt):
    R=np.zeros(N)

    for i in range(N):
        R[i] = volume[i]*(h[i] - h0[i])/dt

    for i in range(Nbound):
        k = CellbyFaceBound[i]
        R[k] += fbound[i]

    # Fluxes aux faces intérieures
    for i in range(Nint):
        k1 = CellsbyFaceInt[i,0]
        k2 = CellsbyFaceInt[i,1]
        flux = Tint[i]*(f_psi(b[k2])-f_psi(b[k1]))
        R[k1] += flux
        R[k2] -= flux
    return R


def Jacobian(b,dt):
    A = np.zeros([N,N])

    for i in range(N):
        A[i,i] = volume[i]/dt

    for i in range(Nint):
        k1 = CellsbyFaceInt[i,0]
        k2 = CellsbyFaceInt[i,1]
        dflux1 = Tint[i]*f_psip(b[k1])
        dflux2 = - Tint[i]*f_psip(b[k2])
        A[k1,k1] += dflux1
        A[k1,k2] += dflux2
        A[k2,k1] -= dflux1
        A[k2,k2] -= dflux2

    return A 



#data structure for the uniform 1D mesh of size dx of the domain (0,Lx)
# cells m = 0:N-1

X = np.linspace(dx/2,Lx-dx/2,N)
volume = dx*np.ones(N)
 

#Interior faces: i = 0:Nint-1
CellsbyFaceInt = np.zeros([Nint,2],dtype=int)
surfaceint = np.ones(Nint)  # surface des faces intérieures
for i in range(Nint):
    CellsbyFaceInt[i,0] = i      # cellule gauche
    CellsbyFaceInt[i,1] = i + 1  # cellule droite 
  

#Boundary
CellbyFaceBound = np.zeros(Nbound,dtype=int)
fbound = np.zeros(Nbound)
CellbyFaceBound[0] = 0        # face gauche, cellule 0
CellbyFaceBound[1] = N-1      # face droite, cellule N-1
fbound[0] = g0  # flux à x=0
fbound[1] = g1  # flux à x=Lx 


#transmissibilities of interior faces
Tint = np.zeros(Nint)
for i in range(Nint):
    k1 = CellsbyFaceInt[i,0]
    k2 = CellsbyFaceInt[i,1]
    Tint[i] = surfaceint[i]/np.abs(X[k2]-X[k1])
    

#  simulation 

#initialization 
h0 = f_hinit(X)
h = h0 
t = 0

hs = np.zeros([N,ndt+1])
hs[:,0] = h0

plt.figure(1)
plt.title('h')
plt.plot (X,h0,'-r')  
plt.show()

for n in range(ndt):  # time loop 
    t = t + dt    
    hsea = f_hsea(t)
    
    b = np.ones(N)*hsea - h # bathymetry
    R = residual(h,b,h0,dt) # initial newton residual 
    normR = np.linalg.norm(R)
    normR0 = normR

    itn = 1
    while ((normR/normR0>eps)and(itn<Newtmax)): #Newton loop 
        itn = itn + 1
        A = Jacobian(b,dt) #Jacobian matrix A = dR/dh(h)
        dh = -np.linalg.solve(A,R) #solution dh 
        h = h + dh #Newton update         
        b = np.ones(N)*hsea - h
        R = residual(h,b,h0,dt) #residual  
        normR = np.linalg.norm(R) #residual norm  

    if (itn % 10 == 0):
        print('it newton',itn,' normR ',normR)
    if (itn>Newtmax-1): # if Newton not converged
        print('it newton',itn)
        print('newton non converge',itn,normR)
        sys.exit()

    h0 = h 
    hs[:,n+1] = h    
 
 

 #plot h 
    plt.figure(1)
    plt.title('h') 
    plt.plot (X,h,'-r')  
    
plt.show()
# computation of hs taking out erosions 
for j in range(ndt-1,-1,-1):
    for k in range(N):
        if j == ndt-1:
            hs[k,j] = hs[k,j+1]
        else:
            hs[k,j] = min(hs[k,j],hs[k,j+1]) 
        

# plt.figure(2)
# plt.title('hs')
# for it in range(ndt+1):
#     plt.plot(X,hs[:,it],'-b')





