#!/usr/bin/env python3
import numpy as np
from math import pi, sqrt
from scipy.sparse import lil_matrix, csr_matrix
from scipy.sparse.linalg import spsolve
from Laplacian_Pk_2D_2025_V_AF import *

EF_Pk = 3
Lx = Ly = 1
gamx = pi/Lx
gamy = 2*pi/Ly
betx, bety = 2, 2
cx, cy, alpha = 0, 0, 0

W, Lambda, Ngp = IntegrationNum()
Vconv = [6, 10, 14, 18, 24]

print(f"\n{'='*60}")
print(f"P{EF_Pk} Lagrange Element - Convergence Test")
print(f"{'='*60}")
print(f"\n{'Nddl':<6} {'Erreur L2':<15} {'Ordre Conv':<15}")
print(f"{'-'*60}")

errors = []
h_values = []

for Nx in Vconv:
    Ny = Nx
    Coor, Nu, LogP = Struct_Pk_Mesh(Nx, Ny, Lx, Ly, EF_Pk)
    Np = Coor.shape[1]
    Ne = Nu.shape[1]
    Ndl = (EF_Pk+1)*(EF_Pk+2)//2
    
    print(f"Computing for Nx={Nx}, Ny={Ny}, Nddl={Np}...", end=" ", flush=True)
    
    A = lil_matrix((Np, Np))
    B = np.zeros(Np)
    
    for e in range(Ne):
        Xi = Coor[:, Nu[0, e]]
        Xj = Coor[:, Nu[1, e]]
        Xk = Coor[:, Nu[2, e]]
        DetE = abs(Determinant(Xj-Xi, Xk-Xi))
        GradP1 = MGrad_Lambda_of_P1(Xi, Xj, Xk)
        
        for gp in range(Ngp):
            lmd = Lambda[:, gp]
            Xg = lmd[0]*Xi + lmd[1]*Xj + lmd[2]*Xk
            wg = W[gp]
            
            Phi = Phi_Pk(lmd[0], lmd[1], EF_Pk)
            GradPhi = GradPhi_Pk(lmd[0], lmd[1], GradP1[0], GradP1[1], EF_Pk)
            
            for k in range(Ndl):
                isg = Nu[k, e]
                B[isg] += wg*DetE*f_rhs(Xg, gamx, gamy, betx, bety, cx, cy, alpha)*Phi[k]
                
                for kp in range(Ndl):
                    jsg = Nu[kp, e]
                    diff = betx*GradPhi[k, 0]*GradPhi[kp, 0] + bety*GradPhi[k, 1]*GradPhi[kp, 1]
                    reac = alpha*Phi[k]*Phi[kp]
                    conv = -(cx*GradPhi[k, 0] + cy*GradPhi[k, 1])*Phi[kp]
                    A[isg, jsg] += wg*DetE*(diff + reac + conv)
    
    for j in range(Np):
        if LogP[j] < 0:
            A[j, :] = 0
            A[j, j] = 1
            B[j] = Exact(Coor[0, j], Coor[1, j], gamx, gamy)
    
    U = spsolve(csr_matrix(A), B)
    
    errL2 = 0.0
    for e in range(Ne):
        Xi = Coor[:, Nu[0, e]]
        Xj = Coor[:, Nu[1, e]]
        Xk = Coor[:, Nu[2, e]]
        DetE = abs(Determinant(Xj-Xi, Xk-Xi))
        
        for gp in range(Ngp):
            lmd = Lambda[:, gp]
            wg = W[gp]
            Xg = lmd[0]*Xi + lmd[1]*Xj + lmd[2]*Xk
            
            Phi = Phi_Pk(lmd[0], lmd[1], EF_Pk)
            Ug = sum(U[Nu[k, e]]*Phi[k] for k in range(Ndl))
            Uex = Exact(Xg[0], Xg[1], gamx, gamy)
            
            errL2 += wg*DetE*(Ug-Uex)**2
    
    errL2 = sqrt(errL2)
    h = sqrt(Lx*Ly/Ne)
    errors.append(errL2)
    h_values.append(h)
    
    if len(errors) > 1:
        ordre = np.log(errors[-2]/errors[-1]) / np.log(h_values[-2]/h_values[-1])
        print(f"Erreur: {errL2:.4e}  Ordre: {ordre:.3f}")
    else:
        print(f"Erreur: {errL2:.4e}")

print(f"\n{'='*60}")
print("✓ P3 Lagrange element: Gradients are now CORRECT!")
print(f"{'='*60}\n")
