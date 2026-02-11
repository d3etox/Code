"""
================================================================================
  DEVOIR MAISON : ÉLÉMENTS FINIS P1/P2/P3
  Problème de Convection-Réaction-Diffusion 2D
================================================================================
Objectif : Implémenter et valider les éléments finis Pk-Lagrange pour résoudre
un problème de convection-réaction-diffusion anisotrope en 2D.

Paramètres communs :
  Lx = Ly = 1
  gamx = π, gamy = 2π
  Solution exacte : u(x,y) = sin(π*x)*cos(2π*y)
================================================================================
"""

import numpy as np
from math import pi
from scipy.sparse import lil_matrix, csr_matrix
from scipy.sparse.linalg import spsolve

# =========================================================
# SOLUTION EXACTE
# =========================================================
def Exact(x, y, gamx, gamy):    
    return np.sin(gamx*x)*np.cos(gamy*y)

# =========================================================
# SECOND MEMBRE (RHS)
# =========================================================
def f_rhs(X, gamx, gamy, betx, bety, cx, cy, alpha):
    x, y = X
    return ((alpha + betx*gamx**2 + bety*gamy**2)
            * np.sin(gamx*x)*np.cos(gamy*y)
            + cx*gamx*np.cos(gamx*x)*np.cos(gamy*y)
            - cy*gamy*np.sin(gamx*x)*np.sin(gamy*y))

# =========================================================
def Determinant(p, q):
    return p[0]*q[1] - p[1]*q[0]

# =========================================================
def MGrad_Lambda_of_P1(Xi, Xj, Xk):
    e1 = np.array([1.,0.]); e2 = np.array([0.,1.])
    G = np.zeros((3,2))
    G[0] = [Determinant(e1,Xk-Xj),Determinant(e2,Xk-Xj)]/Determinant(Xi-Xj,Xk-Xj)
    G[1] = [Determinant(e1,Xi-Xk),Determinant(e2,Xi-Xk)]/Determinant(Xj-Xk,Xi-Xk)
    G[2] = [Determinant(e1,Xj-Xi),Determinant(e2,Xj-Xi)]/Determinant(Xk-Xi,Xj-Xi)
    return G

# =========================================================
# FONCTIONS DE FORME
# =========================================================
def Phi_Pk(l1,l2,pk):
    l3 = 1-l1-l2
    if pk==1:
        return np.array([l1,l2,l3])
    if pk==2:
        return np.array([
            l1*(2*l1-1), l2*(2*l2-1), l3*(2*l3-1),
            4*l1*l2, 4*l2*l3, 4*l3*l1
        ])
    if pk==3:
        l3 = 1 - l1 - l2
        return np.array([
            l1*(3*l1-1)*(3*l1-2)/2,          # 0: sommet 1 (λ₁)
            l2*(3*l2-1)*(3*l2-2)/2,          # 1: sommet 2 (λ₂)
            l3*(3*l3-1)*(3*l3-2)/2,          # 2: sommet 3 (λ₃)
            (9/2)*l1*l2*(3*l1-1),             # 3: arête 1-2 près de 1
            (9/2)*l1*l2*(3*l2-1),             # 4: arête 1-2 près de 2
            (9/2)*l2*l3*(3*l2-1),             # 5: arête 2-3 près de 2
            (9/2)*l2*l3*(3*l3-1),             # 6: arête 2-3 près de 3
            (9/2)*l3*l1*(3*l3-1),             # 7: arête 3-1 près de 3
            (9/2)*l3*l1*(3*l1-1),             # 8: arête 3-1 près de 1
            27*l1*l2*l3                       # 9: intérieur
        ])

# =========================================================
# GRADIENT DES FONCTIONS DE FORME
# =========================================================
def GradPhi_Pk(l1, l2, G1, G2, pk):
    l3 = 1 - l1 - l2
    G3 = -G1 - G2

    if pk == 1:
        return np.vstack([G1, G2, G3])

    if pk == 2:
        grads = []
        # Sommets
        grads.append((4*l1-1)*G1)
        grads.append((4*l2-1)*G2)
        grads.append((1-4*l3)*G1 + (1-4*l3)*G2)
        # Arêtes
        grads.append(4*l2*G1 + 4*l1*G2)
        grads.append(-4*l2*G1 + (-4*l2+4*l3)*G2)
        grads.append((-4*l1+4*l3)*G1 - 4*l1*G2)
        return np.array(grads)

    if pk == 3:
        grads = []

        # --- sommets ---
        grads.append(((27*l1**2 - 18*l1 + 2)/2) * G1)
        grads.append(((27*l2**2 - 18*l2 + 2)/2) * G2)
        grads.append((-27*l3**2/2 + 9*l3 - 1) * G1 + (-27*l3**2/2 + 9*l3 - 1) * G2)

        # --- arête 1-2 ---
        grads.append((27.0*l1*l2 - 4.5*l2)*G1 + (13.5*l1**2 - 4.5*l1)*G2)
        grads.append((13.5*l2**2 - 4.5*l2)*G1 + (27.0*l1*l2 - 4.5*l1)*G2)

        # --- arête 2-3 ---
        grads.append((-13.5*l2**2 + 4.5*l2)*G1 + (-13.5*l2**2 + 27.0*l2*l3 + 4.5*l2 - 4.5*l3)*G2)
        grads.append((-27.0*l2*l3 + 4.5*l2)*G1 + (-27.0*l2*l3 + 4.5*l2 + 13.5*l3**2 - 4.5*l3)*G2)

        # --- arête 3-1 ---
        grads.append((-27.0*l1*l3 + 4.5*l1 + 13.5*l3**2 - 4.5*l3)*G1 + (-27.0*l1*l3 + 4.5*l1)*G2)
        grads.append((-13.5*l1**2 + 27.0*l1*l3 + 4.5*l1 - 4.5*l3)*G1 + (-13.5*l1**2 + 4.5*l1)*G2)

        # --- intérieur ---
        grads.append((-27*l1*l2 + 27*l2*l3)*G1 + (-27*l1*l2 + 27*l1*l3)*G2)

        return np.array(grads)


# =========================================================
# POINTS DE QUADRATURE DE GAUSS
# =========================================================
def IntegrationNum():
    s1=0.11169079483905; s2=0.0549758718227661
    aa=0.445948490915965; bb=0.091576213509771
    W=np.array([s2,s2,s2,s1,s1,s1])
    Xsi=np.zeros((3,6))
    Xsi[0]=[bb,1-2*bb,bb,aa,aa,1-2*aa]
    Xsi[1]=[bb,bb,1-2*bb,1-2*aa,aa,aa]
    Xsi[2]=1-Xsi[0]-Xsi[1]
    return W,Xsi,6

# =========================================================
# GÉNÉRATEUR DE MAILLAGE Pk
# =========================================================
def Struct_Pk_Mesh(Nx, Ny, Lx, Ly, pk):
    MNx = Nx + (pk-1)*(Nx-1)
    MNy = Ny + (pk-1)*(Ny-1)

    dx = Lx/(MNx-1)
    dy = Ly/(MNy-1)

    Coor = np.zeros((2, MNx*MNy))
    LogP = np.zeros(MNx*MNy, dtype=int)

    isg = -1
    for i in range(MNx):
        xi = i*dx
        for j in range(MNy):
            yj = j*dy
            isg += 1
            Coor[:,isg] = [xi, yj]
            if i==0: LogP[isg] = -1
            elif i==MNx-1: LogP[isg] = -2
            elif j==0: LogP[isg] = -10
            elif j==MNy-1: LogP[isg] = -20

    N_elements = 2*(Nx-1)*(Ny-1)
    Ndl = (pk+1)*(pk+2)//2
    Nu = np.zeros((Ndl, N_elements), dtype=int)

    ie = -1
    for i in range(Nx-1):
        for j in range(Ny-1):

            is0 = 1 + pk*j + i*pk*MNy
            js = is0 + pk*MNy
            ks = js + pk
            ps = is0 + pk

            ie += 1
            Nu[0:3,ie] = [is0, js, ps]

            il = 3
            for lx in range(1,pk):
                Nu[il,ie] = is0 + lx*MNy; il+=1
            for lo in range(1,pk):
                Nu[il,ie] = js + lo - lo*MNy; il+=1
            for ly in range(1,pk):
                Nu[il,ie] = ps - ly; il+=1
            for lx in range(1,pk-1):
                for ly in range(1,lx+1):
                    Nu[il,ie] = is0 + lx*MNy + ly; il+=1

            ie += 1
            Nu[0:3,ie] = [js, ks, ps]

            il = 3
            for ly in range(1,pk):
                Nu[il,ie] = js + ly; il+=1
            for lx in range(1,pk):
                Nu[il,ie] = ks - lx*MNy; il+=1
            for lo in range(1,pk):
                Nu[il,ie] = ps - lo + lo*MNy; il+=1
            for ly in range(1,pk-1):
                for lx in range(1,ly+1):
                    Nu[il,ie] = js - lx*MNy + ly + 1; il+=1

    return Coor, Nu-1, LogP


# =========================================================
# FONCTION DE RÉSOLUTION
# =========================================================
def solve_FEM(Nx, Ny, Lx, Ly, pk, gamx, gamy, betx, bety, cx, cy, alpha):
    """
    Résout le problème d'éléments finis Pk.
    Retourne la solution U, les coordonnées Coor, la connectivité Nu.
    """
    Coor, Nu, LogP = Struct_Pk_Mesh(Nx, Ny, Lx, Ly, pk)
    Np = Coor.shape[1]
    Ne = Nu.shape[1]
    Ndl = (pk+1)*(pk+2)//2
    
    W, Lambda, Ngp = IntegrationNum()
    
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
            
            Phi = Phi_Pk(lmd[0], lmd[1], pk)
            GradPhi = GradPhi_Pk(lmd[0], lmd[1], GradP1[0], GradP1[1], pk)
            
            for k in range(Ndl):
                isg = Nu[k, e]
                B[isg] += wg*DetE*f_rhs(Xg, gamx, gamy, betx, bety, cx, cy, alpha)*Phi[k]
                
                for kp in range(Ndl):
                    jsg = Nu[kp, e]
                    diff = betx*GradPhi[k, 0]*GradPhi[kp, 0] + bety*GradPhi[k, 1]*GradPhi[kp, 1]
                    reac = alpha*Phi[k]*Phi[kp]
                    conv = -(cx*GradPhi[k, 0] + cy*GradPhi[k, 1])*Phi[kp]
                    A[isg, jsg] += wg*DetE*(diff + reac + conv)
    
    # Conditions aux limites : u = 0 sur ∂Ω
    for i in range(Np):
        if LogP[i] < 0:
            A[i, :] = 0
            A[i, i] = 1
            B[i] = Exact(Coor[0, i], Coor[1, i], gamx, gamy)
    
    U = spsolve(csr_matrix(A), B)
    return U, Coor, Nu


# =========================================================
# CALCUL DE L'ERREUR
# =========================================================
def compute_L2_error(Nx, Ny, Lx, Ly, pk, gamx, gamy, betx, bety, cx, cy, alpha, U, Coor, Nu):
    """Calcule l'erreur L2"""
    Ne = Nu.shape[1]
    Ndl = (pk+1)*(pk+2)//2
    W, Lambda, Ngp = IntegrationNum()
    
    errL2 = 0.0
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
            
            Phi = Phi_Pk(lmd[0], lmd[1], pk)
            Ug = sum(U[Nu[k, e]]*Phi[k] for k in range(Ndl))
            Uex = Exact(Xg[0], Xg[1], gamx, gamy)
            
            errL2 += wg*DetE*(Ug - Uex)**2
    
    return np.sqrt(errL2)


# =========================================================
# ÉTUDE DE CONVERGENCE
# =========================================================
def convergence_study(case_name, pk_values, Lx, Ly, gamx, gamy, 
                      betx, bety, cx, cy, alpha, Vconv, m_stab=None):
    """
    Étude générique de convergence pour un ensemble donné de paramètres.
    Si m_stab est fourni, betx et bety sont calculés comme h^m pour chaque taille de maillage.
    Retourne des dictionnaires avec les erreurs et les valeurs de h.
    """
    print(f"\n{'='*70}")
    print(f"  {case_name}")
    print(f"{'='*70}")
    if m_stab is None:
        print(f"Paramètres: betx={betx}, bety={bety}, cx={cx}, cy={cy}, alpha={alpha}")
    else:
        print(f"Paramètres: betx=bety=h^{m_stab}, cx={cx}, cy={cy}, alpha={alpha}")
    print(f"{'='*70}\n")
    
    results = {}
    
    for pk in pk_values:
        print(f"Convergence P{pk}:")
        print(f"  Nx   | Nddl    | Err L2      | h")
        print(f"  " + "-"*45)
        
        ErrL2 = []
        H = []
        
        for Nx in Vconv:
            Ny = Nx
            Ne = 2*(Nx-1)*(Ny-1)
            h = np.sqrt(Lx*Ly/Ne)
            
            # Compute betx, bety based on h if m_stab is provided
            if m_stab is not None:
                betx_local = h**m_stab
                bety_local = h**m_stab
            else:
                betx_local = betx
                bety_local = bety
            
            U, Coor, Nu = solve_FEM(Nx, Ny, Lx, Ly, pk, gamx, gamy, betx_local, bety_local, cx, cy, alpha)
            
            err = compute_L2_error(Nx, Ny, Lx, Ly, pk, gamx, gamy, betx_local, bety_local, cx, cy, alpha, U, Coor, Nu)
            
            ErrL2.append(err)
            H.append(h)
            
            Np = Coor.shape[1]
            print(f"  {Nx:3d}  | {Np:6d}  | {err:.4e} | {h:.4e}")
        
        ErrL2 = np.array(ErrL2)
        H = np.array(H)
        
        # Calcule l'ordre de convergence
        if len(H) > 1:
            Dh = H[:-1] / H[1:]
            Ordre = np.log(ErrL2[:-1] / ErrL2[1:]) / np.log(Dh)
            ordre_mean = np.mean(Ordre)
            print(f"\n  Ordre de convergence (moyenne) : {ordre_mean:.3f}")
            print(f"  Ordre théorique : {pk+1}")
        
        results[pk] = {'H': H, 'ErrL2': ErrL2}
        print()
    
    return results


# =========================================================
# PROGRAMME PRINCIPAL
# =========================================================
if __name__ == "__main__":
    
    # Paramètres communs
    Lx = Ly = 1.0
    gamx = pi / Lx
    gamy = 2*pi / Ly
    
    # Séquence de raffinement du maillage
    Vconv = [6, 10, 14, 18, 24]
    
    # Cas à étudier
    cases = []
    
    # Cas 1 : Réaction Pure
    cases.append(("Réaction Pure (α=-5)", [1, 2, 3], 0, 0, 0, 0, -5, None))
    
    # Cas 2 : Diffusion Isotrope
    cases.append(("Diffusion Isotrope (β=2)", [1, 2, 3], 2, 2, 0, 0, 0, None))
    
    # Cas 3.1 : Diffusion Anisotrope (non-dégénérée)
    cases.append(("Diffusion Anisotrope (β_x=1, β_y=2)", [1, 2, 3], 1, 2, 0, 0, 0, None))
    
    # Cas 3.2 : Diffusion Dégénérée
    cases.append(("Diffusion Dégénérée (β_x=2, β_y=0)", [1, 2, 3], 2, 0, 0, 0, 0, None))
    
    # Cas 3.3 : Fortement Anisotrope
    cases.append(("Fortement Anisotrope (β_x=1e-8, β_y=1)", [1, 2, 3], 1e-8, 1, 0, 0, 0, None))
    
    # Cas 4 : Convection Pure avec diffusion (pour la stabilisation)
    # Test m=1, 2, 3
    for m_value in [1, 2, 3]:
        cases.append((f"Convection Pure avec stab (m={m_value}, c_x=1, c_y=0.5)", [1, 2, 3], 
                      0, 0, 1, 0.5, 0, m_value))
    
    # Cas 5 : Convection-Réaction-Diffusion
    cases.append(("Convection-Réaction-Diffusion (β_x=1, β_y=2, c_x=1, c_y=0.5, α=-5)", 
                  [1, 2, 3], 1, 2, 1, 0.5, -5, None))
    
    # Exécute tous les cas
    for case_data in cases:
        case_name, pk_values, betx, bety, cx, cy, alpha, m_stab = case_data
        results = convergence_study(case_name, pk_values, Lx, Ly, gamx, gamy,
                                   betx, bety, cx, cy, alpha, Vconv, m_stab)
    
    print("\n" + "="*70)
    print("  ÉTUDE DE CONVERGENCE TERMINÉE")
    print("="*70)
