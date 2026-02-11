import numpy as np
from math import pi, sqrt
from scipy.sparse import lil_matrix, csr_matrix
from scipy.sparse.linalg import spsolve
import matplotlib.pyplot as plt
import matplotlib.tri as mtri

# =========================================================
# EXACT SOLUTION
# =========================================================
def Exact(x, y, gamx, gamy):    
    return np.sin(gamx*x)*np.cos(gamy*y)

# =========================================================
# RHS
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
# BASIS FUNCTIONS
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
# GRAD PHI
# =========================================================

def GradPhi_Pk(l1, l2, G1, G2, pk):

    l3 = 1 - l1 - l2
    G3 = -G1 - G2

    if pk == 1:
        return np.vstack([G1, G2, G3])

    if pk == 2:
        grads = []
        # Vertices
        grads.append((4*l1-1)*G1)
        grads.append((4*l2-1)*G2)
        grads.append((1-4*l3)*G1 + (1-4*l3)*G2)
        # Edges
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
# GAUSS POINTS
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
# MAILLEUR Pk (SCILAB -> PYTHON)
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
# ================= CONVERGENCE STUDY =====================
# =========================================================

EF_Pk = 3

Lx=Ly=1
gamx=pi/Lx; gamy=2*pi/Ly

betx=1;bety=2
cx=1; cy=0.5
alpha=-5

W,Lambda,Ngp = IntegrationNum()

Vconv = [6, 10, 14, 18, 24]

ErrL1 = []
ErrL2 = []
ErrLinf = []
H = []
NddLC = []   # nombre de ddl
VNx = []     # Nx (optionnel, comme Scilab)

for Nx in Vconv:

    Ny = Nx
    Coor,Nu,LogP = Struct_Pk_Mesh(Nx,Ny,Lx,Ly,EF_Pk)
    print(f"Mesh generated with Nx={Nx}, Ny={Ny}")
    Np = Coor.shape[1]
    NddLC.append(Np)
    VNx.append(Nx)
    Ne = Nu.shape[1]
    Ndl = (EF_Pk+1)*(EF_Pk+2)//2
    print("Nu[:,0] =", Nu[:,0])
    print("Np =", Np)
    A = lil_matrix((Np,Np))
    B = np.zeros(Np)

    for e in range(Ne):
        Xi=Coor[:,Nu[0,e]]; Xj=Coor[:,Nu[1,e]]; Xk=Coor[:,Nu[2,e]]
        DetE=abs(Determinant(Xj-Xi,Xk-Xi))
        GradP1=MGrad_Lambda_of_P1(Xi,Xj,Xk)

        for gp in range(Ngp):
            lmd=Lambda[:,gp]
            Xg=lmd[0]*Xi+lmd[1]*Xj+lmd[2]*Xk
            wg=W[gp]

            Phi=Phi_Pk(lmd[0],lmd[1],EF_Pk)
            GradPhi=GradPhi_Pk(lmd[0],lmd[1],GradP1[0],GradP1[1],EF_Pk)

            for k in range(Ndl):
                isg=Nu[k,e]
                B[isg]+=wg*DetE*f_rhs(Xg,gamx,gamy,betx,bety,cx,cy,alpha)*Phi[k]

                for kp in range(Ndl):
                    jsg=Nu[kp,e]
                    diff = betx*GradPhi[k,0]*GradPhi[kp,0]+bety*GradPhi[k,1]*GradPhi[kp,1]
                    reac = alpha*Phi[k]*Phi[kp]
                    conv = -(cx*GradPhi[k,0]+cy*GradPhi[k,1])*Phi[kp]
                    A[isg,jsg]+=wg*DetE*(diff+reac+conv)

    for i in range(Np):
        if LogP[i]<0:
            A[i,:]=0; A[i,i]=1
            B[i]=Exact(Coor[0,i],Coor[1,i],gamx,gamy)

    U = spsolve(csr_matrix(A),B)

    # -------- ERRORS L1, L2, Linf --------
    errL1 = 0.0
    errL2 = 0.0
    errLinf = 0.0

    for e in range(Ne):
        Xi=Coor[:,Nu[0,e]]; Xj=Coor[:,Nu[1,e]]; Xk=Coor[:,Nu[2,e]]
        DetE=abs(Determinant(Xj-Xi,Xk-Xi))

        for gp in range(Ngp):
            lmd=Lambda[:,gp]; wg=W[gp]
            Xg=lmd[0]*Xi+lmd[1]*Xj+lmd[2]*Xk

            Phi=Phi_Pk(lmd[0],lmd[1],EF_Pk)
            Ug = sum(U[Nu[k,e]]*Phi[k] for k in range(Ndl))
            Uex = Exact(Xg[0],Xg[1],gamx,gamy)

            diff = abs(Ug - Uex)
            errL1   += wg*DetE*diff
            errL2   += wg*DetE*diff**2
            errLinf = max(errLinf, diff)

    ErrL1.append(errL1)
    ErrL2.append(np.sqrt(errL2))
    ErrLinf.append(errLinf)

    H.append(np.sqrt(Lx*Ly/Ne))

ErrL1 = np.array(ErrL1)
ErrL2 = np.array(ErrL2)
ErrLinf = np.array(ErrLinf)
H = np.array(H)
NddLC = np.array(NddLC)
VNx = np.array(VNx)

Dh = H[:-1] / H[1:]

OrdreL1   = np.log(ErrL1[:-1]   / ErrL1[1:])   / np.log(Dh)
OrdreL2   = np.log(ErrL2[:-1]   / ErrL2[1:])   / np.log(Dh)
OrdreLinf = np.log(ErrLinf[:-1] / ErrLinf[1:]) / np.log(Dh)

print("\nEF_Pk =", EF_Pk)
print("Nddl | Ordre L1 | Ordre L2 | Ordre Linf | Erreur L2")
print("------------------------------------------------------")

for i in range(1, len(H)):
    print(f"{NddLC[i]:5d} | "
          f"{OrdreL1[i-1]:8.3f} | "
          f"{OrdreL2[i-1]:8.3f} | "
          f"{OrdreLinf[i-1]:10.3f} | "
          f"{ErrL2[i]:10.4e}")

# =========================================================
# =================== PLOT LAST MESH ======================
# =========================================================

x = Coor[0]; y = Coor[1]
triangles = Nu[:3,:].T
triang = mtri.Triangulation(x,y,triangles)

Uex = Exact(x,y,gamx,gamy)
Err = np.abs(U-Uex)

plt.figure(figsize=(18,5))

plt.subplot(1,3,1)
plt.tricontourf(triang,U,30); plt.colorbar(); plt.triplot(triang,lw=0.3)
plt.title("Solution EF EF_Pk=" + str(EF_Pk))

plt.subplot(1,3,2)
plt.tricontourf(triang,Uex,30); plt.colorbar(); plt.triplot(triang,lw=0.3)
plt.title("Solution exacte")

plt.subplot(1,3,3)
plt.tricontourf(triang,Err,30); plt.colorbar(); plt.triplot(triang,lw=0.3)
plt.title("Erreur L2 EF_Pk=" + str(EF_Pk))

plt.tight_layout()
plt.show()

p_mes = OrdreL2[-1]

plt.figure()

plt.loglog(H, ErrL2, 'o-', label='Erreur L2')

C = ErrL2[-1]/(H[-1]**p_mes)
plt.loglog(H, C*H**p_mes, '--', label=f'Pente mesurée ≈ {p_mes:.2f}')

plt.gca().invert_xaxis()
plt.xlabel('h')
plt.ylabel('Erreur L2')
plt.title(f'Convergence FEM P{EF_Pk}')
plt.legend()
plt.grid(True, which="both", ls="--")
plt.show()
"""



# liste des puissances à tester
m_list = [1, 2, 3]

ErrL2_m = {}  # dictionnaire pour stocker l'erreur pour chaque m

for m in m_list:
    ErrL2 = []
    H = []

    for Nx in Vconv:

        Ny = Nx
        Coor,Nu,LogP = Struct_Pk_Mesh(Nx,Ny,Lx,Ly,EF_Pk)

        Np = Coor.shape[1]
        Ne = Nu.shape[1]
        # taille caractéristique par élément
        dh = np.sqrt(Lx * Ly / Ne)  # dh élémentaire pour Ne éléments
        Ndl = (EF_Pk+1)*(EF_Pk+2)//2
        betx = dh**m
        bety = dh**m
        A = lil_matrix((Np,Np))
        B = np.zeros(Np)

        for e in range(Ne):
            Xi=Coor[:,Nu[0,e]]; Xj=Coor[:,Nu[1,e]]; Xk=Coor[:,Nu[2,e]]
            DetE=abs(Determinant(Xj-Xi,Xk-Xi))
            GradP1=MGrad_Lambda_of_P1(Xi,Xj,Xk)

            for gp in range(Ngp):
                lmd=Lambda[:,gp]
                Xg=lmd[0]*Xi+lmd[1]*Xj+lmd[2]*Xk
                wg=W[gp]

                Phi=Phi_Pk(lmd[0],lmd[1],EF_Pk)
                GradPhi=GradPhi_Pk(lmd[0],lmd[1],GradP1[0],GradP1[1],EF_Pk)

                for k in range(Ndl):
                    isg=Nu[k,e]
                    B[isg]+=wg*DetE*f_rhs(Xg,gamx,gamy,betx,bety,cx,cy,alpha)*Phi[k]

                    for kp in range(Ndl):
                        jsg=Nu[kp,e]
                        diff = betx*GradPhi[k,0]*GradPhi[kp,0]+bety*GradPhi[k,1]*GradPhi[kp,1]
                        reac = alpha*Phi[k]*Phi[kp]
                        conv = -(cx*GradPhi[k,0]+cy*GradPhi[k,1])*Phi[kp]
                        A[isg,jsg]+=wg*DetE*(diff+reac+conv)

        for i in range(Np):
            if LogP[i]<0:
                A[i,:]=0; A[i,i]=1
                B[i]=Exact(Coor[0,i],Coor[1,i],gamx,gamy)

        U = spsolve(csr_matrix(A),B)

        # -------- L2 ERROR --------
        err = 0
        for e in range(Ne):
            Xi=Coor[:,Nu[0,e]]; Xj=Coor[:,Nu[1,e]]; Xk=Coor[:,Nu[2,e]]
            DetE=abs(Determinant(Xj-Xi,Xk-Xi))
            for gp in range(Ngp):
                lmd=Lambda[:,gp]; wg=W[gp]
                Xg=lmd[0]*Xi+lmd[1]*Xj+lmd[2]*Xk
                Phi=Phi_Pk(lmd[0],lmd[1],EF_Pk)
                Ug = sum(U[Nu[k,e]]*Phi[k] for k in range(Ndl))
                Uex = Exact(Xg[0],Xg[1],gamx,gamy)
                err += wg*DetE*(Ug-Uex)**2

        ErrL2.append(sqrt(err))
        H.append(sqrt(Lx*Ly/Ne))

    ErrL2_m[m] = (np.array(H), np.array(ErrL2))

plt.figure()
for m in m_list:
    H, ErrL2 = ErrL2_m[m]
    ordre = np.log(ErrL2[:-1]/ErrL2[1:]) / np.log(H[:-1]/H[1:])
    print(f"Ordre moyen pour m={m} : {np.mean(ordre):.2f}")
    plt.loglog(H, ErrL2, 'o-', label=f'm = {m}')

plt.gca().invert_xaxis()
plt.xlabel('h')
plt.ylabel('Erreur L2')
plt.title(f'Convergence FEM P{EF_Pk} selon m')
plt.legend()
plt.grid(True, which="both", ls="--")
plt.show()
"""