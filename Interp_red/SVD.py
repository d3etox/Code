import numpy as np

def pui_itere(A,x0,eps):
    x = x0
    lam_, lam = 1,0
    it = 0
    while np.abs(lam-lam_)>eps:
        it += 1
        lam_ = lam
        norm = np.linalg.norm(x)
        if norm==0:
            print('Attention norm(x)=0, division par 10e-6')
            b = x/(10e-6)
        else :
            b = x/norm
        x = np.matmul(A,b)
        lam = np.matmul(b.T,x)
    print(f'convergence en {it} iteration')
    print('lambda =',lam)
    return lam,b

F = np.array([[1,1,0.5],[1,1,0.25],[0.5,0.25,2]])
x0 = np.array([1,0,0])
pui_itere(F,x0,0.01)

def inverse(A,x0,eps):
    x = x0
    it = 0
    lam_, lam = 1,0
    while np.abs(lam-lam_)>eps:
        