import numpy as np
import matplotlib.pyplot as plt

x = 2
y = x
y += 2
print('x =', x)
print('y =', y)

def f(x):
    return (np.exp(-x)-x)

def fprime(x):
    return (-np.exp(-x)-1)

def Newton(f,df,x0,tol,max_iter):
    x = x0
    for i in range (max_iter):
        fx = f(x)
        dfx = df(x)
        if dfx == 0 :
            raise ValueError("Dérivée nulle, échec")
        xnew = x - fx/dfx 
        if np.abs(xnew - x) < tol :
            return xnew
    return xnew

x_eq = Newton(f,fprime,0.5,10e-7,1000)
print(x_eq)
x = np.linspace(-1,1,500)
y = f(x)
plt.plot(x,f(x))
plt.show()