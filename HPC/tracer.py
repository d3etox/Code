import numpy as np
import matplotlib.pyplot as plt

# Lecture des données
data = np.loadtxt('solution_global.dat')
x = data[:,0]
u_num = data[:,1]

# Solution analytique
u_exact = 0.5 * x * (1 - x)

plt.plot(x, u_num, '-', label='Numérique')
plt.plot(x, u_exact, '-', label='Analytique')
plt.xlabel('x')
plt.ylabel('u(x)')
plt.legend()
plt.title('Comparaison solution numérique et analytique')
plt.grid(True)
plt.show()
