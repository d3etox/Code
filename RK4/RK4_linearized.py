import numpy as np
import matplotlib.pyplot as plt
import time
import sys
import math as mt

# Constante modifiable :
Tnode = np.array([50,60,55,45,65,40,0,0,0,0,0,0,12,6,0,-6,18,24,30,35,40,25,20,15,20,28,36,12,4,-2]) # Temperature des noeuds initaux en degré
Tnode_0_deg = Tnode + 40
h = 0.5 # pas de temps de calcul
Np = 300 # pas de temps d'écriture
N = 84600 # temps total en secondes
Cp = 150 # Capacitance entre 2 noeuds (pour tous les couples de noeuds)
Conduc = 50 # Conductance entre 2 noeuds (pour tous les couples de noeuds)
boundary = 6 # Indice du noeud boundary
i_pui = 24 # Indice du noeud auquel on ajoute de la puissance
pui = 2 # Puissance en Watt     

# Constante non modifiable :
node = np.array([0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29])
nn = len(node) # number Node
Tnode_0 = Tnode_0_deg + 273.15 # Temperature des noeuds initiaux en Kelvin
boltz = 5.67*10**(-8) # constante de Boltzman
Nh = N / (60 * 60) # nombre d'heure
Ne = int(N/Np) # nombre d'iteration
eps = 1e-6 

#Capacitance du noeud i
MC = np.full(nn,Cp,dtype=np.float64)
MC[boundary] = 1e+12

#Matrice de conductance :
GL = np.zeros((nn,nn))
GL[0, 1], GL[0, 2], GL[0, 3], GL[0, 4], GL[1, 2], GL[1, 4], GL[1, 5], GL[2, 3], GL[2, 5], GL[3, 4], GL[3, 5], GL[4, 5], \
GL[6, 7], GL[6, 8], GL[6, 9], GL[6, 10], GL[7, 8], GL[7, 10], GL[7, 11], GL[8, 9], GL[8, 11], GL[9, 10], GL[9, 11], \
GL[10, 11], GL[12, 13], GL[12, 14], GL[12, 15], GL[12, 16], GL[13, 14], GL[13, 16], GL[13, 17], GL[14, 15], GL[14, 17], \
GL[15, 16], GL[15, 17], GL[16, 17], GL[18, 19], GL[18, 20], GL[18, 21], GL[18, 22], GL[19, 20], GL[19, 22], GL[19, 23], \
GL[20, 21], GL[20, 23], GL[21, 22], GL[21, 23], GL[22, 23], GL[24, 25], GL[24, 26], GL[24, 27], GL[24, 28], GL[25, 26], \
GL[25, 28], GL[25, 29], GL[26, 27], GL[26, 29], GL[27, 28], GL[27, 29], GL[28, 29] = np.full(60,Conduc)
for i in node :
    for j in node :
        GL[j,i] = GL[i,j]
    
#Matrice d'échange radiatif :
GR = np.zeros((nn,nn))
GR[2, 6], GR[3, 6], GR[2, 7], GR[2, 10], GR[3, 10], GR[2, 11], GR[3, 11], GR[2, 12], GR[3, 12], GR[6, 12], GR[9, 12], \
GR[10, 12], GR[2, 13], GR[3, 13], GR[4, 13], GR[5, 13], GR[9, 13], GR[10, 13], GR[11, 13], GR[2, 14], GR[3, 14], \
GR[5, 14], GR[9, 14], GR[10, 14], GR[11, 14], GR[2, 18], GR[3, 18], GR[5, 18], GR[8, 18], GR[9, 18], GR[10, 18], \
GR[11, 18], GR[13, 18], GR[14, 18], GR[2, 19], GR[5, 19], GR[11, 19], GR[13, 21], GR[14, 21], GR[17, 21], GR[2, 22], \
GR[3, 22], GR[5, 22], GR[8, 22], GR[9, 22], GR[10, 22], GR[11, 22], GR[13, 22], GR[14, 22], GR[17, 22], GR[0, 25], \
GR[2, 25], GR[3, 25], GR[6, 25], GR[9, 25], GR[10, 25], GR[12, 25], GR[13, 25], GR[18, 25], GR[22, 25], GR[6, 26], \
GR[9, 26], GR[18, 26], GR[22, 26], GR[12, 27], GR[13, 27], GR[14, 27], GR[0, 28], GR[2, 28], GR[3, 28], GR[12, 28], \
GR[13, 28], GR[14, 28], GR[2, 29], GR[3, 29], GR[6, 29], GR[9, 29], GR[10, 29], GR[12, 29], GR[13, 29], GR[14, 29], \
GR[18, 29], GR[21, 29], GR[22, 29] = (10**(-6))*np.array((80,4.5833,957.52,2993.88,324.17,160.83,10.83,39,540,2.33,56.67,51.25,135.75,2941.09,75,324.375,79.5,218.75,53.75,105.45,214.69,20.63,265,197.5,77.08,507.8,62.81,64.69,430.83,47.92,22.97,599.58,38.25,107.7,216.56,90,12.5,5.625,182.81,22.5,750.075,146.25,570,55.4167,18.33,27.33,366.165,210.75,543.188,90,36.7708,152.125,151.354,147.291,42.7084,46.3543,5.667,18.1667,6.9375,3.333,71.0418,65.8336,46,15.1667,34.2708,5.52083,33.333,26.3541,33.2291,207.503,74.5834,45.7292,41.9792,61.4586,125.728,42.5,148.646,40.9375,43.5625,41.9792,65.2503,64,1.25,38.4792))
for i in node :
    for j in node :
        GR[j,i] = GR[i,j]

#Puissance dissipée au noeud i
Q = np.zeros(nn)
Q[i_pui] = pui

GR_bis = GR*boltz

def ODE(t,y): #equation de la chaleur seulement avec termes conductif, radiatif et avec une puissance supplémentaire
    cond = np.zeros(nn)
    rad = np.zeros(nn)
    pui = np.zeros(nn)
    result = np.zeros(nn)
    for i in node :
        cond[i] = np.sum(GL[i] * (y - y[i]))
        rad[i] = np.sum(boltz * GR[i] * (y**4 - y[i]**4))
        pui[i] = Q[i]
        result[i]=(cond[i] + rad[i] + pui[i]) / MC[i]
    return result

def RungeKutta4(n):
    TnodeRK4 = np.zeros((nn,Ne+1))
    TnodeRK4[:,0] = Tnode_0_deg.copy()
    Tnodebis = Tnode_0.copy()
    t=0
    start_time = time.time()
    while (t < n ):
        t+=h
        k1 = ODE(t , Tnodebis)
        k2 = ODE(t + 0.5 * h , Tnodebis + h * k1 * 0.5)
        k3 = ODE(t + 0.5 * h , Tnodebis + h * k2 * 0.5)
        k4 = ODE(t + h , Tnodebis + h * k3)
        Tnodebis += (h / 6.) * (k1 + 2 * k2 + 2 * k3 + k4)
        if (abs(t % Np) < eps) :  #Np = pas d'écriture 
            TnodeRK4[:, mt.ceil(t/Np)] = Tnodebis - 273.15
        # Progression et temps estimé restant
        elapsed_time = time.time() - start_time
        progress = (t / n) * 100
        time_per_iteration = elapsed_time / (t / h)
        remaining_iterations = (n - t) / h
        estimated_time_remaining = time_per_iteration * remaining_iterations
        minutes, seconds = divmod(int(estimated_time_remaining), 60)
        sys.stdout.write('\rProgression: [{:.0f}%] Temps restant : \
                {} minutes {} secondes'.format(progress, minutes, seconds))
        sys.stdout.flush()
    return TnodeRK4

def algo(n): #l'équation s'écrit donc ay' + by = c 
    start_time = time.time()
    TnodeAlgo = np.zeros((nn,Ne+1))
    TnodeAlgo[:,0] = Tnode_0_deg.copy()
    Tnew = Tnode_0.copy()
    Tnode_bis = Tnode_0.copy()
    T0 = np.zeros((nn,nn))
    PGL = np.zeros(nn)
    PGR = np.zeros(nn)
    Pfin = np.zeros(nn)
    err_rad = np.zeros((nn,nn))
    t=0
    while (t < n ):
        t+=h
        # DL de T**4 en T0 = (T[i] + T[j])/2 pour tout i,j
        T0 = (Tnode_bis[:, np.newaxis] + Tnode_bis[np.newaxis, :]) / 2
        a = MC
        c = np.array([np.sum((GL[i] + GR[i] * boltz * 4 * T0[i]**3)*Tnode_bis) + Q[i]*np.abs(np.cos(t)) for i in node])
        b = np.array([np.sum(GL[i] + GR[i] * boltz * 4 * T0[i]**3) for i in node])
        Cste = Tnode_0 - c/b
        Tnode_bis = c/b + Cste * np.exp(-(b/a)*t)  #y(t) = c/b + Cste*exp(-b/a*t)
        if (np.abs(t % Np) < eps) :
            TnodeAlgo[:, mt.ceil(t/Np)] = Tnode_bis - 273.15
        # Progression et temps estimé restant
        elapsed_time = time.time() - start_time
        progress = (t / n) * 100
        time_per_iteration = elapsed_time / (t / h)
        remaining_iterations = (n - t) / h
        estimated_time_remaining = time_per_iteration * remaining_iterations
        minutes, seconds = divmod(int(estimated_time_remaining), 60)
        sys.stdout.write('\rProgression: [{:.0f}%] Temps restant : {} minutes {} secondes'.format(progress, minutes, seconds))
        sys.stdout.flush()
    return TnodeAlgo

start_time_RK4 = time.time()    
solution_RK4= RungeKutta4(N)
end_time_RK4 = time.time()

start_time_algo = time.time()
solution_algo = algo(N)
end_time_algo = time.time()

# Calcul du temps d'exécution pour chaque méthode
temps_execution_algo = end_time_algo - start_time_algo
temps_execution_RK4 = end_time_RK4 - start_time_RK4

# Convertir le temps d'exécution en minutes et secondes pour la méthode linéarisé
minutes_algo, secondes_algo = divmod(temps_execution_algo, 60)
# Convertir le temps d'exécution en minutes et secondes pour la méthode RungeKutta4
minutes_RK4, secondes_RK4 = divmod(temps_execution_RK4, 60)

# Afficher le temps de calcul pour chaque méthode
print(f"Temps d'exécution pour la méthode linéaire : {int(minutes_algo)} minutes {int(secondes_algo)} secondes")
print(f"Temps d'exécution pour la méthode RungeKutta4 : {int(minutes_RK4)} minutes {int(secondes_RK4)} secondes")

# Calculer l'erreur maximale entre les deux méthodes pour chaque nœud
erreur_max_par_noeud = np.max(np.abs(solution_RK4 - solution_algo) , axis=1)
erreur_mean_par_noeud = np.mean(np.abs(solution_RK4 - solution_algo) , axis=1)
for i, erreur_max in enumerate(erreur_max_par_noeud):
    print(f"Erreur (maximale,moyenne) pour le nœud {i+1} : {np.round(erreur_max,2), np.round(erreur_mean_par_noeud[i],2)} au temps : {int((Np*h))*np.argmax(np.abs(solution_RK4[i] - solution_algo[i]))}")
print(f"Erreur maximal moyenne pour tous les noeuds : {np.mean(erreur_max_par_noeud)}")

# Temps total en heure
temps_heure = np.linspace(0, Nh, Ne+1)

# Créer deux sous-graphiques verticaux
plt.subplot(2, 1, 1)  # 2 lignes, 1 colonne, premier sous-graphique
# Tracer les températures de chaque nœud pour chaque méthode (Algo)
for i in range(0,25,6) :
    plt.plot(temps_heure, solution_algo[i], label=f"Node {i+1}", linestyle='-')

# Étiquettes et légendes
plt.xlabel('Time (hours)')
plt.ylabel('Temperature')
plt.title('Temperature evolution for each node - linearized'"      Temps d'exécution :"f"{int(minutes_algo)} minutes {int(secondes_algo)} secondes")
plt.legend(loc='upper right')

plt.subplot(2, 1, 2)  # 2 lignes, 1 colonne, deuxième sous-graphique
# Tracer les températures de chaque nœud pour chaque méthode (RK4)
for i in range(0,25,6):
    plt.plot(temps_heure, solution_RK4[i], label=f"Node {i+1}", linestyle='-')

# Ajouter les erreurs maximales sur le graphique RK4
for i in range(0,25,6):
    erreur_max_index = np.argmax(np.abs(solution_algo[i] - solution_RK4[i]))
    plt.plot(temps_heure[erreur_max_index], solution_RK4[i, erreur_max_index], 'ro')
    plt.text(temps_heure[erreur_max_index], solution_RK4[i, erreur_max_index], f'Err {i+1}: {erreur_max_par_noeud[i]:.2f}', fontsize=12, color='red')

# Étiquettes et légendes
plt.xlabel('Time (hours)')
plt.ylabel('Temperature')
plt.title('Temperature evolution for each node - RK4'"             Temps d'exécution :"f"{int(minutes_RK4)} minutes {int(secondes_RK4)} secondes")
plt.legend(loc='upper right')
plt.show()