import numpy as np
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans, SpectralClustering, MeanShift, Birch
import matplotlib.cm as cm
from mpl_toolkits.mplot3d import Axes3D

# Charger les données
data = np.loadtxt('data_tpmam5inum.csv', delimiter=',')

# Coordonnées
x = data[:, 0]
y = data[:, 1]
z = data[:, 2]

# Surfaces
surface = data[:, 6]
surface_scaled = surface / np.max(surface) * 1000

# Températures en °C
temperatures = data[:, 3:6].mean(axis=1)

# --------------------------
# Paramètres
# --------------------------
k = 500
X = np.column_stack((x, y, z, temperatures))
cmap = plt.cm.hot
norm = plt.Normalize(temperatures.min(), temperatures.max())

# --------------------------
# 1ère figure : KMeans
# --------------------------
kmeans = KMeans(n_clusters=k, random_state=0, n_init="auto").fit(X)
labels_kmeans = kmeans.labels_

fig1 = plt.figure()
ax1 = fig1.add_subplot(111, projection='3d')
cluster_temp_means = []
cluster_errors1 = []

for i in range(k):
    idx = labels_kmeans == i
    area = np.sum(surface[idx])
    cx = np.sum(x[idx] * surface[idx]) / area
    cy = np.sum(y[idx] * surface[idx]) / area
    cz = np.sum(z[idx] * surface[idx]) / area
    temp_mean = np.mean(temperatures[idx])
    cluster_temp_means.append(temp_mean)
    temp_error_mean = np.mean(np.abs(temperatures[idx] - temp_mean))
    cluster_errors1.append(temp_error_mean)
    ax1.scatter(cx, cy, cz,
                c=cmap(norm(temp_mean)),
                s=area / np.max(surface) * 1000,
                marker='o', alpha=0.8)



mappable = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
mappable.set_array([])
cbar = fig1.colorbar(mappable, ax=ax1, shrink=0.5, aspect=10)
cbar.set_label('Température moyenne du cluster (°C)')
ax1.set_xlabel('X'); ax1.set_ylabel('Y'); ax1.set_zlabel('Z')
ax1.set_title('KMeans (taille ~ aire, couleur ~ temp moyenne)')

# --------------------------
# 2ème figure : géométrie complète
# --------------------------
fig2 = plt.figure()
ax2 = fig2.add_subplot(111, projection='3d')
sc = ax2.scatter(x, y, z,
                 c=temperatures, cmap='hot',
                 s=surface_scaled, alpha=0.8)
ax2.set_xlabel('X'); ax2.set_ylabel('Y'); ax2.set_zlabel('Z')
cbar = fig2.colorbar(sc, ax=ax2, shrink=0.5, aspect=10)
cbar.set_label('Température moyenne (°C)')
ax2.set_title('Géométrie complète (tous les nœuds)')

# --------------------------
# 3ème figure : SpectralClustering
# --------------------------
spectral = SpectralClustering(n_clusters=k, affinity='nearest_neighbors', random_state=0)
labels_spectral = spectral.fit_predict(X)

fig3 = plt.figure()
ax3 = fig3.add_subplot(111, projection='3d')
cluster_temp_means = []
cluster_errors2 = []

for i in range(k):
    idx = labels_spectral == i
    area = np.sum(surface[idx])
    cx = np.sum(x[idx] * surface[idx]) / area
    cy = np.sum(y[idx] * surface[idx]) / area
    cz = np.sum(z[idx] * surface[idx]) / area
    temp_mean = np.mean(temperatures[idx])
    cluster_temp_means.append(temp_mean)
    temp_error_mean = np.mean(np.abs(temperatures[idx] - temp_mean))
    cluster_errors2.append(temp_error_mean)
    ax3.scatter(cx, cy, cz,
                c=cmap(norm(temp_mean)),
                s=area / np.max(surface) * 1000,
                marker='o', alpha=0.8)



mappable3 = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
mappable3.set_array([])
cbar3 = fig3.colorbar(mappable3, ax=ax3, shrink=0.5, aspect=10)
cbar3.set_label('Température moyenne du cluster (°C)')
ax3.set_xlabel('X'); ax3.set_ylabel('Y'); ax3.set_zlabel('Z')
ax3.set_title('SpectralClustering (taille ~ aire, couleur ~ temp moyenne)')

# --------------------------
# 4ème figure : MeanShift
# --------------------------
meanshift = MeanShift()
labels_ms = meanshift.fit_predict(X)
unique_ms = np.unique(labels_ms)
fig4 = plt.figure()
ax4 = fig4.add_subplot(111, projection='3d')
cluster_temp_means = []
cluster_errors3 = []

for i in unique_ms:
    idx = labels_ms == i
    area = np.sum(surface[idx])
    cx = np.sum(x[idx] * surface[idx]) / area
    cy = np.sum(y[idx] * surface[idx]) / area
    cz = np.sum(z[idx] * surface[idx]) / area
    temp_mean = np.mean(temperatures[idx])
    cluster_temp_means.append(temp_mean)
    temp_error_mean = np.mean(np.abs(temperatures[idx] - temp_mean))
    cluster_errors3.append(temp_error_mean)
    ax4.scatter(cx, cy, cz,
                c=cmap(norm(temp_mean)),
                s=area / np.max(surface) * 1000,
                marker='o', alpha=0.8)



mappable4 = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
mappable4.set_array([])
cbar4 = fig4.colorbar(mappable4, ax=ax4, shrink=0.5, aspect=10)
cbar4.set_label('Température moyenne du cluster (°C)')
ax4.set_xlabel('X'); ax4.set_ylabel('Y'); ax4.set_zlabel('Z')
ax4.set_title('MeanShift (taille ~ aire, couleur ~ temp moyenne)')

# --------------------------
# 5ème figure : Birch
# --------------------------
birch = Birch(n_clusters=k)
labels_birch = birch.fit_predict(X)
fig5 = plt.figure()
ax5 = fig5.add_subplot(111, projection='3d')
cluster_temp_means = []
cluster_errors4 = []

for i in np.unique(labels_birch):
    idx = labels_birch == i
    area = np.sum(surface[idx])
    cx = np.sum(x[idx] * surface[idx]) / area
    cy = np.sum(y[idx] * surface[idx]) / area
    cz = np.sum(z[idx] * surface[idx]) / area
    temp_mean = np.mean(temperatures[idx])
    cluster_temp_means.append(temp_mean)
    temp_error_mean = np.mean(np.abs(temperatures[idx] - temp_mean))
    cluster_errors4.append(temp_error_mean)
    ax5.scatter(cx, cy, cz,
                c=cmap(norm(temp_mean)),
                s=area / np.max(surface) * 1000,
                marker='o', alpha=0.8)



mappable5 = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
mappable5.set_array([])
cbar5 = fig5.colorbar(mappable5, ax=ax5, shrink=0.5, aspect=10)
cbar5.set_label('Température moyenne du cluster (°C)')
ax5.set_xlabel('X'); ax5.set_ylabel('Y'); ax5.set_zlabel('Z')
ax5.set_title('Birch (taille ~ aire, couleur ~ temp moyenne)')

print(f"KMeans - Erreur moyenne maximale = {max(cluster_errors1):.2f} °C")
print(f"SpectralClustering - Erreur moyenne maximale = {max(cluster_errors2):.2f} °C")
print(f"MeanShift - Erreur moyenne maximale = {max(cluster_errors3):.2f} °C")
print(f"Birch - Erreur moyenne maximale = {max(cluster_errors4):.2f} °C")

plt.show()

