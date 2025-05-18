import numpy as np
import h5py
import matplotlib.pyplot as plt
import os
from tqdm import tqdm

#me ubico en la carpeta donde se encuentra el archivo hdf5
os.chdir('Ek/')

N_E = 10000 #defino el número de puntos en el rango de energia
E = np.linspace(-15, 15, N_E)  #defino el rango de energia
E = E.reshape(1, N_E) #reshape para que sea un vector fila
sigma = 0.1  #defino el ancho de la gaussiana
N_k = len(os.listdir('.')) #defino el número de k-points


rho_nu= np.zeros((10, N_E)) #defino la densidad de estados
for file in tqdm(os.listdir('.')):
    #para un k dado
    #importo el archivo hdf5
    with h5py.File(file, 'r') as f:
        # Get the first group
        eigenEnergy = f['eigE'][:]
        eigenV2 = f['eigv2'][:]

    for mu in range(10):
        ev2 = eigenV2[:, mu] #eigenvector
        ev2 = ev2.reshape(10,1)
        M = np.matmul(ev2, np.exp(-(E - eigenEnergy[mu])**2 / (2*sigma**2)))
        rho_nu += 1/np.sqrt(2*np.pi*sigma**2) *  M #calculo la densidad de estados
rho_nu /= N_k #normalizo la densidad de estados
rho = np.sum(rho_nu, axis=0) #suma de la densidad de estados

os.chdir('..')
with h5py.File('dos.h5', 'w') as f:
    f.create_dataset('rho', data=rho)
    f.create_dataset('E', data=E[0])

markers = ['o', 's', 'D', '^', 'v', '<', '>', 'p', '*', 'h']
# Mismos grosores para todas las curvas
linewidth = 1.8

# Separación razonable entre marcadores (cada 500 puntos aprox)

# Etiquetas para la leyenda
labels = [
    r'$s_a$', r'$s_c$', r'$X_a$', r'$Y_a$', r'$Z_a$',
    r'$X_c$', r'$Y_c$', r'$Z_c$', r'$s_a^*$', r'$s_c^*$'
]

plt.figure(figsize=(10, 6))
# Graficar DOS total (sin marcador)
plt.plot(E[0], rho, 'k', label='Total DoS', linewidth=1.8)

base_markevery = 500

# Graficar proyecciones con marcadores distintos y desfasados
for i in range(10):
    plt.plot(
        E[0], rho_nu[i],
        label=labels[i],
        linewidth=1.8,
        marker=markers[i],
        markevery=(i * 200, base_markevery),  # desfase para cada curva
        markersize=6,
        alpha=0.7
    )

plt.xlabel("Energy (eV)")
plt.ylabel(r"Density of States $(eV-atom)^{-1}$")
plt.title("Density of States (DoS) of GaP")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig('dos.png', dpi=300)
plt.show()


