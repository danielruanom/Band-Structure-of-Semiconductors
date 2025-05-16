import numpy as np
import h5py
import matplotlib.pyplot as plt
import os
from tqdm import tqdm

#me ubico en la carpeta donde se encuentra el archivo hdf5
os.chdir('Ek_GaP/')

N_E = 1000 #defino el número de puntos en el rango de energia
E = np.linspace(-13, 12, N_E)  #defino el rango de energia
E = E.reshape(1, N_E) #reshape para que sea un vector fila
sigma = 0.05  #defino el ancho de la gaussiana
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

E_plot = np.linspace(-13, 12, N_E) #defino el rango de energia para el plot
plt.plot(E_plot, rho, label='Total DOS')
plt.xlabel('Energy (eV)')
plt.ylabel('Density of States')
plt.ylim(0, 2.5)
plt.grid()
plt.show()


