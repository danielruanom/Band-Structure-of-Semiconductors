import numpy as np
from ase.dft.kpoints import monkhorst_pack
import matplotlib.pyplot as plt
import h5py
from tqdm import tqdm

a=1
#vectores de la red reciproca a la FCC (diamante) sin considerar base atomica
b1 = 2*np.pi/a * np.array([-1, 1, 1])
b2 = 2*np.pi/a * np.array([1, -1, 1])
b3 = 2*np.pi/a * np.array([1, 1, -1])

#genero malla homogenea en la red reciproca NxNxN
N = 37

KMESH = monkhorst_pack((N, N, N))
kzb = []
#represento 3d la malla

for k in KMESH:
    kzb.append(k[0]*b1 + k[1]*b2 + k[2]*b3)

kzb = np.array(kzb)

# values for GaP
E_sa=-8.1124
E_sc=-2.1976
E_pa=1.1250
E_pc=4.1150
E_soa=8.5150
E_soc=7.1850
E_ss=-7.4709/4.0
E_sapc=4.2771/4.0
E_xx=2.1516/4.0
E_pasc=-6.3190/4.0
E_soapc=4.6541/4.0
E_xy=5.1369/4.0
E_pasoc=-5.0950/4.0
Delta_c=0.00 
delta_c=Delta_c/3.
Delta_a=0.00
delta_a=Delta_a/3.
# Define the d vectors
d1=np.array([0,a/2.,a/2.])
d2=np.array([a/2.,0.,a/2.])
d3=np.array([a/2.,a/2.,0.])

E_dim=10
Ea=np.zeros((E_dim,E_dim), dtype=complex) #sea n el número de eigenvalues
def E_ev(Ka):
        #Ka=np.array([kx,ky,kz])  

        g0=(1+np.exp(-1j*np.dot(d1,Ka))+np.exp(-1j*np.dot(d2,Ka))+np.exp(-1j*np.dot(d3,Ka)))/4.0
        g1=(1+np.exp(-1j*np.dot(d1,Ka))-np.exp(-1j*np.dot(d2,Ka))-np.exp(-1j*np.dot(d3,Ka)))/4.0
        g2=(1-np.exp(-1j*np.dot(d1,Ka))+np.exp(-1j*np.dot(d2,Ka))-np.exp(-1j*np.dot(d3,Ka)))/4.0
        g3=(1-np.exp(-1j*np.dot(d1,Ka))-np.exp(-1j*np.dot(d2,Ka))+np.exp(-1j*np.dot(d3,Ka)))/4.0

        Ea[0,0]=E_sa; Ea[1,1]=E_sc; Ea[2,2]=E_pa; Ea[3,3]=E_pa; Ea[4,4]=E_pa; Ea[5,5]=E_pc; Ea[6,6]=E_pc; Ea[7,7]=E_pc; Ea[8,8]=E_soa; Ea[9,9]=E_soc #diagonal
        
        Ea[0,1]=4.*E_ss*g0;  Ea[0,5]=4.*E_sapc*g1; Ea[0,6]=4.*E_sapc*g2; Ea[0,7]=4.*E_sapc*g3 #row 1
        Ea[1,0]=4.*E_ss*np.conj(g0); Ea[1,2]=4.*E_pasc*np.conj(g1); Ea[1,3]=4.*E_pasc*np.conj(g2); Ea[1,4]=4.*E_pasc*np.conj(g3) #row 2 
        Ea[2,1]=4.*E_pasc*g1; Ea[2,5]=4.*E_xx*g0; Ea[2,6]=4.*E_xy*g3; Ea[2,7]=4.*E_xy*g2; Ea[2,9]=4.*E_pasoc*g1
        Ea[3,1]=4.*E_pasc*g2; Ea[3,5]=4.*E_xy*g3; Ea[3,6]=4.*E_xx*g0; Ea[3,7]=4.*E_xy*g1; Ea[3,9]=4.*E_pasoc*g2
        Ea[4,1]=4.*E_pasc*g3; Ea[4,5]=4.*E_xy*g2; Ea[4,6]=4.*E_xy*g1; Ea[4,7]=4.*E_xx*g0; Ea[4,9]=4.*E_pasoc*g3
        Ea[5,0]=4.*E_sapc*np.conj(g1); Ea[5,2]=4.*E_xx*np.conj(g0); Ea[5,3]=4.*E_xy*np.conj(g3); Ea[5,4]=4.*E_xy*np.conj(g2); Ea[5,8]=4.*E_soapc*np.conj(g1)
        Ea[6,0]=4.*E_sapc*np.conj(g2); Ea[6,2]=4.*E_xy*np.conj(g3); Ea[6,3]=4.*E_xx*np.conj(g0); Ea[6,4]=4.*E_xy*np.conj(g1); Ea[6,8]=4.*E_soapc*np.conj(g2) 
        Ea[7,0]=4.*E_sapc*np.conj(g3); Ea[7,2]=4.*E_xy*np.conj(g2); Ea[7,3]=4.*E_xy*np.conj(g1); Ea[7,4]=4.*E_xx*np.conj(g0); Ea[7,8]=4.*E_soapc*np.conj(g3)
        Ea[8,5]=4.*E_soapc*g1; Ea[8,6]=4.*E_soapc*g2; Ea[8,7]=4.*E_soapc*g3
        Ea[9,2]=4.*E_pasoc*np.conj(g1); Ea[9,3]=4.*E_pasoc*np.conj(g2); Ea[9,4]=4.*E_pasoc*np.conj(g3);    
    
        return Ea

#creo directorio Ek si no existe
import os
if not os.path.exists('Ek'):
    os.makedirs('Ek')
#calculo eigenvalores y eigenvectores para cada k
i=0
for k in tqdm(kzb):
    E = E_ev(k)
    eigE, eigv = np.linalg.eig(E)
    eigE = np.real(eigE)
    eigv2 = np.real(np.conjugate(eigv)*eigv)
    with h5py.File(f'Ek/k_{i}.h5', 'w') as f:
        f.create_dataset('eigE', data=eigE)
        f.create_dataset('eigv2', data=eigv2)
    i+=1


      