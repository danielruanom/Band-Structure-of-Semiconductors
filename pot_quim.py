import numpy as np
import h5py
import matplotlib.pyplot as plt
import os
from tqdm import tqdm
from scipy.optimize import minimize

def f_FD(E, nu, T): # Distribucion de Fermi-Dirac, E [eV], nu [eV], T [K]
    k = 8.617333262e-5 # eV/K
    return 1/(np.exp((E - nu)/(k*T)) + 1)

def O_min(nu, T, E_c, rho_c, E_v, rho_v): # Funcion a minimizar para encontrar el potencial quimico
    n_c = np.trapz(rho_c * f_FD(E_c, nu, T), E_c)
    n_v = np.trapz(rho_v* (1 - f_FD(E_v, nu, T)), E_v)
    return (n_c - n_v)**2

with h5py.File('dos.h5', 'r') as f:
    E = f['E'][:]
    rho = f['rho'][:]

#divido la densidad de estados en dos partes, valencia y conduccion
rho_c = rho[E>1]
E_c = E[E>1]

rho_v = rho[E<0.5]
E_v = E[E<0.5]
T = 5000
nu0 = 2.349/2
nu = minimize(O_min, nu0, args=(T, E_c, rho_c, E_v, rho_v), method='Nelder-Mead', options={'disp': True}).x[0]
print('Potencial quimico:', nu)
#calculo la densidad de carga
n_c = np.trapz(rho_c * f_FD(E_c, nu, T), E_c)
n_v = np.trapz(rho_v* (1 - f_FD(E_v, nu, T)), E_v)
print('Densidad de carga:', n_c, n_v)

#plt.figure(figsize=(10, 6))
E_c = E[E>nu]
rho_c = rho[E>nu]
E_v = E[E<nu]
rho_v = rho[E<nu]
plt.plot(E_c, rho_c, label=r'$\rho_c$')
#cojo color de esta linea
c1 = plt.gca().lines[-1].get_color()
plt.plot(E_v, rho_v, label=r'$\rho_v$')
c2 = plt.gca().lines[-1].get_color()
plt.plot(E, f_FD(E, nu, T),color = c1, label=r'$f_{FD}$')
plt.plot(E, 1 - f_FD(E, nu, T),color = c2, label=r'$1 - f_{FD}$')
plt.axvline(nu, color='k', linestyle='--', label=r'$\mu$')
plt.xlabel('Energy (eV)')
plt.grid()
plt.legend()
plt.savefig('a.png', dpi=300)
plt.show()

plt.figure(figsize=(10, 6))
# Plot rho_c * f_FD y sombrear
line_c, = plt.plot(E_c, rho_c * f_FD(E_c, nu, T), label=r'$\rho_c \cdot f_{FD}$')
plt.fill_between(E_c, 0, rho_c * f_FD(E_c, nu, T), color=line_c.get_color(), alpha=0.3, label=r'$n_c$')
# Plot rho_v * (1 - f_FD) y sombrear
line_v, = plt.plot(E_v, rho_v * (1 - f_FD(E_v, nu, T)), label=r'$\rho_v \cdot (1 - f_{FD})$')
plt.fill_between(E_v, 0, rho_v * (1 - f_FD(E_v, nu, T)), color=line_v.get_color(), alpha=0.3, label=r'$p_v$')
# Línea del potencial químico
plt.axvline(nu, color='k', linestyle='--', label=r'$\mu$')
plt.xlabel('Energy (eV)')
plt.grid()
plt.legend()
plt.tight_layout()
plt.savefig('b.png', dpi=300)
plt.show()

#estudio como cambia con T
T_list = np.linspace(390, 500, 15)
nu_list = []
nc_list = []
nv_list = []
for T in tqdm(T_list):
    nu = minimize(O_min, nu0, args=(T, E_c, rho_c, E_v, rho_v), method='Nelder-Mead', options={'disp': True}).x[0]
    nu_list.append(nu)
    n_c = np.trapz(rho_c * f_FD(E_c, nu, T), E_c)
    n_v = np.trapz(rho_v* (1 - f_FD(E_v, nu, T)), E_v)
    nc_list.append(n_c)
    nv_list.append(n_v)

nu_list = np.array(nu_list)


def nu_anal(T, a, b): 
    return a + b * T

#fiteo nu_anal a nu_list
from scipy.optimize import curve_fit

# Ajuste a los datos numéricos
popt, pcov = curve_fit(nu_anal, T_list, nu_list)

# Resultado:
a = popt[0]
b = popt[1]
T_ajust = np.linspace(0, 500, 100)
print(a,b)
plt.figure(figsize=(10, 6))
plt.axhline(nu_anal(0,a,b), color='k', linestyle='--', label=r'$E_{gap}/2$')
plt.plot(T_ajust, nu_anal(T_ajust, a, b),'r-', label='Ajuste lineal')
plt.plot(T_list, nu_list,'ko')
plt.xlabel('Temperature (K)')
plt.ylabel(r'$\mu$ (eV)')
plt.legend()
plt.xlim(0, 500)
plt.grid()
plt.savefig('nu.png', dpi=300)
plt.show()