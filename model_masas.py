import sympy as sp
import numpy as np
import matplotlib.pyplot as plt
import cmath 

hbarc2=(197*10**(-5))**2    #MeV*A^2   
me=0.511                #MeV

for cond in {'Si'}:#,C'Sn','SiC','GaAs','AlP','GaP','AlAs','GaSb','AlSb','InP','InAs','InSb','ZnSe','ZnTe'}: # 'GaAs',

    if cond == 'Si':
        # values for GaAs
        E_sa=-4.2000
        E_sc=-4.2000
        E_ss=-8.300/4.0
        E_sapc=5.7292/4.0
        E_xx=1.7150/4.0
        E_pa=1.7150
        E_pc=1.7150
        E_pasc=-5.7292/4.0
        E_soapc=5.3749/4.0
        E_xy=4.5750/4.0
        E_soa=6.6850
        E_soc=6.6850
        E_pasoc=-5.3749/4.0
        Delta_c=0.00
        delta_c=Delta_c/3.
        Delta_a=0.00
        delta_a=Delta_a/3. 

        a=5.431

        p=6

    if cond == 'Ge':
        # values for GaAs
        E_sa=-5.8800
        E_sc=-5.8800
        E_ss=-6.7800/4.0
        E_sapc=5.4649/4.0
        E_xx=1.6100/4.0
        E_pa=1.6100
        E_pc=1.6100
        E_pasc=-5.4649/4.0
        E_soapc=5.2191/4.0
        E_xy=4.9000/4.0
        E_soa=6.3900
        E_soc=6.3900
        E_pasoc=-5.2191/4.0
        Delta_c=0.00
        delta_c=Delta_c/3.
        Delta_a=0.00
        delta_a=Delta_a/3. 

        a=5.64

        p=4

    # dimensions hamiltonian
    E_dim=10

    m_LG=400
    m_GX= 400
    Ka_x = np.concatenate([
      np.linspace(np.pi/a, 0, m_LG, endpoint=False),               # de L a Gamma en x
      np.linspace(0, 2*np.pi/a, m_GX),             # de Gamma a X en x
    ])

    Ka_y = np.concatenate([
       np.linspace(np.pi/a, 0, m_LG, endpoint=False),              # de L a Gamma en y
       np.linspace(0, 0, m_GX),                    # de Gamma a X en y
    ])

    Ka_z = np.concatenate([
        np.linspace(np.pi/a, 0, m_LG, endpoint=False),             # de L a Gamma en z
        np.linspace(0, 0, m_GX),                   # de Gamma a X en z
    ])          

    n=m_LG+m_GX

    k_path=np.zeros((n))
    ka = np.vstack([Ka_x, Ka_y, Ka_z]).T  # shape (n, 3)

    # Arrays to store eigenvalues
    E_eigenvalues=np.zeros((n,E_dim), dtype=complex)
    E_eigenvectors_valencia=np.zeros((n,E_dim), dtype=complex)
    E_eigenvectors_conduccion=np.zeros((n,E_dim), dtype=complex)
    
    dE_valencia=np.zeros((n))
    d2E_valencia=np.zeros((n))
    m_valencia=np.zeros((n))

    dE_conducion=np.zeros((n))
    d2E_conducion=np.zeros((n))
    m_conducion=np.zeros((n))

        # Define the energy matrix
    Ea=np.zeros((E_dim,E_dim), dtype=complex) #sea n el número de eigenvalues


    def E_ev(Kx,Ky,Kz):
        
        Ka=np.array([Kx,Ky,Kz])        
         # Define the d vectors
        d1=np.array([0,a/2.,a/2.])
        d2=np.array([a/2.,0.,a/2.])
        d3=np.array([a/2.,a/2.,0.])

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
    


    for j in range(n): #hay que modificar el for para los rangos que miramos
        
        Ka=np.array([Ka_x[j],Ka_y[j],Ka_z[j]])  

        Ea=E_ev(Ka_x[j],Ka_y[j],Ka_z[j])

        # Store the eigenvalues for this value of k
        eigenvalues, eigenvectors = np.linalg.eigh(Ea)
        E_eigenvalues[j, :] = eigenvalues   
        E_eigenvectors_valencia[j,:] = np.abs(eigenvectors[:,3])  #revisar si asi se cogen bien o hay que poner traspuesta o algo
        E_eigenvectors_conduccion[j,:] = np.abs(eigenvectors[:,4])

        #para masa efectiva
        if j!=0:
            dk = np.linalg.norm(ka[j] - ka[j - 1])
            k_path[j] = k_path[j - 1] + dk
      
    capa_valencia=E_eigenvalues[:,3].real #la de menor energia en el gap
    capa_conducion=E_eigenvalues[:,4].real #la de mayor energia en el gap  

    # Find the index and energy value of the minimum in the conduction band (capa_conducion)
    min_conducion_index = np.argmin(capa_conducion)
    min_conducion_energy = capa_conducion[min_conducion_index]
    
    # Find the index and energy value of the maximum in the valence band (capa_valencia)
    max_valencia_index = np.argmax(capa_valencia)
    max_valencia_energy = capa_valencia[max_valencia_index]

    energy_gap=abs(max_valencia_energy-min_conducion_energy)
    
    ### MASAS EFECTIVAS ###

    
    #funcion para la derivada^2 para la masa efectiva
    def compute_effective_d2E_tensor(H_k_function, k0, band_index, delta=1e-3):
        """
        Computes the 3x3 effective mass tensor at k0 for a given band.
        """

        d2E_tensor = np.zeros((3, 3))
        directions = np.identity(3)

        for i in range(3):
            for j in range(3):

                if i == j:
                    # Second derivative ∂²E/∂k_i²
                    E = np.linalg.eigvalsh(H_k_function(k0))[band_index]
                    E_plus = np.linalg.eigvalsh(H_k_function(k0 + delta * directions[i]))[band_index]
                    E_minus = np.linalg.eigvalsh(H_k_function(k0 - delta * directions[i]))[band_index]
                    d2E_tensor[i, i] = (E_plus - 2 * E + E_minus) / (delta ** 2)
                else:
                    # Mixed second derivative (symmetric approximation)
                    E_pp = np.linalg.eigvalsh(H_k_function(k0 + delta * (directions[i] + directions[j])))[band_index]
                    E_pm = np.linalg.eigvalsh(H_k_function(k0 + delta * (directions[i] - directions[j])))[band_index]
                    E_mp = np.linalg.eigvalsh(H_k_function(k0 - delta * (directions[i] - directions[j])))[band_index]
                    E_mm = np.linalg.eigvalsh(H_k_function(k0 - delta * (directions[i] + directions[j])))[band_index]
                    d2E_tensor[i, j] = (E_pp - E_pm - E_mp + E_mm) / (4 * delta ** 2)

        # Invert to get effective mass tensor (if you include hbar^2, do it here)
        return np.linalg.inv(d2E_tensor)
    

    # Get the nth eigenvalue at (kx, ky, kz)
    def E_n(kx,ky,kz, band_index):
        H = E_ev(kx,ky,kz)
        evals = np.linalg.eigh(H)[0]  # sorted eigenvalues
        return np.real(evals[band_index])  # or take absolute value if needed
    

    # Finite difference Hessian
    def effective_inv_mass_tensor(k0, band_index, delta=1e-4):
        k0 = np.array(k0)
        dim = 3
        hessian = np.zeros((dim, dim))

        for i in range(dim):
            for j in range(dim):
                # Base vectors
                dk_i = np.zeros(3)
                dk_j = np.zeros(3)
                dk_i[i] = delta
                dk_j[j] = delta

                fpp = E_n(*(k0 + dk_i + dk_j), band_index)
                fpm = E_n(*(k0 + dk_i - dk_j), band_index)
                fmp = E_n(*(k0 - dk_i + dk_j), band_index)
                fmm = E_n(*(k0 - dk_i - dk_j), band_index)

                hessian[i, j] = (fpp - fpm - fmp + fmm) / (4 * delta**2)

        # Effective mass tensor: m_eff^{-1} = (1 / ħ²) * Hessian
        return hessian 
 
    pSi=6

    aux1_v=effective_inv_mass_tensor(ka[max_valencia_index,:],4)
    auxm1_v,auxm2_v,auxm3_v=hbarc2/me/10**(-6)/np.linalg.eigvalsh(aux1_v)
    print(f"c  auxm1_v={auxm1_v}, auxm2_v={auxm2_v}, auxm3_v={auxm3_v}, ")    
    auxm_dos_v = p**(2/3) * (np.abs(auxm1_v * auxm2_v * auxm3_v))**(1/3) 
    print(f"c m_dos_v_{cond}={auxm_dos_v}"+r" $m_e$")

    ev=np.linalg.eigh(aux1_v)[1]
    print(ev)

    aux1_c=effective_inv_mass_tensor(ka[min_conducion_index,:],4)
    auxm1_c,auxm2_c,auxm3_c=hbarc2/me/10**(-6)/np.linalg.eigvalsh(aux1_c)
    print(f"v auxm1_c={auxm1_c}, auxm2_c={auxm2_c}, auxm3_v={auxm3_c}, ")    
    auxm_dos_c = p**(2/3) * (np.abs(auxm1_c * auxm2_c * auxm3_c))**(1/3) 
    print(f"v m_dos_c_{cond}={auxm_dos_c}"+r" $m_e$")
    '''
    d2E_inv_v = compute_effective_d2E_tensor(E_ev, ka[max_valencia_index,:], band_index=3)
    d2E_inv_c = compute_effective_d2E_tensor(E_ev, ka[min_conducion_index,:], band_index=4)
    d2E_eigvals_v = np.linalg.eigvalsh(d2E_inv_v)
    d2E_eigvals_c = np.linalg.eigvalsh(d2E_inv_c)
    m1_v, m2_v, m3_v = hbarc2/me/d2E_eigvals_v/10**(-6)
    m1_c, m2_c, m3_c = hbarc2/me/d2E_eigvals_c/10**(-6)
    print(f" m1_v={m1_v}, m2_v={m2_v}, m3_v={m3_v}, ")
    print(f" m1_c={m1_c}, m2_v={m2_c}, m3_v={m3_c}, ")
    if cond=='Si':
        p = 6  # Si conduction band
    if cond=='Ge':
        p=4
    m_dos_v = p**(2/3) * (np.abs(m1_v * m2_v * m3_v))**(1/3)
    m_dos_c = p**(2/3) * (np.abs(m1_c * m2_c * m3_c))**(1/3)
    print(f"m_dos_v_{cond}={m_dos_v}"+r" $m_e$")
    print(f"m_dos_c_{cond}={m_dos_c}"+r" $m_e$")

    dE_valencia=np.gradient(E_eigenvalues[:,3],k_path)
    d2E_valencia=np.gradient(dE_valencia,k_path)
    #mask_v = np.isfinite(d2E_valencia)
    #d2E_valencia_clean = d2E_valencia[mask_v]
    #k_path_v_clean = k_path[mask_v]
    m_valencia=hbarc2/me/d2E_valencia/10**(-6)#está en unidades de m_e

    dE_conducion=np.gradient(E_eigenvalues[:,4],k_path)
    d2E_conducion=np.gradient(dE_conducion,k_path)
    #mask_c = np.isfinite(d2E_conducion)
    #d2E_conducion_clean = d2E_conducion[mask_c]
    #k_path_c_clean = k_path[mask_c]
    m_conducion=hbarc2/me/d2E_conducion/10**(-6) #está en unidades de m_e

    ###



    # Plot 2)     masas efectivas
    k_norm_2=np.linspace(-1,1,m_LG+m_GX)    
    plt.plot(k_norm_2, m_valencia[0:m_LG+m_GX],color='r', linewidth=1, label='valencia')  # Take real part
    plt.plot(k_norm_2, m_conducion[0:m_LG+m_GX],color='b', linewidth=1, label='conduccion')  # Take real part

    # **Set x and y limits**
    plt.xlim(-1, 1)  # Set x-axis limits
    #plt.ylim(-15, 15)  # Adjust y-axis as needed

    # Add vertical lines
    plt.axvline(x=k_norm_2[m_LG-1], color='k', linestyle='--', linewidth=1)   # L -> Γ
    plt.axvline(x=k_norm_2[m_LG+m_GX-1], color='k', linestyle='--', linewidth=1) # Γ -> X

    # Remove x-axis numbers
    plt.xticks([], [])  # Empty labels
    # Add custom letters at certain x positions
    plt.text(-1.0, -16.5, 'L', fontsize=12, ha='center', va='bottom')
    plt.text(k_norm_2[m_LG-1], -16.5, r'$\Gamma$', fontsize=12, ha='center', va='bottom')
    plt.text(k_norm_2[m_LG+m_GX-1], -16.5, 'X', fontsize=12, ha='center', va='bottom')
    # New letters in between
    plt.text((-1.0 + k_norm_2[m_LG-1]) / 2, -16.5, r'$\Lambda$', fontsize=12, ha='center', va='bottom')
    plt.text((k_norm_2[m_LG-1] + k_norm_2[m_LG+m_GX-1]) / 2, -16.5, r'$\Delta$', fontsize=12, ha='center', va='bottom')
    plt.text(-0.925, 15 -0.25, cond, fontsize=12, ha='left', va='top')
    plt.legend()
    #plt.xlabel(r'vector de onda $k$')
    plt.ylabel(r'masa efectiva $(m_e)$')
    plt.gcf().set_size_inches(4.5,4.5)
    plt.tight_layout()
    plt.savefig(f'Figure_mass_{cond}.pdf')
    plt.show()

    print(r'La masa efectiva en la banda de valencia del '+cond+'es '+str(np.abs(m_valencia[max_valencia_index]))+r'$m_e$')
    print(r'La masa efectiva en la banda de conducción del '+cond+'es '+str(np.abs(m_conducion[min_conducion_index]))+r'$m_e$')

    '''