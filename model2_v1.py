import sympy as sp
import numpy as np
import matplotlib.pyplot as plt
import cmath 

hbarc2=(197*10**(-5))**2    #MeV*A^2   
me=0.511                #MeV
A=10**(-10)             #Amstrong
c=3*10**8               #m/s

for cond in {'Si','Ge'}:#,C'Sn','SiC','GaAs','AlP','GaP','AlAs','GaSb','AlSb','InP','InAs','InSb','ZnSe','ZnTe'}: # 'GaAs',

    if cond == 'C':
        # values for GaAs
        E_sa=-4.5450
        E_sc=-4.545
        E_ss=-22.7250/4.0
        E_sapc=15.2206/4.0
        E_xx=3.8400/4.0
        E_pa=3.8400
        E_pc=3.8400
        E_pasc=-15.2206/4.0
        E_soapc=8.2109/4.0
        E_xy=11.6700/4.0
        E_soa=11.3700
        E_soc=11.3700
        E_pasoc=-8.2109/4.0
        Delta_c=0.00
        delta_c=Delta_c/3.
        Delta_a=0.00
        delta_a=Delta_a/3. 

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

    if cond == 'Sn':
        # values for GaAs
        E_sa=-5.6700
        E_sc=-5.6700
        E_ss=-5.6700/4.0
        E_sapc=4.5116/4.0
        E_xx=1.3300/4.0
        E_pa=1.3300
        E_pc=1.3300
        E_pasc=-4.5116/4.0
        E_soapc=5.8939/4.0
        E_xy=4.0800/4.0
        E_soa=5.9000
        E_soc=5.9000
        E_pasoc=-5.8939/4.0
        Delta_c=0.00
        delta_c=Delta_c/3.
        Delta_a=0.00
        delta_a=Delta_a/3. 

    if cond == 'SiC':
        # values for GaAs
        E_sa=-8.4537
        E_sc=-4.8463
        E_ss=-12.4197/4.0
        E_sapc=9.4900/4.0
        E_xx=3.0380/4.0
        E_pa=2.1234
        E_pc=4.3466
        E_pasc=-9.2007/4.0
        E_soapc=8.7138/4.0
        E_xy=5.9216/4.0
        E_soa=9.6534
        E_soc=9.3166
        E_pasoc=-4.4051/4.0
        Delta_c=0.00
        delta_c=Delta_c/3.
        Delta_a=0.00
        delta_a=Delta_a/3. 

    if cond == 'AlP':
        # values for GaAs
        E_sa=-7.8466
        E_sc=-1.2534
        E_ss=-7.4535/4.0
        E_sapc=5.2451/4.0
        E_xx=2.3749/4.0
        E_pa=1.3169
        E_pc=4.2831
        E_pasc=-5.2775/4.0
        E_soapc=5.2508/4.0
        E_xy=4.8378/4.0
        E_soa=8.7069
        E_soc=7.4231
        E_pasoc=-4.6388/4.0
        Delta_c=0.00
        delta_c=Delta_c/3.
        Delta_a=0.00
        delta_a=Delta_a/3.

    if cond == 'AlAs':
        # values for AlAs
        E_sa=-7.5273
        E_sc=-1.1627
        E_pa=0.9833
        E_pc=3.5867
        E_soa=7.4833
        E_soc=6.7267
        E_ss=-6.6642/4.0
        E_sapc=5.1106/4.0
        E_xx=1.8780/4.0
        E_pasc=-5.4965/4.0
        E_soapc=4.5216/4.0
        E_xy=4.2919/4.0
        E_pasoc=-4.9950/4.0
        Delta_c=0.00 
        delta_c=Delta_c/3.
        Delta_a=0.00
        delta_a=Delta_a/3.

    if cond == 'AlSb':
        # values for AlSb
        E_sa=-6.1714
        E_sc=-2.0716
        E_pa=0.9807
        E_pc=3.0163
        E_soa=6.7607
        E_soc=6.1543
        E_ss=-5.6448/4.0
        E_sapc=4.9121/4.0
        E_xx=1.7199/4.0
        E_pasc=-4.2137/4.0
        E_soapc=4.3662/4.0
        E_xy=3.6648/4.0
        E_pasoc=-3.0739/4.0
        Delta_c=0.00 
        delta_c=Delta_c/3.
        Delta_a=0.00
        delta_a=Delta_a/3.

    if cond == 'GaP':
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

    if cond == 'GaAs':
        # values for GaAs
        E_sa=-8.3431
        E_sc=-2.6569
        E_ss=-6.4513/4.0
        E_sapc=4.48/4.0
        E_xx=1.9546/4.0
        E_pa=1.0414
        E_pc=3.36686
        E_pasc=-5.7839/4.0
        E_soapc=4.8422/4.0
        E_xy=5.0779/4.0
        E_soa=8.5914
        E_soc=6.7386
        E_pasoc=-4.8077/4.0
        Delta_c=0.013
        delta_c=Delta_c/3.
        Delta_a=0.38
        delta_a=Delta_a/3.

    if cond == 'GaSb':
        # values for GaSb
        E_sa=-7.3207
        E_sc=-3.8993
        E_pa=0.8554
        E_pc=2.9146
        E_soa=6.6354
        E_soc=5.9846
        E_ss=-6.1567/4.0
        E_sapc=4.9601/4.0
        E_xx=1.5789/4.0
        E_pasc=-4.6675/4.0
        E_soapc=4.9895/4.0
        E_xy=4.1285/4.0
        E_pasoc=-4.2180/4.0
        Delta_c=0.00 
        delta_c=Delta_c/3.
        Delta_a=0.00
        delta_a=Delta_a/3.

    if cond == 'InP':
        # values for InP
        E_sa=-8.5274
        E_sc=-1.4826
        E_pa=0.8735
        E_pc=4.0465
        E_soa=8.2635
        E_soc=7.0665
        E_ss=-5.3614/4.0
        E_sapc=2.2265/4.0
        E_xx=1.8801/4.0
        E_pasc=-5.5825/4.0
        E_soapc=3.4623/4.0
        E_xy=4.2324/4.0
        E_pasoc=-4.4814/4.0
        Delta_c=0.00 
        delta_c=Delta_c/3.
        Delta_a=0.00
        delta_a=Delta_a/3.

    if cond == 'InAs':
        # values for InAs
        E_sa=-9.5381
        E_sc=-2.7219
        E_pa=0.9099
        E_pc=3.7201
        E_soa=7.4099
        E_soc=6.7401
        E_ss=-5.6052/4.0
        E_sapc=3.0354/4.0
        E_xx=1.8398/4.0
        E_pasc=-5.4389/4.0
        E_soapc=3.3744/4.0
        E_xy=4.4693/4.0
        E_pasoc=-3.9097/4.0
        Delta_c=0.00 
        delta_c=Delta_c/3.
        Delta_a=0.00
        delta_a=Delta_a/3.

    if cond == 'InSb':
        # values for InSb
        E_sa=-8.0157
        E_sc=-3.4643
        E_pa=0.6738
        E_pc=2.9126
        E_soa=6.4530
        E_soc=5.9362
        E_ss=-5.5193/4.0
        E_sapc=3.7880/4.0
        E_xx=1.4018/4.0
        E_pasc=-4.5900/4.0
        E_soapc=3.5666/4.0
        E_xy=3.8761/4.0
        E_pasoc=-3.4048/4.0
        Delta_c=0.00 
        delta_c=Delta_c/3.
        Delta_a=0.00
        delta_a=Delta_a/3.

    if cond == 'ZnSe':
        # values for ZnSe
        E_sa=-11.8383
        E_sc=0.0183
        E_pa=1.5072
        E_pc=5.9928
        E_soa=7.5872
        E_soc=8.9928
        E_ss=-6.2163/4.0
        E_sapc=3.4980/4.0
        E_xx=3.0054/4.0
        E_pasc=-6.3191/4.0
        E_soapc=2.5891/4.0
        E_xy=5.9942/4.0
        E_pasoc=-3.9533/4.0
        Delta_c=0.00 
        delta_c=Delta_c/3.
        Delta_a=0.00
        delta_a=Delta_a/3.

    if cond == 'ZnTe':
        # values for ZnTe
        E_sa=-9.8150
        E_sc=0.9350
        E_pa=1.4834
        E_pc=5.2666
        E_soa=7.0834
        E_soc=8.2666
        E_ss=-6.5765/4.0
        E_sapc=5.9827/4.0
        E_xx=2.7951/4.0
        E_pasc=-5.8199/4.0
        E_soapc=1.3196/4.0
        E_xy=5.4670/4.0
        E_pasoc=0.0000
        Delta_c=0.00 
        delta_c=Delta_c/3.
        Delta_a=0.00
        delta_a=Delta_a/3.


    # Define the lattice constant (este valor no afecta al calculo del gap)
    a=1.
    if cond=='Si':
        a=5.43
    if cond=='Ge':
        a=5.64

    # dimensions hamiltonian
    E_dim=10


    # el rango va a ser de a) k=0 a k=(2pi/a,0,0), y b) de k=0 a k=pi/a(1,1,1)
    # Define the wave vector
    '''
    kminx_a=0
    kmaxx_a=2*np.pi/a
    kx_a=np.linspace(kminx_a,kmaxx_a,m)

    kminx_b=0
    kmaxx_b=np.pi/a
    kx_b=np.linspace(kminx_b,kmaxx_b,m)
    kminy=0
    kmaxy=np.pi/a
    ky=np.linspace(kminy,kmaxy,m)
    kminz=0
    kmaxz=np.pi/a
    kz=np.linspace(kminz,kmaxz,m)
    '''
    #para punto K, 3/4*2*np.pi/a*(1,1,0)

    m_LG=100
    m_GX= 100
    m_XK=30
    m_KG=100
    Ka_x = np.concatenate([
      np.linspace(np.pi/a, 0, m_LG, endpoint=False),               # de L a Gamma en x
      np.linspace(0, 2*np.pi/a, m_GX, endpoint=False),             # de Gamma a X en x
     np.linspace(2*np.pi/a, 1*2*np.pi/a, m_XK, endpoint=False),  # de X a U en x              ###### esta parte hay que revisarla
      np.linspace(3/4*2*np.pi/a, 0, m_KG)          # de K a Gamma en x
    ])

    Ka_y = np.concatenate([
       np.linspace(np.pi/a, 0, m_LG, endpoint=False),              # de L a Gamma en y
       np.linspace(0, 0, m_GX, endpoint=False),                    # de Gamma a X en y
      np.linspace(0, 1/4*2*np.pi/a, m_XK, endpoint=False),         # de X a U en y              ###### esta parte hay que revisarla
      np.linspace(3/4*2*np.pi/a, 0, m_KG)          # de K a Gamma en y
    ])

    Ka_z = np.concatenate([
        np.linspace(np.pi/a, 0, m_LG, endpoint=False),             # de L a Gamma en z
        np.linspace(0, 0, m_GX, endpoint=False),                   # de Gamma a X en z
       np.linspace(0, 1/4*2*np.pi/a, m_XK, endpoint=False),                    # de X a U en z              ###### esta parte hay que revisarla
       np.linspace(0, 0, m_KG)                     # de K a Gamma en z
    ])          

    n=m_LG+m_GX+m_XK+m_KG

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


    E_eigenvectors_valencia_SO=np.zeros((n,2*E_dim), dtype=complex)
    E_eigenvectors_conduccion_SO=np.zeros((n,2*E_dim), dtype=complex)

    # matrices para el calculo de SO
    Hso_a=np.zeros((6,6), dtype=complex) #sea n el número de eigenvalues
    Hso_c=np.zeros((6,6), dtype=complex) #sea n el número de eigenvalues

        # Define the energy matrix
    Ea=np.zeros((E_dim,E_dim), dtype=complex) #sea n el número de eigenvalues
    E_SO_eigenvalues=np.zeros((n,2*E_dim), dtype=complex) #para SO


     # Define the d vectors
    d1=np.array([0,a/2.,a/2.])
    d2=np.array([a/2.,0.,a/2.])
    d3=np.array([a/2.,a/2.,0.])


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




    for j in range(n): #hay que modificar el for para los rangos que miramos

        
        Ka=np.array([Ka_x[j],Ka_y[j],Ka_z[j]])  

        Ea=E_ev(Ka)

        # Store the eigenvalues for this value of k
        eigenvalues, eigenvectors = np.linalg.eigh(Ea)

        E_eigenvalues[j, :] = eigenvalues
            
        E_eigenvectors_valencia[j,:] = np.abs(eigenvectors[:,3])  #revisar si asi se cogen bien o hay que poner traspuesta o algo
        E_eigenvectors_conduccion[j,:] = np.abs(eigenvectors[:,4])


        #para masa efectiva
        if j!=0:
            dk = np.linalg.norm(ka[j] - ka[j - 1])
            k_path[j] = k_path[j - 1] + dk


        ###aqui el calculo de la parte de SO

        #voy a hacer un estudio de eigenvectors para GaAs en banda valencia y conduccion
        if cond == 'GaAs':
            H=np.zeros((2*E_dim,2*E_dim), dtype=complex) #sea n el número de eigenvalues    
            H[:10, :10] = Ea  # Top-left block
            H[10:20, 10:20] = Ea  # Bottom-right block

            #ahora, la parte de S-0
            Hso_a[0,1] =-1j*delta_a;        Hso_a[0,5] =delta_a
            Hso_a[1,0] =1j*delta_a;         Hso_a[1,5] =-1j*delta_a              
            Hso_a[2,3] =-delta_a;           Hso_a[2,4] =1j*delta_a
            Hso_a[3,2] =-delta_a;           Hso_a[3,4] =1j*delta_a
            Hso_a[4,2] =-1j*delta_a;        Hso_a[4,3] =-1j*delta_a
            Hso_a[5,0] =delta_a;            Hso_a[5,1] =1j*delta_a

            Hso_c[0,1] =-1j*delta_c;        Hso_c[0,5] =delta_c
            Hso_c[1,0] =1j*delta_c;         Hso_c[1,5] =-1j*delta_c              
            Hso_c[2,3] =-delta_c;           Hso_c[2,4] =1j*delta_c
            Hso_c[3,2] =-delta_c;           Hso_c[3,4] =1j*delta_c
            Hso_c[4,2] =-1j*delta_c;        Hso_c[4,3] =-1j*delta_c
            Hso_c[5,0] =delta_c;            Hso_c[5,1] =1j*delta_c

            # doing a:b means a, a+1, ... b-1
            rows_a = np.r_[2:5, 12:15]  # Selects rows 3-4 and 8-9
            cols_a = np.r_[2:5, 12:15]  # Selects cols 3-4 and 8-9

            rows_c = np.r_[5:8, 15:18]  # Selects cols 3-4 and 8-9
            cols_c = np.r_[5:8, 15:18]  # Selects cols 3-4 and 8-9

            H[np.ix_(rows_a, cols_a)] += Hso_a
            H[np.ix_(rows_c, cols_c)] += Hso_c

            # Store the eigenvalues for this value of k
            eigenvalues_SO, eigenvectors_SO = np.linalg.eig(H)

            E_SO_eigenvalues[j, :] = np.sort(eigenvalues_SO)

            idx_SO = np.argsort(eigenvalues_SO)
            E_eigenvectors_SO = eigenvectors_SO[:, idx_SO] 
            E_eigenvectors_valencia_SO[j,:] = np.abs(E_eigenvectors_SO[:,7].T)  #revisar si asi se cogen bien o hay que poner traspuesta o algo
            E_eigenvectors_conduccion_SO[j,:] = np.abs(E_eigenvectors_SO[:,10].T)

    if cond != 'Sn':        
        capa_valencia=E_eigenvalues[:,3].real #la de menor energia en el gap
        capa_conducion=E_eigenvalues[:,4].real #la de mayor energia en el gap  
    else:
        capa_valencia=E_eigenvalues[:,7].real #la de menor energia en el gap
        capa_conducion=E_eigenvalues[:,8].real #la de mayor energia en el gap  

    # Find the index and energy value of the minimum in the conduction band (capa_conducion)
    min_conducion_index = np.argmin(capa_conducion)
    min_conducion_energy = capa_conducion[min_conducion_index]
    
    # Find the index and energy value of the maximum in the valence band (capa_valencia)
    max_valencia_index = np.argmax(capa_valencia)
    max_valencia_energy = capa_valencia[max_valencia_index]

    energy_gap=abs(max_valencia_energy-min_conducion_energy)

    if min_conducion_index == max_valencia_index:
        print(f"{cond} is a direct semiconductor with a energy gap of {energy_gap:.3f} eV")
    else:
        print(f"{cond} is a indirect semiconductor with a energy gap of {energy_gap:.3f} eV")        


    ### MASAS EFECTIVAS ###
   
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


    # Plot 1)        # this is without the SO coupling for the GaAs
    k_norm=np.linspace(-1,1,n)
    plt.plot(k_norm, E_eigenvalues[:,0].real,color='k', linewidth=1, label='E1')  # Take real part
    plt.plot(k_norm, E_eigenvalues[:,1].real,color='k', linewidth=1, label='E2')  # Take real part
    plt.plot(k_norm, E_eigenvalues[:,2].real,color='k', linewidth=1, label='E3')  # Take real part
    plt.plot(k_norm, E_eigenvalues[:,3].real,color='k', linewidth=1, label='E4')  # Take real part
    plt.plot(k_norm, E_eigenvalues[:,4].real,color='k', linewidth=1, label='E5')  # Take real part
    plt.plot(k_norm, E_eigenvalues[:,5].real,color='k', linewidth=1, label='E6')  # Take real part
    plt.plot(k_norm, E_eigenvalues[:,6].real,color='k', linewidth=1, label='E7')  # Take real part
    plt.plot(k_norm, E_eigenvalues[:,7].real,color='k', linewidth=1, label='E8')  # Take real part
    plt.plot(k_norm, E_eigenvalues[:,8].real,color='k', linewidth=1, label='E9')  # Take real part
    plt.plot(k_norm, E_eigenvalues[:,9].real,color='k', linewidth=1, label='E10')  # Take real part

    # **Set x and y limits**
    plt.xlim(-1, 1)  # Set x-axis limits
    plt.ylim(-15, 15)  # Adjust y-axis as needed

    # Add vertical lines
    plt.axvline(x=k_norm[m_LG-1], color='k', linestyle='--', linewidth=1)   # L -> Γ
    plt.axvline(x=k_norm[m_LG+m_GX-1], color='k', linestyle='--', linewidth=1) # Γ -> X
    plt.axvline(x=k_norm[m_LG+m_GX+m_XK-1], color='k', linestyle='--', linewidth=1) # X -> K

    # Add horizontal lines at the energy values where the min and max energies occur
    plt.axhline(y=min_conducion_energy, color='k', linestyle='--', linewidth=1)
    plt.axhline(y=max_valencia_energy, color='k', linestyle='--', linewidth=1)

    # Shaded area (between the horizontal lines)
    plt.fill_between(k_norm, min_conducion_energy, max_valencia_energy, where=(capa_conducion > min_conducion_energy) & (capa_valencia < max_valencia_energy), 
                 color='gray', alpha=0.3, label='Energy Gap')

    # Remove x-axis numbers
    plt.xticks([], [])  # Empty labels
    # Add custom letters at certain x positions
    plt.text(-1.0, -16.5, 'L', fontsize=12, ha='center', va='bottom')
    plt.text(k_norm[m_LG-1], -16.5, r'$\Gamma$', fontsize=12, ha='center', va='bottom')
    plt.text(k_norm[m_LG+m_GX-1], -16.5, 'X', fontsize=12, ha='center', va='bottom')
    plt.text(k_norm[m_LG+m_GX+m_XK-1], -16.5, 'U, K', fontsize=12, ha='center', va='bottom')
    plt.text(1.0, -16.5, r'$\Gamma$', fontsize=12, ha='center', va='bottom')
    # New letters in between
    plt.text((-1.0 + k_norm[m_LG-1]) / 2, -16.5, r'$\Lambda$', fontsize=12, ha='center', va='bottom')
    plt.text((k_norm[m_LG-1] + k_norm[m_LG+m_GX-1]) / 2, -16.5, r'$\Delta$', fontsize=12, ha='center', va='bottom')
    plt.text((k_norm[m_LG+m_GX+m_XK-1] + 1.0) / 2, -16.5, r'$\Sigma$', fontsize=12, ha='center', va='bottom')
    
    plt.text(-0.925, 15 -0.25, cond, fontsize=12, ha='left', va='top')

    #plt.xlabel(r'vector de onda $k$')
    plt.ylabel('E (eV)')
    plt.gcf().set_size_inches(4.5,4.5)
    plt.tight_layout()
    plt.savefig(f'Figure_1_{cond}.pdf')
    plt.show()



    # Plot 2)     masas efectivas
    k_norm_2=np.linspace(-1,1,m_LG+m_GX)    
    plt.plot(k_norm_2, np.abs(m_valencia[0:m_LG+m_GX]),color='r',marker='o', linewidth=1, label='valencia')  # Take real part
    plt.plot(k_norm_2, np.abs(m_conducion[0:m_LG+m_GX]),color='b', linewidth=1, label='conduccion')  # Take real part

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



### Aquí las figuras de las parte de autofunciones
    if cond in ['C','Si','Ge','Sn','SiC','GaAs','AlP','GaP','AlAs','GaSb','AlSb','InP','InAs','InSb','ZnSe','ZnTe']:
        
        # plot 2) Eigenvectors para las bandas de conduccion y valencia en GaAs
        colors = ['r', 'b', 'g', 'm', 'c', 'y', 'k', 'orange', 'darkred', 'navy']  # Up to 8 unique colors
        labels = [f'$s_a$', f'$s_c$', f'$p_xa$', f'$p_ya$', f'$p_za$', f'$p_xc$', f'$p_yc$', f'$p_zc$', f'$s^*a$', rf'$s^*c$']  # orbital labels
        linestyles = [':','-','-','-','-','-','-','-',':','-']
        markers =['.','.','None','None','None','None','None','None','.','.']
        #for i in range(8):
        #    plt.plot(k_norm, E_eigenvectors_valencia[:, i].real,color=colors[i],linestyle='--',linewidth=1)  # Valence bands E1–E4

        for i in range(10):
            plt.plot(k_norm, np.abs(E_eigenvectors_conduccion[:, i]),color=colors[i],linestyle='-',linewidth=1, label=labels[i])  # Conduction bands E5–E8
        
        # **Set x and y limits**
        plt.xlim(-1, 1)  # Set x-axis limits
        #plt.ylim(-15, 15)  # Adjust y-axis as needed
        
        # Add vertical lines
        plt.axvline(x=k_norm[m_LG-1], color='k', linestyle='--', linewidth=1)   # L -> Γ
        plt.axvline(x=k_norm[m_LG+m_GX-1], color='k', linestyle='--', linewidth=1) # Γ -> X
        plt.axvline(x=k_norm[m_LG+m_GX+m_XK-1], color='k', linestyle='--', linewidth=1) # X -> K

        # Remove x-axis numbers
        plt.xticks([], [])  # Empty labels
        # Add custom letters at certain x positions
        plt.text(-1.0, -0.1, 'L', fontsize=12, ha='center', va='bottom')
        plt.text(k_norm[m_LG-1], -0.1, r'$\Gamma$', fontsize=12, ha='center', va='bottom')
        plt.text(k_norm[m_LG+m_GX-1], -0.1, 'X', fontsize=12, ha='center', va='bottom')
        plt.text(k_norm[m_LG+m_GX+m_XK-1], -0.1, 'U, K', fontsize=12, ha='center', va='bottom')
        plt.text(1.0, -0.1, r'$\Gamma$', fontsize=12, ha='center', va='bottom')
        # New letters in between
        plt.text((-1.0 + k_norm[m_LG-1]) / 2, -0.1, r'$\Lambda$', fontsize=12, ha='center', va='bottom')
        plt.text((k_norm[m_LG-1] + k_norm[m_LG+m_GX-1]) / 2, -0.1, r'$\Delta$', fontsize=12, ha='center', va='bottom')
        #plt.text((k_norm[m_LG+m_GX-1] + k_norm[m_LG+m_GX+m_XK-1]) / 2, -0.15, 'U', fontsize=12, ha='center', va='top')
        plt.text((k_norm[m_LG+m_GX+m_XK-1] + 1.0) / 2, -0.1, r'$\Sigma$', fontsize=12, ha='center', va='bottom')

        y_max = max(np.max(np.abs(E_eigenvectors_conduccion[:, i])) for i in range(10))
        plt.text(-0.925, y_max*0.95, cond, fontsize=12, ha='left', va='top')
        plt.legend(loc='center left', bbox_to_anchor=(1.02, 0.5), fontsize=9)
        #plt.xlabel(r'vector de onda $k$')
        #plt.gcf().set_size_inches(4.5,4.5)
        plt.tight_layout()
        plt.savefig(f'Figure_autof_conduccion_{cond}.pdf')
        plt.show()


        for i in [0, 1, 8, 9]: #loop en los s and s*
            plt.plot(k_norm, np.abs(E_eigenvectors_valencia[:, i]),color=colors[i],linestyle='-',linewidth=1, label=labels[i])  # Conduction bands E5–E8
        
        # **Set x and y limits**
        plt.xlim(-1, 1)  # Set x-axis limits
        #plt.ylim(-15, 15)  # Adjust y-axis as needed
        
        # Add vertical lines
        plt.axvline(x=k_norm[m_LG-1], color='k', linestyle='--', linewidth=1)   # L -> Γ
        plt.axvline(x=k_norm[m_LG+m_GX-1], color='k', linestyle='--', linewidth=1) # Γ -> X
        plt.axvline(x=k_norm[m_LG+m_GX+m_XK-1], color='k', linestyle='--', linewidth=1) # X -> K

        # Remove x-axis numbers
        plt.xticks([], [])  # Empty labels
        # Add custom letters at certain x positions
        y_max = max(np.max(np.abs(E_eigenvectors_valencia[:, i])) for i in [0, 1, 8, 9])
        y_min = max(np.min(np.abs(E_eigenvectors_valencia[:, i])) for i in [0, 1, 8, 9])
        plt.text(-1.0, y_min*0.1, 'L', fontsize=12, ha='center', va='bottom')
        plt.text(k_norm[m_LG-1], y_min*0.1, r'$\Gamma$', fontsize=12, ha='center', va='bottom')
        plt.text(k_norm[m_LG+m_GX-1], y_min*0.1, 'X', fontsize=12, ha='center', va='bottom')
        plt.text(k_norm[m_LG+m_GX+m_XK-1], y_min*0.1, 'U, K', fontsize=12, ha='center', va='bottom')
        plt.text(1.0, y_min*0.1, r'$\Gamma$', fontsize=12, ha='center', va='bottom')
        # New letters in between
        plt.text((-1.0 + k_norm[m_LG-1]) / 2, y_min*0.1, r'$\Lambda$', fontsize=12, ha='center', va='bottom')
        plt.text((k_norm[m_LG-1] + k_norm[m_LG+m_GX-1]) / 2, y_min*0.1, r'$\Delta$', fontsize=12, ha='center', va='bottom')
        plt.text((k_norm[m_LG+m_GX+m_XK-1] + 1.0) / 2, y_min*0.1, r'$\Sigma$', fontsize=12, ha='center', va='bottom')
        
        # compute y_max for bands 0,1,8,9 inline

        plt.text(-0.925, y_max*0.95, cond, fontsize=12, ha='left', va='top')
        plt.legend(loc='center left', bbox_to_anchor=(1.02, 0.5), fontsize=9)
        #plt.xlabel(r'vector de onda $k$')
        #plt.gcf().set_size_inches(4.5,4.5)
        plt.tight_layout()
        plt.savefig(f'Figure_autof_valencias_{cond}.pdf')
        plt.show()

        for i in range(2, 8): #loop en los p
            plt.plot(k_norm, np.abs(E_eigenvectors_valencia[:, i]),color=colors[i],linestyle='-',linewidth=1, label=labels[i])  # Conduction bands E5–E8
        
        # **Set x and y limits**
        plt.xlim(-0.5, -0.3)  # Set x-axis limits
        #plt.ylim(-15, 15)  # Adjust y-axis as needed
        
        # Add vertical lines
        plt.axvline(x=k_norm[m_LG-1], color='k', linestyle='--', linewidth=1)   # L -> Γ
        plt.axvline(x=k_norm[m_LG+m_GX-1], color='k', linestyle='--', linewidth=1) # Γ -> X
        plt.axvline(x=k_norm[m_LG+m_GX+m_XK-1], color='k', linestyle='--', linewidth=1) # X -> K

        # Remove x-axis numbers
        plt.xticks([], [])  # Empty labels
        # Add custom letters at certain x positions

        y_max = max(np.max(np.abs(E_eigenvectors_valencia[:, i])) for i in range(2, 8))
        y_min = max(np.min(np.abs(E_eigenvectors_valencia[:, i])) for i in range(2, 8))

        plt.text(-1.0, y_min*0.5, 'L', fontsize=12, ha='center', va='bottom')
        plt.text(k_norm[m_LG-1], y_min*0.5, r'$\Gamma$', fontsize=12, ha='center', va='bottom')
        plt.text(k_norm[m_LG+m_GX-1], y_min*0.5, 'X', fontsize=12, ha='center', va='bottom')
        plt.text(k_norm[m_LG+m_GX+m_XK-1], y_min*0.5, 'U, K', fontsize=12, ha='center', va='bottom')
        plt.text(1.0, y_min*0.5, r'$\Gamma$', fontsize=12, ha='center', va='bottom')
        # New letters in between
        plt.text((-1.0 + k_norm[m_LG-1]) / 2, y_min*0.5, r'$\Lambda$', fontsize=12, ha='center', va='bottom')
        plt.text((k_norm[m_LG-1] + k_norm[m_LG+m_GX-1]) / 2, y_min*0.5, r'$\Delta$', fontsize=12, ha='center', va='bottom')
        plt.text((k_norm[m_LG+m_GX+m_XK-1] + 1.0) / 2, y_min*0.5, r'$\Sigma$', fontsize=12, ha='center', va='bottom')
        
        plt.text(-0.49, y_max*0.95, cond, fontsize=12, ha='left', va='top')
        plt.legend(loc='center left', bbox_to_anchor=(1.02, 0.5), fontsize=9)
        #plt.xlabel(r'vector de onda $k$')
        #plt.gcf().set_size_inches(4.5,4.5)
        plt.tight_layout()
        plt.savefig(f'Figure_autof_valenciap_close_{cond}.pdf')
        plt.show()


        for i in range(2, 8): #loop en los p
            plt.plot(k_norm, np.abs(E_eigenvectors_valencia[:, i]),color=colors[i],linestyle='-',linewidth=1, label=labels[i])  # Conduction bands E5–E8
        
        # **Set x and y limits**
        plt.xlim(-1, 1)  # Set x-axis limits
        #plt.ylim(-15, 15)  # Adjust y-axis as needed
        
        # Add vertical lines
        plt.axvline(x=k_norm[m_LG-1], color='k', linestyle='--', linewidth=1)   # L -> Γ
        plt.axvline(x=k_norm[m_LG+m_GX-1], color='k', linestyle='--', linewidth=1) # Γ -> X
        plt.axvline(x=k_norm[m_LG+m_GX+m_XK-1], color='k', linestyle='--', linewidth=1) # X -> K

        # Remove x-axis numbers
        plt.xticks([], [])  # Empty labels
        # Add custom letters at certain x positions
        plt.text(-1.0, -0.1, 'L', fontsize=12, ha='center', va='bottom')
        plt.text(k_norm[m_LG-1], -0.1, r'$\Gamma$', fontsize=12, ha='center', va='bottom')
        plt.text(k_norm[m_LG+m_GX-1], -0.1, 'X', fontsize=12, ha='center', va='bottom')
        plt.text(k_norm[m_LG+m_GX+m_XK-1], -0.1, 'U, K', fontsize=12, ha='center', va='bottom')
        plt.text(1.0, -0.1, r'$\Gamma$', fontsize=12, ha='center', va='bottom')
        # New letters in between
        plt.text((-1.0 + k_norm[m_LG-1]) / 2, -0.1, r'$\Lambda$', fontsize=12, ha='center', va='bottom')
        plt.text((k_norm[m_LG-1] + k_norm[m_LG+m_GX-1]) / 2, -0.1, r'$\Delta$', fontsize=12, ha='center', va='bottom')
        plt.text((k_norm[m_LG+m_GX+m_XK-1] + 1.0) / 2, -0.1, r'$\Sigma$', fontsize=12, ha='center', va='bottom')
        
        y_max = max(np.max(np.abs(E_eigenvectors_valencia[:, i])) for i in range(2, 8))
        plt.text(-0.925, y_max*0.95, cond, fontsize=12, ha='left', va='top')
        plt.legend(loc='center left', bbox_to_anchor=(1.02, 0.5), fontsize=9)
        #plt.xlabel(r'vector de onda $k$')
        #plt.gcf().set_size_inches(4.5,4.5)
        plt.tight_layout()
        plt.savefig(f'Figure_autof_valenciap_{cond}.pdf')
        plt.show()



        # De los dos gráficas anteriores vemos que la capa de valencia no tiene prácticamente contribuciones de los
        # orbitales s y s*, y son combinaciones lineales de los p


        #lo mismo, pero en un rango más pequeño de k
        # plot 2) Eigenvectors para las bandas de conduccion y valencia en GaAs
        colors = ['r', 'b', 'g', 'm', 'c', 'y', 'k', 'orange', 'darkred', 'navy']  # Up to 8 unique colors
        labels = [f'$s_a$', f'$s_c$', f'$p_xa$', f'$p_ya$', f'$p_za$', f'$p_xc$', f'$p_yc$', f'$p_zc$', f'$s^*a$', rf'$s^*c$']  # orbital labels

        #for i in range(8):
        #    plt.plot(k_norm, E_eigenvectors_valencia[:, i].real,color=colors[i],linestyle='--',linewidth=1)  # Valence bands E1–E4

        for i in range(10):
            plt.plot(k_norm, np.abs(E_eigenvectors_conduccion[:, i]),color=colors[i],linestyle='-',linewidth=1, label=labels[i])  # Conduction bands E5–E8
        
        # **Set x and y limits**
        plt.xlim(-1, 0.20212121212121)  # Set x-axis limits, approx
        #plt.ylim(-15, 15)  # Adjust y-axis as needed
        
        # Add vertical lines
        plt.axvline(x=k_norm[m_LG-1], color='k', linestyle='--', linewidth=1)   # L -> Γ
        plt.axvline(x=k_norm[m_LG+m_GX-1], color='k', linestyle='--', linewidth=1) # Γ -> X
        plt.axvline(x=k_norm[m_LG+m_GX+m_XK-1], color='k', linestyle='--', linewidth=1) # X -> K

        # Remove x-axis numbers
        plt.xticks([], [])  # Empty labels
        # Add custom letters at certain x positions
        plt.text(-1.0, -0.1, 'L', fontsize=12, ha='center', va='bottom')
        plt.text(k_norm[m_LG-1], -0.1, r'$\Gamma$', fontsize=12, ha='center', va='bottom')
        plt.text(k_norm[m_LG+m_GX-1], -0.1, 'X', fontsize=12, ha='center', va='bottom')
        #plt.text(k_norm[m_LG+m_GX+m_XK-1], -0.1, 'U, K', fontsize=12, ha='center', va='bottom')
        #plt.text(1.0, -0.1, r'$\Gamma$', fontsize=12, ha='center', va='bottom')
        # New letters in between
        #plt.text((-1.0 + k_norm[m_LG-1]) / 2, -0.1, r'$\Lambda$', fontsize=12, ha='center', va='bottom')
        #plt.text((k_norm[m_LG-1] + k_norm[m_LG+m_GX-1]) / 2, -0.1, r'$\Delta$', fontsize=12, ha='center', va='bottom')
        #plt.text((k_norm[m_LG+m_GX-1] + k_norm[m_LG+m_GX+m_XK-1]) / 2, -0.15, 'U', fontsize=12, ha='center', va='top')
        #plt.text((k_norm[m_LG+m_GX+m_XK-1] + 1.0) / 2, -0.1, r'$\Sigma$', fontsize=12, ha='center', va='bottom')

        y_max = max(np.max(np.abs(E_eigenvectors_conduccion[:, i])) for i in range(10))
        plt.text(-0.925, y_max*0.95, cond, fontsize=12, ha='left', va='top')
        plt.legend(loc='center left', bbox_to_anchor=(1, 0.5),fontsize=9)
        #plt.xlabel(r'vector de onda $k$')
        #plt.gcf().set_size_inches(4.5,4.5)
        plt.tight_layout()
        plt.savefig(f'Figure_autof_conduccion_close_{cond}.pdf')
        plt.show()


        ### SO ###

        if cond == 'GaAs':

            #me falta poner autofunciones en SO

            capa_valencia_SO=E_SO_eigenvalues[:,7].real #la de menor energia en el gap
            capa_conducion_SO=E_SO_eigenvalues[:,9].real #la de mayor energia en el gap    

            min_conducion_index_SO = np.argmin(capa_conducion_SO)
            min_conducion_energy_SO = capa_conducion_SO[min_conducion_index_SO]

            max_valencia_index_SO = np.argmax(capa_valencia_SO)
            max_valencia_energy_SO = capa_valencia_SO[max_valencia_index_SO]

            energy_gap_SO=abs(max_valencia_energy_SO-min_conducion_energy_SO)

            if min_conducion_index_SO == max_valencia_index_SO:
                print(f"{cond} is a direct semiconductor in SO with a energy gap of {energy_gap_SO:.3f} eV")
            else:
                print(f"{cond} is a indirect semiconductor in SO with a energy gap of {energy_gap_SO:.3f} eV") 


            ### Aqui la figura de la parte de SO
            plt.plot(k_norm, E_SO_eigenvalues[:,0].real,color='k', linewidth=1, label='E1')  # Take real part
            plt.plot(k_norm, E_SO_eigenvalues[:,1].real,color='k', linewidth=1, label='E2')  # Take real part
            plt.plot(k_norm, E_SO_eigenvalues[:,2].real,color='k', linewidth=1, label='E3')  # Take real part
            plt.plot(k_norm, E_SO_eigenvalues[:,3].real,color='k', linewidth=1, label='E4')  # Take real part
            plt.plot(k_norm, E_SO_eigenvalues[:,4].real,color='k', linewidth=1, label='E5')  # Take real part
            plt.plot(k_norm, E_SO_eigenvalues[:,5].real,color='k', linewidth=1, label='E6')  # Take real part
            plt.plot(k_norm, E_SO_eigenvalues[:,6].real,color='k', linewidth=1, label='E7')  # Take real part
            plt.plot(k_norm, E_SO_eigenvalues[:,7].real,color='k', linewidth=1, label='E8')  # Take real part
            plt.plot(k_norm, E_SO_eigenvalues[:,8].real,color='k', linewidth=1, label='E9')  # Take real part
            plt.plot(k_norm, E_SO_eigenvalues[:,9].real,color='k', linewidth=1, label='E10')  # Take real part
            plt.plot(k_norm, E_SO_eigenvalues[:,10].real,color='k', linewidth=1, label='E11')  # Take real part
            plt.plot(k_norm, E_SO_eigenvalues[:,11].real,color='k', linewidth=1, label='E12')  # Take real part
            plt.plot(k_norm, E_SO_eigenvalues[:,12].real,color='k', linewidth=1, label='E13')  # Take real part
            plt.plot(k_norm, E_SO_eigenvalues[:,13].real,color='k', linewidth=1, label='E14')  # Take real part
            plt.plot(k_norm, E_SO_eigenvalues[:,14].real,color='k', linewidth=1, label='E15')  # Take real part
            plt.plot(k_norm, E_SO_eigenvalues[:,15].real,color='k', linewidth=1, label='E16')  # Take real part
            plt.plot(k_norm, E_SO_eigenvalues[:,16].real,color='k', linewidth=1, label='E17')  # Take real part
            plt.plot(k_norm, E_SO_eigenvalues[:,17].real,color='k', linewidth=1, label='E18')  # Take real part
            plt.plot(k_norm, E_SO_eigenvalues[:,18].real,color='k', linewidth=1, label='E19')  # Take real part
            plt.plot(k_norm, E_SO_eigenvalues[:,19].real,color='k', linewidth=1, label='E20')  # Take real part

            # **Set x and y limits**
            plt.xlim(-1, 1)  # Set x-axis limits
            plt.ylim(-15, 15)  # Adjust y-axis as needed

            # Add vertical lines
            plt.axvline(x=k_norm[m_LG-1], color='k', linestyle='--', linewidth=1)   # L -> Γ
            plt.axvline(x=k_norm[m_LG+m_GX-1], color='k', linestyle='--', linewidth=1) # Γ -> X
            plt.axvline(x=k_norm[m_LG+m_GX+m_XK-1], color='k', linestyle='--', linewidth=1) # X -> K

            # Add horizontal lines at the energy values where the min and max energies occur
            plt.axhline(y=min_conducion_energy, color='k', linestyle='--', linewidth=1)
            plt.axhline(y=max_valencia_energy, color='k', linestyle='--', linewidth=1)

            # Shaded area (between the horizontal lines)
            plt.fill_between(k_norm, min_conducion_energy, max_valencia_energy, where=(capa_conducion > min_conducion_energy) & (capa_valencia < max_valencia_energy), 
                         color='gray', alpha=0.3, label='Energy Gap')


            # Remove x-axis numbers
            plt.xticks([], [])  # Empty labels
            # Add custom letters at certain x positions
            plt.text(-1.0, -16.5, 'L', fontsize=12, ha='center', va='bottom')
            plt.text(k_norm[m_LG-1], -16.5, r'$\Gamma$', fontsize=12, ha='center', va='bottom')
            plt.text(k_norm[m_LG+m_GX-1], -16.5, 'X', fontsize=12, ha='center', va='bottom')
            plt.text(k_norm[m_LG+m_GX+m_XK-1], -16.5, 'U, K', fontsize=12, ha='center', va='bottom')
            plt.text(1.0, -16.5, r'$\Gamma$', fontsize=12, ha='center', va='bottom')
            # New letters in between
            plt.text((-1.0 + k_norm[m_LG-1]) / 2, -16.5, r'$\Lambda$', fontsize=12, ha='center', va='bottom')
            plt.text((k_norm[m_LG-1] + k_norm[m_LG+m_GX-1]) / 2, -16.5, r'$\Delta$', fontsize=12, ha='center', va='bottom')
            #plt.text((k_norm[m_LG+m_GX-1] + k_norm[m_LG+m_GX+m_XK-1]) / 2, -16.5, 'U', fontsize=12, ha='center', va='top')
            plt.text((k_norm[m_LG+m_GX+m_XK-1] + 1.0) / 2, -16.5, r'$\Sigma$', fontsize=12, ha='center', va='bottom')

            plt.text(-0.925, 15 -0.25, cond, fontsize=12, ha='left', va='top')

            #plt.xlabel(r'vector de onda $k$')
            plt.ylabel('E (eV)')
            plt.gcf().set_size_inches(4.5,4.5)
            plt.tight_layout()
            plt.savefig(f'Figure_SO_{cond}.pdf')
            plt.show()

            ### segunda grafica para SO, donde ahora el tomo el rango de K más pequeño

            plt.plot(k_norm, E_SO_eigenvalues[:,0].real,color='k', linewidth=1, label='E1')  # Take real part
            plt.plot(k_norm, E_SO_eigenvalues[:,1].real,color='k', linewidth=1, label='E2')  # Take real part
            plt.plot(k_norm, E_SO_eigenvalues[:,2].real,color='k', linewidth=1, label='E3')  # Take real part
            plt.plot(k_norm, E_SO_eigenvalues[:,3].real,color='k', linewidth=1, label='E4')  # Take real part
            plt.plot(k_norm, E_SO_eigenvalues[:,4].real,color='k', linewidth=1, label='E5')  # Take real part
            plt.plot(k_norm, E_SO_eigenvalues[:,5].real,color='k', linewidth=1, label='E6')  # Take real part
            plt.plot(k_norm, E_SO_eigenvalues[:,6].real,color='k', linewidth=1, label='E7')  # Take real part
            plt.plot(k_norm, E_SO_eigenvalues[:,7].real,color='k', linewidth=1, label='E8')  # Take real part
            plt.plot(k_norm, E_SO_eigenvalues[:,8].real,color='k', linewidth=1, label='E9')  # Take real part
            plt.plot(k_norm, E_SO_eigenvalues[:,9].real,color='k', linewidth=1, label='E10')  # Take real part
            plt.plot(k_norm, E_SO_eigenvalues[:,10].real,color='k', linewidth=1, label='E11')  # Take real part
            plt.plot(k_norm, E_SO_eigenvalues[:,11].real,color='k', linewidth=1, label='E12')  # Take real part
            plt.plot(k_norm, E_SO_eigenvalues[:,12].real,color='k', linewidth=1, label='E13')  # Take real part
            plt.plot(k_norm, E_SO_eigenvalues[:,13].real,color='k', linewidth=1, label='E14')  # Take real part
            plt.plot(k_norm, E_SO_eigenvalues[:,14].real,color='k', linewidth=1, label='E15')  # Take real part
            plt.plot(k_norm, E_SO_eigenvalues[:,15].real,color='k', linewidth=1, label='E16')  # Take real part
            plt.plot(k_norm, E_SO_eigenvalues[:,16].real,color='k', linewidth=1, label='E17')  # Take real part
            plt.plot(k_norm, E_SO_eigenvalues[:,17].real,color='k', linewidth=1, label='E18')  # Take real part
            plt.plot(k_norm, E_SO_eigenvalues[:,18].real,color='k', linewidth=1, label='E19')  # Take real part
            plt.plot(k_norm, E_SO_eigenvalues[:,19].real,color='k', linewidth=1, label='E20')  # Take real part

            plt.legend(fontsize=6, loc='upper right', ncol=2)  # Adjust legend size and placement

            # **Set x and y limits**
            plt.xlim(-1, 0.1212121212121)  # Set x-axis limits
            plt.ylim(-15, 15)  # Adjust y-axis as needed

            # Add vertical lines
            plt.axvline(x=k_norm[m_LG-1], color='k', linestyle='--', linewidth=1)   # L -> Γ
            plt.axvline(x=k_norm[m_LG+m_GX-1], color='k', linestyle='--', linewidth=1) # Γ -> X
            plt.axvline(x=k_norm[m_LG+m_GX+m_XK-1], color='k', linestyle='--', linewidth=1) # X -> K
    
            # Add horizontal lines at the energy values where the min and max energies occur
            plt.axhline(y=min_conducion_energy, color='k', linestyle='--', linewidth=1)
            plt.axhline(y=max_valencia_energy, color='k', linestyle='--', linewidth=1)

            # Shaded area (between the horizontal lines)
            plt.fill_between(k_norm, min_conducion_energy, max_valencia_energy, where=(capa_conducion > min_conducion_energy) & (capa_valencia < max_valencia_energy), 
                         color='gray', alpha=0.3, label='Energy Gap')


            # Remove x-axis numbers
            plt.xticks([], [])  # Empty labels
            # Add custom letters at certain x positions
            plt.text(-1.0, -16.5, 'L', fontsize=12, ha='center', va='bottom')
            plt.text(k_norm[m_LG-1], -16.5, r'$\Gamma$', fontsize=12, ha='center', va='bottom')
            plt.text(k_norm[m_LG+m_GX-1], -16.5, 'X', fontsize=12, ha='center', va='bottom')
            #plt.text(k_norm[m_LG+m_GX+m_XK-1], -16.5, 'U, K', fontsize=12, ha='center', va='bottom')
            #plt.text(1.0, -16.5, r'$\Gamma$', fontsize=12, ha='center', va='bottom')
            # New letters in between
            #plt.text((-1.0 + k_norm[m_LG-1]) / 2, -16.5, r'$\Lambda$', fontsize=12, ha='center', va='bottom')
            #plt.text((k_norm[m_LG-1] + k_norm[m_LG+m_GX-1]) / 2, -16.5, r'$\Delta$', fontsize=12, ha='center', va='bottom')
            #plt.text((k_norm[m_LG+m_GX-1] + k_norm[m_LG+m_GX+m_XK-1]) / 2, -16.5, 'U', fontsize=12, ha='center', va='top')
            #plt.text((k_norm[m_LG+m_GX+m_XK-1] + 1.0) / 2, -16.5, r'$\Sigma$', fontsize=12, ha='center', va='bottom')

            plt.text(-0.925, 15 -0.25, cond, fontsize=12, ha='left', va='top')

            #plt.xlabel(r'vector de onda $k$')
            plt.ylabel('E (eV)')
            plt.gcf().set_size_inches(4.5,4.5)
            plt.tight_layout()
            plt.savefig(f'Figure_SO_{cond}.pdf')
            plt.show()


