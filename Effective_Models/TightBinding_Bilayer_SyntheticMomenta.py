import numpy as np
import scipy
import cmath
import scipy.linalg as sla
import matplotlib.pyplot as plt

# Tight-binding Minh model for monolayer
t1 = -0.81 # intracell hopping integrals along edge
t2 = -0.74 # intercell hopping integrals
t3 = -0.33 # intracell hopping integrals along diagonals
alpha = -0.1  # Symmetry-breaking factor
a = 0.2   # lattice parameter (nm)
D = 0.25   # interlayer distance (nm) 
delta_x = 0 # shift (nm) 
delta_y = 0 # shift (nm)
dxy = np.sqrt(delta_x*delta_x + delta_y*delta_y) # length of the shift (nm)

# The list of momenta
Gamma_x = 0
Gamma_y = 0
X_x = np.pi/a
X_y = 0
M_x = np.pi/a
M_y = np.pi/a

Nk = 101

kx_array = np.concatenate( (np.linspace(Gamma_x,X_x,Nk),np.linspace(X_x,M_x,Nk)[1:Nk],np.linspace(M_x,Gamma_x,Nk)[1:Nk]) )
ky_array = np.concatenate( (np.linspace(Gamma_y,X_y,Nk),np.linspace(X_y,M_y,Nk)[1:Nk],np.linspace(M_y,Gamma_y,Nk)[1:Nk]) )

print(kx_array)
print(ky_array)

# Arrays of energy
E_array = np.zeros((3*Nk-2,8))

# We scan the momentum arrays
for i in range(3*Nk-2):
    # The momenta
    kx = kx_array[i]
    ky = ky_array[i]

    # The Hamiltonian
    Hamiltonian = np.zeros((8,8),dtype=complex)

    t = np.exp(-(dxy*dxy+D*D)/(a*a)) * (2.5*D*D/(dxy*dxy+D*D) - 2)  

    Hamiltonian[0:8,0:8] = t 

    # Block (0,0) 
    Hamiltonian[0,0] = 0.0 
    Hamiltonian[0,1] = t1+t2*cmath.exp(-1j*kx*a)
    Hamiltonian[0,2] = t3*(1-alpha)
    Hamiltonian[0,3] = t1+t2*cmath.exp(-1j*ky*a)

    Hamiltonian[1,0] = np.conjugate(Hamiltonian[0,1])
    Hamiltonian[1,1] = 0.0  
    Hamiltonian[1,2] = t1+t2*cmath.exp(-1j*ky*a)
    Hamiltonian[1,3] = t3*(1+alpha)

    Hamiltonian[2,0] = np.conjugate(Hamiltonian[0,2])
    Hamiltonian[2,1] = np.conjugate(Hamiltonian[1,2])
    Hamiltonian[2,2] = 0.0 
    Hamiltonian[2,3] = t1+t2*cmath.exp(1j*kx*a)

    Hamiltonian[3,0] = np.conjugate(Hamiltonian[0,3])
    Hamiltonian[3,1] = np.conjugate(Hamiltonian[1,3])
    Hamiltonian[3,2] = np.conjugate(Hamiltonian[2,3])
    Hamiltonian[3,3] = 0.0 

    # Block (0,1) 
    Hamiltonian[0,4] = t
    Hamiltonian[1,5] = t 
    Hamiltonian[2,6] = t 
    Hamiltonian[3,7] = t 

    # Block (1,0) 
    Hamiltonian[4,0] = np.conjugate(Hamiltonian[0,4]) 
    Hamiltonian[5,1] = np.conjugate(Hamiltonian[1,5]) 
    Hamiltonian[6,2] = np.conjugate(Hamiltonian[2,6]) 
    Hamiltonian[7,3] = np.conjugate(Hamiltonian[3,7])

    # Block (1,1) 
    Hamiltonian[4,4] = 0.0 
    Hamiltonian[4,5] = t1+t2*cmath.exp(-1j*kx*a)  
    Hamiltonian[4,6] = t3*(1-alpha)
    Hamiltonian[4,7] = t1+t2*cmath.exp(-1j*ky*a)

    Hamiltonian[5,4] = np.conjugate(Hamiltonian[4,5]) 
    Hamiltonian[5,5] = 0.0 
    Hamiltonian[5,6] = t1+t2*cmath.exp(-1j*ky*a)
    Hamiltonian[5,7] = t3*(1+alpha)

    Hamiltonian[6,4] = np.conjugate(Hamiltonian[4,6])
    Hamiltonian[6,5] = np.conjugate(Hamiltonian[5,6])
    Hamiltonian[6,6] = 0.0 
    Hamiltonian[6,7] = t1+t2*cmath.exp(1j*kx*a)

    Hamiltonian[7,4] = np.conjugate(Hamiltonian[4,7]) 
    Hamiltonian[7,5] = np.conjugate(Hamiltonian[5,7])
    Hamiltonian[7,6] = np.conjugate(Hamiltonian[6,7])
    Hamiltonian[7,7] = 0.0

    # Eigenvalues and eigenvectors
    Evals,Evecs = sla.eigh(Hamiltonian)

    E_array[i,:] = Evals

    print(Hamiltonian)
    print('t = '+str(t))

# Plot the band structure
fig,ax = plt.subplots()
ax.plot(E_array[:,0:8])
ax.vlines(0,-8,8,color='black')
ax.vlines(Nk-1,-8,8,color='black')
ax.vlines(2*(Nk-1),-8,8,color='black')
ax.vlines(3*(Nk-1),-8,8,color='black')
tick_locs = [i*(Nk-1) for i in range(4)]
tick_labs = [r'$\Gamma$','X','M',r'$\Gamma$']
ax.set_xticks(tick_locs)
ax.set_xticklabels(tick_labs,size=16)
ax.set_xlim(0,3*(Nk-1))
ax.set_ylim(-4,4)
ax.set_ylabel('Energy (eV)')
ax.set_title(r'$\alpha = $'+str(alpha))
plt.savefig('alpha_'+str(alpha)+'.png')
plt.show()
