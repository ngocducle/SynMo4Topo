import numpy as np
import scipy
import cmath
import matplotlib.pyplot

# Tight-binding Minh model for monolayer
t1 = -0.81 # intracell hopping integrals along edge
t2 = -0.74 # intercell hopping integrals
t3 = -0.33 # intracell hopping integrals along diagonals
a = 0.2    # lattice parameter (nm)

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

# We scan the momentum arrays


# The Hamiltonian
Hamiltonian = np.zeros((4,4),dtype=complex)

#Hamiltonian[0,1] = t1+2*t2*cmath.exp(-1j*)

print(Hamiltonian)
