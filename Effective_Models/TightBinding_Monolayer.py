import numpy as np 
import scipy 
import cmath 
import scipy.linalg as sla 
import matplotlib.pyplot as plt 

# The geometrical parameters 
a = 0.2 # [nm] lattice parameter 
d = a/2.2 # Size of basis 
alpha = -0.1   # C4 symmetry breaking parameter 
D = 0.25    # interlayer distance 

# Intralayer hopping 
t1 = -np.exp(-d**2/a**2)
t2 = -np.exp(-(a-d)**2/a**2)
t3 = -0.5*np.exp(-2*d**2/a**2)
print('t1 = '+str(t1))
print('t2 = '+str(t2))
print('t3 = '+str(t3))

# The list of momenta 
Gamma_x = 0 
Gamma_y = 0 
X_x = np.pi/a 
X_y = 0
M_x = np.pi/a 
M_y = np.pi/a 

Nk = 101 

kx_array = np.concatenate( (np.linspace(Gamma_x,X_x,Nk), np.linspace(X_x,M_x,Nk)[1:Nk], np.linspace(M_x,Gamma_x,Nk)[1:Nk]) ) 
ky_array = np.concatenate( (np.linspace(Gamma_y,X_y,Nk), np.linspace(X_y,M_y,Nk)[1:Nk], np.linspace(M_y,Gamma_y,Nk)[1:Nk]) ) 

#print(kx_array)
#print(ky_array)

# Arrays of energy 
E_array = np.zeros((3*Nk-2,4)) 

# We scan the momentum arrays 
for i in range(3*Nk-2):
    # The momenta 
    kx = kx_array[i]
    ky = ky_array[i]

    # The Hamiltonian 
    Hamiltonian = np.zeros((4,4),dtype=complex) 

    Hamiltonian[0,0] = 0 
    Hamiltonian[0,1] = t1*cmath.exp(1j*kx*d) + t2*cmath.exp(1j*kx*(d-a))
    Hamiltonian[0,2] = t3*(1-alpha)*cmath.exp(1j*(kx*d+ky*d)) 
    Hamiltonian[0,3] = t1*cmath.exp(1j*ky*d) + t2*cmath.exp(1j*ky*(d-a)) 

    Hamiltonian[1,0] = np.conjugate(Hamiltonian[0,1]) 
    Hamiltonian[1,1] = 0 
    Hamiltonian[1,2] = t1*cmath.exp(1j*ky*d) + t2*cmath.exp(1j*ky*(d-a)) 
    Hamiltonian[1,3] = t3*(1+alpha)*cmath.exp(1j*(-kx*d+ky*d))

    Hamiltonian[2,0] = np.conjugate(Hamiltonian[0,2]) 
    Hamiltonian[2,1] = np.conjugate(Hamiltonian[1,2])
    Hamiltonian[2,2] = 0 
    Hamiltonian[2,3] = t1*cmath.exp(-1j*kx*d) + t2*cmath.exp(-1j*kx*(d-a)) 

    Hamiltonian[3,0] = np.conjugate(Hamiltonian[0,3]) 
    Hamiltonian[3,1] = np.conjugate(Hamiltonian[1,3])
    Hamiltonian[3,2] = np.conjugate(Hamiltonian[2,3]) 
    Hamiltonian[3,3] = 0 

    Evals,Evecs = sla.eigh(Hamiltonian)

    E_array[i,:] = Evals 

# Plot the band structure 
fig,ax = plt.subplots()
ax.plot(E_array[:,0:4])
ax.vlines(0,-8,8,color='black')
ax.vlines(Nk-1,-8,8,color='black')
ax.vlines(2*(Nk-1),-8,8,color='black')
ax.vlines(3*(Nk-1),-8,8,color='black')
tick_locs = [i*(Nk-1) for i in range(4)]
tick_labs = [r'$\Gamma$','X','M',r'$\Gamma$']
ax.set_xticks(tick_locs)
ax.set_xticklabels(tick_labs,size=16)
ax.set_xlim(0,3*(Nk-1))
ax.set_ylim(-6,6)
ax.set_ylabel('Energy (eV)')
ax.set_title(r'$\alpha = $'+str(alpha))
plt.show()
