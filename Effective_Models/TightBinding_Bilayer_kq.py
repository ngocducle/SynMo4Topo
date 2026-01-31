import numpy as np 
import scipy 
import cmath 
import scipy.linalg as sla 
import matplotlib.pyplot as plt

# The geometrical parameters 
a = 0.2 # lattice parameter (nm)
d = a/2.2 # Size of basis 
alpha = -0.1 # C4 symmetry breaking parameter 
D = 0.25 # interlayer distance (nm) 
xAtoms = [0,d,d,0] # arrays of x-coordinates of atoms 
yAtoms = [0,0,d,d] # arrays of y-coordiantes of atoms 

# Tight-binding Minh model for monolayer 
t1 = -np.exp(-d**2/a**2)       # intracell hopping integrals along edge 
t2 = -np.exp(-(a-d)**2/a**2)   # intercell hopping integrals 
t3 = -0.5*np.exp(-2*d**2/a**2) # intracell hopping integrals along diagonals 

print('t1 = '+str(t1))
print('t2 = '+str(t2)) 
print('t3 = '+str(t3))

delta_x = 0 # shift (nm)
delta_y = 0 # shift (nm) 

# The list of momenta 
Nk = 21 
k_array = np.linspace(0,0.2,Nk)

Nq = 51 
q_array = np.linspace(0,0.5,Nq)

##### The bilayer Hamiltonian 
def TightBinding_Bilayer_Hamiltonian(t1,t2,t3,alpha,a,d,D,kx,ky,delta_x,delta_y):
    Hamiltonian = np.zeros((8,8),dtype=complex)

    xAtoms = [0,d,d,0]
    yAtoms = [0,0,d,d] 

    # Block (0,0) 
    Hamiltonian[0,0] = 0.0 
    Hamiltonian[0,1] = t1*cmath.exp(1j*kx*d) + t2*cmath.exp(1j*kx*(d-a))
    Hamiltonian[0,2] = t3*(1-alpha)*cmath.exp(1j*(kx*d+ky*d))
    Hamiltonian[0,3] = t1*cmath.exp(1j*ky*d) + t2*cmath.exp(1j*ky*(d-a))

    Hamiltonian[1,0] = np.conjugate(Hamiltonian[0,1])
    Hamiltonian[1,1] = 0.0
    Hamiltonian[1,2] = t1*cmath.exp(1j*ky*d) + t2*cmath.exp(1j*ky*(d-a))
    Hamiltonian[1,3] = t3*(1+alpha)*cmath.exp(1j*(-kx*d+ky*d))

    Hamiltonian[2,0] = np.conjugate(Hamiltonian[0,2])
    Hamiltonian[2,1] = np.conjugate(Hamiltonian[1,2])
    Hamiltonian[2,2] = 0.0 
    Hamiltonian[2,3] = t1*cmath.exp(-1j*kx*d) + t2*cmath.exp(-1j*kx*(d-a))

    Hamiltonian[3,0] = np.conjugate(Hamiltonian[0,3])
    Hamiltonian[3,1] = np.conjugate(Hamiltonian[1,3])
    Hamiltonian[3,2] = np.conjugate(Hamiltonian[2,3])
    Hamiltonian[3,3] = 0.0 

    # Block (0,1) 
    for j in range(0,4):
        for k in range(0,4):
            Hamiltonian[j,4+k] = 0 

            # The sum over R 
            for nx in range(-3,4):
                for ny in range(-3,4):
                    #disp_x = nx*a + a/2 + xAtoms[j] + delta_x - xAtoms[k] 
                    #disp_y = ny*a + a/2 + yAtoms[j] + delta_y - yAtoms[k] 

                    disp_x = xAtoms[k] - xAtoms[j] - a/2 - delta_x + nx*a 
                    disp_y = yAtoms[k] - yAtoms[j] - a/2 - delta_y + ny*a 

                    c2 = D**2 / (disp_x**2 + disp_y**2 + D**2) 

                    if np.sqrt(disp_x**2 + disp_y**2) < 4*a:
                        t = np.exp(-(disp_x**2 + disp_y**2 + D**2)/a**2) * (2.5*c2-2) 
                    else:
                        t = 0

                    Hamiltonian[j,4+k] = Hamiltonian[j,4+k] + t*cmath.exp(1j*(kx*disp_x+ky*disp_y))

    # Block (1,0)
    for j in range(0,4):
        for k in range(0,4):
            Hamiltonian[4+j,k] = np.conjugate(Hamiltonian[k,4+j]) 

    # Block (1,1)
    Hamiltonian[4,4] = 0.0
    Hamiltonian[4,5] = t1*cmath.exp(1j*kx*d) + t2*cmath.exp(1j*kx*(d-a))
    Hamiltonian[4,6] = t3*(1-alpha)*cmath.exp(1j*(kx*d+ky*d))
    Hamiltonian[4,7] = t1*cmath.exp(1j*ky*d) + t2*cmath.exp(1j*ky*(d-a))

    Hamiltonian[5,4] = np.conjugate(Hamiltonian[4,5])
    Hamiltonian[5,5] = 0.0 
    Hamiltonian[5,6] = t1*cmath.exp(1j*ky*d) + t2*cmath.exp(1j*ky*(d-a))
    Hamiltonian[5,7] = t3*(1+alpha)*cmath.exp(1j*(-kx*d+ky*d))

    Hamiltonian[6,4] = np.conjugate(Hamiltonian[4,6])
    Hamiltonian[6,5] = np.conjugate(Hamiltonian[5,6])
    Hamiltonian[6,6] = 0.0 
    Hamiltonian[6,7] = t1*cmath.exp(-1j*kx*d) + t2*cmath.exp(-1j*kx*(d-a))

    Hamiltonian[7,4] = np.conjugate(Hamiltonian[4,7])
    Hamiltonian[7,5] = np.conjugate(Hamiltonian[5,7]) 
    Hamiltonian[7,6] = np.conjugate(Hamiltonian[6,7])
    Hamiltonian[7,7] = 0.0 


    return Hamiltonian 

### ======================================================================########################
# The list of momenta 
Gamma_x = 0 
Gamma_y = 0 
X_x = np.pi/a 
X_y = 0 
M_x = np.pi/a 
M_y = np.pi/a 

Nk = 1001 

kx_array = np.concatenate( (np.linspace(Gamma_x,X_x,Nk), np.linspace(X_x,M_x,Nk)[1:Nk],np.linspace(M_x,Gamma_x,Nk)[1:Nk]) )
ky_array = np.concatenate( (np.linspace(Gamma_y,X_y,Nk), np.linspace(X_y,M_y,Nk)[1:Nk],np.linspace(M_y,Gamma_y,Nk)[1:Nk]) )

# Arrays of energy 
E_array = np.zeros((3*Nk-2,8)) 

# We scan the momentum arrays 
for i in range(0,3*Nk-2):
    # The momenta 
    kx = kx_array[i]
    ky = ky_array[i]

    # The Hamiltonian 
    Hamiltonian = TightBinding_Bilayer_Hamiltonian(t1,t2,t3,alpha,a,d,D,kx,ky,delta_x,delta_y)

    # Eigenvalues and eigenvectors 
    Evals,Evecs = sla.eigh(Hamiltonian)

    E_array[i,:] = Evals 

# Plot the band structure 
fig,ax = plt.subplots()
ax.plot(E_array[:,0:8],color='purple') 
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
plt.show()

### ==============================================================================================
### Plot the 2D band structure 
Nk = 41 
k_array = np.linspace(-0.2,0.2,Nk)

Nq = 51
q_array = np.linspace(-0.5,0.5,Nq) 

E_array = np.zeros((Nk,Nq,8))

for ik in range(0,Nk):
    for iq in range(0,Nq):
        k = k_array[ik] 
        q = q_array[iq] 

        delta_x = q*a
        delta_y = q*a

        # The Hamiltonian 
        Hamiltonian = TightBinding_Bilayer_Hamiltonian(t1,t2,t3,alpha,a,d,D,kx,ky,delta_x,delta_y)

        # Eigenvalues and eigenvectors 
        Evals,Evecs = sla.eigh(Hamiltonian)
        E_array[ik,iq,:] = Evals 

# Plot the band structure 
X,Y = np.meshgrid(k_array,q_array)

fig,ax = plt.subplots(subplot_kw={'projection':'3d'},figsize=(12,10))
ax.plot_surface(X,Y,E_array[:,:,0].T)
ax.plot_surface(X,Y,E_array[:,:,1].T)
plt.show()
