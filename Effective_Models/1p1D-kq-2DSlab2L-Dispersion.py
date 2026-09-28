import numpy as np 
import scipy 
import scipy.linalg as sla 
import cmath 
import matplotlib.pyplot as plt 
from matplotlib import cm,colors 

##### ===========================================================================
##### FUNCTION: Hamiltonian of 2D photonic crystal slab bilayer with 
##### kx = ky = k, qx = qy = q 
def Hamiltonian(k,q,omega1,eta1,v1,U1,W1,alpha1,
                 omega2,eta2,v2,U2,W2,alpha2,
                 V,beta,dist,d0):
    H = np.zeros((8,8),dtype=complex)

    U1p = U1*(1+alpha1)
    U1m = U1*(1-alpha1)
    U2p = U2*(1+alpha2)
    U2m = U2*(1-alpha2)

    K = 2.0*np.pi 

    ### Block (1,1)
    H[0,0] = omega1 + eta1 + v1*k + v1*k*k/np.sqrt(2)
    H[0,1] = W1 
    H[0,2] = W1 
    H[0,3] = U1*(1+alpha1)

    H[1,0] = W1
    H[1,1] = omega1 - eta1 + v1*k*k/np.sqrt(2)
    H[1,2] = U1*(1-alpha1)
    H[1,3] = W1 

    H[2,0] = W1 
    H[2,1] = U1*(1-alpha1)
    H[2,2] = omega1 - eta1 + v1*k*k/np.sqrt(2)
    H[2,3] = W1 

    H[3,0] = U1*(1+alpha1)
    H[3,1] = W1 
    H[3,2] = W1 
    H[3,3] = omega1 + eta1 - v1*k + v1*k*k/np.sqrt(2)

    ### Block (1,2)
    H[0,4] = -(V + beta*k + beta*k*k/np.sqrt(2))*cmath.exp(-1j*K*q)*np.exp(-dist/d0)
    H[1,5] = (V + beta*k*k/np.sqrt(2))*np.exp(-dist/d0) 
    H[2,6] = (V + beta*k*k/np.sqrt(2))*np.exp(-dist/d0)
    H[3,7] = -(V - beta*k + beta*k*k/np.sqrt(2))*cmath.exp(1j*K*q)*np.exp(-dist/d0)

    ### Block (2,1)
    H[4,0] = -(V + beta*k + beta*k*k/np.sqrt(2))*cmath.exp(1j*K*q)*np.exp(-dist/d0)
    H[5,1] = (V + beta*k*k/np.sqrt(2))*np.exp(-dist/d0)
    H[6,2] = (V + beta*k*k/np.sqrt(2))*np.exp(-dist/d0)
    H[7,3] = -(V - beta*k + beta*k*k/np.sqrt(2))*cmath.exp(-1j*K*q)*np.exp(-dist/d0)

    ### Block (2,2)
    H[4,4] = omega2 + eta2 + v2*k + v2*k*k/np.sqrt(2)
    H[4,5] = W2 
    H[4,6] = W2 
    H[4,7] = U2*(1+alpha2)

    H[5,4] = W2 
    H[5,5] = omega2 - eta2 + v2*k*k/np.sqrt(2)
    H[5,6] = U2*(1-alpha2) 
    H[5,7] = W2 

    H[6,4] = W2 
    H[6,5] = U2*(1-alpha2)
    H[6,6] = omega2 - eta2 + v2*k*k/np.sqrt(2)
    H[6,7] = W2 

    H[7,4] = U2*(1+alpha2)
    H[7,5] = W2 
    H[7,6] = W2 
    H[7,7] = omega2 + eta2 - v2*k + v2*k*k/np.sqrt(2)

    return H

##### ============================================================================
##### The parameters of the calculations 
#omega = 0.29780940 
#eta = 0.0005
#v = 0.317
#U = -0.01536996
#W = 0.00146639
#alpha = 0.05

#pomega = 0.0 
#omega1 = omega*(1 + pomega)
#omega2 = omega*(1 - pomega)

#v1 = v 
#v2 = v 

#pU = -0.1
#U1 = U*(1+pU)
#U2 = U*(1-pU)

#pW = 0.1
#W1 = W*(1+pW)
#W2 = W*(1-pW)

#eta1 = eta 
#eta2 = eta 
#alpha1 = alpha 
#alpha2 = alpha 

shift = 0.0012
omega1 = 0.31005326 + shift    
v1 = 0.32164670   
U1 = -0.01938546   
W1 = 0.00289494

omega2 = 0.28832678 + shift   
v2 = 0.30520880   
U2 = -0.01112112   
W2 = 0.00062349

alpha = -0.062
eta = 0.0032
eta1 = eta 
eta2 = eta 
alpha1 = alpha 
alpha2 = alpha 

d0 = 0.35 
dist = 0.1
V = 0.038
beta = -0.3 


##### ============================================================================
##### The band structure in 3D plot 
### The array of intrinsic momenta k
Nk = 201
Kmax = 0.04
k_array = np.linspace(0.0,Kmax,Nk)
dk = (k_array.max() - k_array.min())/(Nk-1)

### The array of synthetic momenta q 
Nq = 201
q_array = np.linspace(0.0,0.2,Nq)
dq = (q_array.max() - q_array.min())/(Nq-1)

### The array of energy 
Energy_array = np.zeros((Nk,Nq,8))

##### We scan over the intrinsic and the synthetic momenta 
for ik in range(Nk):
    ### The genuine momentum 
    k = k_array[ik] 

    for iq in range(Nq):
        ### The synthetic momentum 
        q = q_array[iq]

        ### The Hamiltonian 
        H = Hamiltonian(k,q,omega1,eta1,v1,U1,W1,alpha1,
                 omega2,eta2,v2,U2,W2,alpha2,
                 V,beta,dist,d0) 

        ### Diagonalize the Hamiltonian
        E,states = sla.eigh(H)

        ### We save the energy eigenvalues to the array Energy_array
        Energy_array[ik,iq,:] = E 

### Plot the energy band structures (3D)
# The arrays of domains and colormap
X,Y = np.meshgrid(k_array,q_array)

fig,ax = plt.subplots(subplot_kw={'projection':'3d'},figsize=(10,8))
ax.plot_surface(X,Y,Energy_array[:,:,0].T,cmap='autumn')
ax.plot_surface(X,Y,Energy_array[:,:,1].T,cmap='winter')
#ax.plot_surface(X,Y,Energy_array[:,:,2].T)
#ax.plot_surface(X,Y,Energy_array[:,:,3].T)
ax.set_xlabel('k',fontsize=14)
ax.set_ylabel('q',fontsize=14)
ax.set_title(f'alpha = {alpha:.2f}, eta = {eta:.4f}',fontsize=14)
ax.view_init(elev=5, azim=225, roll=0)
plt.savefig(f'Dispersion-alpha_{alpha:.2f}-eta_{eta:.4f}.png')
#plt.show()
#plt.close() 


##### ==============================================================================
##### The band structure along k for q = 0.5
### The value of q 
q = 0.5 

### The array of k 
Nk = 201 
Kmax = 0.05
k_array = np.linspace(-Kmax,Kmax,Nk)

### The array of energy 
Energy_array = np.zeros((Nk,8))

### We scan over the range of k 
for ik in range(Nk):
    # The genuine momentum k 
    k = k_array[ik]

    # The Hamiltonian 
    H = Hamiltonian(k,q,omega1,eta1,v1,U1,W1,alpha1,
                 omega2,eta2,v2,U2,W2,alpha2,
                 V,beta,dist,d0)

    # The energy eigenvalues 
    E,states = sla.eigh(H)

    # The array of energy 
    Energy_array[ik,:] = E 

### Plot the dispersion curve at q = 0.5
fig,ax = plt.subplots()
ax.plot(k_array,Energy_array[:,0])
ax.plot(k_array,Energy_array[:,1])
ax.plot(k_array,Energy_array[:,2])
ax.plot(k_array,Energy_array[:,3])
ax.set_xlabel('k',fontsize=14)
ax.set_ylabel('E',fontsize=14)
ax.set_title('q = '+str(q),fontsize=14)
#plt.show()
#plt.close()

##### =============================================================================
##### The band structure along q for fixed k 
### The value of k 
k = 0.00

### The array of q 
Nq = 1001 
q_array = np.linspace(-0.3,0.3,Nq)

### The array of energy 
Energy_array = np.zeros((Nq,8))

### We scan over the values of q 
for iq in range(Nq):
    # The synthetic momentum 
    q = q_array[iq]

    # The Hamiltonian 
    H = Hamiltonian(k,q,omega1,eta1,v1,U1,W1,alpha1,
                 omega2,eta2,v2,U2,W2,alpha2,
                 V,beta,dist,d0)

    # The energy eigenvalues 
    E,states = sla.eigh(H)

    # The array of energy 
    Energy_array[iq,:] = E 

### Plot the dispersion curve at k 
fig,ax = plt.subplots()
ax.plot(q_array,Energy_array[:,0])
ax.plot(q_array,Energy_array[:,1])
#ax.plot(q_array,Energy_array[:,2])
#ax.plot(q_array,Energy_array[:,3])
ax.set_xlabel('q',fontsize=14)
ax.set_ylabel('E',fontsize=14)
ax.set_title(f'alpha = {alpha:.2f}, eta = {eta:.4f}, k = {k:.4f}',fontsize=14)
plt.savefig(f'Dispersion-k_0-alpha_{alpha:.2f}-eta_{eta:.4f}.png')
plt.show()
plt.close()