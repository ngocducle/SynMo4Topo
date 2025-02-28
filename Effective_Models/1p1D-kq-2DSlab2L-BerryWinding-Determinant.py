import numpy as np 
import scipy 
import scipy.linalg as sla 
import cmath 
import matplotlib.pyplot as plt 
from matplotlib import cm,colors 

##### ======================================================================
##### FUNCTION: Hamiltonian of 2D photonic crystal slab bilayer with 
##### kx = ky = k, qx = qy = q 
##### The synthetic momentum: q = delta/(sqrt(2)*a) - 0.5
##### that means q = 0 at shift = sqrt(2)*a/2 
def Hamiltonian(k,q,omega1,eta1,v1,U1,W1,alpha1,
                 omega2,eta2,v2,U2,W2,alpha2,
                 V,beta,dist,d0):
    
    Hamiltonian = np.zeros((8,8),dtype=complex)

    K = 2.0*np.pi 

    ### Block (1,1)
    Hamiltonian[0,0] = omega1 + eta1 + v1*k + v1*k*k/np.sqrt(2)
    Hamiltonian[0,1] = W1 
    Hamiltonian[0,2] = W1 
    Hamiltonian[0,3] = U1*(1+alpha1)

    Hamiltonian[1,0] = W1
    Hamiltonian[1,1] = omega1 - eta1 + v1*k*k/np.sqrt(2)
    Hamiltonian[1,2] = U1*(1-alpha1)
    Hamiltonian[1,3] = W1 

    Hamiltonian[2,0] = W1 
    Hamiltonian[2,1] = U1*(1-alpha1)
    Hamiltonian[2,2] = omega1 - eta1 + v1*k*k/np.sqrt(2)
    Hamiltonian[2,3] = W1 

    Hamiltonian[3,0] = U1*(1+alpha1)
    Hamiltonian[3,1] = W1 
    Hamiltonian[3,2] = W1 
    Hamiltonian[3,3] = omega1 + eta1 - v1*k + v1*k*k/np.sqrt(2)

    ### Block (1,2)
    Hamiltonian[0,4] = -(V + beta*k + beta*k*k/np.sqrt(2))*cmath.exp(-1j*K*q)*np.exp(-dist/d0)
    Hamiltonian[1,5] = (V + beta*k*k/np.sqrt(2))*np.exp(-dist/d0) 
    Hamiltonian[2,6] = (V + beta*k*k/np.sqrt(2))*np.exp(-dist/d0)
    Hamiltonian[3,7] = -(V - beta*k + beta*k*k/np.sqrt(2))*cmath.exp(1j*K*q)*np.exp(-dist/d0)

    ### Block (2,1)
    Hamiltonian[4,0] = -(V + beta*k + beta*k*k/np.sqrt(2))*cmath.exp(1j*K*q)*np.exp(-dist/d0)
    Hamiltonian[5,1] = (V + beta*k*k/np.sqrt(2))*np.exp(-dist/d0)
    Hamiltonian[6,2] = (V + beta*k*k/np.sqrt(2))*np.exp(-dist/d0)
    Hamiltonian[7,3] = -(V - beta*k + beta*k*k/np.sqrt(2))*cmath.exp(-1j*K*q)*np.exp(-dist/d0)

    ### Block (2,2)
    Hamiltonian[4,4] = omega2 + eta2 + v2*k + v2*k*k/np.sqrt(2)
    Hamiltonian[4,5] = W2 
    Hamiltonian[4,6] = W2 
    Hamiltonian[4,7] = U2*(1+alpha2)

    Hamiltonian[5,4] = W2 
    Hamiltonian[5,5] = omega2 - eta2 + v2*k*k/np.sqrt(2)
    Hamiltonian[5,6] = U2*(1-alpha2) 
    Hamiltonian[5,7] = W2 

    Hamiltonian[6,4] = W2 
    Hamiltonian[6,5] = U2*(1-alpha2)
    Hamiltonian[6,6] = omega2 - eta2 + v2*k*k/np.sqrt(2)
    Hamiltonian[6,7] = W2 

    Hamiltonian[7,4] = U2*(1+alpha2)
    Hamiltonian[7,5] = W2 
    Hamiltonian[7,6] = W2 
    Hamiltonian[7,7] = omega2 + eta2 - v2*k + v2*k*k/np.sqrt(2)

    return Hamiltonian

##### ===================================================================
##### The 2x2 Hamiltonian 
def Dispersion_A2x2(A11,A12,A21,A22,k,q):
    Ep = np.sqrt((A11*q+A12*k)*(A11*q+A12*k)+(A21*q+A22*k)*(A21*q+A22*k))
    Em = -Ep 

    return Ep, Em 

##### ====================================================================
##### The arrays of K and Q
Nk = 21 
k_array = np.linspace(-0.1,0.1,Nk)

Nq = 51
q_array = np.linspace(-0.2,0.2,Nq)

X,Y = np.meshgrid(k_array,q_array)

##### The coefficients 
A11 = 1 
A12 = 0.5
A21 = -0.2
A22 = -1

##### The dispersion 
Ep_array,Em_array = Dispersion_A2x2(A11,A12,A21,A22,X,Y)

#print(Ep_array)
#print(Em_array)

### Plot the figure 
fig,ax = plt.subplots(subplot_kw={'projection':'3d'})
ax.plot_surface(X,Y,Ep_array)
ax.plot_surface(X,Y,Em_array)
ax.set_xlabel('k',fontsize=14)
ax.set_ylabel('q',fontsize=14)
plt.show()

##### ======================================================================
##### The parameters of the effective model 
omega = 0.29780940 
eta = 0.0032
v = 0.317 
U = -0.01503117
W = 0.001466 
alpha = -0.062

pomega = 0.0
omega1 = omega*(1+pomega)
omega2 = omega*(1-pomega)

v1 = v
v2 = v 

pU = 0
U1 = U*(1+pU)
U2 = U*(1-pU)

pW = 0
W1 = W*(1+pW)
W2 = W*(1-pW)

eta1 = eta 
eta2 = eta 
alpha1 = alpha 
alpha2 = alpha 

d0 = 0.35 
dist = 0.1 
V0 = 0.038
V = V0*np.exp(-dist/d0)
beta = -0.3

##### ================================================================================
##### The arrays of k and q 
### The array of intrinsic momenta k
Nk = 201 
Kmax = 0.075
k_array = np.linspace(0.030,0.034,Nk)
dk = (k_array.max()-k_array.min())/(Nk-1)

### The array of intrinsic momenta q 
Nq = 201
Qmax = 0.005
q_array = np.linspace(-Qmax,Qmax,Nq)
dq = (q_array.max()-q_array.min())/(Nq-1)

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

X,Y = np.meshgrid(k_array,q_array)

fig,ax = plt.subplots(subplot_kw={'projection':'3d'},figsize=(12,10))
ax.plot_surface(X,Y,Energy_array[:,:,0].T)
ax.plot_surface(X,Y,Energy_array[:,:,1].T)
#ax.plot_surface(X,Y,Energy_array[:,:,2].T)
#ax.plot_surface(X,Y,Energy_array[:,:,3].T)
ax.set_xlabel('k',fontsize=14)
ax.set_ylabel('q',fontsize=14)
plt.show()