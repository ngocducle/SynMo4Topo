import numpy as np 
import scipy 
from scipy.linalg import eigh
import cmath 
import matplotlib.pyplot as plt 
from matplotlib import cm,colors 

##### ===========================================================================
##### FUNCTION: Hamiltonian of 2D photonic crystal slab bilayer with 
##### kx = ky = k, qx = qy = q 
##### The synthetic momentum: q = delta/(sqrt(2)*a) - 0.5
##### that means q = 0 at shift = sqrt(2)*a/2
##### The model is derived directly from the Maxwell's equation 
def Hamiltonian(k,q,omega1,eta1,v1,U1,W1,alpha1,
                 omega2,eta2,v2,U2,W2,alpha2,
                 V,beta,dist,d0):
    
    H = np.zeros((8,8),dtype=complex)

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
##### FUNCTION: the derivative of the Hamiltonian with respect to 
###   the genuine momentum k 
def dH_k(k,q,v1,v2,beta):
    dHk = np.zeros((8,8),dtype=complex)

    K = 2.0*np.pi 

    dHk[0,0] = v1*(1+np.sqrt(2)*k)
    dHk[1,1] = v1*np.sqrt(2)*k 
    dHk[2,2] = v1*np.sqrt(2)*k 
    dHk[3,3] = v1*(-1+np.sqrt(2)*k)

    dHk[0,4] = -np.exp(-dist/d0)*beta*(1+np.sqrt(2)*k)*cmath.exp(-1j*K*q)
    dHk[1,5] = np.exp(-dist/d0)*np.sqrt(2)*beta*k 
    dHk[2,6] = np.exp(-dist/d0)*np.sqrt(2)*beta*k 
    dHk[3,7] = np.exp(-dist/d0)*beta*(1-np.sqrt(2)*k)*cmath.exp(1j*K*q)

    dHk[4,0] = -np.exp(-dist/d0)*beta*(1+np.sqrt(2)*k)*cmath.exp(1j*K*q)
    dHk[5,1] = np.exp(-dist/d0)*np.sqrt(2)*beta*k 
    dHk[6,2] = np.exp(-dist/d0)*np.sqrt(2)*beta*k 
    dHk[7,3] = np.exp(-dist/d0)*beta*(1-np.sqrt(2)*k)*cmath.exp(-1j*K*q)

    dHk[4,4] = v2*(1+np.sqrt(2)*k)
    dHk[5,5] = v2*np.sqrt(2)*k 
    dHk[6,6] = v2*np.sqrt(2)*k 
    dHk[7,7] = v2*(-1+np.sqrt(2)*k)

    return dHk  

##### ============================================================================
##### FUNCTION: the derivative of the Hamiltonian with respect to 
###   the synthetic momentum q
def dH_q(k,q,dist,V,d0,beta):
    dHq = np.zeros((8,8),dtype=complex)

    K = 2.0*np.pi 

    dHq[0,4] = np.exp(-dist/d0)*1j*K*(V+beta*k+beta*k*k/np.sqrt(2))*cmath.exp(-1j*K*q)
    dHq[3,7] = -np.exp(-dist/d0)*1j*K*(V-beta*k+beta*k*k/np.sqrt(2))*cmath.exp(1j*K*q)
    dHq[4,0] = -np.exp(-dist/d0)*1j*K*(V+beta*k+beta*k*k/np.sqrt(2))*cmath.exp(1j*K*q)
    dHq[7,3] = np.exp(-dist/d0)*1j*K*(V-beta*k+beta*k*k/np.sqrt(2))*cmath.exp(-1j*K*q)

    return dHq 

##### ==============================================================================
##### Range of alpha and eta 
Na = 41
alpha_array = np.linspace(-0.08,0.08,Na)

Neta = 51 
eta_array = np.linspace(-0.005,0.005,Neta)

Chern = np.zeros((Na,Neta,8))

##### ============================================================================
##### The parameters of the calculations 
omega = 0.29780940 
v = 0.317
U = -0.01536996
W = 0.00146639

pomega = 0.0 
omega1 = omega*(1 + pomega)
omega2 = omega*(1 - pomega)

v1 = v 
v2 = v 

pU = 0.1
U1 = U*(1+pU)
U2 = U*(1-pU)

pW = -0.1 
W1 = W*(1+pW)
W2 = W*(1-pW)

d0 = 0.35 
dist = 0.1
V = 0.038
beta = -0.3 

##### ================================================================================
##### The arrays of k and q 
### The array of intrinsic momenta k
Nk = 101 
Kmax = 0.05
k_array = np.linspace(-Kmax,Kmax,Nk)
dk = (k_array.max()-k_array.min())/(Nk-1)

### The array of intrinsic momenta q 
Nq = 101
Qmax = 0.5
q_array = np.linspace(-Qmax,Qmax,Nq)
dq = (q_array.max()-q_array.min())/(Nq-1)

### The array of energy 
#Energy_array = np.zeros((Nk,Nq,8))

######## ===============================================================================
######## Scan over alpha and eta 
for ia in range(Na):
    for ieta in range(Neta):
        ### Take the values of alpha and eta 
        alpha = alpha_array[ia] 
        eta = eta_array[ieta]

        ### Values of alpha and eta for the 2 layers 
        eta1 = eta 
        eta2 = eta 
        alpha1 = alpha 
        alpha2 = alpha 

        ### (Re)initialize the array of Berry curvature 
        F_array = np.zeros((Nk,Nq,8))

        ##### We scan over the intrinsic and the synthetic momenta 
        for ik in range(Nk):
            ### The genuine momentum 
            k = k_array[ik] 

            for iq in range(Nq):
                ### The synthetic momentum 
                q = q_array[iq]

                ### The Hamiltonian 
                H = Hamiltonian(k,q,omega1,eta1,v1,U1,W1,alpha1,omega2,eta2,v2,U2,W2,alpha2,
                     V,beta,dist,d0) 
        
                ### The derivative dH/dk 
                dHk = dH_k(k,q,v1,v2,beta)

                ### The derivative dH/ddelta 
                dHq = dH_q(k,q,dist,V,d0,beta)

                ### Diagonalize the Hamiltonian
                E,states = eigh(H)

                ### We save the energy eigenvalues to the array Energy_array
                #Energy_array[ik,iq,:] = E 

                ### ATTENTION! The formula to evaluate the Berry curvature is:
                #
                # F_{kq}^n = \sum_{m \ne n} (-2)*<n| dHk |m><m | dHq |n>/(En-Em)^2
                #
                # In fact: <n| dHk |m> and <m| dHq |n> are the matrix elements of 
                # the operators dHk and dHq in the basis of the energy eigenstates 
                # of the Hamiltonian
                #
                # Therefore, we reexpress the matrices dHk and dHq in the basis of 
                # the eigenstates. The transformation is done by the formula:
                #
                #   A' = states^{\dagger}*A*states 
                #
                # here A = dHk or dHq 
                # and the j-th column of states is the eigenvector corresponding
                # to the j-th eigenvalue

                dHke = np.matmul((states.conjugate()).transpose(),np.matmul(dHk,states))
                dHqe = np.matmul((states.conjugate()).transpose(),np.matmul(dHq,states))

                for n in range(8):
                    for m in range(8):
                        if (m != n):
                            val = -2.0*np.imag(dHke[n,m]*dHqe[m,n]) / (E[n]-E[m])**2
                            F_array[ik,iq,n] = F_array[ik,iq,n] + val

        ##### The Chern number 
        for band in range(8):
            Chern[ia,ieta,band] = np.sum(F_array[:,:,band])*dk*dq/(2*np.pi)

        #Chern[ia,ieta,:] = np.sum(F_array,axis=(0,1))/(2*np.pi)

##### Print the Chern number 
#print(Chern)

##### Save the data 
with open('alpha.npy','wb') as f:
    np.save(f,alpha_array)

with open('eta.npy','wb') as f:
    np.save(f,eta_array)

with open('ChernMap.npy','wb') as f:
    np.save(f,Chern) 

##### Plot the map of Chern number 
X,Y = np.meshgrid(alpha_array,eta_array)
norm = colors.Normalize(vmin=-2,vmax=2)
cmap = 'coolwarm'

fig,ax = plt.subplots(2,4,sharex=True,sharey=True,figsize=(24,12))
ax[0,0].pcolormesh(X,Y,Chern[:,:,0].T,shading='gouraud',cmap=cmap,norm=norm)
ax[0,1].pcolormesh(X,Y,Chern[:,:,1].T,shading='gouraud',cmap=cmap,norm=norm)
ax[0,2].pcolormesh(X,Y,Chern[:,:,2].T,shading='gouraud',cmap=cmap,norm=norm)
ax[0,3].pcolormesh(X,Y,Chern[:,:,3].T,shading='gouraud',cmap=cmap,norm=norm)
ax[1,0].pcolormesh(X,Y,Chern[:,:,4].T,shading='gouraud',cmap=cmap,norm=norm)
ax[1,1].pcolormesh(X,Y,Chern[:,:,5].T,shading='gouraud',cmap=cmap,norm=norm)
ax[1,2].pcolormesh(X,Y,Chern[:,:,6].T,shading='gouraud',cmap=cmap,norm=norm)
ax[1,3].pcolormesh(X,Y,Chern[:,:,7].T,shading='gouraud',cmap=cmap,norm=norm)
ax[1,0].set_xlabel('alpha',fontsize=14)
ax[1,1].set_xlabel('alpha',fontsize=14)
ax[1,2].set_xlabel('alpha',fontsize=14)
ax[1,3].set_xlabel('alpha',fontsize=14)
ax[0,0].set_ylabel('eta',fontsize=14)
ax[1,0].set_ylabel('eta',fontsize=14)
ax[0,0].set_title('Band 1',fontsize=14)
ax[0,1].set_title('Band 2',fontsize=14)
ax[0,2].set_title('Band 3',fontsize=14)
ax[0,3].set_title('Band 4',fontsize=14)
ax[1,0].set_title('Band 5',fontsize=14)
ax[1,1].set_title('Band 6',fontsize=14)
ax[1,2].set_title('Band 7',fontsize=14)
ax[1,3].set_title('Band 8',fontsize=14)
fig.colorbar(cm.ScalarMappable(norm=norm,cmap=cmap),
             ax = ax) 
plt.savefig('ChernMap.png')
plt.show()