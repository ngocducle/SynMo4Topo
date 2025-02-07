import numpy as np 
import scipy 
import scipy.linalg as sla
import cmath 
import matplotlib.pyplot as plt 
from matplotlib import cm,colors 

##### ==================================================================================
#####
#####           Calculate the Berry winding for 2D slab bilayer system 
#####           with 1+1D dimension
#####
##### ===================================================================================

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

##### ============================================================================
##### FUNCTION: the derivative of the Hamiltonian with respect to 
###   the genuine momentum k 
def dH_k(k,q,v1,v2,beta,dist,d0):
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

pU = 0.0 
U1 = U*(1+pU)
U2 = U*(1-pU)

pW = 0.0 
W1 = W*(1+pW)
W2 = W*(1-pW)

eta1 = eta 
eta2 = eta 
alpha1 = alpha 
alpha2 = alpha 

d0 = 0.35 
dist = 0.1 
V = 0.038
beta = -0.3

##### ================================================================================
##### The arrays of k and q 
### The array of intrinsic momenta k
Nk = 201 
Kmax = 0.10
k_array = np.linspace(-Kmax,Kmax,Nk)
dk = (k_array.max()-k_array.min())/(Nk-1)

### The array of intrinsic momenta q 
Nq = 201
Qmax = 0.5
q_array = np.linspace(-Qmax,Qmax,Nq)
dq = (q_array.max()-q_array.min())/(Nq-1)

### The array of energy 
Energy_array = np.zeros((Nk,Nq,8))

### The array of Berry curvature 
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
        dHk = dH_k(k,q,v1,v2,beta,dist,d0)

        ### The derivative dH/ddelta 
        dHq = dH_q(k,q,dist,V,d0,beta)

        ### Diagonalize the Hamiltonian
        E,states = sla.eigh(H)

        ### We save the energy eigenvalues to the array Energy_array
        Energy_array[ik,iq,:] = E 

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

#with open('Berry_Curvature.txt','w') as file:
#    np.savetxt(file,F_array[:,:,0],fmt='%.8f')

##### ========================================================================
#####       Calculate the Berry winding
##### ========================================================================
### The center of the circular contour 
center_k = 0.0 
center_q = 0.0 

### The radius of the contour
contour_radius = 0.025

### The array of angles 
Ntheta = 50 
theta_array = np.arange(0,2*np.pi,Ntheta)



##### ========================================================================
###         Plot the Berry curvature map for bands 1 
##### ========================================================================
### The arrays of domains and colormaps 
X,Y = np.meshgrid(k_array,q_array)

norm = colors.Normalize(vmin = -np.max(abs(F_array[:,:,0:2])),
                        vmax = np.max(abs(F_array[:,:,0:2])))

fig,ax = plt.subplots(figsize=(8,8))
ax.pcolormesh(X,Y,F_array[:,:,0].T,shading='gouraud',cmap='coolwarm',norm=norm)
ax.set_xlabel('k',fontsize=20)
ax.set_ylabel('q',fontsize=20)
ax.set_title(r'Band 1, $\alpha = $'+str(alpha)+r', $\eta = $'+str(eta),fontsize=20)
fig.colorbar(cm.ScalarMappable(norm=norm,cmap='coolwarm'),
             orientation='vertical',
             shrink=1.0,
             ax=ax)
plt.savefig('Band1.png')
plt.show()
