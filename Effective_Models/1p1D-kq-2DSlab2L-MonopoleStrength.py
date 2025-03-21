import numpy as np 
import scipy 
import scipy.linalg as sla 
import cmath 
import matplotlib.pyplot as plt 
from matplotlib import cm,colors 

##### =====================================================================
##### FUNCTION: Hamiltonian 
def Hamiltonian(k,q,m,omega,eta,v,U,W,alpha,
                ro,rv,rW,V,beta,dist,d0):
    
    H = np.zeros((8,8),dtype=complex)

    K = 2.0*np.pi 

    omega1 = omega*(1+ro*m)
    eta1 = eta 
    alpha1 = alpha 
    v1 = v*(1+rv*m)
    U1 = U*(1+m)
    W1 = W*(1+rW*m)

    omega2 = omega*(1-ro*m)
    eta2 = eta 
    alpha2 = alpha 
    v2 = v*(1-rv*m)
    U2 = U*(1-m)
    W2 = W*(1-rW*m)

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

##### =====================================================================
##### FUNCTION: The derivative of the Hamiltonian with respect to 
##### the genuine momentum k 
def dH_k(k,q,m,v,rv,beta,dist,d0):
    dHk = np.zeros((8,8),dtype = complex)

    K = 2.0*np.pi 

    v1 = v*(1+rv*m)
    v2 = v*(1-rv*m)

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

##### =====================================================================
##### FUNCTION: The derivative of the Hamiltonian with respect to 
##### the synthetic momentum q 
def dH_q(k,q,dist,V,d0,beta):
    dHq = np.zeros((8,8),dtype=complex)

    K = 2.0*np.pi 

    dHq[0,4] = np.exp(-dist/d0)*1j*K*(V+beta*k+beta*k*k/np.sqrt(2))*cmath.exp(-1j*K*q)
    dHq[3,7] = -np.exp(-dist/d0)*1j*K*(V-beta*k+beta*k*k/np.sqrt(2))*cmath.exp(1j*K*q)
    dHq[4,0] = -np.exp(-dist/d0)*1j*K*(V+beta*k+beta*k*k/np.sqrt(2))*cmath.exp(1j*K*q)
    dHq[7,3] = np.exp(-dist/d0)*1j*K*(V-beta*k+beta*k*k/np.sqrt(2))*cmath.exp(-1j*K*q)

    return dHq

##### =====================================================================
##### FUNCTION: The derivative of the Hamiltonian with respect to 
##### the mass m 
def dH_m(k,omega,v,W,U,alpha,ro,rv,rW):
    dHm = np.zeros((8,8),dtype=complex)

    ### Block (1,1) 
    dHm[0,0] = omega*ro + v*rv*k*k/np.sqrt(2) + v*rv*k
    dHm[0,1] = W*rW 
    dHm[0,2] = W*rW 
    dHm[0,3] = U*(1+alpha) 

    dHm[1,0] = W*rW 
    dHm[1,1] = omega*ro + v*rv*k*k/np.sqrt(2)
    dHm[1,2] = U*(1-alpha)
    dHm[1,3] = W*rW 
    
    dHm[2,0] = W*rW 
    dHm[2,1] = U*(1-alpha)
    dHm[2,2] = omega*ro + v*rv*k*k/np.sqrt(2)
    dHm[2,3] = W*rW 
    
    dHm[3,0] = U*(1+alpha)
    dHm[3,1] = W*rW 
    dHm[3,2] = W*rW 
    dHm[3,3] = omega*ro + v*rv*k*k/np.sqrt(2) - v*rv*k 

    ### Block (2,2)
    dHm[4,4] = -dHm[0,0]
    dHm[4,5] = -dHm[0,1]
    dHm[4,6] = -dHm[0,2]
    dHm[4,7] = -dHm[0,3]

    dHm[5,4] = -dHm[1,0]
    dHm[5,5] = -dHm[1,1]
    dHm[5,6] = -dHm[1,2]
    dHm[5,7] = -dHm[1,3]

    dHm[6,4] = -dHm[2,0]
    dHm[6,5] = -dHm[2,1]
    dHm[6,6] = -dHm[2,2]
    dHm[6,7] = -dHm[2,3]

    dHm[7,4] = -dHm[3,0]
    dHm[7,5] = -dHm[3,1]
    dHm[7,6] = -dHm[3,2]
    dHm[7,7] = -dHm[3,3]

    return dHm   

##### =====================================================================
##### The MAIN program goes here 
#m = -0.10164457   
shift = 0.0012 
omega = 0.29793535 + shift    
v = 0.31429072   
U = -0.01563435   
W = 0.00148791

eta = 0.0032 
alpha = -0.062 

ro = 0.1353 
rv = 0.0892 
rW = 2.4215 

d0 = 0.35 
dist = 0.1 
V = 0.039 
beta = -0.3

##### ========================================================================
##### Plot the Berry curvature F_{kq} in the plane m = const 
### The value of m 
m = 0.1

### The array of genuine momenta k 
Nk = 201 
Kmax = 0.1 
k_array = np.linspace(-Kmax,Kmax,Nk)
dk = (k_array.max()-k_array.min())/(Nk-1)

### The array of synthetic momenta q 
Nq = 201 
Qmax = 0.5 
q_array = np.linspace(-Qmax,Qmax,Nq)
dq = (q_array.max()-q_array.min())/(Nq-1)

### The array of energy 
Energy_array = np.zeros((Nk,Nq,8))

### The arrays of Berry curvature 
Fkq_array = np.zeros((Nk,Nq,8))
Fqm_array = np.zeros((Nk,Nq,8))
Fmk_array = np.zeros((Nk,Nq,8))

##### We scan over the intrinsic and the synthetic momenta 
for ik in range(Nk):
    # The genuine momentum 
    k = k_array[ik]

    for iq in range(Nq):
        # The synthetic momentum 
        q = q_array[iq]

        # The Hamiltonian 
        H = Hamiltonian(k,q,m,omega,eta,v,U,W,alpha,ro,rv,rW,V,beta,dist,d0)

        # The derivative dH/dk 
        dHk = dH_k(k,q,m,v,rv,beta,dist,d0)

        # The derivative dH/dq 
        dHq = dH_q(k,q,dist,V,d0,beta)

        # The derivative dH/dm 
        dHm = dH_m(k,omega,v,W,U,alpha,ro,rv,rW)

        # Diagonalize the Hamiltonian 
        E,states = sla.eigh(H)

        # Convert dHk, dHq, dHm to the basis of energy eigenstates 
        dHke = np.matmul((states.conjugate()).transpose(),np.matmul(dHk,states))
        dHqe = np.matmul((states.conjugate()).transpose(),np.matmul(dHq,states))
        dHme = np.matmul((states.conjugate()).transpose(),np.matmul(dHm,states))

        for n in range(8):
            for l in range(8):
                if (l != n):
                    Fkq_array[ik,iq,n] = Fkq_array[ik,iq,n] -2.0*np.imag(dHke[n,l]*dHqe[l,n]) / (E[n]-E[l])**2
                    Fqm_array[ik,iq,n] = Fqm_array[ik,iq,n] -2.0*np.imag(dHqe[n,l]*dHme[l,n]) / (E[n]-E[l])**2
                    Fmk_array[ik,iq,n] = Fmk_array[ik,iq,n] -2.0*np.imag(dHme[n,l]*dHke[l,n]) / (E[n]-E[l])**2 

##### =======================================================================
#####           Plot the 2D maps of the Berry curvature of bands 1 and 2
##### =======================================================================

### The arrays of domains and colormaps 
X,Y = np.meshgrid(k_array,q_array)
cmap = 'RdBu'
maxabs_kq = abs(Fkq_array[:,:,0]).max()
maxabs_qm = abs(Fqm_array[:,:,0]).max()
maxabs_mk = abs(Fmk_array[:,:,0]).max()

### Plot the Berry curvature 
fig,axs = plt.subplots(1,3,sharey=True,figsize=(15,5))

norm_kq = colors.Normalize(vmin = -np.max(abs(Fkq_array[:,:,0:2])),
                           vmax = np.max(abs(Fkq_array[:,:,0:2])))

norm_qm = colors.Normalize(vmin = -np.max(abs(Fqm_array[:,:,0:2])),
                           vmax = np.max(abs(Fqm_array[:,:,0:2])))

norm_mk = colors.Normalize(vmin = -np.max(abs(Fmk_array[:,:,0:2])),
                           vmax = np.max(abs(Fmk_array[:,:,0:2])))

images = []

images.append(axs[0].imshow(np.flipud(Fkq_array[:,:,0].T),
                            cmap='coolwarm',
                            norm=norm_kq,
                            aspect=1.25))
images.append(axs[1].imshow(np.flipud(Fqm_array[:,:,1].T),
                            cmap='coolwarm',
                            norm=norm_qm,
                            aspect=1.25))
images.append(axs[2].imshow(np.flipud(Fmk_array[:,:,1].T),
                            cmap='coolwarm',
                            norm=norm_mk,
                            aspect=1.25))

axs[0].set_xticks([0,50,100,150,200])
axs[0].set_xticklabels([-Kmax,-0.5*Kmax,0,0.5*Kmax,Kmax],fontsize=15)
axs[1].set_xticks([0,50,100,150,200])
axs[1].set_xticklabels([-Kmax,-0.5*Kmax,0,0.5*Kmax,Kmax],fontsize=15)
axs[2].set_xticks([0,50,100,150,200])
axs[2].set_xticklabels([-Kmax,-0.5*Kmax,0,0.5*Kmax,Kmax],fontsize=15)
axs[0].set_yticks([0,50,100,150,200])
axs[0].set_yticklabels([Qmax,0.5*Qmax,0,-0.5*Qmax,-Qmax],fontsize=15)

plt.show()
