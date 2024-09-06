import numpy as np 
import scipy 
import scipy.linalg as sla 
import cmath 
import matplotlib.pyplot as plt 
from matplotlib import cm,colors 

##### ===========================================================================
##### FUNCTION: Hamiltonian of 2D photonic crystal slab bilayer with 
##### kx = ky = k, qx = qy = q 
##### The synthetic momentum: q = delta/(sqrt(2)*a) - 0.5
##### that means q = 0 at shift = sqrt(2)*a/2
def Hamiltonian(k,q,omega1,eta1,v1,U1,W1,alpha1,omega2,eta2,v2,U2,W2,alpha2,V):
    Hamiltonian = np.zeros((8,8),dtype=complex)

    U1p = U1*(1+alpha1)
    U1m = U1*(1-alpha1)
    U2p = U2*(1+alpha2)
    U2m = U2*(1-alpha2)

    K = 2.0*np.pi 

    ### Block (1,1)
    Hamiltonian[0,0] = omega1 + v1*k + eta1 + v1*k**2/(2*K)
    Hamiltonian[0,1] = W1 
    Hamiltonian[0,2] = W1 
    Hamiltonian[0,3] = U1*(1+alpha1) 

    Hamiltonian[1,0] = W1
    Hamiltonian[1,1] = omega1 - eta1 + v1*k**2/(2*K)
    Hamiltonian[1,2] = U1*(1-alpha1) 
    Hamiltonian[1,3] = W1 

    Hamiltonian[2,0] = W1 
    Hamiltonian[2,1] = U1*(1-alpha1)
    Hamiltonian[2,2] = omega1 - eta1 + v1*k**2/(2*K)
    Hamiltonian[2,3] = W1 

    Hamiltonian[3,0] = U1*(1+alpha1) 
    Hamiltonian[3,1] = W1 
    Hamiltonian[3,2] = W1 
    Hamiltonian[3,3] = omega1 - v1*k + eta1 + v1*k**2/(2*K)

    ### Block (1,2)
    Hamiltonian[0,4] = V
    Hamiltonian[1,5] = V
    Hamiltonian[2,6] = V 
    Hamiltonian[3,7] = V  

    ### Block (2,1)
    Hamiltonian[4,0] = V 
    Hamiltonian[5,1] = V 
    Hamiltonian[6,2] = V
    Hamiltonian[7,3] = V

    ### Block (2,2)
    Hamiltonian[4,4] = omega2 + v2*k + eta2 + v2*k**2/(2*K)
    Hamiltonian[4,5] = -W2*cmath.exp(-1j*K*q) 
    Hamiltonian[4,6] = -W2*cmath.exp(-1j*K*q) 
    Hamiltonian[4,7] = U2*(1+alpha2)*cmath.exp(-1j*2*K*q) 

    Hamiltonian[5,4] = -W2*cmath.exp(1j*K*q)
    Hamiltonian[5,5] = omega2 - eta2 + v2*k**2/(2*K) 
    Hamiltonian[5,6] = U2*(1-alpha2) 
    Hamiltonian[5,7] = -W2*cmath.exp(-1j*K*q) 

    Hamiltonian[6,4] = -W2*cmath.exp(1j*K*q) 
    Hamiltonian[6,5] = U2*(1-alpha2) 
    Hamiltonian[6,6] = omega2 - eta2 + v2*k**2/(2*K) 
    Hamiltonian[6,7] = -W2*cmath.exp(-1j*K*q) 

    Hamiltonian[7,4] = U2*(1+alpha2)*cmath.exp(1j*2*K*q) 
    Hamiltonian[7,5] = -W2*cmath.exp(1j*K*q) 
    Hamiltonian[7,6] = -W2*cmath.exp(1j*K*q) 
    Hamiltonian[7,7] = omega2 - v2*k + eta2 + v2*k**2/(2*K)


    return Hamiltonian 

##### ============================================================================
##### FUNCTION: the derivative of the Hamiltonian with respect to 
###   the genuine momentum k 
def dH_k(k,v1,v2):
    dHk = np.zeros((8,8),dtype=complex)

    K = 2.0*np.pi 

    dHk[0,0] = v1 + v1*k/K 
    dHk[1,1] = v1*k/K 
    dHk[2,2] = v1*k/K    
    dHk[3,3] = -v1 + v1*k/K  

    dHk[4,4] = v2 + v2*k/K 
    dHk[5,5] = v2*k/K 
    dHk[6,6] = v2*k/K 
    dHk[7,7] = -v2 + v2*k/K  

    return dHk  

##### ============================================================================
##### FUNCTION: the derivative of the Hamiltonian with respect to 
###   the synthetic momentum q
def dH_q(k,q,W2,alpha2,V):
    dHq = np.zeros((8,8),dtype=complex)

    K = 2.0*np.pi 

    dHq[4,5] = 1j*K*W2*cmath.exp(-1j*K*q)
    dHq[4,6] = 1j*K*W2*cmath.exp(-1j*K*q)
    dHq[4,7] = -1j*2*K*U2*(1+alpha2)*cmath.exp(-1j*2*K*q)

    dHq[5,4] = -1j*K*W2*cmath.exp(1j*K*q)
    dHq[5,7] = 1j*K*W2*cmath.exp(-1j*K*q) 

    dHq[6,4] = -1j*K*W2*cmath.exp(1j*K*q)
    dHq[6,7] = 1j*K*W2*cmath.exp(-1j*K*q)

    dHq[7,4] = 1j*2*K*U2*(1+alpha2)*cmath.exp(1j*2*K*q)
    dHq[7,5] = -1j*K*W2*cmath.exp(1j*K*q)
    dHq[7,6] = -1j*K*W2*cmath.exp(1j*K*q)

    return dHq 

##### ============================================================================
##### The parameters of the calculations 
omega = 0.29780940 
eta = -0.003 
v = 0.317
U = -0.01536996
W = 0.00146639
alpha = 0.01

pomega = 0.0 
omega1 = omega*(1 + pomega)
omega2 = omega*(1 - pomega)

v1 = v 
v2 = v 

pU = 0.01
U1 = U*(1+pU)
U2 = U*(1-pU)

pW = -0.01 
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

##### ================================================================================
##### The arrays of k and q 
### The array of intrinsic momenta k
Nk = 201 
Kmax = 0.05
k_array = np.linspace(-Kmax,Kmax,Nk)
dk = (k_array.max()-k_array.min())/(Nk-1)

### The array of intrinsic momenta q 
Nq = 201
q_array = np.linspace(-0.1,0.1,Nq)
dq = (q_array.max()-q_array.min())/(Nq-1)

### The array of energy 
Energy_array = np.zeros((Nk,Nq,8))

### The array of Berry curvature 
F_array = np.zeros((Nk,Nq,8))

##### We scan over the intrinsic and the synthetic momenta 
for ik in range(Nk):
    ### The genuine momentum 
    k = k_array[ik] 

    ### The derivative dH/dk 
    dHk = dH_k(k,v1,v2)

    for iq in range(Nq):
        ### The synthetic momentum 
        q = q_array[iq]

        ### The Hamiltonian 
        H = Hamiltonian(k,q,omega1,eta1,v1,U1,W1,alpha1,omega2,eta2,v2,U2,W2,alpha2,V) 

        ### The derivative dH/ddelta 
        dHq = dH_q(k,q,W2,alpha2,V)

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

with open('Berry_Curvature.txt','w') as file:
    np.savetxt(file,F_array[:,:,0],fmt='%.8f')
##### =================================================================================
#####       Plot the 2D maps of the Berry curvature of bands 1 and 2 
##### =================================================================================
 
### The arrays of domains and colormaps 
X,Y = np.meshgrid(k_array,q_array)
cmap = 'RdBu'
maxabs = abs(F_array[:,:,0]).max()
#F_array = abs(F_array)
namesave = "Berry_curvature_alpha_{:.4f}.png".format(alpha)

fig,ax = plt.subplots(2,1,sharex=True,figsize=(5,12))
vmin,vmax = F_array[:,:,0].min(),F_array[:,:,0].max()  
norm0 = colors.Normalize(vmin=vmin,vmax=vmax)
vmin,vmax = F_array[:,:,1].min(),F_array[:,:,1].max()
norm1 = colors.Normalize(vmin=vmin,vmax=vmax)  
#vmin,vmax = F_array.min(),F_array.max()
#norm = colors.Normalize(vmin=vmin,vmax=vmax)

ax[0].pcolormesh(X,Y,F_array[:,:,0].T,shading='gouraud',cmap='Blues')
ax[0].set_title('Band 1',fontsize=14)
#ax[0].set_xticks([-0.02,-0.01,0,0.01,0.02])
#ax[0].set_yticks([-0.04,-0.02,0.0,0.02,0.04])
ax[0].set_xlabel('k',fontsize=16)
#ax[0].set_ylabel('q',fontsize=14)

ax[1].pcolormesh(X,Y,F_array[:,:,1].T,shading='gouraud',cmap='Reds')
#ax[1].set_xticks([-0.02,-0.01,0,0.01,0.02])
#ax[1].set_yticks([-0.04,-0.02,0.0,0.02,0.04])
ax[1].set_title('Band 2',fontsize=14)
ax[1].set_xlabel('k',fontsize=16)
#ax[1].set_ylabel('q',fontsize=14)

ax[0].tick_params(axis='y',labelsize=14)
ax[1].tick_params(axis='y',labelsize=14)
ax[1].tick_params(axis='x',labelsize=14)

fig.colorbar(cm.ScalarMappable(norm=norm0,cmap='Blues'),
             orientation='vertical',
             shrink=1.0,
             ax=ax[0])

fig.colorbar(cm.ScalarMappable(norm=norm1,cmap='Reds'),
             orientation='vertical',
             shrink=1.0,
             ax=ax[1])

#plt.colorbar()

plt.savefig(namesave)
plt.show()
plt.close()

##### ========================================================================
###         Plot the dispersion surfaces with Berry curvature 
##### ========================================================================
##### For band 1 and band 2 only 
F_array_3D = F_array 
maxabs = abs(F_array_3D[:,:,0]).max()
vmin, vmax = -maxabs,maxabs
cmap = 'coolwarm'
norm = colors.Normalize(vmin=vmin,vmax=vmax)
scamap = plt.cm.ScalarMappable(norm=norm,cmap=cmap)
linewidth = 0 

fig,ax = plt.subplots(subplot_kw = {'projection':'3d'},
                      figsize=(12,10))

fcolors = scamap.to_rgba(F_array_3D[:,:,0].T)
ax.plot_surface(X,Y,Energy_array[:,:,0].T,
                linewidth=linewidth,
                antialiased='False',
                rstride=1,
                cstride=1,
                facecolors=fcolors,
                cmap=cmap)

fcolors = scamap.to_rgba(F_array_3D[:,:,1].T)
ax.plot_surface(X,Y,Energy_array[:,:,1].T,
                linewidth=linewidth,
                antialiased='False',
                rstride=1,
                cstride=1,
                facecolors=fcolors,
                cmap=cmap)

#ax.set_xticks([-0.02,-0.01,0.0,0.01,0.02])
#ax.set_yticks([-0.04,-0.02,0.0,0.02,0.04])
ax.set_xlabel('k',fontsize=14)
ax.set_ylabel('q',fontsize=14)
ax.set_zlabel(r'$\omega a / (2 \pi c)$',fontsize=14)
fig.colorbar(scamap,
             orientation='vertical',
             shrink=0.4,
             location = 'right',
             ax=ax)
ax.view_init(elev=10,azim=135,roll=0)
plt.savefig('Dis-BCur-'+namesave)
plt.show()
plt.close()

##### ========================================================================
###         Plot the Berry curvature 
##### ========================================================================
##### For band 1 and band 2 only 
fig,ax = plt.subplots(subplot_kw = {'projection':'3d'},
                      figsize=(12,10))

ax.plot_surface(X,Y,F_array[:,:,0].T,
                linewidth=linewidth,
                antialiased='False',
                rstride=1,
                cstride=1)
#ax.set_xticks([-0.02,-0.01,0.0,0.01,0.02])
#ax.set_yticks([-0.04,-0.02,0.0,0.02,0.04])
ax.set_xlabel('k',fontsize=14)
ax.set_ylabel('q',fontsize=14)
ax.set_zlabel(r'Berry curvature (a.u.)',fontsize=14)
plt.savefig('BCur-'+namesave)
plt.show()