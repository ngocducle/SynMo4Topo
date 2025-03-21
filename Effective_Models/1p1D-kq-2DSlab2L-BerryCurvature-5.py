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
##### The model is derived directly from the Maxwell's equation 
##### We include the parameter m describing the degree of inversion symmetry
##### breaking 
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

##### ============================================================================
##### The parameters of the calculations
data = np.loadtxt('2DSlab1L-SHole-ParameterRetrieval-Fit_m-MATLAB.txt')

alpha = 0.062
eta = -0.0032
eta1 = eta 
eta2 = eta 
alpha1 = alpha 
alpha2 = alpha 

for i0 in range(21):

    omega = data[i0,2]
    v = data[i0,3]   
    U = data[i0,4]   
    W = data[i0,5]

    ro = 0.1353 
    rv = 0.0892  
    rW = 2.4215

    m = data[i0,1]
    mo = ro*m 
    mv = rv*m 
    mU = m 
    mW = rW*m 

    shift = 0.0012
    omega1 = (omega + shift)*(1 + mo)  
    v1 = v*(1 + mv)    
    U1 = U*(1 + mU)    
    W1 = W*(1 + mW) 

    omega2 = (omega + shift)*(1 - mo)    
    v2 = v*(1 - mv)    
    U2 = U*(1 - mU)    
    W2 = W*(1 - mW) 

    print('m = '+str(m))
    print('omega1 = '+str(omega1))
    print('v1 = '+str(v1))
    print('U1 = '+str(U1))
    print('W1 = '+str(W1))

    print('omega2 = '+str(omega2))
    print('v2 = '+str(v2))
    print('U2 = '+str(U2))
    print('W2 = '+str(W2))

    d0 = 0.35 
    dist = 0.1
    V = 0.039
    beta = -0.3 

    ##### ================================================================================
    ##### The arrays of k and q 
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
                for l in range(8):
                    if (l != n):
                        val = -2.0*np.imag(dHke[n,l]*dHqe[l,n]) / (E[n]-E[l])**2
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
    namesave = "Berry_curvature_m_{:.4f}.png".format(m)

    ### The arrays of domains and colormaps 
    fig,axs = plt.subplots(2,1,sharex=True,figsize=(5,12))

    #F_array = np.tanh(F_array) 

    norm = colors.Normalize(vmin = -0.04*np.max(abs(F_array[:,:,0:2])),
                        vmax = 0.04*np.max(abs(F_array[:,:,0:2])))

    images = []

    images.append(axs[0].imshow(np.flipud(F_array[:,:,0].T),
                            cmap='coolwarm',
                            norm=norm,
                            aspect=1.25))
    images.append(axs[1].imshow(np.flipud(F_array[:,:,1].T),
                            cmap='coolwarm',
                            norm=norm,
                            aspect=1.25))

    axs[1].set_xticks([0,50,100,150,200])
    axs[1].set_xticklabels([-Kmax,-0.5*Kmax,0,0.5*Kmax,Kmax],fontsize=15)
    axs[0].set_yticks([0,50,100,150,200])
    axs[0].set_yticklabels([Qmax,0.5*Qmax,0,-0.5*Qmax,-Qmax],fontsize=15)
    axs[1].set_yticks([0,50,100,150,200])
    axs[1].set_yticklabels([Qmax,0.5*Qmax,0,-0.5*Qmax,-Qmax],fontsize=15)
    #axs[0].set_title(r'$\alpha=$'+str(alpha)+r', $\eta = $'+str(eta),fontsize=15)
    axs[0].set_title('m = {:.4f}'.format(m),fontsize=24)

    plt.xlabel('k',fontsize=16)
    #plt.ylabel('q',fontsize=16)

    fig.colorbar(images[0],
             ax=axs,
             orientation='vertical',
             shrink=1.0,
             location='right')

    plt.savefig('imshow-'+namesave)
    #plt.show()

    ##### ========================================================================
    ###         Plot the Berry curvature map for bands 1 and 2 separately
    ##### ========================================================================
    norm = colors.Normalize(vmin = -np.max(abs(F_array[:,:,0:2])),
                        vmax = np.max(abs(F_array[:,:,0:2])))

    fig,ax = plt.subplots(figsize=(8,8))
    ax.pcolormesh(X,Y,F_array[:,:,0].T,shading='gouraud',cmap='coolwarm',norm=norm)
    ax.set_xticks([-Kmax,-0.5*Kmax,0,0.5*Kmax,Kmax])
    ax.set_xticklabels([-Kmax,-0.5*Kmax,0,0.5*Kmax,Kmax],fontsize=15)
    ax.set_yticks([Qmax,0.5*Qmax,0,-0.5*Qmax,-Qmax])
    ax.set_yticklabels([Qmax,0.5*Qmax,0,-0.5*Qmax,-Qmax],fontsize=15)
    ax.set_xlabel('k',fontsize=20)
    ax.set_ylabel('q',fontsize=20)
    #ax.set_title(r'Band 1, $\alpha = $'+str(alpha)+r', $\eta = $'+str(eta),fontsize=20)
    ax.set_title('Band 1: m = {:.4f}'.format(m),fontsize=20)
    fig.colorbar(cm.ScalarMappable(norm=norm,cmap='coolwarm'),
             orientation='vertical',
             shrink=1.0,
             ax=ax)
    #plt.savefig('Band1.png')
    #plt.show()

    fig,ax = plt.subplots(figsize=(8,8))
    ax.pcolormesh(X,Y,F_array[:,:,1].T,shading='gouraud',cmap='coolwarm',norm=norm)
    ax.set_xticks([-Kmax,-0.5*Kmax,0,0.5*Kmax,Kmax])
    ax.set_xticklabels([-Kmax,-0.5*Kmax,0,0.5*Kmax,Kmax],fontsize=15)
    ax.set_yticks([Qmax,0.5*Qmax,0,-0.5*Qmax,-Qmax])
    ax.set_yticklabels([Qmax,0.5*Qmax,0,-0.5*Qmax,-Qmax],fontsize=15)
    ax.set_xlabel('k',fontsize=20)
    ax.set_ylabel('q',fontsize=20)
    #ax.set_title(r'Band 2, $\alpha = $'+str(alpha)+r', $\eta = $'+str(eta),fontsize=20)
    ax.set_title('Band 2: m = {:.4f}'.format(m),fontsize=20)
    fig.colorbar(cm.ScalarMappable(norm=norm,cmap='coolwarm'),
             orientation='vertical',
             shrink=1.0,
             ax=ax)
    #plt.savefig('Band2.png')
    #plt.show()