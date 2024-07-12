import numpy as np 
import scipy 
import scipy.linalg as sla 
import cmath 
import matplotlib.pyplot as plt 
from matplotlib import cm,colors 

##### ================================================================================
##### FUNCTION: Hamiltonian of 2D photonic crystal slab bilayer with 
###   (kx,ky), qx = qy = q 

def Hamiltonian(kx,ky,q,omega1,U1,W1,v1,alpha1,omega2,U2,W2,v2,alpha2,V):
    Hamiltonian = np.zeros((8,8),dtype=complex)

    K = np.pi 

    U1D = U1*(1+alpha1)
    U1A = U1*(1-alpha1)
    U2D = U2*(1+alpha2)
    U2A = U2*(1-alpha2)

    ### Block (1,1)
    Hamiltonian[0,0] = omega1 + v1*(kx+ky)/np.sqrt(2.0)
    Hamiltonian[0,1] = W1 
    Hamiltonian[0,2] = W1 
    Hamiltonian[0,3] = U1*(1+alpha1) 

    Hamiltonian[1,0] = W1 
    Hamiltonian[1,1] = omega1 + v1*(kx-ky)/np.sqrt(2.0)
    Hamiltonian[1,2] = U1*(1-alpha1) 
    Hamiltonian[1,3] = W1 

    Hamiltonian[2,0] = W1 
    Hamiltonian[2,1] = U1*(1-alpha1) 
    Hamiltonian[2,2] = omega1 + v1*(-kx+ky)/np.sqrt(2.0)
    Hamiltonian[2,3] = W1 

    Hamiltonian[3,0] = U1*(1+alpha1) 
    Hamiltonian[3,1] = W1 
    Hamiltonian[3,2] = W1 
    Hamiltonian[3,3] = omega1 + v1*(-kx-ky)/np.sqrt(2.0) 

    ### Block (1,2)
    Hamiltonian[0,4] = V*cmath.exp(-1j*2*q*K)
    Hamiltonian[1,5] = V
    Hamiltonian[2,6] = V
    Hamiltonian[3,7] = V*cmath.exp(1j*2*q*K) 

    ### Block (2,1)
    Hamiltonian[4,0] = V*cmath.exp(1j*2*q*K)
    Hamiltonian[5,1] = V
    Hamiltonian[6,2] = V
    Hamiltonian[7,3] = V*cmath.exp(-1j*2*q*K)

    ### Block (2,2)
    Hamiltonian[4,4] = omega2 + v2*(kx+ky)/np.sqrt(2.0)
    Hamiltonian[4,5] = W2 
    Hamiltonian[4,6] = W2 
    Hamiltonian[4,7] = U2*(1+alpha2) 

    Hamiltonian[5,4] = W2 
    Hamiltonian[5,5] = omega2 + v2*(kx-ky)/np.sqrt(2.0)
    Hamiltonian[5,6] = U2*(1-alpha2) 
    Hamiltonian[5,7] = W2 

    Hamiltonian[6,4] = W2 
    Hamiltonian[6,5] = U2*(1-alpha2) 
    Hamiltonian[6,6] = omega2 + v2*(-kx+ky)/np.sqrt(2.0)
    Hamiltonian[6,7] = W2 

    Hamiltonian[7,4] = U2*(1+alpha2) 
    Hamiltonian[7,5] = W2 
    Hamiltonian[7,6] = W2 
    Hamiltonian[7,7] = omega2 + v2*(-kx-ky)/np.sqrt(2.0)

    return Hamiltonian 

##### ==================================================================================
###   FUNCTION: the derivative of the Hamiltonian with respect to 
###   the intrinsic momentum kx 
def dH_kx(v1,v2):
    dHkx = np.zeros((8,8),dtype=complex)

    dHkx[0,0] = v1 / np.sqrt(2.0) 
    dHkx[1,1] = v1 / np.sqrt(2.0)
    dHkx[2,2] = - v1 / np.sqrt(2.0)
    dHkx[3,3] = - v1 / np.sqrt(2.0) 
    dHkx[4,4] = v2 / np.sqrt(2.0)
    dHkx[5,5] = v2 / np.sqrt(2.0)
    dHkx[6,6] = - v2 / np.sqrt(2.0)
    dHkx[7,7] = - v2 / np.sqrt(2.0) 

    return dHkx 

##### ==================================================================================
###   FUNCTION: the derivative of the Hamiltonian with respect to 
###   the intrinsic momentum ky 
def dH_ky(v1,v2):
    dHky = np.zeros((8,8),dtype=complex)

    dHky[0,0] = v1 / np.sqrt(2.0) 
    dHky[1,1] = - v1 / np.sqrt(2.0)
    dHky[2,2] = v1 / np.sqrt(2.0)
    dHky[3,3] = - v1 / np.sqrt(2.0) 
    dHky[4,4] = v2 / np.sqrt(2.0)
    dHky[5,5] = - v2 / np.sqrt(2.0)
    dHky[6,6] = v2 / np.sqrt(2.0)
    dHky[7,7] = - v2 / np.sqrt(2.0) 

    return dHky 

##### =============================================================================
##### FUNCTION: the derivative of the Hamiltonian with respect to 
###   the synthetic momentum q 
def dH_q(V,q):
    dHq = np.zeros((8,8),dtype=complex)

    K = 2.0*np.pi 

    dHq[0,4] = -1j*K*V*cmath.exp(-1j*q*K)
    dHq[3,7] = 1j*K*V*cmath.exp(1j*q*K)
    dHq[4,0] = 1j*K*V*cmath.exp(1j*q*K)
    dHq[7,3] = -1j*K*V*cmath.exp(-1j*q*K)

    return dHq 

##### ============================================================================= #
#####                                                                               #
#####                       The MAIN program goes here                              #
#####                                                                               #
##### ============================================================================= # 

def main():
    ### The arrays of genuine momenta 
    Nk = 21
    
    kx_array = np.linspace(-0.1,0.1,Nk)
    dkx = (kx_array.max() - kx_array.min())/(Nk-1)

    ky_array = np.linspace(-0.1,0.1,Nk)
    dky = (ky_array.max() - ky_array.min())/(Nk-1)

    ### The arrays of synthetic momenta 
    Nq = 1
    q_array = np.linspace(0.0,0.0,Nq)
    dq = (q_array.max() - q_array.min())/(Nq-1)

    ### The parameters 
    U = -0.01537 
    dU = 0.01*U 

    W = 0.001466 
    dW = -0.01*W  

    omega1 = 0.29780940
    v1 = 0.317
    U1 = U + dU 
    W1 = W + dW 
    alpha1 = -0.1

    omega2 = 0.29780940
    v2 = 0.317
    U2 = U - dU 
    W2 = W - dW 
    alpha2 = -0.1

    V = 0.025

    ### Arrays of energy 
    Energy_array = np.zeros((Nk,Nk,Nq,8))

    ### Arrays of Berry curvature 
    Fxy_array = np.zeros((Nk,Nk,Nq,8))
    Fyq_array = np.zeros((Nk,Nk,Nq,8))
    Fqx_array = np.zeros((Nk,Nk,Nq,8))

    ### Array of Chern number 
    Chern_number = np.zeros(8)

    ### The derivatives of the Hamiltonian
    dHkx = dH_kx(v1,v2)
    dHky = dH_ky(v1,v2)

    ##### We scan over the intrinsic and synthetic momenta 
    for ix in range(Nk):
        kx = kx_array[ix] 

        for iy in range(Nk):
            ky = ky_array[iy]

            for iq in range(Nq):
                q = q_array[iq]

                ### The Hamiltonian 
                H = Hamiltonian(kx,ky,q,
                                omega1,U1,W1,v1,alpha1,
                                omega2,U2,W2,v2,alpha2,V)
                
                ### The derivative dH/dq 
                dHq = dH_q(V,q)

                ### Diagonalize the Hamiltonian 
                E,states = sla.eigh(H)

                ### We save the energy eigenvalues to the array Energy_array 
                Energy_array[ix,iy,iq,:] = E 

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

                dHkxe = np.matmul((states.conjugate()).transpose(),np.matmul(dHkx,states))
                dHkye = np.matmul((states.conjugate()).transpose(),np.matmul(dHky,states))
                dHqe = np.matmul((states.conjugate()).transpose(),np.matmul(dHq,states))

                for n in range(8):
                    for m in range(8):
                        if (m != n):
                            val = -2.0*np.imag(dHkxe[n,m]*dHkye[m,n])/(E[n]-E[m])**2
                            Fxy_array[ix,iy,iq,n] = Fxy_array[ix,iy,iq,n] + val 

                for n in range(8):
                    for m in range(8):
                        if (m != n):
                            val = -2.0*np.imag(dHkye[n,m]*dHqe[m,n])/(E[n]-E[m])**2
                            Fyq_array[ix,iy,iq,n] = Fyq_array[ix,iy,iq,n] + val 

                for n in range(8):
                    for m in range(8):
                        if (m != n):
                            val = -2.0*np.imag(dHqe[n,m]*dHkxe[m,n])/(E[n]-E[m])**2
                            Fqx_array[ix,iy,iq,n] = Fqx_array[ix,iy,iq,n] + val 

    print(Fxy_array)
    print(Fyq_array)
    print(Fqx_array)

    ##### Plot the Berry flux 
    ### Make the grid 
    X,Y,Q = np.meshgrid(kx_array,ky_array,q_array)

    ### Plot the figure
    ax = plt.figure().add_subplot(projection='3d')
    ax.quiver(X,Y,Q,
              Fyq_array[:,:,:,0],Fqx_array[:,:,:,0],Fxy_array[:,:,:,0],
              length=0.01,normalize=True)
    ax.set_xlabel('kx')
    ax.set_ylabel('ky')
    ax.set_zlabel('q')
    plt.show()

    #for iq in range(Nq):
    #    fig,ax = plt.subplots()
    #    ax.quiver(X,Y,Fyq_array[:,:,iq,0],Fqx_array[:,:,iq,0])
    #    ax.set_xlabel('kx')
    #    ax.set_ylabel('ky')
    #    ax.set_title('q')
    #    plt.show()


##### Run the MAIN program 
if __name__ == '__main__':
    main()

