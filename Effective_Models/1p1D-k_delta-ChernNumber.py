import numpy as np 
import scipy 
import scipy.linalg as sla 
import cmath 
import matplotlib.pyplot as plt 
from matplotlib import cm, colors 

##### FUNCTION: Hamiltonian of 2D photonic crystal slab bilayer at M-point
###   The slabs are different from each other
def _1p1D_2DSlab2L_M_k_delta(q,delta,omega1,v1,U1,W1,omega2,v2,U2,W2,V):
    omega1pp = omega1 + v1*q*np.sqrt(2) 
    omega1pm = omega1 
    omega1mp = omega1 
    omega1mm = omega1 - v1*q*np.sqrt(2)

    omega2pp = omega2 + v2*q*np.sqrt(2) 
    omega2pm = omega2 
    omega2mp = omega2  
    omega2mm = omega2 - v2*q*np.sqrt(2) 

    K = 2.0*np.pi 

    Vpp = V * cmath.exp(-1j*K*delta)
    Vpm = V 
    Vmp = V 
    Vmm = V * cmath.exp(1j*K*delta) 

    Hamiltonian = np.array(
        [
            [omega1pp,W1,W1,U1,Vpp,0.0,0.0,0.0],
            [W1,omega1pm,U1,W1,0.0,Vpm,0.0,0.0],
            [W1,U1,omega1mp,W1,0.0,0.0,Vmp,0.0],
            [U1, W1, W1, omega1mm,0.0,0.0,0.0,Vmm],
            [Vmm,0.0,0.0,0.0,omega2pp,W2,W2,U2],
            [0.0,Vmp,0.0,0.0,W2,omega2pm,U2,W2],
            [0.0,0.0,Vpm,0.0,W2,U2,omega2mp,W2],
            [0.0,0.0,0.0,Vpp,U2,W2,W2,omega2mm] 
        ]
    )

    evalues,evectors = sla.eigh(Hamiltonian,eigvals_only=False)

    return evalues, evectors 

##### FUNCTION: The derivative of the Hamiltonian with respect to q 
def dq_1p1D_2DSlab2L_M_k_delta(v1,v2):
    dHdq = np.array(
        [
        [v1*np.sqrt(2.0),0.0,0.0,0.0,0.0,0.0,0.0,0.0],
        [0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0], 
        [0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0],
        [0.0,0.0,0.0,-v1*np.sqrt(2.0),0.0,0.0,0.0,0.0],
        [0.0,0.0,0.0,0.0,v2*np.sqrt(2.0),0.0,0.0,0.0],
        [0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0],
        [0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0],
        [0.0,0.0,0.0,0.0,0.0,0.0,0.0,-v2*np.sqrt(2.0)]
        ]
    )

    return dHdq 

##### FUNCTION: The derivative of the Hamiltonian with respect to delta 
def ddelta_1p1D_2DSlab2L_M_k_delta(delta,V):
    K = 2.0*np.pi 

    Vpp = V * cmath.exp(-1j*K*delta)
    Vmm = V * cmath.exp(1j*K*delta) 

    dHddelta = np.array(
        [
            [0.0,0.0,0.0,0.0, -1j*K*Vpp,0.0,0.0,0.0],
            [0.0,0.0,0.0,0.0,  0.0,0.0,0.0,0.0],
            [0.0,0.0,0.0,0.0,  0.0,0.0,0.0,0.0],
            [0.0,0.0,0.0,0.0,  0.0,0.0,0.0,1j*K*Vmm],
            [1j*K*Vmm,0.0,0.0,0.0, 0.0,0.0,0.0,0.0],
            [0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0],
            [0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0],
            [0.0,0.0,0.0,-1j*K*Vpp, 0.0,0.0,0.0,0.0] 
        ]
    )

    return dHddelta

##### FUNCTION: program to calculate the Berry curvature 
def Berry_curvature_q_delta(q,delta,omega1,v1,U1,W1,omega2,v2,U2,W2,V):
    dHdq = dq_1p1D_2DSlab2L_M_k_delta(v1,v2)
    #print(dHdq) 

    dHddelta = ddelta_1p1D_2DSlab2L_M_k_delta(delta,V)
    #print(dHddelta) 

    evalues,evectors = _1p1D_2DSlab2L_M_k_delta(q,delta,
                                                omega1,v1,U1,W1,
                                                omega2,v2,U2,W2,
                                                V) 

    curvature = np.zeros(8)

    for i in range(8):
        Vector_i = evectors[:,i] 

        for n in range(8):
            if (n != i):
                Vector_n = evectors[:,n]

                S1 = np.matmul( (Vector_i.conjugate()).transpose(),
                                np.matmul(dHdq,Vector_n))
                S2 = np.matmul( (Vector_n.conjugate()).transpose(),
                                np.matmul(dHddelta,Vector_i))
                
                S = -2.0*np.imag(S1*S2)/(evalues[i]-evalues[n])**2 

                curvature[i] = curvature[i] + S 

    return curvature 

##### The MAIN program goes here
def main():
    omega1 = 0.300464408
    omega2 = 0.300464408
    v1 = 0.35
    v2 = 0.35

    U = -0.016615673 
    DeltaU = 0.02
    U1 = U+DeltaU 
    U2 = U-DeltaU 

    W = 0.001744918
    DeltaW = 0.001
    W1 = W+DeltaW 
    W2 = W-DeltaW 

    d0 = 0.35 
    dist = 0.05 
    V0 = 0.044
    V = V0*np.exp(-dist/d0)

    ### The array of q
    Nq = 201 
    q_array = np.linspace(-0.5,0.5,Nq)

    ### The array of delta 
    Ndelta = 201 
    delta_array = np.linspace(-0.5,0.5,Ndelta)

    ### Initialization of the array of Berry curvature
    Curvature_array = np.zeros((Nq,Ndelta,8)) 

    ### We scan over the array of Berry curvature
    for i in range(Nq):
        for j in range(Ndelta):
            q = q_array[i] 
            delta = delta_array[j]

            Curvature_array[i,j,:] = Berry_curvature_q_delta(q,delta,
                                                   omega1,v1,U1,W1,
                                                   omega2,v2,U2,W2,
                                                   V) 


    ### Plot the Berry curvature
    X,Y = np.meshgrid(q_array+0.5,delta_array)

    for i in range(8):
        fig,ax = plt.subplots()
        ax.pcolormesh(X,Y,Curvature_array[:,:,i].T,shading='gouraud',cmap='RdBu')
        ax.set_xlabel('q',fontsize=14)
        ax.set_ylabel(r'$\delta$',fontsize=14)
        ax.set_title('Band '+str(i+1))

        maxabs = abs(Curvature_array[:,:,i]).max()
        vmin = -maxabs 
        vmax = maxabs 
        norm = colors.Normalize(vmin=vmin,vmax=vmax)
        fig.colorbar(cm.ScalarMappable(norm=norm,cmap='RdBu'),
                 ax = ax)
    
    plt.show()

if __name__ == "__main__":
    main()