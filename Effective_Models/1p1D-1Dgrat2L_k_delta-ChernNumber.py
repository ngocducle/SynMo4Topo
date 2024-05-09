import numpy as np 
import scipy 
import scipy.linalg as sla
import cmath 
import matplotlib.pyplot as plt
from matplotlib import cm, colors 

##### FUNCTION: Hamiltonian of 1D grating bilayer at X-point
###   The slabs are different from each other
def _1p1D_1Dgrat2L_X_k_delta(q,delta,v,U1,U2,V):
    Vp = V*cmath.exp(-1j*np.pi*delta) 
    Vm = V*cmath.exp(1j*np.pi*delta)

    Hamiltonian = np.array(
        [
            [v*q,U1,Vp,0.0],
            [U1,-v*q,0.0,Vm],
            [Vm,0.0,v*q,U2],
            [0.0,Vp,U2,-v*q]
        ]
    )

    return Hamiltonian 

##### FUNCTION: The derivative of the Hamiltonian with respect to q 
def dq_1p1D_1Dgrat2L_X_k_delta(v):
    dHdq = np.array(
        [
            [v,0.0,0.0,0.0],
            [0.0,-v,0.0,0.0],
            [0.0,0.0,v,0.0],
            [0.0,0.0,0.0,-v]
        ]
    )

    return dHdq 

##### FUNCTION: The derivative of the Hamiltonian with respect to delta
def ddelta_1p1D_1Dgrat2L_X_k_delta(delta,V):
    Vp = V*cmath.exp(-1j*np.pi*delta) 
    Vm = V*cmath.exp(1j*np.pi*delta)

    dHddelta = np.array(
        [
            [0.0,0.0,-1j*np.pi*Vp,0.0],
            [0.0,0.0,0.0,-1j*np.pi*Vm],
            [1j*np.pi*Vm,0.0,0.0,0.0],
            [0.0,-1j*np.pi*Vp,0.0,0.0] 
        ]
    )

    return dHddelta 

##### FUNCTION: function to calculate the Berry curvature 
def Berry_curvature_q_delta(q,delta,v,U1,U2,V):
    H = _1p1D_1Dgrat2L_X_k_delta(q,delta,v,U1,U2,V) 
    dHdq = dq_1p1D_1Dgrat2L_X_k_delta(v)
    dHddelta = ddelta_1p1D_1Dgrat2L_X_k_delta(delta,V) 

    evalues,evectors = sla.eigh(H,eigvals_only=False) 

    curvature = np.zeros(4)

    for i in range(4):
        Vector_i = evectors[:,i] 

        for n in range(4):
            if (n != i):
                Vector_n = evectors[:,n]

                S1 = np.matmul( (Vector_i.conjugate()).transpose(),
                                np.matmul(dHdq,Vector_n))
                S2 = np.matmul( (Vector_n.conjugate()).transpose(),
                                np.matmul(dHddelta,Vector_i))
                
                S = -2.0*np.imag(S1*S2)/(evalues[i]-evalues[n])**2 

                curvature[i] = curvature[i] + S 

    return evalues,curvature 

##### The MAIN program goes here
def main():
    U = 0.0208
    Delta = -0.05*U
    Gamma = -0.05*U  
    U1 = U+Delta
    U2 = U-Delta  
    V = np.sqrt(U**2-Delta**2)-Gamma 
    ng = 3.1024 
    v = 1.0/(2.0*np.pi*ng)
    

    """
    q = 0.1 
    delta = 0.2 

    F = Berry_curvature_q_delta(q,delta,v,U1,U2,V) 
    print(F)
    """ 

    ### Array of genuine momenta 
    Nk = 501
    k_array = np.linspace(-0.5,0.5,Nk) 

    ### Array of intrinsic momenta
    Ndelta = 201
    delta_array = np.linspace(-0.5,0.5,Ndelta) 

    ### Arrays of energy 
    Energy = np.zeros((Nk,Ndelta,4))

    ### Array of Berry curvature 
    Berry_curvature = np.zeros((Nk,Ndelta,4))

    ### Calculate the Berry curvature 
    for i in range(Nk):
        for j in range(Ndelta):
            q = k_array[i]
            delta = delta_array[j] 
            evalues,Berry_curvature[i,j,:] = Berry_curvature_q_delta(q,delta,v,U1,U2,V) 
            Energy[i,j,:] = evalues 

    print(Berry_curvature)

    ### Plot the figure 
    X,Y = np.meshgrid(k_array,delta_array) 

    for i in range(4):
        fig,ax = plt.subplots()
        ax.pcolormesh(X,Y,Berry_curvature[:,:,i].T,shading='gouraud',cmap='RdBu')
        ax.set_xlabel('q',fontsize=14)
        ax.set_ylabel(r'$\delta$',fontsize=14)
        ax.set_title('Band '+str(i+1))

        maxabs = abs(Berry_curvature[:,:,i]).max()
        vmin = -maxabs 
        vmax = maxabs 
        norm = colors.Normalize(vmin=vmin,vmax=vmax)
        fig.colorbar(cm.ScalarMappable(norm=norm,cmap='RdBu'),
                    ax = ax)
        
        ax = plt.figure(figsize=(12,10)).add_subplot(projection='3d')
        band1 = ax.plot_surface(X,Y,Energy[:,:,0].T)
        band2 = ax.plot_surface(X,Y,Energy[:,:,1].T)
        band3 = ax.plot_surface(X,Y,Energy[:,:,2].T)
        band4 = ax.plot_surface(X,Y,Energy[:,:,3].T)
    
    plt.show()



    

if __name__ == '__main__':
    main()
