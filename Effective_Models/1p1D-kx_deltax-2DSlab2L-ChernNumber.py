import numpy as np 
import scipy 
import scipy.linalg as sla
import cmath 
import matplotlib.pyplot as plt 
from matplotlib import cm,colors 

##### FUNCTION: Hamiltonian of 2D photonic crystal slab bilayer with 
###   ky = 0, qy = 0 => (kx,qx) 
def Hamiltonian(kx,qx,omega1,U1,W1,v1,omega2,U2,W2,v2,V):
    Hamiltonian = np.zeros((8,8),dtype=complex)

    K = 2.0*np.pi 

    ### Block (1,1)
    Hamiltonian[0,0] = omega1 + v1*kx/np.sqrt(2.0)
    Hamiltonian[0,1] = W1 
    Hamiltonian[0,2] = W1 
    Hamiltonian[0,3] = U1 

    Hamiltonian[1,0] = W1 
    Hamiltonian[1,1] = omega1 + v1*kx/np.sqrt(2.0)
    Hamiltonian[1,2] = U1 
    Hamiltonian[1,3] = W1 

    Hamiltonian[2,0] = W1 
    Hamiltonian[2,1] = U1 
    Hamiltonian[2,2] = omega1 - v1*kx/np.sqrt(2.0) 
    Hamiltonian[2,3] = W1  

    Hamiltonian[3,0] = U1 
    Hamiltonian[3,1] = W1 
    Hamiltonian[3,2] = W1 
    Hamiltonian[3,3] = omega1 -v1*kx/np.sqrt(2.0) 

    ### Block (1,2)
    Hamiltonian[0,4] = V*cmath.exp(-1j*0.5*qx*K)
    Hamiltonian[1,5] = V*cmath.exp(-1j*0.5*qx*K) 
    Hamiltonian[2,6] = V*cmath.exp(1j*0.5*qx*K) 
    Hamiltonian[3,7] = V*cmath.exp(1j*0.5*qx*K) 

    ### Block (2,1)
    Hamiltonian[4,0] = V*cmath.exp(1j*0.5*qx*K)
    Hamiltonian[5,1] = V*cmath.exp(1j*0.5*qx*K) 
    Hamiltonian[6,2] = V*cmath.exp(-1j*0.5*qx*K)
    Hamiltonian[7,3] = V*cmath.exp(-1j*0.5*qx*K)

    ### Block (2,2)
    Hamiltonian[4,4] = omega2 + v2*kx/np.sqrt(2.0) 
    Hamiltonian[4,5] = W2 
    Hamiltonian[4,6] = W2 
    Hamiltonian[4,7] = U2 

    Hamiltonian[5,4] = W2 
    Hamiltonian[5,5] = omega2 + v2*kx/np.sqrt(2.0)
    Hamiltonian[5,6] = U2 
    Hamiltonian[5,7] = W2 

    Hamiltonian[6,4] = W2 
    Hamiltonian[6,5] = U2 
    Hamiltonian[6,6] = omega2 - v2*kx/np.sqrt(2.0) 
    Hamiltonian[6,7] = W2 
    
    Hamiltonian[7,4] = U2 
    Hamiltonian[7,5] = W2 
    Hamiltonian[7,6] = W2 
    Hamiltonian[7,7] = omega2 - v2*kx/np.sqrt(2.0) 

    return Hamiltonian

##### FUNCTION: the derivative of the Hamiltonian with respect to 
###   the intrinsic momentum k
def dH_kx(v1,v2):
    dHkx = np.zeros((8,8),dtype=complex)

    dHkx[0,0] = v1/np.sqrt(2.0)
    dHkx[1,1] = v1/np.sqrt(2.0)
    dHkx[2,2] = -v1/np.sqrt(2.0)
    dHkx[3,3] = -v1/np.sqrt(2.0)   
    dHkx[4,4] = v2/np.sqrt(2.0)
    dHkx[5,5] = v2/np.sqrt(2.0)
    dHkx[6,6] = -v2/np.sqrt(2.0)  
    dHkx[7,7] = -v2/np.sqrt(2.0) 

    return dHkx  

##### FUNCTION: the derivative of the Hamiltonian with respect to 
###   the synthetic momentum delta 
def dH_qx(V,qx):
    dHqx = np.zeros((8,8),dtype=complex)

    K = 2.0*np.pi 

    dHqx[0,4] = -1j*0.5*K*V*cmath.exp(-1j*0.5*qx*K)
    dHqx[1,5] = -1j*0.5*K*V*cmath.exp(-1j*0.5*qx*K)
    dHqx[2,6] = 1j*0.5*K*V*cmath.exp(1j*0.5*qx*K) 
    dHqx[3,7] = 1j*0.5*K*V*cmath.exp(1j*0.5*qx*K)
    dHqx[4,0] = 1j*0.5*K*V*cmath.exp(1j*0.5*qx*K) 
    dHqx[5,1] = 1j*0.5*K*V*cmath.exp(1j*0.5*qx*K) 
    dHqx[6,2] = -1j*0.5*K*V*cmath.exp(-1j*0.5*qx*K)
    dHqx[7,3] = -1j*0.5*K*V*cmath.exp(-1j*0.5*qx*K)

    return dHqx 

##### The MAIN program goes here 
def main():
    omega1 = 0.300464408
    omega2 = 0.300464408 
    v1 = 0.35
    v2 = 0.35 

    U = -0.016615673 
    pU = -0.01
    DeltaU = pU*U
    U1 = U+DeltaU 
    U2 = U-DeltaU 

    W = 0.001744918
    pW = -0.01 
    DeltaW = pW*W 
    W1 = W+DeltaW 
    W2 = W-DeltaW 

    d0 = 0.35 
    dist = 0.05
    V0 = 0.044
    V = V0*np.exp(-dist/d0)

    ### The array of intrinsic momenta k
    Nk = 201
    Kmax = 0.011
    kx_array = np.linspace(-Kmax,Kmax,Nk)
    dkx = (kx_array.max() - kx_array.min())/(Nk-1)

    ### The array of synthetic momenta q 
    Nq = 201 
    q_r = 0.05
    qx_array = np.linspace(0.5-q_r,0.5+q_r,Nq)
    dqx = (qx_array.max() - qx_array.min())/(Nq-1)

    ### The derivative dH/dkx 
    dHkx = dH_kx(v1,v2)

    ### Arrays of energy 
    Energy_array = np.zeros((Nk,Nq,8))

    ### Array of Berry curvature 
    F_array = np.zeros((Nk,Nq,8))

    ### Array of Chern number 
    Chern_number = np.zeros(8)

    ##### We scan over the intrinsic and synthetic momenta 
    for i in range(Nk):
        kx = kx_array[i]
        
        for j in range(Nq):
            qx = qx_array[j]
        
            ### The Hamiltonian 
            H = Hamiltonian(kx,qx,omega1,U1,W1,v1,omega2,U2,W2,v2,V)
            #print('H = ')
            #print(H)

            ### The derivative dH/dqx  
            dHqx = dH_qx(V,qx)

            ### Diagonalize the Hamitonian
            E,states = sla.eigh(H)

            ### We save the energy eigenvalues to the array Energy_array 
            Energy_array[i,j,:] = E 

            ### ATTENTION! The formula to evaluate the Berry curvature is:
            #
            # F_{q delta}^n = \sum_{m \ne n} (-2)*<n| dHq |m><m | dHdelta |n>/(En-Em)^2
            #
            # In fact: <n| dHq |m> and <m| dHdelta |n> are the matrix elements of 
            # the operators dHq and dHdelta in the basis of the energy eigenstates 
            # of the Hamiltonian
            #
            # Therefore, we reexpress the matrices dHq and dHdelta in the basis of 
            # the eigenstates. The transformation is done by the formula:
            #
            #   A' = states^{\dagger}*A*states 
            #
            # here A = dHq or dHdelta
            # and the j-th column of states is the eigenvector corresponding
            # to the j-th eigenvalue
            dHkxe = np.matmul((states.conjugate()).transpose(),np.matmul(dHkx,states))
            dHqxe = np.matmul((states.conjugate()).transpose(),np.matmul(dHqx,states))

            for n in range(8):
                for m in range(8):
                    if (m != n):
                        val = -2.0*np.imag(dHkxe[n,m]*dHqxe[m,n]) / (E[n]-E[m])**2
                        F_array[i,j,n] = F_array[i,j,n] + val 

    ### Calculate the Chern numbers 
    Chern_number = np.sum(F_array,axis=(0,1))*dkx*dqx/(2.0*np.pi)

    print('# Chern number C1 = '+str(Chern_number[0]))
    print('# Chern number C2 = '+str(Chern_number[1]))
    print('# Chern number C3 = '+str(Chern_number[2]))
    print('# Chern number C4 = '+str(Chern_number[3]))
    print('# Chern number C5 = '+str(Chern_number[4]))
    print('# Chern number C6 = '+str(Chern_number[5]))
    print('# Chern number C7 = '+str(Chern_number[6]))
    print('# Chern number C8 = '+str(Chern_number[7]))

    ### Print the results to a file 
    for n in range(8):
        with open('1p1D-k_delta_2DSlab2L-ChernNumber-'+str(n+1)+'.txt','w') as file:
            for i in range(Nk):
                for j in range(Nq):
                    file.write('%.8f    ' % kx_array[i])
                    file.writelines('%.8f   ' % qx_array[j])
                    file.writelines('%.8f   ' % F_array[i,j,n])
                    file.write('\n')

    ### The maximal absolute value of the Berry curvature of all the 8 bands 
    maxabs = abs(F_array).max()

    ### Plot the 2D maps of the Berry curvature of the 8 bands 
    X,Y = np.meshgrid(kx_array+0.5,qx_array)
    cmap = 'coolwarm'

    fig,ax = plt.subplots(2,4,sharex=True,sharey=True,figsize=(16,10))
    vmin,vmax = -maxabs,maxabs 
    norm = colors.Normalize(vmin=vmin,vmax=vmax)

    for i in range(2):
        for j in range(4):
            ax[i,j].pcolormesh(X,Y,F_array[:,:,4*i+j].T,shading='gouraud',cmap=cmap)
            ax[i,j].set_xlabel(r'$k_x$',fontsize=14)
            ax[i,j].set_title('Band '+str(4*i+j+1))
            ax[i,j].set_aspect('equal')

    ax[0,0].set_ylabel(r'$q_x$',fontsize=14)
    ax[1,0].set_ylabel(r'$q_x$',fontsize=14)
    
    fig.colorbar(cm.ScalarMappable(norm=norm,cmap=cmap),
        orientation='vertical',
        shrink=0.5,
        ax=ax)
    
    ### Plot the 2D maps of the absolute value of the Berry curvature of the 8 bands 
    X,Y = np.meshgrid(kx_array+0.5,qx_array)
    cmap = 'bone'

    fig,ax = plt.subplots(2,4,sharex=True,sharey=True,figsize=(16,10))
    vmin,vmax = 0.0,maxabs 
    norm = colors.Normalize(vmin=vmin,vmax=vmax)

    for i in range(2):
        for j in range(4):
            ax[i,j].pcolormesh(X,Y,abs(F_array[:,:,4*i+j]).T,shading='gouraud',cmap=cmap)
            ax[i,j].set_xlabel(r'$k_x$',fontsize=14)
            ax[i,j].set_title('Band '+str(4*i+j+1))
            ax[i,j].set_aspect('equal')

    ax[0,0].set_ylabel(r'$q_x$',fontsize=14)
    ax[1,0].set_ylabel(r'$q_x$',fontsize=14)
    
    fig.colorbar(cm.ScalarMappable(norm=norm,cmap=cmap),
        orientation='vertical',
        shrink=0.5,
        ax=ax)
    
    ### Plot the dispersion surfaces with Berry curvature 
    ### All the 8 bands 
    F_array_3D = F_array 
    maxabs = abs(F_array_3D).max() 
    maxabs = 0.1*maxabs 
    vmin, vmax = -maxabs,maxabs 
    cmap = 'coolwarm'
    norm = colors.Normalize(vmin=vmin,vmax=vmax)
    scamap = plt.cm.ScalarMappable(norm=norm,cmap=cmap)
    linewidth=0 

    xmin = 0  # int(Nk/2) - 20
    xmax = Nk # int(Nk/2) + 20 

    fig,ax = plt.subplots(subplot_kw = {'projection':'3d'},
                          figsize=(12,10)) 

    for i in range(8):
        fcolors = scamap.to_rgba(F_array_3D[:,:,i].T) 
        ax.plot_surface(X[:,xmin:xmax],
                        Y[:,xmin:xmax],
                        Energy_array[xmin:xmax,:,i].T,
                        linewidth=linewidth,
                        antialiased='True',
                        rstride=10,
                        cstride=10,
                        facecolors=fcolors,
                        cmap=cmap)
        
    ax.set_xlabel(r'$k_x$',fontsize=14)
    ax.set_ylabel(r'$q_x$',fontsize=14)
    ax.set_title(r'$\Delta U/U = $'+str(pU)+', $\Delta W/W = $'+str(pW),fontsize=14)
    fig.colorbar(scamap,
                 orientation='vertical',
                 shrink=0.4,
                 ax = ax)
    ax.view_init(elev=10,azim=75,roll=0)
    plt.savefig('Bands'+'_pU'+str(pU)+'_pW'+str(pW)+'.png')        

    ### Plot the dispersion surfaces of couples of bands 
    for i in range(0,8,4):
        fig,ax = plt.subplots(subplot_kw = {'projection':'3d'},
                          figsize=(12,10)) 
        fcolors1 = scamap.to_rgba(F_array_3D[:,:,i].T) 
        ax.plot_surface(X[:,xmin:xmax],
                        Y[:,xmin:xmax],
                        Energy_array[xmin:xmax,:,i].T,
                        linewidth=linewidth,
                        antialiased='True',
                        rstride=1,
                        cstride=1,
                        facecolors=fcolors1,
                        cmap=cmap)
        fcolors2 = scamap.to_rgba(F_array_3D[:,:,i+1].T) 
        ax.plot_surface(X[:,xmin:xmax],
                        Y[:,xmin:xmax],
                        Energy_array[xmin:xmax,:,i+1].T,
                        linewidth=linewidth,
                        antialiased='True',
                        rstride=1,
                        cstride=1,
                        facecolors=fcolors2,
                        cmap=cmap)
        fcolors3 = scamap.to_rgba(F_array_3D[:,:,i+2].T) 
        ax.plot_surface(X[:,xmin:xmax],
                        Y[:,xmin:xmax],
                        Energy_array[xmin:xmax,:,i+2].T,
                        linewidth=linewidth,
                        antialiased='True',
                        rstride=1,
                        cstride=1,
                        facecolors=fcolors3,
                        cmap=cmap)
        fcolors4 = scamap.to_rgba(F_array_3D[:,:,i+3].T) 
        ax.plot_surface(X[:,xmin:xmax],
                        Y[:,xmin:xmax],
                        Energy_array[xmin:xmax,:,i+3].T,
                        linewidth=linewidth,
                        antialiased='True',
                        rstride=1,
                        cstride=1,
                        facecolors=fcolors4,
                        cmap=cmap)
        ax.set_xlabel(r'$k_x$',fontsize=14)
        ax.set_ylabel(r'$q_x$',fontsize=14)
        ax.set_title('Bands-'+str(i+1)+'-'+str(i+2)+'-'+str(i+3)+'-'+str(i+4)+'; '+
                     r'$\Delta U/U = $'+str(pU)+', $\Delta W/W = $'+str(pW),
                     fontsize=14)
        fig.colorbar(scamap,
                     orientation='vertical',
                     shrink=0.4,
                     ax = ax)
        ax.view_init(elev=10,azim=75,roll=0)
        plt.savefig('Bands-'+str(i+1)+'-'+str(i+2)+'-'+str(i+3)+'-'+str(i+4)
                    +'_pU'+str(pU)+'_pW'+str(pW)+'.png')     

    ### Show the figures        
    plt.show()

if __name__ == "__main__":
    main()