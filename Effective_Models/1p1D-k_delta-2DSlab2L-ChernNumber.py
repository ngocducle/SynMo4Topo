import numpy as np 
import scipy 
import scipy.linalg as sla
import cmath 
import matplotlib.pyplot as plt 
from matplotlib import cm,colors 

##### FUNCTION: Hamiltonian of 2D photonic crystal slab bilayer with 
###   kx = ky = k, delta_x = delta_y = delta 
def Hamiltonian(q,delta,omega1,U1,W1,v1,omega2,U2,W2,v2,V):
    Hamiltonian = np.zeros((8,8),dtype=complex)

    K = 2.0*np.pi 

    ### Block (1,1)
    Hamiltonian[0,0] = omega1 + np.sqrt(2.0)*v1*q
    Hamiltonian[0,1] = W1 
    Hamiltonian[0,2] = W1 
    Hamiltonian[0,3] = U1 

    Hamiltonian[1,0] = W1 
    Hamiltonian[1,1] = omega1 
    Hamiltonian[1,2] = U1 
    Hamiltonian[1,3] = W1 

    Hamiltonian[2,0] = W1 
    Hamiltonian[2,1] = U1 
    Hamiltonian[2,2] = omega1  
    Hamiltonian[2,3] = W1  

    Hamiltonian[3,0] = U1 
    Hamiltonian[3,1] = W1 
    Hamiltonian[3,2] = W1 
    Hamiltonian[3,3] = omega1 -np.sqrt(2.0)*v1*q 

    ### Block (1,2)
    Hamiltonian[0,4] = V*cmath.exp(-1j*delta*K)
    Hamiltonian[1,5] = V 
    Hamiltonian[2,6] = V 
    Hamiltonian[3,7] = V*cmath.exp(1j*delta*K) 

    ### Block (2,1)
    Hamiltonian[4,0] = V*cmath.exp(1j*delta*K)
    Hamiltonian[5,1] = V 
    Hamiltonian[6,2] = V 
    Hamiltonian[7,3] = V*cmath.exp(-1j*delta*K)

    ### Block (2,2)
    Hamiltonian[4,4] = omega2 + np.sqrt(2.0)*v2*q 
    Hamiltonian[4,5] = W2 
    Hamiltonian[4,6] = W2 
    Hamiltonian[4,7] = U2 

    Hamiltonian[5,4] = W2 
    Hamiltonian[5,5] = omega2
    Hamiltonian[5,6] = U2 
    Hamiltonian[5,7] = W2 

    Hamiltonian[6,4] = W2 
    Hamiltonian[6,5] = U2 
    Hamiltonian[6,6] = omega2  
    Hamiltonian[6,7] = W2 
    
    Hamiltonian[7,4] = U2 
    Hamiltonian[7,5] = W2 
    Hamiltonian[7,6] = W2 
    Hamiltonian[7,7] = omega2 - np.sqrt(2.0)*v2*q 

    return Hamiltonian

##### FUNCTION: the derivative of the Hamiltonian with respect to 
###   the intrinsic momentum k
def dH_q(v1,v2):
    dHq = np.zeros((8,8),dtype=complex)

    dHq[0,0] = np.sqrt(2.0)*v1 
    dHq[3,3] = -np.sqrt(2.0)*v1 
    dHq[4,4] = np.sqrt(2.0)*v2 
    dHq[7,7] = -np.sqrt(2.0)*v2 

    return dHq 

##### FUNCTION: the derivative of the Hamiltonian with respect to 
###   the synthetic momentum delta 
def dH_delta(V,delta):
    dHdelta = np.zeros((8,8),dtype=complex)

    K = 2.0*np.pi 

    dHdelta[0,4] = -1j*K*V*cmath.exp(-1j*delta*K)
    dHdelta[3,7] = 1j*K*V*cmath.exp(1j*delta*K)
    dHdelta[4,0] = 1j*K*V*cmath.exp(1j*delta*K)
    dHdelta[7,3] = -1j*K*V*cmath.exp(-1j*delta*K)

    return dHdelta 

##### The MAIN program goes here 
def main():
    omega1 = 0.300464408
    omega2 = 0.300464408 
    v1 = 0.35
    v2 = 0.35 

    U = -0.016615673
    pU = -0.00038
    DeltaU = pU*U
    U1 = U+DeltaU 
    U2 = U-DeltaU 

    W = 0.001744918
    pW = -0.01 
    DeltaW = pW*W 
    W1 = W+DeltaW 
    W2 = W-DeltaW 

    d0 = 0.35 
    dist = 0.4
    V0 = 0.044
    V = V0*np.exp(-dist/d0)

    ### The array of intrinsic momenta k
    Nk = 201
    Kmax = 0.1
    k_array = np.linspace(-Kmax,Kmax,Nk)
    dk = (k_array.max() - k_array.min())/(Nk-1)

    ### The array of synthetic momenta delta 
    Ndelta = 201 
    delta_array = np.linspace(-0.5,0.5,Ndelta)
    ddelta = (delta_array.max() - delta_array.min())/(Ndelta-1)

    ### The derivative dH/dk 
    dHq = dH_q(v1,v2)

    ### Arrays of energy 
    Energy_array = np.zeros((Nk,Ndelta,8))

    ### Array of Berry curvature 
    F_array = np.zeros((Nk,Ndelta,8))

    ### Array of Chern number 
    Chern_number = np.zeros(8)

    ##### We scan over the intrinsic and synthetic momenta 
    for i in range(Nk):
        k = k_array[i]
        
        for j in range(Ndelta):
            delta = delta_array[j]
        
            ### The Hamiltonian 
            H = Hamiltonian(k,delta,omega1,U1,W1,v1,omega2,U2,W2,v2,V)
            #print('H = ')
            #print(H)

            ### The derivative dH/ddelta 
            dHdelta = dH_delta(V,delta)

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
            dHqe = np.matmul((states.conjugate()).transpose(),np.matmul(dHq,states))
            dHdeltae = np.matmul((states.conjugate()).transpose(),np.matmul(dHdelta,states))

            for n in range(8):
                for m in range(8):
                    if (m != n):
                        val = -2.0*np.imag(dHqe[n,m]*dHdeltae[m,n]) / (E[n]-E[m])**2
                        F_array[i,j,n] = F_array[i,j,n] + val 

    ### Calculate the Chern numbers 
    Chern_number = np.sum(F_array,axis=(0,1))*dk*ddelta/(2.0*np.pi)

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
                for j in range(Ndelta):
                    file.write('%.8f    ' % k_array[i])
                    file.writelines('%.8f   ' % delta_array[j])
                    file.writelines('%.8f   ' % F_array[i,j,n])
                    file.write('\n')

    ### The maximal absolute value of the Berry curvature of all the 8 bands 
    maxabs = abs(F_array).max()

    ### Plot the 2D maps of the Berry curvature of the 8 bands 
    X,Y = np.meshgrid(k_array+0.5,delta_array)
    cmap = 'coolwarm'

    fig,ax = plt.subplots(2,4,sharex=True,sharey=True,figsize=(16,10))
    vmin,vmax = -maxabs,maxabs 
    norm = colors.Normalize(vmin=vmin,vmax=vmax)

    for i in range(2):
        for j in range(4):
            ax[i,j].pcolormesh(X,Y,F_array[:,:,4*i+j].T,shading='gouraud',cmap=cmap)
            ax[i,j].set_xlabel('k',fontsize=14)
            ax[i,j].set_title('Band '+str(4*i+j+1))
            ax[i,j].set_aspect('equal')

    ax[0,0].set_ylabel(r'$\delta$',fontsize=14)
    ax[1,0].set_ylabel(r'$\delta$',fontsize=14)
    
    fig.colorbar(cm.ScalarMappable(norm=norm,cmap=cmap),
        orientation='vertical',
        shrink=0.5,
        ax=ax)
    
    plt.savefig('Berry_curvature_maps_pU'+str(pU)+'_pW'+str(pW)+'.png')
    
    ### Plot the 2D maps of the absolute value of the Berry curvature of the 8 bands 
    X,Y = np.meshgrid(k_array+0.5,delta_array)
    cmap = 'bone'

    fig,ax = plt.subplots(2,4,sharex=True,sharey=True,figsize=(16,10))
    vmin,vmax = 0.0,maxabs 
    norm = colors.Normalize(vmin=vmin,vmax=vmax)

    for i in range(2):
        for j in range(4):
            ax[i,j].pcolormesh(X,Y,abs(F_array[:,:,4*i+j]).T,shading='gouraud',cmap=cmap)
            ax[i,j].set_xlabel('k',fontsize=14)
            ax[i,j].set_title('Band '+str(4*i+j+1))
            ax[i,j].set_aspect('equal')

    ax[0,0].set_ylabel(r'$\delta$',fontsize=14)
    ax[1,0].set_ylabel(r'$\delta$',fontsize=14)
    
    fig.colorbar(cm.ScalarMappable(norm=norm,cmap=cmap),
        orientation='vertical',
        shrink=0.5,
        ax=ax)
    plt.savefig('Abs_Berry_curvature_maps_pU'+str(pU)+'_pW'+str(pW)+'.png')
    
    
    ### Plot the dispersion surfaces with Berry curvature 
    ### All the 8 bands 
    F_array_3D = F_array 
    maxabs = abs(F_array_3D).max() 
    maxabs = 300 
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
                        rstride=1,
                        cstride=1,
                        facecolors=fcolors,
                        cmap=cmap)
        
    ax.set_xlabel(r'$k a / (2 \ pi)$',fontsize=14)
    ax.set_ylabel(r'$\delta$',fontsize=14)
    ax.set_title(r'$\Delta U/U = $'+str(pU)+', $\Delta W/W = $'+str(pW),fontsize=14)
    fig.colorbar(scamap,
                 orientation='vertical',
                 shrink=0.4,
                 ax = ax)
    ax.view_init(elev=5,azim=75,roll=0)
    plt.savefig('Bands_pU'+str(pU)+'_pW'+str(pW)+'.png')        

    ### Plot the dispersion surfaces of couples of bands 
    for i in range(6):
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
        ax.set_xlabel(r'$k a / (2 \ pi)$',fontsize=14)
        ax.set_ylabel(r'$\delta$',fontsize=14)
        #ax.set_title('Bands '+str(i+1)+'+'+str(i+2),fontsize=14)
        ax.set_title(r'$\Delta U/U = $'+str(pU)+', $\Delta W/W = $'+str(pW),fontsize=14)
        fig.colorbar(scamap,
                     orientation='vertical',
                     shrink=0.4,
                     ax = ax)
        ax.view_init(elev=5,azim=75,roll=0)
        plt.savefig('Bands_'+str(i+1)+'_'+str(i+2)+'_pU'+str(pU)+'_pW'+str(pW)+'.png')     

    ### Show the figures        
    plt.show()

if __name__ == "__main__":
    main()