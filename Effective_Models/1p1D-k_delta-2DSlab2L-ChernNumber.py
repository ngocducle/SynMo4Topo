import numpy as np 
import scipy 
import scipy.linalg as sla
import cmath 
import matplotlib.pyplot as plt 
from matplotlib import cm,colors 

##### ============================================================================
##### FUNCTION: Hamiltonian of 2D photonic crystal slab bilayer with 
###   kx = ky = k, q_x = q_y = q  
def Hamiltonian(k,q,omega1,U1,W1,v1,omega2,U2,W2,v2,V,alpha):
    Hamiltonian = np.zeros((8,8),dtype=complex)

    K = 2.0*np.pi 

    U1D = U1*(1+alpha)
    U1A = U1*(1-alpha)
    U2D = U2*(1+alpha)
    U2A = U2*(1-alpha)

    ### Block (1,1)
    Hamiltonian[0,0] = omega1 + np.sqrt(2.0)*v1*k
    Hamiltonian[0,1] = W1 
    Hamiltonian[0,2] = W1 
    Hamiltonian[0,3] = U1D 

    Hamiltonian[1,0] = W1 
    Hamiltonian[1,1] = omega1 
    Hamiltonian[1,2] = U1A 
    Hamiltonian[1,3] = W1 

    Hamiltonian[2,0] = W1 
    Hamiltonian[2,1] = U1A 
    Hamiltonian[2,2] = omega1  
    Hamiltonian[2,3] = W1  

    Hamiltonian[3,0] = U1D 
    Hamiltonian[3,1] = W1 
    Hamiltonian[3,2] = W1 
    Hamiltonian[3,3] = omega1 -np.sqrt(2.0)*v1*k 

    ### Block (1,2)
    Hamiltonian[0,4] = V*cmath.exp(-1j*q*K)
    Hamiltonian[1,5] = V 
    Hamiltonian[2,6] = V 
    Hamiltonian[3,7] = V*cmath.exp(1j*q*K) 

    ### Block (2,1)
    Hamiltonian[4,0] = V*cmath.exp(1j*q*K)
    Hamiltonian[5,1] = V 
    Hamiltonian[6,2] = V 
    Hamiltonian[7,3] = V*cmath.exp(-1j*q*K)

    ### Block (2,2)
    Hamiltonian[4,4] = omega2 + np.sqrt(2.0)*v2*k 
    Hamiltonian[4,5] = W2 
    Hamiltonian[4,6] = W2 
    Hamiltonian[4,7] = U2D 

    Hamiltonian[5,4] = W2 
    Hamiltonian[5,5] = omega2
    Hamiltonian[5,6] = U2A 
    Hamiltonian[5,7] = W2 

    Hamiltonian[6,4] = W2 
    Hamiltonian[6,5] = U2A 
    Hamiltonian[6,6] = omega2  
    Hamiltonian[6,7] = W2 
    
    Hamiltonian[7,4] = U2D  
    Hamiltonian[7,5] = W2 
    Hamiltonian[7,6] = W2 
    Hamiltonian[7,7] = omega2 - np.sqrt(2.0)*v2*k 

    return Hamiltonian

##### =============================================================================
##### FUNCTION: the derivative of the Hamiltonian with respect to 
###   the intrinsic momentum k
def dH_k(v1,v2):
    dHk = np.zeros((8,8),dtype=complex)

    dHk[0,0] = np.sqrt(2.0)*v1 
    dHk[3,3] = -np.sqrt(2.0)*v1 
    dHk[4,4] = np.sqrt(2.0)*v2 
    dHk[7,7] = -np.sqrt(2.0)*v2 

    return dHk  

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
    omega1 = 0.300464408
    omega2 = 0.300464408 
    v1 = 0.35
    v2 = 0.35 

    # Anisotropic coefficient along the diagonals 
    alpha = -0.10

    U = -0.016615673
    pU = 0.00 
    DeltaU = pU*U
    U1 = U+DeltaU 
    U2 = U-DeltaU 

    W = 0.001744918
    pW = 0.00 
    DeltaW = pW*W 
    W1 = W+DeltaW 
    W2 = W-DeltaW 

    d0 = 0.35 
    dist = 0.4
    V0 = 0.044
    V = V0*np.exp(-dist/d0)

    ### The array of intrinsic momenta k
    Nk = 201
    Kmax = 0.025
    k_array = np.linspace(-Kmax,Kmax,Nk)
    dk = (k_array.max() - k_array.min())/(Nk-1)

    ### The array of synthetic momenta delta 
    Nq = 201 
    q_array = np.linspace(0.50,0.55,Nq)
    dq = (q_array.max() - q_array.min())/(Nq-1)

    ### The derivative dH/dk 
    dHk = dH_k(v1,v2)

    ### Arrays of energy 
    Energy_array = np.zeros((Nk,Nq,8))

    ### Array of Berry curvature 
    F_array = np.zeros((Nk,Nq,8))

    ### Array of Chern number 
    Chern_number = np.zeros(8)

    ##### We scan over the intrinsic and synthetic momenta 
    for i in range(Nk):
        k = k_array[i]
        
        for j in range(Nq):
            q = q_array[j]
        
            ### The Hamiltonian 
            H = Hamiltonian(k,q,omega1,U1,W1,v1,omega2,U2,W2,v2,V,alpha)
            #print('H = ')
            #print(H)

            ### The derivative dH/ddelta 
            dHq = dH_q(V,q)

            ### Diagonalize the Hamitonian
            E,states = sla.eigh(H)

            ### We save the energy eigenvalues to the array Energy_array 
            Energy_array[i,j,:] = E 

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
                        F_array[i,j,n] = F_array[i,j,n] + val 



    ##### ============================================================================
    ###                         Calculate the Chern numbers 
    ##### ============================================================================
    Chern_number = np.sum(F_array,axis=(0,1))*dk*dq/(2.0*np.pi)

    print('# Chern number C1 = '+str(Chern_number[0]))
    print('# Chern number C2 = '+str(Chern_number[1]))
    print('# Chern number C3 = '+str(Chern_number[2]))
    print('# Chern number C4 = '+str(Chern_number[3]))
    print('# Chern number C5 = '+str(Chern_number[4]))
    print('# Chern number C6 = '+str(Chern_number[5]))
    print('# Chern number C7 = '+str(Chern_number[6]))
    print('# Chern number C8 = '+str(Chern_number[7]))

    ##### =============================================================================
    ###                         Print the results to a file 
    ##### =============================================================================
    for n in range(8):
        with open('1p1D-k_delta_2DSlab2L-ChernNumber-'+str(n+1)+'.txt','w') as file:
            for i in range(Nk):
                for j in range(Nq):
                    file.write('%.8f    ' % k_array[i])
                    file.writelines('%.8f   ' % q_array[j])
                    file.writelines('%.8f   ' % F_array[i,j,n])
                    file.write('\n')


    ##### =============================================================================
    #####           Plot the 2D maps of the Berry curvature of the 8 bands 
    ##### ============================================================================= 
    ### The maximal absolute value of the Berry curvature of all the 8 bands 
    maxabs = abs(F_array).max()

    ### The arrays of domains and colormap 
    X,Y = np.meshgrid(k_array,q_array)
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

    ax[0,0].set_ylabel('q',fontsize=14)
    ax[1,0].set_ylabel('q',fontsize=14)
    
    fig.colorbar(cm.ScalarMappable(norm=norm,cmap=cmap),
        orientation='vertical',
        shrink=0.5,
        ax=ax)
    
    plt.savefig('Berry_curvature_maps_pU'+str(pU)+'_pW'+str(pW)+
                '_alpha'+str(alpha)+'.png')
    
    ##### ================================================================================
    ###    Plot the 2D maps of the absolute value of the Berry curvature of the 8 bands 
    ##### ================================================================================
    X,Y = np.meshgrid(k_array,q_array)
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

    ax[0,0].set_ylabel('q',fontsize=14)
    ax[1,0].set_ylabel('q',fontsize=14)
    
    fig.colorbar(cm.ScalarMappable(norm=norm,cmap=cmap),
        orientation='vertical',
        shrink=0.5,
        ax=ax)
    plt.savefig('Abs_Berry_curvature_maps_pU'+str(pU)+'_pW'+str(pW)+
                '_alpha'+str(alpha)+'.png')
    
    ##### ================================================================================
    ###                Plot the dispersion surfaces with Berry curvature 
    ##### ================================================================================
    ### All the 8 bands 
    ### ----------------------------------------------------------------------------------
    F_array_3D = F_array 
    maxabs = abs(F_array_3D).max() 
    maxabs = maxabs 
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
        
    ax.set_xlabel('k',fontsize=14)
    ax.set_ylabel('q',fontsize=14)
    ax.set_zlabel('E',fontsize=14)
    ax.set_title(r'$\Delta U/U = $'+str(pU)+', $\Delta W/W = $'+str(pW)
                 +r', $\alpha = $'+str(alpha),fontsize=14)
    fig.colorbar(scamap,
                 orientation='vertical',
                 shrink=0.4,
                 ax = ax)
    ax.view_init(elev=5,azim=75,roll=0)
    plt.savefig('Bands_pU'+str(pU)+'_pW'+str(pW)+'_alpha'+str(alpha)+'.png')        
    
    
    ### -------------------------------------------------------------------------------
    ### Plot the dispersion surfaces of couples of bands 
    ### -------------------------------------------------------------------------------
    for i in range(2,3):
        fig,ax = plt.subplots(subplot_kw = {'projection':'3d'},
                          figsize=(12,10)) 
        """fcolors1 = scamap.to_rgba(F_array_3D[:,:,i].T) 
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
                        cmap=cmap)"""
        ax.plot_surface(X,Y,Energy_array[:,:,i].T,cmap='spring') 
        ax.plot_surface(X,Y,Energy_array[:,:,i+1].T,cmap='summer')
        ax.set_xlabel('k',fontsize=14)
        ax.set_ylabel('q',fontsize=14)
        ax.set_zlabel('E',fontsize=14)
        ax.set_title('Bands '+str(i+1)+'+'+str(i+2)+r'; $\Delta U/U = $'+str(pU)+
                     ', $\Delta W/W = $'+str(pW)+r', $\alpha = $'+str(alpha),
                     fontsize=14)
        fig.colorbar(scamap,
                     orientation='vertical',
                     shrink=0.4,
                     ax = ax)
        ax.view_init(elev=5,azim=75,roll=0)
        plt.savefig('Bands_'+str(i+1)+'_'+str(i+2)+'_pU'+str(pU)+'_pW'+str(pW)+
                    '_alpha'+str(alpha)+'.png')     

    ### ----------------------------------------------------------------------------------
    ### Select a band select = n (2<= n <= 7) and plot it together with bands n-1 and n+1 
    ### ----------------------------------------------------------------------------------
    select = 3
    fig,ax = plt.subplots(subplot_kw = {'projection':'3d'},
                          figsize=(12,10)) 
    fcolors1 = scamap.to_rgba(F_array_3D[:,:,select-2].T) 
    ax.plot_surface(X[:,xmin:xmax],
                    Y[:,xmin:xmax],
                    Energy_array[xmin:xmax,:,select-2].T,
                    linewidth=linewidth,
                    antialiased='True',
                    rstride=1,
                    cstride=1,
                    facecolors=fcolors1,
                    cmap=cmap)
    fcolors2 = scamap.to_rgba(F_array_3D[:,:,select-1].T) 
    ax.plot_surface(X[:,xmin:xmax],
                    Y[:,xmin:xmax],
                    Energy_array[xmin:xmax,:,select-1].T,
                    linewidth=linewidth,
                    antialiased='True',
                    rstride=1,
                    cstride=1,
                    facecolors=fcolors2,
                    cmap=cmap)
    fcolors3 = scamap.to_rgba(F_array_3D[:,:,select].T) 
    ax.plot_surface(X[:,xmin:xmax],
                    Y[:,xmin:xmax],
                    Energy_array[xmin:xmax,:,select].T,
                    linewidth=linewidth,
                    antialiased='True',
                    rstride=1,
                    cstride=1,
                    facecolors=fcolors3,
                    cmap=cmap)
    ax.set_xlabel('k',fontsize=14)
    ax.set_ylabel('q',fontsize=14)
    ax.set_zlabel('E',fontsize=14)
    ax.set_title('Bands '+str(select-1)+'+'+str(select)+'+'+str(select+1)+';'
                 +r'$\Delta U/U = $'+str(pU)+', $\Delta W/W = $'+str(pW)
                 +r', $\alpha = $'+str(alpha),fontsize=14)
    fig.colorbar(scamap,
                 orientation='vertical',
                 shrink=0.4,
                 ax = ax)
    ax.view_init(elev=5,azim=75,roll=0)
    plt.savefig('Bands_'+str(select-1)+'_'+str(select)+'_'+str(select+1)
                +'_pU'+str(pU)+'_pW'+str(pW)+'_alpha'+str(alpha)+'.png')     
    

    ### Check for the dispersion of couple of bands select and select+1  
    ### (1 <= select <= 7)
    ### k-direction with varying value of q 
    select = 3  

    fig,ax = plt.subplots(1,2,figsize=(10,8),sharey=True)
    ax[0].plot(k_array,Energy_array[:,100,select-1],color='red',label='q = '+str(round(q_array[100],6))) 
    ax[0].plot(k_array,Energy_array[:,100,select],color='red')
    #ax[0].plot(k_array,Energy_array[:,145,select-1],color='green',label='q = '+str(round(q_array[145],6))) 
    #ax[0].plot(k_array,Energy_array[:,145,select],color='green')
    #ax[0].plot(k_array,Energy_array[:,180,select-1],color='blue',label='q = '+str(round(q_array[180],6)))
    #ax[0].plot(k_array,Energy_array[:,180,select],color='blue')
    ax[0].set_xlabel('k',fontsize=14)
    ax[0].set_ylabel('E',fontsize=14) 
    ax[0].legend(fontsize=14,loc='center',bbox_to_anchor=(0.275,-0.625,0.5,1))
    ax[0].set_title('Bands '+str(select)+'+'+str(select+1)+r', $\alpha = $'+str(alpha),fontsize=14)

    ax[1].plot(q_array,Energy_array[100,:,select-1],color=[0.5,0,0],label='k = '+str(round(k_array[100],6))) 
    ax[1].plot(q_array,Energy_array[100,:,select],color=[0.5,0,0])
    #ax[1].plot(q_array,Energy_array[145,:,select-1],color=[0,0.5,0],label='k = '+str(round(k_array[145],6))) 
    #ax[1].plot(q_array,Energy_array[145,:,select],color=[0,0.5,0])
    #ax[1].plot(q_array,Energy_array[180,:,select-1],color=[0,0,0.5],label='k = '+str(round(k_array[180],6)))
    #ax[1].plot(q_array,Energy_array[180,:,select],color=[0,0,0.5])
    ax[1].set_xlabel('q',fontsize=14)
    ax[1].legend(fontsize=14,loc='center',bbox_to_anchor=(0.275,-0.625,0.5,1)) 
    ax[1].set_title('Bands '+str(select)+'+'+str(select+1)+r', $\alpha = $'+str(alpha),fontsize=14)
    

    plt.yticks(fontsize=14)
    plt.savefig('Bands_'+str(select)+'_'+str(select+1)+'_dispersion_k_q.png')

    ##### ============================================================================
    ### Show the figures        
    plt.show()
##### ================================================================================ 

if __name__ == "__main__":
    main()