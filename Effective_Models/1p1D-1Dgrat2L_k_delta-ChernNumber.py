import numpy as np 
import scipy 
import scipy.linalg as sla
import cmath 
import matplotlib.pyplot as plt
from matplotlib import cm,colors 
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.colors import LightSource 

##### FUNCTION: Hamiltonian of 1D grating bilayer at X-point
###   The slabs are different from each other
def Hamiltonian(q,delta,v,U1,U2,V):
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
def dH_q(v):
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
def dH_delta(delta,V):
    Vp = V*cmath.exp(-1j*np.pi*delta) 
    Vm = V*cmath.exp(1j*np.pi*delta)

    dHddelta = np.array(
        [
            [0.0,0.0,-1j*np.pi*Vp,0.0],
            [0.0,0.0,0.0,1j*np.pi*Vm],
            [1j*np.pi*Vm,0.0,0.0,0.0],
            [0.0,-1j*np.pi*Vp,0.0,0.0] 
        ]
    )

    return dHddelta 

##### The MAIN program goes here
def main():
    U = 0.03
    Delta = -0.05*U
    Gamma = -0.05*U  
    U1 = U+Delta
    U2 = U-Delta  
    V = np.sqrt(U**2-Delta**2)-Gamma 
    ng = 3.1024 
    v = 1.0/(2.0*np.pi*ng)

    ### Array of genuine momenta 
    Nk = 201
    Kmax = 0.5
    k_array = np.linspace(-Kmax,Kmax,Nk)
    dk = (k_array.max() - k_array.min())/(Nk-1)

    ### Array of intrinsic momenta
    Ndelta = 201
    delta_array = np.linspace(-0.5,0.5,Ndelta) 
    ddelta = (delta_array.max() - delta_array.min())/(Ndelta-1)

    ### Arrays of energy 
    Energy_array = np.zeros((Nk,Ndelta,4))

    ### Array of Berry curvature 
    F_array = np.zeros((Nk,Ndelta,4))

    ### Array of Chern number 
    Chern_number = np.zeros(4)

    ### Calculate the Berry curvature 
    for i in range(Nk):
        for j in range(Ndelta):
            q = k_array[i]
            delta = delta_array[j] 

            H = Hamiltonian(q,delta,v,U1,U2,V)
            dHq = dH_q(v)
            dHdelta = dH_delta(delta,V) 

            # The i-th element of E is the i-th energy eigenvalue
            # The j-th column of states is the eigenstate corresponding to E[j] 
            E,states = sla.eigh(H)

            # We save the energy eigenvalues to the array Energy_array 
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
            dHqe = np.matmul( (states.conjugate()).transpose(), np.matmul(dHq,states) )
            dHdeltae = np.matmul((states.conjugate()).transpose(), np.matmul(dHdelta,states)) 

            #print(dHqe)
            #print(dHdeltae)

            for n in range(4):
                for m in range(4):
                    if (m != n):
                        val = -2*np.imag(dHqe[n,m]*dHdeltae[m,n]) / (E[n]-E[m])**2 
                        F_array[i,j,n] = F_array[i,j,n] + val 

    #print(Energy_array)
    #print(F_array)

    ### Calculate the Chern numbers 
    Chern_number = np.sum(F_array,axis = (0,1))*dk*ddelta/(2.0*np.pi)

    print('# Chern number C1 = '+str(Chern_number[0]))
    print('# Chern number C2 = '+str(Chern_number[1]))
    print('# Chern number C3 = '+str(Chern_number[2]))
    print('# Chern number C4 = '+str(Chern_number[3]))    

    ### Print the results to a file
    for n in range(4): 
        with open('1p1D-1Dgrat2L-k-delta-Farray-'+str(n+1)+'.txt','w') as file:
            for i in range(Nk):
                for j in range(Ndelta):
                    file.write('%.8f    ' % k_array[i])
                    file.writelines('%.8f    ' % delta_array[j])
                    file.writelines('%.8f    ' % F_array[i,j,n])
                    file.write('\n')

    ### Plot the 2D maps of the Berry curvature of the 4 bands 
    X,Y = np.meshgrid(k_array+0.5,delta_array) 
    cmap = 'coolwarm'

    maxabs = abs(F_array).max() 
    vmin, vmax = -maxabs, maxabs 
    norm = colors.Normalize(vmin=vmin,vmax=vmax)

    for i in range(4):
        fig,ax = plt.subplots()
        ax.pcolormesh(X,Y,F_array[:,:,i].T,shading='gouraud',cmap=cmap)
        ax.set_xlabel('q',fontsize=14)
        ax.set_ylabel(r'$\delta$',fontsize=14)
        ax.set_title('Band '+str(i+1))
        fig.colorbar(cm.ScalarMappable(norm=norm,cmap=cmap),ax = ax)
        plt.savefig('Berry_curvature_Band_'+str(i+1)+'.png')

    
    ### Plot the 2D maps of the absolute value of the Berry curvature of the 4 bands 
    maxabs = abs(F_array).max() 
    vmin, vmax = 0.0, maxabs 
    norm = colors.Normalize(vmin=vmin,vmax=vmax)

    for i in range(4):
        fig,ax = plt.subplots()
        ax.pcolormesh(X,Y,abs(F_array[:,:,i]).T,shading='gouraud',cmap='bone')
        ax.set_xlabel('q',fontsize=14)
        ax.set_ylabel(r'$\delta$',fontsize=14)
        ax.set_title('Band '+str(i+1))
        fig.colorbar(cm.ScalarMappable(norm=norm,cmap='bone'),ax = ax)
        plt.savefig('Abs_Berry_curvature_Band_'+str(i+1)+'.png')

    ### Plot the dispersion surface with Berry curvature     
    F_array_3D = F_array
    maxabs = abs(F_array_3D).max()
    vmin, vmax = -maxabs,maxabs  
    norm = colors.Normalize(vmin=vmin,vmax=vmax)
    scamap = plt.cm.ScalarMappable(norm=norm,cmap=cmap)

    ### Check for the symmetry of F_array_3D
    #print(F_array_3D.max())
    #print(F_array_3D.min())

    #for i in range(Nk):
    #    print(str(F_array_3D[i,99,2])+','
    #         +str(F_array_3D[i,100,2])+','
    #         +str(F_array_3D[i,101,2]))

    linewidth=0 

    fig,ax = plt.subplots(subplot_kw={'projection':'3d'})
    
    fcolors1 = scamap.to_rgba(F_array_3D[:,:,0].T)   
    ax.plot_surface(X,Y,Energy_array[:,:,0].T,linewidth=linewidth,
                    antialiased='True',rstride=5,cstride=5,
                    facecolors=fcolors1,cmap=cmap)

    fcolors2 = scamap.to_rgba(F_array_3D[:,:,1].T)  
    print(np.shape(fcolors2)) 
    ax.plot_surface(X,Y,Energy_array[:,:,1].T,linewidth=linewidth,
                    antialiased='True',rstride=1,cstride=5,
                    facecolors=fcolors2,cmap=cmap)
    
    fcolors3 = scamap.to_rgba(F_array_3D[:,:,2].T)  
    ax.plot_surface(X,Y,Energy_array[:,:,2].T,linewidth=linewidth,
                    antialiased='True',rstride=1,cstride=5,
                    facecolors=fcolors3,cmap=cmap)
    
    fcolors4 = scamap.to_rgba(F_array_3D[:,:,3].T)  
    ax.plot_surface(X,Y,Energy_array[:,:,3].T,linewidth=linewidth,
                    antialiased='True',rstride=5,cstride=5,
                    facecolors=fcolors4,cmap=cmap)
    
    ax.set_xlabel(r'$k a / (2 \ pi)$',fontsize=14)
    ax.set_ylabel(r'$\delta$',fontsize=14)
    fig.colorbar(scamap,
                 orientation='vertical',
                 shrink=0.4,
                 ax = ax)
    ax.view_init(elev=5,azim=45,roll=0)
    plt.savefig('Bands.png')    

    plt.show()
    
if __name__ == '__main__':
    main()