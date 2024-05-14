import numpy as np 
import scipy 
import scipy.linalg as sla 
import matplotlib.pyplot as plt 
from matplotlib import cm,colors 
from mpl_toolkits.mplot3d import Axes3D 
from matplotlib.colors import LightSource 

def Hamiltonian(kx,ky,U,V,Delta,v):
    diag = U - V + 0.5*v**2 * kx**2/U 
    xi = -v*Delta/U 
    Hamiltonian = np.array(
        [
            [diag, xi*kx-1j*V*ky],
            [xi*kx+1j*V*ky,-diag]  
        ]
    )

    return Hamiltonian 

def dH_x(kx,U,Delta,v):
    diag = v**2 * kx/U 
    xi = -v*Delta/U 

    dHx = np.array(
        [
            [diag,xi],
            [xi,-diag] 
        ]
    )

    return dHx 

def dH_y(V):
    dHy = np.array(
        [
            [0,-1j*V],
            [1j*V,0]
        ]
    )

    return dHy 

def main():
    v = 0.5
    U = 0.02 
    V = 1.5*U 
    Delta = 0.1*U 

    Nx = 201
    Ny = 201
    Kx = np.linspace(-0.5,0.5,Nx)
    Ky = np.linspace(-0.5,0.5,Ny)
    dx = (Kx.max()-Kx.min())/(Nx-1)
    dy = (Ky.max()-Ky.min())/(Ny-1)

    Energy_array = np.zeros((Nx,Ny,2))
    F_array = np.zeros((Nx,Ny,2))

    for i in range(Nx):
        for j in range(Ny):
            kx = Kx[i]
            ky = Ky[j]

            H = Hamiltonian(kx,ky,U,V,Delta,v) 
            dHx = dH_x(kx,U,Delta,v) 
            dHy = dH_y(V) 

            # The i-th element of E is the i-th energy eigenvalue
            # The j-th column of states is the eigenstate corresponding to E[j]
            E,states = sla.eigh(H) 

            # We save the energy eigenvalues to the array Energy_array 
            Energy_array[i,j,:] = E  

            ### ATTENTION! The formula to evaluate the Berry curvature is:
            #
            # F_{xy}^n = \sum_{m \ne n} (-2)*<n| dHx |m><m| dHy |n>/(En-Em)^2
            #
            # In fact: <n| dHx |m> and <m| dHy |n> are the matrix elements of 
            # the operators dHx and dHy in the basis of the energy eigenstates
            # of the Hamiltonian
            # 
            # Therefore, we reexpress the matrices dHx and dHy in the basis of 
            # the eigenstates. The transformation is done by the formula:
            #
            #   A' = states^{\dagger}*A*states
            #
            # here A = dHx or dHy
            # and the j-th column of states is the eigenvector coresponding 
            # to the j-th eigenvalue (j is not the loop index in this program)
            dHxe = np.matmul((states.conjugate()).transpose(), np.matmul(dHx,states))
            dHye = np.matmul((states.conjugate()).transpose(), np.matmul(dHy,states))

            for n in range(2):
                for m in range(2):
                    if (m != n):
                        val = -2*np.imag(dHxe[n,m]*dHye[m,n]) / (E[n]-E[m])**2 
                        F_array[i,j,n] = F_array[i,j,n] + val 
            
    C1 = np.sum(np.sum(F_array[:,:,0]))*dx*dy/(2.0*np.pi)
    print('Chern number C1 = '+str(C1)) 

    C2 = np.sum(np.sum(F_array[:,:,1]))*dx*dy/(2.0*np.pi)
    print('Chern number C2 = '+str(C2))

    """
    ### Plot the 2D maps of the Berry curvature of the 2 bands 
    X,Y = np.meshgrid(Kx,Ky)
    #print('X = ')
    #print(X)

    #print('Y = ')
    #print(Y)
 
    cmap = 'coolwarm'

    maxabs = abs(F_array).max()
    vmin,vmax = -maxabs,maxabs 
    norm = colors.Normalize(vmin=vmin,vmax=vmax)

    for i in range(2):
        fig,ax = plt.subplots()
        ax.pcolormesh(X,Y,F_array[:,:,i].T,shading='gouraud',cmap=cmap) 
        ax.set_xlabel(r'$k_x$',fontsize=14)
        ax.set_ylabel(r'$k_y$',fontsize=14)
        ax.set_title('Band '+str(i+1))
        fig.colorbar(cm.ScalarMappable(norm=norm,cmap=cmap),ax=ax)
        plt.savefig('Berry_curvature_Band_'+str(i+1)+'.png')

    ### Plot the dispersion surface with Berry curvature 
    maxabs = abs(F_array).max()
    F_array_3D = np.sinh(5*F_array/maxabs)
    vmin, vmax = -1, 1
    norm = colors.Normalize(vmin=vmin,vmax=vmax)
    scamap = plt.cm.ScalarMappable(norm=norm,cmap=cmap)

    linewidth=0

    fig,ax = plt.subplots(subplot_kw={'projection':'3d'})

    fcolors1 = scamap.to_rgba(F_array_3D[:,:,0].T)
    ax.plot_surface(X,Y,Energy_array[:,:,0].T,linewidth=linewidth,
                    facecolors=fcolors1,cmap=cmap)
    
    fcolors2 = scamap.to_rgba(F_array_3D[:,:,1].T)
    ax.plot_surface(X,Y,Energy_array[:,:,1].T,linewidth=linewidth,
                    facecolors=fcolors2,cmap=cmap) 
    
    ax.set_xlabel(r'$k_x a / (2 \pi)$',fontsize=14)
    ax.set_ylabel(r'$k_y a / (2 \pi)$',fontsize=14)
    fig.colorbar(scamap,
                 orientation='vertical',
                 shrink=0.4,
                 ax = ax)
    ax.view_init(elev=5,azim=120,roll=0)
    plt.savefig('Bands.png')

    plt.show()
    """
    
if __name__ == '__main__':
    main()