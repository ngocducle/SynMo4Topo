import numpy as np 
import scipy 
import scipy.linalg as sla
import cmath 
import matplotlib.pyplot as plt
from matplotlib import cm, colors 
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

##### FUNCTION: calculate the product <V1|M|V2> 
def mul(V1,M,V2):
    Prod = np.matmul((V1.conjugate()).transpose(),np.matmul(M,V2))

    return Prod 

##### FUNCTION: function to calculate the Berry curvature 
def Berry_curvature_q_delta(dH1,dH2,E,V):
    F = np.zeros(4)

    for n in range(4):
        for m in range(4):
            if (m != n):
                F[n] = F[n] -2*np.imag(mul(V[:,n],dH1,V[:,m])*mul(V[:,m],dH2,V[:,n])) \
                    / (E[n]-E[m])**2

    return F 

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

    ### Array of genuine momenta 
    Nk = 251
    Kmax = 0.5
    k_array = np.linspace(-Kmax,Kmax,Nk)

    ### Array of intrinsic momenta
    Ndelta = 201
    delta_array = np.linspace(-0.5,0.5,Ndelta) 

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

            energy,states = sla.eigh(H)

            Energy_array[i,j,:] = energy 
            F_array[i,j,:] = Berry_curvature_q_delta(dHq,dHdelta,energy,states)

    #print(F_array)

    ### Calculate the Chern numbers 
    Chern_number = np.sum(F_array,axis = (0,1))*2.0*Kmax/(Nk*Ndelta*2.0*np.pi)

    print('# Chern numbers = ')    
    print(Chern_number)

    ### Plot the figure 
    X,Y = np.meshgrid(k_array+0.5,delta_array) 

    maxabs = abs(F_array).max() 
    vmin, vmax = -maxabs, maxabs 
    norm = colors.Normalize(vmin=vmin,vmax=vmax)

    cmap = 'rainbow'

    """
    for i in range(4):
        fig,ax = plt.subplots()
        ax.pcolormesh(X,Y,Berry_curvature[:,:,i].T,shading='gouraud',cmap=cmap)
        ax.set_xlabel('q',fontsize=14)
        ax.set_ylabel(r'$\delta$',fontsize=14)
        ax.set_title('Band '+str(i+1))
        fig.colorbar(cm.ScalarMappable(norm=norm,cmap=cmap),ax = ax)
        plt.savefig('Berry_curvature_Band_'+str(i+1)+'.png')
    """
    fig,ax = plt.subplots()
    ax.pcolormesh(X,Y,F_array[:,:,0].T,shading='gouraud',cmap=cmap)
    ax.set_xlabel('q',fontsize=14)
    ax.set_ylabel(r'$\delta$',fontsize=14)
    ax.set_title('Band 1')
    fig.colorbar(cm.ScalarMappable(norm=norm,cmap=cmap),ax = ax)
    plt.savefig('Berry_curvature_Band_1.png')

    fig,ax = plt.subplots()
    ax.pcolormesh(X,Y,F_array[:,:,1].T,shading='gouraud',cmap=cmap)
    ax.set_xlabel('q',fontsize=14)
    ax.set_ylabel(r'$\delta$',fontsize=14)
    ax.set_title('Band 2')
    fig.colorbar(cm.ScalarMappable(norm=norm,cmap=cmap),ax = ax)
    plt.savefig('Berry_curvature_Band_2.png')

    fig,ax = plt.subplots()
    ax.pcolormesh(X,Y,F_array[:,:,2].T,shading='gouraud',cmap=cmap)
    ax.set_xlabel('q',fontsize=14)
    ax.set_ylabel(r'$\delta$',fontsize=14)
    ax.set_title('Band 3')
    fig.colorbar(cm.ScalarMappable(norm=norm,cmap=cmap),ax = ax)
    plt.savefig('Berry_curvature_Band_3.png')

    fig,ax = plt.subplots()
    ax.pcolormesh(X,Y,F_array[:,:,3].T,shading='gouraud',cmap=cmap)
    ax.set_xlabel('q',fontsize=14)
    ax.set_ylabel(r'$\delta$',fontsize=14)
    ax.set_title('Band 4')
    fig.colorbar(cm.ScalarMappable(norm=norm,cmap=cmap),ax = ax)
    plt.savefig('Berry_curvature_Band_4.png')

    ### Plot the dispersion surface with Berry curvature 
    rgb = np.tile([1.0,1.0,1.0],(Ndelta,Nk,1))

    """
    ls = LightSource(azdeg=225,altdeg=30)
    illuminated_surface1 = ls.shade_rgb(rgb,Berry_curvature[:,:,0].T)
    illuminated_surface2 = ls.shade_rgb(rgb,Berry_curvature[:,:,1].T)
    illuminated_surface3 = ls.shade_rgb(rgb,Berry_curvature[:,:,2].T)
    illuminated_surface4 = ls.shade_rgb(rgb,Berry_curvature[:,:,3].T)

    fig,ax = plt.subplots(subplot_kw={'projection':'3d'})    
    ax.plot_surface(X,Y,Energy[:,:,0].T,linewidth=0.01,
                    facecolors=illuminated_surface1,cmap=cmap)
    ax.plot_surface(X,Y,Energy[:,:,1].T,linewidth=0.01,
                    facecolors=illuminated_surface2,cmap=cmap)
    ax.plot_surface(X,Y,Energy[:,:,2].T,linewidth=0.01,
                    facecolors=illuminated_surface3,cmap=cmap)
    ax.plot_surface(X,Y,Energy[:,:,3].T,linewidth=0.01,
                    facecolors=illuminated_surface4,cmap=cmap)
    ax.set_xlabel(r'$k a / (2 \ pi)$',fontsize=14)
    ax.set_ylabel(r'$\delta$',fontsize=14)
    fig.colorbar(cm.ScalarMappable(norm=norm,cmap=cmap),
                 orientation='vertical',
                 shrink=0.4,
                 ax = ax)
    plt.savefig('Bands.png')    
    plt.show()
    """

    fig,ax = plt.subplots(subplot_kw={'projection':'3d'}) 
    scamap = plt.cm.ScalarMappable(norm=norm,cmap=cmap)
    fcolors1 = scamap.to_rgba(F_array[:,:,0].T)   
    ax.plot_surface(X,Y,Energy_array[:,:,0].T,linewidth=1e-7,
                    facecolors=fcolors1,cmap=cmap)

    fcolors2 = scamap.to_rgba(F_array[:,:,1].T)  
    ax.plot_surface(X,Y,Energy_array[:,:,1].T,linewidth=1e-7,
                    facecolors=fcolors2,cmap=cmap)
    
    fcolors3 = scamap.to_rgba(F_array[:,:,2].T)  
    ax.plot_surface(X,Y,Energy_array[:,:,2].T,linewidth=1e-7,
                    facecolors=fcolors3,cmap=cmap)
    
    fcolors4 = scamap.to_rgba(F_array[:,:,3].T)  
    ax.plot_surface(X,Y,Energy_array[:,:,3].T,linewidth=1e-7,
                    facecolors=fcolors4,cmap=cmap)
    
    ax.set_xlabel(r'$k a / (2 \ pi)$',fontsize=14)
    ax.set_ylabel(r'$\delta$',fontsize=14)
    fig.colorbar(scamap,
                 orientation='vertical',
                 shrink=0.4,
                 ax = ax)
    plt.savefig('Bands.png')    
    plt.show()

if __name__ == '__main__':
    main()