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
            [0.0,0.0,0.0,1j*np.pi*Vm],
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

    ### Array of genuine momenta 
    Nk = 251
    Kmax = 0.5
    k_array = np.linspace(-Kmax,Kmax,Nk)

    ### Array of intrinsic momenta
    Ndelta = 201
    delta_array = np.linspace(-0.5,0.5,Ndelta) 

    ### Arrays of energy 
    Energy = np.zeros((Nk,Ndelta,4))

    ### Array of Berry curvature 
    Berry_curvature = np.zeros((Nk,Ndelta,4))

    ### Array of Chern number 
    Chern_number = np.zeros(4)

    ### Calculate the Berry curvature 
    for i in range(Nk):
        for j in range(Ndelta):
            q = k_array[i]
            delta = delta_array[j] 
            evalues,Berry_curvature[i,j,:] = Berry_curvature_q_delta(q,delta,v,U1,U2,V) 
            Energy[i,j,:] = evalues 

    print(Berry_curvature)

    ### Calculate the Chern numbers 
    Chern_number = np.sum(Berry_curvature,axis = (0,1))*2.0*Kmax/(Nk*Ndelta*2.0*np.pi)

    print('# Chern numbers = ')    
    print(Chern_number)

    ### Plot the figure 
    X,Y = np.meshgrid(k_array+0.5,delta_array) 

    maxabs = abs(Berry_curvature).max() 
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
    ax.pcolormesh(X,Y,Berry_curvature[:,:,0].T,shading='gouraud',cmap=cmap)
    ax.set_xlabel('q',fontsize=14)
    ax.set_ylabel(r'$\delta$',fontsize=14)
    ax.set_title('Band 1')
    fig.colorbar(cm.ScalarMappable(norm=norm,cmap=cmap),ax = ax)
    plt.savefig('Berry_curvature_Band_1.png')

    fig,ax = plt.subplots()
    ax.pcolormesh(X,Y,Berry_curvature[:,:,1].T,shading='gouraud',cmap=cmap)
    ax.set_xlabel('q',fontsize=14)
    ax.set_ylabel(r'$\delta$',fontsize=14)
    ax.set_title('Band 2')
    fig.colorbar(cm.ScalarMappable(norm=norm,cmap=cmap),ax = ax)
    plt.savefig('Berry_curvature_Band_2.png')

    fig,ax = plt.subplots()
    ax.pcolormesh(X,Y,Berry_curvature[:,:,2].T,shading='gouraud',cmap=cmap)
    ax.set_xlabel('q',fontsize=14)
    ax.set_ylabel(r'$\delta$',fontsize=14)
    ax.set_title('Band 3')
    fig.colorbar(cm.ScalarMappable(norm=norm,cmap=cmap),ax = ax)
    plt.savefig('Berry_curvature_Band_3.png')

    fig,ax = plt.subplots()
    ax.pcolormesh(X,Y,Berry_curvature[:,:,3].T,shading='gouraud',cmap=cmap)
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
    fcolors1 = scamap.to_rgba(Berry_curvature[:,:,0].T)   
    ax.plot_surface(X,Y,Energy[:,:,0].T,linewidth=1e-7,
                    facecolors=fcolors1,cmap=cmap)

    fcolors2 = scamap.to_rgba(Berry_curvature[:,:,1].T)  
    ax.plot_surface(X,Y,Energy[:,:,1].T,linewidth=1e-7,
                    facecolors=fcolors2,cmap=cmap)
    
    fcolors3 = scamap.to_rgba(Berry_curvature[:,:,2].T)  
    ax.plot_surface(X,Y,Energy[:,:,2].T,linewidth=1e-7,
                    facecolors=fcolors3,cmap=cmap)
    
    fcolors4 = scamap.to_rgba(Berry_curvature[:,:,3].T)  
    ax.plot_surface(X,Y,Energy[:,:,3].T,linewidth=1e-7,
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