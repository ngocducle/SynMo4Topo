import numpy as np 
import scipy 
import scipy.linalg as sla 
import matplotlib.pyplot as plt 

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

def mul(V1,M,V2):
    Prod = np.matmul((V1.conjugate()).transpose(),np.matmul(M,V2))

    return Prod 

def Berry_curvature(dH1,dH2,energy,states):
    E1 = energy[0]
    E2 = energy[1] 

    V1 = states[:,0]
    V2 = states[:,1]

    F = np.zeros(2) 
 
    F[0] = -2*np.imag(mul(V1,dH1,V2)*mul(V2,dH2,V1))/(E2-E1)**2 
    F[1] = -2*np.imag(mul(V2,dH1,V1)*mul(V1,dH2,V2))/(E1-E2)**2 

    return F  

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

    F_array = np.zeros((Nx,Ny,2))

    
    for i in range(Nx):
        for j in range(Ny):
            kx = Kx[i]
            ky = Ky[j]

            H = Hamiltonian(kx,ky,U,V,Delta,v) 
            dHx = dH_x(kx,U,Delta,v) 
            dHy = dH_y(V) 

            energy,states = sla.eigh(H) 

            F_array[i,j,:] = Berry_curvature(dHx,dHy,energy,states)
            
    C1 = np.sum(np.sum(F_array[:,:,0]))*dx*dy/(2.0*np.pi)
    print('Chern number C1 = '+str(C1)) 

    C2 = np.sum(np.sum(F_array[:,:,1]))*dx*dy/(2.0*np.pi)
    print('Chern number C2 = '+str(C2))

    """
    X,Y = np.meshgrid(Kx,Ky)
    #print('X = ')
    #print(X)

    #print('Y = ')
    #print(Y)
 
    fig,ax = plt.subplots(subplot_kw={'projection':'3d'})
    ax.plot_surface(X,Y,F_array[:,:,0].T)
    ax.plot_surface(X,Y,F_array[:,:,1].T)
    ax.set_xlabel('kx',fontsize=14)
    ax.set_ylabel('ky',fontsize=14)
    plt.show()
    """

if __name__ == '__main__':
    main()