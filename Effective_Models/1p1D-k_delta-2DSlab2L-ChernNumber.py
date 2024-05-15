import numpy as np 
import scipy 
import scipy.linalg as sla 
import cmath 
import matplotlib.pyplot as plt 
from matplotlib import cm, colors 

##### FUNCTION: Hamiltonian of 2D photonic crystal slab bilayer at M-point
###   The slabs are different from each other
def Hamiltonian(q,delta,omega1,v1,U1,W1,omega2,v2,U2,W2,V):
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
            [omega1pp, W1, W1, U1, Vpp, 0.0, 0.0, 0.0],
            [W1, omega1pm, U1, W1, 0.0, Vpm, 0.0, 0.0],
            [W1, U1, omega1mp, W1, 0.0, 0.0, Vmp, 0.0],
            [U1, W1, W1, omega1mm, 0.0, 0.0, 0.0, Vmm],
            [Vmm, 0.0, 0.0, 0.0, omega2pp, W2, W2, U2],
            [0.0, Vmp, 0.0, 0.0, W2, omega2pm, U2, W2],
            [0.0, 0.0, Vpm, 0.0, W2, U2, omega2mp, W2],
            [0.0, 0.0, 0.0, Vpp, U2, W2, W2, omega2mm] 
        ]
    )

    return Hamiltonian  

##### FUNCTION: The derivative of the Hamiltonian with respect to q 
def dH_q(v1,v2):
    dHdq = np.array(
        [
        [v1*np.sqrt(2.0), 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [0.0,             0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 
        [0.0,             0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0,-v1*np.sqrt(2.0), 0.0, 0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0,   0.0, v2*np.sqrt(2.0),0.0,0.0, 0.0],
        [0.0, 0.0, 0.0,   0.0,           0.0, 0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0,   0.0,           0.0, 0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0,   0.0, 0.0, 0.0, 0.0,-v2*np.sqrt(2.0)]
        ]
    )

    return dHdq 

##### FUNCTION: The derivative of the Hamiltonian with respect to delta 
def dH_delta(delta,V):
    K = 2.0*np.pi 

    Vpp = V * cmath.exp(-1j*K*delta)
    Vmm = V * cmath.exp(1j*K*delta) 

    dHddelta = np.array(
        [
            [0.0, 0.0, 0.0, 0.0, -1j*K*Vpp, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0,      0.0,  0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0,      0.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0,  0.0, 0.0, 0.0, 1j*K*Vmm],
            [1j*K*Vmm, 0.0, 0.0,  0.0, 0.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0,       0.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0,       0.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0,-1j*K*Vpp,  0.0, 0.0, 0.0, 0.0] 
        ]
    )

    return dHddelta

    
##### The MAIN program goes here
def main():
    omega1 = 0.300464408
    omega2 = 0.300464408
    v1 = 0.35
    v2 = 0.35

    U = -0.016615673 
    DeltaU = 0.001
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

    ### The array of intrinsic momenta
    Nk = 201
    Kmax = 0.5
    k_array = np.linspace(-Kmax,Kmax,Nk)
    dk = (k_array.max() - k_array.min())/(Nk-1)

    ### The array of delta 
    Ndelta = 201 
    delta_array = np.linspace(-0.5,0.5,Ndelta)
    ddelta = (delta_array.max() - delta_array.min())/(Ndelta-1)

    ### The arrays of energy 
    Energy1 = np.zeros((Nk,Ndelta))
    Energy2 = np.zeros((Nk,Ndelta))
    Energy3 = np.zeros((Nk,Ndelta))
    Energy4 = np.zeros((Nk,Ndelta))
    Energy5 = np.zeros((Nk,Ndelta))
    Energy6 = np.zeros((Nk,Ndelta))
    Energy7 = np.zeros((Nk,Ndelta))
    Energy8 = np.zeros((Nk,Ndelta))

    ### The arrays of Berry curvature
    F_array1 = np.zeros((Nk,Ndelta))
    F_array2 = np.zeros((Nk,Ndelta))
    F_array3 = np.zeros((Nk,Ndelta))
    F_array4 = np.zeros((Nk,Ndelta))
    F_array5 = np.zeros((Nk,Ndelta))
    F_array6 = np.zeros((Nk,Ndelta))
    F_array7 = np.zeros((Nk,Ndelta))
    F_array8 = np.zeros((Nk,Ndelta)) 

    ### We scan over the arrays of genuine and synthetic momenta 
    for i in range(Nk):
        for j in range(Ndelta):
            q = k_array[i]
            delta = delta_array[j]

            H = Hamiltonian(q,delta,omega1,v1,U1,W1,omega2,v2,U2,W2,V)
            dHq = dH_q(v1,v2)
            dHdelta = dH_delta(delta,V)  

            # The i-th element of E is the i-th energy eigenvalue 
            # The j-th element of states is the eigenstate corresponding to E[j]
            E,states = sla.eigh(H) 

            # We save the energy eigenvalues to the arrays Energy1,2,3,4,5,6,7,8 
            Energy1[i,j] = E[0]
            Energy2[i,j] = E[1]
            Energy3[i,j] = E[2]
            Energy4[i,j] = E[3]
            Energy5[i,j] = E[4]
            Energy6[i,j] = E[5]
            Energy7[i,j] = E[6]
            Energy8[i,j] = E[7]

            ### ATTENTION! The formula to evaluate the Berry curvature for band n (1<=n<=8) is:
            #
            #  F_{q delta}^n = \sum_{m \ne n} (-2)*<n| dHq |m><m | dHdelta |n>/(En-Em)^2
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
            dHqe = np.matmul((states.conjugate()).transpose(), np.matmul(dHq,states))
            dHdeltae = np.matmul((states.conjugate()).transpose(),np.matmul(dHdelta,states)) 
   
            for m in range(8):
                if (m != 0):
                    val = -2*np.imag(dHqe[0,m]*dHdeltae[m,0]) / (E[0]-E[m])**2
                    F_array1[i,j] = F_array1[i,j] + val 

                if (m != 1):
                    val = -2*np.imag(dHqe[1,m]*dHdeltae[m,1]) / (E[1]-E[m])**2
                    F_array2[i,j] = F_array2[i,j] + val 

                if (m != 2):
                    val = -2*np.imag(dHqe[2,m]*dHdeltae[m,2]) / (E[2]-E[m])**2
                    F_array3[i,j] = F_array3[i,j] + val 

                if (m != 3):
                    val = -2*np.imag(dHqe[3,m]*dHdeltae[m,3]) / (E[3]-E[m])**2
                    F_array4[i,j] = F_array4[i,j] + val

                if (m != 4):
                    val = -2*np.imag(dHqe[4,m]*dHdeltae[m,4]) / (E[4]-E[m])**2
                    F_array5[i,j] = F_array5[i,j] + val 

                if (m != 5):
                    val = -2*np.imag(dHqe[5,m]*dHdeltae[m,5]) / (E[5]-E[m])**2
                    F_array6[i,j] = F_array6[i,j] + val 

                if (m != 6):
                    val = -2*np.imag(dHqe[6,m]*dHdeltae[m,6]) / (E[6]-E[m])**2
                    F_array7[i,j] = F_array7[i,j] + val 

                if (m != 7):
                    val = -2*np.imag(dHqe[7,m]*dHdeltae[m,7]) / (E[7]-E[m])**2
                    F_array8[i,j] = F_array8[i,j] + val 

    ### Calculate the Chern numbers 
    C1 = np.sum(F_array1)*dk*ddelta/(2.0*np.pi)
    C2 = np.sum(F_array2)*dk*ddelta/(2.0*np.pi)
    C3 = np.sum(F_array3)*dk*ddelta/(2.0*np.pi)
    C4 = np.sum(F_array4)*dk*ddelta/(2.0*np.pi)
    C5 = np.sum(F_array5)*dk*ddelta/(2.0*np.pi)
    C6 = np.sum(F_array6)*dk*ddelta/(2.0*np.pi)
    C7 = np.sum(F_array7)*dk*ddelta/(2.0*np.pi)
    C8 = np.sum(F_array8)*dk*ddelta/(2.0*np.pi)

    print('# Chern number C1 = '+str(C1))
    print('# Chern number C2 = '+str(C2))
    print('# Chern number C3 = '+str(C3))
    print('# Chern number C4 = '+str(C4))
    print('# Chern number C5 = '+str(C5))
    print('# Chern number C6 = '+str(C6))
    print('# Chern number C7 = '+str(C7))
    print('# Chern number C8 = '+str(C8))

    ### Plot the 2D maps of the Berry curvature of the 8 bands 
    X,Y = np.meshgrid(k_array+0.5,delta_array)
    cmap = 'coolwarm' 

    maxabs = np.max( [abs(F_array1).max(), abs(F_array2).max(),
                      abs(F_array3).max(), abs(F_array4).max(),
                      abs(F_array5).max(), abs(F_array6).max(),
                      abs(F_array7).max(), abs(F_array8).max(),
                      ] )
    
    vmin,vmax = -maxabs,maxabs   
    norm = colors.Normalize(vmin=vmin,vmax=vmax)

    fig,ax = plt.subplots()
    ax.pcolormesh(X,Y,F_array1.T,shading='gouraud',cmap=cmap)
    ax.set_xlabel(r'$ka/(2\pi)$',fontsize=14)
    ax.set_ylabel(r'$\delta$',fontsize=14)
    ax.set_title('Band 1')
    fig.colorbar(cm.ScalarMappable(norm=norm,cmap=cmap),ax = ax)
    plt.savefig('Berry_curvature_Band_1.png')

    fig,ax = plt.subplots()
    ax.pcolormesh(X,Y,F_array2.T,shading='gouraud',cmap=cmap)
    ax.set_xlabel(r'$ka/(2\pi)$',fontsize=14)
    ax.set_ylabel(r'$\delta$',fontsize=14)
    ax.set_title('Band 2')
    fig.colorbar(cm.ScalarMappable(norm=norm,cmap=cmap),ax = ax)
    plt.savefig('Berry_curvature_Band_2.png')

    fig,ax = plt.subplots()
    ax.pcolormesh(X,Y,F_array3.T,shading='gouraud',cmap=cmap)
    ax.set_xlabel(r'$ka/(2\pi)$',fontsize=14)
    ax.set_ylabel(r'$\delta$',fontsize=14)
    ax.set_title('Band 3')
    fig.colorbar(cm.ScalarMappable(norm=norm,cmap=cmap),ax = ax)
    plt.savefig('Berry_curvature_Band_3.png')

    fig,ax = plt.subplots()
    ax.pcolormesh(X,Y,F_array4.T,shading='gouraud',cmap=cmap)
    ax.set_xlabel(r'$ka/(2\pi)$',fontsize=14)
    ax.set_ylabel(r'$\delta$',fontsize=14)
    ax.set_title('Band 4')
    fig.colorbar(cm.ScalarMappable(norm=norm,cmap=cmap),ax = ax)
    plt.savefig('Berry_curvature_Band_4.png')

    fig,ax = plt.subplots()
    ax.pcolormesh(X,Y,F_array5.T,shading='gouraud',cmap=cmap)
    ax.set_xlabel(r'$ka/(2\pi)$',fontsize=14)
    ax.set_ylabel(r'$\delta$',fontsize=14)
    ax.set_title('Band 5')
    fig.colorbar(cm.ScalarMappable(norm=norm,cmap=cmap),ax = ax)
    plt.savefig('Berry_curvature_Band_5.png')

    fig,ax = plt.subplots()
    ax.pcolormesh(X,Y,F_array6.T,shading='gouraud',cmap=cmap)
    ax.set_xlabel(r'$ka/(2\pi)$',fontsize=14)
    ax.set_ylabel(r'$\delta$',fontsize=14)
    ax.set_title('Band 6')
    fig.colorbar(cm.ScalarMappable(norm=norm,cmap=cmap),ax = ax)
    plt.savefig('Berry_curvature_Band_6.png')

    fig,ax = plt.subplots()
    ax.pcolormesh(X,Y,F_array7.T,shading='gouraud',cmap=cmap)
    ax.set_xlabel(r'$ka/(2\pi)$',fontsize=14)
    ax.set_ylabel(r'$\delta$',fontsize=14)
    ax.set_title('Band 7')
    fig.colorbar(cm.ScalarMappable(norm=norm,cmap=cmap),ax = ax)
    plt.savefig('Berry_curvature_Band_7.png')

    fig,ax = plt.subplots()
    ax.pcolormesh(X,Y,F_array8.T,shading='gouraud',cmap=cmap)
    ax.set_xlabel(r'$ka/(2\pi)$',fontsize=14)
    ax.set_ylabel(r'$\delta$',fontsize=14)
    ax.set_title('Band 8')
    fig.colorbar(cm.ScalarMappable(norm=norm,cmap=cmap),ax = ax)
    plt.savefig('Berry_curvature_Band_8.png')

    ### Plot the 2D maps of the absolute value of the Berry curvature of the 8 bands
    vmin,vmax = 0.0,maxabs 
    norm = colors.Normalize(vmin=vmin,vmax=vmax)
    cmap = 'bone'

    fig,ax = plt.subplots()
    ax.pcolormesh(X,Y,abs(F_array1.T),shading='gouraud',cmap=cmap)
    ax.set_xlabel(r'$ka/(2\pi)$',fontsize=14)
    ax.set_ylabel(r'$\delta$',fontsize=14)
    ax.set_title('Band 1')
    fig.colorbar(cm.ScalarMappable(norm=norm,cmap=cmap),ax = ax)
    plt.savefig('Abs_Berry_curvature_Band_1.png')

    fig,ax = plt.subplots()
    ax.pcolormesh(X,Y,abs(F_array2.T),shading='gouraud',cmap=cmap)
    ax.set_xlabel(r'$ka/(2\pi)$',fontsize=14)
    ax.set_ylabel(r'$\delta$',fontsize=14)
    ax.set_title('Band 2')
    fig.colorbar(cm.ScalarMappable(norm=norm,cmap=cmap),ax = ax)
    plt.savefig('Abs_Berry_curvature_Band_2.png')

    fig,ax = plt.subplots()
    ax.pcolormesh(X,Y,abs(F_array3.T),shading='gouraud',cmap=cmap)
    ax.set_xlabel(r'$ka/(2\pi)$',fontsize=14)
    ax.set_ylabel(r'$\delta$',fontsize=14)
    ax.set_title('Band 3')
    fig.colorbar(cm.ScalarMappable(norm=norm,cmap=cmap),ax = ax)
    plt.savefig('Abs_Berry_curvature_Band_3.png')

    fig,ax = plt.subplots()
    ax.pcolormesh(X,Y,abs(F_array4.T),shading='gouraud',cmap=cmap)
    ax.set_xlabel(r'$ka/(2\pi)$',fontsize=14)
    ax.set_ylabel(r'$\delta$',fontsize=14)
    ax.set_title('Band 4')
    fig.colorbar(cm.ScalarMappable(norm=norm,cmap=cmap),ax = ax)
    plt.savefig('Abs_Berry_curvature_Band_4.png')

    fig,ax = plt.subplots()
    ax.pcolormesh(X,Y,abs(F_array5.T),shading='gouraud',cmap=cmap)
    ax.set_xlabel(r'$ka/(2\pi)$',fontsize=14)
    ax.set_ylabel(r'$\delta$',fontsize=14)
    ax.set_title('Band 5')
    fig.colorbar(cm.ScalarMappable(norm=norm,cmap=cmap),ax = ax)
    plt.savefig('Abs_Berry_curvature_Band_5.png')

    fig,ax = plt.subplots()
    ax.pcolormesh(X,Y,abs(F_array6.T),shading='gouraud',cmap=cmap)
    ax.set_xlabel(r'$ka/(2\pi)$',fontsize=14)
    ax.set_ylabel(r'$\delta$',fontsize=14)
    ax.set_title('Band 6')
    fig.colorbar(cm.ScalarMappable(norm=norm,cmap=cmap),ax = ax)
    plt.savefig('Abs_Berry_curvature_Band_6.png')

    fig,ax = plt.subplots()
    ax.pcolormesh(X,Y,abs(F_array7.T),shading='gouraud',cmap=cmap)
    ax.set_xlabel(r'$ka/(2\pi)$',fontsize=14)
    ax.set_ylabel(r'$\delta$',fontsize=14)
    ax.set_title('Band 7')
    fig.colorbar(cm.ScalarMappable(norm=norm,cmap=cmap),ax = ax)
    plt.savefig('Abs_Berry_curvature_Band_7.png')

    fig,ax = plt.subplots()
    ax.pcolormesh(X,Y,abs(F_array8.T),shading='gouraud',cmap=cmap)
    ax.set_xlabel(r'$ka/(2\pi)$',fontsize=14)
    ax.set_ylabel(r'$\delta$',fontsize=14)
    ax.set_title('Band 8')
    fig.colorbar(cm.ScalarMappable(norm=norm,cmap=cmap),ax = ax)
    plt.savefig('Abs_Berry_curvature_Band_8.png')
    
    ### Plot the dispersion surface with Berry curvature 
    F_array_1_3D = np.sinh(20.0*F_array1/maxabs) 
    F_array_2_3D = np.sinh(20.0*F_array2/maxabs)
    F_array_3_3D = np.sinh(20.0*F_array3/maxabs) 
    F_array_4_3D = np.sinh(20.0*F_array4/maxabs)
    F_array_5_3D = np.sinh(20.0*F_array5/maxabs) 
    F_array_6_3D = np.sinh(20.0*F_array6/maxabs)
    F_array_7_3D = np.sinh(20.0*F_array7/maxabs) 
    F_array_8_3D = np.sinh(20.0*F_array8/maxabs)
    vmin,vmax = -1,1 
    norm = colors.Normalize(vmin=vmin,vmax=vmax)
    cmap = 'coolwarm'
    scamap = plt.cm.ScalarMappable(norm=norm,cmap=cmap)

    linewidth = 0 

    fig,ax = plt.subplots(subplot_kw={'projection':'3d'})

    fcolors1 = scamap.to_rgba(F_array_1_3D.T) 
    ax.plot_surface(X,Y,Energy1.T,linewidth=linewidth,
                    facecolors=fcolors1,cmap=cmap)
    
    fcolors2 = scamap.to_rgba(F_array_2_3D.T) 
    ax.plot_surface(X,Y,Energy2.T,linewidth=linewidth,
                    facecolors=fcolors2,cmap=cmap)
    
    fcolors3 = scamap.to_rgba(F_array_3_3D.T)
    ax.plot_surface(X,Y,Energy3.T,linewidth=linewidth,
                    facecolors=fcolors3,cmap=cmap)
    
    fcolors4 = scamap.to_rgba(F_array_4_3D.T)
    ax.plot_surface(X,Y,Energy4.T,linewidth=linewidth,
                    facecolors=fcolors4,cmap=cmap)
    
    fcolors5 = scamap.to_rgba(F_array_5_3D.T) 
    ax.plot_surface(X,Y,Energy5.T,linewidth=linewidth,
                    facecolors=fcolors5,cmap=cmap)
    
    fcolors6 = scamap.to_rgba(F_array_6_3D.T) 
    ax.plot_surface(X,Y,Energy6.T,linewidth=linewidth,
                    facecolors=fcolors6,cmap=cmap) 
    
    fcolors7 = scamap.to_rgba(F_array_7_3D.T)
    ax.plot_surface(X,Y,Energy7.T,linewidth=linewidth,
                    facecolors=fcolors7,cmap=cmap)
    
    fcolors8 = scamap.to_rgba(F_array_8_3D.T) 
    ax.plot_surface(X,Y,Energy8.T,linewidth=linewidth,
                    facecolors=fcolors8,cmap=cmap)

    ax.set_xlabel(r'$k a / (2\pi)$',fontsize=14)
    ax.set_ylabel(r'$\delta$',fontsize=14)
    fig.colorbar(scamap,
                 orientation='vertical',
                 shrink=0.4,
                 ax = ax)
    ax.view_init(elev=5,azim=45,roll=0)
    plt.savefig('Bands.png')

    plt.show()
    
if __name__ == "__main__":
    main()